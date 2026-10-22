"""ChEMBL ligands for a UniProt target (get_ligands).

Written against responses validated on 2026-10-17 for KRAS / P01116 (target,
activity, molecule) and 2026-10-21 for EGFR (assay, document, mechanism).
Fields are read defensively with .get(); if the API does not return them the
output says so instead of filling gaps.

No project imports on purpose (avoids cycles). Only httpx (dependency of mcp).
"""
import asyncio
import logging
import re
import statistics
from collections import defaultdict

import httpx

log = logging.getLogger("protein-mcp.chembl")  # never print() to stdout

BASE = "https://www.ebi.ac.uk/chembl/api/data"
PAGE = 1000                 # ChEMBL page size used when paging activities
MAX_ACTIVITIES = 6000       # hard cap on activities fetched per call
MAX_ASSAYS = 20000          # EGFR has 5966 assays (2026-10-21); a 5000 cap lost a third of them
MAX_MECHANISMS = 2000       # EGFR has 90 mechanism rows (2026-10-21)
MAX_LIGANDS = 15            # default number of ligands returned
MIN_CONFIDENCE = 8          # assay confidence_score threshold (client-side)
STANDARD_TYPES = "IC50,Ki,Kd"
IN_BATCH = 50               # ids per __in query
TIMEOUT = 60.0
MAX_DOCS_PER_GROUP = 5
MAX_DESCS_PER_GROUP = 2

ASSAY_TYPES = {
    "B": "binding", "F": "functional", "A": "ADME",
    "T": "toxicity", "P": "physicochemical", "U": "unassigned",
}
NOT_STATED = "not stated"


# ---------------------------------------------------------------- helpers
def _num(x):
    try:
        return None if x is None else float(x)
    except (TypeError, ValueError):
        return None


def _short(text, n=160):
    text = " ".join(str(text or "").split())
    return text if len(text) <= n else text[: n - 3] + "..."


async def _get(client, endpoint, params):
    url = f"{BASE}/{endpoint}.json"
    last = None
    for attempt in range(3):
        try:
            r = await client.get(url, params=params)
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
            else:
                r.raise_for_status()
                return r.json()
        except httpx.HTTPError as e:
            last = f"{type(e).__name__}: {e}"
        await asyncio.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"ChEMBL {endpoint} failed: {last}")


async def _get_all(client, endpoint, params, key, cap):
    """Page through an endpoint. Returns (items, total_count)."""
    first = await _get(client, endpoint, {**params, "limit": PAGE, "offset": 0})
    items = list(first.get(key) or [])
    total = int((first.get("page_meta") or {}).get("total_count") or len(items))
    offsets = list(range(PAGE, min(total, cap), PAGE))
    sem = asyncio.Semaphore(4)

    async def one(off):
        async with sem:
            d = await _get(client, endpoint, {**params, "limit": PAGE, "offset": off})
            return d.get(key) or []

    for chunk in await asyncio.gather(*(one(o) for o in offsets)):
        items.extend(chunk)
    return items, total


async def _safe(coro, label):
    try:
        return (await coro), None
    except Exception as e:  # isolated failure per source
        log.warning("%s failed: %s", label, e)
        return None, f"{label}: {e}"


def _batches(ids, n=IN_BATCH):
    ids = list(ids)
    return [ids[i:i + n] for i in range(0, len(ids), n)]


def _mut_tokens(value):
    return [t for t in re.split(r"[/,;+\s]+", (value or "").upper()) if t]


def _mut_matches(requested, value):
    """Exact variant token, or a bare position ("G12") matching any substitution."""
    req = requested.strip().upper()
    toks = _mut_tokens(value)
    if re.fullmatch(r"[A-Z]\d+", req):
        return any(re.fullmatch(re.escape(req) + r"[A-Z]", t) for t in toks)
    return req in toks


# ------------------------------------------------------------ target step
def _pick_target(targets, uniprot_id):
    chosen, ignored = None, []
    for t in targets:
        comps = t.get("target_components") or []
        is_single = (
            t.get("target_type") == "SINGLE PROTEIN"
            and len(comps) == 1
            and comps[0].get("accession") == uniprot_id
        )
        if is_single and chosen is None:
            chosen = t
        else:
            ignored.append({
                "chembl_id": t.get("target_chembl_id"),
                "type": t.get("target_type"),
                "name": t.get("pref_name"),
            })
    return chosen, ignored


# ------------------------------------------------------- detail fetchers
async def _fetch_assays(client, target_id):
    """Assay records of the target, keyed by assay id. Returns (map, total_count)."""
    items, total = await _get_all(
        client, "assay", {"target_chembl_id": target_id}, "assays", MAX_ASSAYS)
    return {
        a["assay_chembl_id"]: {
            "confidence": a.get("confidence_score"),
            "cell_type": a.get("assay_cell_type"),
        }
        for a in items if a.get("assay_chembl_id")
    }, total


async def _fetch_mechanisms(client, target_id):
    """Mechanisms keyed by PARENT molecule. Verified 2026-10-21 (EGFR): 15 of 90
    rows are recorded on a salt (erlotinib HCl CHEMBL1079742) while activities use
    the parent (CHEMBL553); keying by the salt hid approved drugs.
    Returns (map, rows_read, rows_total)."""
    items, total = await _get_all(client, "mechanism", {"target_chembl_id": target_id},
                                  "mechanisms", MAX_MECHANISMS)
    out = {}
    for m in items:
        mid = m.get("parent_molecule_chembl_id") or m.get("molecule_chembl_id")
        if not mid:
            continue
        cur = out.setdefault(mid, {"mechanism_of_action": None, "action_type": None,
                                   "max_phase": None})
        cur["mechanism_of_action"] = cur["mechanism_of_action"] or m.get("mechanism_of_action")
        cur["action_type"] = cur["action_type"] or m.get("action_type")
        phase = _num(m.get("max_phase"))
        if phase is not None and (cur["max_phase"] is None or phase > cur["max_phase"]):
            cur["max_phase"] = phase
    return out, len(items), total


async def _fetch_activities_for(client, base_params, parent_ids):
    """All qualifying activities of the given parent molecules (no cap per batch
    beyond MAX_ACTIVITIES), so drugs keep full statistics when the ranked
    download is truncated. Returns (activities, batches_capped): a batch with more
    than MAX_ACTIVITIES rows is cut and counted."""
    out, capped = [], 0
    for batch in _batches(parent_ids):
        items, total = await _get_all(
            client, "activity",
            {**base_params, "parent_molecule_chembl_id__in": ",".join(batch)},
            "activities", MAX_ACTIVITIES)
        out.extend(items)
        capped += len(items) < total
    return out, capped


async def _fetch_in(client, endpoint, key, id_field, ids):
    out = {}
    for batch in _batches(ids):
        d = await _get(client, endpoint,
                       {f"{id_field}__in": ",".join(batch), "limit": len(batch)})
        for item in d.get(key) or []:
            out[item[id_field]] = item
    return out


# ------------------------------------------------------------ aggregation
_DESC_MUT_RE = re.compile(r"\b[A-Z]\d{1,3}[A-Z]\b")
MAX_FLAGGED = 10


def _label_conflicts_description(label, description):
    """True when the assay text names variants and none is the curated label.

    ChEMBL's assay_variant_mutation can disagree with the assay description
    (seen on KRAS: label G12C, text 'KRAS G12D mutant'). Such rows are excluded
    and listed, never silently trusted.
    """
    named = set(_DESC_MUT_RE.findall((description or "").upper()))
    if not named or not label:
        return False
    return not (named & set(_mut_tokens(label)))


def _assay_class(assay, assay_type):
    if assay.get("cell_type"):
        return "cellular"
    return ASSAY_TYPES.get(assay_type, assay_type or "unknown")


def _aggregate(acts, uniprot_id, mutation, assays, min_conf):
    excluded = defaultdict(int)
    conf_dist = defaultdict(int)
    flagged = []
    groups = {}
    seen = set()
    for a in acts:
        aid = a.get("activity_id")
        if aid in seen:
            continue
        seen.add(aid)

        mut_raw = a.get("assay_variant_mutation")
        if mutation and not _mut_matches(mutation, mut_raw):
            excluded["other_or_unstated_mutation"] += 1
            continue
        if a.get("potential_duplicate") in (1, "1", True):
            excluded["potential_duplicate"] += 1
            continue
        if a.get("data_validity_comment"):
            excluded["data_validity_flagged"] += 1
            continue
        if a.get("standard_units") != "nM":
            excluded["units_not_nM"] += 1
            continue
        value, pchembl = _num(a.get("standard_value")), _num(a.get("pchembl_value"))
        if value is None or pchembl is None:
            excluded["missing_value_or_pchembl"] += 1
            continue

        assay = (assays or {}).get(a.get("assay_chembl_id"))
        if assays and assay is None:
            # Without its assay record the confidence filter and the cellular /
            # biochemical split cannot be applied: excluded and counted, not guessed.
            excluded["assay_record_missing"] += 1
            continue
        assay = assay or {}
        conf = assay.get("confidence")
        conf_dist[str(conf)] += 1
        if assays and min_conf and conf is not None and conf < min_conf:
            excluded["below_min_confidence"] += 1
            continue
        if a.get("standard_relation") != "=":
            excluded["censored_relation_not_equal"] += 1
            continue

        mol = a.get("parent_molecule_chembl_id") or a.get("molecule_chembl_id")
        if not mol:
            continue
        if _label_conflicts_description(mut_raw, a.get("assay_description")):
            excluded["mutation_label_conflicts_description"] += 1
            if len(flagged) < MAX_FLAGGED:
                flagged.append({
                    "activity_id": aid, "molecule_chembl_id": mol,
                    "labelled_mutation": mut_raw,
                    "assay_chembl_id": a.get("assay_chembl_id"),
                    "assay_description": _short(a.get("assay_description")),
                })
            continue
        mut = (mut_raw or "").strip() or NOT_STATED
        # Cellular and biochemical potencies are never pooled in one range.
        aclass = _assay_class(assay, a.get("assay_type"))
        key = (mol, mut, a.get("standard_type"), aclass)
        g = groups.setdefault(key, {
            "nM": [], "pchembl": [], "assays": set(), "docs": {}, "types": set(),
            "cells": set(), "descs": [], "smiles": a.get("canonical_smiles"),
        })
        g["nM"].append(value)
        g["pchembl"].append(pchembl)
        g["assays"].add(a.get("assay_chembl_id"))
        if assay.get("cell_type"):
            g["cells"].add(assay["cell_type"])
        d = a.get("assay_description")
        if d and d not in g["descs"] and len(g["descs"]) < MAX_DESCS_PER_GROUP:
            g["descs"].append(d)
        doc = a.get("document_chembl_id")
        if doc and doc not in g["docs"]:
            g["docs"][doc] = (a.get("document_journal"), a.get("document_year"))
    return groups, dict(excluded), dict(conf_dist), flagged


def _rank(groups):
    by_mol = defaultdict(lambda: {"best": 0.0, "n": 0})
    for (mol, *_rest), g in groups.items():
        by_mol[mol]["best"] = max(by_mol[mol]["best"], max(g["pchembl"]))
        by_mol[mol]["n"] += len(g["nM"])
    return sorted(by_mol.items(), key=lambda kv: (-kv[1]["best"], -kv[1]["n"]))


def _phase_rank(phase):
    """Sort key for ChEMBL max_phase: 4, 3, 2, 1, 0.5 (early phase 1), then -1
    (clinical, phase unknown per ChEMBL's schema), then None (preclinical)."""
    if phase is None:
        return 0.0
    return 0.1 if phase < 0 else phase


def _molecule_block(mol_id, m, mech):
    retrieved = m is not None
    m = m or {}
    props = m.get("molecule_properties") or {}
    structs = m.get("molecule_structures") or {}
    return {
        "chembl_id": mol_id,
        "name": m.get("pref_name"),   # often None for research compounds: never invented
        # False when the molecule endpoint failed: name/type/phase are then unknown,
        # not absent from ChEMBL.
        "molecule_record_retrieved": retrieved,
        "link": f"https://www.ebi.ac.uk/chembl/compound_report_card/{mol_id}/",
        "molecule_type": m.get("molecule_type"),
        "max_phase": _num(m.get("max_phase")),   # None = no phase recorded
        "first_approval": m.get("first_approval"),
        "chemical_probe_flag": m.get("chemical_probe"),
        "mechanism": mech,
        "identity": {
            "smiles": structs.get("canonical_smiles"),
            "inchi_key": structs.get("standard_inchi_key"),
            "molecular_weight": _num(props.get("mw_freebase") or props.get("full_mwt")),
        },
        "properties": {   # None when ChEMBL has no value (e.g. macrocyclic peptides)
            "alogp": _num(props.get("alogp")),
            "psa": _num(props.get("psa")),
            "hba": props.get("hba"),
            "hbd": props.get("hbd"),
            "ro5_violations": props.get("num_ro5_violations"),
            "aromatic_rings": props.get("aromatic_rings"),
            "qed_weighted": _num(props.get("qed_weighted")),
        },
    }


def _potency_rows(mol_id, groups, docs_info):
    rows = []
    for (mol, mut, stype, aclass), g in groups.items():
        if mol != mol_id:
            continue
        docs = []
        for doc_id, (journal, year) in list(g["docs"].items())[:MAX_DOCS_PER_GROUP]:
            info = docs_info.get(doc_id) or {}
            links = {"chembl_document": f"https://www.ebi.ac.uk/chembl/document_report_card/{doc_id}/"}
            if info.get("pubmed_id"):
                links["pubmed"] = f"https://pubmed.ncbi.nlm.nih.gov/{info['pubmed_id']}/"
            if info.get("doi"):
                links["doi"] = f"https://doi.org/{info['doi']}"
            docs.append({
                "document_chembl_id": doc_id,
                "pubmed_id": info.get("pubmed_id"),
                "doi": info.get("doi"),
                "journal": info.get("journal") or journal,
                "year": info.get("year") or year,
                "links": links,
            })
        rows.append({
            "mutation": mut,
            "standard_type": stype,
            "assay_class": aclass,   # cellular / binding / functional ...
            "n_measurements": len(g["nM"]),
            "best_nM": round(min(g["nM"]), 3),
            "median_nM": round(statistics.median(g["nM"]), 3),
            "worst_nM": round(max(g["nM"]), 3),
            "best_pchembl": round(max(g["pchembl"]), 2),
            "n_assays": len(g["assays"]),
            "cell_types": sorted(g["cells"]),
            "assay_examples": [_short(d) for d in g["descs"]],
            "documents": docs,
        })
    rows.sort(key=lambda r: -r["best_pchembl"])
    return rows


# ------------------------------------------------------------ public API
_POTENCY_KEEP = ("mutation", "standard_type", "assay_class", "n_measurements",
                 "best_nM", "median_nM", "worst_nM", "best_pchembl", "cell_types")
COMPACT_POTENCY_ROWS = 3    # EGFR: 51 ligands x up to 23 groups = 94 KB without this cap


def compact(result):
    """Chat-sized view of get_ligands. Drops SMILES, properties, assay texts and the
    document list of each potency row (one source link is kept), and keeps only the
    COMPACT_POTENCY_ROWS best potency groups per ligand (potency_groups_total says
    how many exist); the full result goes to the dossier."""
    if not result.get("ok"):
        return result

    def source(p):
        for d in p.get("documents") or []:
            links = d.get("links") or {}
            if links.get("pubmed") or links.get("doi"):
                return links.get("pubmed") or links.get("doi")
        return None

    ligands = []
    for lig in result.get("ligands") or []:
        ligands.append({
            "chembl_id": lig.get("chembl_id"),
            "name": lig.get("name"),
            **({"molecule_record_retrieved": False} if lig.get("molecule_record_retrieved") is False else {}),
            "link": lig.get("link"),
            "max_phase": lig.get("max_phase"),
            "first_approval": lig.get("first_approval"),
            "mechanism_of_action": (lig.get("mechanism") or {}).get("mechanism_of_action"),
            "inchi_key": (lig.get("identity") or {}).get("inchi_key"),
            "best_pchembl": lig.get("best_pchembl"),
            "n_measurements": lig.get("n_measurements"),
            "potency_groups_total": len(lig.get("potency") or []),
            "potency": [   # best groups first; the rest stay in the dossier
                {**{k: p.get(k) for k in _POTENCY_KEEP},
                 "source": source(p), "n_documents": len(p.get("documents") or [])}
                for p in (lig.get("potency") or [])[:COMPACT_POTENCY_ROWS]
            ],
        })
    flagged = result.get("label_conflicts_flagged") or []
    out = {k: v for k, v in result.items() if k not in ("ligands", "label_conflicts_flagged")}
    out["label_conflicts_flagged"] = {   # total excluded: activities.excluded
        "listed_in_dossier": len(flagged),
        "examples": [{k: f.get(k) for k in ("molecule_chembl_id", "labelled_mutation", "assay_chembl_id")}
                     for f in flagged[:3]],
    }
    out["ligands"] = ligands
    return out


async def get_ligands(uniprot_id, mutation="", max_ligands=MAX_LIGANDS,
                      min_confidence=MIN_CONFIDENCE, _client=None):
    uid = (uniprot_id or "").strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{6,10}", uid):
        return {"ok": False, "error": f"'{uniprot_id}' is not a valid UniProt accession."}
    max_ligands = max(1, min(int(max_ligands), 50))
    mutation = (mutation or "").strip()

    own = _client is None
    client = _client or httpx.AsyncClient(timeout=TIMEOUT, follow_redirects=True)
    try:
        try:
            td = await _get(client, "target",
                            {"target_components__accession": uid, "limit": 50})
        except Exception as e:
            return {"ok": False, "uniprot_id": uid, "error": str(e)}
        target, ignored = _pick_target(td.get("targets") or [], uid)
        if not target:
            return {"ok": True, "uniprot_id": uid, "target": None,
                    "other_targets_ignored": ignored,
                    "message": "No single-protein ChEMBL target for this accession: sin registro en ChEMBL."}
        tid = target["target_chembl_id"]

        act_params = {
            "target_chembl_id": tid,
            "standard_type__in": STANDARD_TYPES,
            "pchembl_value__isnull": "false",
        }
        # Most potent first, activity_id as tie-breaker: ordering by pChEMBL alone
        # repeated and skipped rows across pages (verified on EGFR, 2026-10-21).
        ranked_params = {**act_params, "order_by": ["-pchembl_value", "activity_id"]}
        (acts_res, assays_res, mech_res) = await asyncio.gather(
            _safe(_get_all(client, "activity", ranked_params, "activities", MAX_ACTIVITIES),
                  "activities"),
            _safe(_fetch_assays(client, tid), "assays"),
            _safe(_fetch_mechanisms(client, tid), "mechanisms"),
        )
        if acts_res[0] is None:
            return {"ok": False, "uniprot_id": uid, "error": acts_res[1]}
        acts, total = acts_res[0]
        truncated = len(acts) < total
        cutoff = min((_num(a.get("pchembl_value")) for a in acts
                      if _num(a.get("pchembl_value")) is not None), default=None)
        (assays, assays_total), assays_err = (assays_res[0] or (None, None)), assays_res[1]
        mech_data, mech_err = mech_res
        mechs, mech_rows_read, mech_rows_total = mech_data or ({}, None, None)
        # Failed lookup: drugs without IC50/Ki/Kd are unknown, never "none".
        mech_ok = mech_err is None

        mech_acts, mech_acts_err, mech_capped = [], None, 0
        if truncated and mechs:
            mech_fetch, mech_acts_err = await _safe(
                _fetch_activities_for(client, act_params, list(mechs)),
                "activities of mechanism compounds")
            mech_acts, mech_capped = mech_fetch or ([], 0)
        fetched_ids = {a.get("activity_id") for a in acts}
        extra = [a for a in mech_acts if a.get("activity_id") not in fetched_ids]

        groups, excluded, conf_dist, flagged = _aggregate(
            acts + extra, uid, mutation, assays, min_confidence)
        ranked = _rank(groups)
        present = {m for m, _ in ranked}
        # Same filters without the mutation one: a drug with data only on other
        # variants must not be presented as having no IC50/Ki/Kd at all.
        present_any = present if not mutation else {
            k[0] for k in _aggregate(acts + extra, uid, "", assays, min_confidence)[0]}
        mech_unchecked = truncated and mech_acts_err is not None
        top_ids = [m for m, _ in ranked[:max_ligands]]
        clinical_ids = [m for m in mechs if m in present and m not in top_ids]
        detail_ids = top_ids + clinical_ids
        missing_ids = sorted(m for m in mechs if m not in present)
        other_variant_ids = [m for m in missing_ids if m in present_any]
        # Truncated download + failed extra download: below-cutoff data may exist.
        unchecked_ids = [m for m in missing_ids if m not in present_any and mech_unchecked]
        no_activity_ids = [m for m in missing_ids if m not in present_any and not mech_unchecked]

        mol_ids = detail_ids + missing_ids
        mol_map, mol_err = (await _safe(
            _fetch_in(client, "molecule", "molecules", "molecule_chembl_id", mol_ids),
            "molecules")) if mol_ids else ({}, None)
        mol_map = mol_map or {}

        doc_ids = {d for (m, *_rest), g in groups.items() if m in detail_ids
                   for d in list(g["docs"])[:MAX_DOCS_PER_GROUP]}
        docs_raw, doc_err = (await _safe(
            _fetch_in(client, "document", "documents", "document_chembl_id", doc_ids),
            "documents")) if doc_ids else ({}, None)
        docs_info = {k: {"pubmed_id": v.get("pubmed_id"), "doi": v.get("doi"),
                         "journal": v.get("journal"), "year": v.get("year")}
                     for k, v in (docs_raw or {}).items()}

        stats = {m: s for m, s in ranked}
        ligands = []
        for mol_id in detail_ids:
            block = _molecule_block(mol_id, mol_map.get(mol_id), mechs.get(mol_id))
            if block["max_phase"] is None and mechs.get(mol_id):
                block["max_phase"] = mechs[mol_id].get("max_phase")
            block["best_pchembl"] = round(stats[mol_id]["best"], 2)
            block["n_measurements"] = stats[mol_id]["n"]
            block["potency"] = _potency_rows(mol_id, groups, docs_info)
            ligands.append(block)
        # Clinical-stage compounds first (what a discovery analyst checks first),
        # then by best pChEMBL.
        ligands.sort(key=lambda l: (-_phase_rank(l["max_phase"]), -l["best_pchembl"]))

        # Compounds with a recorded mechanism that are not in `ligands`, listed with
        # name and phase so an analyst never reads "no approved drugs".
        def mechanism_entries(ids):
            out = []
            for mol_id in ids:
                m, mech = mol_map.get(mol_id), mechs[mol_id]
                phase = _num((m or {}).get("max_phase"))
                out.append({
                    "chembl_id": mol_id,
                    "name": (m or {}).get("pref_name"),
                    "molecule_record_retrieved": m is not None,
                    "link": f"https://www.ebi.ac.uk/chembl/compound_report_card/{mol_id}/",
                    "molecule_type": (m or {}).get("molecule_type"),
                    "max_phase": phase if phase is not None else mech.get("max_phase"),
                    "first_approval": (m or {}).get("first_approval"),
                    "mechanism_of_action": mech.get("mechanism_of_action"),
                    "action_type": mech.get("action_type"),
                })
            out.sort(key=lambda d: (-_phase_rank(d["max_phase"]), d["chembl_id"]))
            return out

        # No IC50/Ki/Kd at all (antibodies, proteins: evolocumab on PCSK9).
        mechanism_only = mechanism_entries(no_activity_ids)
        # Data exist, but only for other variants than `mutation` (EGFR C797S: lapatinib).
        other_variants_only = mechanism_entries(other_variant_ids)
        # Not checked: their full activity download failed.
        unchecked = mechanism_entries(unchecked_ids)
        if not mech_ok:
            mechanism_only = other_variants_only = unchecked = None

        errors = [e for e in (assays_err, mech_err, mech_acts_err, mol_err, doc_err) if e]
        conf_applied = min_confidence if (assays and min_confidence) else None
        assays_incomplete = assays is not None and assays_total is not None and len(assays) < assays_total

        notes = [
            "Ligands are listed clinical-stage first (max_phase), then by best pChEMBL: potency is not selectivity or developability, and many top values are single measurements (see n_measurements).",
            "ChEMBL max_phase: 4 approved, 3/2/1 clinical phase, 0.5 early phase 1, -1 clinical but phase unknown to ChEMBL, None preclinical.",
            "Potency is grouped by (mutation, standard_type, assay_class); cellular and biochemical values and different mutations are never pooled. 'not stated' does NOT prove wild type.",
            "Activities whose curated mutation label disagrees with the mutation named in the assay text are excluded and listed in label_conflicts_flagged.",
            "Only compounds with a recorded mechanism or inside the top-N are checked for max_phase: other clinical-stage compounds may exist among the ligands not returned.",
            "Only exact '=' relations enter the statistics; censored values (>, <) are counted in activities.excluded.",
            "name is ChEMBL's preferred name and is None for most research compounds; no names are inferred.",
            "Linking PDB ligands to ChEMBL (by InChIKey) is not done here: no match means 'sin registro', never a guess.",
        ]
        if mech_ok:
            notes.insert(2, "Compounds with a recorded mechanism on this target but no qualifying IC50/Ki/Kd (typically antibodies and other proteins) are listed in mechanism_compounds_without_qualifying_activity with name, type and phase: they are drugs too.")
        # State of the extra download of mechanism compounds (only needed when truncated).
        if not truncated:
            mech_download = None
        elif not mech_ok:
            mech_download = "mechanism_lookup_failed"
        elif not mechs:
            mech_download = "no_mechanism_compounds"
        elif mech_acts_err:
            mech_download = "failed"
        elif mech_capped:
            mech_download = "partial"
        else:
            mech_download = "complete"
        if other_variants_only:
            notes.append(
                f"{len(other_variants_only)} compounds with a recorded mechanism have qualifying data only for "
                f"other variants than '{mutation}' (mechanism_compounds_other_variants_only): they are not "
                "untested; they have no qualifying measurement on this variant in ChEMBL (one may exist but "
                "be excluded: censored relation, low confidence, label conflict).")
        if not mech_ok:
            notes.append(
                "The ChEMBL mechanism lookup FAILED: drugs without IC50/Ki/Kd (e.g. antibodies), mechanisms "
                "of action and the phase of compounds outside the top-N were NOT checked. The "
                "mechanism_compounds_* lists are None (not checked), which does not mean there are none.")
        elif mech_rows_total is not None and mech_rows_read < mech_rows_total:
            notes.append(f"Only {mech_rows_read} of {mech_rows_total} mechanism records were read: "
                         "some drugs with a recorded mechanism may be missing.")
        if truncated:
            extra_text = {
                "mechanism_lookup_failed": "; the mechanism lookup failed, so drugs got no extra download",
                "failed": (
                    f"; the extra download of every activity of the {len(mechs)} compounds with a recorded "
                    f"mechanism FAILED: {len(unchecked or [])} of them could not be checked "
                    "(mechanism_compounds_unchecked) and the statistics of the others only count "
                    "measurements above the cutoff"),
                "no_mechanism_compounds": "; no compound has a recorded mechanism on this target",
                "partial": (f", plus the activities of the {len(mechs)} compounds with a recorded mechanism "
                            f"({mech_capped} batch(es) cut at {MAX_ACTIVITIES} activities)"),
                "complete": f", plus every activity of the {len(mechs)} compounds with a recorded mechanism",
            }[mech_download]
            notes.append(
                f"Truncated: the {len(acts)} most potent of {total} activities were fetched "
                f"(lowest pChEMBL fetched {cutoff}; activities tied at that value may be partly missing)"
                f"{extra_text}. ligands_with_qualifying_data is therefore a lower bound. For other compounds, "
                "n_measurements, median and worst values only count measurements above that cutoff, and "
                "weaker ligands are not seen."
                + (" With a mutation filter, activities of that variant below the cutoff are also missed."
                   if mutation else ""))
        if errors:
            notes.append("Partial data, these sources failed: " + "; ".join(errors))
        if mol_err:
            notes.append("Molecule records could not be read: name, type and phase are unknown where "
                         "molecule_record_retrieved is false (not absent from ChEMBL).")
        if assays is None:
            notes.append("Assay confidence could not be read: no confidence filter was applied.")
        if assays_incomplete:
            notes.append(f"Only {len(assays)} of {assays_total} assay records were read; activities of the "
                         "other assays are excluded (activities.excluded.assay_record_missing).")

        return {
            "ok": True,
            "uniprot_id": uid,
            "target": {
                "chembl_id": tid,
                "name": target.get("pref_name"),
                "type": target.get("target_type"),
                "organism": target.get("organism"),
                "link": f"https://www.ebi.ac.uk/chembl/target_report_card/{tid}/",
            },
            "other_targets_ignored": ignored,
            "mutation_filter": mutation or None,
            "activities": {
                "total_in_chembl": total,
                "fetched": len(acts) + len(extra),
                "fetched_for_mechanism_compounds": len(extra),
                "mechanism_compounds_download": mech_download,
                "truncated": truncated,
                "ranked_download_cutoff_pchembl": cutoff if truncated else None,
                "used": sum(len(g["nM"]) for g in groups.values()),
                "excluded": excluded,
            },
            "confidence": {"min_applied": conf_applied, "distribution_before_filter": conf_dist,
                           "assay_records_read": len(assays) if assays is not None else None,
                           "assay_records_total": assays_total},
            "mechanisms": {"retrieved": mech_ok, "records_read": mech_rows_read,
                           "records_total": mech_rows_total},
            "label_conflicts_flagged": flagged,
            "selection": (
                ("At least " if truncated else "")
                + f"{len(ranked)} ligands have qualifying data"
                + (f" for '{mutation}'" if mutation else "")
                + (" (download truncated: only the most potent activities were seen)" if truncated else "")
                + f"; returned the top {len(top_ids)} by best pChEMBL (listed clinical-stage first)"
                + (f" plus {len(clinical_ids)} with a recorded mechanism" if clinical_ids else "")
                + "."
                + (f" {len(mechanism_only)} more compounds with a recorded mechanism have no qualifying"
                   " IC50/Ki/Kd (see mechanism_compounds_without_qualifying_activity)."
                   if mechanism_only else "")
                + (f" {len(other_variants_only)} have data only for other variants"
                   " (see mechanism_compounds_other_variants_only)." if other_variants_only else "")
                + (f" {len(unchecked)} could not be checked (see mechanism_compounds_unchecked)."
                   if unchecked else "")
                + (" WARNING: the ChEMBL mechanism lookup failed, so drugs without IC50/Ki/Kd"
                   " (e.g. antibodies) were NOT checked." if not mech_ok else "")
            ),
            "ligands_with_qualifying_data": len(ranked),
            "ligands_with_qualifying_data_is_lower_bound": truncated,
            "ligands_returned": len(ligands),
            "ligands": ligands,
            "mechanism_compounds_without_qualifying_activity": mechanism_only,
            "mechanism_compounds_other_variants_only": other_variants_only,
            "mechanism_compounds_unchecked": unchecked,
            "notes": notes,
        }
    finally:
        if own:
            await client.aclose()