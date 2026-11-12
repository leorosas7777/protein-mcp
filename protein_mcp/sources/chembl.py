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
import time
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
# MEASURED 2026-11-12, three probes of the same molecule URL: HTTP 200 in 1.8 s, HTTP 500 in
# 63.2 s, HTTP 200 in 1.4 s. A healthy answer takes under 2 s and a failing one can take over a
# minute, so one unlucky request outlived Claude Desktop's own limit on a tool call; with three
# attempts and the back-off in _get, the old flat 60 s allowed 3*60+9 = 189 s per endpoint, and
# analyze_admet calls two or three of them. Measured consequence (Desktop test of 2026-11-12):
# analyze_admet timed out twice with nothing to show. The UNKNOWN branches already existed and
# were honest; what was missing was reaching them quickly.
# The short read timeout belongs to the ADME path, which answers one chat message and has a
# budget. get_ligands keeps the old 60 s on purpose: it pages with limit=1000, one failed page
# kills the whole tool, and pages that END WELL above 20 s were measured on 2026-11-12 (offset
# 5000 of CHEMBL203 returned HTTP 200 after 71.2 s, while offsets 0 and 1000 took 6-7 s). A
# 20 s cap there would turn a slow success into a failure, which is a regression, so the two
# callers get two timeouts and the reason is written here rather than shared by accident.
READ_TIMEOUT = 20.0        # ADME path: ten times a healthy answer, well under a caller's window
CONNECT_TIMEOUT = 10.0
TIMEOUT = 60.0                                                  # get_ligands and its paging
ADME_TIMEOUT = httpx.Timeout(READ_TIMEOUT, connect=CONNECT_TIMEOUT)
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


class Deadline:
    """A wall-clock budget shared by every request of one call.

    httpx bounds a single request; this bounds the TOTAL, retries and back-off included, so a
    slow source does not outlive the caller that is waiting. When it runs out the answer is an
    honest UNKNOWN inside the caller's window instead of a hang with nothing to show.

    The bound is not exact, and the slack is stated rather than implied: a request already in
    flight is given at least 1 s of read time, and CONNECT_TIMEOUT is not reduced by what is
    left, so a call can overshoot the budget by roughly CONNECT_TIMEOUT. It cannot overshoot by
    a multiple of it, which is what the 3 x 60 s of retries did before.
    """

    def __init__(self, budget):
        self.budget = float(budget)
        self._start = time.monotonic()

    def spent(self):
        return round(time.monotonic() - self._start, 1)

    def remaining(self):
        return self.budget - (time.monotonic() - self._start)

    def expired(self):
        return self.remaining() <= 0


async def _get(client, endpoint, params, deadline=None):
    url = f"{BASE}/{endpoint}.json"
    last = None
    for attempt in range(3):
        if deadline is not None and deadline.expired():
            raise RuntimeError(
                f"ChEMBL {endpoint} gave up after {deadline.spent()} s of a "
                f"{deadline.budget} s budget: {last or 'the first attempt did not complete'}")
        kwargs = {}
        if deadline is not None:
            # never wait past the budget for one request either
            kwargs["timeout"] = httpx.Timeout(
                min(READ_TIMEOUT, max(1.0, deadline.remaining())), connect=CONNECT_TIMEOUT)
        try:
            r = await client.get(url, params=params, **kwargs)
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
            else:
                r.raise_for_status()
                return r.json()
        except httpx.HTTPError as e:
            last = f"{type(e).__name__}: {e}"
        nap = 1.5 * (attempt + 1)
        if deadline is not None:
            nap = min(nap, max(0.0, deadline.remaining()))
        await asyncio.sleep(nap)
    raise RuntimeError(f"ChEMBL {endpoint} failed: {last}")


async def _get_all(client, endpoint, params, key, cap, deadline=None):
    """Page through an endpoint. Returns (items, total_count)."""
    first = await _get(client, endpoint, {**params, "limit": PAGE, "offset": 0}, deadline)
    items = list(first.get(key) or [])
    total = int((first.get("page_meta") or {}).get("total_count") or len(items))
    offsets = list(range(PAGE, min(total, cap), PAGE))
    sem = asyncio.Semaphore(4)

    async def one(off):
        async with sem:
            d = await _get(client, endpoint, {**params, "limit": PAGE, "offset": off}, deadline)
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
                    "message": "No single-protein ChEMBL target for this accession: no record in ChEMBL."}
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
            "Linking PDB ligands to ChEMBL (by InChIKey) is not done here: no match means 'no record', never a guess.",
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

# ---------------------------------------------------------------- ADME (Phase 3 step 7)
# Measured against the real API on 2026-11-11 with the two ligands of this project
# (ESTADO 8.5). Three findings shape everything below:
#   1. assay_type=A is NOT ADMET. It holds 1,168,142 rows whose standard_type vocabulary
#      includes 'Saluretic potency', 'Hammett constant', 'Tail vein toxicity' and
#      'No. of animals protected/animals tested'. Of vemurafenib's 17 "A" rows only ONE is
#      ADME (CLint in human hepatocytes); the rest are HepG2 cytotoxicity and a B-Raf
#      western blot. So the filter here is a POSITIVE allowlist, never ChEMBL's own label.
#   2. Units are not homogeneous inside one standard_type: 'CL' arrives as mL.min-1.kg-1,
#      mL.min-1.g-1 and 'ml kg-1' (not a clearance unit at all) in the same 3000-row sample,
#      and 'Kd' arrives in min-1. Nothing here is averaged, aggregated or converted: every
#      row keeps its own units verbatim, with its relation and its assay description.
#   3. The vocabulary is case-inconsistent ('Papp' and 'permeability' in the same molecule),
#      so matching is case-folded while the output prints what ChEMBL returned.
# hERG is deliberately NOT here: it lives in assay_type=F against target CHEMBL240 and is
# off-target safety, not ADME (decision of 2026-11-11; ESTADO 6).
ADME_TYPES = {
    # absorption / permeability / physicochemical
    "papp": "permeability", "permeability": "permeability", "logpapp": "permeability",
    "absorption": "absorption", "solubility": "solubility",
    "logd": "lipophilicity", "logp": "lipophilicity", "pka": "ionisation",
    # distribution
    "vd": "volume of distribution", "vdss": "volume of distribution",
    "ppb": "plasma protein binding", "fu": "fraction unbound",
    "bpr": "tissue-to-plasma ratio",
    # metabolism / clearance
    "cl": "clearance", "clint": "clearance", "clearance": "clearance", "clh": "clearance",
    "t1/2": "half-life", "% remaining": "metabolic stability",
    "remaining": "metabolic stability", "stability": "metabolic stability",
    "drug metabolism": "metabolic stability",
    # exposure (in vivo pharmacokinetics)
    "auc": "exposure", "cmax": "exposure", "tmax": "exposure", "cp": "exposure",
    "auc/dose": "exposure", "cmax/dose": "exposure", "plasma level": "exposure",
    "f": "oral bioavailability", "bioavailability": "oral bioavailability",
    "mrt": "mean residence time",
}
# MEASURED 2026-11-11, second round, after the reviewer found the allowlist too narrow. The nine
# types added above were checked one by one against the real API first, and two candidates were
# REJECTED because their standard_type lies about their content:
#   - 'pKa1'/'pKa2' (57 and 66 rows) carry "Partition coefficient at pH 7.4 between chloroform
#     and water as a measure of lipid solubility", which is a logP, not an acid constant;
#   - 'Solubility ratio' (365 rows) is "[HbS+drug]/[HbS-drug]", the solubility of sickle
#     haemoglobin, not of the compound.
# Those two are why the allowlist is read as a list of types whose CONTENT was checked, not a
# list of names that sound right. Added with their measured evidence: 'CLH' (14,593 rows,
# uL/min, "Hepatic clearance in rat liver microsomes"), 'logPapp' (722, Caco-2 apparent
# permeability), 'Cp' (10,210, ug ml-1, plasma concentration), 'Plasma level' (396),
# 'AUC/dose' (2,290), 'Cmax/dose' (408), 'Stability' (50,706, %, "First pass metabolism and
# metabolic bioavailability using rat hepatic microsomes"), 'Drug metabolism' (27,774,
# "Metabolic stability ... assessed as disappearance of compound") and 'BPR' (3,489, "Ratio of
# drug level in brain against plasma"; named tissue-to-plasma, not blood-to-plasma, because
# that is what the descriptions say).
# STILL NOT ON THE LIST, and this is a declared limitation, not an oversight: the generic types
# 'Activity' and 'Ratio' carry real ADME in some rows (P-gp and BCRP efflux ratios, CYP450
# induction for vemurafenib) mixed with cytotoxicity and target-engagement in others. They
# cannot be separated without reading each assay description, so admitting them would bring
# back the exact defect the allowlist exists to prevent. The excluded rows are therefore
# reported as "not on this project's allowlist", NEVER as "not ADME".
# Words in an assay description that say the depositor CALCULATED the value rather than
# measuring it. Real rows exist: 'pKa' rows read "Acid dissociation constant was calculated
# based on carbon atom CH2". A flagged row is still reported, with the flag, because the claim
# comes from the source's own text.
_CALCULATED_WORDS = ("calculated", "calcd", "predicted", "estimated", "computed", "in silico")
MAX_ADME_ROWS = 300      # vemurafenib has 17 rows of assay_type=A; CHEMBL53463 has 1143
MAX_EXCLUDED_SHOWN = 40  # the excluded rows travel with their description (F3), so they are capped
# Budget for ALL the ChEMBL work of one analyze_admet call (molecule lookup plus activities), so
# the tool reaches its UNKNOWN branch inside the caller's window instead of hanging. A CHOICE,
# not a measurement: Claude Desktop's limit on a tool call is not published, and the measured
# facts behind the number are the ones in the READ_TIMEOUT comment (a healthy answer in under
# 2 s, a failing one in 63.2 s).
ADME_BUDGET = 40.0
# MEASURED 2026-11-11: the page size decides whether this query works at all. With
# assay_type=A and one molecule, limit=1000 returned HTTP 500 three times out of three and
# limit=200 three out of three, while limit=50 and limit=20 returned 200 (with an occasional
# 500, which _get retries). So _get_all's PAGE=1000 cannot be used here.
ADME_PAGE = 20
# The assay types queried for ADME. The loop and the output both read THIS tuple, so the
# output cannot claim a type that was not asked for (reviewer finding, 2026-11-11).
ADME_ASSAY_TYPES = ("A", "P")
# fields of `molecule` that are records of the ChEMBL entry: not calculations, not measurements
RECORD_FIELDS = ("max_phase", "first_approval", "withdrawn_flag", "black_box_warning", "oral",
                 "parenteral", "topical", "molecule_type", "natural_product", "prodrug",
                 "therapeutic_flag", "availability_type", "chirality", "dosed_ingredient",
                 "indication_class", "usan_year", "usan_stem", "first_in_class", "structure_type",
                 "inorganic_flag", "polymer_flag", "orphan", "chemical_probe", "veterinary")
_INCHIKEY_RE = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]")


def adme_class(standard_type):
    """The ADME class of a standard_type, or None if it is not on the allowlist."""
    return ADME_TYPES.get(" ".join(str(standard_type or "").split()).lower())


async def fetch_molecule_by_inchikey(inchi_key: str, client=None, deadline=None) -> dict:
    """Exact 27-character InChIKey -> the one ChEMBL molecule, or found: False.

    Exact match only (decision of 2026-11-11): a connectivity-block match would also hit a
    different salt, tautomer or stereoisomer, and this is the link that lets ADMET data be
    attributed to the compound that was docked and simulated. No compound name is ever used
    to search, and none is ever inferred.
    """
    key = (inchi_key or "").strip().upper()
    if not _INCHIKEY_RE.fullmatch(key):
        return {"found": False, "inchi_key": inchi_key,
                "error": "Not a 27-character standard InChIKey, so ChEMBL was not queried."}
    own = client is None
    client = client or httpx.AsyncClient(timeout=ADME_TIMEOUT, follow_redirects=True)
    try:
        d = await _get(client, "molecule",
                       {"molecule_structures__standard_inchi_key": key}, deadline)
        mols = d.get("molecules") or []
        total = int((d.get("page_meta") or {}).get("total_count") or len(mols))
        if not mols:
            return {"found": False, "inchi_key": key, "matches": 0,
                    "note": "No ChEMBL molecule has this exact InChIKey (no record)."}
        if total > 1:
            # not observed in the two real cases (both total_count=1); never guess which one
            return {"found": False, "inchi_key": key, "matches": total,
                    "candidates": [m.get("molecule_chembl_id") for m in mols[:10]],
                    "note": (f"{total} ChEMBL molecules share this exact InChIKey, so none was "
                             "chosen. The candidate ids are listed.")}
        m = mols[0]
        structs = m.get("molecule_structures") or {}
        cid = m.get("molecule_chembl_id")
        return {
            "found": True,
            "inchi_key": key,
            "matches": 1,
            "chembl_id": cid,
            "link": f"https://www.ebi.ac.uk/chembl/compound_report_card/{cid}/" if cid else None,
            # None for most research compounds; never substituted by a synonym or a guess
            "pref_name": m.get("pref_name"),
            "inchi_key_in_chembl": structs.get("standard_inchi_key"),
            "canonical_smiles_in_chembl": structs.get("canonical_smiles"),
            "record": {f: m.get(f) for f in RECORD_FIELDS},
            # MEASURED 2026-11-11: TXV (CHEMBL5440746) comes back with prodrug, chirality,
            # availability_type, first_in_class, inorganic_flag, orphan and veterinary all at
            # -1, while vemurafenib has 0/1 there (and chirality 2, availability_type 1, so
            # those two are NOT binary flags -- an earlier wording here claimed they were and
            # that claim came from no source). All this says is that ChEMBL returned -1, which
            # is named and NOT interpreted: this project has not read ChEMBL's schema, so it
            # cannot say what -1 means, and least of all read it as "no".
            "record_fields_minus_one": sorted(f for f in RECORD_FIELDS if m.get(f) == -1),
            # The other half of the same asymmetry (reviewer finding): these flags come back
            # False on an entry with nothing annotated (TXV: max_phase None, molecule_type
            # None, no ATC, no cross-references, seven fields at -1 -- and oral/parenteral/
            # topical/therapeutic_flag/dosed_ingredient/withdrawn_flag all False). The API does
            # not distinguish "annotated as no" from "not annotated", so a False here must not
            # be read as a negative finding. Named, not reinterpreted.
            "record_boolean_fields_false": sorted(
                f for f in RECORD_FIELDS if m.get(f) is False),
            "atc_classifications": list(m.get("atc_classifications") or []),
            "cross_references": list(m.get("cross_references") or []),
            # calculated BY ChEMBL with its own engine; kept apart from our RDKit values on purpose
            "chembl_calculated_properties": m.get("molecule_properties"),
            "molecule_hierarchy": m.get("molecule_hierarchy"),
        }
    finally:
        if own:
            await client.aclose()


async def fetch_adme_activities(chembl_id: str, client=None, deadline=None) -> dict:
    """ADME rows for one molecule: assay_type A and P fetched, then the positive allowlist applied.

    Both types are queried because MEASURED 2026-11-11: vemurafenib's only aqueous solubility
    (2.0 ug.mL-1 in FaSSIF at pH 6.5) is filed under assay_type=P, not A, so querying A alone
    silently loses real ADME measurements. B (binding), F (functional), T (toxicity) and U are
    NOT queried: that is where off-target safety lives (hERG included), which is out of this
    step's scope, so a measurement filed only there is missed and the note says so.

    Every surviving row keeps its own units and its assay_type, and the rows that did not survive
    are counted by standard_type, so "few measurements" is never silently "we filtered them out".
    """
    own = client is None
    client = client or httpx.AsyncClient(timeout=ADME_TIMEOUT, follow_redirects=True)
    try:
        rows, totals, out_of_time, stopped_by_clock = [], {}, [], []
        for at in ADME_ASSAY_TYPES:
            params = {"molecule_chembl_id": chembl_id, "assay_type": at}
            total, offset, read_here = None, 0, 0
            while True:
                # MAX_ADME_ROWS is a cap on ROWS READ, not on the offset: with the offset form
                # a molecule could come back with MAX+PAGE-1 rows (reviewer finding)
                room = MAX_ADME_ROWS - len(rows)
                if room <= 0:
                    break
                # out of time is NOT out of rows: a type left unread here is UNKNOWN for its own
                # reason, and the two reasons are reported separately on purpose
                if deadline is not None and deadline.expired():
                    if total is None:
                        out_of_time.append(at)       # never asked about at all
                    else:
                        stopped_by_clock.append(at)  # asked, then cut off mid-type
                    break
                d = await _get(client, "activity",
                               {**params, "limit": min(ADME_PAGE, room), "offset": offset},
                               deadline)
                page = d.get("activities") or []
                if total is None:
                    total = int((d.get("page_meta") or {}).get("total_count") or len(page))
                for a in page:
                    a["_assay_type"] = at
                rows.extend(page)
                offset += len(page)
                read_here += len(page)
                if not page or read_here >= total:
                    break
            # None, not 0, when the cap was already spent and this type was never asked about:
            # "ChEMBL has no rows of type P" and "we never asked for type P" are different
            # statements, and only the first is a fact (found by the truncation test)
            totals[at] = total
        known = [v for v in totals.values() if v is not None]
        total = sum(known)
        unqueried = sorted(at for at, v in totals.items()
                           if v is None and at not in out_of_time)
        kept, dropped, excluded = [], defaultdict(int), []
        for a in rows:
            cls = adme_class(a.get("standard_type"))
            if not cls:
                dropped[a.get("standard_type")] += 1
                # F3 (Desktop test of 2026-11-12): the counts alone left the analyst unable to
                # see WHAT was excluded, and the whole point of the allowlist being narrower than
                # ADME is that these rows may hold real measurements (P-gp and BCRP efflux ratios,
                # CYP450 induction). The row travels with its own assay description so the
                # analyst can see it; it is NOT classified and must never be added to the
                # allowlisted rows or counted as an ADME measurement.
                if len(excluded) < MAX_EXCLUDED_SHOWN:
                    ex_assay = a.get("assay_chembl_id")
                    excluded.append({
                        "on_adme_allowlist": False,
                        "assay_type": a.get("_assay_type"),
                        "assay_type_meaning": ASSAY_TYPES.get(a.get("_assay_type")),
                        "standard_type": a.get("standard_type"),
                        "relation": a.get("standard_relation"),
                        "value": a.get("standard_value"),
                        "units": a.get("standard_units"),
                        "assay_description": _short(a.get("assay_description"), 300),
                        "assay_chembl_id": ex_assay,
                        "assay_link": (f"https://www.ebi.ac.uk/chembl/assay_report_card/{ex_assay}/"
                                       if ex_assay else None),
                    })
                continue
            doc, assay = a.get("document_chembl_id"), a.get("assay_chembl_id")
            desc = _short(a.get("assay_description"), 300)
            # the source's own text, not our judgement: 'pKa' rows exist that read "Acid
            # dissociation constant was calculated based on carbon atom CH2"
            said = [w for w in _CALCULATED_WORDS if w in (desc or "").lower()]
            kept.append({
                "adme_class": cls,
                "description_says_value_was_calculated": bool(said),
                "calculated_wording_found": said or None,
                "assay_type": a.get("_assay_type"),
                "assay_type_meaning": ASSAY_TYPES.get(a.get("_assay_type")),
                "standard_type": a.get("standard_type"),      # as ChEMBL returned it
                "relation": a.get("standard_relation"),       # a censored value stays censored
                "value": a.get("standard_value"),
                "units": a.get("standard_units"),             # never converted
                "assay_description": _short(a.get("assay_description"), 300),
                "assay_chembl_id": assay,
                "assay_link": (f"https://www.ebi.ac.uk/chembl/assay_report_card/{assay}/"
                               if assay else None),
                "target_organism": a.get("target_organism"),
                "bao_label": a.get("bao_label"),
                "data_validity_comment": a.get("data_validity_comment"),
                "potential_duplicate": a.get("potential_duplicate"),
                "activity_comment": a.get("activity_comment"),
                "document_chembl_id": doc,
                "document_link": (f"https://www.ebi.ac.uk/chembl/document_report_card/{doc}/"
                                  if doc else None),
            })
        truncated = total > len(rows) or bool(unqueried) or bool(out_of_time)
        asked = [at for at in ADME_ASSAY_TYPES
                 if at not in unqueried and at not in out_of_time]
        # REVIEWER FINDING 2026-11-12, introduced by F1 the same day: when the budget was spent
        # before any type could be asked about, this returned retrieved=True with rows_total 0,
        # which reads as "ChEMBL holds no ADME measurements" -- the exact shape of TXV's real
        # zero. Nothing was retrieved, so say so and let the caller's UNKNOWN branch handle it.
        if not asked:
            return {
                "retrieved": False,
                "lookup_reached_chembl": True,
                "assay_types_skipped_at_the_time_budget": out_of_time or None,
                "assay_types_skipped_at_the_row_cap": unqueried or None,
                "note": ("No assay type could be queried: the "
                         + ("time budget ran out" if out_of_time else "row cap was already spent")
                         + f" first ({', '.join(out_of_time or unqueried)}). Whether ADME "
                         "measurements exist is UNKNOWN. This is NOT 'no data' and NOT zero "
                         "rows: the source was not asked."),
            }
        # the cap and the clock are different reasons and the sentence must name the real one
        cap_reached = len(rows) >= MAX_ADME_ROWS
        queried = ", ".join(asked)
        not_queried = tuple(at for at in ASSAY_TYPES if at not in asked)
        return {
            "retrieved": True,
            "assay_types_queried": {at: ASSAY_TYPES[at] for at in asked},
            "assay_types_not_queried": {at: ASSAY_TYPES[at] for at in not_queried},
            "assay_types_skipped_at_the_row_cap": unqueried or None,
            "assay_types_skipped_at_the_time_budget": out_of_time or None,
            "assay_types_cut_off_at_the_time_budget": stopped_by_clock or None,
            # a value of None means "never asked", not "ChEMBL has none"
            "rows_total_by_assay_type": totals,
            "rows_total": total,
            "rows_total_is_lower_bound": bool(unqueried) or bool(out_of_time),
            "rows_read": len(rows),
            "truncated": truncated,
            "rows_not_read": total - len(rows),
            "row_cap": MAX_ADME_ROWS,
            "adme_rows": kept,
            "rows_excluded_by_allowlist_by_standard_type":
                dict(sorted(dropped.items(), key=lambda kv: -kv[1])),
            "rows_excluded_by_allowlist": excluded,
            "rows_excluded_shown": (f"{len(excluded)} of {sum(dropped.values())}"
                                    if dropped else "0 of 0"),
            "rows_excluded_note": (
                "These rows were read from the same query and are NOT on this project's ADME "
                "allowlist, which is narrower than ADME on purpose. That is a statement about "
                "the allowlist, NOT about the data: some of these rows do hold real ADME "
                "measurements. They are shown with their own assay description so they can be "
                "read one by one, and each is unclassified: do not add them to the allowlisted "
                "rows, do not count them as ADME measurements, and do not call them 'not ADME'."
            ) if dropped else None,
            "note": (f"{len(kept)} of {len(rows)} rows read from assay types {queried} are on "
                     "this project's ADME allowlist. The rest are listed by standard_type: they "
                     "are NOT 'not ADME', they are not on the allowlist, which is narrower than "
                     "ADME on purpose. The generic types 'Activity' and 'Ratio' in particular "
                     "carry real ADME in some rows (efflux ratios, CYP induction) mixed with "
                     "cytotoxicity in others, and are excluded because the two cannot be told "
                     f"apart without reading each assay description. Assay types "
                     f"{', '.join(not_queried)} were not queried, so a measurement filed only "
                     "there is not in this result."
                     + (f" TRUNCATED: ChEMBL has at least {total} rows in the types queried and "
                        f"only {len(rows)} were read, so {total - len(rows)} were never examined "
                        "and this is a LOWER BOUND. Reason: "
                        + (f"the row cap of {MAX_ADME_ROWS} was reached." if cap_reached
                           else "the time budget ran out, NOT the row cap (which is "
                                f"{MAX_ADME_ROWS} and was not reached): ChEMBL was slow, not "
                                "empty." if out_of_time or stopped_by_clock
                           else "a page came back short of the total ChEMBL reported.")
                        if truncated else "")
                     + (" The row cap was reached before assay type(s) "
                        f"{', '.join(unqueried)} were queried at all, so how many rows they hold "
                        "is UNKNOWN, not zero." if unqueried else "")
                     + (" The time budget ran out before assay type(s) "
                        f"{', '.join(out_of_time)} were queried at all, so how many rows they "
                        "hold is UNKNOWN, not zero: ChEMBL was slow, not empty."
                        if out_of_time else "")),
        }
    finally:
        if own:
            await client.aclose()
