"""RCSB PDB: experimental structures for a UniProt protein.

(Named rcsb.py, not pdb.py, to avoid shadowing Python's stdlib `pdb` module.)
"""
import re

import httpx

from common import clean_text, log, short

# >>> EDIT THIS LINE to change how many experimental PDB structures are returned <<<
MAX_PDB_MODELS = 7

CANDIDATE_POOL = 100  # best-resolution PDB entries examined before choosing

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"

_MUTATION_RE = re.compile(r"[A-Z]\d{1,4}[A-Z]")                 # e.g. G12C
_TITLE_MUTATION_RE = re.compile(r"(?<![A-Za-z0-9])([A-Z]\d{1,4}[A-Z])(?![A-Za-z0-9])")
_NOT_MUTATIONS = {"H2O"}

# Components that appear in many crystals but are not the ligand of interest.
_NUCLEOTIDES = {"GDP", "GTP", "GNP", "GCP", "GSP", "ADP", "ATP", "ANP", "ACP"}
_ADDITIVES = {
    "MG", "ZN", "CA", "NA", "K", "CL", "MN", "CO", "NI", "FE", "CD", "IOD", "BR",
    "SO4", "PO4", "NO3", "GOL", "EDO", "PEG", "PGE", "PG4", "DMS", "ACT", "TRS",
    "MES", "BME", "FMT", "EPE", "IMD", "CIT", "TLA", "MPD",
}
# Cofactors / prosthetic groups: part of the protein's own machinery, not a
# drug-like ligand. Listed separately so they do not pollute the ligand table.
_COFACTORS = {
    "HEM", "HEC", "HEA", "HEB", "FAD", "FMN", "NAD", "NAP", "NDP", "NAI",
    "SAM", "SAH", "COA", "PLP", "TPP", "BTN", "CLA", "BCL", "SF4", "FES",
}

_ENTRIES_QUERY = """
query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info { resolution_combined }
    rcsb_accession_info { initial_release_date }
    rcsb_primary_citation { pdbx_database_id_PubMed pdbx_database_id_DOI }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description pdbx_mutation }
      rcsb_polymer_entity_container_identifiers { uniprot_ids }
    }
    nonpolymer_entities {
      nonpolymer_comp { chem_comp { id name } }
    }
  }
}
"""


def _pdb_search_payload(accession: str, rows: int, mutation: str = "") -> dict:
    nodes = [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_name",
            "operator": "exact_match", "value": "UniProt"}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": accession}},
    ]
    if mutation:
        # Server-side text filter, so the 100-entry pool is built from entries that
        # already match (otherwise a rare variant never enters the best-resolution pool).
        nodes.append({"type": "terminal", "service": "full_text",
                      "parameters": {"value": mutation}})
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": nodes,
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            "results_content_type": ["experimental"],
        },
    }


def _res_sort_key(resolution):
    """Best resolution first; entries without resolution (e.g. NMR) go last."""
    return (resolution is None, resolution if resolution is not None else 0.0)


def _parse_pdb_entry(e: dict, accession: str) -> dict:
    """Flatten one RCSB entry. Links are built ONLY from returned IDs."""
    pdb_id = e.get("rcsb_id")
    res = (e.get("rcsb_entry_info") or {}).get("resolution_combined")
    resolution = res[0] if isinstance(res, list) and res else None
    if isinstance(resolution, (int, float)):
        resolution = round(resolution, 3)  # RCSB sometimes returns 1.82202726164
    methods = [m.get("method") for m in (e.get("exptl") or []) if m and m.get("method")]
    title = clean_text((e.get("struct") or {}).get("title"))

    # `pdbx_mutation` is free text typed by the depositor: it can hold "YES", be
    # empty although the title names the mutation, or list real mutations. Keep
    # the three cases separate and never treat "empty" as proof of wild type.
    mutations, unparsed, partners = [], [], []
    for ent in e.get("polymer_entities") or []:
        info = (ent or {}).get("rcsb_polymer_entity") or {}
        ids = ((ent or {}).get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []
        if accession in ids:
            raw = (info.get("pdbx_mutation") or "").strip()
            for token in (t.strip() for t in raw.split(",") if t.strip()):
                (mutations if _MUTATION_RE.fullmatch(token) else unparsed).append(token)
        elif info.get("pdbx_description"):
            partners.append(short(info["pdbx_description"], 60))
    title_muts = {m for m in _TITLE_MUTATION_RE.findall(title) if m not in _NOT_MUTATIONS}
    title_muts -= set(mutations)

    comps = {}
    for ent in e.get("nonpolymer_entities") or []:
        chem = (((ent or {}).get("nonpolymer_comp") or {}).get("chem_comp")) or {}
        if chem.get("id"):
            comps[chem["id"]] = chem.get("name")
    ligands = [
        {"id": i, "name": short(n, 60)}
        for i, n in sorted(comps.items())
        if i not in _NUCLEOTIDES and i not in _ADDITIVES and i not in _COFACTORS
    ]

    citation = e.get("rcsb_primary_citation") or {}
    pubmed, doi = citation.get("pdbx_database_id_PubMed"), citation.get("pdbx_database_id_DOI")
    links = {}
    if pdb_id:
        links["rcsb"] = f"https://www.rcsb.org/structure/{pdb_id}"
        links["download_cif"] = f"https://files.rcsb.org/download/{pdb_id}.cif"
    if pubmed:
        links["pubmed"] = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed}/"
    if doi:
        links["doi"] = f"https://doi.org/{doi}"
    return {
        "pdb_id": pdb_id,
        "title": title if len(title) <= 140 else title[:137] + "...",
        "method": "; ".join(methods) or None,
        "resolution_angstrom": resolution,
        "release_date": ((e.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10] or None,
        "mutations": sorted(set(mutations)),
        "mutations_from_title": sorted(title_muts),  # mentioned in the title, not in the mutation field
        "mutations_unparsed": sorted(set(unparsed)),  # depositor text that is not a mutation (e.g. "YES")
        "nucleotides": sorted(i for i in comps if i in _NUCLEOTIDES),
        "ligands": ligands,
        "cofactors": sorted(i for i in comps if i in _COFACTORS),
        "partner_chains": partners[:3],
        "links": links,
    }


def _coarse_key(s: dict):
    muts = tuple(sorted(set(s["mutations"]) | set(s["mutations_from_title"]) | set(s["mutations_unparsed"])))
    return (muts, tuple(s["nucleotides"]), bool(s["ligands"]))


def _select_diverse(entries: list, limit: int) -> list:
    """Pick up to `limit` entries. `entries` must be sorted best resolution first.

    Pass 1: the best-resolution entry for each distinct combination of
    (mutations, nucleotide, has-ligand), so the picks cover different states
    instead of ten near-identical structures. Pass 2: fill remaining slots with
    the next best resolution.
    """
    chosen, seen_keys, chosen_ids = [], set(), set()
    for s in entries:
        key = _coarse_key(s)
        if key not in seen_keys:
            seen_keys.add(key)
            chosen.append(s)
            chosen_ids.add(s["pdb_id"])
            if len(chosen) == limit:
                return chosen
    for s in entries:
        if s["pdb_id"] not in chosen_ids:
            chosen.append(s)
            chosen_ids.add(s["pdb_id"])
            if len(chosen) == limit:
                break
    return chosen


def _mutation_matches(s: dict, mutation: str) -> bool:
    """True if the requested variant is NAMED in the mutation field or the title.

    RCSB's full-text search also hits entries that only mention the text somewhere
    (e.g. a G12N structure whose record cites G12D), so each hit is confirmed here.
    A bare position such as "G12" matches any substitution at that position.
    """
    wanted = mutation.strip().lower()
    if not wanted:
        return True
    pos_only = re.fullmatch(r"[a-z]\d{1,4}", wanted) is not None
    tokens = {t.lower() for t in s["mutations"] + s["mutations_from_title"] + s["mutations_unparsed"]}
    for t in tokens:
        if t == wanted or (pos_only and re.fullmatch(re.escape(wanted) + r"[a-z]", t)):
            return True
    suffix = r"[A-Za-z]?" if pos_only else ""
    pattern = r"(?<![A-Za-z0-9])" + re.escape(mutation.strip()) + suffix + r"(?![A-Za-z0-9])"
    return re.search(pattern, s["title"], re.I) is not None


async def fetch_pdb(client, accession: str, limit: int, mutation: str = "") -> dict:
    # 1) Which experimental entries contain this UniProt protein (best resolution first).
    empty_note = (
        f"No experimental PDB structures found for this accession matching '{mutation}'."
        if mutation else "No experimental PDB structures found for this accession."
    )
    try:
        resp = await client.post(
            RCSB_SEARCH, json=_pdb_search_payload(accession, CANDIDATE_POOL, mutation)
        )
        if resp.status_code == 204:  # RCSB answers 204 No Content when nothing matches
            return {"ok": True, "total_entries": 0, "returned": 0, "structures": [],
                    "mutation_filter": mutation or None, "note": empty_note}
        resp.raise_for_status()
        sdata = resp.json()
    except (httpx.HTTPError, ValueError) as exc:
        log.error("RCSB search failed: %s", exc)
        return {"ok": False, "error": f"RCSB search failed: {exc}"}
    if not isinstance(sdata, dict) or "result_set" not in sdata:
        return {"ok": False, "error": "Unexpected RCSB search response format.",
                "response_keys": sorted(sdata.keys()) if isinstance(sdata, dict) else None}
    ids = [r.get("identifier") for r in sdata["result_set"] if isinstance(r, dict) and r.get("identifier")]
    total = sdata.get("total_count")
    if not ids:
        return {"ok": True, "total_entries": total or 0, "returned": 0, "structures": [],
                "mutation_filter": mutation or None, "note": empty_note}

    # 2) Details for the candidates in ONE GraphQL request.
    try:
        resp = await client.post(RCSB_GRAPHQL, json={"query": _ENTRIES_QUERY, "variables": {"ids": ids}})
        resp.raise_for_status()
        gdata = resp.json()
    except (httpx.HTTPError, ValueError) as exc:
        log.error("RCSB GraphQL failed: %s", exc)
        return {"ok": False, "error": f"RCSB details request failed: {exc}"}
    raw = ((gdata or {}).get("data") or {}).get("entries") if isinstance(gdata, dict) else None
    if raw is None:
        return {"ok": False, "error": "Unexpected RCSB GraphQL response format.",
                "graphql_errors": str((gdata or {}).get("errors"))[:300] if isinstance(gdata, dict) else None}

    parsed = [_parse_pdb_entry(e, accession) for e in raw if isinstance(e, dict)]
    parsed.sort(key=lambda s: _res_sort_key(s["resolution_angstrom"]))
    examined = len(parsed)

    # Confirm the variant: RCSB full-text hits can merely mention it (false positives).
    discarded = []
    if mutation:
        confirmed = [s for s in parsed if _mutation_matches(s, mutation)]
        discarded = [{"pdb_id": s["pdb_id"], "title": s["title"]} for s in parsed if s not in confirmed]
        parsed = confirmed
        if not parsed:
            return {
                "ok": True, "total_entries": total, "mutation_filter": mutation, "examined": examined,
                "confirmed_in_examined": 0, "discarded_count": len(discarded),
                "discarded_not_confirmed": discarded[:10], "returned": 0, "structures": [],
                "note": (f"RCSB text search found {total} entries, but none of the {examined} examined "
                         f"names '{mutation}' in its mutation field or title."),
            }

    chosen = _select_diverse(parsed, limit)
    chosen.sort(key=lambda s: _res_sort_key(s["resolution_angstrom"]))
    result = {
        "ok": True,
        "total_entries": total,
        "mutation_filter": mutation or None,
        "examined": examined,
        "returned": len(chosen),
        "selection": (
            (f"RCSB text search for '{mutation}' found {total} entries; examined the {examined} "
             f"best-resolution ones, {len(parsed)} of which name the variant in the mutation field or "
             f"title ({len(discarded)} discarded). Picked up to {limit} from the confirmed ones: "
             if mutation else
             f"Examined the {examined} best-resolution entries of {total}; picked up to {limit}: ")
            + "first the best-resolution entry for each distinct (mutations, nucleotide, ligand yes/no) "
            "combination, then the next best resolution."
        ),
        "structures": chosen,
    }
    if mutation:
        result["confirmed_in_examined"] = len(parsed)
        result["discarded_count"] = len(discarded)
        result["discarded_not_confirmed"] = discarded[:10]
    return result