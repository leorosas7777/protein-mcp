"""Protein research MCP server.

Tools:
  - resolve_protein: protein/gene name -> UniProt accession candidates.
  - search_literature: UniProt accession -> papers from Europe PMC
    (PubMed/MEDLINE records, PMC full-text records and, optionally, preprints).
  - get_structures: UniProt accession -> a few experimental PDB structures
    (RCSB) plus the AlphaFold DB model with its confidence.

IMPORTANT: this server talks to Claude Desktop over stdio.
Never use print() here: anything written to stdout corrupts the protocol.
Use logging (goes to stderr) instead.
"""
import asyncio
import html
import logging
import re

import httpx
from mcp.server.fastmcp import FastMCP

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("protein-mcp")

mcp = FastMCP("protein-research")

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,id,protein_name,gene_names,organism_name,length,reviewed"


def _parse_entry(entry: dict) -> dict:
    """Flatten one UniProt JSON entry into the fields we care about."""
    desc = entry.get("proteinDescription", {})
    name = (
        desc.get("recommendedName", {}).get("fullName", {}).get("value")
        or desc.get("submissionNames", [{}])[0].get("fullName", {}).get("value")
        or "unknown"
    )
    genes = entry.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value") if genes else None
    accession = entry.get("primaryAccession")
    return {
        "uniprot_id": accession,
        "entry_name": entry.get("uniProtkbId"),
        "protein_name": name,
        "gene": gene,
        "organism": entry.get("organism", {}).get("scientificName"),
        "length_aa": entry.get("sequence", {}).get("length"),
        "reviewed": entry.get("entryType", "").startswith("UniProtKB reviewed"),
        "source_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
    }


@mcp.tool()
async def resolve_protein(
    name: str, organism_id: int = 9606, max_results: int = 5
) -> dict:
    """Resolve a protein or gene name to UniProt accession candidates.

    Args:
        name: Protein or gene name, e.g. "KRAS", "BACE1", "insulin receptor".
        organism_id: NCBI taxonomy ID. Default 9606 (human). Use 10090 for mouse.
        max_results: How many candidates to return (1-10).

    Returns a list of candidates. If more than one plausible match exists,
    ask the user which one they mean before continuing. Never guess.
    """
    max_results = max(1, min(max_results, 10))
    base = f"(organism_id:{organism_id}) AND (reviewed:true)"
    clean = name.strip()
    # 1) Exact gene-symbol match (precise, avoids entries that merely mention
    #    the name in their text). 2) Fall back to free text if nothing matches.
    queries = []
    if clean and " " not in clean:
        queries.append(f"(gene_exact:{clean}) AND {base}")
    queries.append(f"({clean}) AND {base}")

    candidates: list[dict] = []
    query = queries[-1]
    async with httpx.AsyncClient(timeout=20.0) as client:
        for query in queries:
            params = {"query": query, "fields": FIELDS, "format": "json", "size": max_results}
            try:
                resp = await client.get(UNIPROT_SEARCH, params=params)
                resp.raise_for_status()
            except httpx.HTTPError as exc:
                log.error("UniProt request failed: %s", exc)
                return {"ok": False, "error": f"UniProt request failed: {exc}"}
            candidates = [_parse_entry(e) for e in resp.json().get("results", [])]
            if candidates:
                break

    if not candidates:
        return {
            "ok": True,
            "candidates": [],
            "note": (
                "No reviewed (Swiss-Prot) entry found for this organism. "
                "Try another spelling, the gene symbol, or a different organism_id."
            ),
        }
    return {
        "ok": True,
        "query_used": query,
        "candidates": candidates,
        "needs_user_confirmation": len(candidates) > 1,
    }


# ---------------------------------------------------------------------------
# Step 2: literature search (Europe PMC)
# ---------------------------------------------------------------------------
EUROPEPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}"
_ACCESSION_RE = re.compile(r"[A-Z0-9]{6,10}")
_INLINE_TAG_RE = re.compile(r"</?(?:sup|sub|i|b|em|strong|u|span)\b[^>]*>", re.I)
_TAG_RE = re.compile(r"</?[a-zA-Z][^>]*>")
_SORT_MAP = {"relevance": None, "recent": "P_PDATE_D desc", "cited": "CITED desc"}


def _clean_text(value) -> str:
    """Clean HTML that Europe PMC leaves in titles and abstracts.

    Entities are decoded FIRST: Europe PMC often sends escaped tags such as
    "&lt;sup&gt;", which would turn into real tags if decoded after stripping.
    Inline tags (sup/sub/i/b...) are removed without a space (KRAS<sup>G13C</sup>
    -> KRASG13C); block tags (h4, p...) become a space.
    """
    if not value:
        return ""
    text = html.unescape(str(value))
    text = _INLINE_TAG_RE.sub("", text)
    text = _TAG_RE.sub(" ", text)
    return re.sub(r"\s+", " ", text).strip()


def _topic_clause(topic: str) -> str:
    """Turn the optional topic into a Europe PMC clause, or "" if empty.

    Plain words are restricted to title/abstract and ALL must appear, so papers
    that only mention the topic in their full text do not add noise. Advanced
    syntax (contains ":" or parentheses or AND/OR/NOT) is passed through as is.
    """
    topic = (topic or "").strip()
    if not topic:
        return ""
    if re.search(r"[:()]|\b(?:AND|OR|NOT)\b", topic):
        return f"({topic})"
    parts = []
    for word in re.findall(r"[A-Za-z0-9][A-Za-z0-9\-]*", topic):
        parts.append(f'TITLE_ABS:"{word}"' if "-" in word else f"TITLE_ABS:{word}")
    return "(" + " AND ".join(parts) + ")" if parts else ""


def _build_query(gene, protein_name, topic: str, include_preprints: bool) -> str:
    """Build a Europe PMC query. Returns "" if there is nothing to search for."""
    terms = []
    for term in (gene, protein_name):
        if term:
            terms.append('TITLE_ABS:"%s"' % str(term).replace('"', "'"))
    if not terms:
        return ""
    query = "(" + " OR ".join(terms) + ")"
    clause = _topic_clause(topic)
    if clause:
        query += f" AND {clause}"
    if not include_preprints:
        query += " NOT SRC:PPR"
    return query


def _parse_paper(r: dict) -> dict:
    """Flatten one Europe PMC result. Links are built ONLY from returned IDs."""
    source, rec_id = r.get("source"), r.get("id")
    pmid, pmcid, doi = r.get("pmid"), r.get("pmcid"), r.get("doi")
    links = {}
    if pmid:
        links["pubmed"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    if doi:
        links["doi"] = f"https://doi.org/{doi}"
    if pmcid:
        links["pmc"] = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
    if source and rec_id:
        links["europepmc"] = f"https://europepmc.org/article/{source}/{rec_id}"
    authors = _clean_text(r.get("authorString"))
    if len(authors) > 150:
        authors = authors[:147] + "..."
    abstract = _clean_text(r.get("abstractText"))
    # In "core" results the journal and publication types are nested objects.
    journal_obj = (r.get("journalInfo") or {}).get("journal") or {}
    journal = _clean_text(journal_obj.get("title") or r.get("journalTitle")) or None
    raw_types = (r.get("pubTypeList") or {}).get("pubType")
    pub_types = raw_types if isinstance(raw_types, list) else ([raw_types] if raw_types else [])
    return {
        "title": _clean_text(r.get("title")),
        "authors": authors,
        "journal": journal,
        "year": r.get("pubYear"),
        "pub_types": pub_types,  # e.g. "Journal Article", "Review"
        "source_db": source,  # MED = PubMed/MEDLINE, PPR = preprint, PMC = PMC-only
        "is_preprint": source == "PPR",
        "pmid": pmid,
        "pmcid": pmcid,
        "doi": doi,
        "open_access": r.get("isOpenAccess") == "Y",
        "full_text_in_europepmc": r.get("inEPMC") == "Y",
        "has_pdf": r.get("hasPDF") == "Y",
        "cited_by_count": r.get("citedByCount"),
        "abstract": abstract or None,
        "links": links,
    }


@mcp.tool()
async def search_literature(
    uniprot_id: str,
    topic: str = "",
    max_results: int = 10,
    sort: str = "relevance",
    include_preprints: bool = False,
) -> dict:
    """Search the scientific literature (Europe PMC) for a protein.

    Args:
        uniprot_id: UniProt accession from resolve_protein, e.g. "P01116".
        topic: Optional extra words to narrow the search, e.g. "inhibitor",
            "crystal structure", "allosteric". Plain words must ALL appear in the
            title or abstract. Advanced Europe PMC syntax (with ":", "(" or
            AND/OR/NOT) is passed through unchanged.
        max_results: Papers to return (1-25). Keep it small: abstracts are long.
        sort: "relevance" (default), "recent" or "cited".
        include_preprints: Include bioRxiv/medRxiv etc. Default False. Preprints
            are NOT peer reviewed; they are flagged with is_preprint=true.

    Rules for using the result:
    - Base every statement ONLY on the returned titles and abstracts, and cite
      each one with the links provided. Never add facts that are not in them.
    - Only abstracts are returned, not full text. open_access/has_pdf only say
      whether a legal free full text exists.
    - If total_hits is much larger than returned, say the review is partial.
    """
    accession = uniprot_id.strip().upper()
    if not _ACCESSION_RE.fullmatch(accession):
        return {
            "ok": False,
            "error": (
                f"'{uniprot_id}' does not look like a UniProt accession "
                "(e.g. P01116). Use resolve_protein first."
            ),
        }
    if sort not in _SORT_MAP:
        return {"ok": False, "error": f"sort must be one of {sorted(_SORT_MAP)}"}
    max_results = max(1, min(max_results, 25))

    async with httpx.AsyncClient(timeout=30.0, follow_redirects=True) as client:
        # 1) Protein identity (gene symbol + recommended name) from UniProt.
        try:
            resp = await client.get(
                UNIPROT_ENTRY.format(accession=accession), params={"format": "json"}
            )
            if resp.status_code == 404:
                return {"ok": False, "error": f"UniProt accession {accession} not found."}
            resp.raise_for_status()
            entry = _parse_entry(resp.json())
        except (httpx.HTTPError, ValueError) as exc:
            log.error("UniProt entry request failed: %s", exc)
            return {"ok": False, "error": f"UniProt request failed: {exc}"}

        gene = entry["gene"]
        protein_name = None if entry["protein_name"] == "unknown" else entry["protein_name"]
        query = _build_query(gene, protein_name, topic, include_preprints)
        if not query:
            return {
                "ok": False,
                "error": f"No gene symbol or protein name found for {accession}; cannot search.",
            }

        # 2) Literature search in Europe PMC.
        params = {
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": max_results,
        }
        if _SORT_MAP[sort]:
            params["sort"] = _SORT_MAP[sort]
        try:
            resp = await client.get(EUROPEPMC_SEARCH, params=params)
            resp.raise_for_status()
            data = resp.json()
        except (httpx.HTTPError, ValueError) as exc:
            log.error("Europe PMC request failed: %s", exc)
            return {
                "ok": False,
                "error": f"Europe PMC request failed: {exc}",
                "query_used": query,
            }

    results = (data.get("resultList") or {}).get("result") if isinstance(data, dict) else None
    if results is None:
        return {
            "ok": False,
            "error": "Unexpected Europe PMC response format (no resultList.result).",
            "query_used": query,
            "response_keys": sorted(data.keys()) if isinstance(data, dict) else None,
        }

    papers = [_parse_paper(p) for p in results]
    return {
        "ok": True,
        "protein": {
            "uniprot_id": accession,
            "gene": gene,
            "protein_name": protein_name,
            "source_url": entry["source_url"],
        },
        "query_used": query,
        "total_hits": data.get("hitCount"),
        "returned": len(papers),
        "sort": sort,
        "include_preprints": include_preprints,
        "papers": papers,
        "notes": [
            "Search matches the gene symbol or protein name in title/abstract; "
            "ambiguous symbols can bring unrelated papers, so check titles.",
            "Only abstracts are returned, not full text.",
            "All links are built from IDs returned by Europe PMC.",
        ],
    }


# ---------------------------------------------------------------------------
# Step 3: 3D structures (RCSB PDB experimental + AlphaFold DB prediction)
# ---------------------------------------------------------------------------
# >>> EDIT THIS LINE to change how many experimental PDB structures are returned <<<
MAX_PDB_MODELS = 7

# >>> EDIT THIS LINE: AlphaFold returns the canonical model and sometimes an isoform model.
# The canonical one is preferred unless the isoform's mean pLDDT is higher by at least this many points. <<<
AF_ISOFORM_MARGIN = 5.0

CANDIDATE_POOL = 100  # best-resolution PDB entries examined before choosing

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

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


def _short(value, limit: int = 80):
    text = _clean_text(value)
    if len(text) > limit:
        return text[: limit - 3] + "..."
    return text or None


def _plddt_band(score) -> str:
    """Standard AlphaFold confidence bands for a pLDDT score."""
    if not isinstance(score, (int, float)):
        return "unknown"
    if score >= 90:
        return "very high"
    if score >= 70:
        return "confident"
    if score >= 50:
        return "low"
    return "very low"


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
    title = _clean_text((e.get("struct") or {}).get("title"))

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
            partners.append(_short(info["pdbx_description"], 60))
    title_muts = {m for m in _TITLE_MUTATION_RE.findall(title) if m not in _NOT_MUTATIONS}
    title_muts -= set(mutations)

    comps = {}
    for ent in e.get("nonpolymer_entities") or []:
        chem = (((ent or {}).get("nonpolymer_comp") or {}).get("chem_comp")) or {}
        if chem.get("id"):
            comps[chem["id"]] = chem.get("name")
    ligands = [
        {"id": i, "name": _short(n, 60)}
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


def _state_key(s: dict):
    """Grouping used when the user already fixed the variant: nucleotide + ligand yes/no."""
    return (tuple(s["nucleotides"]), bool(s["ligands"]))


def _select_diverse(entries: list, limit: int, group_by_mutation: bool = True,
                    prefer_ligand_fill: bool = False) -> list:
    """Pick up to `limit` entries. `entries` must be sorted best resolution first.

    Pass 1: the best-resolution entry for each distinct group, so the picks cover
    different states instead of ten near-identical structures. The group is
    (mutations, nucleotide, has-ligand), or just (nucleotide, has-ligand) when
    `group_by_mutation` is False (a mutation filter is active: second-site
    mutations should not eat the slots).
    Pass 2: fill the remaining slots with the next best resolution; with
    `prefer_ligand_fill`, ligand-bound entries go first (still by resolution).
    """
    key_fn = _coarse_key if group_by_mutation else _state_key
    chosen, seen_keys, chosen_ids = [], set(), set()
    for s in entries:
        key = key_fn(s)
        if key not in seen_keys:
            seen_keys.add(key)
            chosen.append(s)
            chosen_ids.add(s["pdb_id"])
            if len(chosen) == limit:
                return chosen
    rest = [s for s in entries if s["pdb_id"] not in chosen_ids]
    if prefer_ligand_fill:
        rest.sort(key=lambda s: not s["ligands"])  # stable: keeps resolution order inside each half
    for s in rest:
        chosen.append(s)
        if len(chosen) == limit:
            break
    return chosen


async def _fetch_pdb(client, accession: str, limit: int, mutation: str = "") -> dict:
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

    raw_entries = [e for e in raw if isinstance(e, dict)]
    parsed = [_parse_pdb_entry(e, accession) for e in raw_entries]
    examined = len(parsed)

    dropped = 0
    if mutation:
        # RCSB's full-text search also matches entries that merely mention the text
        # somewhere (e.g. a G12N structure for a G12D query). Keep only entries that
        # name the mutation in the mutation field or in the (full) title.
        needle = mutation.lower()
        kept = []
        for e, s in zip(raw_entries, parsed):
            full_title = _clean_text((e.get("struct") or {}).get("title")).lower()
            named = {m.lower() for m in s["mutations"] + s["mutations_from_title"]}
            if needle in named or needle in full_title:
                kept.append(s)
        dropped = examined - len(kept)
        parsed = kept
        if not parsed:
            return {"ok": True, "total_entries": total, "mutation_filter": mutation, "examined": examined,
                    "returned": 0, "structures": [],
                    "note": (f"RCSB text-matched {total} entries, but none names '{mutation}' in its "
                             "mutation field or title.")}

    parsed.sort(key=lambda s: _res_sort_key(s["resolution_angstrom"]))
    chosen = _select_diverse(parsed, limit, group_by_mutation=not mutation, prefer_ligand_fill=bool(mutation))
    chosen.sort(key=lambda s: _res_sort_key(s["resolution_angstrom"]))
    if mutation:
        selection = (
            f"RCSB text-matched {total} entries for '{mutation}'; examined the {examined} best-resolution ones "
            f"and kept the {len(parsed)} that name it in the mutation field or title ({dropped} discarded as "
            f"false matches). Picked up to {limit}: first the best-resolution entry for each "
            "(nucleotide, ligand yes/no) combination, then ligand-bound entries by resolution, then the rest."
        )
    else:
        selection = (
            f"Examined the {examined} best-resolution entries of {total}; picked up to {limit}: "
            "first the best-resolution entry for each distinct (mutations, nucleotide, ligand yes/no) "
            "combination, then the next best resolution."
        )
    return {
        "ok": True,
        "total_entries": total,
        "mutation_filter": mutation or None,
        "examined": examined,
        "confirmed_matching_filter": len(parsed) if mutation else None,
        "returned": len(chosen),
        "selection": selection,
        "structures": chosen,
    }


async def _fetch_alphafold(client, accession: str) -> dict:
    try:
        resp = await client.get(ALPHAFOLD_API.format(accession=accession))
        if resp.status_code == 404:
            return {"ok": True, "available": False, "note": "No AlphaFold DB model for this accession."}
        resp.raise_for_status()
        data = resp.json()
    except (httpx.HTTPError, ValueError) as exc:
        log.error("AlphaFold request failed: %s", exc)
        return {"ok": False, "error": f"AlphaFold DB request failed: {exc}"}
    models = [m for m in (data if isinstance(data, list) else [data]) if isinstance(m, dict)]
    if not models:
        return {"ok": True, "available": False, "note": "AlphaFold DB returned no model."}

    def score(m):
        value = m.get("globalMetricValue")
        return value if isinstance(value, (int, float)) else 0

    def isoform_of(m):
        """'AF-P01116-2-F1' -> 'P01116-2'; canonical entries ('AF-P01116-F1') -> None."""
        mo = re.fullmatch(r"AF-([A-Z0-9]+)-(\d+)-F\d+", str(m.get("entryId") or ""))
        return f"{mo.group(1)}-{mo.group(2)}" if mo else None

    # Prefer the canonical-sequence model unless an isoform model is clearly more
    # confident (AF_ISOFORM_MARGIN pLDDT points): a ~1 point gap is not meaningful,
    # while mixing up sequences hurts traceability.
    top = max(models, key=score)
    canonical = [m for m in models if isoform_of(m) is None]
    best, reason = top, "Highest mean pLDDT."
    if canonical and isoform_of(top) is not None:
        top_canonical = max(canonical, key=score)
        gap = round(score(top) - score(top_canonical), 2)
        if gap < AF_ISOFORM_MARGIN:
            best = top_canonical
            reason = (f"Canonical-sequence model preferred: the best isoform model ({top.get('entryId')}) "
                      f"is only {gap} pLDDT points higher (margin required: {AF_ISOFORM_MARGIN}).")
        else:
            reason = (f"Isoform model chosen: it is {gap} pLDDT points above the canonical one "
                      f"(margin required: {AF_ISOFORM_MARGIN}).")
    plddt = best.get("globalMetricValue")
    start, end = best.get("uniprotStart"), best.get("uniprotEnd")
    best_isoform = isoform_of(best)
    return {
        "ok": True,
        "available": True,
        "entry_id": best.get("entryId"),
        "selection_reason": reason,
        "isoform": best_isoform,  # None = the canonical UniProt sequence
        "is_canonical_sequence": best_isoform is None,
        "other_models": [
            {"entry_id": m.get("entryId"), "isoform": isoform_of(m), "mean_plddt": m.get("globalMetricValue")}
            for m in models if m is not best
        ],
        "mean_plddt": plddt,
        "confidence_band": _plddt_band(plddt),
        "fraction_of_residues": {
            "very_high_over_90": best.get("fractionPlddtVeryHigh"),
            "confident_70_90": best.get("fractionPlddtConfident"),
            "low_50_70": best.get("fractionPlddtLow"),
            "very_low_under_50": best.get("fractionPlddtVeryLow"),
        },
        "residues": f"{start}-{end}" if start and end else None,
        "model_version": best.get("latestVersion"),
        "created": best.get("modelCreatedDate"),
        "models_considered": len(models),
        "links": {
            "alphafold": f"https://alphafold.ebi.ac.uk/entry/{accession}",
            "pdb": best.get("pdbUrl"),
            "cif": best.get("cifUrl"),
            "pae_image": best.get("paeImageUrl"),
        },
    }


async def _safe(coro, label: str) -> dict:
    """The API shapes are partly unverified: never let one source crash the tool."""
    try:
        return await coro
    except Exception as exc:  # noqa: BLE001
        log.exception("%s failed", label)
        return {"ok": False, "error": f"{label} failed unexpectedly: {exc!r}"}


@mcp.tool()
async def get_structures(uniprot_id: str, max_models: int = 0, mutation: str = "") -> dict:
    """Get 3D structures for a protein: experimental PDB entries + AlphaFold model.

    Args:
        uniprot_id: UniProt accession from resolve_protein, e.g. "P01116".
        max_models: How many experimental PDB structures to return (1-25).
            0 (default) uses the MAX_PDB_MODELS constant in server.py.
        mutation: Optional text to restrict the PDB entries to one variant,
            e.g. "G12D". Use it when the user cares about a specific mutation;
            without it the picks come from the best-resolution entries overall,
            which may not include rare variants. RCSB matches it as free text and
            the tool then keeps only entries that name it in the mutation field or
            title (see `selection` for how many were discarded).

    The PDB structures are a SELECTION, not the full list: state "N of
    total_entries" and read `selection` to explain how they were chosen.
    The AlphaFold model is the canonical-sequence one unless an isoform model is
    clearly more confident; `selection_reason` says which case applied.

    Rules for using the result:
    - Report only what is in the result; never invent ligands, mutations or
      resolutions. Link every structure with its `links`.
    - Mutation data is depositor-written and unreliable. `mutations` is the
      official field, `mutations_from_title` are mutations named only in the
      title, `mutations_unparsed` is non-mutation text. An empty `mutations`
      does NOT prove wild type: check the title and say it is unconfirmed.
    - Resolution is not binding-site quality: check the ligand density before
      trusting a pose.
    - The AlphaFold model may be a non-canonical isoform (see `isoform`,
      `is_canonical_sequence`, `other_models`): always say which one was chosen.
    - AlphaFold gives one static prediction with no ligands, nucleotides or ions;
      do not use it to describe inhibitor-bound conformations.
    """
    accession = uniprot_id.strip().upper()
    if not _ACCESSION_RE.fullmatch(accession):
        return {
            "ok": False,
            "error": (
                f"'{uniprot_id}' does not look like a UniProt accession "
                "(e.g. P01116). Use resolve_protein first."
            ),
        }
    limit = max_models if max_models > 0 else MAX_PDB_MODELS
    limit = max(1, min(limit, 25))
    mutation = (mutation or "").strip()
    if mutation and not re.fullmatch(r"[A-Za-z0-9 _\-/.]{1,40}", mutation):
        return {
            "ok": False,
            "error": f"Invalid mutation filter '{mutation}': use letters, digits and - _ / . only (max 40 chars).",
        }

    async with httpx.AsyncClient(timeout=30.0, follow_redirects=True) as client:
        pdb, af = await asyncio.gather(
            _safe(_fetch_pdb(client, accession, limit, mutation), "RCSB PDB"),
            _safe(_fetch_alphafold(client, accession), "AlphaFold DB"),
        )

    return {
        "ok": bool(pdb.get("ok") or af.get("ok")),
        "uniprot_id": accession,
        "experimental": pdb,
        "alphafold": af,
        "notes": [
            "PDB structures are a diverse selection, not the complete list.",
            "AlphaFold is a single static prediction without ligands, nucleotides or ions.",
            "All links are built from IDs/URLs returned by RCSB and AlphaFold DB.",
        ],
    }


if __name__ == "__main__":
    mcp.run()  # stdio transport