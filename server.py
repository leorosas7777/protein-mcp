"""Protein research MCP server.

Tools:
  - resolve_protein: protein/gene name -> UniProt accession candidates.
  - search_literature: UniProt accession -> papers from Europe PMC
    (PubMed/MEDLINE records, PMC full-text records and, optionally, preprints).

IMPORTANT: this server talks to Claude Desktop over stdio.
Never use print() here: anything written to stdout corrupts the protocol.
Use logging (goes to stderr) instead.
"""
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


if __name__ == "__main__":
    mcp.run()  # stdio transport
