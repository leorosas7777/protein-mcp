"""Europe PMC: literature search for a UniProt protein."""
import re

import httpx

from common import clean_text, log, parse_accession
from uniprot import fetch_entry

EUROPEPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
_SORT_MAP = {"relevance": None, "recent": "P_PDATE_D desc", "cited": "CITED desc"}


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
    authors = clean_text(r.get("authorString"))
    if len(authors) > 150:
        authors = authors[:147] + "..."
    abstract = clean_text(r.get("abstractText"))
    # In "core" results the journal and publication types are nested objects.
    journal_obj = (r.get("journalInfo") or {}).get("journal") or {}
    journal = clean_text(journal_obj.get("title") or r.get("journalTitle")) or None
    raw_types = (r.get("pubTypeList") or {}).get("pubType")
    pub_types = raw_types if isinstance(raw_types, list) else ([raw_types] if raw_types else [])
    return {
        "title": clean_text(r.get("title")),
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


async def search(
    uniprot_id: str, topic: str, max_results: int, sort: str, include_preprints: bool
) -> dict:
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    if sort not in _SORT_MAP:
        return {"ok": False, "error": f"sort must be one of {sorted(_SORT_MAP)}"}
    max_results = max(1, min(max_results, 25))

    async with httpx.AsyncClient(timeout=30.0, follow_redirects=True) as client:
        # 1) Protein identity (gene symbol + recommended name) from UniProt.
        got = await fetch_entry(client, accession)
        if not got["ok"]:
            return got
        entry = got["entry"]

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
