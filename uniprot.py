"""UniProt: protein/gene name -> accession candidates, and entry lookup."""
import httpx

from common import log

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}"
FIELDS = "accession,id,protein_name,gene_names,organism_name,length,reviewed"


def parse_entry(entry: dict) -> dict:
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


async def resolve(name: str, organism_id: int, max_results: int) -> dict:
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
            candidates = [parse_entry(e) for e in resp.json().get("results", [])]
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


async def fetch_entry(client: httpx.AsyncClient, accession: str) -> dict:
    """Return {"ok": True, "entry": parsed} or {"ok": False, "error": ...}."""
    try:
        resp = await client.get(
            UNIPROT_ENTRY.format(accession=accession), params={"format": "json"}
        )
        if resp.status_code == 404:
            return {"ok": False, "error": f"UniProt accession {accession} not found."}
        resp.raise_for_status()
        return {"ok": True, "entry": parse_entry(resp.json())}
    except (httpx.HTTPError, ValueError) as exc:
        log.error("UniProt entry request failed: %s", exc)
        return {"ok": False, "error": f"UniProt request failed: {exc}"}
