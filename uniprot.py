"""UniProt: protein/gene name -> accession candidates, entry lookup and profile annotation."""
import re

import httpx

from common import clean_text, log

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}"
FIELDS = "accession,id,protein_name,gene_names,organism_name,length,reviewed"

# Evidence codes seen in the real KRAS entry (2026-10-20). Others pass through as the raw code.
ECO_LABELS = {
    "ECO:0000269": "experimental",
    "ECO:0000250": "by similarity",       # inferred from another protein (source: UniProtKB)
    "ECO:0000305": "curator inference",
}


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


async def fetch_entry_json(client: httpx.AsyncClient, accession: str) -> dict:
    """Return {"ok": True, "json": raw entry} or {"ok": False, "error": ...}.

    A deleted entry answers HTTP 200 with entryType "Inactive" (verified with
    Q0ZZZ9), not 404: it is reported as an error, never parsed as a protein.
    """
    try:
        resp = await client.get(
            UNIPROT_ENTRY.format(accession=accession), params={"format": "json"}
        )
        if resp.status_code == 404:
            return {"ok": False, "error": f"UniProt accession {accession} not found."}
        resp.raise_for_status()
        data = resp.json()
    except (httpx.HTTPError, ValueError) as exc:
        log.error("UniProt entry request failed: %s", exc)
        return {"ok": False, "error": f"UniProt request failed: {exc}"}
    if data.get("entryType") == "Inactive":
        reason = data.get("inactiveReason") or {}
        detail = "; ".join(x for x in (reason.get("inactiveReasonType"), reason.get("deletedReason")) if x)
        return {"ok": False, "error": f"UniProt accession {accession} is inactive ({detail or 'no reason given'})."}
    return {"ok": True, "json": data}


async def fetch_entry(client: httpx.AsyncClient, accession: str) -> dict:
    """Return {"ok": True, "entry": parsed} or {"ok": False, "error": ...}."""
    result = await fetch_entry_json(client, accession)
    if not result["ok"]:
        return result
    return {"ok": True, "entry": parse_entry(result["json"])}


def _evidence(evidences) -> dict:
    """PubMed IDs and evidence types behind one statement."""
    evidences = evidences or []
    codes = sorted({e.get("evidenceCode") for e in evidences if e.get("evidenceCode")})
    return {
        "pubmed": sorted({e["id"] for e in evidences if e.get("source") == "PubMed" and e.get("id")}),
        "evidence": [ECO_LABELS.get(c, c) for c in codes],
    }


def _sequence_feature(f: dict) -> dict:
    """One 'Natural variant' or 'Mutagenesis' feature. `change` is only built for
    single-residue substitutions; anything else (e.g. the G10 -> GG insertion)
    keeps change=None and the raw fields."""
    loc = f.get("location") or {}
    start, end = (loc.get("start") or {}).get("value"), (loc.get("end") or {}).get("value")
    alt = f.get("alternativeSequence") or {}
    orig, alts = alt.get("originalSequence"), alt.get("alternativeSequences") or []
    change = None
    if start is not None and start == end and orig and len(orig) == 1 and alts and all(len(a) == 1 for a in alts):
        change = f"{orig}{start}{'/'.join(alts)}"
    dbsnp = [x["id"] for x in f.get("featureCrossReferences") or [] if x.get("database") == "dbSNP" and x.get("id")]
    return {
        "change": change,
        "start": start,
        "end": end,
        "original": orig,
        "alternatives": alts,
        "description": f.get("description"),
        "feature_id": f.get("featureId"),
        "dbsnp": dbsnp,
        **_evidence(f.get("evidences")),
        "links": {"dbsnp": [f"https://www.ncbi.nlm.nih.gov/snp/{r}" for r in dbsnp]},
    }


def parse_profile(entry: dict) -> dict:
    """Profile fields of one UniProt entry, as returned (validated on P01116, 2026-10-20)."""
    comments = entry.get("comments") or []

    def of_type(t):
        return [c for c in comments if c.get("commentType") == t]

    function = [
        {"text": clean_text(t.get("value")), **_evidence(t.get("evidences"))}
        for c in of_type("FUNCTION") for t in c.get("texts") or [] if t.get("value")
    ]

    locations = []
    for c in of_type("SUBCELLULAR LOCATION"):
        for sl in c.get("subcellularLocations") or []:
            loc = sl.get("location") or {}
            locations.append({
                "molecule": c.get("molecule"),  # e.g. "Isoform 2B"; None = the entry as a whole
                "location": loc.get("value"),
                "topology": (sl.get("topology") or {}).get("value"),
                "orientation": (sl.get("orientation") or {}).get("value"),
                **_evidence(loc.get("evidences")),
            })

    # Some DISEASE comments have no `disease` block, only a note (2 of 10 in KRAS):
    # name stays None and the note is kept verbatim.
    diseases = []
    for c in of_type("DISEASE"):
        d = c.get("disease") or {}
        xref = d.get("diseaseCrossReference") or {}
        mim = xref.get("id") if xref.get("database") == "MIM" else None
        note_texts = (c.get("note") or {}).get("texts") or []
        links = {}
        if d.get("diseaseAccession"):
            links["uniprot_disease"] = f"https://www.uniprot.org/diseases/{d['diseaseAccession']}"
        if mim:
            links["omim"] = f"https://omim.org/entry/{mim}"
        diseases.append({
            "name": d.get("diseaseId"),
            "acronym": d.get("acronym"),
            "uniprot_disease_id": d.get("diseaseAccession"),
            "mim": mim,
            "description": d.get("description"),
            "note": " ".join(clean_text(t.get("value")) for t in note_texts if t.get("value")) or None,
            **_evidence((d.get("evidences") or []) + [e for t in note_texts for e in t.get("evidences") or []]),
            "links": links,
        })

    features = entry.get("features") or []
    xrefs = entry.get("uniProtKBCrossReferences") or []
    # UniProt gives versioned gene IDs (ENSG00000133703.15), one per transcript:
    # strip the version (Open Targets needs it bare) and deduplicate.
    ensembl_genes = sorted({
        re.sub(r"\.\d+$", "", p["value"])
        for x in xrefs if x.get("database") == "Ensembl"
        for p in x.get("properties") or [] if p.get("key") == "GeneId" and p.get("value")
    })
    return {
        **parse_entry(entry),
        "taxon_id": (entry.get("organism") or {}).get("taxonId"),
        "annotation_score": entry.get("annotationScore"),
        "protein_existence": entry.get("proteinExistence"),
        "last_annotation_update": (entry.get("entryAudit") or {}).get("lastAnnotationUpdateDate"),
        "function": function,
        "subcellular_locations": locations,
        "diseases": diseases,
        "natural_variants": [_sequence_feature(f) for f in features if f.get("type") == "Natural variant"],
        "mutagenesis": [_sequence_feature(f) for f in features if f.get("type") == "Mutagenesis"],
        "ensembl_gene_ids": ensembl_genes,
        "chembl_target_ids": [x["id"] for x in xrefs if x.get("database") == "ChEMBL" and x.get("id")],
    }
