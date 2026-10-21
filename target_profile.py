"""get_protein_profile: UniProt annotation + Open Targets target profile.

(Named target_profile.py, not profile.py, to avoid shadowing Python's stdlib `profile` module.)

The two sources stay separate in the output. get_protein_profile returns the FULL
result (server.py saves it to dossiers/<accession>/profile.json) and compact()
builds the small summary that goes to the chat.
"""
import httpx

import opentargets
import uniprot
from common import parse_accession, safe, short

TOP_DISEASES = 10      # associations listed per block (overall / oncology) in the summary
FUNCTION_CHARS = 600   # function text kept in the summary; the full text is in the dossier
HUMAN = 9606


async def _open_targets(client, prof: dict) -> dict:
    if prof["taxon_id"] != HUMAN:
        return {"ok": True, "found": False, "note": "Open Targets covers human targets only: sin registro."}
    genes = prof["ensembl_gene_ids"]
    if not genes:
        return {"ok": True, "found": False, "note": "UniProt lists no Ensembl gene: sin registro en Open Targets."}
    ot = await safe(opentargets.fetch_target(client, genes[0], prof["uniprot_id"]), "Open Targets")
    if len(genes) > 1:
        ot["other_ensembl_genes_not_queried"] = genes[1:]
    return ot


def _location_text(loc: dict) -> str:
    text = "; ".join(x for x in (loc["location"], loc["topology"], loc["orientation"]) if x)
    return f"{text} ({loc['molecule']})" if loc["molecule"] else text


def _compact_uniprot(p: dict) -> dict:
    text = " ".join(f["text"] for f in p["function"])
    return {
        "gene": p["gene"],
        "protein_name": p["protein_name"],
        "organism": p["organism"],
        "reviewed": p["reviewed"],
        "annotation_score": p["annotation_score"],
        "last_annotation_update": p["last_annotation_update"],
        "links": {"uniprot": p["source_url"]},
        "function": {
            "text": text[:FUNCTION_CHARS] + ("..." if len(text) > FUNCTION_CHARS else "") if text else None,
            "truncated": len(text) > FUNCTION_CHARS,
            "pubmed_count": len({x for f in p["function"] for x in f["pubmed"]}),
            "evidence": sorted({e for f in p["function"] for e in f["evidence"]}),
        },
        "subcellular_locations": [_location_text(loc) for loc in p["subcellular_locations"]],
        "diseases": {
            "count": len(p["diseases"]),
            "items": [
                {"name": d["name"], "acronym": d["acronym"], "mim": d["mim"],
                 "note": None if d["name"] else short(d["note"], 160),  # unnamed: only a note exists
                 "links": d["links"]}
                for d in p["diseases"]
            ],
        },
        "natural_variants": {
            "count": len(p["natural_variants"]),
            "changes": [v["change"] or v["feature_id"] for v in p["natural_variants"]],
        },
        "mutagenesis_count": len(p["mutagenesis"]),
        "ensembl_gene_ids": p["ensembl_gene_ids"],
        "chembl_target_ids": p["chembl_target_ids"],
    }


def _compact_ot(ot: dict, n: int) -> dict:
    if not ot.get("ok") or not ot.get("found"):
        keys = ("ok", "found", "note", "error", "ensembl_gene_id", "data_version")
        return {k: ot[k] for k in keys if ot.get(k) is not None}

    def slim(r):
        return {"name": r["name"], "disease_id": r["disease_id"], "score": r["score"],
                "datatype_scores": r["datatype_scores"], "link": r["links"].get("evidence")}

    assoc = ot["associations"]
    rows = assoc["rows"]
    oncology = [r for r in rows if r["is_oncology"]]
    names = ot["modality_names"]
    return {
        "ok": True,
        "found": True,
        "ensembl_gene_id": ot["ensembl_gene_id"],
        "data_version": ot["data_version"],
        "symbol": ot["symbol"],
        "uniprot_crosscheck": ot["uniprot_crosscheck"],
        "tractability": {
            names.get(m, m): {"yes": [x["label"] for x in flags if x["value"] is True],
                              "no": [x["label"] for x in flags if x["value"] is False]}
            for m, flags in ot["tractability"].items()
        },
        "associations": {
            "total": assoc["total"],
            "examined": assoc["examined"],
            "top_overall": [slim(r) for r in rows[:n]],
            "oncology_in_examined": len(oncology),
            "top_oncology": [slim(r) for r in oncology[:n]],
        },
        "links": ot["links"],
        **({"other_ensembl_genes_not_queried": ot["other_ensembl_genes_not_queried"]}
           if ot.get("other_ensembl_genes_not_queried") else {}),
    }


def compact(full: dict, max_diseases: int = TOP_DISEASES) -> dict:
    """Chat-sized summary of the full profile (the full one goes to the dossier)."""
    if not full.get("ok"):
        return full
    n = max(1, min(int(max_diseases), 50))
    return {
        "ok": True,
        "uniprot_id": full["uniprot_id"],
        "summary_note": (
            "Compact summary. The full result (all variants with descriptions and PubMed IDs, "
            "disease descriptions, every examined association with its datatype scores, all "
            "tractability flags) is in dossier_path."
        ),
        "uniprot": _compact_uniprot(full["uniprot"]),
        "open_targets": _compact_ot(full["open_targets"], n),
    }


async def get_protein_profile(uniprot_id: str, _client=None) -> dict:
    """FULL profile; server.py saves it to the dossier and sends compact() to the chat."""
    accession, err = parse_accession(uniprot_id)
    if err:
        return err

    own = _client is None
    client = _client or httpx.AsyncClient(timeout=60.0, follow_redirects=True)
    try:
        up = await uniprot.fetch_entry_json(client, accession)
        if not up["ok"]:
            return {**up, "uniprot_id": accession}
        prof = uniprot.parse_profile(up["json"])
        ot = await _open_targets(client, prof)
    finally:
        if own:
            await client.aclose()

    return {"ok": True, "uniprot_id": accession, "uniprot": prof, "open_targets": ot}
