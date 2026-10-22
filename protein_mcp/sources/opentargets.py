"""Open Targets Platform: target-disease associations and tractability (human targets only).

Validated against the real GraphQL API on 2026-10-20 (KRAS, ENSG00000133703, data 26.09):
- a VERSIONED Ensembl ID ("ENSG00000133703.15", as UniProt gives it) returns
  target: null with HTTP 200 and no errors, exactly like a missing target;
- associatedDiseases rows come ordered by score (descending); page size must be 0-3000;
- tractability lists every flag, true and false, per modality.
"""
import re

import httpx

from protein_mcp.common import log

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
OT_WEB = "https://platform.opentargets.org"

ASSOCIATION_POOL = 500          # top associations examined, by overall score (KRAS has 2262)
ONCOLOGY_AREA = "MONDO_0045024"  # Open Targets therapeutic area "cancer or benign tumor"
MODALITIES = {"SM": "small molecule", "AB": "antibody", "PR": "PROTAC", "OC": "other clinical"}

_ENSG_RE = re.compile(r"ENSG\d{11}")

_TARGET_QUERY = """
query ($id: String!, $size: Int!) {
  meta { dataVersion { year month iteration } }
  target(ensemblId: $id) {
    id
    approvedSymbol
    approvedName
    proteinIds { id source }
    tractability { label modality value }
    associatedDiseases(page: {index: 0, size: $size}) {
      count
      rows {
        score
        disease { id name therapeuticAreas { id name } }
        datatypeScores { id score }
      }
    }
  }
}
"""


def _num(x):
    return round(x, 3) if isinstance(x, (int, float)) else None


def _parse_row(row: dict, gene_id: str) -> dict:
    d = row.get("disease") or {}
    did = d.get("id")
    areas = d.get("therapeuticAreas") or []
    return {
        "disease_id": did,
        "name": d.get("name"),
        "score": _num(row.get("score")),
        "therapeutic_areas": [a.get("name") for a in areas if a.get("name")],
        "is_oncology": any(a.get("id") == ONCOLOGY_AREA for a in areas),
        # only the evidence types that have data are returned
        "datatype_scores": {
            s["id"]: _num(s.get("score")) for s in row.get("datatypeScores") or [] if s.get("id")
        },
        "links": {
            "disease": f"{OT_WEB}/disease/{did}",
            "evidence": f"{OT_WEB}/evidence/{gene_id}/{did}",
        } if did else {},
    }


async def fetch_target(client, ensembl_id: str, accession: str, pool: int = ASSOCIATION_POOL) -> dict:
    gene_id = re.sub(r"\.\d+$", "", (ensembl_id or "").strip())
    if not _ENSG_RE.fullmatch(gene_id):
        return {"ok": False, "error": f"'{ensembl_id}' is not an Ensembl human gene ID (ENSG + 11 digits)."}
    try:
        resp = await client.post(
            OT_GRAPHQL, json={"query": _TARGET_QUERY, "variables": {"id": gene_id, "size": pool}}
        )
        resp.raise_for_status()
        body = resp.json()
    except (httpx.HTTPError, ValueError) as exc:
        log.error("Open Targets request failed: %s", exc)
        return {"ok": False, "error": f"Open Targets request failed: {exc}"}
    if body.get("errors"):
        message = (body["errors"][0] or {}).get("message", "")
        return {"ok": False, "error": f"Open Targets error: {message[:200]}"}

    data = body.get("data") or {}
    dv = (data.get("meta") or {}).get("dataVersion") or {}
    version = ".".join(str(v) for v in (dv.get("year"), dv.get("month")) if v) or None
    target = data.get("target")
    if target is None:
        return {
            "ok": True, "found": False, "ensembl_gene_id": gene_id, "data_version": version,
            "note": "sin registro en Open Targets para este gen.",
        }

    tractability = {}
    for t in target.get("tractability") or []:
        tractability.setdefault(t.get("modality"), []).append({"label": t.get("label"), "value": t.get("value")})
    assoc = target.get("associatedDiseases") or {}
    rows = [_parse_row(r, gene_id) for r in assoc.get("rows") or []]
    swissprot = [p.get("id") for p in target.get("proteinIds") or [] if p.get("source") == "uniprot_swissprot"]
    return {
        "ok": True,
        "found": True,
        "ensembl_gene_id": gene_id,
        "data_version": version,
        "symbol": target.get("approvedSymbol"),
        "name": target.get("approvedName"),
        "swissprot_ids": swissprot,
        "uniprot_crosscheck": accession in swissprot,  # False = Open Targets maps this gene to another protein
        "tractability": tractability,
        "modality_names": MODALITIES,
        "associations": {
            "total": assoc.get("count"),
            "examined": len(rows),
            "selection": f"top {len(rows)} of {assoc.get('count')} by overall association score",
            "rows": rows,
        },
        "links": {"target": f"{OT_WEB}/target/{gene_id}"},
    }
