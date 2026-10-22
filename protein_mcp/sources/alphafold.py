"""AlphaFold DB: predicted model and confidence for a UniProt protein."""
import re

import httpx

from protein_mcp.common import log

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

# >>> EDIT THIS LINE: an isoform model replaces the canonical one only if its mean
# >>> pLDDT is higher by at least this many points. <<<
AF_ISOFORM_MARGIN = 5.0


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


async def fetch_alphafold(client, accession: str) -> dict:
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

    # Prefer the canonical UniProt sequence: structures, papers and ligand data all
    # refer to it. An isoform model wins only if it is clearly better.
    best_overall = max(models, key=score)
    canonical = [m for m in models if isoform_of(m) is None]
    if not canonical:
        best = best_overall
        reason = (
            "No canonical-sequence model exists in AlphaFold DB for this accession; "
            "chose the model with the highest mean pLDDT."
        )
    else:
        best_canonical = max(canonical, key=score)
        if best_overall is best_canonical:
            best = best_canonical
            reason = "Canonical sequence model; it also has the highest mean pLDDT."
        elif score(best_overall) - score(best_canonical) >= AF_ISOFORM_MARGIN:
            best = best_overall
            reason = (
                f"Isoform model {best_overall.get('entryId')} chosen: its mean pLDDT "
                f"({score(best_overall)}) beats the canonical model's ({score(best_canonical)}) "
                f"by at least {AF_ISOFORM_MARGIN}."
            )
        else:
            best = best_canonical
            reason = (
                f"Canonical model preferred: the best isoform model ({best_overall.get('entryId')}, "
                f"pLDDT {score(best_overall)}) is better by less than {AF_ISOFORM_MARGIN} "
                f"than the canonical one ({score(best_canonical)})."
            )
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