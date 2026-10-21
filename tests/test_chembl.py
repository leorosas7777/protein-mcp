"""Offline test for chembl.py using shapes copied from the REAL KRAS responses
(2026-10-17): string numbers, null pref_name, null molecule_properties, a PPI
target next to the single-protein target, assay_variant_mutation set.
Run: uv run python tests/test_chembl.py
"""
import asyncio
import json
import sys
from pathlib import Path

import httpx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # project root
import chembl  # noqa: E402


def act(aid, mol, val, pchembl, mut="G12V", rel="=", units="nM", assay="CHEMBL3223233", **kw):
    d = {
        "activity_id": aid, "assay_chembl_id": assay, "assay_type": "B",
        "assay_variant_mutation": mut, "assay_variant_accession": "P01116",
        "assay_description": f"Inhibition of recombinant K-Ras {mut or 'protein'} mutant by HTRF assay",
        "document_chembl_id": "CHEMBL3217533", "document_journal": "Medchemcomm",
        "document_year": 2013, "molecule_chembl_id": mol, "parent_molecule_chembl_id": mol,
        "pchembl_value": pchembl, "potential_duplicate": 0, "standard_relation": rel,
        "standard_type": "IC50", "standard_units": units, "standard_value": val,
        "data_validity_comment": None, "canonical_smiles": "CC",
    }
    d.update(kw)
    return d


ACTS = [
    act(1, "CHEMBL3218635", "1200.0", "5.92"),
    act(2, "CHEMBL3218635", "3000.0", "5.52"),
    act(3, "CHEMBL3218635", "800.0", "6.10", mut="G12D"),
    act(4, "CHEMBL3218636", "10000.0", "5.00"),
    act(5, "CHEMBL3218636", "9000.0", "5.05", rel=">"),            # censored
    act(6, "CHEMBL3218636", "50.0", "7.3", potential_duplicate=1),  # duplicate
    act(7, "CHEMBL3218637", "20.0", "7.7", mut=None, assay="CHEMBL_LOWCONF"),
    act(8, "CHEMBL3218637", "5.0", "8.3", mut=None),
    act(8, "CHEMBL3218637", "5.0", "8.3", mut=None),                # repeated id
    # Real KRAS defect: label G12C but the assay text says G12D (top-ranked row in the live run)
    act(9, "CHEMBL3218638", "0.16", "9.8", mut="G12C", standard_type="Kd",
        assay_description="Binding affinity to KRAS G12D mutant (unknown origin) by SPR analysis"),
    # Real KRAS defect: cell assay typed 'B' (binding); must not be pooled with biochemical rows
    act(10, "CHEMBL3218635", "50.0", "7.3", assay="CHEMBL_CELL",
        assay_description="Inhibition of KRAS G12V mutant in human NCI-H358 cells by in-cell western"),
]
MOLECULE = {
    "molecule_chembl_id": "CHEMBL3218635", "pref_name": None, "max_phase": None,
    "molecule_type": "Small molecule", "first_approval": None, "chemical_probe": 0,
    "molecule_properties": {"alogp": None, "full_mwt": "1396.59", "mw_freebase": "1396.59",
                            "psa": None, "hba": None, "hbd": None, "num_ro5_violations": None},
    "molecule_structures": {"canonical_smiles": "CC", "standard_inchi_key": "MOGVORJDAHDGBV-GTTJADIHSA-N"},
}
TARGETS = {"targets": [
    {"target_chembl_id": "CHEMBL2189121", "target_type": "SINGLE PROTEIN", "pref_name": "GTPase KRas",
     "organism": "Homo sapiens", "target_components": [{"accession": "P01116"}]},
    {"target_chembl_id": "CHEMBL5291977", "target_type": "PROTEIN-PROTEIN INTERACTION",
     "pref_name": "GTPase KRas/RAF1", "organism": "Homo sapiens",
     "target_components": [{"accession": "P01116"}, {"accession": "P04049"}]},
]}


def handler(req: httpx.Request):
    p, q = req.url.path, dict(req.url.params)
    if p.endswith("/target.json"):
        return httpx.Response(200, json=TARGETS)
    if p.endswith("/activity.json"):
        assert q["target_chembl_id"] == "CHEMBL2189121", q   # never the PPI target
        return httpx.Response(200, json={"activities": ACTS, "page_meta": {"total_count": len(ACTS)}})
    if p.endswith("/assay.json"):
        return httpx.Response(200, json={"assays": [
            {"assay_chembl_id": "CHEMBL3223233", "confidence_score": 9, "assay_cell_type": None},
            {"assay_chembl_id": "CHEMBL_LOWCONF", "confidence_score": 4, "assay_cell_type": "HCT116"},
            {"assay_chembl_id": "CHEMBL_CELL", "confidence_score": 9, "assay_cell_type": "NCI-H358"},
        ], "page_meta": {"total_count": 3}})
    if p.endswith("/mechanism.json"):
        return httpx.Response(200, json={"mechanisms": [
            {"molecule_chembl_id": "CHEMBL_DRUG_NO_DATA", "mechanism_of_action": "KRAS G12C inhibitor",
             "action_type": "INHIBITOR", "max_phase": 4}]})
    if p.endswith("/molecule.json"):
        ids = q["molecule_chembl_id__in"].split(",")
        return httpx.Response(200, json={"molecules": [{**MOLECULE, "molecule_chembl_id": i} for i in ids]})
    if p.endswith("/document.json"):
        return httpx.Response(200, json={"documents": [
            {"document_chembl_id": "CHEMBL3217533", "pubmed_id": 24000000, "doi": "10.1039/x", "journal": "Medchemcomm", "year": 2013}]})
    return httpx.Response(404)


async def run(**kw):
    async with httpx.AsyncClient(transport=httpx.MockTransport(handler)) as c:
        return await chembl.get_ligands("P01116", _client=c, **kw)


async def main():
    r = await run()
    a = r["activities"]
    assert r["target"]["chembl_id"] == "CHEMBL2189121" and len(r["other_targets_ignored"]) == 1
    # used: ids 1,2,3,4,8,10 (5 censored, 6 duplicate, 7 low confidence, 9 label conflict, repeated 8 deduped)
    assert a["used"] == 6, a
    assert a["excluded"] == {"potential_duplicate": 1, "below_min_confidence": 1,
                             "censored_relation_not_equal": 1,
                             "mutation_label_conflicts_description": 1}, a["excluded"]
    assert [f["molecule_chembl_id"] for f in r["label_conflicts_flagged"]] == ["CHEMBL3218638"]
    ids = [l["chembl_id"] for l in r["ligands"]]
    assert "CHEMBL3218638" not in ids, ids                          # conflicting row never ranks
    assert ids[0] == "CHEMBL_DRUG_NO_DATA" or ids[0] == "CHEMBL3218637", ids
    k = next(l for l in r["ligands"] if l["chembl_id"] == "CHEMBL3218635")
    muts = {(p["mutation"], p["assay_class"]): p for p in k["potency"]}
    # never merged across mutations NOR across cellular/biochemical
    assert set(muts) == {("G12V", "binding"), ("G12V", "cellular"), ("G12D", "binding")}, muts
    assert muts[("G12V", "cellular")]["cell_types"] == ["NCI-H358"]
    g = muts[("G12V", "binding")]
    assert (g["n_measurements"], g["best_nM"], g["worst_nM"]) == (2, 1200.0, 3000.0)
    assert k["name"] is None and k["properties"]["alogp"] is None    # nulls stay null
    assert g["documents"][0]["links"]["pubmed"].endswith("/24000000/")
    assert r["mechanism_compounds_without_qualifying_activity"] == ["CHEMBL_DRUG_NO_DATA"]

    r2 = await run(mutation="G12")                                   # bare position
    assert {p["mutation"] for l in r2["ligands"] for p in l["potency"]} == {"G12V", "G12D"}
    r3 = await run(mutation="G12D")
    assert [l["chembl_id"] for l in r3["ligands"]] == ["CHEMBL3218635"]
    r4 = await run(mutation="G12C")
    assert r4["ligands_returned"] == 0
    bad = await chembl.get_ligands("not an id")
    assert bad["ok"] is False
    print(json.dumps({"selection": r["selection"], "excluded": a["excluded"],
                      "confidence": r["confidence"]}, indent=1))
    print("ALL OK")


asyncio.run(main())