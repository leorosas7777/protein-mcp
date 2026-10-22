"""Offline test for chembl.py using shapes copied from the REAL KRAS responses
(2026-10-17): string numbers, null pref_name, null molecule_properties, a PPI
target next to the single-protein target, assay_variant_mutation set.
Run: uv run python tests/test_chembl.py
"""
import asyncio
import json
from pathlib import Path

import httpx

from protein_mcp.sources import chembl


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
# Real shapes (ChEMBL mechanism endpoint, 2026-10-21): EGFR records erlotinib on its
# salt (molecule CHEMBL1079742, parent CHEMBL553); PCSK9 antibodies have no activity
# and come back with molecule_structures/molecule_properties = null.
MECHANISMS = [
    {"molecule_chembl_id": "CHEMBL_SALT", "parent_molecule_chembl_id": "CHEMBL3218636",
     "mechanism_of_action": "GTPase KRas inhibitor", "action_type": "INHIBITOR", "max_phase": 4},
    {"molecule_chembl_id": "CHEMBL_ANTIBODY", "parent_molecule_chembl_id": "CHEMBL_ANTIBODY",
     "mechanism_of_action": "GTPase KRas inhibitor", "action_type": "INHIBITOR", "max_phase": 4},
]
ANTIBODY = {"molecule_chembl_id": "CHEMBL_ANTIBODY", "pref_name": "EVOLOCUMAB", "max_phase": "4.0",
            "molecule_type": "Antibody", "first_approval": 2015, "chemical_probe": 0,
            "molecule_properties": None, "molecule_structures": None}
TARGETS = {"targets": [
    {"target_chembl_id": "CHEMBL2189121", "target_type": "SINGLE PROTEIN", "pref_name": "GTPase KRas",
     "organism": "Homo sapiens", "target_components": [{"accession": "P01116"}]},
    {"target_chembl_id": "CHEMBL5291977", "target_type": "PROTEIN-PROTEIN INTERACTION",
     "pref_name": "GTPase KRas/RAF1", "organism": "Homo sapiens",
     "target_components": [{"accession": "P01116"}, {"accession": "P04049"}]},
]}


ASSAYS = [
    {"assay_chembl_id": "CHEMBL3223233", "confidence_score": 9, "assay_cell_type": None},
    {"assay_chembl_id": "CHEMBL_LOWCONF", "confidence_score": 4, "assay_cell_type": "HCT116"},
    {"assay_chembl_id": "CHEMBL_CELL", "confidence_score": 9, "assay_cell_type": "NCI-H358"},
]
FAIL = set()   # endpoints that answer HTTP 503: "mechanism", "mechanism_activities", "molecule"


def handler(req: httpx.Request):
    p, q = req.url.path, dict(req.url.params)
    if p.endswith("/target.json"):
        return httpx.Response(200, json=TARGETS)
    if p.endswith("/activity.json"):
        assert q["target_chembl_id"] == "CHEMBL2189121", q   # never the PPI target
        if "parent_molecule_chembl_id__in" in q:            # activities of mechanism compounds
            if "mechanism_activities" in FAIL:
                return httpx.Response(503)
            ids = q["parent_molecule_chembl_id__in"].split(",")
            rows = [a for a in ACTS if a["parent_molecule_chembl_id"] in ids]
            return httpx.Response(200, json={"activities": rows, "page_meta": {"total_count": len(rows)}})
        # ranked download: most potent first, activity_id as tie-breaker (real API order)
        assert req.url.params.get_list("order_by") == ["-pchembl_value", "activity_id"], req.url
        rows = sorted(ACTS, key=lambda a: (-float(a["pchembl_value"]), a["activity_id"]))
        off, lim = int(q["offset"]), int(q["limit"])
        return httpx.Response(200, json={"activities": rows[off:off + lim],
                                         "page_meta": {"total_count": len(ACTS)}})
    if p.endswith("/assay.json"):   # paged like the real endpoint (EGFR: 5966 assays)
        off, lim = int(q["offset"]), int(q["limit"])
        return httpx.Response(200, json={"assays": ASSAYS[off:off + lim],
                                         "page_meta": {"total_count": len(ASSAYS)}})
    if p.endswith("/mechanism.json"):
        if "mechanism" in FAIL:
            return httpx.Response(503)
        return httpx.Response(200, json={"mechanisms": MECHANISMS, "page_meta": {"total_count": len(MECHANISMS)}})
    if p.endswith("/molecule.json"):
        if "molecule" in FAIL:
            return httpx.Response(503)
        ids = q["molecule_chembl_id__in"].split(",")
        return httpx.Response(200, json={"molecules": [
            ANTIBODY if i == "CHEMBL_ANTIBODY" else {**MOLECULE, "molecule_chembl_id": i} for i in ids]})
    if p.endswith("/document.json"):
        return httpx.Response(200, json={"documents": [
            {"document_chembl_id": "CHEMBL3217533", "pubmed_id": 24000000, "doi": "10.1039/x", "journal": "Medchemcomm", "year": 2013}]})
    return httpx.Response(404)


async def no_sleep(*_):
    pass


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
    # mechanism recorded on a salt still reaches the parent's ligand row, and ranks first (phase 4)
    assert ids[0] == "CHEMBL3218636", ids
    salt = r["ligands"][0]
    assert salt["max_phase"] == 4 and salt["mechanism"]["mechanism_of_action"] == "GTPase KRas inhibitor"
    k = next(l for l in r["ligands"] if l["chembl_id"] == "CHEMBL3218635")
    muts = {(p["mutation"], p["assay_class"]): p for p in k["potency"]}
    # never merged across mutations NOR across cellular/biochemical
    assert set(muts) == {("G12V", "binding"), ("G12V", "cellular"), ("G12D", "binding")}, muts
    assert muts[("G12V", "cellular")]["cell_types"] == ["NCI-H358"]
    g = muts[("G12V", "binding")]
    assert (g["n_measurements"], g["best_nM"], g["worst_nM"]) == (2, 1200.0, 3000.0)
    assert k["name"] is None and k["properties"]["alogp"] is None    # nulls stay null
    assert g["documents"][0]["links"]["pubmed"].endswith("/24000000/")
    # antibody without IC50/Ki/Kd: listed with name, type and phase, never as a bare id
    assert r["mechanism_compounds_without_qualifying_activity"] == [{
        "chembl_id": "CHEMBL_ANTIBODY", "name": "EVOLOCUMAB", "molecule_record_retrieved": True,
        "link": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL_ANTIBODY/",
        "molecule_type": "Antibody", "max_phase": 4.0, "first_approval": 2015,
        "mechanism_of_action": "GTPase KRas inhibitor", "action_type": "INHIBITOR"}]
    assert "1 more compounds with a recorded mechanism" in r["selection"]
    assert a["truncated"] is False and a["fetched_for_mechanism_compounds"] == 0
    assert a["ranked_download_cutoff_pchembl"] is None

    # truncated download: keeps the most potent rows, then fetches every activity of
    # mechanism compounds so the phase-4 drug below the cutoff is not lost
    page, cap = chembl.PAGE, chembl.MAX_ACTIVITIES
    chembl.PAGE = chembl.MAX_ACTIVITIES = 4
    try:
        rt = await run()
    finally:
        chembl.PAGE, chembl.MAX_ACTIVITIES = page, cap
    at = rt["activities"]
    assert at["truncated"] and at["total_in_chembl"] == len(ACTS), at
    assert at["ranked_download_cutoff_pchembl"] == 7.7, at          # top 4: ids 9, 8, 8, 7
    assert (at["fetched"], at["fetched_for_mechanism_compounds"]) == (7, 3), at   # + ids 4, 5, 6
    tids = [l["chembl_id"] for l in rt["ligands"]]
    assert tids[0] == "CHEMBL3218636" and "CHEMBL3218635" not in tids, tids   # 35 is below the cutoff
    assert any(n.startswith("Truncated: the 4 most potent of 11") for n in rt["notes"]), rt["notes"]
    assert rt["mechanism_compounds_without_qualifying_activity"][0]["name"] == "EVOLOCUMAB"
    assert rt["mechanism_compounds_unchecked"] == []
    assert at["mechanism_compounds_download"] == "complete" and r["activities"]["mechanism_compounds_download"] is None
    assert rt["ligands_with_qualifying_data_is_lower_bound"] is True
    assert rt["selection"].startswith("At least "), rt["selection"]
    assert r["ligands_with_qualifying_data_is_lower_bound"] is False and not r["selection"].startswith("At least")
    assert r["mechanisms"] == {"retrieved": True, "records_read": len(MECHANISMS), "records_total": len(MECHANISMS)}

    # truncated AND the extra download fails: below-cutoff drugs are "unchecked",
    # never "no qualifying activity", and the note does not claim the download
    orig_sleep = asyncio.sleep
    asyncio.sleep = no_sleep   # _get retries with sleeps
    chembl.PAGE = chembl.MAX_ACTIVITIES = 4
    FAIL.add("mechanism_activities")
    try:
        rf = await run()
    finally:
        chembl.PAGE, chembl.MAX_ACTIVITIES = page, cap
        FAIL.clear()
        asyncio.sleep = orig_sleep
    assert rf["mechanism_compounds_without_qualifying_activity"] == [], rf["selection"]
    assert {x["chembl_id"] for x in rf["mechanism_compounds_unchecked"]} == {"CHEMBL3218636", "CHEMBL_ANTIBODY"}
    trunc_note = next(n for n in rf["notes"] if n.startswith("Truncated"))
    assert "FAILED" in trunc_note and "plus every activity" not in trunc_note, trunc_note
    assert rf["activities"]["mechanism_compounds_download"] == "failed"

    # same failure when every mechanism compound is already above the cutoff: unchecked
    # is empty, but the failure is still declared (the statistics are partial)
    asyncio.sleep = no_sleep
    chembl.PAGE = chembl.MAX_ACTIVITIES = 4
    FAIL.add("mechanism_activities")
    saved_mechs = MECHANISMS[:]
    MECHANISMS[:] = [{**MECHANISMS[0], "molecule_chembl_id": "CHEMBL3218637",
                      "parent_molecule_chembl_id": "CHEMBL3218637"}]   # pChEMBL 8.3: inside the top 4
    try:
        rin = await run()
    finally:
        chembl.PAGE, chembl.MAX_ACTIVITIES = page, cap
        FAIL.clear()
        MECHANISMS[:] = saved_mechs
        asyncio.sleep = orig_sleep
    assert rin["mechanism_compounds_unchecked"] == [] and rin["activities"]["mechanism_compounds_download"] == "failed"
    in_note = next(n for n in rin["notes"] if n.startswith("Truncated"))
    assert "FAILED" in in_note and "statistics of the others only count" in in_note, in_note

    # mechanism endpoint fails: the drug lists are None ("not checked"), never [] ("none"),
    # and neither selection nor the truncation note claim the mechanism download
    asyncio.sleep = no_sleep
    chembl.PAGE = chembl.MAX_ACTIVITIES = 4
    FAIL.add("mechanism")
    try:
        rmf = await run()
    finally:
        chembl.PAGE, chembl.MAX_ACTIVITIES = page, cap
        FAIL.clear()
        asyncio.sleep = orig_sleep
    assert rmf["ok"] is True and rmf["mechanisms"]["retrieved"] is False, rmf["mechanisms"]
    for field in ("mechanism_compounds_without_qualifying_activity", "mechanism_compounds_other_variants_only",
                  "mechanism_compounds_unchecked"):
        assert rmf[field] is None, (field, rmf[field])
    assert "WARNING: the ChEMBL mechanism lookup failed" in rmf["selection"], rmf["selection"]
    assert any(n.startswith("The ChEMBL mechanism lookup FAILED") for n in rmf["notes"]), rmf["notes"]
    mf_note = next(n for n in rmf["notes"] if n.startswith("Truncated"))
    assert "plus every activity" not in mf_note and "mechanism lookup failed" in mf_note, mf_note
    assert rmf["activities"]["mechanism_compounds_download"] == "mechanism_lookup_failed"
    assert not any("are listed in mechanism_compounds_without_qualifying_activity" in n for n in rmf["notes"])
    assert chembl.compact(rmf)["mechanism_compounds_without_qualifying_activity"] is None

    # molecule endpoint fails: names are "not retrieved", not "absent from ChEMBL"
    asyncio.sleep = no_sleep
    FAIL.add("molecule")
    try:
        rm = await run()
    finally:
        FAIL.clear()
        asyncio.sleep = orig_sleep
    assert all(l["molecule_record_retrieved"] is False and l["name"] is None for l in rm["ligands"])
    assert rm["mechanism_compounds_without_qualifying_activity"][0]["molecule_record_retrieved"] is False
    assert any("molecule_record_retrieved is false" in n for n in rm["notes"]), rm["notes"]
    assert chembl.compact(rm)["ligands"][0]["molecule_record_retrieved"] is False
    assert "molecule_record_retrieved" not in chembl.compact(r)["ligands"][0]   # only when false

    # assay records beyond what was read: those activities are excluded and counted,
    # never passed unfiltered (EGFR: 966 of 5966 assays lost to a 5000 cap)
    cap_assays = chembl.MAX_ASSAYS
    chembl.PAGE = chembl.MAX_ASSAYS = 2
    try:
        ra = await run()
    finally:
        chembl.PAGE, chembl.MAX_ASSAYS = page, cap_assays
    assert ra["activities"]["excluded"]["assay_record_missing"] == 1, ra["activities"]   # id 10 (CHEMBL_CELL)
    assert ra["confidence"]["assay_records_read"] == 2 and ra["confidence"]["assay_records_total"] == 3
    assert any(n.startswith("Only 2 of 3 assay records") for n in ra["notes"]), ra["notes"]
    assert r["confidence"]["assay_records_read"] == 3 and not any("assay records were read" in n for n in r["notes"])

    # ChEMBL max_phase -1 = clinical, phase unknown: below phase 0.5, above preclinical
    order = sorted([None, -1.0, 0.5, 4.0, 1.0], key=lambda x: -chembl._phase_rank(x))
    assert order == [4.0, 1.0, 0.5, -1.0, None], order

    r2 = await run(mutation="G12")                                   # bare position
    assert {p["mutation"] for l in r2["ligands"] for p in l["potency"]} == {"G12V", "G12D"}
    r3 = await run(mutation="G12D")
    assert [l["chembl_id"] for l in r3["ligands"]] == ["CHEMBL3218635"]
    # phase-4 drug with data only on G12V: "other variants only", NOT "no qualifying activity"
    # (live EGFR C797S listed lapatinib and dacomitinib as having no IC50)
    assert [x["chembl_id"] for x in r3["mechanism_compounds_other_variants_only"]] == ["CHEMBL3218636"]
    assert [x["chembl_id"] for x in r3["mechanism_compounds_without_qualifying_activity"]] == ["CHEMBL_ANTIBODY"]
    assert "1 have data only for other variants" in r3["selection"], r3["selection"]
    assert r["mechanism_compounds_other_variants_only"] == []                  # no filter: never used
    r4 = await run(mutation="G12C")
    assert r4["ligands_returned"] == 0
    bad = await chembl.get_ligands("not an id")
    assert bad["ok"] is False
    assert chembl.compact(bad) is bad

    # compact view for the chat: no SMILES/assay texts, one source link per potency row
    c = chembl.compact(r)
    text = json.dumps(c)
    assert "canonical_smiles" not in text and "smiles" not in text and "assay_examples" not in text
    ck = next(l for l in c["ligands"] if l["chembl_id"] == "CHEMBL3218635")
    cg = next(p for p in ck["potency"] if (p["mutation"], p["assay_class"]) == ("G12V", "binding"))
    assert cg["source"].endswith("/24000000/") and cg["n_documents"] == 1 and cg["best_nM"] == 1200.0
    assert ck["potency_groups_total"] == 3 and len(ck["potency"]) == 3
    big = {**r, "ligands": [{**k, "potency": k["potency"] * 3}]}       # 9 groups -> 3 in chat
    bk = chembl.compact(big)["ligands"][0]
    assert bk["potency_groups_total"] == 9 and len(bk["potency"]) == chembl.COMPACT_POTENCY_ROWS
    assert c["label_conflicts_flagged"]["listed_in_dossier"] == 1
    assert c["activities"] == r["activities"] and c["selection"] == r["selection"]
    assert len(json.dumps(c)) < len(json.dumps(r))
    print(json.dumps({"selection": r["selection"], "excluded": a["excluded"],
                      "confidence": r["confidence"]}, indent=1))
    print("ALL OK")


asyncio.run(main())