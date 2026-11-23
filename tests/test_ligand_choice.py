"""Offline test for select_ligands (protein_mcp/files/ligand_choice.py) and for the compound path
of dock_ligand (Phase 3, step 6b). Same fixtures as test_dock: the real cleaned manifest, the real
mmCIFs, real Meeko, and Vina replaced by the stand-in that returns the real 1.2.7 output.

Every check has a control that must fail if the code stops doing the thing (verified by hand with
the mutations named in the comments). Run: uv run python tests/test_ligand_choice.py
"""
import asyncio
import json
import sys
from pathlib import Path

from protein_mcp.files import choice, dock, ligand_choice

sys.path.insert(0, str(Path(__file__).parent))
from test_choice import ACC  # noqa: E402
from test_dock import calls, make_fake_vina, make_prep  # noqa: E402

# the real canonical SMILES of the 6 compounds of step 6b, from ChEMBL on 2026-11-17
VEM = "CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(-c4ccc(Cl)cc4)cc23)c1F"
SORA = "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
LONG = "C" * 40                      # ~45 A: does not fit any box of this project
SALT = "CC(=O)Oc1ccccc1C(=O)O.[Na+]"

SEL = [{"id": "CHEMBL1229517", "name": "VEMURAFENIB", "level": "T1", "smiles": VEM},
       {"id": "CHEMBL1336", "name": "SORAFENIB", "level": "T3", "smiles": SORA}]


def manifest(prep: Path) -> dict:
    return json.loads((prep / "manifest.json").read_text(encoding="utf-8"))


def test_box_clearance_matches_the_measurement():
    """The clearance of the 6 compounds against the real box of 8C7Y TXV U, as measured on
    2026-11-17 with the separate script. This ties the tool to that measurement: if the count
    changes, these numbers move and the test fails.
    Control that must fail: any of the six figures changed by 0.01 (checked by hand)."""
    box = {"center": [0.02, 110.789, 15.175], "size": [21.8, 16.27, 22.0]}
    fit = dock.box_clearance(VEM, box)
    assert fit["min_clearance_angstrom"] == 4.2, fit
    assert fit["conformers_tried"] == 5 and fit["seeds_without_conformer"] == []
    assert fit["meets_threshold"] and not fit["box_limited"], fit
    fit = dock.box_clearance(SORA, box)
    assert fit["min_clearance_angstrom"] == 1.99, fit
    # 1.99 is above the 1.5 A minimum the user set on 2026-11-18 and below the 3.0 A flag
    assert fit["meets_threshold"] and fit["box_limited"], fit
    assert fit["limiting_box_axis_angstrom"] == 22.0 and fit["limiting_ligand_extent_angstrom"] == 20.01, fit
    # the long chain does not fit at all
    fit = dock.box_clearance(LONG, box)
    assert not fit["fits"] and not fit["meets_threshold"] and fit["min_clearance_angstrom"] < 0, fit
    # a box smaller than the compound in one axis only: the limiting axis is named, not the first one
    fit = dock.box_clearance(ASPIRIN, {"center": [0, 0, 0], "size": [30.0, 30.0, 5.0]})
    assert fit["limiting_box_axis_angstrom"] == 5.0, fit


def test_refusals_write_nothing():
    """A refused selection must leave the manifest byte for byte as it was.
    Control that must fail: writing the manifest before validating the rows."""
    prep, _ = make_prep()
    before = (prep / "manifest.json").read_bytes()
    bad = [
        ([{"id": "X", "smiles": SALT, "level": "T3"}], "2 components"),
        ([{"id": "X", "smiles": "not a smiles", "level": "T3"}], "could not read"),
        ([{"id": "X", "smiles": VEM, "level": "T2"}], "must be one of T1, T3"),
        ([{"id": "", "smiles": VEM, "level": "T3"}], "'id' is required"),
        ([{"id": "X", "level": "T3"}], "'smiles' is required"),
        ([{"id": "X", "smiles": VEM, "level": "T3"},
          {"id": "x", "smiles": SORA, "level": "T3"}], "already used by another compound"),
        ([{"id": "A", "smiles": VEM, "level": "T3"},
          {"id": "B", "smiles": SALT, "level": "T3"}], "components"),
    ]
    for compounds, text in bad:
        r = ligand_choice.select_ligands(ACC, compounds, reason="a reason")
        assert not r["ok"], (compounds, r)
        assert any(text in e for e in r["refused"]), (text, r["refused"])
        assert (prep / "manifest.json").read_bytes() == before, text
    # one bad row refuses the whole selection, and says which ones had passed so far
    r = ligand_choice.select_ligands(ACC, bad[-1][0], reason="a reason")
    assert r["accepted_so_far"] == ["A"] and "1 of 2" in r["error"], r
    # no reason: the valid selection is shown and nothing is written (as in select_structure)
    r = ligand_choice.select_ligands(ACC, SEL)
    assert not r["ok"] and "reason is required" in r["error"]
    assert r["valid_selection"] == ["CHEMBL1229517 (VEMURAFENIB, T1)", "CHEMBL1336 (SORAFENIB, T3)"], r
    assert (prep / "manifest.json").read_bytes() == before
    assert manifest(prep).get("ligand_choice") is None
    # a reason longer than the limit, and too many compounds
    assert "longer than" in ligand_choice.select_ligands(ACC, SEL, reason="x" * 1001)["error"]
    many = [{"id": f"C{i}", "smiles": ASPIRIN, "level": "T3"} for i in range(51)]
    assert "more than the 50" in ligand_choice.select_ligands(ACC, many, reason="r")["error"]
    assert (prep / "manifest.json").read_bytes() == before
    return prep


def test_records_the_selection():
    prep, _ = make_prep()
    r = ligand_choice.select_ligands(ACC, SEL, reason="the researcher picked them")
    assert r["ok"] and r["recorded"] == "2 of 2 compounds given", r
    assert r["levels"] == {"T1": ["CHEMBL1229517"], "T3": ["CHEMBL1336"]}, r["levels"]
    assert r["heavy_atom_range"] == [32, 33], r["heavy_atom_range"]   # sorafenib 32, vemurafenib 33
    rec = manifest(prep)["ligand_choice"]
    assert rec["reason"] == "the researcher picked them" and rec["recorded_by"] == "select_ligands"
    assert [c["slug"] for c in rec["compounds"]] == ["chembl1229517", "chembl1336"]
    vem = rec["compounds"][0]
    assert len(vem["inchikey"]) == 27 and vem["inchikey"].count("-") == 2, vem["inchikey"]
    assert vem["heavy_atoms"] == 33 and vem["values"] == "calculated"
    assert vem["canonical_smiles"] and vem["smiles"] == VEM
    # a second selection supersedes the first and keeps it in the history with the date
    r2 = ligand_choice.select_ligands(ACC, [SEL[0]], reason="solo el ancla")
    assert r2["ok"] and r2["previous_selections"] == 1 and "2 compounds recorded" in r2["replaced_selection"]
    m = manifest(prep)
    assert len(m["ligand_choice"]["compounds"]) == 1
    assert m["ligand_choice_history"][0]["reason"] == "the researcher picked them"
    assert m["ligand_choice_history"][0]["superseded_at"]
    # carry_over: prepare_structures replaces the manifest, the selection does not depend on it
    new = {"created_at": "2026-11-19T00:00:00Z"}
    status = ligand_choice.carry_over(m, new)
    assert status["kept"] == "1 compounds" and new["ligand_choice"]["reason"] == "solo el ancla"
    assert new["ligand_choice"]["revalidated_at"] == "2026-11-19T00:00:00Z"
    assert len(new["ligand_choice_history"]) == 1
    # control that must fail: carry_over on a manifest without a selection must not invent one
    empty = {}
    assert ligand_choice.carry_over({}, empty) is None and empty["ligand_choice"] is None
    return prep


def test_dock_refuses_without_a_selection(record: Path):
    """Nothing is docked from a string typed in the chat (ESTADO 7.4, decision (b) 1).
    Control that must fail: resolving `compound` as a SMILES when nothing is recorded."""
    prep, _ = make_prep()
    before = len(calls(record))
    r = asyncio.run(dock.start_dock(ACC, "3OG7", "032 C", compound=VEM))
    assert not r["ok"] and "select_ligands" in r["error"], r
    assert "no compound selection is recorded" in r["error"]
    ligand_choice.select_ligands(ACC, SEL, reason="r")
    # a compound that is not in the selection: refused, and the recorded ones are listed
    r = asyncio.run(dock.start_dock(ACC, "3OG7", "032 C", compound="CHEMBL25"))
    assert not r["ok"] and "not in the recorded selection" in r["error"], r
    assert r["recorded_compounds"] == ["CHEMBL1229517 (VEMURAFENIB)", "CHEMBL1336 (SORAFENIB)"], r
    # conformer 'crystal' with a compound: there is no crystal pose of it
    r = asyncio.run(dock.start_dock(ACC, "3OG7", "032 C", conformer="crystal", compound="CHEMBL1336"))
    assert not r["ok"] and "has no crystal pose" in r["error"], r
    # a compound that does not fit the box: refused naming the box axis and the figure
    ligand_choice.select_ligands(ACC, [{"id": "LONG", "smiles": LONG, "level": "T3"}], reason="r")
    r = asyncio.run(dock.start_dock(ACC, "3OG7", "032 C", compound="LONG"))
    assert not r["ok"] and "below the declared minimum" in r["error"], r
    assert "box axis" in r["error"] and "is NOT enlarged" in r["error"]
    assert r["box_fit"]["min_clearance_angstrom"] < 0
    assert len(calls(record)) == before, "Vina was started for a refused compound"
    assert not (prep / "dock").exists(), "a refused compound left files in the dossier"
    return prep


def test_two_compounds_do_not_overwrite_each_other(record: Path):
    """The point of the slug in the run id (step 6b): two compounds in the SAME pocket must keep
    two separate runs, two sets of files and two state files.
    Control that must fail: run_id_of ignoring the slug (verified by hand: with `return
    f"{comp}_{conformer}"` the manifest keeps one run and this test fails)."""
    prep, _ = make_prep()
    ligand_choice.select_ligands(ACC, SEL, reason="r")
    a = asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C", compound="CHEMBL1229517"))
    b = asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C", compound="CHEMBL1336"))
    assert a["ok"] and b["ok"], (a.get("error"), b.get("error"))
    m = manifest(prep)
    runs = next(s for s in m["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]
    assert sorted(runs) == ["032_regenerated_chembl1229517", "032_regenerated_chembl1336"], sorted(runs)
    assert a["replaced_run_computed_at"] is None and b["replaced_run_computed_at"] is None
    files = sorted(p.name for p in (prep / "dock" / "3OG7_032_C").iterdir())
    for slug in ("chembl1229517", "chembl1336"):
        assert f"poses_032_regenerated_{slug}.sdf" in files, files
        assert f"ligand_032_regenerated_{slug}.pdbqt" in files
    # the crystal re-docking of the same pocket keeps the id it had before step 6b existed
    c = asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C"))
    assert c["ok"] and dock.run_id_of("032", "regenerated") == "032_regenerated"
    runs = next(s for s in manifest(prep)["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]
    assert sorted(runs) == ["032_regenerated", "032_regenerated_chembl1229517", "032_regenerated_chembl1336"]
    assert "redocking" in runs["032_regenerated"] and "compound_docking" not in runs["032_regenerated"]
    return prep


def test_compound_run_has_no_rmsd_and_no_criterion(record: Path):
    """A score without a reference pose must never arrive with an RMSD, a 0 or a verdict.
    Control that must fail: `rmsd_angstrom: 0.0` or a `within_criterion` in the block."""
    prep, _ = make_prep()
    ligand_choice.select_ligands(ACC, SEL, reason="r")
    r = asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C", compound="CHEMBL1336"))
    assert r["ok"], r
    cd = r["compound_docking"]
    assert cd["rmsd_angstrom"] is None and "no reference pose" in cd["rmsd_not_applicable"]
    assert "within_criterion" not in cd and "criterion" not in cd
    assert all("RMSD not applicable" in p for p in r["poses"])
    assert cd["top_pose_score_kcal_mol"] == -11.354 and cd["heavy_atoms"] == 32   # sorafenib
    assert cd["ligand_efficiency_kcal_mol_per_heavy_atom"] == round(-11.354 / 32, 4)
    # in the box of 3OG7 032 C sorafenib is NOT box_limited: that box is 26.3 A in its long axis,
    # so the flag is a property of the pocket's box and not of the compound
    assert not cd["box_limited"] and "box_limited_note" not in cd, cd
    assert r["box_fit"]["min_clearance_angstrom"] > 3.0, r["box_fit"]
    run = next(s for s in manifest(prep)["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"][
        "032_regenerated_chembl1336"]
    assert "redocking" not in run, "a compound run must not carry the re-docking block"
    assert run["compound_docking"]["is_preregistered_control"] is False
    assert run["compound_docking"]["rmsd_negative_control"] is None
    assert run["ligand"]["source"].startswith("SMILES recorded by select_ligands")
    assert run["ligand"]["level"] == "T3" and run["ligand"]["inchikey"]
    assert run["ligand"]["canonical_smiles"] != VEM, "the compound docked was not the one asked for"
    assert "No input coordinates" in run["ligand"]["note"], run["ligand"]["note"]
    # T1 says the reference pose exists but is not computed yet; T3 says there is none
    t1 = asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C", compound="CHEMBL1229517"))
    assert "not implemented yet" in t1["compound_docking"]["rmsd_not_applicable"], t1["compound_docking"]
    assert not t1["compound_docking"]["box_limited"]
    # parameters off the pre-registration: the score stops being comparable with the set
    off = asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C", compound="CHEMBL1336", exhaustiveness=8))
    assert "not comparable with the other compounds" in off["compound_docking"]["note"], off["compound_docking"]
    # the same compound in the real box of 8C7Y TXV U (21.8 x 16.27 x 22.0) IS box_limited, and the
    # flag and its note travel with the score, in the reply and in the manifest
    tight = asyncio.run(dock.dock_ligand(ACC, "8C7Y", "TXV U", compound="CHEMBL1336"))
    assert tight["ok"], tight
    assert tight["box_fit"]["min_clearance_angstrom"] == 1.99, tight["box_fit"]
    assert tight["compound_docking"]["box_limited"]
    assert "NOT comparable" in tight["compound_docking"]["box_limited_note"]
    rec = next(s for s in manifest(prep)["structures"] if s["pdb_id"] == "8C7Y")["docking"]["TXV U"]["runs"][
        "TXV_regenerated_chembl1336"]
    assert rec["compound_docking"]["box_limited"] and rec["compound_docking"]["box_limited_note"]
    assert rec["compound_docking"]["box_fit"]["limiting_box_axis_angstrom"] == 22.0
    return prep


def test_dock_status_finds_a_compound_run(record: Path):
    """The conformer filter used endswith on the run id, which stopped matching once the slug was
    appended. Control that must fail: going back to endswith (verified by hand)."""
    prep, _ = make_prep()
    ligand_choice.select_ligands(ACC, SEL, reason="r")
    asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C", compound="CHEMBL1336"))
    asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C"))
    out = asyncio.run(dock.dock_status(ACC, conformer="regenerated"))
    assert out["ok"] and "2 of 2" in out["selection"], out["selection"]
    out = asyncio.run(dock.dock_status(ACC, compound="chembl1336"))
    assert out["ok"] and "1 of 2" in out["selection"], out["selection"]
    out = asyncio.run(dock.dock_status(ACC, compound="chembl9999"))
    assert not out["ok"] and "No docking run matches" in out["error"], out


def main():
    import tempfile
    record = make_fake_vina(Path(tempfile.mkdtemp()))
    test_box_clearance_matches_the_measurement()
    test_refusals_write_nothing()
    test_records_the_selection()
    test_dock_refuses_without_a_selection(record)
    test_two_compounds_do_not_overwrite_each_other(record)
    test_compound_run_has_no_rmsd_and_no_criterion(record)
    test_dock_status_finds_a_compound_run(record)
    print("ALL OK")


if __name__ == "__main__":
    main()
