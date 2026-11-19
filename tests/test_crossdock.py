"""Offline test for cross_dock_rmsd (protein_mcp/files/crossdock.py), the RMSD of a level-T1
compound against its own crystallographic pose in another structure (Phase 3, step 6b).

Same fixtures as test_dock and test_ligand_choice: the real cleaned manifest, the real mmCIFs of
3OG7 and 8C7Y, the real ModelServer SDF of 032 (vemurafenib) and Vina replaced by the stand-in
that returns the real 1.2.7 output.

The measurement this file ties the code to was made on 2026-11-17 with a separate script
(dossiers/_step6b_prep/retrieval/anchor_check.py) and is quoted in ESTADO 7.4.

Every check has a control that must fail if the code stops doing the thing (verified by hand with
the mutations named in the comments). Run: uv run python tests/test_crossdock.py
"""
import asyncio
import json
import sys
from pathlib import Path

from protein_mcp.files import crossdock, dock, ligand_choice

sys.path.insert(0, str(Path(__file__).parent))
from test_choice import ACC  # noqa: E402
from test_dock import REAL_RMSD, SDF, make_fake_vina, make_prep  # noqa: E402
from test_ligand_choice import SORA, VEM  # noqa: E402

# measured 2026-11-17, 3OG7 chain A on 8C7Y chain A, gemmi SupSelect.CaP
SUP_RMSD, SUP_PAIRS = 1.671, 241
ALL_PAIRS = [(1.671, 241), (1.691, 243), (1.932, 248), (2.06, 249)]

SEL = [{"id": "CHEMBL1229517", "name": "VEMURAFENIB", "level": "T1", "smiles": VEM},
       {"id": "CHEMBL1336", "name": "SORAFENIB", "level": "T3", "smiles": SORA}]


def manifest(prep: Path) -> dict:
    return json.loads((prep / "manifest.json").read_text(encoding="utf-8"))


def prepared() -> Path:
    """A prep folder with the selection recorded and VEMURAFENIB docked into 8C7Y TXV U."""
    prep, _ = make_prep()
    out = ligand_choice.select_ligands(ACC, SEL, "the researcher's own words, for the test")
    assert out["ok"], out
    run = asyncio.run(dock.dock_ligand(ACC, "8C7Y", "TXV U", compound="CHEMBL1229517"))
    assert run["ok"], run.get("error")
    return prep


def test_superposition_matches_the_measurement():
    """gemmi on the two real mmCIFs must give the figures of 2026-11-17, and the pair CHOSEN must
    be the lowest of the four, not the first one tried.
    Control that must fail: returning tried[0] instead of the minimum (checked by hand: the
    iteration order gives 3OG7 A -> 8C7Y A first, so the assertion on the full list is what
    catches it), or any of the four figures moving by 0.01."""
    prep = prepared()
    (sup, mob, tgt), tried = crossdock.best_superposition(prep / "raw" / "8C7Y.cif",
                                                          prep / "raw" / "3OG7.cif")
    assert (mob, tgt) == ("A", "A"), (mob, tgt)
    assert round(float(sup.rmsd), 3) == SUP_RMSD and int(sup.count) == SUP_PAIRS, (sup.rmsd, sup.count)
    assert [(t["rmsd_angstrom"], t["pairs"]) for t in tried] == ALL_PAIRS, tried
    assert tried == sorted(tried, key=lambda d: d["rmsd_angstrom"]), tried


def test_the_moved_pose_lands_in_the_box():
    """The reference pose carried into 8C7Y must fall where it was measured to fall: 33 of 33
    heavy atoms in the box, centroid 3.64 A from the box centre.
    Control that must fail: applying the inverse transform, or not applying it at all (then the
    pose stays in the 3OG7 frame, 0 of 33 atoms in the box). Checked by hand both ways."""
    from rdkit import Chem

    prep = prepared()
    (sup, _m, _t), _tried = crossdock.best_superposition(prep / "raw" / "8C7Y.cif",
                                                         prep / "raw" / "3OG7.cif")
    with open(SDF, "rb") as fh:
        crystal = next(Chem.ForwardSDMolSupplier(fh, removeHs=False))
    box = {"center": [0.02, 110.789, 15.175], "size": [21.8, 16.27, 22.0]}
    moved = crossdock.in_box(crossdock.heavy_xyz(crossdock.moved_copy(crystal, sup.transform)), box)
    assert moved["heavy_atoms_in_box"] == "33 of 33", moved
    assert moved["centroid_in_box"] and moved["centroid_to_box_centre_angstrom"] == 3.64, moved
    # the control: the same pose NOT moved is nowhere near this box
    still = crossdock.in_box(crossdock.heavy_xyz(crystal), box)
    assert still["heavy_atoms_in_box"] == "0 of 33" and not still["centroid_in_box"], still


def test_rmsd_reproduces_the_independently_measured_values():
    """With the identity superposition (a structure on itself) the moved reference IS the crystal
    pose, so the RMSDs of the real poses must be the list measured in the prototype of 2026-11-05.
    This is what ties the measurement to dock.py's own function.
    Control that must fail: rdMolAlign.GetBestRMS instead of CalcRMS (it superposes the ligands
    and returns smaller numbers: 0.44 instead of 0.5396 for mode 1). Checked by hand."""
    from rdkit import Chem

    prep = prepared()
    (sup, mob, tgt), _tried = crossdock.best_superposition(prep / "raw" / "3OG7.cif",
                                                           prep / "raw" / "3OG7.cif")
    # any chain on itself is the identity; which pair wins among equals does not matter here
    assert mob == tgt and round(float(sup.rmsd), 3) == 0.0, (mob, tgt, sup.rmsd)
    with open(SDF, "rb") as fh:
        crystal = next(Chem.ForwardSDMolSupplier(fh, removeHs=False))
    moved = crossdock.moved_copy(crystal, sup.transform)
    run = asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C"))
    assert run["ok"], run.get("error")
    rec = next(s for s in manifest(prep)["structures"]
               if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]["032_regenerated"]
    poses_path = prep / next(v["path"] for k, v in rec["files"].items() if k.endswith(".sdf"))
    with open(poses_path, "rb") as fh:
        poses = [m for m in Chem.ForwardSDMolSupplier(fh, removeHs=False)]
    assert crossdock._rmsds(poses, moved) == REAL_RMSD, crossdock._rmsds(poses, moved)


def test_refuses_a_compound_without_a_reference_pose():
    """T3 is not a gap in the data: no crystal of this target holds the compound, so the RMSD does
    not exist. The refusal must say so and must not fall back on another molecule.
    Control that must fail: dropping the level check (then SORAFENIB reaches the InChIKey search
    and is refused for a different reason, so the message is asserted, not just the failure)."""
    prepared()
    out = crossdock.cross_dock_rmsd(ACC, "CHEMBL1336")
    assert not out["ok"] and out.get("level") == "T3", out
    assert "not applicable" in out["error"] and "T1" in out["error"], out["error"]
    # and a compound that was never recorded is refused by the selection gate, not by this one
    out = crossdock.cross_dock_rmsd(ACC, "CHEMBL999999")
    assert not out["ok"] and "CHEMBL999999" in out["error"], out["error"]


def test_the_reference_is_found_by_inchikey_and_not_by_name():
    """ESTADO 7.1: an identity is a key, not a name. The reference must still be found when the
    recorded name is wrong, and must NOT be found when the key does not match.
    Control that must fail: matching on comp_id or on `name` (then the wrong-name case finds
    nothing and the wrong-key case still finds 3OG7 032 C). Checked by hand, both directions."""
    prep = prepared()
    path = prep / "manifest.json"
    m = manifest(prep)
    entry = next(e for e in m["ligand_choice"]["compounds"] if e["id"] == "CHEMBL1229517")
    entry["name"] = "NOT THE NAME OF THIS COMPOUND"
    path.write_text(json.dumps(m, ensure_ascii=False), encoding="utf-8")
    out = crossdock.cross_dock_rmsd(ACC, "CHEMBL1229517", "8C7Y", "TXV U")
    assert out["ok"], out.get("error")
    assert out["cross_docking"]["reference_pose"]["from"] == "3OG7 032 C", out["cross_docking"]

    m = manifest(prep)
    entry = next(e for e in m["ligand_choice"]["compounds"] if e["id"] == "CHEMBL1229517")
    entry["inchikey"] = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    path.write_text(json.dumps(m, ensure_ascii=False), encoding="utf-8")
    out = crossdock.cross_dock_rmsd(ACC, "CHEMBL1229517", "8C7Y", "TXV U")
    assert not out["ok"], f"a wrong InChIKey still found a reference pose: {out.get('cross_docking')}"
    assert "crystallographic instance" in out["error"].lower(), out["error"]
    assert "not applicable" in out["error"], out["error"]


def test_the_number_never_travels_without_its_bar():
    """The RMSD, the shift bar, the superposition it came from and whether the criterion decides
    must all be in the result and in the chat view. The bar is the point of the whole module.
    Control that must fail: removing shift_angstrom from the block, or reporting within_criterion
    without decides (both checked by hand)."""
    prepared()
    out = crossdock.cross_dock_rmsd(ACC, "CHEMBL1229517", "8C7Y", "TXV U")
    assert out["ok"], out.get("error")
    b = out["cross_docking"]
    assert b["superposition"]["rmsd_angstrom"] == SUP_RMSD, b["superposition"]
    assert b["rmsd_negative_control"]["passed"], b["rmsd_negative_control"]
    assert len(b["rmsd_angstrom_by_pose"]) == 9, b["rmsd_angstrom_by_pose"]
    bar = b["reference_pose_uncertainty"]["shift_angstrom"]
    assert isinstance(bar, float) and bar > 0, b["reference_pose_uncertainty"]
    lo, hi = b["within_criterion_with_bar"]["range_angstrom"]
    assert abs((hi - lo) - 2 * bar) < 0.01, (lo, hi, bar)
    assert b["within_criterion_with_bar"]["decides"] in (True, False)
    small = crossdock.compact(out)
    for key in ("top_pose_rmsd_angstrom", "shift_bar_angstrom", "superposition_rmsd_angstrom",
                "criterion_angstrom", "within_criterion", "decides", "negative_control"):
        assert key in small, key
    # the criterion note travels with it: a number read alone is the failure this prevents
    assert any("edge of the noise" in n.lower() for n in small["notes"]), small["notes"]


def test_it_is_written_into_the_run_it_measures():
    """The RMSD belongs to the docking run, not to a separate place where it could be read without
    the run's parameters. The per-pose values must be filled in too.
    Control that must fail: writing the block at the pocket level instead of inside the run."""
    prep = prepared()
    out = crossdock.cross_dock_rmsd(ACC, "CHEMBL1229517", "8C7Y", "TXV U")
    assert out["ok"] and out["recorded_in_manifest"], out
    s = next(x for x in manifest(prep)["structures"] if x["pdb_id"] == "8C7Y")
    run = s["docking"]["TXV U"]["runs"]["TXV_regenerated_chembl1229517"]
    assert run["cross_docking"]["top_pose_rmsd_angstrom"] == out["cross_docking"]["top_pose_rmsd_angstrom"]
    assert run["poses"][0]["rmsd_to_reference_pose_angstrom"] is not None, run["poses"][0]
    assert run["poses"][0]["rmsd_reference"] == "3OG7 032 C", run["poses"][0]
    # the field of the re-docking kind stays null: this run never had a pose of its own to compare
    assert run["poses"][0]["rmsd_to_crystal_pose_angstrom"] is None, run["poses"][0]


def main():
    import tempfile
    make_fake_vina(Path(tempfile.mkdtemp()))
    test_superposition_matches_the_measurement()
    test_the_moved_pose_lands_in_the_box()
    test_rmsd_reproduces_the_independently_measured_values()
    test_refuses_a_compound_without_a_reference_pose()
    test_the_reference_is_found_by_inchikey_and_not_by_name()
    test_the_number_never_travels_without_its_bar()
    test_it_is_written_into_the_run_it_measures()
    print("ALL OK")


if __name__ == "__main__":
    main()
