"""Offline test for open_in_pymol (protein_mcp/files/viewer.py): scene from the manifest, ASCII
.pml, PyMOL lookup and launch check. The manifest is the real BRAF V600E one (2026-10-26), cut to
the fields the viewer reads: 3OG7 (vemurafenib on chain A, GLU 600 only in chain B), 5JSM (6NB
pocket on chains A and B), 7ZR0 (BRAF is auth K; ATP sits on HSP90 only), plus made-up cases
for IDs PyMOL cannot select. The PyMOL defects it guards against (non-ASCII paths, stdout, no
Qt) were seen on this machine: ESTADO_PROYECTO.md section 8.2.
Run: uv run python tests/test_viewer.py
"""
import asyncio
import json
import os
import sys
import tempfile
from pathlib import Path

from protein_mcp import dossier
from protein_mcp.files import viewer

viewer.LAUNCH_CHECK_SECONDS = 1.5


def res(chain, seq, name="ALA"):
    return {"auth_asym_id": chain, "resname": name, "auth_seq_id": seq, "uniprot": None}


def structure(pdb, receptor, chains, instances, pockets, ok=True):
    return {"pdb_id": pdb, "nonpolymer_instances": instances,
            "prepared": {"ok": ok, "receptor": {"path": receptor}, "metrics": {"chains": chains},
                         "pockets": pockets} if ok else {"ok": False, "error": "Receptor is empty"}}


BOX_3OG7 = {"center": [2.643, -2.28, -19.403], "size": [26.3, 15.13, 16.4]}
MANIFEST = {
    "uniprot_id": "P15056", "created_at": "2026-10-26T23:32:07Z",
    "structures": [
        structure("3OG7", "clean/3OG7_receptor.pdb",
                  [{"auth_asym_id": "A", "variant_check": {"observed": None, "status": "not modelled"}},
                   {"auth_asym_id": "B", "variant_check": {"observed": "GLU 600", "status": "variant"}}],
                  [{"comp_id": "032", "label_asym_id": "C", "auth_asym_id": "A", "name": "Vemurafenib",
                    "neighbor_target_chains": {"A": {"auth_asym_id": "A", "closest_angstrom": 2.59}},
                    "sdf": {"ok": True, "path": "ligands/3OG7_032_C.sdf"}}],
                  [{"comp_id": "032", "label_asym_id": "C", "auth_asym_id": "A", "box": BOX_3OG7,
                    "residues": [res("A", "463", "ILE"), res("A", "483", "LYS"), res("A", "594", "ASP")]}]),
        structure("5JSM", "clean/5JSM_receptor.pdb",
                  [{"auth_asym_id": "A", "variant_check": {"observed": "GLU 600", "status": "variant"}},
                   {"auth_asym_id": "B", "variant_check": {"observed": "GLU 600", "status": "variant"}}],
                  [{"comp_id": "6NB", "label_asym_id": "J", "auth_asym_id": "A",
                    "neighbor_target_chains": {"A": {}, "B": {}}, "sdf": {"ok": True, "path": "ligands/5JSM_6NB_J.sdf"}}],
                  [{"comp_id": "6NB", "label_asym_id": "J", "auth_asym_id": "A", "spans_target_chains": True,
                    "box": {"center": [1, 2, 3], "size": [40, 20, 20]},
                    "residues": [res("A", "461"), res("A", "463"), res("B", "461"), res("B", "463")]}]),
        structure("7ZR0", "clean/7ZR0_receptor.pdb",
                  [{"auth_asym_id": "K", "variant_check": {"observed": None, "status": "not modelled"}}],
                  [{"comp_id": "ATP", "label_asym_id": "E", "auth_asym_id": "A", "neighbor_target_chains": {},
                    "neighbor_other_chains": {"A": {"auth_asym_id": "A", "closest_angstrom": 2.33}},
                    "sdf": {"ok": True, "path": "ligands/7ZR0_ATP_E.sdf"}}],
                  []),
        # made-up: negative number, insertion code, a chain name PyMOL cannot select, a pocket error,
        # a path outside prep/ and a ligand without SDF
        structure("9XYZ", "clean/9XYZ_receptor.pdb",
                  [{"auth_asym_id": "A", "variant_check": {"observed": "GLU 100A", "status": "variant"}}],
                  [{"comp_id": "LIG", "label_asym_id": "B", "sdf": {"ok": True, "path": "../../outside.sdf"}},
                   {"comp_id": "ZZZ", "label_asym_id": "C", "sdf": None},
                   {"comp_id": "EDO", "label_asym_id": "E", "category": "additive", "sdf": None},
                   # shape of the real failure of 6V34 QOP C (2026-10-26, before the retries)
                   {"comp_id": "QOP", "label_asym_id": "F", "category": "ligand",
                    "sdf": {"ok": False, "url": "https://models.rcsb.org/v1/6V34/ligand?label_asym_id=C&encoding=sdf",
                            "error": "RemoteProtocolError: Server disconnected without sending a response."}}],
                  [{"comp_id": "LIG", "label_asym_id": "B", "box": {"center": [0, 0], "size": [1, 1, 1]},
                    "residues": [res("A", "-5"), res("A", "100A"), res("A B", "7")]},
                   {"comp_id": "LIG", "label_asym_id": "D", "error": "gemmi failed"}]),
        structure("1BAD", "clean/1BAD_receptor.pdb", [], [], [], ok=False),
        structure("2MIS", "clean/2MIS_receptor.pdb", [], [], []),          # receptor file not on disk
    ],
    "alphafold": {"entry_id": "AF-P15056-F1", "file": {"ok": True, "path": "raw/AF-P15056-F1-model_v6.cif"}},
}
# 7ZR0 partner chains as in the real manifest: HSP90-beta is A and B, CDC37 is C (the title lists
# "HSP90-CDC37-BRAF"; Claude once read it as A = HSP90, B = CDC37 in the Desktop test)
MANIFEST["structures"][2]["other_polymer_entities"] = [
    {"chains": [{"label_asym_id": "A", "auth_asym_id": "A"}, {"label_asym_id": "B", "auth_asym_id": "B"}],
     "description": "Heat shock protein HSP 90-beta", "uniprot_ids": ["P08238"]},
    {"chains": [{"label_asym_id": "C", "auth_asym_id": "C"}], "description": "Hsp90 co-chaperone Cdc37",
     "uniprot_ids": ["Q16543"]},
]


def make_prep(root: Path) -> Path:
    """dossiers/P15056/prep under a folder with a space and non-ASCII letters (a Windows user name
    with accents): PyMOL cannot `load` such a path, so the .pml must escape it."""
    dossier.DOSSIER_DIR = root / "dossiers ñá"
    prep = dossier.DOSSIER_DIR / "P15056" / "prep"
    for rel in ["clean/3OG7_receptor.pdb", "clean/5JSM_receptor.pdb", "clean/7ZR0_receptor.pdb",
                "clean/9XYZ_receptor.pdb", "ligands/3OG7_032_C.sdf", "ligands/5JSM_6NB_J.sdf",
                "ligands/7ZR0_ATP_E.sdf", "raw/AF-P15056-F1-model_v6.cif"]:
        (prep / rel).parent.mkdir(parents=True, exist_ok=True)
        (prep / rel).write_text("stub\n", encoding="utf-8")
    (dossier.DOSSIER_DIR / "P15056" / "outside.sdf").write_text("stub\n", encoding="utf-8")
    (prep / "manifest.json").write_text(json.dumps(MANIFEST), encoding="utf-8")
    return prep


def test_scene_and_pml(prep: Path):
    r = asyncio.run(viewer.open_in_pymol("P15056", start=False))
    assert r["ok"], r
    assert r["selection"] == "4 of 6 structures in the manifest + the AlphaFold model"
    assert {s["pdb_id"]: s["reason"] for s in r["skipped"]} == {
        "1BAD": "Receptor is empty", "2MIS": "receptor file missing from prep/"}
    by = {s["pdb_id"]: s for s in r["structures"]}
    assert [s["enabled"] for s in r["structures"]] == [True, False, False, False]
    assert by["3OG7"]["variant_position"] == {"selection": "variant_pos_3OG7", "chains": ["B"], "by_chain": [
        {"auth_asym_id": "A", "observed": None, "status": "not modelled"},
        {"auth_asym_id": "B", "observed": "GLU 600", "status": "variant"}]}
    assert by["7ZR0"]["variant_position"] is None
    assert by["7ZR0"]["ligands"][0]["on_target_chain"] is False and by["3OG7"]["ligands"][0]["on_target_chain"]
    assert by["3OG7"]["ligands"][0]["has_pocket"] and not by["7ZR0"]["ligands"][0]["has_pocket"]
    assert by["7ZR0"]["ligands"][0]["partner_chains"] == ["A: Heat shock protein HSP 90-beta"]
    assert by["3OG7"]["ligands"][0]["partner_chains"] == []
    assert any(n.startswith("Loaded but not on a target chain")
               and "7ZR0 ATP E (on A: Heat shock protein HSP 90-beta)" in n for n in r["notes"])
    assert [l["comp_id"] for l in by["9XYZ"]["ligands"]] == []          # outside prep/ and no SDF: not loaded
    assert by["9XYZ"]["ligands_not_loaded"] == [               # the additive EDO is left out on purpose
        {"comp_id": "LIG", "label_asym_id": "B", "reason": "SDF file not found in prep/"},
        {"comp_id": "ZZZ", "label_asym_id": "C", "reason": "no SDF in the manifest"},
        {"comp_id": "QOP", "label_asym_id": "F", "reason": "SDF download failed"}]
    assert by["9XYZ"]["pockets_not_shown"] == [
        {"comp_id": "LIG", "label_asym_id": "D", "reason": "pocket error in the manifest: gemmi failed"}]
    p9 = by["9XYZ"]["pockets"]
    assert len(p9) == 1 and p9[0]["box_object"] is None and p9[0]["selection"] == "pocket_9XYZ_LIG_B"
    assert p9[0]["n_residues"] == 3 and p9[0]["n_residues_selected"] == 2 and p9[0]["residues_not_selectable"] == ["A B:7"]
    assert by["3OG7"]["pockets"][0]["n_residues_selected"] == 3 and "residues_not_selectable" not in by["3OG7"]["pockets"][0]
    assert r["alphafold"]["aligned_to"] == "rec_3OG7" and r["pymol"]["launched"] is False

    pml = Path(r["pml_path"])
    raw = pml.read_bytes()
    raw.decode("ascii")                                    # pure ASCII, no BOM
    text = raw.decode("ascii")
    assert "dossiers \\xf1\\xe1/P15056/prep/clean/3OG7_receptor.pdb" in text
    assert "\nload " not in text and text.startswith("# ")
    lines = text.splitlines()
    assert lines[lines.index("python") + 1] == "from pymol import cmd" and lines[-1] == "python end"
    assert "_do('select', 'pocket_3OG7_032_C', 'rec_3OG7 and ((chain A and resi 463+483+594))')" in text
    assert "_do('select', 'pocket_5JSM_6NB_J', 'rec_5JSM and ((chain A and resi 461+463) or (chain B and resi 461+463))')" in text
    assert "'rec_9XYZ and ((chain A and resi \\\\-5+100A))'" in text      # "A B" dropped, -5 escaped
    assert "_do('select', 'variant_pos_5JSM', 'rec_5JSM and ((chain A and resi 600) or (chain B and resi 600))')" in text
    assert "_do('group', 's_3OG7', 'rec_3OG7 lig_3OG7_032_C box_3OG7_032_C')" in text
    assert "_do('align', 'af_P15056', 'rec_3OG7')" in text
    assert "_do('disable', 's_5JSM')" in text and "_do('disable', 's_3OG7')" not in text
    # the scene ends with deselect; then the chat-control receiver (step 4b), then the loaded marker
    ops_end = max(i for i, l in enumerate(lines) if l.startswith("_do("))
    assert "_do('orient', 'pocket_3OG7_032_C')" in text and lines[ops_end] == "_do('deselect')"
    assert "    _pmcp_receiver(" in "\n".join(lines[ops_end:])
    assert lines[-2] == "print(' protein-mcp: session loaded')"
    assert "pocket_9XYZ_LIG_D" not in text and "outside.sdf" not in text
    # the 3OG7 box: 12 edges, corners at center +- size/2
    box_line = next(l for l in lines if l.startswith("_do('load_cgo'"))
    cgo = eval(box_line[len("_do("):-1])[1]                    # noqa: S307 (our own output)
    assert cgo.count(viewer._CGO_VERTEX) == 24
    xs = {round(v, 3) for v in cgo[7::4][:24]}
    assert xs == {round(2.643 - 13.15, 3), round(2.643 + 13.15, 3)}
    print("scene and pml OK")


def test_subset(prep: Path):
    r = asyncio.run(viewer.open_in_pymol("P15056", pdb_ids=["5jsm", "1ABC"], include_alphafold=False, start=False))
    assert r["ok"] and [s["pdb_id"] for s in r["structures"]] == ["5JSM"] and r["structures"][0]["enabled"]
    assert r["pdb_ids_not_in_manifest"] == ["1ABC"] and r["alphafold"] is None
    assert "af_P15056" not in Path(r["pml_path"]).read_text(encoding="ascii")
    none = asyncio.run(viewer.open_in_pymol("P15056", pdb_ids=["1BAD"], include_alphafold=False, start=False))
    assert not none["ok"] and none["skipped"][0]["pdb_id"] == "1BAD"
    missing = asyncio.run(viewer.open_in_pymol("Q99999", start=False))
    assert not missing["ok"] and "run prepare_structures first" in missing["error"]
    print("subset and errors OK")


def test_render_rejects():
    try:
        viewer.render_pml([("do", ["python import os"])], [])
    except ValueError as exc:
        assert "'do'" in str(exc)
    else:
        raise AssertionError("cmd.do must be rejected")


def test_find_and_launch(prep: Path):
    old = os.environ.pop(viewer.PYMOL_ENV, None)
    path_before = os.environ.get("PATH", "")
    try:
        os.environ["PATH"] = ""
        assert viewer.find_pymol()["found"] is False
        os.environ[viewer.PYMOL_ENV] = str(prep / "no_such_pymol.exe")
        f = viewer.find_pymol()
        assert not f["found"] and "not a file" in f["error"]
        os.environ[viewer.PYMOL_ENV] = f'"{sys.executable}"'        # quoted, as pasted from Explorer
        f = viewer.find_pymol()
        windowless = Path(sys.executable).with_name("pythonw.exe")     # used when it exists (step 4b)
        exe = str(windowless) if Path(sys.executable).name.lower() == "python.exe" and windowless.is_file() else sys.executable
        assert f["found"] and f["command"] == [exe, "-m", "pymol"], f
    finally:
        os.environ["PATH"] = path_before
        os.environ.pop(viewer.PYMOL_ENV, None)
        if old is not None:
            os.environ[viewer.PYMOL_ENV] = old

    pml = prep / "pymol" / "session.pml"
    logf = prep / "pymol" / "test_launch.log"
    # the real output of PyMOL without Qt (2026-10-26), then its exit
    no_qt = [sys.executable, "-c", "import sys; print('Qt not available (pymol.Qt), using GLUT/Tk interface'); "
             "print('NotImplementedError: compile with --glut'); sys.exit(1)"]
    r = asyncio.run(viewer.launch(no_qt, pml, logf))
    assert r["launched"] is False and r["exit_code"] == 1 and "PyQt5" in r["hint"]
    assert any("compile with --glut" in line for line in r["log_tail"])
    # still running: the session messages already in the log are read; PyMOL's echo of the source is not
    stays = [sys.executable, "-c", "import time; "
             "print(\"    6:        print(' protein-mcp: cmd.%s%r failed: %s' % (fn, args[:2], exc))\"); "
             "print(\" protein-mcp: cmd.load('C:/x/7ZR0_receptor.pdb', 'rec_7ZR0') failed:  Error: failed to open file\"); "
             "print(' protein-mcp: session loaded'); time.sleep(30)"]
    r = asyncio.run(viewer.launch(stays, pml, logf))
    assert r["launched"] is True and r["pid"] and r["session_loaded"] is True
    assert r["session_errors"] == ["protein-mcp: cmd.load('C:/x/7ZR0_receptor.pdb', 'rec_7ZR0') failed:  Error: failed to open file"]
    os.kill(r["pid"], 9)
    slow = [sys.executable, "-c", "import time; time.sleep(30)"]
    r = asyncio.run(viewer.launch(slow, pml, logf))
    assert r["launched"] is True and r["session_loaded"] is False and "later errors" in r["note"]
    os.kill(r["pid"], 9)
    bad = asyncio.run(viewer.launch([str(prep / "no_such.exe")], pml, logf))
    assert bad["launched"] is False and "Could not start PyMOL" in bad["error"]
    print("find and launch OK")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        prep = make_prep(Path(tmp))
        test_scene_and_pml(prep)
        test_subset(prep)
        test_render_rejects()
        test_find_and_launch(prep)
    print("ALL OK")


main()
