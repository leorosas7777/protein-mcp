"""Offline test for dock_ligand (protein_mcp/files/dock.py, Phase 3 step 2).
The manifest is the one of test_choice (real 3OG7, 5JSM, 8C7Y mmCIFs cleaned with clean.py, real
RCSB chain replies), with the real mmCIFs in raw/ and the real ModelServer SDF of 3OG7 032 C.
Meeko really prepares the receptor (about 7 s). Vina is replaced by a stand-in (.cmd) that returns
the REAL Vina 1.2.7 output of 2026-11-05 for 3OG7 032 C, regenerated conformer, seed 42,
exhaustiveness 32 (tests/fixtures/dock/vina_3OG7_032_C_regenerated.*; in the .log only the two
local paths were replaced by <work>, the poses file is untouched). Meeko messages are the
real outputs for 5JSM 6NB J (refused) and 8C7Y TXV U (3 residues deleted).
Run: uv run python tests/test_dock.py
"""
import asyncio
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from protein_mcp import dossier
from protein_mcp.files import choice, dock

sys.path.insert(0, str(Path(__file__).parent))
from test_choice import ACC, FIX, build_manifest  # noqa: E402

DFIX = Path(__file__).parent / "fixtures" / "dock"
SDF = FIX / "modelserver_3OG7_C.sdf"
# CalcRMS of the real poses against the real crystal SDF (prototype run, 2026-11-05)
REAL_RMSD = [0.5396, 2.8565, 6.0683, 8.3492, 8.7793, 2.104, 6.6953, 9.3754, 8.2993]
REAL_SCORES = [-11.354, -10.965, -10.936, -10.614, -10.612, -10.52, -10.512, -10.438, -10.345]

FAKE_VINA = r'''
import json, os, shutil, sys
args = sys.argv[1:]
if "--version" in args:
    print("AutoDock Vina v1.2.7")
    sys.exit(0)
with open(os.environ["FAKE_VINA_RECORD"], "a", encoding="utf-8") as fh:
    fh.write(json.dumps(args) + "\n")
touch = os.environ.get("FAKE_VINA_TOUCH_MANIFEST")
if touch:                      # prepare_structures rebuilding the manifest during docking
    m = json.load(open(touch, encoding="utf-8")); m["created_at"] = "2026-11-05T00:00:00Z"
    json.dump(m, open(touch, "w", encoding="utf-8"))
if os.environ.get("FAKE_VINA_FAIL"):
    print("ERROR: something broke"); sys.exit(1)
nap = os.environ.get("FAKE_VINA_SLEEP")        # stay alive, so a background run can be seen and stopped
if nap:
    import time
    time.sleep(float(nap))
shutil.copyfile(os.environ["FAKE_VINA_POSES"], args[args.index("--out") + 1])
print(open(os.environ["FAKE_VINA_LOG"], encoding="utf-8").read())
'''


def make_fake_vina(tmp: Path) -> Path:
    (tmp / "fake_vina.py").write_text(FAKE_VINA, encoding="utf-8")
    cmd = tmp / "vina_fake.cmd"
    cmd.write_text(f'@"{sys.executable}" -I "{tmp / "fake_vina.py"}" %*\r\n', encoding="ascii")
    os.environ.update(PROTEIN_MCP_VINA=str(cmd), FAKE_VINA_RECORD=str(tmp / "calls.jsonl"),
                      FAKE_VINA_POSES=str(DFIX / "vina_3OG7_032_C_regenerated.poses.pdbqt"),
                      FAKE_VINA_LOG=str(DFIX / "vina_3OG7_032_C_regenerated.log"))
    return tmp / "calls.jsonl"


def calls(record: Path) -> list:
    return [json.loads(l) for l in record.read_text(encoding="utf-8").splitlines()] if record.exists() else []


def make_prep(folder_name="dossiers") -> tuple[Path, dict]:
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp()) / folder_name
    prep = dossier.DOSSIER_DIR / ACC / "prep"
    m = build_manifest(prep / "clean")
    (prep / "raw").mkdir(parents=True)
    (prep / "ligands").mkdir()
    shutil.copy(SDF, prep / "ligands" / "3OG7_032_C.sdf")
    for s in m["structures"]:
        if s["file"].get("ok"):
            data = gzip.decompress((FIX / f"{s['pdb_id']}.cif.gz").read_bytes())
            (prep / s["file"]["path"]).write_bytes(data)
            s["file"]["sha256"] = hashlib.sha256(data).hexdigest()
        for i in s.get("nonpolymer_instances") or []:
            if (i.get("sdf") or {}).get("path") == "ligands/3OG7_032_C.sdf":
                i["sdf"]["sha256"] = hashlib.sha256(SDF.read_bytes()).hexdigest()
    choice.write_manifest(prep / "manifest.json", m)
    return prep, m


def run(*args, **kw):
    return asyncio.run(dock.dock_ligand(ACC, *args, **kw))


def test_parse_meeko():
    out = (DFIX / "meeko_5jsm_6NB_J_refused.stdout.txt").read_text(encoding="utf-8")
    err = (DFIX / "meeko_5jsm_6NB_J_refused.stderr.txt").read_text(encoding="utf-8")
    p = dock.parse_meeko(out, err)
    assert [r["residue"] for r in p["refused"]] == ["A:468", "A:475", "A:487", "A:488", "A:575", "A:600", "B:488"], p["refused"]
    assert next(r for r in p["refused"] if r["residue"] == "A:600")["distance_to_box_edge_angstrom"] == 2.456
    assert len(p["deleted"]) == 34 and p["deleted"][0] == {"residue": "A:662", "outside_box_angstrom": 15.064}
    # 34 deleted + 7 refused = the 41 residues without a template of step 1 (ESTADO 8.3)
    assert len(p["residues_without_template"]) == 41, len(p["residues_without_template"])
    p = dock.parse_meeko((DFIX / "meeko_8c7y_TXV_U_ok.stdout.txt").read_text(encoding="utf-8"),
                         (DFIX / "meeko_8c7y_TXV_U_ok.stderr.txt").read_text(encoding="utf-8"))
    assert p["refused"] == [] and p["residues_without_template"] == ["A:607", "A:608", "A:609"]
    assert [(d["residue"], d["outside_box_angstrom"]) for d in p["deleted"]] == [
        ("A:607", 17.487), ("A:608", 14.981), ("A:609", 11.859)]
    assert dock.unaccounted_warning(p) is None
    # a changed message format: a residue without template that is not in the deleted list is reported
    p["deleted"] = p["deleted"][:2]
    assert "A:609" in dock.unaccounted_warning(p)


def test_poses_and_rmsd():
    crystal = dock.read_sdf(SDF)
    text = (DFIX / "vina_3OG7_032_C_regenerated.poses.pdbqt").read_text(encoding="utf-8")
    poses, scores = dock.read_poses(text)
    assert poses.GetNumConformers() == 9 and scores == REAL_SCORES, scores
    assert dock.rmsd_to_crystal(poses, crystal) == REAL_RMSD
    # the $$$$ trap of step 1: the whole SDF must read back with all 9 poses
    tmp = Path(tempfile.mkdtemp()) / "carpeta ñá" / "poses.sdf"
    tmp.parent.mkdir()
    tmp.write_text(dock.poses_sdf(poses, scores), encoding="utf-8")
    assert dock.count_sdf_records(tmp) == 9
    # a truncated file (one MODEL lost) is refused, not silently shortened
    cut = text[: text.rindex("MODEL")]
    try:
        dock.read_poses(cut + "MODEL 9\n")
        raise AssertionError("pose count mismatch not detected")
    except dock.DockError:
        pass
    nc = dock.negative_control(crystal)
    assert nc["passed"] and nc["crystal_moved_5A_angstrom"] == 5.0 and nc["crystal_vs_itself_angstrom"] == 0.0
    # the control goes through the pose-measuring function: a superposing RMSD makes it fail (reviewer H3)
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    def best(poses, ref):
        r, p = Chem.RemoveHs(ref), Chem.RemoveHs(poses)
        return [round(rdMolAlign.GetBestRMS(p, r, prbId=c.GetId()), 4) for c in p.GetConformers()]
    try:
        dock.negative_control(crystal, best)
        raise AssertionError("negative control passed with GetBestRMS")
    except dock.DockError as exc:
        assert exc.extra["rmsd_negative_control"]["crystal_moved_5A_angstrom"] == 0.0


def test_guard():
    crystal = dock.read_sdf(SDF)
    box = {"center": [2.643, -2.28, -19.403], "size": [26.3, 15.13, 16.4]}   # manifest box of 3OG7 032 C
    near = np.array([[2.643, -2.28, -15.0]])
    mol, rec = dock.build_ligand(crystal, "crystal", box, 42)
    g = dock.guard(mol, box, near, "crystal")
    assert g["ligand_atoms_in_box"] == f"{mol.GetNumAtoms()} of {mol.GetNumAtoms()}" and "by construction" in g["note"]
    mol, rec = dock.build_ligand(crystal, "regenerated", box, 42)
    assert rec["translation_angstrom"] and rec["smiles"] == "CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(-c4ccc(Cl)cc4)cc23)c1F"
    pos = mol.GetConformer().GetPositions()
    assert np.allclose(pos.mean(axis=0), box["center"], atol=1e-3)
    dock.guard(mol, box, near, "regenerated")
    # another frame (step 1: 19.7 Å away, Vina docked it anyway): refused here
    far = {"center": [22.3, -2.28, -19.403], "size": box["size"]}
    for args, text in (((mol, far, near, "regenerated"), "Ligand outside the box"),
                       ((mol, box, np.zeros((0, 3)), "regenerated"), "No receptor atom")):
        try:
            dock.guard(*args)
            raise AssertionError("guard passed")
        except dock.DockError as exc:
            assert text in str(exc), exc


def test_ascii_root():
    real = tempfile.gettempdir
    try:
        tempfile.gettempdir = lambda: r"C:\Users\Jesús\AppData\Local\Temp"
        root = dock.ascii_work_root(str(Path(tempfile.mkdtemp(dir=real())) / "vina.exe"))
        assert str(root).isascii() and root.name == "tmp"
        try:
            dock.ascii_work_root(r"C:\Programas ñ\vina.exe")
            raise AssertionError("non-ASCII accepted")
        except dock.DockError as exc:
            assert "ASCII" in str(exc)
    finally:
        tempfile.gettempdir = real


def test_refusals(record: Path):
    prep, m = make_prep()
    before = len(calls(record))
    cases = [
        ({"structure": "5JSM", "pocket": "6NB J"}, "large_box and spans_target_chains"),
        ({"structure": "8C7Y", "pocket": "TXV"}, "matches 2 pockets"),
        ({"structure": "3OG7", "pocket": ""}, "pocket is required"),
        ({}, "no user_choice"),
        ({"structure": "AlphaFold"}, "pocket is required"),
        ({"structure": "3OG7", "pocket": "032 C", "conformer": "smiles"}, "conformer must be"),
        ({"structure": "3OG7", "pocket": "032 C", "exhaustiveness": 0}, "exhaustiveness"),
        ({"structure": "3OG7", "pocket": "032 C", "seed": True}, "seed"),
        ({"structure": "8C7Y", "pocket": "TXV X"}, "no crystal ligand SDF"),   # real 6V34-shape failure
    ]
    for kw, text in cases:
        r = run(**kw)
        assert not r["ok"] and text in r["error"], (kw, r)
    assert "options" in run(structure="8C7Y", pocket="TXV")
    # user_choice without pocket
    m2 = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    m2["user_choice"] = {"structure": "3OG7", "chains": [{"auth_asym_id": "A"}], "pocket": None}
    choice.write_manifest(prep / "manifest.json", m2)
    assert "has no pocket" in run()["error"]
    # receptor that kept a nucleotide/cofactor (another instance of the entry): refused naming it, never deleted
    m3 = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    s = next(x for x in m3["structures"] if x["pdb_id"] == "3OG7")
    # 3OG7 has no such instance: an invented GDP instance in label Z, as clean.py would keep it (KRAS case)
    other = {"label_asym_id": "Z", "comp_id": "GDP", "category": "nucleotide"}
    s["nonpolymer_instances"].append(other)
    s["prepared"]["receptor"]["kept_nonpolymer_label_asym_ids"] = [other["label_asym_id"]]
    choice.write_manifest(prep / "manifest.json", m3)
    r = run(structure="3OG7", pocket="032 C")
    assert f"keeps {other['comp_id']} ({other['label_asym_id']})" in r["error"], r
    # receptor written as mmCIF (chain names that do not fit in PDB): only --read_pdb is validated
    s["prepared"]["receptor"]["kept_nonpolymer_label_asym_ids"] = []
    s["prepared"]["receptor"]["format"] = "mmcif"
    choice.write_manifest(prep / "manifest.json", m3)
    r = run(structure="3OG7", pocket="032 C")
    assert "is mmcif, not PDB" in r["error"], r
    choice.write_manifest(prep / "manifest.json", m)
    # changed files
    (prep / "ligands" / "3OG7_032_C.sdf").write_text(SDF.read_text() + "\n", encoding="utf-8")
    assert "changed since the manifest" in run(structure="3OG7", pocket="032 C")["error"]
    shutil.copy(SDF, prep / "ligands" / "3OG7_032_C.sdf")
    vina = os.environ.pop("PROTEIN_MCP_VINA")
    try:
        assert "PROTEIN_MCP_VINA" in run(structure="3OG7", pocket="032 C")["error"]
        os.environ["PROTEIN_MCP_VINA"] = r"C:\no\such\vina.exe"
        assert "not a file" in run(structure="3OG7", pocket="032 C")["error"]
    finally:
        os.environ["PROTEIN_MCP_VINA"] = vina
    assert len(calls(record)) == before, "Vina was started for a refused request"
    assert not (prep / "dock").exists()


def test_full_run(record: Path, folder="dossiers"):
    prep, m = make_prep(folder)
    progress = []

    async def report(elapsed):
        progress.append(elapsed)
    r = asyncio.run(dock.dock_ligand(ACC, "3OG7", "032 C", progress=report))
    assert r["ok"], r
    assert r["redocking"]["top_pose_rmsd_angstrom"] == 0.5396 and r["redocking"]["within_criterion"]
    assert r["redocking"]["is_preregistered_control"] and r["poses_returned"].startswith("9 of num_modes 9")
    assert r["poses"][0] == "mode 1: -11.35 kcal/mol (calculated), RMSD to crystal pose 0.54 Å"
    assert r["receptor_preparation"]["deleted_residues"] == []      # 3OG7 passes Meeko at once (ESTADO 8.3)
    assert r["ligand"]["rotatable_torsions"] == 7 and r["ligand"]["heavy_atoms"] == 33
    args = calls(record)[-1]
    for flag, value in (("--seed", "42"), ("--exhaustiveness", "32"), ("--num_modes", "9"), ("--energy_range", "3"),
                        ("--min_rmsd", "1"), ("--scoring", "vina"), ("--center_x", "2.643"), ("--size_y", "15.13")):
        assert args[args.index(flag) + 1] == value, (flag, args)
    d = prep / "dock" / "3OG7_032_C"
    assert sorted(p.name for p in d.iterdir()) == sorted([
        "receptor.pdbqt", "receptor_meeko.log", "ligand_032_regenerated.pdbqt", "poses_032_regenerated.pdbqt",
        "poses_032_regenerated.sdf", "vina_032_regenerated.log"])
    assert dock.count_sdf_records(d / "poses_032_regenerated.sdf") == 9
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    s = next(x for x in man["structures"] if x["pdb_id"] == "3OG7")
    run_rec = s["docking"]["032 C"]["runs"]["032_regenerated"]
    assert run_rec["values"] == "calculated" and run_rec["vina"]["version"] == "1.2.7"
    assert [p["rmsd_to_crystal_pose_angstrom"] for p in run_rec["poses"]] == REAL_RMSD
    assert run_rec["guard"]["receptor_atoms_within_8A_of_box_centre"] > 0
    assert all(hashlib.sha256((prep / f["path"]).read_bytes()).hexdigest() == f["sha256"] for f in run_rec["files"].values())
    assert s["docking"]["032 C"]["box"] == next(p for p in s["prepared"]["pockets"] if p["label_asym_id"] == "C")["box"]
    # repeat: replaced, and the crystal conformer is recorded beside it, labelled biased
    r2 = run("3OG7", "032 C")
    assert r2["replaced_run_computed_at"] == run_rec["computed_at"]
    r3 = run("3OG7", "032 C", conformer="crystal")
    assert r3["ok"] and not r3["redocking"]["is_preregistered_control"] and "biased" in r3["redocking"]["note"]
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    assert sorted(next(x for x in man["structures"] if x["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]) == [
        "032_crystal", "032_regenerated"]
    # other exhaustiveness/seed: same threshold, but not labelled as the pre-registered control (reviewer H2)
    r4 = run("3OG7", "032 C", exhaustiveness=8, seed=7)
    assert r4["ok"] and not r4["redocking"]["is_preregistered_control"] and "not the pre-registered control" in r4["redocking"]["note"]
    assert "exhaustiveness 32, seed 42" in r4["redocking"]["criterion"]
    # a superposing RMSD in production is caught by the runtime negative control (reviewer H3)
    real = dock.rmsd_to_crystal
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    dock.rmsd_to_crystal = lambda poses, ref: [round(rdMolAlign.GetBestRMS(Chem.RemoveHs(poses), Chem.RemoveHs(ref),
                                                                           prbId=c.GetId()), 4) for c in poses.GetConformers()]
    try:
        bad = run("3OG7", "032 C")
    finally:
        dock.rmsd_to_crystal = real
    assert not bad["ok"] and "negative control failed" in bad["error"], bad
    # from the user's choice, which names chain A only: the receptor has A and B, and the reply says so (H6)
    choice.select_structure(ACC, "3OG7", "032 C", chains=["A"], reason="test")
    r5 = run()
    assert "(from user_choice)" in r5["docked"] and r5["receptor_preparation"]["receptor_chains"] == ["A", "B"]
    assert "chains A, B" in r5["receptor_preparation"]["chains_note"], r5["receptor_preparation"]
    assert "chains_note" not in r["receptor_preparation"]          # no user_choice: no note
    return prep


def test_manifest_rebuilt(record: Path):
    prep, m = make_prep()
    os.environ["FAKE_VINA_TOUCH_MANIFEST"] = str(prep / "manifest.json")
    try:
        r = run("3OG7", "032 C")
    finally:
        del os.environ["FAKE_VINA_TOUCH_MANIFEST"]
    assert not r["ok"] and "rewrote the manifest" in r["error"] and "nothing was written" in r["error"], r
    assert not (prep / "dock").exists(), "files copied to the dossier although the manifest was rebuilt (reviewer H8)"
    os.environ["FAKE_VINA_FAIL"] = "1"
    try:
        r = run("3OG7", "032 C")
    finally:
        del os.environ["FAKE_VINA_FAIL"]
    assert not r["ok"] and "Vina failed (exit code 1)" in r["error"] and "something broke" in r["error"], r


def _alive(pid: int) -> bool:
    out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"], capture_output=True, text=True).stdout
    return str(pid) in out


def test_cancel_kills_vina():
    """A cancelled request (MCP notifications/cancelled -> CancelledError) must not leave Vina running (reviewer H1)."""
    work = Path(tempfile.mkdtemp())
    args = [sys.executable, "-c", "import os, time; open('pid', 'w').write(str(os.getpid())); time.sleep(120)"]

    async def scenario():
        task = asyncio.create_task(dock.run_vina(args, work))
        for _ in range(100):
            await asyncio.sleep(0.1)
            if (work / "pid").exists() and (work / "pid").read_text():
                break
        await asyncio.sleep(0.5)
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass
    asyncio.run(scenario())
    pid = int((work / "pid").read_text())
    assert not _alive(pid), f"child {pid} still running after the request was cancelled"


def test_stdout_clean():
    """The tool must write nothing to stdout (the MCP stdio): run it in a child process."""
    prep, _ = make_prep()
    code = ("import asyncio; from pathlib import Path; from protein_mcp import dossier; "
            "from protein_mcp.files import dock; "
            f"dossier.DOSSIER_DIR = Path({str(dossier.DOSSIER_DIR)!r}); "
            f"r = asyncio.run(dock.dock_ligand({ACC!r}, '3OG7', '032 C')); assert r['ok'], r")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, env=os.environ.copy())
    assert p.returncode == 0, p.stderr.decode(errors="replace")[-2000:]
    assert p.stdout == b"", p.stdout[:500]


def test_carry_over():
    box = {"center": [1, 2, 3], "size": [20, 20, 20]}
    old = {"structures": [{"pdb_id": "3OG7", "docking": {"032 C": {"raw_sha256": "a", "box": box, "runs": {}}}},
                          {"pdb_id": "8C7Y", "docking": {"TXV U": {"raw_sha256": "b", "box": box, "runs": {}}}},
                          {"pdb_id": "4XV1", "docking": {"904 C": {"raw_sha256": "c", "box": box, "runs": {}}}}]}
    pocket = lambda c, l, b: {"comp_id": c, "label_asym_id": l, "box": b}
    new = {"structures": [
        {"pdb_id": "3OG7", "file": {"ok": True, "sha256": "a"}, "prepared": {"pockets": [pocket("032", "C", box)]}},
        {"pdb_id": "8C7Y", "file": {"ok": True, "sha256": "b"},
         "prepared": {"pockets": [pocket("TXV", "U", {**box, "size": [21, 20, 20]})]}}]}
    st = dock.carry_over(old, new)
    assert st["kept"] == ["3OG7 032 C"], st
    assert st["dropped"] == ["8C7Y TXV U: structure file or box changed, or pocket gone", "4XV1 904 C: not in the new manifest"]
    assert "docking" in new["structures"][0] and "docking" not in new["structures"][1]
    assert dock.carry_over(None, {"structures": []}) == {"kept": [], "dropped": []}


# ---------------------------------------------------------------- step 2.1: background runs

def _wait_state(state_file: Path, wanted, timeout=120.0):
    """Wait for the state file to reach one of `wanted` (what dock_status would report)."""
    import time
    deadline = time.monotonic() + timeout
    last = None
    while time.monotonic() < deadline:
        last = dock.effective_state(dock._read_json(state_file) or {})
        if last in wanted:
            return last
        time.sleep(0.4)
    log = state_file.parent / dock.WORKER_LOG
    raise AssertionError(
        f"state stayed {last!r}, never reached {wanted}.\nstate file: "
        f"{state_file.read_text(encoding='utf-8') if state_file.exists() else 'none'}\n"
        f"{dock.WORKER_LOG}: "
        f"{log.read_text(encoding='utf-8', errors='replace')[-3000:] if log.exists() else 'none'}")


def _status(**kw):
    return asyncio.run(dock.dock_status(ACC, **kw))


def _start(**kw):
    return asyncio.run(dock.start_dock(ACC, **kw))


def test_state_and_stale(tmp_path=None):
    """A recorded state means what it says, except that a run without a heartbeat is vanished."""
    tmp = Path(tempfile.mkdtemp()) / "run_regenerated.state.json"
    st = dock.RunState(tmp, {"state": "queued", "launched_at": dock.utc_now()})
    st.update(state="running", pid=os.getpid(), started_at=dock.utc_now())
    data = dock._read_json(tmp)
    assert data["state"] == "running" and data["pid"] == os.getpid() and data["heartbeat"], data
    assert dock.effective_state(data) == "running"
    for final in ("finished", "failed", "stopped"):
        assert dock.effective_state({**data, "state": final, "heartbeat": "2020-01-01T00:00:00Z"}) == final
    # a late heartbeat on a process that is still alive is NOT death: calling it vanished would
    # hand its lock to a second Vina (the state file keeps the worker's own pid for this)
    stale_but_alive = {**data, "heartbeat": "2020-01-01T00:00:00Z", "pid": os.getpid()}
    assert dock.effective_state(stale_but_alive) == "running", stale_but_alive
    old = {**data, "heartbeat": "2020-01-01T00:00:00Z", "pid": 999999}
    assert dock.effective_state(old) == "vanished", "a run whose process is gone must not read as running"
    assert dock.effective_state({}) == "vanished" and dock.effective_state({"state": "running"}) == "vanished"
    assert dock._age_seconds("nonsense") is None and dock._age_seconds(None) is None
    # the stop marker is what the worker polls
    assert dock.stop_marker(tmp).name == "run_regenerated.stop"
    assert not st.stop_requested()
    dock.stop_marker(tmp).touch()
    assert st.stop_requested()


def test_lock():
    """One Vina at a time on this machine: the second request is refused naming the first, and a
    lock whose run left no heartbeat is taken over."""
    prep, _ = make_prep()
    state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
    live = dock.RunState(state_file, {"state": "running", "launched_at": dock.utc_now(),
                                      "started_at": dock.utc_now(), "elapsed_seconds": 12.0})
    live.update()
    token = dock.acquire_lock(ACC, {"label": "P15056 3OG7 032 C (regenerated)", "state_file": str(state_file)})
    assert dock.lock_path(ACC).exists()
    try:
        dock.acquire_lock(ACC, {"label": "other", "state_file": str(state_file)})
        raise AssertionError("a second lock was granted while a run was live")
    except dock.DockError as exc:
        assert "already in progress" in str(exc) and "3OG7 032 C" in str(exc), exc
        assert exc.extra.get("running_run") == "P15056 3OG7 032 C (regenerated)", exc.extra
    # someone else's lock is never removed
    dock.release_lock(ACC, "not-my-token")
    assert dock.lock_path(ACC).exists(), "release_lock removed a lock that was not its own"
    # the holder stops reporting (shutdown, kill): the lock is taken over, not held for ever.
    # RunState.update always stamps a fresh heartbeat, so the stale one is written to the file.
    stale = dock._read_json(state_file)
    stale["heartbeat"] = "2020-01-01T00:00:00Z"
    state_file.write_text(json.dumps(stale), encoding="utf-8")
    # and the lock itself has to be older than the launch window, or a state file older than the
    # lock reads as "the launcher has not written it yet" and the holder counts as live
    held = dock._read_json(dock.lock_path(ACC))
    held["locked_at"] = "2020-01-01T00:00:01Z"
    dock.lock_path(ACC).write_text(json.dumps(held), encoding="utf-8")
    second = dock.acquire_lock(ACC, {"label": "takeover", "state_file": str(state_file)})
    assert second != token and dock._read_json(dock.lock_path(ACC))["label"] == "takeover"
    dock.release_lock(ACC, second)
    assert not dock.lock_path(ACC).exists()


def test_cli_takes_the_lock(record: Path):
    """A run started by hand goes through the same machine-wide lock, so the command line and the
    chat cannot both have Vina going (the race ESTADO listed as open)."""
    prep, _ = make_prep()
    argv = [ACC, "--structure", "3OG7", "--pocket", "032 C"]
    # a live background run blocks it
    state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
    live = dock.RunState(state_file, {"state": "running", "launched_at": dock.utc_now(),
                                      "started_at": dock.utc_now()})
    live.update()
    token = dock.acquire_lock(ACC, {"label": "background run", "state_file": str(state_file)})
    assert dock.main(argv) == 1, "the CLI ran although a background run held the lock"
    dock.release_lock(ACC, token)
    # a lock held by a process that is gone does not block it, and the CLI leaves no lock behind
    dock.acquire_lock(ACC, {"label": "dead run by hand", "pid": 999999})
    assert dock.main(argv) == 0
    assert not dock.lock_path(ACC).exists(), "the CLI did not release the lock"
    m = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    assert "032_regenerated" in next(s for s in m["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]


def test_two_step_run(record: Path):
    """dock_ligand launches a real detached worker (with the Vina stand-in) and dock_status reads
    the result from the manifest; a second request is refused while it runs."""
    prep, _ = make_prep()
    os.environ["PROTEIN_MCP_DOSSIER_DIR"] = str(dossier.DOSSIER_DIR)
    os.environ["FAKE_VINA_SLEEP"] = "6"
    try:
        empty = _status()
        assert not empty["ok"] and "No docking run is on record" in empty["error"], empty

        r = _start(structure="3OG7", pocket="032 C")
        assert r["ok"] and r["launched"] and r["run"]["state"] in ("queued", "running"), r
        assert r["run"]["launch_flags"] in ("detached+breakaway", "detached"), r["run"]
        assert r["box"]["size"] and "dock_status" in r["next"], r
        assert not (prep / "dock" / "3OG7_032_C" / "poses_032_regenerated.sdf").exists(), \
            "the launching call must not have written results yet"
        state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
        launcher_pid = r["run"]["launcher_pid"]          # must exist: the comparison below needs it
        worker_pid = dock._read_json(state_file)["pid"]
        assert isinstance(launcher_pid, int) and isinstance(worker_pid, int), (launcher_pid, worker_pid)
        assert worker_pid != launcher_pid, \
            ("the worker must record its OWN pid: python.exe of the .venv is a trampoline, so the "
             f"launcher's pid ({launcher_pid}) is not the process that docks")

        # while it is live: refused, and the running run is named
        again = _start(structure="3OG7", pocket="032 C")
        assert not again["ok"] and "already in progress" in again["error"], again
        # another pocket/conformer is refused too: the lock is for the whole machine, not the run.
        # (validation comes first, so an invalid request is told what is wrong, not "busy")
        other = _start(structure="3OG7", pocket="032 C", conformer="crystal")
        assert not other["ok"] and "already in progress" in other["error"], other
        assert not other["launched"]
        invalid = _start(structure="8C7Y", pocket="TXV U")
        assert not invalid["ok"] and "is missing or outside prep/" in invalid["error"], invalid

        running = _status(structure="3OG7")["runs"][0]
        assert running["state"] in ("queued", "running") and running["worker_pid"], running
        assert "running" in (running.get("note") or ""), running

        assert _wait_state(state_file, ("finished",)) == "finished"
        st = _status()
        assert st["ok"] and st["live_runs"] == 0 and "1 of 1 docking runs" in st["selection"], st
        done = st["runs"][0]
        assert done["state"] == "finished", done
        assert done["redocking"]["top_pose_rmsd_angstrom"] == REAL_RMSD[0], done["redocking"]
        assert done["redocking"]["is_preregistered_control"] is True
        assert done["run"]["structure"] == "3OG7" and done["run"]["conformer"] == "regenerated"
        # the durable record is the manifest, and the result file is next to the output
        m = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
        run = next(s for s in m["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]["032_regenerated"]
        assert run["redocking"]["top_pose_rmsd_angstrom"] == REAL_RMSD[0]
        result_file = prep / "dock" / "3OG7_032_C" / "run_regenerated.result.json"
        assert json.loads(result_file.read_text(encoding="utf-8"))["ok"]
        assert not dock.lock_path(ACC).exists(), "the lock was not released after the run finished"
        # with the run over, a new one is allowed again
        third = _start(structure="3OG7", pocket="032 C")
        assert third["ok"], third
        _wait_state(state_file, ("finished", "failed"))
    finally:
        os.environ.pop("FAKE_VINA_SLEEP", None)
        os.environ.pop("PROTEIN_MCP_DOSSIER_DIR", None)


def test_stop_background_run(record: Path):
    """dock_status(stop=True) ends the live run: nothing is written and the lock is freed."""
    prep, _ = make_prep()
    os.environ["PROTEIN_MCP_DOSSIER_DIR"] = str(dossier.DOSSIER_DIR)
    os.environ["FAKE_VINA_SLEEP"] = "90"
    work_root = Path(dock.ascii_work_root())
    before = set(work_root.glob("protein_mcp_dock_*"))    # ignore leftovers of earlier work
    try:
        r = _start(structure="3OG7", pocket="032 C")
        assert r["ok"], r
        state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
        _wait_state(state_file, ("running",))
        for _ in range(60):                 # wait until Vina itself was started
            if (dock._read_json(state_file) or {}).get("vina_pid"):
                break
            __import__("time").sleep(0.4)
        worker_pid = dock._read_json(state_file)["pid"]
        out = _status(stop=True)
        assert out["ok"] and out["stopped"] and out["run"]["state"] == "stopped", out
        assert "Vina was killed" in out["note"], out
        assert "Stopped on request" in (out["run"]["error"] or ""), out["run"]
        assert not _alive(worker_pid), f"the worker {worker_pid} is still running after stop"
        m = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
        assert "docking" not in next(s for s in m["structures"] if s["pdb_id"] == "3OG7"), \
            "a stopped run must write nothing to the manifest"
        assert not dock.lock_path(ACC).exists(), "a stopped run must release the lock"
        # a grandchild of Vina (here the stand-in behind the .cmd) would keep the folder locked
        assert not set(work_root.glob("protein_mcp_dock_*")) - before, \
            "the temporary folder of a stopped run was left behind"
        # stopping again when nothing is live says so instead of failing
        nothing = _status(stop=True)
        assert nothing["ok"] and nothing["stopped"] is False and "Nothing to stop" in nothing["note"], nothing
    finally:
        os.environ.pop("FAKE_VINA_SLEEP", None)
        os.environ.pop("PROTEIN_MCP_DOSSIER_DIR", None)


def test_status_vanished():
    """A run whose process died without writing anything is reported as vanished, not running."""
    prep, _ = make_prep()
    state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
    # as the state file is left when the computer is shut down in the middle of a run
    state_file.parent.mkdir(parents=True, exist_ok=True)
    state_file.write_text(json.dumps(
        {"state": "running", "uniprot_id": ACC, "structure": "3OG7", "pocket": "032 C",
         "pocket_key": "032 C", "run_id": "032_regenerated", "conformer": "regenerated",
         "pid": 999999, "phase": "Vina running",
         "launched_at": dock.utc_now(), "heartbeat": "2020-01-01T00:00:00Z"}), encoding="utf-8")
    st = _status()
    run = st["runs"][0]
    assert run["state"] == "vanished", run
    # the claim is now backed by looking: the manifest has no block from this run
    assert "shut down" in run["error"] and "no block from this run" in run["error"], run
    assert run["last_phase"] == "Vina running", run
    assert _status(structure="8C7Y")["ok"] is False, "a filter that matches nothing must say so"

    # a state file that does not say which run it was cannot claim anything about the dossier
    nameless = json.loads(state_file.read_text(encoding="utf-8"))
    nameless.pop("pocket_key"), nameless.pop("run_id")
    state_file.write_text(json.dumps(nameless), encoding="utf-8")
    run = _status()["runs"][0]
    assert run["state"] == "vanished" and "cannot tell whether it wrote anything" in run["error"], run
    assert "nothing was kept" not in run["error"], run["error"]

    # a worker that never reported is a different fact, and gets a longer grace (STARTUP_GRACE)
    queued = {"state": "queued", "uniprot_id": ACC, "structure": "3OG7", "pocket": "032 C",
              "pocket_key": "032 C", "run_id": "032_regenerated", "conformer": "regenerated",
              "launched_at": dock.utc_now(), "heartbeat": dock.utc_now()}
    state_file.write_text(json.dumps(queued), encoding="utf-8")
    assert dock.effective_state(queued) == "queued"
    old_queued = {**queued, "heartbeat": "2020-01-01T00:00:00Z"}
    state_file.write_text(json.dumps(old_queued), encoding="utf-8")
    run = _status()["runs"][0]
    assert run["state"] == "vanished" and "never started" in run["error"], run
    assert dock.WORKER_LOG in run["error"], run
    # the same age while it was running is reported as a run that stopped, not one that never began
    state_file.write_text(json.dumps({**old_queued, "state": "running", "pid": 999999}), encoding="utf-8")
    assert "stopped without a result" in _status()["runs"][0]["error"]


def test_results_without_a_state_file(record: Path):
    """A run recorded in the manifest with no state file (command line, or an older version) is
    reported as the finished run it is, not as "nothing was ever launched" (reviewer H2)."""
    prep, _ = make_prep()
    assert dock.main([ACC, "--structure", "3OG7", "--pocket", "032 C"]) == 0
    assert not list((prep / "dock").glob("*/run_*.state.json")), "the CLI should write no state file"
    st = _status()
    assert st["ok"] and st["live_runs"] == 0, st
    assert "1 of 1 docking runs on record" in st["selection"], st["selection"]
    run = st["runs"][0]
    assert run["state"] == "finished" and run["state_file"] is None, run
    assert run["redocking"]["top_pose_rmsd_angstrom"] == REAL_RMSD[0], run["redocking"]
    assert "not started by dock_ligand" in run["note_on_state"], run["note_on_state"]
    assert _status(structure="3OG7", conformer="regenerated")["runs"][0]["state"] == "finished"
    assert _status(conformer="crystal")["ok"] is False, "a filter that matches nothing must still say so"


def test_vanished_but_written(record: Path):
    """A worker that died after writing the manifest and before recording "finished" must not be
    told to run again as if nothing had been kept (reviewer H5)."""
    prep, _ = make_prep()
    launched_at = dock.utc_now()               # the run starts BEFORE it writes the manifest
    state = {"state": "running", "uniprot_id": ACC, "structure": "3OG7", "pocket": "032 C",
             "pocket_key": "032 C", "run_id": "032_regenerated", "conformer": "regenerated",
             "pid": 999999, "phase": "reading the poses and writing the dossier",
             "launched_at": launched_at, "started_at": launched_at,
             "heartbeat": "2020-01-01T00:00:00Z"}
    assert dock.main([ACC, "--structure", "3OG7", "--pocket", "032 C"]) == 0   # writes the manifest
    state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
    state_file.write_text(json.dumps(state), encoding="utf-8")
    run = _status()["runs"][0]
    assert run["state"] == "finished", run
    assert run["redocking"]["top_pose_rmsd_angstrom"] == REAL_RMSD[0], run["redocking"]
    assert "without recording the end" in run["note_on_state"], run["note_on_state"]
    assert "error" not in run, run.get("error")

    # the block of an EARLIER run must NOT be handed to a later run that died: same pocket and
    # conformer overwrite the same key, so only `computed_at` tells them apart (reviewer F2)
    later = dict(state, launched_at=dock.utc_now(), started_at=dock.utc_now())
    m = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    block = next(s for s in m["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]
    block["032_regenerated"]["computed_at"] = "2020-01-01T00:00:00Z"      # written long before
    (prep / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    state_file.write_text(json.dumps(later), encoding="utf-8")
    run = _status()["runs"][0]
    assert run["state"] == "vanished", run
    assert "no block from this run" in run["error"], run["error"]
    assert "redocking" not in run, "the RMSD of an earlier run was reported as this run's result"


def test_lock_not_held_for_ever(record: Path):
    """A background lock whose state file is gone is expired, not alive: judging it by the pid in
    the lock would be judging the MCP server and would block every later run (reviewer H3)."""
    prep, _ = make_prep()
    state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
    token = dock.acquire_lock(ACC, {"label": "dead background run", "state_file": str(state_file)})
    assert dock.lock_path(ACC).exists() and not state_file.exists()
    held = dock._read_json(dock.lock_path(ACC))
    # a lock just taken counts as live: the launcher writes the state file right afterwards
    assert dock._lock_holder_live(held) is True, held
    # but one taken long ago whose state file cannot be read is expired, not a live run
    held["locked_at"] = "2020-01-01T00:00:00Z"
    dock.lock_path(ACC).write_text(json.dumps(held), encoding="utf-8")
    assert dock._lock_holder_live(dock._read_json(dock.lock_path(ACC))) is False, \
        "a lock whose state file cannot be read must not count as a live run"
    second = dock.acquire_lock(ACC, {"label": "new run", "state_file": str(state_file)})
    assert second != token
    dock.release_lock(ACC, second)


def test_recycled_pid_is_not_a_live_run():
    """Windows reuses process ids: a dead run whose id went to another program must not be
    reported as running, so the image name is checked too (reviewer H7)."""
    assert dock._pid_alive(os.getpid(), expect="python") is True
    assert dock._pid_alive(os.getpid(), expect="vina_1.2.7_win") is False, \
        "the image name is not being checked: any process with that id would pass"
    assert dock._pid_alive(999999, expect="python") is False
    assert dock._pid_alive(None) is False
    # a run with a late heartbeat whose pid is alive stays "running", but says so
    prep, _ = make_prep()
    state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
    state_file.parent.mkdir(parents=True, exist_ok=True)
    state_file.write_text(json.dumps(
        {"state": "running", "uniprot_id": ACC, "structure": "3OG7", "pocket": "032 C",
         "pocket_key": "032 C", "run_id": "032_regenerated", "conformer": "regenerated",
         "pid": os.getpid(), "phase": "Vina running", "elapsed_seconds": 55.0,
         "launched_at": dock.utc_now(), "heartbeat": "2020-01-01T00:00:00Z"}), encoding="utf-8")
    run = _status()["runs"][0]
    assert run["state"] == "running", run
    assert "has not reported" in run["warning"] and "same id" in run["warning"], run.get("warning")
    assert run["heartbeat_age_seconds"] > dock.STALE_AFTER


def test_stop_after_it_finished(record: Path):
    """stop=True on a run that has just finished must not claim it was stopped, nor tell the user
    to kill it from the Task Manager (reviewer H4)."""
    prep, _ = make_prep()
    assert dock.main([ACC, "--structure", "3OG7", "--pocket", "032 C"]) == 0
    state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
    state_file.write_text(json.dumps(
        {"state": "running", "uniprot_id": ACC, "structure": "3OG7", "pocket": "032 C",
         "pocket_key": "032 C", "run_id": "032_regenerated", "conformer": "regenerated",
         "pid": os.getpid(), "phase": "Vina running", "launched_at": dock.utc_now(),
         "heartbeat": dock.utc_now()}), encoding="utf-8")
    # it is live when stop is called and records "finished" while the stop waits, which is the
    # race the reviewer described: the result is already written
    def finish():
        data = dock._read_json(state_file)
        data.update(state="finished", heartbeat=dock.utc_now(), finished_at=dock.utc_now())
        state_file.write_text(json.dumps(data), encoding="utf-8")

    timer = __import__("threading").Timer(2.0, finish)
    timer.start()
    try:
        out = _status(stop=True)
    finally:
        timer.cancel()
    assert out["ok"] and out["stopped"] is False, out
    assert "NOT stopped" in out["note"] and "Task Manager" not in out["note"], out["note"]
    assert out["run"]["state"] == "finished"
    assert out["run"]["redocking"]["top_pose_rmsd_angstrom"] == REAL_RMSD[0]


def test_unreadable_state_says_so(record: Path):
    """An unreadable state file is not evidence: no story about shutdowns or about the dossier
    being empty, and never ok: true (reviewer F3)."""
    prep, _ = make_prep()
    assert dock.main([ACC, "--structure", "3OG7", "--pocket", "032 C"]) == 0   # a real result exists
    state_file = dock.state_path(prep, "3OG7", "032", "C", "regenerated")
    state_file.parent.mkdir(parents=True, exist_ok=True)
    state_file.write_text("{ truncated", encoding="utf-8")
    run = _status()["runs"][0]
    assert run["state"] == "unknown" and run["ok"] is False, run
    assert "could not be read" in run["error"] and "unknown" in run["error"], run["error"]
    for invented in ("shut down", "nothing was kept", "s ago"):
        assert invented not in run["error"], (invented, run["error"])


def test_selection_count_with_a_filter(record: Path):
    """The "N of total" of dock_status must not count a filtered-out state file twice (F4)."""
    prep, _ = make_prep()
    for pdb, pocket in (("3OG7", "032 C"), ("8C7Y", "TXV U")):
        state_file = dock.state_path(prep, pdb, *pocket.split(), "regenerated")
        dock.RunState(state_file, {"state": "finished", "uniprot_id": ACC, "structure": pdb,
                                   "pocket": pocket, "pocket_key": pocket,
                                   "run_id": f"{pocket.split()[0]}_regenerated",
                                   "conformer": "regenerated"}).update()
    assert "2 of 2 docking runs" in _status()["selection"]
    assert "1 of 2 docking runs" in _status(structure="3OG7")["selection"], _status(structure="3OG7")["selection"]


def test_refusal_keeps_the_previous_result(record: Path):
    """A request refused by the lock must not destroy the previous run's result file (F5)."""
    prep, _ = make_prep()
    os.environ["PROTEIN_MCP_DOSSIER_DIR"] = str(dossier.DOSSIER_DIR)
    try:
        out_dir = prep / "dock" / "3OG7_032_C"
        out_dir.mkdir(parents=True)
        result_file = out_dir / "run_regenerated.result.json"
        result_file.write_text('{"ok": true, "kept": "the previous run"}', encoding="utf-8")
        other_state = dock.state_path(prep, "8C7Y", "TXV", "U", "regenerated")
        dock.RunState(other_state, {"state": "running", "launched_at": dock.utc_now(),
                                    "started_at": dock.utc_now(), "pid": os.getpid()}).update()
        token = dock.acquire_lock(ACC, {"label": "another run", "state_file": str(other_state)})
        try:
            r = _start(structure="3OG7", pocket="032 C")
            assert not r["ok"] and "already in progress" in r["error"], r
            assert result_file.exists(), "a refused request deleted the previous run's result file"
            assert json.loads(result_file.read_text(encoding="utf-8"))["kept"] == "the previous run"
        finally:
            dock.release_lock(ACC, token)
    finally:
        os.environ.pop("PROTEIN_MCP_DOSSIER_DIR", None)


def test_cpu_zero_is_not_a_measurement(record: Path):
    """Vina reached through a wrapper measures 0 s of CPU; reporting 0 would read as a starved
    Vina, so it is reported as not measurable (reviewer H12)."""
    prep, _ = make_prep()
    r = run("3OG7", "032 C")
    assert r["ok"], r.get("error")
    vina = r["vina"]
    # the stand-in is a .cmd, i.e. a trampoline, whatever its wall time: no guard on wall_seconds,
    # which is what made the first version of this test vacuous (reviewer F1)
    assert vina.get("cpu_seconds") is None, vina
    assert "Not measurable" in vina.get("cpu_seconds_note", ""), vina
    assert "cpu_per_wall_second" not in vina, vina
    # and the measurement itself works on a process that is not a trampoline: the .venv python.exe
    # IS one (it reads ~0), which is why Vina's own handle is the one that counts
    burn = "x=0\nfor i in range(8_000_000): x+=i*i\n"
    real = subprocess.Popen([sys._base_executable, "-I", "-c", burn])
    real.wait()
    assert dock._cpu_seconds(real._handle) > 0.1, "GetProcessTimes measured no CPU for a CPU-bound process"
    stub = subprocess.Popen([sys.executable, "-I", "-c", burn])
    stub.wait()
    assert dock._cpu_seconds(stub._handle) < 0.1, \
        "the .venv python.exe no longer trampolines: the note about wrappers needs rechecking"


def main():
    record = make_fake_vina(Path(tempfile.mkdtemp()))
    test_parse_meeko()
    test_poses_and_rmsd()
    test_guard()
    test_ascii_root()
    test_carry_over()
    test_refusals(record)
    test_full_run(record)
    test_full_run(record, "dossiers ñá")      # non-ASCII dossier path (ESTADO 8.3)
    test_manifest_rebuilt(record)
    test_cancel_kills_vina()
    test_stdout_clean()
    test_state_and_stale()
    test_lock()
    test_cli_takes_the_lock(record)
    test_two_step_run(record)
    test_stop_background_run(record)
    test_status_vanished()
    test_results_without_a_state_file(record)
    test_vanished_but_written(record)
    test_lock_not_held_for_ever(record)
    test_recycled_pid_is_not_a_live_run()
    test_stop_after_it_finished(record)
    test_unreadable_state_says_so(record)
    test_selection_count_with_a_filter(record)
    test_refusal_keeps_the_previous_result(record)
    test_cpu_zero_is_not_a_measurement(record)
    print("ALL OK")


if __name__ == "__main__":
    main()
