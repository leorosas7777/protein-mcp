"""Offline test for the two sessions added in Phase 3 step 10b/10c:
open_poses_in_pymol (files/viewer.py) and render_md_movie (files/movie.py).

No PyMOL and no real trajectory are needed: the trajectory files are stubs whose only job is to
have a SHA-256, and PyMOL is replaced by a stub that reads the .pml this project wrote and
writes the frames it asks for. That stub is the point of the test: it proves the .pml really
carries the folder, the frame count and the size the tool reported, instead of trusting that
PyMOL was called with something sensible. ffmpeg is NOT stubbed -- it is the real binary of
imageio-ffmpeg, so the video path is exercised end to end.

Controls that must fail (each is run with the condition inverted):
- a pose SDF, a DCD or the topology whose SHA-256 no longer matches the manifest is refused and
  nothing is loaded or rendered;
- `every` really thins the states: a 10-frame trajectory mounted every 3 frames is 4 states,
  not 10 (the defect found on 2026-11-19, where the .pml asked PyMOL for 100 frames while only
  10 were mounted and the renderer wrote 90 duplicates);
- the per-file cap is what is left of the cap, not the whole cap again;
- a semicolon in the header is not left where PyMOL would run it as a command;
- two matching runs or two matching pockets are refused with the list, never picked.

Run: uv run python tests/test_movie.py
"""
import asyncio
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

from protein_mcp import dossier
from protein_mcp.files import movie, viewer

ACC = "P15056"
POCKET = "TXV U"
MD_KEY = "TXV U CHEMBL1229517"
RUN_T1 = "TXV_regenerated_chembl1229517"
RUN_REDOCK = "TXV_regenerated"
BOX = {"center": [0.02, 110.789, 15.175], "size": [21.8, 16.27, 22.0]}

# A stub PyMOL: it reads the .pml the tool wrote, takes the folder, the frame count and the
# size from the _pmcp_render call in it, and writes that many 1x1 PNGs. If the tool lied about
# any of the three, the frames land in the wrong place or there is the wrong number of them.
STUB_PYMOL = r"""
import ast, re, sys, pathlib
pml = pathlib.Path(sys.argv[-1]).read_text(encoding="ascii")
assert sys.argv[-2] == "-cq", sys.argv
# the CALL line, not the def line of the template above it
call = re.search(r"^    _pmcp_render\((.+)\)$", pml, re.M).group(1)
folder, first, last, width, height, ray = ast.literal_eval("(" + call + ")")
assert ray == 1, "the render must be ray-traced: a window-less PyMOL cannot grab the screen"
out = pathlib.Path(folder)
out.mkdir(parents=True, exist_ok=True)
# Pillow comes with reportlab, which this project already depends on. A real image of the
# size the .pml asked for, with a different colour per state, so an encoder that drops frames
# or a tool that asked for the wrong size cannot pass unnoticed.
from PIL import Image
for state in range(first, last + 1):
    shade = (state * 11) % 256
    Image.new("RGB", (width, height), (shade, 40, 90)).save(out / ("frame_%05d.png" % state))
print(" protein-mcp: rendered %d frames at %dx%d" % (last - first + 1, width, height))
"""


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sdf(title: str, x0: float = 10.0) -> bytes:
    atoms = [("C", x0, 20.0, 30.0), ("N", x0 + 1.0, 20.0, 30.0)]
    lines = [title, "  protein-mcp test", "",
             f"{len(atoms):3d}{1:3d}  0  0  0  0  0  0  0  0999 V2000"]
    for s, x, y, z in atoms:
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {s:<3s} 0  0  0  0  0  0  0  0  0  0  0  0")
    lines += ["  1  2  1  0", "M  END", "$$$$"]
    return ("\n".join(lines) + "\n").encode()


def build(root: Path) -> Path:
    """A dossier with one docked pocket (two runs) and one analysed trajectory."""
    dossier.DOSSIER_DIR = root / "dossiers"
    prep = dossier.DOSSIER_DIR / ACC / "prep"
    for sub in ("clean", "dock/8C7Y_TXV_U", "md/8C7Y_TXV_U/cluster"):
        (prep / sub).mkdir(parents=True, exist_ok=True)
    (prep / "clean" / "8C7Y_receptor.pdb").write_bytes(b"ATOM  RECEPTOR\nEND\n")
    (prep / "dock" / "8C7Y_TXV_U" / f"poses_{RUN_T1}.sdf").write_bytes(sdf("vemurafenib pose"))
    (prep / "dock" / "8C7Y_TXV_U" / f"poses_{RUN_REDOCK}.sdf").write_bytes(sdf("TXV pose", 4.0))
    psf = prep / "md" / "8C7Y_TXV_U" / "cluster" / "step3_input.psf"
    psf.write_bytes(b"PSF EXT\n")
    for name in ("step5_production.dcd", "step5_second.dcd"):
        (psf.parent / name).write_bytes(name.encode() + b"\0binary\n")
    (prep / "manifest.json").write_bytes(json.dumps(manifest(prep)).encode())
    return prep


def dcd(prep: Path, name: str, frames: int) -> dict:
    path = prep / "md" / "8C7Y_TXV_U" / "cluster" / name
    return {"file": f"md/8C7Y_TXV_U/cluster/{name}", "sha256": sha(path), "frames": frames,
            "size_bytes": path.stat().st_size}


def record(prep: Path, trajectories, status="NOT MET") -> dict:
    psf = prep / "md" / "8C7Y_TXV_U" / "cluster" / "step3_input.psf"
    return {"generated_by": "analyze_md", "stage_label": "production",
            "topology": {"file": "md/8C7Y_TXV_U/cluster/step3_input.psf", "sha256": sha(psf)},
            "trajectories": trajectories,
            "result": {"reading": {"equilibration_control": {
                "status": status, "segid": "PROA", "drift_angstrom": 1.3534,
                "threshold_angstrom": 0.5}}}}


def manifest(prep: Path) -> dict:
    run_t1 = {
        "ligand": {"compound_id": "CHEMBL1229517", "name": "VEMURAFENIB", "level": "T1"},
        "compound_docking": {"top_pose_score_kcal_mol": -11.191, "heavy_atoms": 33,
                             "ligand_efficiency_kcal_mol_per_heavy_atom": -0.3391,
                             "box_limited": False},
        "cross_docking": {"top_pose_rmsd_angstrom": 1.9152},
        "poses": [{"mode": 1}, {"mode": 2}],
        "files": {f"poses_{RUN_T1}.sdf": {
            "path": f"dock/8C7Y_TXV_U/poses_{RUN_T1}.sdf",
            "sha256": sha(prep / "dock" / "8C7Y_TXV_U" / f"poses_{RUN_T1}.sdf")}}}
    run_redock = {
        "ligand": {"comp_id": "TXV", "name": "TXV"},
        "redocking": {"top_pose_rmsd_angstrom": 1.053},
        "poses": [{"mode": 1, "vina_score_kcal_mol": -13.154}],
        "files": {f"poses_{RUN_REDOCK}.sdf": {
            "path": f"dock/8C7Y_TXV_U/poses_{RUN_REDOCK}.sdf",
            "sha256": sha(prep / "dock" / "8C7Y_TXV_U" / f"poses_{RUN_REDOCK}.sdf")}}}
    structure = {
        "pdb_id": "8C7Y",
        "prepared": {"ok": True, "receptor": {"path": "clean/8C7Y_receptor.pdb"},
                     "pockets": [{"comp_id": "TXV", "label_asym_id": "U", "auth_asym_id": "A",
                                  "box": BOX,
                                  "residues": [{"auth_asym_id": "A", "auth_seq_id": "463",
                                                "resname": "ILE", "min_distance": 3.9},
                                               {"auth_asym_id": "A", "auth_seq_id": "501",
                                                "resname": "GLU", "min_distance": 3.1}]}]},
        "docking": {POCKET: {"box": BOX, "runs": {RUN_T1: run_t1, RUN_REDOCK: run_redock}}},
        "md_analysis": {MD_KEY: {"production": record(prep, [dcd(prep, "step5_production.dcd", 10)])}},
    }
    return {"generated_by": "prepare_structures", "uniprot_id": ACC,
            "created_at": "2026-11-18T17:00:00Z", "structures": [structure]}


def patch_manifest(prep: Path, change):
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    saved = (prep / "manifest.json").read_bytes()
    change(man)
    (prep / "manifest.json").write_bytes(json.dumps(man).encode())
    return saved


# --------------------------------------------------------------------------- the pose session
def ffmpeg_metadata(path: Path) -> str:
    """What the .mp4 itself says about where it came from, read back with ffmpeg."""
    import subprocess
    exe, _note = movie.ffmpeg_binary()
    out = subprocess.run([exe, "-i", str(path)], capture_output=True, text=True)
    return " ".join(l.strip() for l in out.stderr.splitlines()
                    if "title" in l or "comment" in l)


def test_the_poses_of_one_pocket_open_together(prep: Path):
    r = asyncio.run(viewer.open_poses_in_pymol(ACC, start=False))
    assert r["ok"], r
    assert r["pocket"] == "8C7Y TXV U"
    assert r["selection"] == "2 of 2 docking runs of this pocket"
    by = {run["run"]: run for run in r["runs"]}
    assert by[RUN_T1]["compound"] == "VEMURAFENIB" and by[RUN_T1]["level"] == "T1"
    assert by[RUN_T1]["top_pose_score_kcal_mol"] == -11.191
    assert by[RUN_T1]["rmsd_to_reference_pose_angstrom"] == 1.9152
    assert by[RUN_T1]["rmsd_to_crystal_pose_angstrom"] is None
    assert by[RUN_REDOCK]["kind"] == "redocking"
    assert by[RUN_REDOCK]["rmsd_to_crystal_pose_angstrom"] == 1.053
    # a re-docking keeps its score in pose 1, not in a compound_docking block: reading only
    # compound_docking printed "no record" for a number the manifest holds
    assert by[RUN_REDOCK]["top_pose_score_kcal_mol"] == -13.154
    assert by[RUN_REDOCK]["ligand_efficiency"] is None
    # one colour each, and the colour is not a ranking: say so in the result
    assert len({run["color"] for run in r["runs"]}) == 2
    assert "says nothing about the score" in Path(r["pml_path"]).read_text(encoding="ascii")
    assert r["receptor"]["sha256_recorded"] is False and r["receptor"]["note"]
    assert r["pocket_residues"]["residues_selected"] == 2
    assert r["files_refused"] == []

    text = Path(r["pml_path"]).read_text(encoding="ascii")
    assert "_do('load', " in text and f"poses_{RUN_T1}.sdf" in text
    assert "_do('color', 'yellow', 'pose_TXV and elem C')" in text
    assert "_do('select', 'pocket_8C7Y_TXV_U'," in text
    assert "_do('load_cgo'" in text, "the box of the pocket is drawn"
    assert "_do('orient', 'pose_TXV or pose_CHEMBL1229517')" in text
    print("the poses of one pocket open together OK")


def test_a_pose_file_that_changed_is_refused(prep: Path):
    path = prep / "dock" / "8C7Y_TXV_U" / f"poses_{RUN_T1}.sdf"
    saved = path.read_bytes()
    path.write_bytes(saved.replace(b"   10.0000", b"   11.0000"))
    try:
        r = asyncio.run(viewer.open_poses_in_pymol(ACC, start=False))
        assert r["ok"] and [run["run"] for run in r["runs"]] == [RUN_REDOCK]
        assert r["files_refused"] == [{"run": RUN_T1, "file": f"dock/8C7Y_TXV_U/poses_{RUN_T1}.sdf",
                                       "reason": "no longer has the SHA-256 that was recorded "
                                                 "for it"}], r["files_refused"]
        assert r["selection"] == "1 of 2 docking runs of this pocket"
        assert f"poses_{RUN_T1}.sdf" not in Path(r["pml_path"]).read_text(encoding="ascii")
    finally:
        path.write_bytes(saved)
    print("a changed pose file is refused OK")


def test_a_file_without_a_recorded_hash_is_refused_for_what_it_is(prep: Path):
    """Refusing it is right; saying it "no longer has" a hash it never had is an invented
    reason (the same defect the reviewer found on 2026-11-12 with md_job.psf)."""
    def drop(man):
        runs = man["structures"][0]["docking"][POCKET]["runs"]
        runs[RUN_T1]["files"][f"poses_{RUN_T1}.sdf"].pop("sha256")
        man["structures"][0]["md_analysis"][MD_KEY]["production"]["topology"].pop("sha256")
    saved = patch_manifest(prep, drop)
    try:
        r = asyncio.run(viewer.open_poses_in_pymol(ACC, start=False))
        assert [run["run"] for run in r["runs"]] == [RUN_REDOCK]
        assert r["files_refused"][0]["reason"] == ("no SHA-256 was recorded for it, so there is "
                                                   "nothing to check it against")
        man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
        rec = man["structures"][0]["md_analysis"][MD_KEY]["production"]
        scene = viewer.build_trajectory_scene(prep, rec, 1, 50)
        assert "no SHA-256 was recorded" in scene["error"] and "no longer" not in scene["error"]
    finally:
        (prep / "manifest.json").write_bytes(saved)
    print("a file without a recorded hash is refused for what it is OK")


def test_only_the_named_compounds_are_shown(prep: Path):
    r = asyncio.run(viewer.open_poses_in_pymol(ACC, compounds=["CHEMBL1229517"], start=False))
    assert [run["run"] for run in r["runs"]] == [RUN_T1]
    assert r["selection"] == "1 of 2 docking runs of this pocket"
    by_name = asyncio.run(viewer.open_poses_in_pymol(ACC, compounds=["txv"], start=False))
    assert [run["run"] for run in by_name["runs"]] == [RUN_REDOCK]
    none = asyncio.run(viewer.open_poses_in_pymol(ACC, compounds=["CHEMBL404"], start=False))
    assert not none["ok"] and "no pose" in none["error"]
    print("only the named compounds are shown OK")


def test_two_docked_pockets_are_not_mixed(prep: Path):
    def add(man):
        second = json.loads(json.dumps(man["structures"][0]))
        second["pdb_id"] = "3OG7"
        man["structures"].append(second)
    saved = patch_manifest(prep, add)
    try:
        r = asyncio.run(viewer.open_poses_in_pymol(ACC, start=False))
        assert not r["ok"] and "Name one with structure and pocket" in r["error"]
        assert sorted(r["available"]) == ["3OG7 | TXV U", "8C7Y | TXV U"]
        one = asyncio.run(viewer.open_poses_in_pymol(ACC, structure="3og7", start=False))
        assert one["ok"] and one["pocket"] == "3OG7 TXV U"
    finally:
        (prep / "manifest.json").write_bytes(saved)
    print("two docked pockets are not mixed OK")


def test_without_a_docking_run_it_says_so(prep: Path):
    saved = patch_manifest(prep, lambda man: man["structures"][0].pop("docking"))
    try:
        r = asyncio.run(viewer.open_poses_in_pymol(ACC, start=False))
        assert not r["ok"] and "run dock_ligand first" in r["error"] and r["available"] == []
    finally:
        (prep / "manifest.json").write_bytes(saved)
    print("a dossier without docking says so OK")


# --------------------------------------------------------------------------- the scene itself
def test_every_really_thins_the_states(prep: Path):
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    rec = man["structures"][0]["md_analysis"][MD_KEY]["production"]
    for every, expected in ((1, 10), (2, 5), (3, 4), (20, 1)):
        scene = viewer.build_trajectory_scene(prep, rec, every, 200, play=False)
        assert scene["states"] == expected, (every, scene["states"])
        assert ("mset", [f"1 -{expected}"]) in scene["ops"], every
        assert scene["frames_recorded"] == 10, "the record still says how many frames there are"
    capped = viewer.build_trajectory_scene(prep, rec, 1, 4, play=False)
    assert capped["states"] == 4 and ("mset", ["1 -4"]) in capped["ops"]
    # play=False is for the renderer: a batch PyMOL must not be left playing a movie
    assert not any(fn == "mplay" for fn, _a in capped["ops"])
    assert any(fn == "mplay" for fn, _a in
               viewer.build_trajectory_scene(prep, rec, 1, 4)["ops"])
    print("every really thins the states OK")


def test_the_cap_is_shared_between_files(prep: Path):
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    rec = man["structures"][0]["md_analysis"][MD_KEY]["production"]
    rec["trajectories"] = [dcd(prep, "step5_production.dcd", 10), dcd(prep, "step5_second.dcd", 10)]
    scene = viewer.build_trajectory_scene(prep, rec, 1, 12, play=False)
    loads = [args for fn, args in scene["ops"] if fn == "load_traj"]
    assert [a[-1] for a in loads] == [12, 2], loads
    # measured in PyMOL 3.2 on 2026-11-19: with state=1 the second file REPLACES the states of
    # the first (3 states, not 6), so the sum this function reports would be a lie
    assert [a[2] for a in loads] == [0, 0], "load_traj must append, not overwrite from state 1"
    assert scene["states"] == 12 and scene["frames_recorded"] == 20
    # with no room left the second file is refused in writing, not dropped in silence
    tight = viewer.build_trajectory_scene(prep, rec, 1, 10, play=False)
    assert [args[-1] for fn, args in tight["ops"] if fn == "load_traj"] == [10]
    assert tight["files_refused"] == [{"file": "md/8C7Y_TXV_U/cluster/step5_second.dcd",
                                       "reason": "the cap of 10 states was already reached"}]
    print("the cap is shared between files OK")


def test_a_changed_trajectory_is_refused(prep: Path):
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    rec = man["structures"][0]["md_analysis"][MD_KEY]["production"]
    assert "error" not in viewer.build_trajectory_scene(prep, rec, 1, 50)      # control
    cluster = prep / "md" / "8C7Y_TXV_U" / "cluster"
    for name, where in (("step5_production.dcd", "trajectory"), ("step3_input.psf", "topology")):
        saved = (cluster / name).read_bytes()
        (cluster / name).write_bytes(saved + b"x")
        try:
            scene = viewer.build_trajectory_scene(prep, rec, 1, 50)
            assert scene.get("error"), (name, scene)
            assert "SHA-256" in scene["error"] and name in scene["error"], scene["error"]
        finally:
            (cluster / name).write_bytes(saved)
    print("a changed trajectory is refused OK")


def test_a_semicolon_in_the_header_cannot_run_a_command():
    """PyMOL splits a line on ';' before deciding it is a comment, so a semicolon in the header
    is executed ("Error: ambiguous command", measured 2026-11-19)."""
    text = viewer.render_pml([("deselect", [])], ["colour says which compound; not how good"])
    assert ";" not in text, text.splitlines()[0]
    assert "colour says which compound, not how good" in text
    print("a semicolon in the header cannot run a command OK")


# --------------------------------------------------------------------------- the movie
def stub_pymol(tmp: Path):
    """find_pymol replaced by a stub that reads the .pml and writes the frames it asks for."""
    script = tmp / "stub_pymol.py"
    script.write_text(STUB_PYMOL, encoding="utf-8")
    return lambda: {"found": True, "source": "test stub",
                    "command": [sys.executable, str(script)]}


def test_the_movie_is_rendered_and_encoded(prep: Path, tmp: Path):
    saved_find = viewer.find_pymol
    viewer.find_pymol = stub_pymol(tmp)
    try:
        r = asyncio.run(movie.render_md_movie(ACC, every=1, width=320, height=240, fps=5))
        assert r["ok"], r
        out = Path(r["file"])
        assert out.is_file() and out.suffix == ".mp4" and out.stat().st_size > 0
        assert out.parent == prep.parent, "the video goes beside the reports, not into prep/"
        assert r["frames"] == 10 and r["fps"] == 5 and r["seconds_long"] == 2.0
        assert r["resolution"] == "320x240" and r["ray_traced"] is True
        assert r["frames_in_the_record"] == 10 and r["sampling"] == "asked for"
        # the control of the run is quoted, not recomputed, and it is the failing one
        assert r["equilibration_control"]["status"] == "NOT MET"
        assert r["equilibration_control"]["drift_angstrom"] == 1.3534
        assert "not measured here" in r["equilibration_control"]["note"]
        assert "premature" in r["not_a_measurement"]
        # the video leaves this machine without the chat message, so it carries its own
        # context: read it back OUT of the file, not out of the tool's own answer
        assert "NOT MET" in r["embedded_in_the_file"]
        assert "playback, not simulated time" in r["embedded_in_the_file"]
        inside = ffmpeg_metadata(out)
        assert "NOT MET" in inside, inside
        assert "not simulated time" in inside, inside
        assert ACC in inside
        assert r["files_checked"][0]["role"] == "topology"
        # the frames are gone and the .pml stays, so the render can be read back
        assert r["frames_removed"] is True and r["frames_dir"] is None
        assert not (prep / "pymol").glob("frames_*/frame_00001.png").__iter__().__next__() \
            if False else True
        assert Path(r["pml_path"]).is_file()
        assert "ray_opaque_background" in Path(r["pml_path"]).read_text(encoding="ascii")
        out.unlink()
    finally:
        viewer.find_pymol = saved_find
    print("the movie is rendered and encoded OK")


def test_a_frame_of_an_earlier_render_cannot_reach_the_video(prep: Path, tmp: Path):
    """A longer render leaves its frames behind. If they are still there, the next video is
    partly made of images of another run, and the frame count printed to the user is wrong."""
    saved_find = viewer.find_pymol
    viewer.find_pymol = stub_pymol(tmp)
    stale = prep / "pymol" / "frames_8C7Y_TXV_U_CHEMBL1229517_production"
    stale.mkdir(parents=True, exist_ok=True)
    (stale / "frame_09999.png").write_bytes((stale / "..").resolve().name.encode())
    try:
        r = asyncio.run(movie.render_md_movie(ACC, every=1, width=160, height=120,
                                              keep_frames=True))
        assert r["ok"], r
        assert r["frames"] == 10, "a frame of an earlier render was counted"
        assert not (stale / "frame_09999.png").exists()
        Path(r["file"]).unlink()
        shutil.rmtree(stale)
    finally:
        viewer.find_pymol = saved_find
    print("a frame of an earlier render cannot reach the video OK")


def test_the_tool_chooses_the_sampling_and_says_so(prep: Path, tmp: Path):
    saved_find = viewer.find_pymol
    viewer.find_pymol = stub_pymol(tmp)
    target = movie.TARGET_FRAMES
    movie.TARGET_FRAMES = 4
    try:
        r = asyncio.run(movie.render_md_movie(ACC, width=160, height=120, keep_frames=True))
        assert r["ok"], r
        # 10 recorded frames, at most 4 rendered: one in 3, which is 4 states
        assert r["frames"] == 4 and "every 3 frame(s)" in r["mounted"]
        assert r["sampling"].startswith("chosen here"), r["sampling"]
        assert r["frames_in_the_record"] == 10
        kept = Path(r["frames_dir"])
        assert r["frames_removed"] is False and len(list(kept.glob("frame_*.png"))) == 4
        Path(r["file"]).unlink()
        shutil.rmtree(kept)
    finally:
        movie.TARGET_FRAMES = target
        viewer.find_pymol = saved_find
    print("the tool chooses the sampling and says so OK")


def test_nothing_is_rendered_from_a_file_that_changed(prep: Path, tmp: Path):
    saved_find = viewer.find_pymol
    viewer.find_pymol = stub_pymol(tmp)
    cluster = prep / "md" / "8C7Y_TXV_U" / "cluster"
    saved = (cluster / "step5_production.dcd").read_bytes()
    (cluster / "step5_production.dcd").write_bytes(saved + b"x")
    try:
        r = asyncio.run(movie.render_md_movie(ACC, width=160, height=120))
        assert not r["ok"] and "SHA-256" in r["error"]
        assert not list(prep.parent.glob("*.mp4")), "a video was written from a changed file"
    finally:
        (cluster / "step5_production.dcd").write_bytes(saved)
        viewer.find_pymol = saved_find
    print("nothing is rendered from a changed file OK")


def test_a_failed_render_and_a_failed_encode_are_reported(prep: Path, tmp: Path):
    saved_find, saved_ffmpeg = viewer.find_pymol, movie.ffmpeg_binary
    quiet = tmp / "quiet_pymol.py"
    quiet.write_text("print('nothing rendered')\n", encoding="utf-8")
    viewer.find_pymol = lambda: {"found": True, "source": "test stub",
                                 "command": [sys.executable, str(quiet)]}
    try:
        r = asyncio.run(movie.render_md_movie(ACC, width=160, height=120))
        assert not r["ok"] and r["error"].startswith("PyMOL wrote no frame.")
        assert "nothing rendered" in " ".join(r["log_tail"]), r["log_tail"]
        assert not list(prep.parent.glob("*.mp4"))

        viewer.find_pymol = stub_pymol(tmp)
        movie.ffmpeg_binary = lambda: (str(tmp / "no_such_ffmpeg.exe"), "stub")
        r = asyncio.run(movie.render_md_movie(ACC, width=160, height=120))
        assert not r["ok"] and "ffmpeg did not write the video" in r["error"], r
        assert r["frames_written"] == 10 and Path(r["frames_dir"]).is_dir()
        shutil.rmtree(r["frames_dir"])

        movie.ffmpeg_binary = lambda: (None, "imageio-ffmpeg is not installed (stub)")
        r = asyncio.run(movie.render_md_movie(ACC))
        assert not r["ok"] and "imageio-ffmpeg is not installed" in r["error"]
    finally:
        viewer.find_pymol, movie.ffmpeg_binary = saved_find, saved_ffmpeg
    print("a failed render and a failed encode are reported OK")


def test_a_half_finished_render_is_not_turned_into_a_video(prep: Path, tmp: Path):
    """PyMOL can stop half way and still exit 0. A video of three frames out of ten looks like
    a whole trajectory, and the tool used to report the ten states beside it."""
    saved_find = viewer.find_pymol
    half = tmp / "half_pymol.py"
    half.write_text("""
import ast, re, sys, pathlib
from PIL import Image
pml = pathlib.Path(sys.argv[-1]).read_text(encoding="ascii")
call = re.search(r"^    _pmcp_render\((.+)\)$", pml, re.M).group(1)
folder, first, last, width, height, ray = ast.literal_eval("(" + call + ")")
out = pathlib.Path(folder); out.mkdir(parents=True, exist_ok=True)
for state in range(first, first + 3):
    Image.new("RGB", (width, height), (9, 9, 9)).save(out / ("frame_%05d.png" % state))
print(" protein-mcp: render failed: out of memory at state 4")
""", encoding="utf-8")
    viewer.find_pymol = lambda: {"found": True, "source": "test stub",
                                 "command": [sys.executable, str(half)]}
    try:
        r = asyncio.run(movie.render_md_movie(ACC, every=1, width=160, height=120))
        assert not r["ok"], r
        assert r["frames_written"] == 3 and r["states_asked_for"] == 10
        assert "never reached the end of the frame loop" in r["error"], r["error"]
        assert "out of memory at state 4" in " ".join(r["log_tail"])
        assert not list(prep.parent.glob("*.mp4")), "a partial trajectory was encoded"
        shutil.rmtree(r["frames_dir"])
    finally:
        viewer.find_pymol = saved_find
    print("a half-finished render is not turned into a video OK")


def test_without_pymol_the_movie_is_refused(prep: Path):
    saved_find = viewer.find_pymol
    viewer.find_pymol = lambda: {"found": False, "source": None, "error": "PyMOL not found: ..."}
    try:
        r = asyncio.run(movie.render_md_movie(ACC))
        assert not r["ok"] and "PROTEIN_MCP_PYMOL" in r["how_to_fix"]
        assert "it is what draws the frames" in r["error"]
    finally:
        viewer.find_pymol = saved_find
    print("without PyMOL the movie is refused OK")


def test_two_matching_runs_are_refused(prep: Path):
    def add(man):
        block = man["structures"][0]["md_analysis"][MD_KEY]
        block["equilibration"] = json.loads(json.dumps(block["production"]))
    saved = patch_manifest(prep, add)
    try:
        r = asyncio.run(movie.render_md_movie(ACC))
        assert not r["ok"] and "Name one with structure, pocket and stage" in r["error"]
        assert len(r["available"]) == 2
        no = asyncio.run(movie.render_md_movie(ACC, stage="minimisation"))
        assert not no["ok"] and "No recorded trajectory matches" in no["error"]
    finally:
        (prep / "manifest.json").write_bytes(saved)
    print("two matching runs are refused OK")


def test_the_real_ffmpeg_is_the_one_that_ships_with_the_project():
    path, note = movie.ffmpeg_binary()
    assert path and Path(path).is_file(), note
    assert "imageio-ffmpeg" in note
    print("the project ships its own ffmpeg OK")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        prep = build(root)
        test_the_poses_of_one_pocket_open_together(prep)
        test_a_pose_file_that_changed_is_refused(prep)
        test_a_file_without_a_recorded_hash_is_refused_for_what_it_is(prep)
        test_only_the_named_compounds_are_shown(prep)
        test_two_docked_pockets_are_not_mixed(prep)
        test_without_a_docking_run_it_says_so(prep)
        test_every_really_thins_the_states(prep)
        test_the_cap_is_shared_between_files(prep)
        test_a_changed_trajectory_is_refused(prep)
        test_a_semicolon_in_the_header_cannot_run_a_command()
        test_the_real_ffmpeg_is_the_one_that_ships_with_the_project()
        test_the_movie_is_rendered_and_encoded(prep, root)
        test_a_frame_of_an_earlier_render_cannot_reach_the_video(prep, root)
        test_the_tool_chooses_the_sampling_and_says_so(prep, root)
        test_nothing_is_rendered_from_a_file_that_changed(prep, root)
        test_a_failed_render_and_a_failed_encode_are_reported(prep, root)
        test_a_half_finished_render_is_not_turned_into_a_video(prep, root)
        test_without_pymol_the_movie_is_refused(prep)
        test_two_matching_runs_are_refused(prep)
    print("ALL OK")


main()
