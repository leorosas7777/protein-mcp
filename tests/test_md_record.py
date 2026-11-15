"""Offline test for the recording half of analyze_md (Phase 3 step 9, decided 2026-11-12).

With `uniprot_id`, analyze_md writes one md_analysis block per stage into prep/manifest.json and
the per-frame series into a file of their own. The dossier is the real one of test_md_run (the
3OG7 package of test_md_package, then build_md_job on an archive in the real layout), and the
trajectory is the synthetic one of test_md_analyze, whose PSF is put INTO that archive so the
registered package and the trajectory really share a topology.

Controls that must fail, one per check the user asked to refuse rather than warn:
- a PSF one byte different from the registered one is refused and nothing is written;
- a folder outside prep/md/<pocket folder>/ is refused;
- a DCD that changes while the analysis runs, and a manifest rewritten meanwhile, are not
  recorded;
- with the PSF hash absent from md_job (a job recorded before the field existed) the archive is
  the reference, and with the archive gone the recording is refused.
And the one that must stay as it was: without uniprot_id, nothing is written at all.

Run: uv run python tests/test_md_record.py
"""
import hashlib
import json
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from protein_mcp.files import choice, md_analyze, md_run  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from test_choice import ACC  # noqa: E402
from test_md_analyze import write_fixture, write_psf  # noqa: E402
from test_md_package import POCKET, test_build_writes_what_charmm_gui_needs  # noqa: E402
from test_md_run import SITE_REPORT, archive_with_namd_scheme  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def setup():
    """A dossier with a recorded job whose archive holds the fixture's PSF, and a trajectory
    folder under prep/md/3OG7_032_C/cluster with a production (step5) and an equilibration
    (step4) stage."""
    prep, delivered = test_build_writes_what_charmm_gui_needs()
    package = json.loads((delivered.parent / "package.json").read_text(encoding="utf-8"))
    tmp = Path(tempfile.mkdtemp())
    psf_text = write_psf(tmp / "fixture.psf").read_text()
    archive = archive_with_namd_scheme(tmp, delivered, package["to_model"], psf=psf_text)
    built = md_run.build_md_job(ACC, str(archive), "3OG7", POCKET, site_report=SITE_REPORT)
    assert built["ok"], built
    cluster = prep / "md" / "3OG7_032_C" / "cluster"
    write_fixture(cluster)                                    # step5_1.dcd + step5_production.inp
    write_fixture(cluster, name="step4_equilibration.dcd", inp=False, with_minimization=False,
                  shifts=[0.0, 0.5, 1.0, 1.5])
    # the trajectory's topology is a byte copy of the one build_md_job extracted
    shutil.copy(Path(built["run_dir"]) / "system" / "namd" / "step3_input.psf",
                cluster / "step3_input.psf")
    return prep, cluster, archive


def manifest():
    return json.loads(choice.manifest_path(ACC).read_text(encoding="utf-8"))


def block(m=None):
    s = next(x for x in (m or manifest())["structures"] if x["pdb_id"] == "3OG7")
    return s.get("md_analysis", {}).get(POCKET, {})


def run(cluster, **kw):
    return md_analyze.analyze_md(cluster, include_series=False, uniprot_id=ACC,
                                 structure="3OG7", pocket=POCKET, **kw)


def _series_anywhere(node, frames=(4, 6)) -> bool:
    """A per-frame series under ANY key: `values`/`times_ps`, or any list of numbers as long as
    the analysed or total frame count (that is how `frame_indices` was caught by the reviewer)."""
    if isinstance(node, dict):
        return any(k in ("values", "times_ps") or _series_anywhere(v, frames)
                   for k, v in node.items())
    if isinstance(node, list):
        if len(node) in frames and all(isinstance(x, (int, float)) for x in node):
            return True
        return any(_series_anywhere(v, frames) for v in node)
    return False


def test_production_and_equilibration_are_recorded_apart():
    prep, cluster, _archive = setup()
    r = run(cluster)
    assert r["ok"], r
    rec = r["recording"]
    assert rec["recorded"] and rec["stage_label"] == "production", rec
    assert not _series_anywhere(r), "the chat reply must carry no series"
    b = block()
    prod = b["production"]
    assert set(b) == {"production"}, set(b)
    assert prod["stage"] == "step5" and "md_job.scheme" in prod["stage_label_source"]
    # every DCD with its size, SHA-256 and frame count, computed here independently
    dcd = cluster / "step5_1.dcd"
    assert prod["trajectories"] == [{"file": "md/3OG7_032_C/cluster/step5_1.dcd",
                                     "size_bytes": dcd.stat().st_size, "sha256": sha(dcd),
                                     "frames": 6}], prod["trajectories"]
    assert prod["topology"]["sha256"] == sha(cluster / "step3_input.psf")
    assert prod["topology"]["reference"]["source"].startswith("md_job.psf")
    # the series are in their own file, complete, and the manifest holds its hash, not them
    series = prep / prod["series_file"]["path"]
    assert series.name == "analysis_production.json"
    assert sha(series) == prod["series_file"]["sha256"]
    full = json.loads(series.read_text(encoding="utf-8"))["result"]
    # the fixture's analytic answer: a rigid 0-1-2-3 A translation of the ligand
    assert [round(v, 2) for v in full["ligand_rmsd"]["values"]] == [0.0, 1.0, 2.0, 3.0], full
    assert "frame_indices" in full["hbonds"]["count_per_frame"], "the file keeps every series"
    assert len(full["ligand_rmsd"]["values"]) == 4 and len(full["com_distance"]["times_ps"]) == 4
    assert not _series_anywhere(prod["result"]), "series must not be in the manifest"
    assert prod["result"]["reading"]["equilibration_control"]["status"] in (
        "MET", "NOT MET", "NO EVIDENCE", "CANNOT BE ASSESSED")

    e = run(cluster, stage="step4")
    assert e["ok"] and e["recording"]["stage_label"] == "equilibration", e.get("recording")
    b = block()
    assert set(b) == {"production", "equilibration"}, set(b)
    assert b["production"] == prod, "recording the equilibration must not touch production"
    assert b["equilibration"]["trajectories"][0]["frames"] == 4
    assert (cluster.parent / "analysis_equilibration.json").is_file()

    # repeating a stage replaces that stage only
    eq_before = b["equilibration"]
    m = manifest()
    s = next(x for x in m["structures"] if x["pdb_id"] == "3OG7")
    s["md_analysis"][POCKET]["production"]["sentinel"] = "stale"
    choice.write_manifest(choice.manifest_path(ACC), m)
    again = run(cluster)
    assert again["recording"]["recorded"], again["recording"]
    b = block()
    assert b["equilibration"] == eq_before
    assert "sentinel" not in b["production"], "the production record was not replaced"
    # the stage inputs and the restraint file are tied to the record by hash
    inp = cluster / "step5_production.inp"
    assert {"file": "md/3OG7_032_C/cluster/step5_production.inp", "sha256": sha(inp)} in \
        b["production"]["inputs_read"], b["production"]["inputs_read"]
    print("ok  production and equilibration recorded as two runs; a repeat replaces only its own")
    return prep, cluster


def test_a_different_psf_is_refused_and_nothing_is_written():
    """CONTROL: one byte more in the topology and the recording is refused, not warned about."""
    prep, cluster, _archive = setup()
    before = choice.manifest_path(ACC).read_bytes()
    psf = cluster / "step3_input.psf"
    psf.write_bytes(psf.read_bytes() + b"\n")
    r = run(cluster)
    assert r["ok"] is False, r
    assert "SHA-256" in r["error"] and "registered package" in r["error"], r["error"]
    assert r["psf_sha256"] != r["registered_psf_sha256"]
    assert choice.manifest_path(ACC).read_bytes() == before, "the manifest must be untouched"
    assert not list(cluster.parent.glob("analysis_*.json")), "no series file either"
    # and the same trajectory measures fine when nothing is to be recorded
    loose = md_analyze.analyze_md(cluster, include_series=False)
    assert loose["ok"] and loose["recording"]["recorded"] is False, loose.get("recording")
    assert choice.manifest_path(ACC).read_bytes() == before
    print("ok  a PSF one byte different is refused and nothing is written")


def test_a_folder_outside_the_pocket_folder_is_refused():
    prep, cluster, _archive = setup()
    away = Path(tempfile.mkdtemp()) / "cluster"
    shutil.copytree(cluster, away)
    before = choice.manifest_path(ACC).read_bytes()
    r = run(away)
    assert r["ok"] is False and "is not inside" in r["error"], r
    assert r["expected_under"].endswith("3OG7_032_C"), r
    assert choice.manifest_path(ACC).read_bytes() == before
    print("ok  a trajectory folder outside prep/md/3OG7_032_C is refused")


def test_without_uniprot_id_nothing_is_written():
    """Condition 2: the loose use is as it was, and writes nothing anywhere."""
    prep, cluster, _archive = setup()
    before = choice.manifest_path(ACC).read_bytes()
    listing = sorted(p.name for p in cluster.parent.rglob("*"))
    r = md_analyze.analyze_md(cluster, include_series=False)
    assert r["ok"] and r["recording"]["recorded"] is False, r.get("recording")
    assert choice.manifest_path(ACC).read_bytes() == before
    assert sorted(p.name for p in cluster.parent.rglob("*")) == listing
    bad = md_analyze.analyze_md(cluster, structure="3OG7", pocket=POCKET)
    assert bad["ok"] is False and "uniprot_id" in bad["error"], bad
    print("ok  without uniprot_id nothing is written; structure/pocket alone are refused")


def test_a_stage_the_job_does_not_name_is_refused():
    prep, cluster, _archive = setup()
    write_fixture(cluster, name="step6_extra.dcd", inp=False)
    r = run(cluster, stage="step6")
    assert r["ok"] is False and "neither the production stage" in r["error"], r
    assert "step6" not in json.dumps(block())
    print("ok  a stage the recorded job does not name is refused")


def test_changes_during_the_analysis_are_not_recorded():
    """CONTROL: the checks after the analysis are real. A DCD that changes, and a manifest that
    prepare_structures rewrote, leave the manifest without the block."""
    prep, cluster, _archive = setup()
    real = md_analyze._analyze

    def touch_dcd(*a, **kw):
        out = real(*a, **kw)
        dcd = cluster / "step5_1.dcd"
        dcd.write_bytes(dcd.read_bytes() + b"\0")
        return out

    md_analyze._analyze = touch_dcd
    try:
        r = run(cluster)
    finally:
        md_analyze._analyze = real
    assert r["ok"] and r["recording"]["recorded"] is False, r["recording"]
    assert "step5_1.dcd changed" in r["recording"]["why"], r["recording"]
    assert block() == {}

    prep, cluster, _archive = setup()

    def rewrite_manifest(*a, **kw):
        out = real(*a, **kw)
        m = manifest()
        m["created_at"] = "2099-01-01T00:00:00Z"
        choice.write_manifest(choice.manifest_path(ACC), m)
        return out

    md_analyze._analyze = rewrite_manifest
    try:
        r = run(cluster)
    finally:
        md_analyze._analyze = real
    assert r["recording"]["recorded"] is False and "rewrote the manifest" in r["recording"]["why"]
    assert block() == {}
    print("ok  a DCD changed and a manifest rewritten during the analysis are not recorded")


def test_an_older_job_falls_back_to_the_archive_and_refuses_without_it():
    prep, cluster, archive = setup()
    m = manifest()
    s = next(x for x in m["structures"] if x["pdb_id"] == "3OG7")
    s["md_job"][POCKET].pop("psf")               # as recorded before 2026-11-12
    choice.write_manifest(choice.manifest_path(ACC), m)
    r = run(cluster)
    assert r["recording"]["recorded"], r
    assert "inside job.tgz" in block()["production"]["topology"]["reference"]["source"]
    archive.rename(archive.with_suffix(".moved"))
    r = run(cluster, stage="step4")
    assert r["ok"] is False and "no longer at" in r["error"], r
    assert "equilibration" not in block()
    print("ok  without md_job.psf the archive is the reference; without the archive, refused")


def test_the_equilibration_control_is_a_field():
    """The status the report heads section 2 with, against a held and an unheld backbone."""
    def result(drift, held=False):
        return {
            "protein_rmsd": {"drift": drift, "segid": "PROA", "first_quarter_mean": 1.0,
                             "last_quarter_mean": 1.0 + drift, "per_segment": {}},
            "ligand_rmsd": {"max": 1.0, "last_quarter_mean": 1.0},
            "com_distance": {"range": 0.5, "mean": 3.0, "n": 10},
            "hbonds": {"pairs": [], "pairs_found": 0},
            "trajectory": {"restraints": {"held": {"protein_backbone": {
                "known": True, "held": 100 if held else 0, "note": "100 of 100 held"}}}},
        }
    get = lambda r: md_analyze.reading(r)["equilibration_control"]  # noqa: E731
    assert get(result(1.35))["status"] == "NOT MET"
    assert get(result(0.3))["status"] == "MET"
    assert get(result(0.07, held=True))["status"] == "NO EVIDENCE"
    # control: a held backbone that still drifts beyond the threshold is not hidden as "no evidence"
    assert get(result(1.0, held=True))["status"] == "NOT MET"
    c = get(result(1.35))
    assert (c["drift_angstrom"], c["threshold_angstrom"], c["segid"]) == (1.35, 0.5, "PROA")
    print("ok  the control status is a field: NOT MET / MET / NO EVIDENCE, held or not")


# ---------------------------------------------------------------------------------------------
# Section 2 of build_experiment_report, read from what was recorded above
# ---------------------------------------------------------------------------------------------

def _sec2(md: str) -> str:
    return md.split("## 2. Molecular dynamics")[1].split("## 3. ADMET")[0]


def _plotted(doc, title_part):
    """The (x, y) points of the panel whose title contains `title_part`, from the first figure."""
    from reportlab.graphics.charts.lineplots import LinePlot
    from reportlab.graphics.shapes import String
    fig = next(b for b in doc.blocks if b[0] == "figure")
    drawing = fig[1](500.0)
    items = drawing.contents
    titles = [i for i, x in enumerate(items) if isinstance(x, String) and title_part in x.text]
    plot = next(x for x in items[titles[0]:] if isinstance(x, LinePlot))
    return plot.data


def test_the_report_heads_each_run_with_its_control_and_draws_from_the_file():
    from protein_mcp import experiment_report, report
    prep, cluster = test_production_and_equilibration_are_recorded_apart()
    r = experiment_report.build_experiment_report(ACC)
    assert r["ok"], r
    md = Path(r["markdown_path"]).read_text(encoding="utf-8")
    sec = _sec2(md)
    inv = {s["number"]: s for s in r["sections"]}
    assert inv["2"]["status"] == "has content", inv["2"]
    assert not any("molecular dynamics" in x for x in r["sections_not_run"]), r["sections_not_run"]
    assert len(r["md_runs"]) == 2 and r["md_runs"][0].startswith("3OG7 032 C production"), r["md_runs"]

    # 5.1: the configuration table, quoted from md_package and md_job
    m = manifest()
    job = next(x for x in m["structures"] if x["pdb_id"] == "3OG7")["md_job"][POCKET]
    for fact in ("Starting pose", "of docking run 032_regenerated", job["job_id"],
                 job["psf"]["sha256"], "Production as written", "Equilibration as written",
                 "whether the cluster ran these inputs unchanged is not recorded"):
        assert fact in sec, fact

    # each run: its own heading, the control line FIRST, then the figures, all said to be its own
    prod = sec.split("#### Production run")[1].split("#### Equilibration run")[0]
    equi = sec.split("#### Equilibration run")[1]
    for block_, label in ((prod, "production"), (equi, "equilibration")):
        first = block_.split("\n\n", 2)[1]
        assert first.startswith("EQUILIBRATION CONTROL: "), first[:200]
        assert block_.index("EQUILIBRATION CONTROL") < block_.index("Ligand heavy-atom RMSD")
        assert f"comes from the {label} run only" in block_
        assert "PROPORTION OF ANALYSED FRAMES" in block_ and "never a fraction of time" in block_
        assert "Figure (drawn in the PDF only)" in block_
        assert "READING OFFERED, NOT MEASURED" in block_
        assert "None" not in block_.replace("Nonetheless", ""), "a missing value printed as None"
    rec = block()["production"]
    assert rec["trajectories"][0]["sha256"] in prod and rec["topology"]["sha256"] in prod
    assert "All 4 analysed frames plotted" in prod

    # every table width is a FRACTION of the frame (the defect of step 7: points cut the PDF)
    doc, *_ = experiment_report.build_document(ACC, report.load_manifest(ACC), "t")
    widths = [sum(b[3]) for b in doc.blocks if b[0] == "table"]
    assert all(w <= 1.0001 for w in widths), widths
    print("ok  section 2: configuration, then each run headed by its control, figures said to be its own")
    return prep, cluster


def test_the_plots_come_from_the_series_file_and_nothing_else():
    """CONTROL: edit the series file (and its recorded hash) and the plot follows the FILE while
    the manifest summary stays as it was; edit it without the hash and the plots are omitted."""
    from protein_mcp import experiment_report, report
    prep, cluster = test_the_report_heads_each_run_with_its_control_and_draws_from_the_file()
    doc, *_ = experiment_report.build_document(ACC, report.load_manifest(ACC), "t")
    before = [y for _x, y in _plotted(doc, "Ligand heavy-atom RMSD")[0]]
    assert max(before) > 2.9, before                     # the fixture's 0-1-2-3 A translation

    m = manifest()
    rec = next(x for x in m["structures"] if x["pdb_id"] == "3OG7")["md_analysis"][POCKET]["production"]
    series = prep / rec["series_file"]["path"]
    data = json.loads(series.read_text(encoding="utf-8"))
    data["result"]["ligand_rmsd"]["values"] = [0.0, 7.0, 7.0, 7.0]
    series.write_text(json.dumps(data), encoding="utf-8")
    rec["series_file"]["sha256"] = sha(series)
    choice.write_manifest(choice.manifest_path(ACC), m)
    doc, *_ = experiment_report.build_document(ACC, report.load_manifest(ACC), "t")
    after = [y for _x, y in _plotted(doc, "Ligand heavy-atom RMSD")[0]]
    assert after == [0.0, 7.0, 7.0, 7.0], after
    assert rec["result"]["ligand_rmsd"]["max"] != 7.0, "the summary was not the source of the plot"

    # the file changed and the manifest still holds the old hash: no plot, and the reason
    series.write_text(json.dumps(data) + " ", encoding="utf-8")
    r = experiment_report.build_experiment_report(ACC)
    prod = _sec2(Path(r["markdown_path"]).read_text(encoding="utf-8")).split(
        "#### Equilibration run")[0]
    assert "Plots omitted:" in prod and "no longer has the SHA-256" in prod, prod[-1500:]
    print("ok  the plot follows the series file; a file that no longer matches its hash is not drawn")


def test_a_record_from_another_version_reads_no_record():
    """Fields absent from the manifest read 'no record', never None, and a missing control field
    is said to be missing rather than re-derived."""
    from protein_mcp import experiment_report, report
    prep, cluster = test_the_report_heads_each_run_with_its_control_and_draws_from_the_file()
    m = manifest()
    s = next(x for x in m["structures"] if x["pdb_id"] == "3OG7")
    rec = s["md_analysis"][POCKET]["production"]
    rec["result"]["reading"].pop("equilibration_control")
    for k in ("series_file", "trajectories", "topology", "stage_label_source"):
        rec.pop(k)
    rec["result"]["trajectory"].pop("restraints")
    s["md_job"][POCKET].pop("psf")
    choice.write_manifest(choice.manifest_path(ACC), m)
    doc, *_ = experiment_report.build_document(ACC, report.load_manifest(ACC), "t")
    sec = _sec2(report.render_markdown(doc))
    prod = sec.split("#### Production run")[1].split("#### Equilibration run")[0]
    assert "EQUILIBRATION CONTROL: NOT RECORDED" in prod, prod[:600]
    assert "Plots omitted: no series file is recorded" in prod
    assert "not recorded by the version of build_md_job" in sec
    assert "None" not in sec, [l for l in sec.splitlines() if "None" in l][:3]
    # control: with the 'no record' helper switched off, the same manifest does print None
    real = experiment_report._v
    experiment_report._v = lambda x: x
    try:
        doc, *_ = experiment_report.build_document(ACC, report.load_manifest(ACC), "t")
    finally:
        experiment_report._v = real
    assert "None" in _sec2(report.render_markdown(doc))
    print("ok  an older record reads 'no record' and says the control was not recorded")


def test_the_checks_the_first_round_left_without_a_control():
    """Reviewer H6: each of these checks is removed in thought and the test must then fail."""
    import os
    # (1) the folder is inside, but the topology named is outside the pocket folder
    prep, cluster, archive = setup()
    away = Path(tempfile.mkdtemp()) / "step3_input.psf"
    shutil.copy(cluster / "step3_input.psf", away)
    r = run(cluster, topology=str(away))
    assert r["ok"] is False and "outside" in r["error"] and str(away) in r["error"], r
    # (2) a restraint file named by the input that lies OUTSIDE the pocket folder
    inp = cluster / "step5_production.inp"
    ref = Path(tempfile.mkdtemp()) / "posres.ref"
    ref.write_text("REMARK\nEND\n")
    text = inp.read_text()
    inp.write_text(text + f"constraints on\nconsref {os.path.relpath(ref, cluster)}\n")
    r = run(cluster)
    assert r["ok"] is False and "outside" in r["error"] and "posres.ref" in r["error"], r
    inp.write_text(text)
    # (3) an archive whose bytes changed, on the fallback path: refused
    m = manifest()
    s = next(x for x in m["structures"] if x["pdb_id"] == "3OG7")
    s["md_job"][POCKET].pop("psf")
    choice.write_manifest(choice.manifest_path(ACC), m)
    original = archive.read_bytes()
    archive.write_bytes(original + b"\0")
    r = run(cluster)
    assert r["ok"] is False and "no longer has the SHA-256" in r["error"], r
    archive.write_bytes(original)
    # (4) on the fallback path, a trajectory PSF that differs from the one in the archive
    psf = cluster / "step3_input.psf"
    psf.write_bytes(psf.read_bytes() + b"\n")
    r = run(cluster)
    assert r["ok"] is False and "inside job.tgz" in r["error"], r
    assert block() == {}
    print("ok  topology and restraint file outside, a changed archive and a PSF that differs "
          "from the one in the archive are all refused")


def test_a_restraint_keyword_without_the_file_is_not_called_held():
    """Reviewer H1: the input says `constraints on` but the atoms held are UNKNOWN. The control
    must not say the backbone was held; with a small drift it CANNOT BE ASSESSED."""
    r = {
        "protein_rmsd": {"drift": 0.3, "segid": "PROA", "first_quarter_mean": 1.0,
                         "last_quarter_mean": 1.3, "per_segment": {}},
        "ligand_rmsd": {"max": 1.0, "last_quarter_mean": 1.0},
        "com_distance": {"range": 0.5, "mean": 3.0, "n": 10},
        "hbonds": {"pairs": [], "pairs_found": 0},
        "trajectory": {"restraints": {"atoms_held": True, "held": {"protein_backbone": {
            "known": False, "note": "UNKNOWN: no restraint reference file was read"}}}},
    }
    c = md_analyze.reading(r)["equilibration_control"]
    assert c["status"] == "CANNOT BE ASSESSED", c
    assert c["backbone_held"] is None and "UNKNOWN" in c["why"], c
    assert "was held" not in c["why"]
    # the drift is reported both ways, so a negative one is not printed as positive silently
    r["protein_rmsd"]["drift"] = -0.9
    r["trajectory"] = {}
    c = md_analyze.reading(r)["equilibration_control"]
    assert (c["drift_angstrom"], c["drift_signed_angstrom"], c["status"]) == (0.9, -0.9, "NOT MET")
    r["protein_rmsd"]["drift_warning"] = "too few frames"
    assert md_analyze.reading(r)["equilibration_control"]["status"] == "CANNOT BE ASSESSED"
    print("ok  a restraint keyword alone is not 'held'; signed and absolute drift both recorded")


def test_a_hydrogen_bond_count_of_another_length_is_not_drawn():
    """Reviewer H7: the borrowed time axis is used only when the lengths agree."""
    from protein_mcp import experiment_report
    from reportlab.graphics.shapes import String
    full = {"ligand_rmsd": {"values": [0.0, 1.0, 2.0], "times_ps": [0.0, 10.0, 20.0]},
            "com_distance": {"values": [3.0, 3.0, 3.0], "times_ps": [0.0, 10.0, 20.0]},
            "protein_rmsd": {"segid": "PROA", "values": [0, 1, 1], "times_ps": [0.0, 10.0, 20.0]},
            "hbonds": {"count_per_frame": {"values": [1, 2]}}}
    texts = [x.text for x in experiment_report._md_figure(full, 1)(500.0).contents
             if isinstance(x, String)]
    assert "values recorded, but no time axis of the same length: not drawn" in texts, texts
    full["hbonds"]["count_per_frame"]["values"] = [1, 2, 2]
    texts = [x.text for x in experiment_report._md_figure(full, 1)(500.0).contents
             if isinstance(x, String)]
    assert not any("not drawn" in x for x in texts), texts
    print("ok  an H-bond count of another length is not put on the ligand's time axis")


def test_a_failed_manifest_write_keeps_the_previous_series_file():
    """Reviewer H8: the new series file replaces the old one only after the manifest is written."""
    prep, cluster, _archive = setup()
    assert run(cluster)["recording"]["recorded"]
    rec = block()["production"]
    series = prep / rec["series_file"]["path"]
    real = choice.write_manifest

    def full_disk(*a, **k):
        raise OSError(28, "No space left")

    choice.write_manifest = full_disk
    try:
        r = run(cluster)
    finally:
        choice.write_manifest = real
    assert r["recording"]["recorded"] is False and "No space left" in r["recording"]["why"], r
    assert sha(series) == rec["series_file"]["sha256"], "the old record's file was overwritten"
    assert not series.with_name(series.name + ".part").exists()
    print("ok  a failed manifest write leaves the previous series file matching its record")


def main():
    test_the_equilibration_control_is_a_field()
    test_a_restraint_keyword_without_the_file_is_not_called_held()
    test_a_hydrogen_bond_count_of_another_length_is_not_drawn()
    test_the_checks_the_first_round_left_without_a_control()
    test_a_failed_manifest_write_keeps_the_previous_series_file()
    test_the_plots_come_from_the_series_file_and_nothing_else()
    test_a_record_from_another_version_reads_no_record()
    test_production_and_equilibration_are_recorded_apart()
    test_a_different_psf_is_refused_and_nothing_is_written()
    test_a_folder_outside_the_pocket_folder_is_refused()
    test_without_uniprot_id_nothing_is_written()
    test_a_stage_the_job_does_not_name_is_refused()
    test_changes_during_the_analysis_are_not_recorded()
    test_an_older_job_falls_back_to_the_archive_and_refuses_without_it()
    print("ALL OK")


if __name__ == "__main__":
    main()
