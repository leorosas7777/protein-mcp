"""Offline test for the compound key of build_md_job and analyze_md (step 6b, 2026-11-19).

Until today both tools keyed everything on the pocket alone (`032 C`), while build_md_package,
check_charmm_gui_output and analyze_admet had already been given a `compound`. Two compounds
simulated in the SAME pocket therefore collided: the second build_md_job overwrote the first
pocket's job block and wrote into the same `run/` folder, and the second analyze_md overwrote the
first md_analysis record. That is the collision `dock.run_id_of` was given a slug for; this file
is the control that it cannot come back.

The fixture is the real dossier of test_md_record (the 3OG7 package, build_md_job on an archive in
the real CHARMM-GUI layout, and the synthetic trajectory whose PSF is inside that archive). The
two "compounds" are given that same system on purpose: what is under test is the KEYING, not the
chemistry, so both archives must pass the same checks and the only thing that may tell the records
apart is the compound.

Controls that must fail if the keying goes away (checked by hand with tests/mutate_compound.py):
- build_md_job writing both compounds under one key or into one folder;
- analyze_md recording both compounds under one key;
- analyze_md taking the folder of one compound as the folder of another;
- a compound that is not in the recorded selection being accepted.

Run: uv run python tests/test_md_compound.py
"""
import json
import shutil
import sys
from pathlib import Path

from protein_mcp.files import choice, ligand_choice, md_analyze, md_run

sys.path.insert(0, str(Path(__file__).parent))
from test_choice import ACC  # noqa: E402
from test_ligand_choice import SORA, VEM  # noqa: E402
from test_md_package import POCKET  # noqa: E402
from test_md_record import manifest, setup  # noqa: E402

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
# three compounds: two get a package and a job, the third stays without one on purpose, so the
# refusal that lists the jobs of a pocket has something to refuse
SEL = [{"id": "CHEMBL1229517", "name": "VEMURAFENIB", "level": "T1", "smiles": VEM},
       {"id": "CHEMBL1336", "name": "SORAFENIB", "level": "T3", "smiles": SORA},
       {"id": "CHEMBL25", "name": "ASPIRIN", "level": "T3", "smiles": ASPIRIN}]
SLUGS = {"CHEMBL1229517": "chembl1229517", "CHEMBL1336": "chembl1336", "CHEMBL25": "chembl25"}


def jobs(m=None) -> dict:
    s = next(x for x in (m or manifest())["structures"] if x["pdb_id"] == "3OG7")
    return s.get("md_job") or {}


def analyses(m=None) -> dict:
    s = next(x for x in (m or manifest())["structures"] if x["pdb_id"] == "3OG7")
    return s.get("md_analysis") or {}


def fixture():
    """(prep, archive, {compound id: its cluster folder}). The crystal-ligand job of test_md_record
    is left in place: a per-compound job must not disturb it."""
    from test_md_run import SITE_REPORT

    prep, cluster, archive = setup()
    rec = ligand_choice.select_ligands(ACC, SEL, reason="the researcher's three for this test")
    assert rec["ok"], rec
    base = prep / "md" / "3OG7_032_C"
    folders = {}
    for cid in ("CHEMBL1229517", "CHEMBL1336"):
        pkg = prep / "md" / f"3OG7_032_C_{SLUGS[cid]}"
        if pkg.is_dir():
            shutil.rmtree(pkg)
        # the compound's own package folder, with the delivered files and package.json of the
        # pocket: the archive must keep passing _check_job, so the system has to be the same one
        shutil.copytree(base, pkg, ignore=shutil.ignore_patterns("run", "run.*", "cluster",
                                                                 "analysis_*.json"))
        shutil.copytree(cluster, pkg / "cluster")
        folders[cid] = pkg / "cluster"
    built = {cid: md_run.build_md_job(ACC, str(archive), "3OG7", POCKET,
                                      site_report=SITE_REPORT, compound=cid)
             for cid in folders}
    for cid, r in built.items():
        assert r["ok"], (cid, r)
    return prep, archive, folders, built


def test_each_compound_gets_its_own_job_key_and_folder():
    """Control that must fail: build_md_job keying on the pocket alone."""
    prep, _archive, folders, built = fixture()
    keys = set(jobs())
    assert keys == {POCKET, f"{POCKET} CHEMBL1229517", f"{POCKET} CHEMBL1336"}, keys
    for cid in folders:
        key = f"{POCKET} {cid}"
        assert built[cid]["recorded_under"] == f"structures[3OG7].md_job[{key}]", built[cid]
        block = jobs()[key]
        assert block["compound"]["id"] == cid, block["compound"]
        assert block["compound"]["slug"] == SLUGS[cid], block["compound"]
        assert block["dir"] == f"md/3OG7_032_C_{SLUGS[cid]}/run", block["dir"]
        assert Path(built[cid]["run_dir"]) == prep / block["dir"], built[cid]["run_dir"]
        assert cid in built[cid]["built"] and "crystal ligand" not in built[cid]["built"], built[cid]["built"]
    # the two jobs are two folders, and the crystal-ligand job of the fixture is untouched
    dirs = {jobs()[k]["dir"] for k in jobs()}
    assert len(dirs) == 3, dirs
    assert jobs()[POCKET]["dir"] == "md/3OG7_032_C/run", jobs()[POCKET]
    assert jobs()[POCKET].get("compound") in (None, {}), jobs()[POCKET].get("compound")
    # every job is in run_all.slurm under its own label, so submitting it cannot run one twice
    text = (prep / "md" / "run_all.slurm").read_text(encoding="utf-8")
    for cid in folders:
        assert f"3OG7_032_C_{SLUGS[cid]}/run/system/namd" in text, text[:400]


def test_each_compound_gets_its_own_analysis_record():
    """Control that must fail: analyze_md keying the md_analysis block on the pocket alone."""
    prep, _archive, folders, _built = fixture()
    out = {}
    for cid, cluster in folders.items():
        out[cid] = md_analyze.analyze_md(cluster, include_series=False, uniprot_id=ACC,
                                         structure="3OG7", pocket=POCKET, compound=cid)
        assert out[cid]["ok"], (cid, out[cid].get("error"))
        assert out[cid]["recording"]["recorded"], (cid, out[cid]["recording"])
    keys = set(analyses())
    assert keys == {f"{POCKET} CHEMBL1229517", f"{POCKET} CHEMBL1336"}, keys
    for cid in folders:
        key = f"{POCKET} {cid}"
        rec = analyses()[key]["production"]
        assert rec["compound"]["id"] == cid, rec["compound"]
        assert rec["pocket"] == POCKET, rec["pocket"]
        assert out[cid]["recording"]["key"] == f"structures[3OG7].md_analysis[{key}][production]", \
            out[cid]["recording"]["key"]
        # the series file is the compound's own, inside the compound's folder
        series = rec["series_file"]["path"]
        assert series == f"md/3OG7_032_C_{SLUGS[cid]}/analysis_production.json", series
        assert (prep / series).is_file(), series
        head = json.loads((prep / series).read_text(encoding="utf-8"))
        assert head["pocket"] == key, head["pocket"]
        # and every file it recorded is that compound's, not the other's
        for traj in rec["trajectories"]:
            assert traj["file"].startswith(f"md/3OG7_032_C_{SLUGS[cid]}/"), traj["file"]
        assert rec["folder"].startswith(f"md/3OG7_032_C_{SLUGS[cid]}/"), rec["folder"]
        assert rec["topology"]["file"].startswith(f"md/3OG7_032_C_{SLUGS[cid]}/"), rec["topology"]
    # Two records, two series files on disk. This does NOT prove the keying: the series file
    # carries its own generated_at and folder paths, so the two hashes would differ even if both
    # records shared one key. What distinguishes them is the `pocket`/`compound` each record
    # names, asserted above (reviewer H11, 2026-11-22): keep this as a file-identity check only.
    paths = {analyses()[f"{POCKET} {cid}"]["production"]["series_file"]["path"] for cid in folders}
    both = {analyses()[f"{POCKET} {cid}"]["production"]["series_file"]["sha256"] for cid in folders}
    assert len(paths) == 2 and len(both) == 2, (paths, both)
    for cid in folders:
        rec = analyses()[f"{POCKET} {cid}"]["production"]
        assert rec["compound"]["id"] == cid, rec["compound"]
        assert rec["series_file"]["path"].startswith(f"md/3OG7_032_C_{SLUGS[cid]}/"), rec["series_file"]
    assert POCKET not in analyses(), "the crystal ligand was not analysed and must have no record"


def test_the_folder_of_another_compound_is_refused():
    """Control that must fail: resolving the package folder without the slug. Passing compound A
    with compound B's folder is how one run would be recorded as the other."""
    _prep, _archive, folders, _built = fixture()
    r = md_analyze.analyze_md(folders["CHEMBL1336"], include_series=False, uniprot_id=ACC,
                              structure="3OG7", pocket=POCKET, compound="CHEMBL1229517")
    assert not r["ok"], r
    assert "is not inside" in r["error"] and "chembl1229517" in r["expected_under"], r
    assert not analyses(), analyses()
    # and the same folder without any compound: it is not inside the crystal ligand's folder either
    r = md_analyze.analyze_md(folders["CHEMBL1336"], include_series=False, uniprot_id=ACC,
                              structure="3OG7", pocket=POCKET)
    assert not r["ok"] and "is not inside" in r["error"], r
    assert r["expected_under"].endswith("3OG7_032_C"), r["expected_under"]
    assert not analyses(), analyses()


def test_a_compound_without_a_job_is_refused_with_the_jobs_there_are():
    """Control that must fail: a missing md_job being filled from another compound's block, or the
    refusal naming the pocket's jobs without saying which compound each belongs to."""
    _prep, _archive, folders, _built = fixture()
    # CHEMBL25 is in the recorded selection and has no package and no job
    r = md_analyze.analyze_md(folders["CHEMBL1229517"], include_series=False, uniprot_id=ACC,
                              structure="3OG7", pocket=POCKET, compound="CHEMBL25")
    assert not r["ok"], r
    assert f"no md_job is recorded for 3OG7 {POCKET} CHEMBL25" in r["error"], r["error"]
    assert set(r["md_jobs_recorded_for_this_pocket"]) == set(jobs()), r
    assert not analyses(), analyses()
    # a compound that was never selected is refused by the selection, not by the folder
    r = md_analyze.analyze_md(folders["CHEMBL1229517"], include_series=False, uniprot_id=ACC,
                              structure="3OG7", pocket=POCKET, compound="CHEMBL999999")
    assert not r["ok"] and "not in the recorded selection" in r["error"], r
    assert not analyses(), analyses()
    # and compound alone, without uniprot_id, says what it is for instead of being ignored
    r = md_analyze.analyze_md(folders["CHEMBL1229517"], include_series=False,
                              compound="CHEMBL1229517")
    assert not r["ok"] and "needs uniprot_id" in r["error"], r


def test_build_md_job_without_a_compound_names_the_compound_packages():
    """Control that must fail: falling back to any package found, which is how one compound's
    archive would be verified against another compound's system."""
    prep, archive, _folders, _built = fixture()
    from test_md_run import SITE_REPORT

    pkg = prep / "md" / "3OG7_032_C" / "package.json"
    aside = pkg.with_suffix(".json.aside")
    pkg.rename(aside)
    try:
        r = md_run.build_md_job(ACC, str(archive), "3OG7", POCKET, site_report=SITE_REPORT)
        assert not r["ok"], r
        assert "No package at" in r["error"] and "name one with `compound`" in r["error"], r["error"]
        assert set(r["compound_packages"]) == {"3OG7_032_C_chembl1229517", "3OG7_032_C_chembl1336"}, r
        # nothing was written: the crystal-ligand job block is as it was
        assert set(jobs()) == {POCKET, f"{POCKET} CHEMBL1229517", f"{POCKET} CHEMBL1336"}, set(jobs())
    finally:
        aside.rename(pkg)


def test_a_second_run_of_one_compound_replaces_only_that_compounds_record():
    """Control that must fail: a re-analysis writing over every key of the pocket."""
    _prep, _archive, folders, _built = fixture()
    for cid, cluster in folders.items():
        assert md_analyze.analyze_md(cluster, include_series=False, uniprot_id=ACC,
                                     structure="3OG7", pocket=POCKET,
                                     compound=cid)["recording"]["recorded"]
    first = analyses()[f"{POCKET} CHEMBL1336"]["production"]["generated_at"]
    again = md_analyze.analyze_md(folders["CHEMBL1229517"], include_series=False, uniprot_id=ACC,
                                  structure="3OG7", pocket=POCKET, compound="CHEMBL1229517",
                                  pocket_cutoff=6.0)
    assert again["recording"]["recorded"], again["recording"]
    kept = analyses()[f"{POCKET} CHEMBL1336"]["production"]
    assert kept["generated_at"] == first, "the other compound's record was rewritten"
    assert kept["arguments"]["pocket_cutoff"] == 5.0, kept["arguments"]
    mine = analyses()[f"{POCKET} CHEMBL1229517"]["production"]
    assert mine["arguments"]["pocket_cutoff"] == 6.0, mine["arguments"]
    assert mine["arguments"]["compound"] == "CHEMBL1229517", mine["arguments"]


def main():
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as exc:
                failed += 1
                print(f"FAIL {name}: {exc}")
    print(("All passed" if not failed else f"{failed} failed"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
