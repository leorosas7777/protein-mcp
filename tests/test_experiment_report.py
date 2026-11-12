"""Offline test for build_experiment_report (protein_mcp/experiment_report.py, Phase 3 step 3).

The docking blocks are REAL: the manifest of test_choice (real mmCIFs cleaned with clean.py,
real RCSB chain replies, real ModelServer SDF of 3OG7 032 C) is docked with Meeko really
preparing the receptor and the Vina stand-in of test_dock, which returns the real Vina 1.2.7
output of 2026-11-05. The Meeko messages used for the deleted/refused branches are the real
outputs stored in tests/fixtures/dock/ (8C7Y TXV U and 5JSM 6NB J), parsed by dock.parse_meeko.

Checks that matter here (anti-hallucination):
  * every link in the record is in the manifest (and the checker is shown to catch one that is not);
  * the "affinity (kcal/mol)" header that the Vina log in the dossier really contains is NOT copied;
  * every mention of IC50/Ki/Kd is the caveat that they are not comparable, never a value;
  * the reason recorded by select_structure is quoted only on the run it belongs to;
  * the research report gains one line about the experiment record, with no figures in it.
Run: uv run python tests/test_experiment_report.py
"""
import copy
import json
import re
import sys
import tempfile
from pathlib import Path

from protein_mcp import dossier, experiment_report, report
from protein_mcp.files import choice, dock

sys.path.insert(0, str(Path(__file__).parent))
from test_dock import DFIX, calls, make_fake_vina, make_prep, run  # noqa: E402
from test_choice import ACC  # noqa: E402

# Vina's own table header in the log kept in the dossier ("affinity | ... | (kcal/mol)"): the
# report must not repeat it, or a calculated score would be labelled as a measured affinity.
VINA_LOG_HEADER = "mode |   affinity | dist from best mode"


def build(acc=ACC):
    return experiment_report.build_experiment_report(acc)


def md_of(result):
    return Path(result["markdown_path"]).read_text(encoding="utf-8")


def links_are_sourced(md: str, manifest_text: str):
    urls = re.findall(r"\]\((https?://[^)\s]+)\)", md)
    assert urls, "no link in the record: the check would pass on an empty document"
    for url in urls:
        assert url in manifest_text, f"link not backed by the manifest: {url}"


def test_no_runs():
    """No manifest -> clear error; a manifest without docking -> the section says 'not run'."""
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())
    md_path = dossier.DOSSIER_DIR / ACC / f"experiments_{ACC}.md"
    r = build()
    assert not r["ok"] and "Run prepare_structures first" in r["error"], r
    # no manifest -> nothing written (the same path is shown to exist after the next build)
    assert not md_path.exists() and not md_path.with_suffix(".pdf").exists()
    assert not build("not an id")["ok"]
    prep, _m = make_prep()
    r = build()
    assert r["ok"] and r["docking_runs"] == [] and r["pages"] >= 1, r
    # make_prep moves DOSSIER_DIR to a new folder: compare the path inside it
    assert (Path(r["markdown_path"]) == dossier.DOSSIER_DIR / md_path.relative_to(md_path.parents[1])
            and Path(r["markdown_path"]).exists()), "the 'nothing written' check used another path"
    md = md_of(r)
    assert "Not run. No docking run is recorded in prep/manifest.json" in md
    assert f'dock_ligand("{ACC}", structure, pocket)' in md
    assert "docking (no run in the manifest)" in r["sections_not_run"]
    # the two sections that no tool fills yet are present and say so, with the step that would
    assert "## 2. Molecular dynamics" in md and "Phase 3 step 4" in md
    # Section 3 exists from the start and declares the gap by naming the tool that fills it.
    # It must NOT call anything PREDICTED: nothing analyze_admet writes comes from a trained
    # model (ESTADO 7.3), and the placeholder used to promise "PREDICTED values".
    assert "## 3. ADMET" in md and "analyze_admet" in md
    assert "ADMET (no analyze_admet result in the manifest)" in r["sections_not_run"]
    assert "PREDICTED" not in md and "predicted" not in md
    # H5: the pending line does not contradict the two sections declared not run
    assert "Preparation steps still pending in the manifest:" in md
    assert "Molecular dynamics (section 2) has not been run" in md
    assert "nothing listed" not in md
    return prep


def test_real_run(prep: Path):
    """One real docking run (3OG7 032 C, regenerated, the pre-registered control)."""
    assert run("3OG7", "032 C")["ok"]
    manifest_text = choice.manifest_path(ACC).read_text(encoding="utf-8")
    man = json.loads(manifest_text)
    rec = next(s for s in man["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]["032_regenerated"]
    r = build()
    assert r["ok"] and r["docking_runs"] == [f"3OG7 032 C 032_regenerated ({rec['computed_at']})"], r
    assert "docking (no run in the manifest)" not in r["sections_not_run"]
    assert r["manifest_created_at"] == man["created_at"]
    md = md_of(r)

    # provenance, all of it from the manifest
    assert "1 docking run recorded, on 1 of the 5 pockets computed in the manifest (1 of 4 structures)" in md, md[:2000]
    assert ("mmCIF raw/3OG7.cif, SHA-256 "
            + next(s for s in man["structures"] if s["pdb_id"] == "3OG7")["file"]["sha256"]) in md
    assert rec["ligand"]["sdf_sha256"] in md and "ligands/3OG7_032_C.sdf" in md
    assert "Receptor clean/3OG7_receptor.pdb, chains A, B" in md
    assert "Pocket 032 C: box 26.3 × 15.13 × 16.4 Å, centre (2.643, -2.28, -19.403)" in md
    assert "Construct mutations recorded by get_structures: field: I544A, F667E" in md

    # the ligand is named only by the wwPDB synonyms stored in the manifest
    assert "wwPDB synonyms: Vemurafenib, PLX4032" in md and rec["ligand"]["smiles"] in md
    assert "rotatable torsions: 7" in md

    # calculated values, labelled, with the literal scope and the negative control of the run
    assert "| Vina score (kcal/mol, calculated) |" in md and "RMSD to the crystal pose (Å, calculated)" in md
    assert "| 1 | -11.35 | 0.5396 |" in md and "| 9 | -10.35 | 8.2993 |" in md
    assert "Top-pose RMSD: 0.5396 Å; criterion <= 2 Å; within the criterion: True" in md
    assert "Pre-registered control run: True" in md
    assert rec["redocking"]["scope"] in md and "not the accuracy of the method" in md
    assert ("the crystal pose against itself gave 0 Å and the crystal pose moved 5 Å gave 5 Å; "
            "passed: True") in md
    assert "9 poses returned (num_modes 9)." in md and "the rest were outside" not in md
    # the vocabulary line does not call the experimental resolution a calculated value
    assert "Every number in section 1 is CALCULATED" not in md
    assert "the resolution is the experimental value deposited for the PDB entry" in md

    # Meeko and the guard, as recorded (3OG7 passes at once: nothing deleted)
    assert "Meeko mk_prepare_receptor --read_pdb --delete_bad_res_from_box_radius 8; exit code 0" in md
    assert "Deleted by Meeko (outside the box, so the pocket is unchanged): none" in md
    assert "Refused by Meeko (too close to the box, not deleted): none" in md
    assert "Ligand atoms inside the box: " + rec["guard"]["ligand_atoms_in_box"] in md

    # files with their SHA-256, straight from the manifest
    for name, f in rec["files"].items():
        assert f"| {name} | {f['path']} | {f['sha256']} |" in md, name

    # the Vina log in the dossier really has the header; the report does not copy it
    log = (prep / rec["files"][f"vina_032_regenerated.log"]["path"]).read_text(encoding="utf-8")
    assert VINA_LOG_HEADER in log, "the fixture log lost the header: this check would prove nothing"
    assert VINA_LOG_HEADER not in md and "Estimated Free Energy of Binding" not in md
    assert "affinity" not in md.lower().replace("measured affinity", ""), "the report calls a score an affinity"

    # ChEMBL potency is not in this document: every IC50/Ki/Kd mention is the caveat
    assert "CHEMBL" not in md and "pChEMBL" not in md
    mentions = [l for l in md.splitlines() if "IC50" in l or "Ki/Kd" in l]
    assert len(mentions) == 3 and all("not comparable" in l for l in mentions), mentions

    links_are_sourced(md, manifest_text)
    # control: the checker does catch a link that the manifest does not back
    try:
        links_are_sourced(md + "[x](https://example.org/3OG7)", manifest_text)
        raise AssertionError("the link check passed a fabricated link")
    except AssertionError as exc:
        assert "not backed by the manifest" in str(exc)

    pdf = Path(r["pdf_path"]).read_bytes()
    assert pdf.startswith(b"%PDF") and len(pdf) > 10_000 and r["pages"] >= 2
    return rec


def test_user_choice_is_quoted_only_where_it_applies():
    """The reason belongs to one structure and pocket: it must not travel to another run."""
    md = md_of(build())
    assert "No structure recorded by select_structure in this manifest." in md
    other = choice.select_structure(ACC, "5JSM", "6NB J", chains=["A"], reason="es el vemurafenib dimerizado")
    assert other["ok"], other
    md = md_of(build())
    assert "Not the structure and pocket recorded by select_structure (5JSM chain A, pocket 6NB J" in md
    assert "vemurafenib dimerizado" not in md, "the reason given for 5JSM was quoted on the 3OG7 run"
    assert choice.select_structure(ACC, "3OG7", "032 C", reason="pose of the same compound")["ok"]
    md = md_of(build())
    assert 'reason recorded by select_structure: "pose of the same compound"' in md
    assert "Chosen by the user: 3OG7 chain A, pocket 032 C" in md


def test_more_runs():
    """A biased run and a run that is not the pre-registered control are listed as such."""
    assert run("3OG7", "032 C", conformer="crystal")["ok"]
    md = md_of(build())
    assert "run 032_crystal" in md and "run 032_regenerated" in md
    assert "Started from the crystal conformation: biased in favour" in md
    r2 = build()
    assert len(r2["docking_runs"]) == 2 and "2 docking runs recorded" in md_of(r2)
    assert run("3OG7", "032 C", exhaustiveness=8, seed=7)["ok"]
    md = md_of(build())
    assert "Pre-registered control run: False" in md and "not the pre-registered control" in md
    assert "Parameters: exhaustiveness 8, seed 7" in md
    assert "Replaced an earlier run of the same pocket and conformer, recorded" in md


def test_meeko_messages_and_unmeasured_cpu():
    """Real Meeko output (8C7Y TXV U deleted, 5JSM 6NB J refused) and cpu_seconds that could not be read."""
    man = json.loads(choice.manifest_path(ACC).read_text(encoding="utf-8"))
    block = next(s for s in man["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]
    base = copy.deepcopy(block["runs"]["032_regenerated"])
    deleted = dock.parse_meeko((DFIX / "meeko_8c7y_TXV_U_ok.stdout.txt").read_text(encoding="utf-8"),
                               (DFIX / "meeko_8c7y_TXV_U_ok.stderr.txt").read_text(encoding="utf-8"))
    refused = dock.parse_meeko((DFIX / "meeko_5jsm_6NB_J_refused.stdout.txt").read_text(encoding="utf-8"),
                               (DFIX / "meeko_5jsm_6NB_J_refused.stderr.txt").read_text(encoding="utf-8"))
    edited = copy.deepcopy(base)
    edited["receptor"].update(deleted, warning=dock.unaccounted_warning({**deleted, "deleted": deleted["deleted"][:2]}),
                              chains_note="user_choice names chain A; the receptor has chains A, B")
    edited["receptor"]["refused"] = refused["refused"]
    edited["vina"] = {k: v for k, v in base["vina"].items() if k not in ("cpu_seconds", "cpu_per_wall_second")}
    edited["vina"]["cpu_seconds_note"] = "Not measurable: the executable in PROTEIN_MCP_VINA used no CPU itself."
    block["runs"] = {"032_regenerated": edited}
    choice.write_manifest(choice.manifest_path(ACC), man)
    md = md_of(build())
    assert "A:607 (17.487 Å outside the box), A:608 (14.981 Å outside the box), A:609 (11.859 Å outside the box)" in md
    assert "A:600 (2.456 Å from the box edge)" in md and "B:488 (3.912 Å from the box edge)" in md
    assert "Residues matching no Meeko template: A:607, A:608, A:609" in md
    assert "Warning: Meeko reported residues without a template" in md and "A:609" in md
    assert "user_choice names chain A; the receptor has chains A, B" in md
    # an unmeasurable CPU time is said to be unmeasured: no ratio is printed
    assert "CPU time: not measured." in md and "s of CPU per second of wall clock" not in md
    return man


def test_missing_fields_and_earlier_choice():
    """Fields absent from the manifest read 'no record' (never None); fewer poses than num_modes give no
    invented cause; a run started from a choice that was later replaced says so."""
    man = json.loads(choice.manifest_path(ACC).read_text(encoding="utf-8"))
    s = next(s for s in man["structures"] if s["pdb_id"] == "3OG7")
    real = s["docking"]["032 C"]["runs"]["032_regenerated"]
    bare = copy.deepcopy(man)
    bare.pop("created_at", None)
    bare.pop("generated_by", None)
    b3 = next(x for x in bare["structures"] if x["pdb_id"] == "3OG7")
    b3["docking"]["032 C"] = {"runs": {"032_regenerated": {"poses": [], "vina": {}, "redocking": {},
                                                          "guard": {"ligand_atoms_in_box": "1 of 1"}}}}
    def render(m):
        return report.render_markdown(experiment_report.build_document(ACC, m, "2026-11-06T00:00:00Z")[0])
    md = render(bare)
    for bad in ("None", " Å, centre ()", "conformer ''"):
        assert bad not in md, (bad, [l for l in md.splitlines() if bad in l][:3])
    assert "created no record by no record" in md and "heavy atoms: no record" in md
    assert "within the criterion: no record" in md and "box no record" in md
    # control: with the 'no record' helper switched off, the same manifest does print None
    real_v = experiment_report._v
    experiment_report._v = lambda x: x
    try:
        assert "None" in render(bare), "the None check cannot fail: it proves nothing"
    finally:
        experiment_report._v = real_v

    fewer = copy.deepcopy(man)
    f3 = next(x for x in fewer["structures"] if x["pdb_id"] == "3OG7")
    run = f3["docking"]["032 C"]["runs"]["032_regenerated"]
    run["poses"] = run["poses"][:7]
    md = render(fewer)
    assert "7 poses returned; num_modes 9" in md and "the manifest does not record why fewer were returned" in md
    assert "outside that window" not in md

    # the run was started from user_choice = 3OG7 032 C, and the current choice is then the earlier 5JSM one
    later = copy.deepcopy(man)
    l3 = next(x for x in later["structures"] if x["pdb_id"] == "3OG7")
    l3["docking"]["032 C"]["runs"]["032_regenerated"]["selected_by"] = "user_choice"
    later["user_choice"] = next(h for h in later["user_choice_history"] if h["structure"] == "5JSM")
    md = render(later)
    assert "started from an earlier choice, since replaced" in md and "user_choice_history" in md
    assert "vemurafenib dimerizado" not in md and man["user_choice"]["reason"] not in md
    # control: with the current choice still 3OG7 032 C, its reason is quoted and the extra line is absent
    md = render(man)
    assert "started from an earlier choice" not in md and man["user_choice"]["reason"] in md
    assert real["computed_at"]


def test_choice_recorded_again_after_the_run():
    """Second review: the same structure and pocket chosen again with another reason after a run that
    used the first choice must not get the new reason; a dropped choice is not 'no choice ever'."""
    import time
    assert choice.select_structure(ACC, "3OG7", "032 C", reason="REASON-ONE")["ok"]
    assert run()["ok"]                                       # no arguments: started from user_choice
    man = json.loads(choice.manifest_path(ACC).read_text(encoding="utf-8"))
    rec = next(s for s in man["structures"] if s["pdb_id"] == "3OG7")["docking"]["032 C"]["runs"]["032_regenerated"]
    assert rec["selected_by"] == "user_choice" and rec["user_choice_chosen_at"] == man["user_choice"]["chosen_at"], rec
    md = md_of(build())
    # control: while the choice is the one the run used, its reason is quoted
    assert 'reason recorded by select_structure: "REASON-ONE"' in md
    time.sleep(1.1)                                          # chosen_at has a resolution of one second
    assert choice.select_structure(ACC, "3OG7", "032 C", reason="REASON-TWO")["ok"]
    md = md_of(build())
    assert "REASON-TWO" not in md and "REASON-ONE" not in md, "a reason was quoted on a run that did not use it"
    assert "Its reason is not quoted: the choice was recorded again" in md and "user_choice_history" in md

    dropped = json.loads(choice.manifest_path(ACC).read_text(encoding="utf-8"))
    dropped["user_choice"] = None                            # what choice.carry_over leaves when it drops one
    md = report.render_markdown(experiment_report.build_document(ACC, dropped, "2026-11-06T00:00:00Z")[0])
    assert "No structure recorded by select_structure in this manifest." not in md
    assert "No current choice in this manifest (select_structure); this run was started from an earlier one" in md


def test_research_report_gains_one_line_without_figures():
    """The research report says the experiment record exists and when, and repeats no figures."""
    folder = dossier.DOSSIER_DIR / ACC
    (folder / "structures_v600e.json").write_text(json.dumps({
        "tool": "get_structures", "params": {"mutation": "V600E"}, "retrieved_at": "2026-10-21T21:18:49Z",
        "result": {"ok": True, "experimental": {"ok": True, "structures": []}, "alphafold": {"ok": False}}}),
        encoding="utf-8")
    r = report.build_report(ACC)
    assert r["ok"], r
    sec6 = Path(r["markdown_path"]).read_text(encoding="utf-8").split(
        "## 6. Structure preparation (Phase 2)")[1].split("## 7. Methods")[0]
    assert "Experiment record: 1 docking run recorded in this manifest, the most recent on" in sec6, sec6[-1200:]
    assert f'build_experiment_report("{ACC}")' in sec6 and f"experiments_{ACC}.pdf" in sec6
    for figure in ("0.5396", "-11.35", "kcal/mol", "exhaustiveness", "Vina"):
        assert figure not in sec6, f"section 6 of the research report repeats a docking figure: {figure}"


def test_the_section_inventory_is_read_from_the_document():
    """F4 (Desktop test of 2026-11-12): asked what section 3 said, the answer was that the tool
    output does not number the sections -- true, and it left the analyst unable to know what the
    document holds without opening it. The inventory is read BACK from the blocks that were
    written, so it cannot drift from the document the way a hand-kept list would.

    CONTROL: the status of each section is asserted against what the document itself prints, and
    the count is tied to the number of level-2 headings, so an inventory that silently stopped
    following the document would fail here.
    """
    r = build()
    assert r["ok"], r
    sections = r["sections"]
    doc, _n, _n_admet = experiment_report.build_document(
        ACC, report.load_manifest(ACC), "2026-11-12T00:00:00Z")
    headings = [b[2] for b in doc.blocks if b[0] == "h" and b[1] == 2]
    assert len(sections) == len(headings) == 4, (sections, headings)
    assert [s["number"] for s in sections] == ["1", "2", "3", "4"], sections
    assert [s["heading"] for s in sections] == headings
    by_number = {s["number"]: s for s in sections}
    assert "docking" in by_number["1"]["title"].lower()
    assert by_number["1"]["status"] == "has content", by_number["1"]
    # section 2 is the one declared hole, and it says so in the document
    assert by_number["2"]["title"] == "Molecular dynamics"
    assert by_number["2"]["status"] == "DECLARED HOLE (not run)", by_number["2"]
    assert by_number["2"]["says"].startswith("Not run."), by_number["2"]
    assert by_number["4"]["status"] == "has content"
    # a section with no content must not be reported as having any. This manifest records no
    # ADMET result, so section 3 IS a hole -- asserted flatly, because the first version of this
    # test asked `if "Not run." in md: ... else: assert "has content"`, and that `else` ACCEPTED
    # the lie: dropping the "Not run. " prefix from the placeholder made the inventory call an
    # empty section "has content" and the suite stayed green (reviewer finding 2026-11-12).
    assert by_number["3"]["title"] == "ADMET"
    assert by_number["3"]["status"] == "DECLARED HOLE (not run)", by_number["3"]
    assert "analyze_admet" in by_number["3"]["says"]
    for s in sections:
        assert ("says" in s) == (s["status"] != "has content"), s
    # cross-check between the two fields of the same result: a section the inventory calls full
    # while sections_not_run calls it empty is a contradiction the analyst would be handed, and
    # this is what catches a placeholder whose wording drifts away from the inventory's rule
    not_run = " ".join(r["sections_not_run"]).lower()
    keyword = {"1": "docking", "2": "dynamics", "3": "admet", "4": "methods"}
    for s in sections:
        named = keyword[s["number"]] in not_run
        assert named == (s["status"] != "has content"), (s, r["sections_not_run"])


def test_open_pdf_error():
    md_path = dossier.DOSSIER_DIR / ACC / f"experiments_{ACC}.md"
    before = md_path.read_bytes()
    md_path.write_bytes(before + b"\nmarker of the previous build\n")
    before = md_path.read_bytes()
    real = report.render_pdf
    report.render_pdf = lambda *a, **k: (_ for _ in ()).throw(PermissionError(13, "Permission denied"))
    try:
        r = build()
    finally:
        report.render_pdf = real
    assert not r["ok"] and "close it and try again" in r["error"], r
    assert md_path.read_bytes() == before, "the Markdown was replaced although the PDF could not be written"


def main():
    record = make_fake_vina(Path(tempfile.mkdtemp()))
    prep = test_no_runs()
    rec = test_real_run(prep)
    test_user_choice_is_quoted_only_where_it_applies()
    test_more_runs()
    test_meeko_messages_and_unmeasured_cpu()
    test_missing_fields_and_earlier_choice()
    test_choice_recorded_again_after_the_run()
    test_research_report_gains_one_line_without_figures()
    test_the_section_inventory_is_read_from_the_document()
    test_open_pdf_error()
    assert calls(record), "the Vina stand-in was never called: the run was not real"
    print(json.dumps({"pages": build()["pages"], "runs": len(build()["docking_runs"]),
                      "top_pose_rmsd": rec["redocking"]["top_pose_rmsd_angstrom"]}))
    print("ALL OK")


main()
