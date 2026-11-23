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
    # with no MD block in the manifest, section 2 says so and names the tool that records a run
    assert "## 2. Molecular dynamics" in md
    sec2 = md.split("## 2. Molecular dynamics")[1].split("## 3.")[0]
    assert "Not run. No molecular dynamics package, job or analysis is recorded" in sec2, sec2
    assert f'analyze_md(folder, uniprot_id="{ACC}"' in sec2, sec2
    assert "molecular dynamics (no analyze_md result in the manifest)" in r["sections_not_run"]
    # Section 3 exists from the start and declares the gap by naming the tool that fills it.
    # It must NOT call anything PREDICTED: nothing analyze_admet writes comes from a trained
    # model (ESTADO 7.3), and the placeholder used to promise "PREDICTED values".
    assert "## 3. ADMET" in md and "analyze_admet" in md
    assert "ADMET (no analyze_admet result in the manifest)" in r["sections_not_run"]
    assert "PREDICTED" not in md and "predicted" not in md
    # H5: the pending line does not contradict the two sections declared not run
    assert "Preparation steps still pending in the manifest:" in md
    assert "Molecular dynamics (section 2): no analysed run is recorded." in md
    assert "nothing listed" not in md
    return prep


def test_real_run(prep: Path):
    """One real docking run (3OG7 032 C, regenerated, the pre-registered control)."""
    # the result travels into the assert: a refusal here used to abort the suite with a bare
    # AssertionError and no way to tell a real regression from a flaky launch
    r0 = run("3OG7", "032 C")
    assert r0["ok"], r0
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
    # the number, its criterion and the verdict now travel in the opening paragraph of the run
    # instead of a bullet, but all three must still be there and must still agree
    assert ("The best-scoring pose lands 0.5396 Å from the crystallographic pose against a "
            "criterion of 2 Å, so this set-up reproduces a pose it already knew.") in md
    assert "does NOT reproduce" not in md
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
    assert ("Checked before Vina started: ligand atoms inside the box "
            + rec["guard"]["ligand_atoms_in_box"]) in md

    # files with their SHA-256, straight from the manifest. They moved out of the reading path
    # into the reproducibility block, but not one of them may be dropped: the hash is what proves
    # a file on disk is the file this run wrote.
    audit = md.index("## 6. Reproducibility and audit trail")
    assert md.index("## 4. The compound set") < md.index("## 5. Methods") < audit, \
        "the set table goes after the three result sections, and the audit trail last"
    assert "### Docking · 3OG7 · pocket 032 C · run 032_regenerated" in md[audit:]
    for name, f in rec["files"].items():
        assert f"{name}: {f['path']} (SHA-256 {f['sha256']})" in md[audit:], name
        # and NOT in the reading path, which is the whole point of moving it
        assert f["sha256"] not in md[:audit], f"{name}: a hash leaked into sections 1-5"

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
    r0 = run("3OG7", "032 C", conformer="crystal")
    assert r0["ok"], r0
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
    # an unrecorded verdict must print as unrecorded, and must NOT be rendered as a failure
    assert "Whether it is within the criterion: no record" in md
    assert "does NOT reproduce" not in md, "a missing verdict was turned into a failed one"
    assert "box no record" in md
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
    doc, _n, _n_admet, _n_set = experiment_report.build_document(
        ACC, report.load_manifest(ACC), "2026-11-12T00:00:00Z")
    headings = [b[2] for b in doc.blocks if b[0] == "h" and b[1] == 2]
    assert len(sections) == len(headings) == 6, (sections, headings)
    assert [s["number"] for s in sections] == ["1", "2", "3", "4", "5", "6"], sections
    assert "reproducibility" in sections[-1]["title"].lower(), "the audit trail must be the last section"
    assert [s["heading"] for s in sections] == headings
    by_number = {s["number"]: s for s in sections}
    assert "docking" in by_number["1"]["title"].lower()
    assert by_number["1"]["status"] == "has content", by_number["1"]
    # section 2 is the one declared hole, and it says so in the document
    assert by_number["2"]["title"] == "Molecular dynamics"
    assert by_number["2"]["status"] == "DECLARED HOLE (not run)", by_number["2"]
    assert by_number["2"]["says"].startswith("Not run."), by_number["2"]
    assert by_number["5"]["status"] == "has content"
    # section 4 (the compound set) is a hole here too: this manifest records no
    # select_ligands selection, and a table of the set with nothing in it would be a claim
    assert by_number["4"]["title"] == "The compound set side by side"
    assert by_number["4"]["status"] == "DECLARED HOLE (not run)", by_number["4"]
    assert "select_ligands" in by_number["4"]["says"], by_number["4"]
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
    keyword = {"1": "docking", "2": "dynamics", "3": "admet", "4": "compound set",
               "5": "methods", "6": "reproducibility"}
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


def test_the_cluster_install_path_stays_out_of_the_report():
    """2026-11-15, the user's decision: the shareable report names the NAMD release, not where the
    cluster installs it. The full path stays in the manifest."""
    from protein_mcp import experiment_report as er
    from protein_mcp.report import Doc

    real = "/software/local/apps/NAMD_3.0.3_Linux-x86_64-multicore-CUDA/namd3"
    job = {"site": {"namd_binary": {"value": real, "mark": "CONFIRMED"}}}

    def cluster_row(tail):
        saved = er._binary_tail
        er._binary_tail = tail
        try:
            doc = Doc()
            er._md_config(doc, {}, job, [], "test pocket")
        finally:
            er._binary_tail = saved
        rows = [r for b in doc.blocks if b[0] == "table" for r in b[2]]
        return next(r[1] for r in rows if r[0] == "Where it ran")

    def clean(text):
        return "/software" not in text and "local/apps" not in text and \
            "NAMD_3.0.3_Linux-x86_64-multicore-CUDA/namd3" in text and "full path in the manifest" in text

    row = cluster_row(er._binary_tail)
    assert clean(row), row
    # control: without the function the install path is printed, and the same check must fail
    assert not clean(cluster_row(lambda v: v))
    assert er._binary_tail(None) is None and er._binary_tail("CONFIRM") == "CONFIRM"
    assert er._binary_tail("namd3") == "namd3"
    print("ok  the report names the NAMD release, not the cluster's install path")


SET_REASON = ("five candidates by the agreed criteria; VEMURAFENIB is the control/anchor, not a "
              "candidate: it is the crystal ligand of 3OG7")
# the column order of section 4, so the assertions below name what they read
COMPOUND, SCORE, LE, XRMSD, MDRMSD, CONTROL, ADMET = range(7)


def _set_manifest():
    """A manifest with a selection of 3: one with docking+MD+ADMET and a cross-docking RMSD whose
    error bar decides nothing, one docking-only and box_limited, one with no name in ChEMBL, no
    ChEMBL entry at all, and TWO runs of which the first inserted is the biased crystal one."""
    return {
        "created_at": "2026-11-22T00:00:00Z",
        "ligand_choice": {"reason": SET_REASON, "compounds": [
            {"id": "CHEMBL1229517", "name": "VEMURAFENIB", "level": "T1", "heavy_atoms": 33},
            {"id": "CHEMBL1336", "name": "SORAFENIB", "level": "T3", "heavy_atoms": 32},
            {"id": "CHEMBL525191", "name": None, "level": "T3", "heavy_atoms": 25}]},
        "structures": [{"pdb_id": "8C7Y", "docking": {"TXV U": {"runs": {
            "r1": {"ligand": {"conformer": "regenerated"}, "compound_docking": {
                "compound_id": "CHEMBL1229517", "name": "VEMURAFENIB", "level": "T1",
                "top_pose_score_kcal_mol": -11.191, "heavy_atoms": 33,
                "ligand_efficiency_kcal_mol_per_heavy_atom": -0.3391, "box_limited": False,
                "scope": "...no inactive compounds in this run, so nothing here says whether the "
                         "score separates actives from inactives (ESTADO 7.4)."},
                "cross_docking": {
                    "top_pose_rmsd_angstrom": 1.9152, "criterion_angstrom": 2.0,
                    "within_criterion": True,
                    "reference_pose_uncertainty": {"shift_angstrom": 0.767},
                    "within_criterion_with_bar": {"range_angstrom": [1.148, 2.682],
                                                  "decides": False}}},
            "r2": {"ligand": {"conformer": "regenerated"}, "compound_docking": {
                "compound_id": "CHEMBL1336", "name": "SORAFENIB", "level": "T3",
                "top_pose_score_kcal_mol": -11.337, "heavy_atoms": 32,
                "ligand_efficiency_kcal_mol_per_heavy_atom": -0.3543, "box_limited": True,
                "box_fit": {"min_clearance_angstrom": 1.99, "limiting_box_axis_angstrom": 22.0,
                            "limiting_ligand_extent_angstrom": 20.01, "threshold_angstrom": 1.5}}},
            # the biased run is inserted FIRST, so a table that prints "the first match" prints
            # this one: it must then say which it printed, of how many, and that it is biased
            "r3_crystal": {"ligand": {"conformer": "crystal"}, "compound_docking": {
                "compound_id": "CHEMBL525191", "name": None, "level": "T3",
                "top_pose_score_kcal_mol": -13.5, "heavy_atoms": 25,
                "ligand_efficiency_kcal_mol_per_heavy_atom": -0.54, "box_limited": False}},
            "r3": {"ligand": {"conformer": "regenerated"}, "compound_docking": {
                "compound_id": "CHEMBL525191", "name": None, "level": "T3",
                "top_pose_score_kcal_mol": -8.72, "heavy_atoms": 25,
                "ligand_efficiency_kcal_mol_per_heavy_atom": -0.3488, "box_limited": False}}}}},
            "md_analysis": {"TXV U CHEMBL1229517": {"production": {
                "stage_label": "production", "pocket": "TXV U",
                "compound": {"id": "CHEMBL1229517", "name": "VEMURAFENIB"},
                "result": {"ligand_rmsd": {"mean": 1.5799, "n": 100},
                           "reading": {"equilibration_control": {
                               "status": "NOT MET", "drift_angstrom": 2.9923,
                               "threshold_angstrom": 0.5}}}}}},
            "admet": {"TXV U CHEMBL1229517": {
                "compound_identity": {"compound_id": "CHEMBL1229517", "pdb_id": "8C7Y",
                                      "pocket": "TXV U"},
                "chembl_record": {"status": "RECORD (facts of the ChEMBL entry)"}},
                "TXV U CHEMBL525191": {
                    "compound_identity": {"compound_id": "CHEMBL525191", "pdb_id": "8C7Y",
                                          "pocket": "TXV U"},
                    "chembl_record": {"status": "NO RECORD: no entry for this InChIKey"}}}}],
    }


def test_the_compound_set_table_declares_its_holes_and_ranks_nothing():
    """Section 4 (researcher's request 2026-11-22): the set side by side, WITHOUT a ranking.

    CONTROLS: the rows must stay in the order select_ligands recorded them (a table sorted by
    score would be a ranking even with no verdict written); a compound with no MD must say so
    instead of leaving the cell empty; a box_limited run must carry its warning AND its clearance;
    a compound ChEMBL does not name must print "no record" rather than borrow a name; a figure
    printed out of several runs must say so; and an error bar that straddles the criterion must
    say it decides nothing instead of standing beside the criterion as if it had passed.
    """
    m = _set_manifest()
    doc = experiment_report.Doc()
    n = experiment_report._compound_set(doc, m, ACC)
    assert n == 3, n
    table = next(b for b in doc.blocks if b[0] == "table")
    rows = table[2]
    assert len(table[1]) == len(rows[0]) == 7, (table[1], rows[0])
    # CONTROL: the order is the selection's, not the scores'. Sorted by score the first row
    # would be CHEMBL525191 (-13.5), and by ligand efficiency it would be CHEMBL525191 too.
    assert [r[COMPOUND].split("\n")[1] for r in rows] == [
        "CHEMBL1229517", "CHEMBL1336", "CHEMBL525191"], rows
    md = report.render_markdown(doc)
    assert "THIS TABLE IS NOT A RANKING" in md
    assert "Docking recorded for 3 of 3, molecular dynamics for 1 of 3, ADMET for 2 of 3" in md, md
    # the reason is printed WHOLE: the clipped version cut the sentence that says the anchor is
    # not a candidate, leaving "five candidates" beside a count of three rows
    assert SET_REASON in md, "the recorded reason was clipped"
    assert "is the control/anchor, not a candidate" in md
    # the limits the pre-registration (ESTADO 7.4) obliges this table to print
    assert "NO NEGATIVE CONTROL" in md and "separates actives from inactives" in md
    assert "levels T1, T3" in md and "not ranked together" in md
    # the range of the SET (33, 32, 25 in this fixture), which is what makes a raw score
    # incomparable across it; the figure is the selection's, not one compound's
    assert "heavy-atom count runs from 25 to 33" in md, md
    assert "ONE replica per compound" in md
    # the three provenances are labelled, in the headers and in the closing note
    assert all("calculated" in h for h in table[1][SCORE:CONTROL + 1]), table[1]
    assert "retrieved" in table[1][ADMET].lower()
    assert "THE THREE PROVENANCES IN THIS TABLE DO NOT MEAN THE SAME THING" in md
    # the compound with no dynamics says which tools would fill it, in its own cell
    assert "not run (build_md_package, then analyze_md)" in rows[1][MDRMSD]
    # the one that was run carries the figure, its provenance AND the control that qualifies it
    assert "1.5799 mean over 100 frames" in rows[0][MDRMSD] and "8C7Y TXV U" in rows[0][MDRMSD]
    assert "NOT MET" in rows[0][CONTROL] and "2.9923" in rows[0][CONTROL]
    # the cross-docking RMSD of the T1 anchor does not pose as a pass
    assert "1.9152" in rows[0][XRMSD] and "0.767" in rows[0][XRMSD]
    assert "DECIDES NOTHING" in rows[0][XRMSD] and "1.148-2.682" in rows[0][XRMSD], rows[0][XRMSD]
    # a T3 compound has no reference pose: not applicable, not missing
    assert "not applicable" in rows[1][XRMSD] and "level T3" in rows[1][XRMSD]
    # box_limited is not silently comparable, and the clearance travels with the mark
    assert "box_limited" in rows[1][SCORE] and "1.99 A clearance" in rows[1][SCORE], rows[1][SCORE]
    assert "box_limited" not in rows[0][SCORE], rows[0][SCORE]
    # one of several runs: which one, of how many, and that this one is the biased crystal run
    assert "1 of 2 runs recorded" in rows[2][SCORE], rows[2][SCORE]
    assert "biased in favour of the crystal pose" in rows[2][SCORE], rows[2][SCORE]
    # a compound ChEMBL does not name is not given one
    assert rows[2][COMPOUND].startswith(experiment_report.NO_RECORD), rows[2][COMPOUND]
    # no ChEMBL entry is an absence of DATA: in an ADMET column, silence reads as "no liabilities"
    assert "absence of DATA, not an absence of liabilities" in rows[2][ADMET], rows[2][ADMET]
    assert "ChEMBL entry retrieved" in rows[0][ADMET]
    # and with no selection recorded the section is a declared hole, not an empty table
    empty = experiment_report.Doc()
    assert experiment_report._compound_set(empty, {"structures": []}, ACC) == 0
    text = report.render_markdown(empty)
    assert "Not run." in text and "select_ligands" in text
    assert not [b for b in empty.blocks if b[0] == "table"], "an empty set still drew a table"
    print("ok  the compound set is laid side by side with its holes declared and no ranking")


def test_a_dynamics_record_written_before_step_6b_is_not_called_not_run():
    """A record from before the compound block existed is matched by its key instead.

    "not run" is a claim about what was run. Keyed only on the compound block, a trajectory that
    HAD been analysed was printed as "not run" in section 4 while section 2 of the same document
    described it, and the coverage count said 0 of 1 (reviewer H3, 2026-11-22).
    """
    m = _set_manifest()
    legacy = m["structures"][0]["md_analysis"]["TXV U CHEMBL1229517"]["production"]
    legacy.pop("compound")
    legacy.pop("pocket")
    doc = experiment_report.Doc()
    experiment_report._compound_set(doc, m, ACC)
    md = report.render_markdown(doc)
    rows = next(b for b in doc.blocks if b[0] == "table")[2]
    assert "molecular dynamics for 1 of 3" in md, md
    assert "1.5799 mean over 100 frames" in rows[0][MDRMSD], rows[0][MDRMSD]
    assert "not run" not in rows[0][MDRMSD], rows[0][MDRMSD]
    # and the pocket still travels with the figure, read off the key
    assert "TXV U" in rows[0][MDRMSD], rows[0][MDRMSD]
    # CONTROL: a record of ANOTHER compound's key is not borrowed for this row
    m2 = _set_manifest()
    m2["structures"][0]["md_analysis"] = {"TXV U CHEMBL9999999": {"production": {
        "stage_label": "production",
        "result": {"ligand_rmsd": {"mean": 0.4, "n": 10}, "reading": {}}}}}
    doc2 = experiment_report.Doc()
    experiment_report._compound_set(doc2, m2, ACC)
    md2 = report.render_markdown(doc2)
    assert "molecular dynamics for 0 of 3" in md2, md2
    assert "0.4 mean" not in md2, "a record keyed to another compound was printed in this row"
    print("ok  an MD record written before step 6b counts, and another compound's does not")


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
    test_the_cluster_install_path_stays_out_of_the_report()
    test_the_compound_set_table_declares_its_holes_and_ranks_nothing()
    test_a_dynamics_record_written_before_step_6b_is_not_called_not_run()
    assert calls(record), "the Vina stand-in was never called: the run was not real"
    print(json.dumps({"pages": build()["pages"], "runs": len(build()["docking_runs"]),
                      "top_pose_rmsd": rec["redocking"]["top_pose_rmsd_angstrom"]}))
    print("ALL OK")


main()
