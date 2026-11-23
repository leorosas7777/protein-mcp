"""build_experiment_report: prep/manifest.json -> Markdown + PDF record of the in-silico experiments.

A document of its own, separate from build_report (decided 2026-11-05), with the three
result sections present from the start: docking, molecular dynamics, ADMET, joined on
2026-11-22 by section 4, which lays the selected compounds side by side without ranking
them (methods moved to 5 and the audit trail to 6). Sections with no run
say "not run" and name the tool that would fill them, exactly as build_report does with
missing dossier sections.

Two rules shape this module:
  * self-contained provenance: every run states the PDB entry, resolution, pocket, box,
    the SHA-256 of the mmCIF it came from, the receptor chains and, when the run is the
    structure the user chose, the literal reason recorded by select_structure. The reader
    never has to open the research report;
  * the vocabularies are not mixed (pre-registration, ESTADO 7.1): the docking results here
    are CALCULATED; the provenance quotes what the manifest recorded (experimental
    resolution, box from the preparation step, times measured on this machine). Measured potency (ChEMBL IC50/Ki/Kd) stays in the research report, and
    predicted values (ADMET) get their own section when that step exists.

Values are quoted from prep/manifest.json; nothing is recomputed, re-ranked or inferred,
and no text of a Vina log is copied (its "affinity (kcal/mol)" header would label a
calculated score as a measured affinity). The rendering machinery lives in report.py.
"""
import hashlib
import json
from pathlib import Path

from protein_mcp import dossier, report
from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice
from protein_mcp.report import NO_RECORD, Doc, Lines, Link, _box, _num

# Section 2 is written by _md() from the md_package, md_job and md_analysis blocks of the
# manifest. Rules decided with the user on 2026-11-12: the equilibration control heads every run
# so no ligand figure can be read without it; the four metrics are always printed; the full
# backbone is the headline and the observed/built split is marked as not pre-registered; the
# hydrogen-bond figure is a share of FRAMES; every figure says which run (production or
# equilibration) it comes from; plots are drawn from the recorded series file, never recomputed.
MD_STAGES = ("production", "equilibration")
MAX_PLOT_POINTS = 500
MD_VOCABULARY = (
    "Section 2 holds CALCULATED values: a molecular dynamics simulation run with NAMD and measured "
    "with MDAnalysis by analyze_md, quoted from the manifest. Not one of them is a binding affinity "
    "or a free energy, and none is comparable with the docking scores of section 1 or the "
    "measured potencies of the research report. The thresholds in the control line and in the "
    "offered reading are conventions, not measurements of this system: they were written down "
    "before the first production trajectory existed, but after the restrained equilibration of "
    "the same system had been measured, so they cannot have been tuned to a production result "
    "and could have been to that equilibration.")
# Section 3 is written by _admet() from the `admet` blocks analyze_admet records in the manifest.
# What it must never do, decided when that tool was measured against the real API (ESTADO 7.3):
# keep the three provenances in three tables, never average or convert a measured value, keep a
# censored relation attached to its number, and present the rules as offered conventions.
ADMET_VOCABULARY = (
    "Section 3 holds three different kinds of statement, in three separate tables, and they are "
    "not comparable with each other or with the calculated docking scores of section 1: the "
    "ChEMBL record (facts of a database entry), the retrieved measurements (each with the units "
    "and the assay it was measured in) and the calculated properties (computed from the ligand's "
    "atoms and bonds). NOTHING in this section is a prediction: no value here comes from a "
    "trained model.")

VOCABULARY = ("The docking results in section 1 (Vina scores and RMSD values) are CALCULATED by AutoDock Vina "
              "and RDKit, as recorded in the manifest. The provenance lines quote other kinds of value: the "
              "resolution is the experimental value deposited for the PDB entry, the box was computed when the "
              "structure was prepared, and run times are measurements of this machine. A docking score is not a measured affinity and is not comparable with the IC50/Ki/Kd values that "
              "ChEMBL reports for the same target, which stay in the research report. The poses of one run are "
              "listed in the order Vina returned them; no compound is ranked against another here.")


def _v(x):
    """A recorded value as text, or 'no record' when the manifest has none (never 'None')."""
    return NO_RECORD if x is None or x == "" else x


def _boxtext(box: dict) -> str:
    return _box(box) if (box or {}).get("size") and (box or {}).get("center") else NO_RECORD


def _run_rows(m: dict) -> list:
    """[(structure, pocket_key, block, run_id, run)] for every docking run in the manifest."""
    out = []
    for s in m.get("structures") or []:
        for pocket_key, block in sorted((s.get("docking") or {}).items()):
            for run_id, run in sorted((block.get("runs") or {}).items()):
                out.append((s, pocket_key, block, run_id, run))
    return out


def _choice_line(m: dict, pdb: str, pocket_key: str, run: dict):
    """The user's choice, only where it really applies: the literal reason when this run is the
    chosen structure and pocket, otherwise what was chosen instead (never the other way round)."""
    uc = m.get("user_choice") or {}
    from_choice = run.get("selected_by") == "user_choice"
    earlier = ("earlier choices are kept in user_choice_history of the manifest and their reasons are not "
               "quoted here.")
    if not uc:
        if not m.get("user_choice_history"):
            return "No structure recorded by select_structure in this manifest."
        return ("No current choice in this manifest (select_structure); "
                + ("this run was started from an earlier one; " if from_choice else "") + earlier)
    pk = uc.get("pocket") or {}
    same = uc.get("structure") == pdb and f"{pk.get('comp_id')} {pk.get('label_asym_id')}" == pocket_key
    if not same:
        return (f"Not the structure and pocket recorded by select_structure ({choice.label(uc)}): "
                "the reason the user gave there is not quoted here."
                + (f" This run was started from an earlier choice, since replaced; {earlier}" if from_choice else ""))
    if from_choice and run.get("user_choice_chosen_at") != uc.get("chosen_at"):
        # the same structure and pocket chosen again after the run, or a run that predates the record
        why = ("this run does not record which recording of the choice it used"
               if not run.get("user_choice_chosen_at") else
               f"the choice was recorded again ({uc.get('chosen_at')}) after the one this run used "
               f"({run['user_choice_chosen_at']})")
        return f"Current choice: {choice.label(uc)}. Its reason is not quoted: {why}; {earlier}"
    return (f'Chosen by the user: {choice.label(uc)}; reason recorded by select_structure: "{uc.get("reason")}" '
            f'(recorded {uc.get("chosen_at")})'
            + ("" if from_choice else "; this run was started with explicit arguments"))


def _provenance(doc, s: dict, block: dict, run: dict, m: dict, pocket_key: str):
    f = s.get("file") or {}
    lig, rec = run.get("ligand") or {}, run.get("receptor") or {}
    res = s.get("resolution_angstrom")
    items = [
        [Link(s.get("pdb_id"), (s.get("links") or {}).get("rcsb")),
         f" · {s.get('method') or NO_RECORD} · resolution "
         f"{(_num(res) + ' Å') if res is not None else NO_RECORD} · {s.get('title') or NO_RECORD}"],
        f"mmCIF {f.get('path') or NO_RECORD}, SHA-256 {block.get('raw_sha256') or NO_RECORD}"
        + (f" (retrieved {f['retrieved_at']})" if f.get("retrieved_at") else ""),
        f"Receptor {rec.get('source') or NO_RECORD}, chains "
        + (", ".join(rec.get("chains_auth") or []) or NO_RECORD)
        + (f"; {rec['chains_note']}" if rec.get("chains_note") else ""),
        f"Pocket {pocket_key}: box {_boxtext(block.get('box') or {})}; frame: "
        + ((block.get("box") or {}).get("frame") or NO_RECORD),
        (f"Crystal ligand pose {lig.get('sdf') or NO_RECORD}, SHA-256 {lig.get('sdf_sha256') or NO_RECORD}"
         if "compound_docking" not in run else
         "Crystal ligand pose: not applicable, this run docked a compound of the researcher's selection and "
         "not the crystal ligand of the pocket; the box and the receptor are the pocket's, unchanged"
         + (f". A reference pose of the DOCKED compound does exist, in "
            f"{(run['cross_docking'].get('reference_pose') or {}).get('from') or NO_RECORD}, and is "
            "used for the RMSD below" if run.get("cross_docking") else "")),
        _choice_line(m, s.get("pdb_id"), pocket_key, run),
    ]
    mut = report._mutation_cell(s.get("construct") or {})
    items.append("Construct mutations recorded by get_structures: "
                 + ("; ".join(mut) if isinstance(mut, Lines) else str(mut)))
    doc.bullets(items)


def _ligand_block(doc, run: dict):
    lig = run.get("ligand") or {}
    name = lig.get("name")
    syn = ", ".join(lig.get("synonyms") or []) or "no wwPDB synonym"
    # step 6b: a compound of the user's selection is not an instance of this crystal, so the
    # instance and chain of a pocket ligand do not exist for it and must not print as "no record"
    first = (f"{lig.get('comp_id') or NO_RECORD} (instance {_v(lig.get('label_asym_id'))}, auth chain "
             f"{_v(lig.get('auth_asym_id'))}); wwPDB synonyms: {syn}")
    if "compound_docking" in run:
        first = (f"{_v(lig.get('compound_id'))} ({_v(lig.get('name'))}), level {_v(lig.get('level'))}: a compound "
                 f"of the researcher's selection, not a ligand of this crystal, so it has no instance, no auth "
                 f"chain and no wwPDB synonym. InChIKey {_v(lig.get('inchikey'))}; selection recorded "
                 f"{_v(lig.get('selection_recorded_at'))}")
    items = [
        first,
        f"Chemical name: {report._clip(name, 300) if name else NO_RECORD}",
        f"SMILES: {lig.get('smiles') or NO_RECORD}; heavy atoms: {_v(lig.get('heavy_atoms'))}; "
        f"rotatable torsions: {_v(lig.get('rotatable_torsions'))}",
        f"Source: {lig.get('source') or NO_RECORD}; conformer "
        + (f"'{lig['conformer']}'" if lig.get("conformer") else NO_RECORD)
        + (f"; {lig['method']}" if lig.get("method") else ""),
    ]
    if lig.get("translation_angstrom"):
        items.append("Translation applied to the conformer: ("
                     + ", ".join(_num(x) or "?" for x in lig["translation_angstrom"]) + ") Å")
    if lig.get("note"):
        items.append(lig["note"])
    doc.h(4, "Ligand")
    doc.bullets(items)


def _receptor_items(run: dict) -> list:
    """Receptor preparation and the pre-flight checks, as reproducibility lines.

    These used to open each run with two headings of their own. They say how the receptor was
    protonated and that the ligand really was inside the box before Vina started -- worth keeping,
    never worth reading first (user request 2026-11-18).
    """
    rec = run.get("receptor") or {}
    items = [f"Meeko mk_prepare_receptor {rec.get('options') or NO_RECORD}; exit code {_v(rec.get('exit_code'))}"
             + (f"; {_num(rec['seconds'])} s" if rec.get("seconds") is not None else ""),
             "Residues matching no Meeko template: "
             + (", ".join(rec.get("residues_without_template") or []) or "none reported"),
             "Deleted by Meeko (outside the box, so the pocket is unchanged): "
             + (", ".join(f"{d['residue']} ({_num(d['outside_box_angstrom'])} Å outside the box)"
                          for d in rec.get("deleted") or []) or "none"),
             "Refused by Meeko (too close to the box, not deleted): "
             + (", ".join(f"{d['residue']} ({_num(d['distance_to_box_edge_angstrom'])} Å from the box edge)"
                          for d in rec.get("refused") or []) or "none"),
             f"Residue identifiers: {rec.get('residue_ids') or NO_RECORD}"]
    if rec.get("warning"):
        items.append(f"Warning: {rec['warning']}")
    guard = run.get("guard") or {}
    if guard:
        items += [f"Checked before Vina started: ligand atoms inside the box "
                  f"{guard.get('ligand_atoms_in_box') or NO_RECORD}; ligand centroid inside the box "
                  f"{_v(guard.get('ligand_centroid_in_box'))}; receptor atoms within 8 Å of the box "
                  f"centre {_v(guard.get('receptor_atoms_within_8A_of_box_centre'))}, closest at "
                  f"{_num(guard.get('closest_receptor_atom_to_box_centre_angstrom')) or NO_RECORD} Å"]
        if guard.get("note"):
            items.append(guard["note"])
    return items


def _vina_items(run: dict) -> list:
    v = run.get("vina") or {}
    params = ", ".join(f"{k} {v.get(k)}" for k in ("exhaustiveness", "seed", "num_modes", "energy_range",
                                                   "min_rmsd", "scoring") if k in v)
    items = [f"Vina {v.get('version') or NO_RECORD} (executable SHA-256 {v.get('executable_sha256') or NO_RECORD})",
             f"Parameters: {params or NO_RECORD}; CPUs: {v.get('cpu') or NO_RECORD}",
             f"Search box: {_boxtext(v.get('box') or {})}",
             f"Wall-clock time: {_num(v.get('wall_seconds')) or NO_RECORD} s"]
    if v.get("cpu_seconds_note"):
        items.append(f"CPU time: not measured. {v['cpu_seconds_note']}")
    elif v.get("cpu_seconds") is not None:
        items.append(f"CPU time: {_num(v['cpu_seconds'])} s on {_v(v.get('logical_cpus'))} logical CPUs "
                     f"({_num(v.get('cpu_per_wall_second')) or NO_RECORD} s of CPU per second of wall clock); "
                     "a measurement of this machine, not of the docking")
    return items


def _repro(audit: list, title: str, items: list):
    """Collect one audit entry. It is NOT printed here: section 6 prints them all at the end.

    Nothing measured is dropped from this document -- a hash, a Meeko option or a guard check is
    what lets someone else repeat the run, and removing it would break the audit. But none of it
    helps anyone INTERPRET a result, so the whole audit trail leaves the reading path and becomes
    an appendix (user request 2026-11-18, translated from Spanish: "everything that serves to
    audit the experiments goes to the end... at the beginning, only the data the researcher needs
    in order to interpret"). The entry keeps the heading it had, so a reader who wants the trail of one
    run finds it under that run's name.
    """
    kept = [i for i in items if i]
    if kept:
        audit.append({"title": title, "items": kept})


def _score_span(run: dict) -> str:
    """How far apart the poses of ONE run score. Never a comparison between runs."""
    scores = [p.get("vina_score_kcal_mol") for p in run.get("poses") or []
              if p.get("vina_score_kcal_mol") is not None]
    if not scores:
        return "no pose carries a score"
    if len(scores) == 1:
        return f"a single pose at {_num(min(scores), '{:.2f}')} kcal/mol"
    return (f"{len(scores)} poses spanning {_num(min(scores), '{:.2f}')} to "
            f"{_num(max(scores), '{:.2f}')} kcal/mol")


def _rmsd_sentence(xd: dict) -> str:
    """The T1 anchor in prose: the number, the bar, and whether it settles anything."""
    rng = (xd.get("within_criterion_with_bar") or {}).get("range_angstrom") or []
    decides = (xd.get("within_criterion_with_bar") or {}).get("decides")
    ref = (xd.get("reference_pose") or {}).get("from") or NO_RECORD
    crit = _num(xd.get("criterion_angstrom")) or NO_RECORD
    bar = _num((xd.get("reference_pose_uncertainty") or {}).get("shift_angstrom")) or NO_RECORD
    top = _num(xd.get("top_pose_rmsd_angstrom"), "{:.4f}") or NO_RECORD
    out = (f"Its best-scoring pose sits {top} Å from the crystallographic pose of this same compound "
           f"in {ref}, brought into this receptor by superposing the two proteins. The criterion "
           f"written down in advance is {crit} Å. ")
    if len(rng) == 2 and decides is False:
        out += (f"But that superposition moves the reference pose by ±{bar} Å, so the measurement "
                f"really spans {_num(rng[0])}-{_num(rng[1])} Å, which falls on BOTH sides of the "
                "criterion. This run therefore does NOT decide whether the pose was reproduced, and "
                "must not be reported as a pass.")
    elif len(rng) == 2:
        out += (f"With the ±{bar} Å uncertainty of that superposition the measurement spans "
                f"{_num(rng[0])}-{_num(rng[1])} Å, which stays on one side of the criterion, so it "
                "does settle the question.")
    else:
        out += "The uncertainty of that superposition is not recorded, so the criterion settles nothing."
    return out


def _pose_summary(doc, run: dict):
    """One paragraph a reader can act on, before any table.

    Replaces a list of disconnected bullets (user request 2026-11-18). Every figure in it is
    quoted from the manifest; the paragraph adds sentences, never numbers.
    """
    lig = run.get("ligand") or {}
    cd = run.get("compound_docking") or {}
    rd = run.get("redocking") or {}
    xd = run.get("cross_docking") or {}
    if cd:
        who = f"{_v(lig.get('name'))} ({_v(lig.get('compound_id'))})"
        le = _num(cd.get("ligand_efficiency_kcal_mol_per_heavy_atom"), "{:.4f}") or NO_RECORD
        text = (f"{who}, a compound the researcher selected, was docked into this pocket's box with "
                f"the box, the receptor and the Vina settings left exactly as they were for every "
                f"other compound of the set. Vina returned {_score_span(run)}; the best of them scores "
                f"{_num(cd.get('top_pose_score_kcal_mol'), '{:.3f}') or NO_RECORD} kcal/mol. Over its "
                f"{_v(cd.get('heavy_atoms'))} heavy atoms that is {le} kcal/mol per atom. It is "
                "printed because the raw score grows with the number of atoms, so a bigger molecule "
                "scores better for being bigger; it is NOT used to re-order the set and no threshold "
                "says which of two values is better, because this project has measured none. ")
        text += _rmsd_sentence(xd) if xd else (
            "There is no RMSD for this run, and that is a property of the experiment rather than a "
            "missing value: " + (cd.get("rmsd_not_applicable") or NO_RECORD)
            + " The score is a hypothesis about where this compound would sit, not a reproduced pose.")
        if cd.get("box_limited"):
            fit = cd.get("box_fit") or {}
            text += (" This run is marked box_limited: the compound clears the tightest axis of the box "
                     f"by only {_num(fit.get('min_clearance_angstrom')) or NO_RECORD} Å, so Vina never "
                     "explored the orientations it explored for the others. Its score is not comparable "
                     "with theirs, and it is not ranked with them.")
        doc.p(text)
        return
    crit = _num(rd.get("criterion_angstrom")) or NO_RECORD
    top = _num(rd.get("top_pose_rmsd_angstrom"), "{:.4f}") or NO_RECORD
    text = (f"The crystal ligand {_v(lig.get('comp_id'))} of this pocket was re-docked into its own "
            f"receptor, which is a test of the set-up and not of the method: the answer is already in "
            f"the crystal. Vina returned {_score_span(run)}. ")
    verdict = rd.get("within_criterion")
    if verdict is None:
        # a missing field is a missing field. Reading it as False would turn "the manifest does
        # not say" into "this run failed", which is the one substitution this project forbids.
        text += (f"The best-scoring pose lands {top} Å from the crystallographic pose, against a "
                 f"criterion of {crit} Å. Whether it is within the criterion: {NO_RECORD} -- this "
                 "manifest does not record the verdict, which is not the same as failing it.")
    else:
        text += (f"The best-scoring pose lands {top} Å from the crystallographic pose against a "
                 f"criterion of {crit} Å, so this set-up "
                 + ("reproduces" if verdict else "does NOT reproduce")
                 + " a pose it already knew. That says the box, the receptor and the scoring run as "
                   "intended; it says nothing about a compound the crystal does not contain.")
    doc.p(text)


def _cross_docking_items(xd: dict) -> list:
    """The T1 RMSD, and never the number on its own: the bar, the protein fit it came from, and
    whether the criterion decides at all (ESTADO 7.4)."""
    sup = xd.get("superposition") or {}
    unc = xd.get("reference_pose_uncertainty") or {}
    bar = unc.get("shift_angstrom")
    rng = (xd.get("within_criterion_with_bar") or {}).get("range_angstrom") or []
    nc = xd.get("rmsd_negative_control") or {}
    geo = xd.get("reference_pose_geometry") or {}
    ref = xd.get("reference_pose") or {}
    decides = (xd.get("within_criterion_with_bar") or {}).get("decides")
    return [
        f"Top-pose RMSD: {_num(xd.get('top_pose_rmsd_angstrom'), '{:.4f}') or NO_RECORD} Å "
        f"± {_num(bar) or NO_RECORD} Å (the error bar of the reference pose, not of the docking); "
        f"criterion <= {_num(xd.get('criterion_angstrom')) or NO_RECORD} Å; within the criterion: "
        f"{_v(xd.get('within_criterion'))}",
        (f"With the bar the RMSD spans {_num(rng[0]) or NO_RECORD}-{_num(rng[1]) or NO_RECORD} Å, "
         + ("which falls on BOTH sides of the criterion, so THIS RUN DOES NOT DECIDE whether the "
            "pose was reproduced. It must not be reported as a pass."
            if decides is False else "which falls on one side of the criterion, so it decides.")
         if len(rng) == 2 else f"Bar: {NO_RECORD}"),
        f"Reference pose: the crystallographic pose of this same compound in "
        f"{ref.get('from') or NO_RECORD}, matched by {ref.get('matched_by') or NO_RECORD}",
        f"Protein superposition: {sup.get('method') or NO_RECORD}; {sup.get('mobile') or NO_RECORD} "
        f"on {sup.get('target') or NO_RECORD}, RMSD {_num(sup.get('rmsd_angstrom')) or NO_RECORD} Å "
        f"over {_v(sup.get('pairs'))} CA pairs",
        f"Error bar: {unc.get('why') or NO_RECORD}"
        + ("" if unc.get("measured_here") else " (measured here: no)"),
        unc.get("note") or NO_RECORD,
        f"Where the reference pose lands in this receptor, checked by geometry and not by name: "
        f"{geo.get('heavy_atoms_in_box') or NO_RECORD} heavy atoms inside the box, centroid "
        f"{_num(geo.get('centroid_to_box_centre_angstrom')) or NO_RECORD} Å from the box centre and "
        f"{_num(geo.get('centroid_to_pocket_ligand_centroid_angstrom')) or NO_RECORD} Å from the "
        f"centroid of this pocket's own ligand",
        f"RMSD method: {xd.get('rmsd_method') or NO_RECORD}",
        "Negative control, measured with the same function as the poses: the reference pose "
        f"against itself gave {_num(nc.get('crystal_vs_itself_angstrom')) or NO_RECORD} Å and moved "
        f"5 Å gave {_num(nc.get('crystal_moved_5A_angstrom')) or NO_RECORD} Å; passed: "
        f"{_v(nc.get('passed'))}",
        xd.get("criterion_note") or NO_RECORD,
        f"Scope of the RMSD: {xd.get('scope') or NO_RECORD}",
    ]


def _compound_poses_block(doc, run: dict, count: str):
    """Step 6b: a compound run has scores and NO RMSD. The difference between "no record" (the
    manifest lost it) and "not applicable" (there is no reference pose to measure against) is the
    whole point, so the RMSD column is not printed at all and the reason is given instead."""
    cd = run.get("compound_docking") or {}
    fit = cd.get("box_fit") or {}
    xd = run.get("cross_docking") or {}
    doc.h(4, "Result")
    _pose_summary(doc, run)
    if xd:
        doc.p(f"{count}. Both columns are calculated. The RMSD is NOT against a pose of this "
              f"structure: it is against the crystallographic pose of this same compound in "
              f"{xd.get('reference_pose', {}).get('from') or NO_RECORD}, carried here by a protein "
              f"superposition, and it travels with the error bar given above.", small=True)
        doc.table(["Mode", "Vina score (kcal/mol, calculated)",
                   "RMSD to the reference pose (Å, calculated)"],
                  [[p.get("mode"), _num(p.get("vina_score_kcal_mol"), "{:.2f}") or NO_RECORD,
                    _num(p.get("rmsd_to_reference_pose_angstrom"), "{:.4f}") or NO_RECORD]
                   for p in run.get("poses") or []], (0.08, 0.3, 0.32))
    else:
        doc.p(f"{count}. The scores are calculated values of this run, listed in the order Vina "
              "returned them. There is no RMSD column, for the reason given above.", small=True)
        doc.table(["Mode", "Vina score (kcal/mol, calculated)"],
                  [[p.get("mode"), _num(p.get("vina_score_kcal_mol"), "{:.2f}") or NO_RECORD]
                   for p in run.get("poses") or []], (0.08, 0.3))
    doc.p(f"Scope of this run: {cd.get('scope') or NO_RECORD}"
          + (" That line is the state of the run AT DOCKING TIME, when no pose of this compound had "
             "been brought into this frame; the reference pose was added afterwards by "
             "cross_dock_rmsd, from another crystal." if xd else ""), small=True)
    items = [
        f"Box fit against the pocket's unchanged box: clearance "
        f"{_num(fit.get('min_clearance_angstrom')) or NO_RECORD} Å in the tightest axis (box axis "
        f"{_num(fit.get('limiting_box_axis_angstrom')) or NO_RECORD} Å, compound extent "
        f"{_num(fit.get('limiting_ligand_extent_angstrom')) or NO_RECORD} Å), declared minimum "
        f"{_num(fit.get('threshold_angstrom')) or NO_RECORD} Å; measured as: {fit.get('method') or NO_RECORD}",
        f"Parameters as pre-registered for the set: {_v(cd.get('parameters_as_preregistered'))}",
        cd.get("ligand_efficiency_note") or NO_RECORD,
    ]
    if not xd:
        items.append(f"Pass/fail criterion: none. {cd.get('not_a_control_reason') or NO_RECORD}")
    if cd.get("box_limited"):
        items.insert(0, cd.get("box_limited_note") or "box_limited")
    if cd.get("note"):
        items.append(cd["note"])
    return (_cross_docking_items(xd) if xd else []) + items


def _poses_block(doc, run: dict):
    v, rd = run.get("vina") or {}, run.get("redocking") or {}
    n, asked = len(run.get("poses") or []), v.get("num_modes")
    count = (f"{n} poses returned (num_modes {asked})" if n == asked else
             f"{n} poses returned; num_modes {_v(asked)}, energy_range {_v(v.get('energy_range'))} kcal/mol"
             + ("; the manifest does not record why fewer were returned" if isinstance(asked, int) and n < asked
                else ""))
    if "compound_docking" in run:
        return _compound_poses_block(doc, run, count)
    doc.h(4, "Result")
    _pose_summary(doc, run)
    doc.p(f"{count}. Both columns are calculated values of this run; the RMSD is against the crystal pose of this "
          f"same pocket (method: {rd.get('rmsd_method') or NO_RECORD}).", small=True)
    doc.table(["Mode", "Vina score (kcal/mol, calculated)", "RMSD to the crystal pose (Å, calculated)"],
              [[p.get("mode"), _num(p.get("vina_score_kcal_mol"), "{:.2f}") or NO_RECORD,
                _num(p.get("rmsd_to_crystal_pose_angstrom"), "{:.4f}") or NO_RECORD]
               for p in run.get("poses") or []], (0.08, 0.3, 0.32))
    nc = rd.get("rmsd_negative_control") or {}
    doc.p(f"Scope of this run: {rd.get('scope') or NO_RECORD}", small=True)
    items = [f"Pre-registered control run: {_v(rd.get('is_preregistered_control'))}",
             f"RMSD method: {rd.get('rmsd_method') or NO_RECORD}",
             "Negative control of this run, measured with the same function as the poses: the crystal pose "
             f"against itself gave {_num(nc.get('crystal_vs_itself_angstrom')) or NO_RECORD} Å and the crystal "
             f"pose moved 5 Å gave {_num(nc.get('crystal_moved_5A_angstrom')) or NO_RECORD} Å; "
             f"passed: {_v(nc.get('passed'))}"]
    if rd.get("note"):
        items.append(rd["note"])
    return items


def _files_items(run: dict) -> list:
    """Each written file with its hash, as one line. The hash is for the audit, not the reader.

    A reader interpreting a docking result never needs a SHA-256 (user question 2026-11-18). It is
    kept because it is what proves the file on disk is the file this run wrote -- the manifest, the
    MD package and analyze_md all check against it -- so it moves to the reproducibility block
    instead of taking a table column in the reading path.
    """
    return [f"{name}: {f.get('path') or NO_RECORD} (SHA-256 {f.get('sha256') or NO_RECORD})"
            for name, f in sorted((run.get("files") or {}).items())]


def _docking(doc, m: dict, accession: str, audit: list):
    doc.h(2, "1. Molecular docking (AutoDock Vina)")
    runs = _run_rows(m)
    structures = m.get("structures") or []
    pockets = sum(1 for s in structures for pk in (s.get("prepared") or {}).get("pockets") or [] if not pk.get("error"))
    if not runs:
        doc.p("Not run. No docking run is recorded in prep/manifest.json. Run dock_ligand(\""
              f"{accession}\", structure, pocket) and then dock_status to record one, and build this record again.")
        return 0
    docked_pockets = {(s.get("pdb_id"), key) for s, key, _b, _r, _run in runs}
    n_comp = sum(1 for *_x, run in runs if "compound_docking" in run)
    n_cross = sum(1 for *_x, run in runs if run.get("cross_docking"))
    kinds = (f"{len(runs) - n_comp} re-docking{'' if len(runs) - n_comp == 1 else 's'} of the crystal ligand of a "
             f"pocket into its own receptor (a reference pose in the same structure, and therefore an RMSD)")
    if n_comp:
        kinds += (f" and {n_comp} run{'' if n_comp == 1 else 's'} of a compound of the researcher's selection into "
                  "a pocket's unchanged box, with no pass/fail criterion; the levels T1 and T3 are not comparable "
                  "with each other and are not ranked together")
        kinds += (f". Of those, {n_cross} {'has' if n_cross == 1 else 'have'} an RMSD, measured against the "
                  "crystallographic pose of the SAME compound in another structure of this target, carried over by "
                  "a protein superposition and reported with the error bar of that superposition; the rest have no "
                  "reference pose anywhere, so for them the RMSD is not applicable, which is not a missing value"
                  if n_cross else
                  ". None has a reference pose, so for all of them the RMSD is not applicable, which is not a "
                  "missing value")
    doc.p(f"{len(runs)} docking run{'' if len(runs) == 1 else 's'} recorded, on {len(docked_pockets)} of the "
          f"{pockets} pockets computed in the "
          f"manifest ({len({s.get('pdb_id') for s, *_ in runs})} of {len(structures)} structures): {kinds}. One "
          "ligand per run. " + VOCABULARY, small=True)
    doc.table(["PDB", "Pocket", "Run", "Recorded (UTC)", "Written by", "Started from"],
              [[s.get("pdb_id"), key, run_id, _v(run.get("computed_at")), run.get("generated_by") or NO_RECORD,
                run.get("selected_by") or NO_RECORD] for s, key, _b, run_id, run in runs],
              (0.08, 0.1, 0.22, 0.2, 0.2, 0.2))
    for s, key, block, run_id, run in runs:
        doc.h(3, f"{s.get('pdb_id')} · pocket {key} · run {run_id}")
        if run.get("replaced_run_computed_at"):
            doc.p(f"Replaced an earlier run of the same pocket and conformer, recorded "
                  f"{run['replaced_run_computed_at']}; only this one is kept in the manifest.", small=True)
        # the result first, then what it was run on, then the audit trail. The reader who only
        # wants the answer stops after the first paragraph (user request 2026-11-18).
        extra = _poses_block(doc, run)
        _ligand_block(doc, run)
        doc.h(4, "Provenance")
        _provenance(doc, s, block, run, m, key)
        _repro(audit, f"Docking · {s.get('pdb_id')} · pocket {key} · run {run_id}",
               _vina_items(run) + _receptor_items(run) + (extra or []) + _files_items(run))
        if run.get("note"):
            doc.p(f"{run.get('values') or 'calculated'}: {run['note']}", small=True)
    return len(runs)


def _admet_rows(m: dict) -> list:
    """(structure, pocket_key, admet block) for every pocket with a recorded ADMET result."""
    out = []
    for s in m.get("structures") or []:
        for key, block in sorted((s.get("admet") or {}).items()):
            if isinstance(block, dict):
                out.append((s, key, block))
    return out


def _admet(doc, m: dict, accession: str, audit: list) -> int:
    doc.h(2, "3. ADMET")
    rows = _admet_rows(m)
    if not rows:
        doc.p("Not run. No ADMET result is recorded in prep/manifest.json. Run analyze_admet(\""
              f"{accession}\", structure, pocket) to record one, and build this record again.")
        return 0
    doc.p(f"{len(rows)} compound{'' if len(rows) == 1 else 's'} recorded. " + ADMET_VOCABULARY,
          small=True)
    for s, key, block in rows:
        ident = block.get("compound_identity") or {}
        look = block.get("chembl_lookup") or {}
        # A compound of the researcher's selection has no comp_id: it is not a ligand of any
        # crystal, so that field is NOT APPLICABLE rather than missing. Printing "no record" in
        # the heading said the opposite -- and said it in the most visible place in the section,
        # about compounds whose ADMET was in fact retrieved in full.
        who = ident.get("comp_id") or ident.get("name") or ident.get("compound_id")
        doc.h(3, f"{s.get('pdb_id')} · pocket {key} · {who or NO_RECORD}")
        doc.p(f"Recorded {_v(block.get('generated_at'))} by {_v(block.get('generated_by'))}. "
              f"{block.get('selection') or ''}", small=True)

        # The identity block used to open the compound with an eight-row table of files, hashes
        # and run ids. None of it is read when judging a compound (user request 2026-11-18): it
        # answers "is this the same molecule the docking used", which is a yes/no, so the yes/no
        # stays here and the evidence for it moves to the reproducibility block at the end.
        repro = [
            f"InChIKey (computed here with RDKit): {ident.get('inchikey') or NO_RECORD}",
            f"Read from {ident.get('source') or NO_RECORD}, file {ident.get('source_file') or NO_RECORD} "
            f"(SHA-256 {ident.get('source_sha256') or NO_RECORD})",
            f"Docking run: {ident.get('run_id') or 'none (crystallographic ligand)'}",
            f"Agreement between the files: {ident.get('identity_check') or NO_RECORD}",
            f"Heavy atoms, SDF against the manifest's pocket: "
            f"{(block.get('calculated_properties') or {}).get('rdkit', {}).get('heavy_atoms', NO_RECORD)}"
            f" / {_v(ident.get('heavy_atoms_in_manifest_pocket'))}",
            f"ChEMBL molecules with this exact InChIKey: {_v(look.get('matches'))}",
        ]

        rec = block.get("chembl_record") or {}
        doc.h(4, "What ChEMBL holds on this compound")
        # "no record" and "the lookup never completed" are different statements and only the
        # first is a fact about ChEMBL. Printing "No record." for a failed lookup was the report
        # saying zero where the tool said UNKNOWN (reviewer finding 2026-11-12; real case:
        # 5JSM pocket 6NB J, whose lookup returned HTTP 500).
        if str(rec.get("status", "")).startswith("UNKNOWN"):
            doc.p(f"UNKNOWN, not 'no record'. {rec.get('note') or ''} Whether this compound has a "
                  "ChEMBL entry, a name, a development phase or any measurement was NOT "
                  "established: the lookup did not complete. "
                  + (rec.get("name_status") or "")
                  + " Only the calculated properties and the rule reading below apply.",
                  small=True)
        elif rec.get("status") == "NO RECORD" or not rec.get("chembl_id"):
            doc.p(f"No record. {rec.get('note') or 'No ChEMBL molecule has this exact InChIKey.'} "
                  + (rec.get("name_status") or "")
                  + " Only the calculated properties and the rule reading below apply.",
                  small=True)
        else:
            fields = rec.get("fields") or {}
            minus_one = rec.get("fields_minus_one") or []
            false_flags = rec.get("boolean_fields_false") or []
            # The eleven-field table is an entry dump; what a reader acts on is how far this
            # compound got clinically and whether anything is flagged against it. That becomes a
            # sentence, and the full table moves to the reproducibility block.
            phase, approved = fields.get("max_phase"), fields.get("first_approval")
            how_far = ("an approved drug" if phase is not None and float(phase or 0) >= 4
                       else f"in clinical phase {_v(phase)}" if phase not in (None, "", -1)
                       else "not recorded as having reached any clinical phase")
            flags = [name for name, key in (("withdrawn somewhere", "withdrawn_flag"),
                                            ("carries a black-box warning", "black_box_warning"))
                     if fields.get(key) in (True, 1, "1")]
            atc = ", ".join(rec.get("atc_classifications") or [])
            line = Lines([Link(f"ChEMBL {rec['chembl_id']}", rec.get("link")),
                          f" is {rec.get('pref_name') or 'unnamed in ChEMBL (no pref_name)'}, {how_far}"
                          + (f", first approved {approved}" if approved not in (None, "", -1) else "")
                          + ". "
                          + ("Flagged: " + " and ".join(flags) + ". " if flags else
                             "Nothing is flagged against it in the two safety fields ChEMBL carries "
                             "(withdrawn, black-box warning), which is the absence of a flag and not "
                             "a safety statement. ")
                          + (f"ATC: {atc}." if atc else "No ATC classification.")])
            doc.p(line)
            if minus_one or false_flags:
                doc.p("Some of those fields come back as -1 or False on an entry with nothing "
                      "annotated, so neither value means the compound lacks the property. "
                      + (rec.get("defaults_note") or ""), small=True)
            repro.append("ChEMBL entry fields: " + "; ".join(
                f"{f} {_v(fields.get(f))}"
                + (" (-1, a default)" if f in minus_one else " (False, a default)" if f in false_flags else "")
                for f in ("max_phase", "first_approval", "withdrawn_flag", "black_box_warning",
                          "oral", "parenteral", "topical", "molecule_type", "therapeutic_flag",
                          "usan_stem", "first_in_class")))

        meas = block.get("measured_adme") or {}
        doc.h(4, "Measured ADME (retrieved from ChEMBL)")
        if not meas.get("retrieved"):
            # "the lookup did not complete" was printed even when no query was made at all
            if meas.get("no_chembl_entry"):
                doc.p("No ADME measurement is reported because no ChEMBL molecule has this exact "
                      "InChIKey, so no ADME query was made. ChEMBL holds nothing for a compound it "
                      "does not have; a measurement published elsewhere would not be here.",
                      small=True)
            else:
                # the reason travels in the note: "the lookup failed" and "the time budget ran
                # out before any assay type could be asked about" are different, and a document
                # that prints neither leaves a zero looking like a real zero (reviewer, 11-12)
                doc.p("UNKNOWN, not 'no data' and not zero rows: "
                      + (meas.get("note") or "the ChEMBL lookup did not complete")
                      + (f" ({meas['error']})" if meas.get("error") else ""), small=True)
        else:
            adme = meas.get("adme_rows") or []
            read = meas.get("rows_read")
            held = meas.get("rows_total")
            doc.p(f"{len(adme)} of {_v(read)} rows read from ChEMBL assay types "
                  f"{', '.join(meas.get('assay_types_queried') or {})} are on this project's ADME "
                  f"allowlist. Assay types {', '.join(meas.get('assay_types_not_queried') or {})} "
                  "were not queried, so a measurement filed only there is not here. Values are "
                  "quoted with their own units and are never converted or averaged; a relation "
                  "('>' or '<') is part of the value.", small=True)
            if meas.get("truncated"):
                # the cap was blamed for everything, including a clock that stopped the reading
                # with the cap nowhere near (reviewer 11-12): the reason has to be the real one
                skipped = meas.get("assay_types_skipped_at_the_time_budget") or []
                cut = meas.get("assay_types_cut_off_at_the_time_budget") or []
                capped = meas.get("assay_types_skipped_at_the_row_cap") or []
                reason = ("the row cap of " + str(_v(meas.get("row_cap"))) + " was reached"
                          if (read or 0) >= (meas.get("row_cap") or 10 ** 9)
                          else "the time budget ran out, not the row cap: ChEMBL was slow, "
                               "not empty" if skipped or cut
                          else "reading stopped before the total ChEMBL reported was reached")
                doc.p(f"TRUNCATED: ChEMBL holds at least {_v(held)} rows in the assay types "
                      f"queried and only {_v(read)} were read, so {_v(meas.get('rows_not_read'))} "
                      f"were never examined. Reason: {reason}. Every count in this subsection is "
                      "a LOWER BOUND."
                      + (f" Assay type(s) {', '.join(skipped)} were never queried at all, so how "
                         "many rows they hold is UNKNOWN, not zero." if skipped else "")
                      + (f" Assay type(s) {', '.join(cut)} were cut off part-way." if cut else "")
                      + (f" Assay type(s) {', '.join(capped)} were not reached before the row cap "
                         "was spent." if capped else ""), small=True)
            if adme:
                doc.table(["Property", "Type", "Value", "Units", "Assay"],
                          [[r.get("adme_class") or NO_RECORD, r.get("standard_type") or NO_RECORD,
                            f"{r.get('relation') or ''} {_v(r.get('value'))}".strip()
                            + (" [depositor says calculated]"
                               if r.get("description_says_value_was_calculated") else ""),
                            r.get("units") or NO_RECORD, r.get("assay_description") or NO_RECORD]
                           for r in adme],
                          (0.14, 0.1, 0.14, 0.14, 0.48))
                if any(r.get("description_says_value_was_calculated") for r in adme):
                    doc.p("A row marked [depositor says calculated] has an assay description that "
                          "says the value was calculated or estimated rather than measured. The "
                          "wording is the source's own; such a row is not a measurement.",
                          small=True)
            else:
                doc.p("No row in the assay types queried is on this project's ADME allowlist.",
                      small=True)
            other = meas.get("rows_excluded_by_allowlist_by_standard_type") or {}
            if other:
                # these rows are NOT "not ADME": for vemurafenib the excluded 'Ratio' rows are
                # P-gp and BCRP efflux ratios and one 'Activity' row is CYP450 induction, all of
                # which are ADME. Calling them "not ADME" was a claim about the source's data
                # made from this project's own narrower allowlist (reviewer finding, 2026-11-11).
                doc.p("Rows in those assay types that are NOT on this project's allowlist and were "
                      "therefore excluded: "
                      + ", ".join(f"{k or NO_RECORD} ({v})" for k, v in other.items())
                      + ". This does not mean they are not ADME. The allowlist is narrower than "
                        "ADME on purpose: the generic types 'Activity' and 'Ratio' carry real ADME "
                        "in some rows (efflux ratios, CYP induction) and cytotoxicity or "
                        "target-engagement in others, and the two cannot be separated without "
                        "reading each assay description. An excluded row may hold a measurement "
                        "this subsection does not show.", small=True)
                # F3: the counts alone said WHICH TYPES were excluded but not WHAT they were, so
                # the rows that prompted the paragraph above -- the P-gp and BCRP efflux ratios
                # and the CYP450 induction -- stayed invisible to the reader it warns. They are
                # listed here unclassified, under a heading that says so, and they are NOT added
                # to the table above.
                shown = meas.get("rows_excluded_by_allowlist") or []
                if shown:
                    doc.p(f"The excluded rows themselves ({meas.get('rows_excluded_shown')}), "
                          "unclassified and NOT part of the count above, so each can be read on "
                          "its own terms:", small=True)
                    doc.table(
                        ["Type", "Value", "Units", "Assay (ChEMBL's own description)"],
                        [[r.get("standard_type") or NO_RECORD,
                          " ".join(x for x in (r.get("relation"), _v(r.get("value"))) if x),
                          r.get("units") or "—",
                          r.get("assay_description") or NO_RECORD] for r in shown],
                        # FRACTIONS of the usable width, like every other table in this project
                        # (report.py multiplies them by it). Absolute points here squeezed the
                        # whole table into 40% of the frame and cut off the very descriptions
                        # this table exists to show (reviewer finding 2026-11-12).
                        (0.14, 0.12, 0.12, 0.62))

        calc = block.get("calculated_properties") or {}
        rd, ch = calc.get("rdkit") or {}, calc.get("chembl") or {}
        doc.h(4, "Calculated properties (not measured, not predicted)")
        doc.p(calc.get("note") or "", small=True)
        doc.table(["Property", f"This project ({rd.get('engine') or 'RDKit'})", "ChEMBL's own engine"],
                  [["Molecular weight", _v(rd.get("molecular_weight")), _v(ch.get("full_mwt"))],
                   ["cLogP", _v(rd.get("clogp")), _v(ch.get("alogp"))],
                   ["TPSA", _v(rd.get("tpsa")), _v(ch.get("psa"))],
                   ["H-bond donors", _v(rd.get("hbd")), _v(ch.get("hbd"))],
                   ["H-bond acceptors", _v(rd.get("hba")), _v(ch.get("hba"))],
                   ["Rotatable bonds", _v(rd.get("rotatable_bonds")), _v(ch.get("rtb"))],
                   ["Aromatic rings", _v(rd.get("aromatic_rings")), _v(ch.get("aromatic_rings"))],
                   ["Heavy atoms", _v(rd.get("heavy_atoms")), _v(ch.get("heavy_atoms"))],
                   ["Fraction Csp3", _v(rd.get("fraction_csp3")), NO_RECORD],
                   ["QED", _v(rd.get("qed")), _v(ch.get("qed_weighted"))],
                   ["Formula", rd.get("molecular_formula") or NO_RECORD,
                    ch.get("full_molformula") or NO_RECORD]],
                  (0.3, 0.35, 0.35))

        reading = block.get("rule_reading") or {}
        res = reading.get("results") or {}
        doc.h(4, "Property rules (offered, not measured)")
        doc.p(f"{reading.get('status') or ''}. {reading.get('disclaimer') or ''}", small=True)
        doc.table(["Rule", "Cut-offs", "Violations", "Exceeded"],
                  [[name,
                    "; ".join((reading.get("rules_used") or {}).get(name, {}).get("cutoffs") or []),
                    "UNKNOWN" if r.get("violations") is None else str(r.get("violations")),
                    "; ".join(r.get("violated") or []) or
                    (f"not evaluated: {', '.join(r.get('not_evaluated') or [])}"
                     if r.get("not_evaluated") else "none")]
                   for name, r in res.items()],
                  (0.12, 0.38, 0.12, 0.38))
        for note in block.get("notes") or []:
            doc.p(note, small=True)
        _repro(audit, f"ADMET · {s.get('pdb_id')} · pocket {key} · {who or NO_RECORD}", repro)
    return len(rows)


def _md_rows(m: dict) -> list:
    """(structure, pocket_key, md_package, md_job, md_analysis) for every pocket with any MD block."""
    out = []
    for s in m.get("structures") or []:
        blocks = [s.get(k) if isinstance(s.get(k), dict) else {}
                  for k in ("md_package", "md_job", "md_analysis")]
        for key in sorted(set().union(*blocks)):
            out.append((s, key, *[b.get(key) if isinstance(b.get(key), dict) else {}
                                  for b in blocks]))
    return out


def _md_runs(m: dict) -> list:
    """(structure, pocket_key, stage_label, record) for every recorded analysis, in stage order."""
    return [(s, key, label, ana[label]) for s, key, _p, _j, ana in _md_rows(m)
            for label in MD_STAGES if isinstance(ana.get(label), dict)]


def _compound_name(m: dict, pocket_key: str) -> str:
    """The researcher's name for the compound a pocket key ends in, from ligand_choice.

    A key like 'TXV U CHEMBL1229517' says which compound was simulated only to someone who has
    memorised the ChEMBL ids. The name is in the manifest; it is looked up and never invented.
    """
    for c in ((m.get("ligand_choice") or {}).get("compounds") or []):
        cid = str(c.get("id") or "")
        if cid and pocket_key.upper().endswith(cid.upper()):
            return f"{c.get('name') or cid} ({cid}" + (f", level {c['level']}" if c.get("level") else "") + ")"
    return ""


def _md_subject(m: dict, pocket_key: str) -> str:
    """What was actually simulated here, in words: a selected compound, or the crystal ligand."""
    name = _compound_name(m, pocket_key)
    if name:
        return f"{name}, on its DOCKED pose -- no crystal ever held this compound in this pocket"
    return f"{pocket_key.split()[0]}, the crystallographic ligand of this pocket, on its re-docked pose"


def _site_value(site: dict, name: str):
    field = site.get(name) or {}
    return field.get("gres") or field.get("value") if isinstance(field, dict) else None


def _md_config(doc, package: dict, job: dict, audit: list, where: str):
    """Subtask 5.1: the configuration as recorded by build_md_package and build_md_job, quoted."""
    rows, rows_audit = [], []
    if package:
        pose = package.get("pose") or {}
        rmsd = pose.get("rmsd_to_crystal_pose_angstrom")
        rows.append(["Which pose is being simulated",
                     f"pose {_v(pose.get('index'))} (Vina mode {_v(pose.get('mode'))}) of docking run "
                     f"{_v(package.get('run_id'))}, scoring {_v(pose.get('vina_score_kcal_mol'))} kcal/mol"
                     + (f"; it sits {_v(rmsd)} Å from the crystal pose" if rmsd is not None else
                        "; it has no RMSD, because no crystal holds this compound in this pocket")])
        gaps = package.get("to_model")
        rows.append(["Parts of the protein CHARMM-GUI was asked to build",
                     ", ".join(f"chain {g.get('chain')} {g.get('first')}-{g.get('last')} ({g.get('kind')})"
                               for g in gaps) if gaps else
                     ("none: nothing was missing" if gaps == [] else NO_RECORD)])
        rows_audit.append(f"Gaps read from {_v(package.get('missing_residues_source'))}")
    else:
        rows.append(["Package for CHARMM-GUI", "none recorded for this pocket"])
    if job:
        checks = (job.get("verified") or {}).get("checks") or {}
        closest = checks.get("closest_modelled_run_to_ligand") or {}
        built = checks.get("modelled_runs_present")
        rows.append(["Parts it actually built, in the system it returned",
                     (", ".join(built) if built else ("none" if built == [] else NO_RECORD))
                     + (f"; the closest of them, {closest.get('run')}, sits "
                        f"{_v(closest.get('min_distance_angstrom'))} Å from the ligand" if closest else "")])
        sysm = checks.get("system") or {}
        ions = sysm.get("ions") or {}
        rows.append(["The box that was simulated",
                     f"{_v(sysm.get('box_side_angstrom'))} Å per side, {_v(sysm.get('waters'))} waters, "
                     + ("ions " + ", ".join(f"{k} {v}" for k, v in ions.items()) if ions else "no ions recorded")
                     + f", at {_v(sysm.get('temperature_kelvin'))} K"])
        inp, slurm = job.get("inputs") or {}, job.get("slurm") or {}
        val = lambda k: _v((inp.get(k) or {}).get("value") if isinstance(inp.get(k), dict)  # noqa: E731
                           else inp.get(k))
        rows.append(["Settling the system before measuring (equilibration)",
                     f"{val('equilibration_minimize')} steps of minimization, then "
                     f"{val('equilibration_run')} steps of restrained dynamics "
                     f"({_v(inp.get('equilibration_input'))})"])
        rows.append(["The run the figures come from (production)",
                     f"{_v(slurm.get('nanoseconds_total'))} ns in total, written as "
                     f"{_v(slurm.get('chunks'))} chunk(s) of {_v(inp.get('nanoseconds_per_chunk'))} ns, "
                     f"with a frame stored every {val('dcdfreq')} steps ({_v(inp.get('production_input'))})"])
        rows.append(["Time step and cutoff",
                     f"{_v(sysm.get('timestep_fs'))} fs per step; non-bonded cutoff "
                     f"{_v(sysm.get('cutoff_angstrom'))} Å"])
        site = job.get("site") or {}
        rows.append(["Where it ran",
                     f"partition {_v(_site_value(site, 'partition'))}, GPUs "
                     f"{_v(_site_value(site, 'gpus'))}, time limit {_v(_site_value(site, 'time'))}, "
                     f"NAMD {_v(_binary_tail(_site_value(site, 'namd_binary')))}"])
        arch, psf = job.get("archive") or {}, job.get("psf") or {}
        edited = job.get("inputs_edited")
        rows_audit += [
            f"CHARMM-GUI job {_v(job.get('job_id'))}, folder {_v(job.get('job_folder'))}; archive "
            f"SHA-256 {_v(arch.get('sha256'))}",
            "Topology step3_input.psf SHA-256 "
            + (psf.get("sha256") or (f"none: {psf['source']}" if psf.get("source") else
                                     "not recorded by the version of build_md_job that wrote this job")),
            "Inputs edited by build_md_job: "
            + ("none" if edited == [] else (", ".join(str(e.get("file") or e) for e in edited)
                                            if edited else NO_RECORD)),
            "Each cluster value carries its own mark and source in the manifest",
        ]
    else:
        rows.append(["Simulation job", "none recorded for this pocket"])
    doc.h(4, "How this simulation was set up")
    doc.table(["", "As written"], rows, (0.34, 0.66))
    doc.p("Every row is quoted from what build_md_package and build_md_job WROTE. Whether the cluster "
          "ran those inputs unchanged is not recorded here: each analysed run below names the NAMD input "
          "it actually read, and its time axis comes from the trajectory rather than from this table.",
          small=True)
    if rows_audit:
        _repro(audit, f"MD set-up · {where}", rows_audit)


def _binary_tail(value):
    """The NAMD binary without the cluster's directory layout: its folder and file name only.

    The report is meant to be shareable, and the site's install path says where the user computes
    without adding anything the reader needs. The folder of a NAMD release carries its version
    (NAMD_3.0.3_Linux-x86_64-multicore-CUDA), so that much stays; the full path stays in the
    manifest. Anything that is not a path (None, a CONFIRM mark) is returned unchanged.
    """
    if not isinstance(value, str) or ("/" not in value and "\\" not in value):
        return value
    parts = [p for p in value.replace("\\", "/").split("/") if p]
    return ".../" + "/".join(parts[-2:]) + " (full path in the manifest)"


def _control_line(res: dict) -> str:
    """The equilibration control, quoted from the recorded reading. Never re-derived here."""
    reading = res.get("reading") or {}
    c = reading.get("equilibration_control")
    if not isinstance(c, dict):
        return ("EQUILIBRATION CONTROL: NOT RECORDED. This analysis was recorded by a version of "
                "analyze_md that did not write the control as a field; read the full-segment "
                "backbone rows below before any other figure.")
    head = (f"the full backbone of segment {_v(c.get('segid'))} has an absolute drift of "
            f"{_v(c.get('drift_angstrom'))} Å (|last-quarter mean minus first-quarter mean|; signed "
            f"{_v(c.get('drift_signed_angstrom'))} Å) against a threshold of "
            f"{_v(c.get('threshold_angstrom'))} Å. The threshold is a convention, not a "
            "measurement, written down before any production trajectory existed but after the "
            "restrained equilibration of the same system had been measured")
    status = c.get("status")
    if status == "NOT MET":
        return f"EQUILIBRATION CONTROL: NOT MET. {head[0].upper()}{head[1:]}. Every figure below is premature."
    if status == "NO EVIDENCE":
        return (f"EQUILIBRATION CONTROL: NO EVIDENCE. {head[0].upper()}{head[1:]}, but {c.get('why')}. "
                "A figure below that concerns held atoms carries no evidence either.")
    if status == "MET":
        return (f"EQUILIBRATION CONTROL: MET. {head[0].upper()}{head[1:]}. This says the backbone "
                "settled within that convention; it says nothing about whether the pose is correct.")
    return f"EQUILIBRATION CONTROL: {_v(status)}. {_v(c.get('why'))}"


def _num4(x):
    return _num(x, "{:.4f}") or NO_RECORD


def _metric_row(name, series: dict, selection=None):
    s = series or {}
    return [name, selection or s.get("selection") or NO_RECORD,
            *[_num4(s.get(k)) for k in ("min", "max", "mean", "first_quarter_mean",
                                        "last_quarter_mean", "drift", "last_quarter_range")]]


def _series_file(prep: Path, rec: dict):
    """(full result, None) from the recorded series file, or (None, reason) when it cannot be used."""
    sf = rec.get("series_file") or {}
    if not sf.get("path") or not sf.get("sha256"):
        return None, "no series file is recorded for this run"
    path = (prep / sf["path"]).resolve()
    if prep.resolve() not in path.parents or not path.is_file():
        return None, f"{sf['path']} is missing or outside the prep folder"
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != sf["sha256"]:
        return None, (f"{sf['path']} no longer has the SHA-256 recorded in the manifest, so it is not "
                      "the file this run recorded")
    try:
        return json.loads(data.decode("utf-8")).get("result") or {}, None
    except ValueError as exc:
        return None, f"{sf['path']} could not be read ({exc})"


def _sample(values, times, stride):
    """Every `stride`-th point and the last one; nothing when the axis does not match the values."""
    if not values or len(values) != len(times):
        return []
    idx = list(range(0, len(values), stride))
    if idx[-1] != len(values) - 1:
        idx.append(len(values) - 1)
    return [(times[i] / 1000.0, values[i]) for i in idx]


def _md_figure(full: dict, stride: int, header: str = ""):
    """A builder for four panels drawn from the recorded series (reportlab.graphics, vector)."""
    palette = ["#1f4e8c", "#c0504d", "#2e7d32", "#7b1fa2", "#8a6d1d", "#00838f"]
    prot = full.get("protein_rmsd") or {}
    segments = prot.get("per_segment") or {(prot.get("segid") or "protein"): prot}
    panel_protein = []
    for i, (seg, ser) in enumerate(sorted(segments.items())):
        colour = palette[i % len(palette)]
        panel_protein.append((f"{seg} full", ser, colour, False))
        if isinstance(ser.get("observed_only"), dict):
            panel_protein.append((f"{seg} observed only (post hoc)", ser["observed_only"], colour, True))
    hb = dict((full.get("hbonds") or {}).get("count_per_frame") or {})
    # analyze_md records the per-frame bond count WITHOUT a time axis. It is counted over the same
    # analysed frames as the RMSDs (tested in test_md_analyze), so the ligand series' axis is used
    # when the lengths agree, and the caption says so; otherwise the panel is not drawn.
    lig_times = (full.get("ligand_rmsd") or {}).get("times_ps") or []
    if not hb.get("times_ps"):
        hb["times_ps"] = lig_times if len(hb.get("values") or []) == len(lig_times) else []
    panels = [("Protein backbone RMSD per segment (Å)", panel_protein),
              ("Ligand heavy-atom RMSD, fitted on the pocket (Å)",
               [("ligand", full.get("ligand_rmsd") or {}, palette[0], False)]),
              ("Ligand–pocket centre-of-mass distance (Å)",
               [("COM distance", full.get("com_distance") or {}, palette[0], False)]),
              ("Ligand–protein H-bonds per frame (count)", [("H-bonds", hb, palette[0], False)])]

    def build(width):
        from reportlab.graphics.charts.lineplots import LinePlot
        from reportlab.graphics.shapes import Drawing, Line, String
        from reportlab.lib import colors

        font = report._font()[0]
        cols, ph, top = 2, 150, 34
        # a banner INSIDE the drawing: with one figure per stage, two groups of four panels sat
        # side by side in the document with nothing on them saying which run they came from
        # (user request 2026-11-18). The caption said it; the picture did not.
        head = 16 if header else 0
        pw = width / cols
        d = Drawing(width, 2 * (ph + top) + head)
        if header:
            d.add(String(0, 2 * (ph + top) + 4, header, fontName=font, fontSize=8.5))
            d.add(Line(0, 2 * (ph + top), width, 2 * (ph + top), strokeColor=colors.HexColor("#555555"),
                       strokeWidth=0.6))
        for n, (title, lines) in enumerate(panels):
            x0 = (n % cols) * pw
            y0 = (1 - n // cols) * (ph + top)
            d.add(String(x0 + 30, y0 + ph + top - 10, title, fontName=font, fontSize=7.5))
            for k, (label, _s, colour, dashed) in enumerate(lines):
                lx = x0 + 30 + (k % 2) * (pw / 2 - 15)
                ly = y0 + ph + top - 20 - (k // 2) * 8
                d.add(Line(lx, ly + 2, lx + 14, ly + 2, strokeColor=colors.HexColor(colour),
                           strokeWidth=1, strokeDashArray=[3, 2] if dashed else None))
                d.add(String(lx + 17, ly, label, fontName=font, fontSize=6))
            data = [_sample(s.get("values") or [], s.get("times_ps") or [], stride)
                    for _l, s, _c, _d in lines]
            if not any(data):
                recorded = any(s.get("values") for _l, s, _c, _d in lines)
                d.add(String(x0 + 30, y0 + ph / 2,
                             "values recorded, but no time axis of the same length: not drawn"
                             if recorded else "no per-frame values recorded",
                             fontName=font, fontSize=7))
                continue
            lp = LinePlot()
            lp.x, lp.y, lp.width, lp.height = x0 + 30, y0 + 22, pw - 45, ph - 30
            lp.data = [series or [(0.0, 0.0)] for series in data]
            lp.joinedLines = 1
            for k, (_l, _s, colour, dashed) in enumerate(lines):
                lp.lines[k].strokeColor = colors.HexColor(colour)
                lp.lines[k].strokeWidth = 0.8
                if dashed:
                    lp.lines[k].strokeDashArray = [3, 2]
            for axis in (lp.xValueAxis, lp.yValueAxis):
                axis.labels.fontName, axis.labels.fontSize = font, 6
            lp.xValueAxis.valueMin = 0
            lp.yValueAxis.valueMin = 0
            d.add(lp)
            d.add(String(x0 + 30 + (pw - 45) / 2, y0 + 4, "time from the first analysed frame (ns)",
                         fontName=font, fontSize=6, textAnchor="middle"))
        return d
    return build


STAGE_TITLE = {
    "production": "The run the results come from (production)",
    "equilibration": "The run that settled the system first (equilibration)",
}
STAGE_WHY = {
    "production": ("This is the free simulation. Nothing is held in place, so what moves here is the "
                   "system's own behaviour and these are the figures to read."),
    "equilibration": ("This run exists to let the box settle, and parts of the system are held in "
                      "place by springs while it does. A number measured on a held atom describes "
                      "the spring, not the molecule, so the figures here are kept for completeness "
                      "and are not evidence about the pose."),
}


def _md_stage(doc, prep: Path, label: str, rec: dict, audit: list, subject: str = ""):
    res = rec.get("result") or {}
    traj = res.get("trajectory") or {}
    inputs = ", ".join(traj.get("stage_inputs") or []) or "no NAMD input of this stage was found"
    doc.h(4, STAGE_TITLE.get(label, label.capitalize()))
    doc.p(STAGE_WHY.get(label, ""), small=True)
    doc.p(_control_line(res))
    doc.p(f"{_v(traj.get('frames_note'))}, one every {_v(traj.get('dt_ps'))} ps, covering "
          f"{_v(traj.get('span_ns'))} ns between the first and the last of them. "
          f"Restraints: {_v((traj.get('restraints') or {}).get('note'))}", small=True)

    top, lig, pocket = rec.get("topology") or {}, res.get("ligand") or {}, res.get("pocket") or {}
    restr = traj.get("restraints") or {}
    held = restr.get("held") or {}
    com = res.get("com_distance") or {}
    repro = [f"Trajectory {t.get('file')}: {_v(t.get('frames'))} frames, {_v(t.get('size_bytes'))} bytes, "
             f"SHA-256 {_v(t.get('sha256'))}" for t in rec.get("trajectories") or []]
    repro += [
        f"Recorded {_v(rec.get('generated_at'))} by {_v(rec.get('generated_by'))}; the stage is called "
        f"'{label}' because {_v(rec.get('stage_label_source'))}; NAMD input {inputs}; stage "
        f"{_v(rec.get('stage'))}",
        f"Topology {_v(top.get('file'))}, SHA-256 {_v(top.get('sha256'))}"
        + (f"; equal to the PSF of the registered package ({(top.get('reference') or {}).get('source')})"
           if top.get("matches_registered_package") else ""),
        f"Minimization: {_v(traj.get('minimization_note'))}",
        f"Held atoms among those measured: ligand heavy atoms: "
        f"{_v((held.get('ligand_heavy') or {}).get('note'))}; protein backbone: "
        f"{_v((held.get('protein_backbone') or {}).get('note'))}",
        f"Ligand selection: {_v(lig.get('selection'))}; {_v(lig.get('dropped_note'))}",
        f"Pocket selection: {_v(pocket.get('selection'))} ({_v(pocket.get('criterion'))}): "
        + (", ".join(pocket.get("residues") or []) or NO_RECORD),
        f"Built residues inside that pocket: {_v(pocket.get('modelled_residues_in_pocket_note'))}",
        f"Periodic images: {_v(com.get('pbc_corrected_frames'))} of {_v(com.get('n'))} frames re-imaged; "
        f"non-rectangular boxes {_v(com.get('non_orthorhombic_frames'))}; frames without a box "
        f"{_v(com.get('frames_without_box'))}; frames with the pocket split across the box "
        f"{_v(com.get('pocket_split_frames'))}",
        "Arguments: " + (", ".join(f"{k} {v}" for k, v in (rec.get("arguments") or {}).items()
                                   if v not in (None, "")) or NO_RECORD),
    ]

    doc.h(4, "What the trajectory shows")
    prot = res.get("protein_rmsd") or {}
    rows, split = [], False
    for seg, ser in sorted((prot.get("per_segment") or {}).items()):
        rows.append(_metric_row(f"Protein backbone RMSD, {seg}, full segment (pre-registered control)",
                                ser, ser.get("selection") or f"{_v(ser.get('atoms'))} backbone atoms"))
        for part, name in (("observed_only", "observed residues only"),
                           ("modelled_only", "residues built by CHARMM-GUI only")):
            if isinstance(ser.get(part), dict):
                split = True
                rows.append(_metric_row(f"Protein backbone RMSD, {seg}, {name} (post hoc, NOT "
                                        "pre-registered)", ser[part]))
    if not rows:
        rows.append(_metric_row("Protein backbone RMSD (pre-registered control)", prot))
    rows.append(_metric_row("Ligand heavy-atom RMSD, fitted on the pocket backbone (never on the ligand)",
                            res.get("ligand_rmsd") or {},
                            f"{_v((res.get('ligand_rmsd') or {}).get('selection'))}; superposed on "
                            f"{_v((res.get('ligand_rmsd') or {}).get('superposed_on'))}"))
    rows.append(_metric_row("Ligand–pocket centre-of-mass distance", com))
    hb = res.get("hbonds") or {}
    rows.append(_metric_row("Ligand–protein hydrogen bonds per frame (a count, not Å)",
                            hb.get("count_per_frame") or {}, "pairs meeting the criterion below"))
    doc.table(["Metric", "Selection", "Min", "Max", "Mean", "First-quarter mean", "Last-quarter mean",
               "Drift", "Last-quarter range"], rows, (0.2, 0.2) + (0.6 / 7,) * 7)
    doc.p("Values in Å except the hydrogen-bond count. Drift = last-quarter mean minus first-quarter "
          "mean. Each protein segment is superposed on itself; segments are never fitted together. "
          f"Reference: {_v(prot.get('reference'))}.", small=True)
    if split:
        doc.p("HOW TO READ THE SPLIT ROWS. 'Observed' residues were present in the structure that was "
              "uploaded to CHARMM-GUI, so their coordinates come from the crystal. 'Built' residues "
              "were absent and CHARMM-GUI invented a conformation for them: they can move a long way "
              "without the protein doing anything interesting. 'Post hoc' means the split was looked "
              "at only AFTER the pre-registered control had already failed, so it is a second look at "
              "a result, not a test planned in advance, and it cannot turn a failed control into a "
              "passed one.", small=True)
        doc.p("The 'observed only' and 'built only' rows are an analysis added on 2026-11-12 AFTER the "
              "pre-registered full-backbone control failed on this project's first production run. Its "
              "reason (residues CHARMM-GUI built are a guess, not a measurement) does not depend on the "
              "result, but the decision to look came after seeing it: it was NOT pre-registered, and the "
              "full-segment row stays the control. 'Observed' means present in the structure uploaded to "
              "CHARMM-GUI.", small=True)
    warnings = [f"{name}: {s.get('drift_warning')}" for name, s in
                (("protein", prot), ("ligand", res.get("ligand_rmsd") or {}), ("COM", com))
                if (s or {}).get("drift_warning")]
    if warnings:
        doc.p("Too few frames for a plateau statement: " + "; ".join(warnings), small=True)

    pairs = hb.get("pairs") or []
    doc.h(4, "Which contacts the ligand keeps (hydrogen bonds)")
    doc.p("HOW TO READ THIS TABLE. One row is one donor-acceptor pair that satisfied the geometric "
          "criterion in at least one stored frame. The last column is the share of the stored frames "
          "in which it did: 1.00 means the contact was present every time a frame was written, 0.10 "
          "means it appeared in one frame out of ten. It is NOT a fraction of time -- hydrogen bonds "
          "form and break far faster than frames are stored, so a bond that flickers between two "
          "frames is invisible here. A high share says a contact is persistent in this trajectory; it "
          "does not say the contact is strong, and with no apo simulation it cannot say the contact "
          "holds the ligand in place.", small=True)
    doc.p(f"{_v(hb.get('pairs_found'))} pair(s) met the criterion in at least one frame. "
          f"Criterion: {_v(hb.get('criterion'))}.", small=True)
    doc.table(["Donor", "Hydrogen", "Acceptor", "Frames present / analysed", "Share of analysed frames"],
              [[p.get("donor"), p.get("hydrogen"), p.get("acceptor"),
                f"{_v(p.get('frames_present'))} / {_v(p.get('frames_analysed'))}",
                _num(p.get("occupancy_sampled"), "{:.2f}") or NO_RECORD] for p in pairs],
              (0.26, 0.12, 0.26, 0.18, 0.18))

    rmsf = res.get("rmsf") or {}
    doc.h(4, "How much each pocket residue moved (RMSF)")
    doc.p("HOW TO READ THIS TABLE. RMSF is how far one residue wandered from its own average "
          "position across the trajectory, in Å. A small value means that residue stayed put; a "
          "large one means it moved a lot, which can be a flexible loop, a side chain with nothing "
          "to hold it, or a residue CHARMM-GUI built. Rows are not ranked and there is no threshold: "
          "the useful comparison is between residues of the same pocket in the same run. Because no "
          "simulation of this pocket WITHOUT the ligand was run, a residue that barely moves cannot "
          "be credited to the ligand holding it.", small=True)
    doc.p(f"{_v(rmsf.get('selection'))}. {_v(rmsf.get('criterion'))}.", small=True)
    doc.table(["Residue", "RMSF (Å)"],
              [[r.get("residue"), _num4(r.get("rmsf_angstrom"))] for r in rmsf.get("per_residue") or []],
              (0.3, 0.2))

    full, why = _series_file(prep, rec)
    if full is None:
        doc.p(f"Plots omitted: {why}. They are drawn only from the recorded series file and never "
              "recomputed here.", small=True)
    else:
        n = len((full.get("ligand_rmsd") or {}).get("values") or [])
        stride = max(1, -(-n // MAX_PLOT_POINTS))
        sampling = (f"all {n} analysed frames plotted" if stride == 1 else
                    f"every {stride}th of {n} analysed frames plotted, and the last one")
        held_note = "restrained: parts of the system were held by springs" if held else "nothing held"
        banner = (f"{label.upper()} RUN"
                  + (f"  ·  {subject}" if subject else "")
                  + f"  ·  {_v(traj.get('span_ns'))} ns  ·  {held_note}")
        doc.figure(_md_figure(full, stride, banner),
                   f"These four panels all come from the {label} run named in the banner above them, "
                   f"and from no other run. Drawn from {(rec.get('series_file') or {}).get('path')}, "
                   f"whose SHA-256 matches the manifest; nothing is recomputed. {sampling.capitalize()}. "
                   "In the first panel a SOLID line is a whole protein segment (the pre-registered "
                   "control) and a DASHED line of the same colour is that same segment counting only "
                   "the residues that came from the crystal, leaving out the ones CHARMM-GUI built; "
                   "that dashed line was added after the control had already failed, so it is a second "
                   "look and not a planned test. The hydrogen-bond panel counts bonds per stored frame, "
                   "a sample, not a fraction of time; its count is recorded without times, so it is "
                   "placed on the ligand series' time axis (same analysed frames, same length), and "
                   "left out if the lengths differ.")

    reading = res.get("reading") or {}
    doc.h(4, "What this can and cannot be read as")
    doc.p(f"{_v(reading.get('status'))}. {_v(reading.get('disclaimer'))}", small=True)
    doc.bullets(list(reading.get("offered") or []))
    doc.p("Not measured here: " + " ".join(res.get("not_measured") or [NO_RECORD]), small=True)
    repro += [
        "Thresholds used in the reading above: "
        + (", ".join(f"{k} {v}" for k, v in (reading.get("thresholds_used") or {}).items()) or NO_RECORD),
        f"Where those thresholds come from: {_v(reading.get('thresholds_source'))}",
    ]
    _repro(audit, f"MD {label} · {subject or NO_RECORD}", repro)


def _md(doc, m: dict, accession: str, prep: Path, audit: list) -> int:
    doc.h(2, "2. Molecular dynamics")
    rows = _md_rows(m)
    if not rows:
        doc.p("Not run. No molecular dynamics package, job or analysis is recorded in prep/manifest.json. "
              "The chain is build_md_package, CHARMM-GUI, build_md_job, the cluster, and then "
              f"analyze_md(folder, uniprot_id=\"{accession}\", structure, pocket), which records the "
              "analysis; build this record again after it.")
        return 0
    runs = _md_runs(m)
    # Which compound was simulated has to be the first thing said. With seven pockets carrying a
    # package and one carrying a trajectory, a reader otherwise cannot tell which numbers below
    # describe a simulation and which describe a folder that was written and never run (user
    # request 2026-11-18).
    simulated = {(s.get("pdb_id"), k) for s, k, _l, _r in runs}
    done = [f"{s.get('pdb_id')} pocket {k} -- {_md_subject(m, k)}" for s, k, *_ in rows
            if (s.get("pdb_id"), k) in simulated]
    waiting = [(_compound_name(m, k).split(" (")[0] + f" [{s.get('pdb_id')} {k}]")
               if _compound_name(m, k) else f"{s.get('pdb_id')} {k}"
               for s, k, *_ in rows if (s.get("pdb_id"), k) not in simulated]
    lead = (f"Of the {len(rows)} pocket{'' if len(rows) == 1 else 's'} that carry a molecular dynamics "
            f"record, {len(simulated)} ha{'s' if len(simulated) == 1 else 've'} a trajectory that was "
            f"actually run and measured ({len(runs)} analysed stage{'' if len(runs) == 1 else 's'} in "
            "all, since a simulation is analysed stage by stage).")
    if done:
        lead += " SIMULATED AND MEASURED: " + "; ".join(done) + "."
    if waiting:
        lead += (" Package written but NO trajectory analysed, so nothing below is measured for "
                 + ("it" if len(waiting) == 1 else "them") + ": " + "; ".join(waiting) + ".")
    doc.p(lead)
    doc.p(MD_VOCABULARY, small=True)
    for s, key, package, job, analysis in rows:
        subject = _compound_name(m, key)
        doc.h(3, f"{s.get('pdb_id')} · pocket {key}" + (f" · {subject.split(' (')[0]}" if subject else ""))
        doc.p(f"What is simulated here: {_md_subject(m, key)}.", small=True)
        where = (f"{s.get('pdb_id')} pocket {key}"
                 + (f" · {subject.split(' (')[0]}" if subject else ""))
        _md_config(doc, package, job, audit, where)
        recorded = [lab for lab in MD_STAGES if isinstance(analysis.get(lab), dict)]
        if not recorded:
            doc.p("No trajectory analysis is recorded for this pocket. Run analyze_md on the folder "
                  f"brought back from the cluster with uniprot_id=\"{accession}\" to record one.",
                  small=True)
        for label in recorded:
            _md_stage(doc, prep, label, analysis[label], audit, subject=where)
    return len(runs)


def _set_docking(m: dict, compound_id: str) -> list:
    """Every recorded run of `compound_id`, each with the structure and pocket it belongs to.

    The run is found by the compound_id the run itself recorded, never by the run id, whose slug
    is a naming convention and not evidence of which compound was docked. ALL of them are
    returned, and the caller says which one it prints and how many there were: returning the
    first match silently made a second run of the same compound invisible, and with it the fact
    that one of the two may be the biased `conformer="crystal"` run (reviewer H2, 2026-11-22).
    """
    out = []
    for s in m.get("structures") or []:
        for pocket_key, block in sorted((s.get("docking") or {}).items()):
            for run_id, run in (block.get("runs") or {}).items():
                cd = (run or {}).get("compound_docking") or {}
                if cd.get("compound_id") == compound_id:
                    # the conformer is recorded on the run's ligand block, which is where
                    # dock_ligand writes it; it is not a field of `vina`
                    out.append({"cd": cd, "pdb": s.get("pdb_id"), "pocket": pocket_key,
                                "run_id": run_id,
                                "conformer": ((run or {}).get("ligand") or {}).get("conformer"),
                                "cross_docking": (run or {}).get("cross_docking") or {}})
    return out


def _key_is_of(key: str, compound_id: str) -> bool:
    """Whether a manifest key of the form "<pocket> <CHEMBL id>" belongs to this compound."""
    return bool(key) and str(key).endswith(" " + compound_id)


def _set_md(m: dict, compound_id: str) -> list:
    """Every recorded PRODUCTION analysis of `compound_id`, with the structure and pocket.

    Matched on the `compound` block analyze_md writes since step 6b AND on the "<pocket> <CHEMBL
    id>" key, because a record written before 6b has no compound block: keyed only on the block,
    a trajectory that was in fact analysed was printed as "not run" while section 2 of the same
    document described it (reviewer H3, 2026-11-22). "not run" is a claim about what was run.
    """
    out = []
    for s in m.get("structures") or []:
        for key, analysis in sorted((s.get("md_analysis") or {}).items()):
            if not isinstance(analysis, dict):
                continue
            for stage_name, stage in analysis.items():
                if (stage or {}).get("stage_label") != "production":
                    continue
                comp = ((stage or {}).get("compound") or {})
                if comp.get("id") == compound_id or (not comp and _key_is_of(key, compound_id)):
                    out.append({"stage": stage, "key": key, "stage_name": stage_name,
                                "pdb": s.get("pdb_id"),
                                "pocket": (stage or {}).get("pocket")
                                or str(key).rsplit(" ", 1)[0],
                                "identified_by": ("the compound block of the record" if comp
                                                  else "the manifest key, which ends in this "
                                                       "ChEMBL id (record written before the "
                                                       "compound block existed)")})
    return out


def _set_admet(m: dict, compound_id: str) -> list:
    """Every recorded ADMET block of `compound_id`, with the structure and pocket of its key."""
    out = []
    for s in m.get("structures") or []:
        for key, block in sorted((s.get("admet") or {}).items()):
            ident = (block or {}).get("compound_identity") or {}
            if ident.get("compound_id") == compound_id or (not ident and _key_is_of(key, compound_id)):
                out.append({"block": block, "pdb": ident.get("pdb_id") or s.get("pdb_id"),
                            "pocket": ident.get("pocket") or str(key).rsplit(" ", 1)[0]})
    return out


def _compound_set(doc, m: dict, accession: str) -> int:
    """Section 4: one row per selected compound, with the holes declared and NO ranking.

    The researcher asked for the set side by side (2026-11-22). It is a coverage map, not a
    verdict: the rows are in the order select_ligands recorded them, and the document says in
    the open why the columns must not be read as an order of merit. Three reasons, all of them
    measured and none of them a matter of taste: a raw Vina score favours the larger compound
    (so ligand efficiency travels beside it and neither re-ranks the set); a box_limited run
    scores the fit of the box as much as the fit of the pocket; and an MD column built on one
    replica of one length has no error bar, so a difference between two compounds cannot be
    told from noise. Any compound whose MD or ADMET was not run says so in its own cell.
    """
    doc.h(2, "4. The compound set side by side")
    choice_block = m.get("ligand_choice") or {}
    compounds = choice_block.get("compounds") or []
    if not compounds:
        doc.p("Not run. No compound selection is recorded in prep/manifest.json, so there is no set "
              f"to lay side by side. Run select_ligands(\"{accession}\", compounds, reason) to record "
              "one, and build this record again.")
        return 0
    rows, with_md, with_dock, with_admet = [], 0, 0, 0
    heavy_seen, levels_seen, limited, no_negative_control = [], [], [], False
    for c in compounds:
        cid = c.get("id")
        docks, mds, admets = _set_docking(m, cid), _set_md(m, cid), _set_admet(m, cid)
        if c.get("heavy_atoms") is not None:
            heavy_seen.append(c["heavy_atoms"])
        if c.get("level") and c["level"] not in levels_seen:
            levels_seen.append(c["level"])
        if docks:
            with_dock += 1
            first = docks[0]
            cd = first["cd"]
            where = f"{first['pdb']} {first['pocket']}" + _of_n(docks, "run")
            score, le = (_num(cd.get("top_pose_score_kcal_mol")),
                         _num(cd.get("ligand_efficiency_kcal_mol_per_heavy_atom")))
            dock_cell = f"{score} ({where})" if score else f"{NO_RECORD} ({where})"
            if (first.get("conformer") or "") == "crystal":
                dock_cell += (" — started from the crystal conformation: biased in favour of the "
                              "crystal pose, not comparable with a regenerated run")
            fit = cd.get("box_fit") or {}
            if cd.get("box_limited"):
                limited.append(_v(c.get("name")) or cid)
                dock_cell += (f" — box_limited: {_num(fit.get('min_clearance_angstrom'))} A clearance "
                              f"on the {_num(fit.get('limiting_box_axis_angstrom'))} A axis against a "
                              f"{_num(fit.get('limiting_ligand_extent_angstrom'))} A extent, over a "
                              f"threshold of {_num(fit.get('threshold_angstrom'))} A. The score measures "
                              "the fit of the box as much as the fit of the pocket and is NOT "
                              "comparable with the rest")
            if "no inactive compounds" in str(cd.get("scope") or ""):
                no_negative_control = True
            le_cell = f"{le} over {_v(cd.get('heavy_atoms'))} heavy atoms" if le else NO_RECORD
            xd = first.get("cross_docking") or {}
            if xd.get("top_pose_rmsd_angstrom") is not None:
                bar = (xd.get("reference_pose_uncertainty") or {}).get("shift_angstrom")
                # `within_criterion_with_bar` is a BLOCK, not a flag: its `decides` field is what
                # says whether the criterion settles anything once the superposition bar is
                # applied. Read as a boolean it is always truthy, and the cell then printed the
                # figure beside the criterion as if it had passed.
                with_bar = xd.get("within_criterion_with_bar") or {}
                span = with_bar.get("range_angstrom") or []
                rmsd_cell = (f"{_num(xd['top_pose_rmsd_angstrom'])} ± {_num(bar)} vs a criterion of "
                             f"{_num(xd.get('criterion_angstrom'))} A"
                             + ("" if with_bar.get("decides") is not False else
                                (f" — DECIDES NOTHING: with the bar applied the value spans "
                                 f"{_num(span[0])}-{_num(span[1])} A, both sides of the criterion"
                                 if len(span) == 2 else
                                 " — DECIDES NOTHING once the superposition bar is applied")))
            else:
                rmsd_cell = ("not applicable: no crystallographic pose of this compound in another "
                             "structure of this target (level "
                             f"{_v(c.get('level'))}), so there is nothing to measure against")
        else:
            dock_cell = le_cell = rmsd_cell = "not run (dock_ligand)"
        if mds:
            with_md += 1
            first_md = mds[0]
            result = (first_md["stage"].get("result") or {})
            lig = result.get("ligand_rmsd") or {}
            control = ((result.get("reading") or {}).get("equilibration_control") or {})
            where_md = (f"{_v(first_md.get('pdb'))} {_v(first_md.get('pocket'))}"
                        + _of_n(mds, "production run"))
            md_cell = ((f"{_num(lig.get('mean'))} mean over {_v(lig.get('n'))} frames ({where_md})")
                       if lig.get("mean") is not None else f"{NO_RECORD} ({where_md})")
            control_cell = (f"{_v(control.get('status'))} (drift {_num(control.get('drift_angstrom'))} A "
                            f"vs {_num(control.get('threshold_angstrom'))} A)"
                            if control else NO_RECORD)
        else:
            md_cell = control_cell = "not run (build_md_package, then analyze_md)"
        if admets:
            with_admet += 1
            rec = (admets[0]["block"].get("chembl_record") or {})
            status = str(rec.get("status") or "")
            admet_cell = (("ChEMBL entry retrieved; its ADME records and calculated properties are "
                           f"in section 3 ({_v(admets[0].get('pdb'))} {_v(admets[0].get('pocket'))}"
                           + _of_n(admets, "result") + ")")
                          if status.startswith("RECORD") else
                          ("No entry for this compound in ChEMBL: nothing is recorded there. That is "
                           "an absence of DATA, not an absence of liabilities"
                           if status else NO_RECORD))
        else:
            admet_cell = "not run (analyze_admet)"
        rows.append([f"{_v(c.get('name')) if c.get('name') else NO_RECORD}\n{cid}\n{_v(c.get('level'))}",
                     dock_cell, le_cell, rmsd_cell, md_cell, control_cell, admet_cell])
    doc.p(f"{len(compounds)} compound(s) in the selection select_ligands recorded. "
          f"Docking recorded for {with_dock} of {len(compounds)}, molecular dynamics for "
          f"{with_md} of {len(compounds)}, ADMET for {with_admet} of {len(compounds)}. The rows are "
          "in the order the selection recorded them.")
    # the reason in full, never clipped: the clipped version cut the sentence that says which of
    # these compounds is the anchor rather than a candidate, and the row count then read as a
    # count of candidates (reviewer H7, 2026-11-22)
    if choice_block.get("reason"):
        doc.p("Why these compounds, as the researcher recorded it: "
              + str(choice_block["reason"]), small=True)
    doc.p("THIS TABLE IS NOT A RANKING, and the columns must not be read as an order of merit. "
          + (f"The set spans levels {', '.join(levels_seen)}, which are NOT comparable with each "
             "other and are not ranked together: a level T1 compound has a crystallographic pose "
             "to measure against and a T3 compound has none, so the same column means a different "
             "thing in their rows. " if len(levels_seen) > 1 else "")
          + "A Vina score is approximately additive over atoms, so a raw score favours the larger "
          "compound for being larger"
          + (f"; across this set the heavy-atom count runs from {min(heavy_seen)} to "
             f"{max(heavy_seen)}" if heavy_seen else "")
          + ". Ligand efficiency is printed beside it and neither value re-ranks the set, because "
          "no threshold for an interpretable difference has been measured in this project. "
          + ("NO NEGATIVE CONTROL: the set holds no compound known to be inactive against this "
             "target, so nothing here says whether the score separates actives from inactives, and "
             "no ordering of these scores is supported (pre-registered consequence, ESTADO 7.4). "
             if no_negative_control else "")
          + (f"The run of {', '.join(limited)} is marked box_limited, with its clearance in the cell. "
             if limited else "")
          + "The dynamics column comes from ONE replica per compound, so it carries no error bar "
          "and a difference between two compounds cannot be told from noise; where the equilibration "
          "control says NOT MET, the protein was still drifting over the trajectory and nothing "
          "about the ligand is settled. A compound whose dynamics or ADMET was not run says so in "
          "its own cell: an empty cell is a gap in what was run, never a result about the compound.",
          small=True)
    doc.table(["Compound / ChEMBL id / level", "Top pose score (kcal/mol, calculated)",
               "Ligand efficiency (calculated)", "Cross-docking RMSD (A, calculated)",
               "Ligand RMSD, production (A, calculated)", "Equilibration control (calculated)",
               "ADMET (retrieved)"],
              rows, [0.15, 0.16, 0.12, 0.15, 0.15, 0.14, 0.13])
    doc.p("THE THREE PROVENANCES IN THIS TABLE DO NOT MEAN THE SAME THING. The scores, ligand "
          "efficiencies and RMSD values are CALCULATED by AutoDock Vina and RDKit; the dynamics "
          "columns are CALCULATED by NAMD and MDAnalysis; the ADMET column is RETRIEVED from "
          "ChEMBL. A calculated value is not a measurement and none of them is a measured "
          "affinity: measured potency against the target is in the research report, "
          f"report_{accession}.pdf. Nothing in this table comes from a trained model. Every figure "
          "is quoted from the same manifest blocks the sections above print in full, where each one "
          "travels with its own selection and criterion, and the structure and pocket each figure "
          "belongs to is named in its cell.", small=True)
    return len(compounds)


def _of_n(found: list, what: str) -> str:
    """", 1 of N <what>s recorded" when more than one was found. The table prints the first and
    must say so: silently printing one of several hides the others, and one of them may be the
    biased crystal-conformer run or another pocket entirely (reviewer H2)."""
    return "" if len(found) <= 1 else f", 1 of {len(found)} {what}s recorded — see section 1 to 3"


def _methods(doc, m: dict):
    doc.h(2, "5. Methods and limitations")
    doc.p("This record is generated deterministically from prep/manifest.json: every value and every identifier is "
          f"quoted from the manifest, fields without data are shown as '{NO_RECORD}', and nothing is recomputed or "
          "inferred here. No text of a Vina log is copied into this document.")
    doc.bullets([
        VOCABULARY,
        "Re-docking the crystal ligand of a pocket into its own receptor tests the set-up, not the method: the "
        "scope line of each run states what it covers.",
        "The receptor is rigid and was not completed: missing loops stay missing, side chains are not rebuilt, and "
        "residues that match no Meeko template are listed per run (deleted only when outside the box).",
        "Wall-clock and CPU times are measurements of the machine that ran the docking, not properties of the result.",
    ])
    if m.get("decisions"):
        doc.p("Rules that produced the receptor, the pocket and the box used above "
              "(prep/manifest.json, decisions):")
        doc.bullets(list(m["decisions"]))
    doc.p("Preparation steps still pending in the manifest: "
          + ("; ".join(m.get("pending") or []) or "none listed") + ". "
          + ("Molecular dynamics (section 2): no analysed run is recorded." if not _md_runs(m) else
             "Molecular dynamics (section 2): the simulations ran outside this machine; their "
             "trajectories were measured locally by analyze_md and every figure, file hash and "
             "plotted series is quoted from what it recorded."))


def _audit_section(doc, audit: list):
    """Section 6: the whole audit trail, after everything a reader interprets.

    Collected by _repro while sections 1-3 are written, so the order here follows the order of
    the runs above and each entry keeps the name of the run it belongs to. Nothing is summarised
    and nothing is dropped: this is the part that lets someone else repeat the work.
    """
    doc.h(2, "6. Reproducibility and audit trail")
    if not audit:
        doc.p("Nothing to record: no run above wrote files, hashes or tool options into the manifest.")
        return
    doc.p(f"{len(audit)} block{'' if len(audit) == 1 else 's'}, one per run of sections 1 to 3, in the "
          "order those runs appear. This is the material needed to REPEAT or AUDIT the work: file "
          "paths and their SHA-256, the options each tool was given, the checks that ran before it, "
          "the superposition behind a cross-docking RMSD, the raw database fields and the thresholds "
          "a reading used. Not one line of it is needed to interpret a result, which is why it is "
          "here and not above; not one line of it is optional, which is why it is here at all. "
          "A hash is what proves "
          "the file on disk is the file the run wrote -- analyze_md refuses to measure a trajectory "
          "whose topology does not match it.")
    for entry in audit:
        doc.h(3, entry["title"])
        doc.bullets(entry["items"])


def build_document(accession: str, m: dict, generated: str, prep: Path = None) -> tuple:
    """(Doc, docking runs recorded, ADMET results recorded, compounds in the selection).

    `prep` is the folder the manifest's relative paths start from (the series files of section 2);
    it defaults to the dossier's prep folder."""
    doc = Doc()
    doc.h(1, f"In-silico experiment record: {accession}")
    seq = m.get("sequence") or {}
    doc.p(f"Generated {generated} from prep/manifest.json (created {_v(m.get('created_at'))} by "
          f"{_v(m.get('generated_by'))}). Variant filter: {m.get('mutation_filter') or 'none'}. "
          f"UniProt sequence used for the numbering: {seq.get('header') or NO_RECORD}", small=True)
    doc.p("Results of in-silico experiments (calculated docking scores and RMSD values; ADME records, retrieved "
          "measurements and calculated properties of the ligand), each with the provenance of the structure it came "
          "from, quoted from the manifest. No value in this document is a prediction: nothing here comes from a "
          "trained model. Measured potency against the target (ChEMBL activities) is not part of this document; it "
          f"is in the research report, report_{accession}.pdf. Interpretation is not part of either.", small=True)
    doc.p("HOW THIS DOCUMENT IS ORDERED. Sections 1 to 3 hold what a researcher reads to interpret "
          "the experiments: what was run, on what, and what came out. Section 4 lays the selected "
          "compounds side by side, with the holes declared and no ranking. Section 5 states the methods "
          "and the limits. Everything needed to REPEAT or AUDIT a run -- file paths and hashes, tool "
          "options, pre-flight checks, raw database fields, the thresholds behind a reading -- is in "
          "section 6 at the end, grouped by run. Nothing was deleted to get it out of the way.",
          small=True)
    audit = []
    n = _docking(doc, m, accession, audit)
    _md(doc, m, accession, prep or choice.manifest_path(accession).parent, audit)
    n_admet = _admet(doc, m, accession, audit)
    n_set = _compound_set(doc, m, accession)
    _methods(doc, m)
    _audit_section(doc, audit)
    return doc, n, n_admet, n_set


def _section_inventory(doc) -> list:
    """What the document's sections are and which are declared holes, READ BACK from the blocks.

    F4, from the Desktop test of 2026-11-12: asked what section 3 said, the answer was "the tool
    output does not number the sections", which was true -- so the analyst could not know what
    the document holds without opening it. Derived from the document itself and never from a
    hand-kept list: a parallel list is a claim about a document that drifts from it, which is
    the same mistake as the file-name stage rule of analyze_md (ESTADO 4, tool 18).
    """
    out = []
    for b in doc.blocks:
        if b[0] == "h" and b[1] == 2:
            number, _, title = str(b[2]).partition(". ")
            out.append({"number": number.strip(), "title": title.strip() or str(b[2]),
                        "heading": str(b[2]), "status": "has content"})
        elif out and b[0] == "p" and out[-1]["status"] == "has content":
            text = " ".join(str(p) for p in b[1] if p)
            if text.startswith("Not run."):
                out[-1]["status"] = "DECLARED HOLE (not run)"
                out[-1]["says"] = text if len(text) <= 300 else text[:299] + "…"
    return out


def build_experiment_report(uniprot_id: str) -> dict:
    """prep/manifest.json -> experiments_<accession>.md and .pdf in dossiers/<accession>/."""
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    m = report.load_manifest(accession)
    if m is None:
        return {"ok": False, "uniprot_id": accession,
                "error": f"No preparation manifest for {accession}. Run prepare_structures first: the experiment "
                         "record is built from dossiers/<accession>/prep/manifest.json."}
    if m.get("unreadable"):
        return {"ok": False, "uniprot_id": accession,
                "error": f"prep/manifest.json could not be read ({m['unreadable']})."}
    generated = utc_now()
    doc, n_runs, n_admet, n_set = build_document(accession, m, generated)
    folder = dossier.DOSSIER_DIR / accession
    md_path = folder / f"experiments_{accession}.md"
    pdf_path = folder / f"experiments_{accession}.pdf"
    try:
        # PDF first: when it is open in a viewer the write fails and the Markdown is left as it was,
        # so the two files never disagree
        pages = report.render_pdf(doc, pdf_path, title=f"In-silico experiment record {accession}",
                                  footer=f"{accession} · experiment record · generated {generated} "
                                         "from prep/manifest.json only")
        md_path.write_text(report.render_markdown(doc), encoding="utf-8")
    except OSError as exc:   # typically the PDF is open in a viewer on Windows
        log.exception("Could not write the experiment record")
        return {"ok": False, "uniprot_id": accession,
                "error": f"Could not write the experiment record ({exc}). The PDF is written first, so if it could "
                         "not be written neither file was replaced. If the PDF is open, close it and try again."}
    runs = [f"{s.get('pdb_id')} {key} {run_id} ({run.get('computed_at')})"
            for s, key, _b, run_id, run in _run_rows(m)]
    font_note = report._font()[2]
    return {
        "ok": True,
        "uniprot_id": accession,
        "pdf_path": str(pdf_path),
        "markdown_path": str(md_path),
        "pages": pages,
        "generated": generated,
        "manifest_created_at": m.get("created_at"),
        "docking_runs": runs,
        "admet_compounds": [f"{s.get('pdb_id')} {key} "
                            f"{((b.get('compound_identity') or {}).get('inchikey')) or 'no InChIKey'}"
                            for s, key, b in _admet_rows(m)],
        # the inventory is the authoritative answer to "what does section N say": it is read
        # back from the document that was just written, so it cannot disagree with it
        "sections": _section_inventory(doc),
        "md_runs": [f"{s.get('pdb_id')} {key} {label} ({rec.get('generated_at')}): "
                    + str(((rec.get("result") or {}).get("reading") or {})
                          .get("equilibration_control", {}).get("status", "control not recorded"))
                    for s, key, label, rec in _md_runs(m)],
        "sections_not_run": (["molecular dynamics (no analyze_md result in the manifest)"]
                             if not _md_runs(m) else [])
                            + (["ADMET (no analyze_admet result in the manifest)"] if not n_admet else [])
                            + (["docking (no run in the manifest)"] if not n_runs else [])
                            + (["compound set (no select_ligands selection in the manifest)"]
                               if not n_set else []),
        **({"font_warning": font_note} if font_note else {}),
    }
