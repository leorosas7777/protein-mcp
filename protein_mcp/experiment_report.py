"""build_experiment_report: prep/manifest.json -> Markdown + PDF record of the in-silico experiments.

A document of its own, separate from build_report (decided 2026-11-05), with the three
sections present from the start: docking, molecular dynamics, ADMET. Sections with no run
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
        f"Crystal ligand pose {lig.get('sdf') or NO_RECORD}, SHA-256 {lig.get('sdf_sha256') or NO_RECORD}",
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
    items = [
        f"{lig.get('comp_id') or NO_RECORD} (instance {_v(lig.get('label_asym_id'))}, auth chain "
        f"{_v(lig.get('auth_asym_id'))}); wwPDB synonyms: {syn}",
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


def _receptor_block(doc, run: dict):
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
    doc.h(4, "Receptor preparation")
    doc.bullets(items)
    guard = run.get("guard") or {}
    if guard:
        doc.h(4, "Checks run before Vina was started")
        doc.bullets([f"Ligand atoms inside the box: {guard.get('ligand_atoms_in_box') or NO_RECORD}",
                     f"Ligand centroid inside the box: {_v(guard.get('ligand_centroid_in_box'))}",
                     f"Receptor atoms within 8 Å of the box centre: "
                     f"{_v(guard.get('receptor_atoms_within_8A_of_box_centre'))}; closest at "
                     f"{_num(guard.get('closest_receptor_atom_to_box_centre_angstrom')) or NO_RECORD} Å"]
                    + ([guard["note"]] if guard.get("note") else []))


def _vina_block(doc, run: dict):
    v = run.get("vina") or {}
    doc.h(4, "AutoDock Vina run")
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
    doc.bullets(items)


def _poses_block(doc, run: dict):
    v, rd = run.get("vina") or {}, run.get("redocking") or {}
    doc.h(4, "Poses (calculated)")
    n, asked = len(run.get("poses") or []), v.get("num_modes")
    count = (f"{n} poses returned (num_modes {asked})" if n == asked else
             f"{n} poses returned; num_modes {_v(asked)}, energy_range {_v(v.get('energy_range'))} kcal/mol"
             + ("; the manifest does not record why fewer were returned" if isinstance(asked, int) and n < asked
                else ""))
    doc.p(f"{count}. Both columns are calculated values of this run; the RMSD is against the crystal pose of this "
          f"same pocket (method: {rd.get('rmsd_method') or NO_RECORD}).", small=True)
    doc.table(["Mode", "Vina score (kcal/mol, calculated)", "RMSD to the crystal pose (Å, calculated)"],
              [[p.get("mode"), _num(p.get("vina_score_kcal_mol"), "{:.2f}") or NO_RECORD,
                _num(p.get("rmsd_to_crystal_pose_angstrom"), "{:.4f}") or NO_RECORD]
               for p in run.get("poses") or []], (0.08, 0.3, 0.32))
    nc = rd.get("rmsd_negative_control") or {}
    items = [f"Top-pose RMSD: {_num(rd.get('top_pose_rmsd_angstrom'), '{:.4f}') or NO_RECORD} Å; criterion "
             f"<= {_num(rd.get('criterion_angstrom')) or NO_RECORD} Å; within the criterion: "
             f"{_v(rd.get('within_criterion'))}",
             f"Pre-registered control run: {_v(rd.get('is_preregistered_control'))}",
             f"RMSD method: {rd.get('rmsd_method') or NO_RECORD}",
             "Negative control of this run, measured with the same function as the poses: the crystal pose "
             f"against itself gave {_num(nc.get('crystal_vs_itself_angstrom')) or NO_RECORD} Å and the crystal "
             f"pose moved 5 Å gave {_num(nc.get('crystal_moved_5A_angstrom')) or NO_RECORD} Å; "
             f"passed: {_v(nc.get('passed'))}",
             f"Scope: {rd.get('scope') or NO_RECORD}"]
    if rd.get("note"):
        items.append(rd["note"])
    doc.bullets(items)


def _files_block(doc, run: dict):
    rows = [[name, f.get("path"), f.get("sha256")] for name, f in sorted((run.get("files") or {}).items())]
    doc.h(4, "Files written (paths relative to the prep folder)")
    doc.table(["Name", "Path", "SHA-256"], rows, (0.22, 0.36, 0.42))


def _docking(doc, m: dict, accession: str):
    doc.h(2, "1. Molecular docking (AutoDock Vina)")
    runs = _run_rows(m)
    structures = m.get("structures") or []
    pockets = sum(1 for s in structures for pk in (s.get("prepared") or {}).get("pockets") or [] if not pk.get("error"))
    if not runs:
        doc.p("Not run. No docking run is recorded in prep/manifest.json. Run dock_ligand(\""
              f"{accession}\", structure, pocket) and then dock_status to record one, and build this record again.")
        return 0
    docked_pockets = {(s.get("pdb_id"), key) for s, key, _b, _r, _run in runs}
    doc.p(f"{len(runs)} docking run{'' if len(runs) == 1 else 's'} recorded, on {len(docked_pockets)} of the "
          f"{pockets} pockets computed in the "
          f"manifest ({len({s.get('pdb_id') for s, *_ in runs})} of {len(structures)} structures). Re-docking of the "
          "crystal ligand of a pocket into its own receptor: one ligand per run. " + VOCABULARY, small=True)
    doc.table(["PDB", "Pocket", "Run", "Recorded (UTC)", "Written by", "Started from"],
              [[s.get("pdb_id"), key, run_id, _v(run.get("computed_at")), run.get("generated_by") or NO_RECORD,
                run.get("selected_by") or NO_RECORD] for s, key, _b, run_id, run in runs],
              (0.08, 0.1, 0.22, 0.2, 0.2, 0.2))
    for s, key, block, run_id, run in runs:
        doc.h(3, f"{s.get('pdb_id')} · pocket {key} · run {run_id}")
        if run.get("replaced_run_computed_at"):
            doc.p(f"Replaced an earlier run of the same pocket and conformer, recorded "
                  f"{run['replaced_run_computed_at']}; only this one is kept in the manifest.", small=True)
        doc.h(4, "Provenance")
        _provenance(doc, s, block, run, m, key)
        _ligand_block(doc, run)
        _receptor_block(doc, run)
        _vina_block(doc, run)
        _poses_block(doc, run)
        _files_block(doc, run)
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


def _admet(doc, m: dict, accession: str) -> int:
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
        doc.h(3, f"{s.get('pdb_id')} · pocket {key} · {ident.get('comp_id') or NO_RECORD}")
        doc.p(f"Recorded {_v(block.get('generated_at'))} by {_v(block.get('generated_by'))}. "
              f"{block.get('selection') or ''}", small=True)

        doc.h(4, "Compound identity")
        doc.table(["Field", "Value"], [
            ["InChIKey (computed here with RDKit)", ident.get("inchikey") or NO_RECORD],
            ["Read from", ident.get("source") or NO_RECORD],
            ["File", ident.get("source_file") or NO_RECORD],
            ["SHA-256 of that file", ident.get("source_sha256") or NO_RECORD],
            ["Docking run", ident.get("run_id") or "none (crystallographic ligand)"],
            ["Agreement between the files", ident.get("identity_check") or NO_RECORD],
            ["Heavy atoms (SDF / manifest pocket)",
             f"{(block.get('calculated_properties') or {}).get('rdkit', {}).get('heavy_atoms', NO_RECORD)}"
             f" / {_v(ident.get('heavy_atoms_in_manifest_pocket'))}"],
            ["ChEMBL molecules with this exact InChIKey", _v(look.get("matches"))],
        ], (0.38, 0.62))

        rec = block.get("chembl_record") or {}
        doc.h(4, "ChEMBL record (facts of the entry, not measurements)")
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
            doc.p(Lines([Link(f"ChEMBL {rec['chembl_id']}", rec.get("link")),
                         f" · name in ChEMBL: {rec.get('pref_name') or 'unnamed (no pref_name)'}"]),
                  small=True)
            minus_one = rec.get("fields_minus_one") or []
            false_flags = rec.get("boolean_fields_false") or []
            doc.table(["Field", "Value"],
                      [[f, f"{_v(fields.get(f))} (-1, see below)" if f in minus_one
                        else f"{_v(fields.get(f))} (default, see below)" if f in false_flags
                        else _v(fields.get(f))] for f in
                       ("max_phase", "first_approval", "withdrawn_flag", "black_box_warning",
                        "oral", "parenteral", "topical", "molecule_type", "therapeutic_flag",
                        "usan_stem", "first_in_class")]
                      + [["ATC", ", ".join(rec.get("atc_classifications") or []) or NO_RECORD]],
                      (0.38, 0.62))
            if minus_one or false_flags:
                doc.p(rec.get("defaults_note") or "", small=True)
                if minus_one:
                    doc.p("Returned as -1: " + ", ".join(minus_one) + ".", small=True)
                if false_flags:
                    doc.p("Returned as False: " + ", ".join(false_flags) + ". On an entry with "
                          "nothing annotated every one of these is False, so none of them is "
                          "evidence that the compound lacks that property.", small=True)

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


def _site_value(site: dict, name: str):
    field = site.get(name) or {}
    return field.get("gres") or field.get("value") if isinstance(field, dict) else None


def _md_config(doc, package: dict, job: dict):
    """Subtask 5.1: the configuration as recorded by build_md_package and build_md_job, quoted."""
    rows = []
    if package:
        pose = package.get("pose") or {}
        rows.append(["Starting pose",
                     f"pose {_v(pose.get('index'))} (Vina mode {_v(pose.get('mode'))}) of docking run "
                     f"{_v(package.get('run_id'))}; Vina score {_v(pose.get('vina_score_kcal_mol'))} "
                     f"kcal/mol and RMSD to the crystal pose "
                     f"{_v(pose.get('rmsd_to_crystal_pose_angstrom'))} Å, both calculated",
                     "md_package"])
        gaps = package.get("to_model")
        rows.append(["Gaps left for CHARMM-GUI to build",
                     ", ".join(f"{g.get('chain')} {g.get('first')}-{g.get('last')} ({g.get('kind')})"
                               for g in gaps) if gaps else ("none" if gaps == [] else NO_RECORD),
                     f"md_package (gaps read from {_v(package.get('missing_residues_source'))})"])
    else:
        rows.append(["Package for CHARMM-GUI", "no md_package recorded for this pocket", "—"])
    if job:
        checks = (job.get("verified") or {}).get("checks") or {}
        closest = checks.get("closest_modelled_run_to_ligand") or {}
        built = checks.get("modelled_runs_present")
        rows.append(["Residues built in the system CHARMM-GUI returned",
                     (", ".join(built) if built else ("none" if built == [] else NO_RECORD))
                     + (f"; closest to the ligand: {closest.get('run')} at "
                        f"{_v(closest.get('min_distance_angstrom'))} Å" if closest else ""),
                     "md_job (archive checks)"])
        arch, psf = job.get("archive") or {}, job.get("psf") or {}
        rows.append(["CHARMM-GUI job", f"{_v(job.get('job_id'))} (folder {_v(job.get('job_folder'))}); "
                                       f"archive SHA-256 {_v(arch.get('sha256'))}", "md_job"])
        rows.append(["Topology (step3_input.psf) SHA-256",
                     psf.get("sha256") or (f"none: {psf['source']}" if psf.get("source") else
                                           "not recorded by the version of build_md_job that "
                                           "wrote this job"), "md_job"])
        sysm = checks.get("system") or {}
        ions = sysm.get("ions") or {}
        rows.append(["System",
                     f"box side {_v(sysm.get('box_side_angstrom'))} Å; {_v(sysm.get('waters'))} "
                     "waters; ions " + (", ".join(f"{k} {v}" for k, v in ions.items()) or NO_RECORD)
                     + f"; {_v(sysm.get('temperature_kelvin'))} K; time step "
                       f"{_v(sysm.get('timestep_fs'))} fs; cutoff {_v(sysm.get('cutoff_angstrom'))} Å",
                     "md_job (archive checks)"])
        inp, slurm = job.get("inputs") or {}, job.get("slurm") or {}
        val = lambda k: _v((inp.get(k) or {}).get("value") if isinstance(inp.get(k), dict)  # noqa: E731
                           else inp.get(k))
        rows.append(["Equilibration as written",
                     f"{_v(inp.get('equilibration_input'))}: minimize {val('equilibration_minimize')} "
                     f"steps, then run {val('equilibration_run')} steps", "md_job"])
        rows.append(["Production as written",
                     f"{_v(inp.get('production_input'))}: {_v(slurm.get('chunks'))} chunk(s) × "
                     f"{_v(inp.get('nanoseconds_per_chunk'))} ns = {_v(slurm.get('nanoseconds_total'))} "
                     f"ns; dcdfreq {val('dcdfreq')} steps", "md_job"])
        edited = job.get("inputs_edited")
        rows.append(["Inputs edited by build_md_job",
                     "none" if edited == [] else (", ".join(str(e.get("file") or e) for e in edited)
                                                  if edited else NO_RECORD), "md_job"])
        site = job.get("site") or {}
        rows.append(["Cluster job",
                     f"partition {_v(_site_value(site, 'partition'))}; GPUs "
                     f"{_v(_site_value(site, 'gpus'))}; time limit {_v(_site_value(site, 'time'))}; "
                     f"NAMD {_v(_binary_tail(_site_value(site, 'namd_binary')))}",
                     "md_job (each value's mark and source are in the manifest)"])
    else:
        rows.append(["Simulation job", "no md_job recorded for this pocket", "—"])
    doc.h(4, "Configuration (as recorded when the package and the job were written)")
    doc.table(["Field", "Value", "Recorded by"], rows, (0.24, 0.56, 0.2))
    doc.p("Quoted from the md_package and md_job blocks of the manifest. They record what was WRITTEN; "
          "whether the cluster ran these inputs unchanged is not recorded there. Each analysed run below "
          "names its own NAMD input, and its time axis comes from the trajectory, not from this table.",
          small=True)


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


def _md_figure(full: dict, stride: int):
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
        pw = width / cols
        d = Drawing(width, 2 * (ph + top))
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


def _md_stage(doc, prep: Path, label: str, rec: dict):
    res = rec.get("result") or {}
    traj = res.get("trajectory") or {}
    inputs = ", ".join(traj.get("stage_inputs") or []) or "no NAMD input of this stage was found"
    doc.h(4, f"{label.capitalize()} run (stage {_v(rec.get('stage'))}; input {inputs})")
    doc.p(_control_line(res))
    doc.p(f"Recorded {_v(rec.get('generated_at'))} by {_v(rec.get('generated_by'))}. Every figure under "
          f"this heading comes from the {label} run only. The word '{label}' is the label the recorded "
          f"job gives this stage: {_v(rec.get('stage_label_source'))}.", small=True)

    top, lig, pocket = rec.get("topology") or {}, res.get("ligand") or {}, res.get("pocket") or {}
    restr = traj.get("restraints") or {}
    held = restr.get("held") or {}
    com = res.get("com_distance") or {}
    items = [f"Trajectory {t.get('file')}: {_v(t.get('frames'))} frames, {_v(t.get('size_bytes'))} bytes, "
             f"SHA-256 {_v(t.get('sha256'))}" for t in rec.get("trajectories") or []]
    items += [
        f"Topology {_v(top.get('file'))}, SHA-256 {_v(top.get('sha256'))}"
        + (f"; equal to the PSF of the registered package ({(top.get('reference') or {}).get('source')})"
           if top.get("matches_registered_package") else ""),
        f"{_v(traj.get('frames_note'))}; frame spacing {_v(traj.get('dt_ps'))} ps; span {_v(traj.get('span_ns'))} "
        "ns between the first and the last analysed frame",
        f"Minimization: {_v(traj.get('minimization_note'))}",
        f"Restraints: {_v(restr.get('note'))}",
        f"Held atoms among those measured: ligand heavy atoms: "
        f"{_v((held.get('ligand_heavy') or {}).get('note'))}; protein backbone: "
        f"{_v((held.get('protein_backbone') or {}).get('note'))}",
        f"Ligand: {_v(lig.get('selection'))}; {_v(lig.get('dropped_note'))}",
        f"Pocket: {_v(pocket.get('selection'))} ({_v(pocket.get('criterion'))}): "
        + (", ".join(pocket.get("residues") or []) or NO_RECORD),
        f"Built residues in the pocket: {_v(pocket.get('modelled_residues_in_pocket_note'))}",
        f"Periodic images: {_v(com.get('pbc_corrected_frames'))} of {_v(com.get('n'))} frames re-imaged; "
        f"non-rectangular boxes {_v(com.get('non_orthorhombic_frames'))}; frames without a box "
        f"{_v(com.get('frames_without_box'))}; frames with the pocket split across the box "
        f"{_v(com.get('pocket_split_frames'))}",
        "Arguments: " + (", ".join(f"{k} {v}" for k, v in (rec.get("arguments") or {}).items()
                                   if v not in (None, "")) or NO_RECORD),
    ]
    doc.bullets(items)

    doc.h(4, f"Measured over the {label} run (calculated)")
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
    doc.h(4, f"Hydrogen bonds, {label} run")
    doc.p(f"{_v(hb.get('pairs_found'))} ligand–protein pair(s) met the criterion in at least one analysed "
          f"frame. Criterion: {_v(hb.get('criterion'))}. The share is a PROPORTION OF ANALYSED FRAMES, "
          "never a fraction of time: hydrogen bonds break and re-form faster than a frame is stored.",
          small=True)
    doc.table(["Donor", "Hydrogen", "Acceptor", "Frames present / analysed", "Share of analysed frames"],
              [[p.get("donor"), p.get("hydrogen"), p.get("acceptor"),
                f"{_v(p.get('frames_present'))} / {_v(p.get('frames_analysed'))}",
                _num(p.get("occupancy_sampled"), "{:.2f}") or NO_RECORD] for p in pairs],
              (0.26, 0.12, 0.26, 0.18, 0.18))

    rmsf = res.get("rmsf") or {}
    doc.h(4, f"Pocket flexibility (RMSF), {label} run")
    doc.p(f"{_v(rmsf.get('selection'))}. {_v(rmsf.get('criterion'))}. With no apo simulation, a flexible "
          "residue cannot be attributed to the ligand.", small=True)
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
        doc.figure(_md_figure(full, stride),
                   f"{label.capitalize()} run, stage {_v(rec.get('stage'))}: drawn from "
                   f"{(rec.get('series_file') or {}).get('path')}, whose SHA-256 matches the manifest; "
                   f"nothing is recomputed. {sampling.capitalize()}. Solid lines are full segments, dashed "
                   "lines the observed residues only (post hoc, not pre-registered). The hydrogen-bond "
                   "panel counts bonds per stored frame, a sample, not a fraction of time; its count "
                   "is recorded without times, so it is placed on the ligand series' time axis (same "
                   "analysed frames, same length), and left out if the lengths differ.")

    reading = res.get("reading") or {}
    doc.h(4, f"Offered reading of the {label} run (not measured)")
    doc.p(f"{_v(reading.get('status'))}. {_v(reading.get('disclaimer'))} Thresholds: "
          + (", ".join(f"{k} {v}" for k, v in (reading.get("thresholds_used") or {}).items()) or NO_RECORD)
          + f". Source of the thresholds: {_v(reading.get('thresholds_source'))}.", small=True)
    doc.bullets(list(reading.get("offered") or []))
    doc.p("Not measured: " + " ".join(res.get("not_measured") or [NO_RECORD]), small=True)


def _md(doc, m: dict, accession: str, prep: Path) -> int:
    doc.h(2, "2. Molecular dynamics")
    rows = _md_rows(m)
    if not rows:
        doc.p("Not run. No molecular dynamics package, job or analysis is recorded in prep/manifest.json. "
              "The chain is build_md_package, CHARMM-GUI, build_md_job, the cluster, and then "
              f"analyze_md(folder, uniprot_id=\"{accession}\", structure, pocket), which records the "
              "analysis; build this record again after it.")
        return 0
    runs = _md_runs(m)
    names = ", ".join(f"{s.get('pdb_id')} {k} {label}" for s, k, label, _r in runs) or "none"
    doc.p(f"{len(rows)} pocket{'' if len(rows) == 1 else 's'} with a molecular dynamics record; "
          f"{len(runs)} analysed run{'' if len(runs) == 1 else 's'} ({names}). " + MD_VOCABULARY,
          small=True)
    for s, key, package, job, analysis in rows:
        doc.h(3, f"{s.get('pdb_id')} · pocket {key}")
        _md_config(doc, package, job)
        recorded = [lab for lab in MD_STAGES if isinstance(analysis.get(lab), dict)]
        if not recorded:
            doc.p("No trajectory analysis is recorded for this pocket. Run analyze_md on the folder "
                  f"brought back from the cluster with uniprot_id=\"{accession}\" to record one.",
                  small=True)
        for label in recorded:
            _md_stage(doc, prep, label, analysis[label])
    return len(runs)


def _methods(doc, m: dict):
    doc.h(2, "4. Methods and limitations")
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


def build_document(accession: str, m: dict, generated: str, prep: Path = None) -> tuple:
    """(Doc, number of docking runs recorded, number of ADMET results recorded).

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
    n = _docking(doc, m, accession)
    _md(doc, m, accession, prep or choice.manifest_path(accession).parent)
    n_admet = _admet(doc, m, accession)
    _methods(doc, m)
    return doc, n, n_admet


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
    doc, n_runs, n_admet = build_document(accession, m, generated)
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
                            + (["docking (no run in the manifest)"] if not n_runs else []),
        **({"font_warning": font_note} if font_note else {}),
    }
