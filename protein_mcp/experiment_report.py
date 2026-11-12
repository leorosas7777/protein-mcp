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
from protein_mcp import dossier, report
from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice
from protein_mcp.report import NO_RECORD, Doc, Lines, Link, _box, _num

# Sections that have no producing tool yet: declared gaps, never filled from the conversation.
PENDING_SECTIONS = [
    ("2. Molecular dynamics", "No tool writes this section yet. Phase 3 step 4 is to build a delivery package "
                              "for CHARMM-GUI and a SLURM template (the simulation runs outside this machine and "
                              "comes back later), and this section stays empty until a run is recorded."),
]
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
          + ("; ".join(m.get("pending") or []) or "none listed")
          + ". Molecular dynamics (section 2) has not been run.")


def build_document(accession: str, m: dict, generated: str) -> tuple:
    """(Doc, number of docking runs recorded, number of ADMET results recorded)."""
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
    for title, text in PENDING_SECTIONS:
        doc.h(2, title)
        doc.p("Not run. " + text)
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
        "sections_not_run": ["molecular dynamics (Phase 3 step 4)"]
                            + (["ADMET (no analyze_admet result in the manifest)"] if not n_admet else [])
                            + (["docking (no run in the manifest)"] if not n_runs else []),
        **({"font_warning": font_note} if font_note else {}),
    }
