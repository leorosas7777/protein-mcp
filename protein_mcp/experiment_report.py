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
    ("3. ADMET", "No tool writes this section yet (Phase 3 step 7). ADMET values will be PREDICTED values and will "
                 "be kept in their own section, never in the same table as the calculated docking scores or as "
                 "measured potency."),
]

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
          + ". Molecular dynamics and ADMET (sections 2 and 3) have not been run.")


def build_document(accession: str, m: dict, generated: str) -> tuple:
    """(Doc, number of docking runs recorded)."""
    doc = Doc()
    doc.h(1, f"In-silico experiment record: {accession}")
    seq = m.get("sequence") or {}
    doc.p(f"Generated {generated} from prep/manifest.json (created {_v(m.get('created_at'))} by "
          f"{_v(m.get('generated_by'))}). Variant filter: {m.get('mutation_filter') or 'none'}. "
          f"UniProt sequence used for the numbering: {seq.get('header') or NO_RECORD}", small=True)
    doc.p("Results of in-silico experiments (calculated docking scores and RMSD values; predicted ADMET values once "
          "that step exists), each with the provenance of the structure it came from, quoted from the manifest. "
          "Measured potency (ChEMBL) is not part of this document; it is in the research report, "
          f"report_{accession}.pdf. Interpretation is not part of either.", small=True)
    n = _docking(doc, m, accession)
    for title, text in PENDING_SECTIONS:
        doc.h(2, title)
        doc.p("Not run. " + text)
    _methods(doc, m)
    return doc, n


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
    doc, n_runs = build_document(accession, m, generated)
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
        "sections_not_run": ["molecular dynamics (Phase 3 step 4)", "ADMET (Phase 3 step 7)"]
                            + (["docking (no run in the manifest)"] if not n_runs else []),
        **({"font_warning": font_note} if font_note else {}),
    }
