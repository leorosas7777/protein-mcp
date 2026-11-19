"""build_console_page: the dossier and the experiments as one interactive page (Phase 3, step 10).

The written record stays what it is: report_<acc>.pdf holds the context and experiments_<acc>.pdf
the experiments, both archivable and citable. This page is the companion the PDFs cannot be -- it
holds everything the documents must summarise (all 500 associations, all 47 variants, every pose
of every run) and lets the reader open only what they want. It is written beside the PDFs and
links to them, so the detail is always one click away.

What it does NOT do, on purpose:
  * no network, at build time or at open time. The data is embedded in the file because a page
    opened with file:// cannot fetch a sibling file, and no script is loaded from a CDN;
  * no computing. Every number comes from a dossier JSON or from prep/manifest.json, exactly as
    the tools recorded it. The page re-derives nothing, so it cannot disagree with the PDFs;
  * no ranking and no verdict. Compounds are listed in the order the manifest holds them, with
    the score resolution printed beside them; the reader decides.
  * no launching. A web page cannot start a program on the machine that opens it, so the PyMOL
    buttons print the exact open_in_pymol call instead of pretending to run it.

Imports dossier, report (for load_manifest) and files.choice; never a source module.
"""
import json
from pathlib import Path

from protein_mcp import dossier, report
from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.console_assets import TEMPLATE
from protein_mcp.files import choice


def _pick(d, *path, default=None):
    for key in path:
        if not isinstance(d, dict):
            return default
        d = d.get(key)
    return default if d is None else d


def _load(folder: Path, name: str):
    p = folder / name
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except ValueError as exc:
        log.warning("console: %s could not be read (%s)", name, exc)
        return None


def _tool_file(folder: Path, base: str):
    """(the document, its file name) for one tool of the dossier.

    The default call writes `<base>.json` and a call with a mutation, a topic or an order writes
    `<base>_<something>.json` (dossier.dossier_name). The default call comes first, as in
    report.py: a dossier that holds both would otherwise be shown through the filtered query
    while the page says nothing about which one it read. The file name travels with the data so
    the page can name it.
    """
    default = folder / f"{base}.json"
    if default.is_file():
        return _load(folder, default.name), default.name
    for p in sorted(folder.glob(f"{base}_*.json")):
        return _load(folder, p.name), p.name
    return None, None


# --------------------------------------------------------------------------- phase 1
def _context(folder: Path) -> dict:
    pdoc, pfile = _tool_file(folder, "profile")
    pdoc = pdoc or {}
    prof = _pick(pdoc, "result", default={})
    uni, ot = prof.get("uniprot") or {}, prof.get("open_targets") or {}
    sdoc, sfile = _tool_file(folder, "structures")
    ldoc, lfile = _tool_file(folder, "ligands")
    tdoc, tfile = _tool_file(folder, "literature")
    sdoc, ldoc, tdoc = sdoc or {}, ldoc or {}, tdoc or {}
    exp = _pick(sdoc, "result", "experimental", default={})
    af = _pick(sdoc, "result", "alphafold", default={})
    lig = _pick(ldoc, "result", default={})
    lit = _pick(tdoc, "result", default={})

    ligands = []
    for m in lig.get("ligands") or []:
        ligands.append({
            "id": m.get("chembl_id"), "name": m.get("name"), "link": m.get("link"),
            "phase": m.get("max_phase"), "approval": m.get("first_approval"),
            "moa": _pick(m, "mechanism", "mechanism_of_action"),
            "pchembl": m.get("best_pchembl"), "n": m.get("n_measurements"),
            "mw": _pick(m, "identity", "molecular_weight"),
            "alogp": _pick(m, "properties", "alogp"),
            "psa": _pick(m, "properties", "psa"),
            "ro5": _pick(m, "properties", "ro5_violations"),
            "inchikey": _pick(m, "identity", "inchi_key"),
            "potency": [{"mutation": p.get("mutation"), "type": p.get("standard_type"),
                         "assay": p.get("assay_class"), "n": p.get("n_measurements"),
                         "best": p.get("best_nM"), "median": p.get("median_nM"),
                         "worst": p.get("worst_nM"),
                         "example": (p.get("assay_examples") or [None])[0]}
                        for p in m.get("potency") or []],
        })
    return {
        "identity": {
            "name": uni.get("protein_name"), "gene": uni.get("gene"),
            "organism": uni.get("organism"), "length": uni.get("length_aa"),
            "url": uni.get("source_url"), "entry": uni.get("entry_name"),
            "ensembl": ot.get("ensembl_gene_id"), "ot_release": ot.get("data_version"),
        },
        "retrieved": {
            "profile": pdoc.get("retrieved_at"),
            "structures": sdoc.get("retrieved_at"), "ligands": ldoc.get("retrieved_at"),
            "literature": tdoc.get("retrieved_at"),
        },
        # which file of the dossier each panel was read from: with a mutation filter there is
        # more than one per tool, and they do not hold the same rows
        "files": {"profile": pfile, "structures": sfile, "ligands": lfile, "literature": tfile},
        "function": [f.get("text") for f in uni.get("function") or []],
        "variants": [{"change": v.get("change"), "pos": v.get("start"),
                      "description": v.get("description") or "",
                      "dbsnp": (v.get("dbsnp") or [None])[0],
                      "pdb": None, "chembl": None, "papers": None}
                     for v in uni.get("natural_variants") or []],
        "assoc": {"total": _pick(ot, "associations", "total"),
                  "examined": _pick(ot, "associations", "examined"),
                  "selection": _pick(ot, "associations", "selection"),
                  "rows": [{"id": a.get("disease_id"), "name": a.get("name"),
                            "score": a.get("score"), "onc": bool(a.get("is_oncology")),
                            "dts": a.get("datatype_scores") or {}}
                           for a in _pick(ot, "associations", "rows", default=[])]},
        "tractability": {"flags": ot.get("tractability") or {},
                         "names": ot.get("modality_names") or {}},
        "structures": {"total": exp.get("total_entries"), "examined": exp.get("examined"),
                       "confirmed": exp.get("confirmed_in_examined"),
                       "selection": exp.get("selection"),
                       "rows": [{"pdb": s.get("pdb_id"), "title": s.get("title"),
                                 "res": s.get("resolution_angstrom"), "date": s.get("release_date"),
                                 "method": s.get("method"),
                                 "ligands": [x.get("id") for x in s.get("ligands") or []],
                                 "rcsb": (s.get("links") or {}).get("rcsb")}
                                for s in exp.get("structures") or []],
                       "dropped": [{"pdb": d.get("pdb_id"), "title": d.get("title")}
                                   for d in exp.get("discarded_not_confirmed") or []]},
        "alphafold": {"id": af.get("entry_id"), "plddt": af.get("mean_plddt"),
                      "band": af.get("confidence_band"),
                      "fractions": af.get("fraction_of_residues") or {}},
        "ligands": {"rows": ligands, "selection": lig.get("selection"),
                    "qualifying": lig.get("ligands_with_qualifying_data"),
                    "lower_bound": lig.get("ligands_with_qualifying_data_is_lower_bound"),
                    "activities": lig.get("activities") or {},
                    "target": (lig.get("target") or {}).get("chembl_id")},
        "literature": {"query": lit.get("query_used"), "total": lit.get("total_hits"),
                       "rows": [{"title": p.get("title"), "journal": p.get("journal"),
                                 "year": p.get("year"), "pmid": p.get("pmid"),
                                 "doi": p.get("doi"), "oa": bool(p.get("open_access")),
                                 "abstract": p.get("abstract") or ""}
                                for p in lit.get("papers") or []]},
    }


# --------------------------------------------------------------------------- phase 2/3
def _compound_of(m: dict, pocket_key: str) -> dict:
    """The compound a pocket key ends in, from ligand_choice. Looked up, never invented."""
    for c in (m.get("ligand_choice") or {}).get("compounds") or []:
        cid = str(c.get("id") or "")
        if cid and pocket_key.upper().endswith(cid.upper()):
            return c
    return {}


def _experiments(m: dict, prep: Path) -> dict:
    pockets, runs, admet, md = [], [], [], []
    centres = {}
    for s in m.get("structures") or []:
        pdb = s.get("pdb_id")
        for pk in (s.get("prepared") or {}).get("pockets") or []:
            if pk.get("error"):
                continue
            box = pk.get("box") or {}
            pockets.append({
                "pdb": pdb,
                "key": " ".join(x for x in (pk.get("comp_id"), pk.get("label_asym_id")) if x),
                "chain": pk.get("auth_asym_id"), "heavy": pk.get("heavy_atoms"),
                "residues": [f"{r.get('resname')}{r.get('auth_seq_id')}" for r in
                             sorted(pk.get("residues") or [],
                                    key=lambda r: r.get("min_distance") or 99)],
                "distances": [r.get("min_distance") for r in
                              sorted(pk.get("residues") or [],
                                     key=lambda r: r.get("min_distance") or 99)],
                "box": box.get("size"), "centre": box.get("center"),
                "large_box": pk.get("large_box"), "spans": pk.get("spans_target_chains"),
            })
            centres[(pdb, pockets[-1]["key"])] = box.get("center")
        for key, block in sorted((s.get("docking") or {}).items()):
            for run_id, run in sorted((block.get("runs") or {}).items()):
                cd = run.get("compound_docking") or {}
                rd = run.get("redocking") or {}
                xd = run.get("cross_docking") or {}
                ligand = run.get("ligand") or {}
                fit = cd.get("box_fit") or {}
                bar = _pick(xd, "reference_pose_uncertainty", "shift_angstrom")
                rng = _pick(xd, "within_criterion_with_bar", "range_angstrom", default=[])
                runs.append({
                    "pdb": pdb, "pocket": key, "run": run_id,
                    "kind": "compound" if cd else "redock",
                    "id": ligand.get("compound_id"),
                    "name": ligand.get("name") or ligand.get("comp_id"),
                    "level": ligand.get("level"),
                    "smiles": ligand.get("smiles"), "inchikey": ligand.get("inchikey"),
                    "heavy": cd.get("heavy_atoms") or ligand.get("heavy_atoms"),
                    "score": cd.get("top_pose_score_kcal_mol"),
                    "le": cd.get("ligand_efficiency_kcal_mol_per_heavy_atom"),
                    "box_limited": bool(cd.get("box_limited")),
                    "clearance": fit.get("min_clearance_angstrom"),
                    "limiting_axis": fit.get("limiting_box_axis_angstrom"),
                    "extent": fit.get("limiting_ligand_extent_angstrom"),
                    "rmsd_na": cd.get("rmsd_not_applicable"),
                    "redock_rmsd": rd.get("top_pose_rmsd_angstrom"),
                    "redock_criterion": rd.get("criterion_angstrom"),
                    "redock_within": rd.get("within_criterion"),
                    "anchor": ({"rmsd": xd.get("top_pose_rmsd_angstrom"), "bar": bar,
                                "criterion": xd.get("criterion_angstrom"),
                                "range": rng,
                                "decides": _pick(xd, "within_criterion_with_bar", "decides"),
                                "reference": _pick(xd, "reference_pose", "from"),
                                "matched_by": _pick(xd, "reference_pose", "matched_by"),
                                "sup_rmsd": _pick(xd, "superposition", "rmsd_angstrom"),
                                "sup_pairs": _pick(xd, "superposition", "pairs")} if xd else None),
                    "poses": [{"mode": p.get("mode"), "score": p.get("vina_score_kcal_mol"),
                               "rmsd": (p.get("rmsd_to_reference_pose_angstrom")
                                        if xd else p.get("rmsd_to_crystal_pose_angstrom"))}
                              for p in run.get("poses") or []],
                    "seconds": _pick(run, "vina", "wall_seconds"),
                    "scope": cd.get("scope") or rd.get("scope"),
                    "geometry": _pose_geometry(prep, run, centres.get((pdb, key)),
                                               ligand.get("name") or ligand.get("comp_id") or run_id),
                })
        for key, block in sorted((s.get("admet") or {}).items()):
            if not isinstance(block, dict):
                continue
            ident = block.get("compound_identity") or {}
            rec = block.get("chembl_record") or {}
            calc = _pick(block, "calculated_properties", "rdkit", default={})
            meas = block.get("measured_adme") or {}
            rules = _pick(block, "rule_reading", "results", default={})
            admet.append({
                "pdb": pdb, "pocket": key,
                "who": ident.get("comp_id") or ident.get("name") or ident.get("compound_id"),
                "id": ident.get("compound_id") or rec.get("chembl_id"),
                "chembl_name": rec.get("pref_name"), "link": rec.get("link"),
                "status": rec.get("status"),
                "phase": _pick(rec, "fields", "max_phase"),
                "approval": _pick(rec, "fields", "first_approval"),
                "withdrawn": _pick(rec, "fields", "withdrawn_flag"),
                "blackbox": _pick(rec, "fields", "black_box_warning"),
                "atc": rec.get("atc_classifications") or [],
                "props": {k: calc.get(k) for k in
                          ("molecular_weight", "clogp", "tpsa", "hbd", "hba",
                           "rotatable_bonds", "qed", "heavy_atoms")},
                "rules": {name: {"violations": r.get("violations"),
                                 "violated": r.get("violated") or []}
                          for name, r in rules.items()},
                "adme_rows": [{"cls": r.get("adme_class"), "type": r.get("standard_type"),
                               "value": r.get("value"), "relation": r.get("relation"),
                               "units": r.get("units"), "assay": r.get("assay_description")}
                              for r in meas.get("adme_rows") or []],
                "rows_read": meas.get("rows_read"), "retrieved": bool(meas.get("retrieved")),
                "meas_note": meas.get("note"),
            })
        pkgs = s.get("md_package") if isinstance(s.get("md_package"), dict) else {}
        jobs = s.get("md_job") if isinstance(s.get("md_job"), dict) else {}
        anas = s.get("md_analysis") if isinstance(s.get("md_analysis"), dict) else {}
        for key in sorted(set(pkgs) | set(jobs) | set(anas)):
            pkg, job = pkgs.get(key) or {}, jobs.get(key) or {}
            comp = _compound_of(m, key)
            checks = _pick(job, "verified", "checks", default={})
            sysm = checks.get("system") or {}
            stages = []
            for label, rec in sorted((anas.get(key) or {}).items()):
                if not isinstance(rec, dict):
                    continue
                res = rec.get("result") or {}
                ctrl = _pick(res, "reading", "equilibration_control", default={})
                stages.append({
                    "label": label, "status": ctrl.get("status"),
                    "drift": ctrl.get("drift_angstrom"), "threshold": ctrl.get("threshold_angstrom"),
                    "segid": ctrl.get("segid"),
                    "span_ns": _pick(res, "trajectory", "span_ns"),
                    "frames": _pick(res, "trajectory", "frames_note"),
                    "restraints": _pick(res, "trajectory", "restraints", "note"),
                    "series": _series(prep, rec),
                    "hbonds": [{"donor": p.get("donor"), "acceptor": p.get("acceptor"),
                                "present": p.get("frames_present"),
                                "analysed": p.get("frames_analysed"),
                                "share": p.get("occupancy_sampled")}
                               for p in _pick(res, "hbonds", "pairs", default=[])],
                    "rmsf": [{"residue": r.get("residue"), "rmsf": r.get("rmsf_angstrom")}
                             for r in _pick(res, "rmsf", "per_residue", default=[])],
                    "offered": list(_pick(res, "reading", "offered", default=[])),
                })
            md.append({
                "pdb": pdb, "pocket": key,
                "compound": comp.get("name") or comp.get("id"),
                "id": comp.get("id"), "level": comp.get("level"),
                "run": pkg.get("run_id"),
                "pose": _pick(pkg, "pose", "mode"),
                "pose_score": _pick(pkg, "pose", "vina_score_kcal_mol"),
                "pose_rmsd": _pick(pkg, "pose", "rmsd_to_crystal_pose_angstrom"),
                "box_side": sysm.get("box_side_angstrom"), "waters": sysm.get("waters"),
                "temperature": sysm.get("temperature_kelvin"),
                "timestep": sysm.get("timestep_fs"),
                "ns_total": _pick(job, "slurm", "nanoseconds_total"),
                "stages": stages,
            })
    return {"pockets": pockets, "runs": runs, "admet": admet, "md": md,
            "choice": {"structure": _pick(m, "user_choice", "structure"),
                       "reason": _pick(m, "user_choice", "reason"),
                       "resolution": _pick(m, "user_choice", "resolution_angstrom")},
            "compound_choice": {
                "reason": _pick(m, "ligand_choice", "reason"),
                "compounds": [{"id": c.get("id"), "name": c.get("name"), "level": c.get("level")}
                              for c in _pick(m, "ligand_choice", "compounds", default=[])]},
            "mutation": m.get("mutation_filter"), "created": m.get("created_at")}


def _pose_geometry(prep: Path, run: dict, box_centre, name: str) -> dict:
    """The top pose's heavy atoms and bonds, in Ångström relative to the box centre.

    Drawn so the central panel shows the SELECTED compound instead of the same picture for all
    of them (the user reported on 2026-11-18 that clicking a compound printed nothing). These
    are the coordinates of the pose file this run wrote, read with RDKit and only shifted to the
    box centre -- no fitting, no idealisation. The file is checked against the SHA-256 in the
    manifest first: a pose file that changed is not the pose this run produced.
    """
    import hashlib
    rec = next((f for n, f in (run.get("files") or {}).items() if n.endswith(".sdf")), None)
    if not rec or not rec.get("path"):
        return {"error": "no pose SDF is recorded for this run"}
    path = (prep / rec["path"]).resolve()
    if prep.resolve() not in path.parents or not path.is_file():
        return {"error": f"{rec['path']} is missing"}
    if rec.get("sha256") and hashlib.sha256(path.read_bytes()).hexdigest() != rec["sha256"]:
        return {"error": f"{rec['path']} no longer matches the SHA-256 in the manifest"}
    try:
        from rdkit import Chem
    except ImportError:
        return {"error": "RDKit is not available"}
    supplier = Chem.SDMolSupplier(str(path), removeHs=True, sanitize=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None or not mol.GetNumConformers():
        return {"error": f"{rec['path']} holds no readable pose"}
    conf = mol.GetConformer(0)
    cx, cy, cz = (list(box_centre) + [0, 0, 0])[:3] if box_centre else (0, 0, 0)
    atoms, keep = [], {}
    for a in mol.GetAtoms():
        if a.GetAtomicNum() <= 1:
            continue
        p = conf.GetAtomPosition(a.GetIdx())
        keep[a.GetIdx()] = len(atoms)
        atoms.append({"e": a.GetSymbol(),
                      "x": round(p.x - cx, 3), "y": round(p.y - cy, 3), "z": round(p.z - cz, 3)})
    bonds = [[keep[b.GetBeginAtomIdx()], keep[b.GetEndAtomIdx()]] for b in mol.GetBonds()
             if b.GetBeginAtomIdx() in keep and b.GetEndAtomIdx() in keep]
    shifted = "shifted to the box centre" if box_centre else         "in the frame of the file: this pocket records no box centre to shift them to"
    return {"atoms": atoms, "bonds": bonds, "file": rec["path"], "name": name,
            "note": f"heavy atoms of the top-scoring pose, {shifted} and projected; "
                    "no fitting and no idealisation"}


def _series(prep: Path, rec: dict) -> dict:
    """The per-frame series of one analysed stage, read from the recorded file ONLY if its
    SHA-256 still matches the manifest. A file that changed is not the file that was measured."""
    import hashlib
    sf = rec.get("series_file") or {}
    if not sf.get("path") or not sf.get("sha256"):
        return {"error": "no series file is recorded for this run"}
    path = (prep / sf["path"]).resolve()
    if prep.resolve() not in path.parents or not path.is_file():
        return {"error": f"{sf['path']} is missing or outside the prep folder"}
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != sf["sha256"]:
        return {"error": f"{sf['path']} no longer matches the SHA-256 in the manifest"}
    try:
        full = json.loads(data.decode("utf-8")).get("result") or {}
    except ValueError as exc:
        return {"error": f"{sf['path']} could not be read ({exc})"}
    out = {}
    prot = full.get("protein_rmsd") or {}
    for seg, ser in sorted((prot.get("per_segment") or {}).items()):
        out[f"backbone {seg}"] = _line(ser)
        if isinstance(ser.get("observed_only"), dict):
            out[f"backbone {seg} observed"] = _line(ser["observed_only"])
    axis = (full.get("ligand_rmsd") or {}).get("times_ps") or []
    out["ligand RMSD"] = _line(full.get("ligand_rmsd") or {})
    out["COM distance"] = _line(full.get("com_distance") or {})
    hb = dict(_pick(full, "hbonds", "count_per_frame", default={}))
    if not hb.get("times_ps") and len(hb.get("values") or []) == len(axis):
        hb["times_ps"] = axis
    out["H-bonds"] = _line(hb)
    return {k: v for k, v in out.items() if v}


def _line(ser: dict) -> dict:
    values, times = ser.get("values") or [], ser.get("times_ps") or []
    if not values or len(values) != len(times):
        return {}
    return {"v": [round(float(x), 4) for x in values],
            "t": [round(float(x) / 1000.0, 4) for x in times],
            "min": ser.get("min"), "max": ser.get("max"), "drift": ser.get("drift")}



# --------------------------------------------------------------------------- 3D
VIEWER_ASSET = Path(__file__).with_name("assets") / "3Dmol-min.js"
# Identified by SHA-256 against the published build, not by the file name: see
# protein_mcp/assets/PROVENANCE.md, which also holds the two licence files.
VIEWER_NAME = "3Dmol.js 2.5.5"
VIEWER_LICENSE = "BSD-3-Clause (University of Pittsburgh and contributors)"
VIEWER_NOTICE = """/*! 3Dmol.js 2.5.5 - molecular visualisation library by David Koes and contributors.
 * Copyright (c) 2014, University of Pittsburgh and contributors. All rights reserved.
 * Redistributed here under the BSD-3-Clause licence; its full text and the provenance of this
 * build are in protein_mcp/assets/ (3Dmol-LICENSE.txt, 3Dmol-min.js.LICENSE.txt,
 * PROVENANCE.md) of the protein-research MCP server that wrote this page. */
"""
# a receptor is ~700 KB of text and the page carries it whole, so only the structures that
# actually have a docking run are embedded; the rest are named and not loaded
MAX_EMBEDDED_RECEPTORS = 4


def _first_model(sdf_text: str) -> str:
    """The first record of a multi-pose SDF. The others are in the file; the viewer shows the
    top-scoring one, which is the pose every figure in the report is about."""
    head, sep, _rest = sdf_text.partition("$$$$")
    return (head + sep + "\n") if sep else sdf_text


def _structures_3d(m: dict, prep: Path) -> dict:
    """Receptor and top-pose text for the structures that carry a docking run.

    Embedded as text because a page opened with file:// cannot read a sibling file. Nothing is
    converted: the receptor is the cleaned PDB the preparation step wrote and each pose is the
    first record of the SDF the docking run wrote, both verbatim.
    """
    import hashlib
    receptors, poses, skipped, refused = {}, [], [], []
    with_runs = [s for s in m.get("structures") or [] if s.get("docking")]
    for s in with_runs:
        pdb = s.get("pdb_id")
        rec = (s.get("prepared") or {}).get("receptor") or {}
        path = (prep / rec["path"]).resolve() if rec.get("path") else None
        if path is None or prep.resolve() not in path.parents or not path.is_file():
            skipped.append({"pdb": pdb, "reason": "no cleaned receptor in prep/"})
        elif len(receptors) >= MAX_EMBEDDED_RECEPTORS:
            skipped.append({"pdb": pdb, "reason": f"more than {MAX_EMBEDDED_RECEPTORS} receptors "
                                                  "would be embedded; this one was left out"})
        else:
            receptors[pdb] = path.read_text(encoding="utf-8", errors="replace")
        for key, block in sorted((s.get("docking") or {}).items()):
            for run_id, run in sorted((block.get("runs") or {}).items()):
                sdf = next((f for n, f in (run.get("files") or {}).items() if n.endswith(".sdf")), None)
                if not sdf or not sdf.get("path"):
                    continue
                f = (prep / sdf["path"]).resolve()
                if prep.resolve() not in f.parents or not f.is_file():
                    refused.append({"run": run_id, "file": sdf["path"],
                                    "reason": "missing from prep/"})
                    continue
                data = f.read_bytes()
                # the same rule as the flat panel and as every PyMOL session: a pose file that
                # no longer matches the manifest is not the pose that was scored, and the 3D
                # view is exactly where a wrong molecule would look most convincing
                if not sdf.get("sha256"):
                    refused.append({"run": run_id, "file": sdf["path"],
                                    "reason": "no SHA-256 was recorded for it"})
                    continue
                if hashlib.sha256(data).hexdigest() != sdf["sha256"]:
                    refused.append({"run": run_id, "file": sdf["path"],
                                    "reason": "no longer matches the SHA-256 in the manifest"})
                    continue
                lig = run.get("ligand") or {}
                poses.append({
                    "run": run_id, "pdb": pdb, "pocket": key,
                    "name": lig.get("name") or lig.get("compound_id") or lig.get("comp_id") or run_id,
                    "id": lig.get("compound_id"), "level": lig.get("level"),
                    "kind": "compound" if run.get("compound_docking") else "redock",
                    "sdf": _first_model(data.decode("utf-8", errors="replace")),
                })
    return {"receptors": receptors, "poses": poses, "receptors_not_embedded": skipped,
            "poses_not_embedded": refused,
            "note": "The receptor is the cleaned PDB of the preparation step and each pose is the "
                    "first record of that run's SDF, both embedded verbatim. Every pose file was "
                    "checked against the SHA-256 in the manifest first. The other poses of a run "
                    "are in its file and are not shown here. The cleaned receptor carries no "
                    "SHA-256 in the manifest, so it is the one file here that cannot be checked."}


# --------------------------------------------------------------------------- build
def build_console_page(uniprot_id: str) -> dict:
    """dossier + prep/manifest.json -> console_<accession>.html beside the PDFs."""
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    folder = dossier.DOSSIER_DIR / accession
    if not folder.is_dir():
        return {"ok": False, "uniprot_id": accession,
                "error": f"No dossier for {accession}. Run get_protein_profile first."}
    m = report.load_manifest(accession)
    if m is None:
        m = {}
        log.info("console: no preparation manifest for %s; phase 2/3 panels will say so", accession)
    elif m.get("unreadable"):
        return {"ok": False, "uniprot_id": accession,
                "error": f"prep/manifest.json could not be read ({m['unreadable']})."}
    prep = choice.manifest_path(accession).parent
    data = {
        "accession": accession,
        "built_at": utc_now(),
        "context": _context(folder),
        "experiments": _experiments(m, prep) if m else None,
        "three_d": _structures_3d(m, prep) if m else None,
        "documents": [
            {"label": "Context report (phase 1)", "file": f"report_{accession}.pdf",
             "exists": (folder / f"report_{accession}.pdf").is_file()},
            {"label": "Experiment record (phases 2-3)", "file": f"experiments_{accession}.pdf",
             "exists": (folder / f"experiments_{accession}.pdf").is_file()},
        ],
    }
    try:
        library = VIEWER_ASSET.read_text(encoding="utf-8")
    except OSError as exc:
        library, viewer_note = "", f"3D viewer not embedded: {VIEWER_ASSET.name} could not be read ({exc})"
    else:
        # the page redistributes the bundle, and a BSD-3-Clause redistribution has to carry the
        # copyright notice: it goes in the page, not only in protein_mcp/assets/PROVENANCE.md
        library = VIEWER_NOTICE + library
        viewer_note = (f"{VIEWER_NAME} embedded from {VIEWER_ASSET.name} "
                       f"({len(library) // 1024} KB), {VIEWER_LICENSE}")
    data["viewer_credit"] = (f"3D viewer: {VIEWER_NAME}, {VIEWER_LICENSE}"
                             if library else "no 3D viewer embedded")
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    out = folder / f"console_{accession}.html"
    out.write_text(TEMPLATE.replace("/*__DATA__*/null", payload)
                           .replace("/*__VIEWER_LIB__*/", library), encoding="utf-8")
    exps = data["experiments"] or {}
    return {
        "ok": True, "uniprot_id": accession,
        "file": str(out),
        "open_with": "a web browser; it needs no server and no network",
        "size_kb": round(out.stat().st_size / 1024, 1),
        "holds": {
            "associations": len(_pick(data, "context", "assoc", "rows", default=[])),
            "variants": len(_pick(data, "context", "variants", default=[])),
            "ligands_retrieved": len(_pick(data, "context", "ligands", "rows", default=[])),
            "docking_runs": len(exps.get("runs") or []),
            "admet_results": len(exps.get("admet") or []),
            "md_pockets": len(exps.get("md") or []),
            "md_stages_measured": sum(len(x.get("stages") or []) for x in exps.get("md") or []),
        },
        "viewer": viewer_note,
        "embedded_3d": {
            "receptors": sorted((data["three_d"] or {}).get("receptors") or {}),
            "poses": len((data["three_d"] or {}).get("poses") or []),
            "receptors_not_embedded": (data["three_d"] or {}).get("receptors_not_embedded") or [],
            "poses_not_embedded": (data["three_d"] or {}).get("poses_not_embedded") or [],
        },
        "dossier_files_read": _pick(data, "context", "files", default={}),
        "documents_linked": [d["file"] for d in data["documents"] if d["exists"]],
        "documents_missing": [d["file"] for d in data["documents"] if not d["exists"]],
        "notes": [
            "Every value is quoted from the dossier JSON files and prep/manifest.json; the page "
            "computes nothing, so it cannot disagree with the PDFs.",
            "No network at build time or at open time: the data is embedded because a page opened "
            "with file:// cannot read a sibling file, and no script is loaded from a CDN.",
            "The PyMOL buttons print the open_in_pymol call to run; a web page cannot start a "
            "program on the machine that opens it.",
        ],
    }


if __name__ == "__main__":  # python -m protein_mcp.console P15056
    import sys
    print(json.dumps(build_console_page(sys.argv[1] if len(sys.argv) > 1 else "P15056"),
                     indent=1, ensure_ascii=False))
