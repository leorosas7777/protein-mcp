"""cross_dock_rmsd: the RMSD of a level-T1 compound against its crystallographic pose in ANOTHER
structure of the same target (Phase 3, step 6b).

A compound of the researcher's selection is docked into the pocket of a structure that does not
contain it, so its pose cannot be compared with anything in that structure. For a T1 compound a
comparison still exists: the same compound HAS a crystallographic pose in another entry of this
dossier. Bringing that pose into the frame of the docked receptor is a superposition, and the
superposition is itself a measurement with an error -- which is why every number here travels with
it.

Pre-registered in ESTADO 7.4 (written before this module existed) and enforced here:
- the superposition is the GLOBAL one (`gemmi.calculate_superposition`, `SupSelect.CaP`), not the
  one restricted to the pocket. Measured 2026-11-17 on 3OG7 -> 8C7Y: global 1.657 A over 247 CA,
  pocket-restricted 1.709 A, i.e. the restricted fit does NOT improve on the global one, and the
  reference pose MOVES 0.76 A between the two;
- that 0.76 A travels with the RMSD everywhere, and the reporting of the number without it is the
  failure this module exists to prevent;
- the criterion is <= 2.0 A, the same as 7.1, and it is declared in advance to be AT THE EDGE OF
  THE NOISE here: the pocket itself differs ~1.7 A between the two crystals, so passing is not
  strong evidence and failing is not a defect of the pipeline, but the expected case of
  cross-docking between two conformations of the same kinase;
- the reference instance is found by InChIKey, never by a compound name (anti-hallucination rule),
  and where it lands is then checked BY GEOMETRY against the box and the pocket ligand;
- the RMSD is rdMolAlign.CalcRMS (in place, symmetry-corrected), never GetBestRMS, with the same
  negative control as dock.py: the reference pose moved 5 A must measure 5 A.

Nothing here has a network. The result is written into the docking run it belongs to, as a
`cross_docking` block, so the number cannot be separated from the run that produced it.

Pure function + CLI: python -m protein_mcp.files.crossdock P15056 --compound CHEMBL1229517
Imports common, files.choice and files.dock (+ gemmi, rdkit, numpy).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice, dock

# ESTADO 7.4. The bar is a measurement of this project (3OG7 -> 8C7Y, 2026-11-17), not a constant
# of the method: it is recomputed per pair below and this value is what was pre-registered.
PREREGISTERED_SHIFT_BAR = 0.76      # A, global vs pocket-restricted superposition
CRITERION = 2.0                     # A, the criterion of 7.1
MIN_POLYMER = 50                    # residues: shorter chains are not superposed
POCKET_CUTOFF = 5.0                 # A, heavy atom to heavy atom, for the pocket-restricted fit

NOISE_NOTE = (
    "Declared in advance (ESTADO 7.4): this criterion is AT THE EDGE OF THE NOISE for a "
    "cross-docking. The pocket itself differs by about the same amount between the two crystals, "
    "so a value below the criterion is NOT strong evidence that the pose is right, and a value "
    "above it is NOT a defect of the pipeline: it is the expected case of docking into a different "
    "conformation of the same protein. Read it together with the superposition RMSD and the shift "
    "bar, never alone.")
BAR_NOTE = (
    "The reference pose is not an observation in this receptor's frame: it is the pose observed in "
    "another crystal, moved here by a superposition. Superposing on a different set of atoms moves "
    "it, and the shift measured between the global and the pocket-restricted sets is the bar below. "
    "Any report of the RMSD must print the bar beside it.")
METHOD = ("rdMolAlign.CalcRMS, heavy atoms, symmetry-corrected, in place (no superposition of the "
          "two ligands: the only superposition applied is the protein one, reported here)")


class CrossDockError(Exception):
    def __init__(self, message: str, **extra):
        super().__init__(message)
        self.extra = extra


# ------------------------------------------------------------------ pure pieces

def _polymers(path: Path):
    """The structure with ligands and waters removed, as gemmi reads it."""
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_ligands_and_waters()
    return st


def best_superposition(target_path: Path, mobile_path: Path):
    """The chain pair of (mobile -> target) with the lowest CA RMSD, and every pair tried.

    Every pair is returned, not only the winner: the choice of chain is a decision of this tool
    and the reader has to be able to see what it chose between.
    """
    import gemmi
    tgt, mob = _polymers(target_path), _polymers(mobile_path)
    tried, best = [], None
    for cm in mob[0]:
        pm = cm.get_polymer()
        if len(pm) < MIN_POLYMER:
            continue
        for ct in tgt[0]:
            pt = ct.get_polymer()
            if len(pt) < MIN_POLYMER:
                continue
            sup = gemmi.calculate_superposition(pt, pm, pm.check_polymer_type(), gemmi.SupSelect.CaP)
            tried.append({"mobile_chain": cm.name, "target_chain": ct.name,
                          "rmsd_angstrom": round(float(sup.rmsd), 3), "pairs": int(sup.count)})
            if best is None or sup.rmsd < best[0].rmsd:
                best = (sup, cm.name, ct.name)
    if best is None:
        raise CrossDockError(f"No polymer chain of {MIN_POLYMER} residues or more could be "
                             f"superposed between {target_path.name} and {mobile_path.name}.")
    return best, sorted(tried, key=lambda d: d["rmsd_angstrom"])


def _ca_by_number(path: Path, chain: str) -> dict:
    out = {}
    for ch in _polymers(path)[0]:
        if ch.name != chain:
            continue
        for res in ch:
            a = res.find_atom("CA", "*")
            if a:
                out[res.seqid.num] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return out


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Superpose P on Q. (R, t, rmsd). Used ONLY for the pocket-restricted fit of the error bar;
    the fit that moves the reference pose is gemmi's."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1, 1, float(np.sign(np.linalg.det(Vt.T @ U.T)))])
    R = Vt.T @ D @ U.T
    t = qc - R @ pc
    rms = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean()))
    return R, t, rms


def heavy_xyz(mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                     for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() > 1])


def moved_copy(mol, transform):
    """A copy of `mol` with every atom moved by a gemmi.Transform (the protein superposition)."""
    import gemmi
    from rdkit import Chem
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    for i in range(out.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        q = transform.apply(gemmi.Position(p.x, p.y, p.z))
        conf.SetAtomPosition(i, (q.x, q.y, q.z))
    return out


def in_box(xyz: np.ndarray, box: dict) -> dict:
    """Whether the moved reference pose falls in the docking box. Geometry, never a name."""
    c, s = np.array(box["center"], float), np.array(box["size"], float)
    low, high = c - s / 2, c + s / 2
    inside = int(((xyz >= low) & (xyz <= high)).all(1).sum())
    return {"heavy_atoms_in_box": f"{inside} of {len(xyz)}",
            "centroid_in_box": bool(((xyz.mean(0) >= low) & (xyz.mean(0) <= high)).all()),
            "centroid_to_box_centre_angstrom": round(float(np.linalg.norm(xyz.mean(0) - c)), 2)}


# ------------------------------------------------------------------ the tool

def _reference_instance(manifest: dict, prep: Path, entry: dict, exclude_pdb: str, named: str = ""):
    """The crystallographic instance of THIS compound in another structure, found by InChIKey.

    Never by comp_id and never by name: a three-letter code is not an identity (ESTADO 7.1). Every
    instance that could be read is reported, so a second candidate cannot be chosen silently.
    """
    from rdkit import Chem

    want = (entry.get("inchikey") or "").strip()
    if not want:
        raise CrossDockError(f"select_ligands recorded no InChIKey for {entry['id']}, so the "
                             "crystallographic instance of this compound cannot be identified. "
                             "Record the selection again.")
    hits, unreadable = [], []
    for s in manifest.get("structures") or []:
        pdb = (s.get("pdb_id") or "").upper()
        if pdb == exclude_pdb.upper():
            continue                     # the pose must come from ANOTHER crystal
        if named and pdb != named.upper():
            continue
        for ni in s.get("nonpolymer_instances") or []:
            sdf = ni.get("sdf") or {}
            if not sdf.get("ok") or not sdf.get("path"):
                continue
            path = (prep / sdf["path"]).resolve()
            if prep.resolve() not in path.parents or not path.is_file():
                continue
            try:
                with open(path, "rb") as fh:
                    mol = next(Chem.ForwardSDMolSupplier(fh, removeHs=True), None)
                key = Chem.MolToInchiKey(mol) if mol is not None else None
            except Exception as exc:                       # a file RDKit cannot read is not a match
                unreadable.append(f"{pdb} {ni.get('comp_id')} {ni.get('label_asym_id')}: {exc!r}")
                continue
            if key and key == want:
                hits.append({"pdb_id": pdb, "comp_id": ni.get("comp_id"),
                             "label_asym_id": ni.get("label_asym_id"),
                             "sdf": sdf["path"], "sha256": sdf.get("sha256"),
                             "inchikey": key, "raw_file": (s.get("file") or {}).get("path")})
    if not hits:
        raise CrossDockError(
            f"No structure of this dossier holds a crystallographic instance of {entry['id']} "
            f"(InChIKey {want}) outside {exclude_pdb}, so there is no reference pose and no RMSD. "
            "This is 'not applicable', not a missing measurement."
            + (f" {len(unreadable)} ligand file(s) could not be read: {'; '.join(unreadable)}"
               if unreadable else ""),
            candidates_read=0)
    if len(hits) > 1 and not named:
        raise CrossDockError(
            f"{len(hits)} crystallographic instances of {entry['id']} exist in this dossier "
            f"({', '.join(h['pdb_id'] + ' ' + h['comp_id'] + ' ' + h['label_asym_id'] for h in hits)}): "
            "name one with reference, because which one is used changes the RMSD.",
            candidates=hits)
    return hits[0], hits, unreadable


def cross_dock_rmsd(uniprot_id: str, compound: str, structure: str = "", pocket: str = "",
                    reference: str = "", pose: int = 0) -> dict:
    """The RMSD of a T1 compound's docked poses against its crystallographic pose, brought here by
    a protein superposition. See the module docstring; pre-registration in ESTADO 7.4."""
    from rdkit import Chem

    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    path = choice.manifest_path(accession)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"ok": False, "uniprot_id": accession,
                "error": f"No readable manifest at {path}: run prepare_structures first."}
    prep = path.parent
    try:
        if not (compound or "").strip():
            raise CrossDockError("No compound: this tool measures one compound of the recorded "
                                 "selection against its own crystallographic pose.")
        # the level is checked BEFORE the structure is resolved: "this compound has no reference
        # pose" is true of the compound, whatever structure is named, and it is the answer the
        # user needs rather than a complaint about a missing argument
        entry = dock._compound(manifest, compound)
        if (entry.get("level") or "").upper() != "T1":
            raise CrossDockError(
                f"{entry['id']} is level {entry.get('level')}, not T1. An RMSD needs a "
                "crystallographic pose of the SAME compound to compare against, and only T1 has "
                "one. For T3 the RMSD is not applicable, which is a property of the experiment, "
                "not a gap in the data.", level=entry.get("level"))

        s, snap, origin = dock._resolve(manifest, structure, pocket)
        pocket_key = f"{snap['pocket']['comp_id']} {snap['pocket']['label_asym_id']}"
        pdb_id = (s.get("pdb_id") or "").upper()
        runs = ((s.get("docking") or {}).get(pocket_key) or {}).get("runs") or {}
        run_id = dock.run_id_of(pocket_key.split()[0], "regenerated", entry["slug"])
        if run_id not in runs:
            raise CrossDockError(f"No run '{run_id}': dock {entry['id']} into {pdb_id} {pocket_key} "
                                 f"first. Available: {', '.join(sorted(runs)) or 'none'}.")
        run = runs[run_id]

        poses_rel = next((v["path"] for k, v in (run.get("files") or {}).items()
                          if k.endswith(".sdf")), None)
        if not poses_rel:
            raise CrossDockError(f"Run {run_id} records no poses SDF.")
        poses_path = (prep / poses_rel).resolve()
        recorded = (run["files"][Path(poses_rel).name] or {}).get("sha256")
        if recorded and dock.sha256(poses_path) != recorded:
            raise CrossDockError(f"{Path(poses_rel).name} changed since the manifest (SHA-256): "
                                 "dock again before measuring an RMSD against it.")
        with open(poses_path, "rb") as fh:
            poses = [m for m in Chem.ForwardSDMolSupplier(fh, removeHs=False)]
        if any(m is None for m in poses):
            raise CrossDockError(f"{Path(poses_rel).name}: {sum(m is None for m in poses)} poses "
                                 "could not be read by RDKit.")

        ref, candidates, unreadable = _reference_instance(manifest, prep, entry, pdb_id, reference)
        ref_struct = next(x for x in manifest["structures"]
                          if (x.get("pdb_id") or "").upper() == ref["pdb_id"])
        target_raw = (prep / (s.get("file") or {}).get("path")).resolve()
        mobile_raw = (prep / (ref_struct.get("file") or {}).get("path")).resolve()
        for raw, st in ((target_raw, s), (mobile_raw, ref_struct)):
            if dock.sha256(raw) != (st.get("file") or {}).get("sha256"):
                raise CrossDockError(f"{raw.name} changed since the manifest (SHA-256): run "
                                     "prepare_structures again before superposing it.")

        (sup, mob_chain, tgt_chain), tried = best_superposition(target_raw, mobile_raw)
        with open((prep / ref["sdf"]).resolve(), "rb") as fh:
            crystal = next(Chem.ForwardSDMolSupplier(fh, removeHs=False, sanitize=True), None)
        if crystal is None:
            raise CrossDockError(f"RDKit could not read the reference ligand {ref['sdf']}.")
        moved = moved_copy(crystal, sup.transform)

        # the error bar: the SAME reference pose, superposed on the pocket residues instead of the
        # whole chain. The difference is how much the "reference" depends on a choice of this tool
        bar = _shift_bar(prep, target_raw, mobile_raw, tgt_chain, mob_chain, snap, crystal, moved)

        # geometry, never a name: where the reference pose lands in the receptor that was docked
        geometry = in_box(heavy_xyz(moved), snap["pocket"]["box"])
        pocket_sdf_rel = next((ni.get("sdf", {}).get("path") for ni in s.get("nonpolymer_instances") or []
                               if ni.get("comp_id") == snap["pocket"]["comp_id"]
                               and ni.get("label_asym_id") == snap["pocket"]["label_asym_id"]
                               and (ni.get("sdf") or {}).get("ok")), None)
        if pocket_sdf_rel and (prep / pocket_sdf_rel).is_file():
            with open((prep / pocket_sdf_rel).resolve(), "rb") as fh:
                pocket_lig = next(Chem.ForwardSDMolSupplier(fh, removeHs=True, sanitize=False), None)
            if pocket_lig is not None:
                geometry["centroid_to_pocket_ligand_centroid_angstrom"] = round(float(np.linalg.norm(
                    heavy_xyz(moved).mean(0) - heavy_xyz(pocket_lig).mean(0))), 2)
        else:
            # absence of the file is not "the pose is nowhere near the pocket ligand"
            geometry["centroid_to_pocket_ligand_centroid_angstrom"] = None
            geometry["pocket_ligand_not_compared"] = (
                "the SDF of this pocket's own ligand is not in the prep folder, so the distance to "
                "it was not measured; the box check above still applies")

        control = dock.negative_control(moved)
        rmsds = _rmsds(poses, moved)
        top = rmsds[0]
        best_i = int(np.argmin(rmsds))
    except (CrossDockError, dock.DockError) as exc:
        out = {"ok": False, "uniprot_id": accession, "error": str(exc)}
        out.update(getattr(exc, "extra", {}) or {})
        return out
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        log.exception("cross_dock_rmsd failed")
        return {"ok": False, "uniprot_id": accession, "error": f"No RMSD measured: {exc!r}"}

    block = {
        "generated_by": "cross_dock_rmsd",
        "generated_at": utc_now(),
        "level": "T1",
        "compound_id": entry["id"],
        "name": entry.get("name"),
        "inchikey": entry.get("inchikey"),
        "docked_into": f"{pdb_id} {pocket_key}",
        "reference_pose": {
            "from": f"{ref['pdb_id']} {ref['comp_id']} {ref['label_asym_id']}",
            "matched_by": "InChIKey of the crystallographic SDF, not by compound name or comp_id",
            "inchikey": ref["inchikey"], "sdf": ref["sdf"], "sha256": ref["sha256"],
            "other_instances_found": [h for h in candidates if h is not ref] or None,
            "ligand_files_unreadable": unreadable or None,
        },
        "superposition": {
            "method": "gemmi.calculate_superposition, SupSelect.CaP, global (whole polymer chain)",
            "mobile": f"{ref['pdb_id']} chain {mob_chain}",
            "target": f"{pdb_id} chain {tgt_chain}",
            "rmsd_angstrom": round(float(sup.rmsd), 3),
            "pairs": int(sup.count),
            "chain_pairs_tried": tried,
            "set_fixed_in_advance": "ESTADO 7.4: the global set, not the pocket-restricted one",
        },
        "reference_pose_uncertainty": {**bar, "preregistered_bar_angstrom": PREREGISTERED_SHIFT_BAR,
                                       "note": BAR_NOTE},
        "reference_pose_geometry": geometry,
        "rmsd_angstrom_by_pose": rmsds,
        "top_pose_rmsd_angstrom": top,
        "best_pose": {"index": best_i + 1, "rmsd_angstrom": rmsds[best_i],
                      "note": "The best of the poses is reported, and it is NOT the result: the "
                              "result is the top-scored pose, because at prediction time the score "
                              "is all there is to choose with."},
        "criterion_angstrom": CRITERION,
        "within_criterion": top <= CRITERION,
        "within_criterion_with_bar": {
            "range_angstrom": [round(top - bar["shift_angstrom"], 3), round(top + bar["shift_angstrom"], 3)],
            "note": "The criterion applied to the RMSD plus and minus the shift bar. If the two "
                    "ends fall on different sides of the criterion, this run does not decide.",
            "decides": (top + bar["shift_angstrom"] <= CRITERION) or (top - bar["shift_angstrom"] > CRITERION),
        },
        "rmsd_method": METHOD,
        "rmsd_negative_control": control,
        "criterion_note": NOISE_NOTE,
        "scope": (f"One cross-docking of one compound ({entry['id']}) from one crystal "
                  f"({ref['pdb_id']}) into one other ({pdb_id} {pocket_key}). It is not the accuracy "
                  "of the method, and it says nothing about the T3 compounds of this set, which have "
                  "no reference pose at all."),
        "values": "calculated",
    }
    warn = _record(path, manifest, pdb_id, pocket_key, run_id, block)
    out = {"ok": True, "uniprot_id": accession, "structure_origin": origin,
           "run_id": run_id, "cross_docking": block, "recorded_in_manifest": warn is None}
    if warn:
        out["manifest_warning"] = warn
    return out


def _rmsds(poses, reference) -> list:
    """Every pose against the moved reference, with dock.py's own measurement."""
    from rdkit import Chem
    ref = Chem.RemoveHs(reference)
    out = []
    for m in poses:
        try:
            out.append(dock.rmsd_to_crystal(m, ref)[0])
        except Exception as exc:
            raise CrossDockError(
                "The docked pose and the crystallographic pose could not be matched atom for atom "
                f"({exc!r}), so no RMSD was measured. They are the same InChIKey, so the usual "
                "cause is a different protonation or bond order in one of the two files.")
    return out


def _shift_bar(prep: Path, target_raw: Path, mobile_raw: Path, tgt_chain: str, mob_chain: str,
               snap: dict, crystal, moved_global) -> dict:
    """How far the reference pose moves when the superposition set changes (ESTADO 7.4).

    The pocket set is derived by geometry: residues of the target chain with a heavy atom within
    POCKET_CUTOFF of the pocket's own ligand, present in both crystals. If it cannot be built the
    bar is not invented: the pre-registered figure is reported as the only one available.
    """
    from rdkit import Chem
    try:
        tgt_ca, mob_ca = _ca_by_number(target_raw, tgt_chain), _ca_by_number(mobile_raw, mob_chain)
        common = sorted(set(tgt_ca) & set(mob_ca))
        # the pocket ligand's coordinates come from the prepared SDF of the pocket instance
        lig_xyz = None
        sdf = (snap.get("pocket") or {}).get("sdf")
        rel = sdf.get("path") if isinstance(sdf, dict) else sdf
        if rel and (prep / rel).is_file():
            with open((prep / rel).resolve(), "rb") as fh:
                m = next(Chem.ForwardSDMolSupplier(fh, removeHs=True, sanitize=False), None)
            if m is not None:
                lig_xyz = heavy_xyz(m)
        if lig_xyz is None:
            # fall back on the box centre: a sphere of the box's own half-diagonal is not the
            # pocket, so the set is NOT built from it -- the bar is reported as not measured
            return {"shift_angstrom": PREREGISTERED_SHIFT_BAR, "measured_here": False,
                    "why": "the pocket ligand's SDF was not available, so the pocket-restricted "
                           "superposition could not be built; the pre-registered bar is used"}
        nums = []
        stp = _polymers(target_raw)
        for ch in stp[0]:
            if ch.name != tgt_chain:
                continue
            for res in ch:
                xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name != "H"])
                if len(xyz) and np.sqrt(((xyz[:, None, :] - lig_xyz[None, :, :]) ** 2).sum(-1)).min() <= POCKET_CUTOFF:
                    nums.append(res.seqid.num)
        nums = sorted(set(nums) & set(common))
        if len(nums) < 4:
            return {"shift_angstrom": PREREGISTERED_SHIFT_BAR, "measured_here": False,
                    "why": f"only {len(nums)} pocket residues are present in both crystals, too few "
                           "to superpose on; the pre-registered bar is used"}
        R, t, rms = kabsch(np.array([mob_ca[n] for n in nums]), np.array([tgt_ca[n] for n in nums]))
        local = (R @ heavy_xyz(crystal).T).T + t
        shift = float(np.sqrt(((local - heavy_xyz(moved_global)) ** 2).sum(1).mean()))
        return {"shift_angstrom": round(shift, 3), "measured_here": True,
                "pocket_superposition_rmsd_angstrom": round(rms, 3),
                "pocket_residues_used": len(nums),
                "why": ("the same reference pose superposed on the pocket residues instead of the "
                        "whole chain, atom-to-atom RMSD between the two placements")}
    except Exception as exc:                       # an error bar that fails is reported, not hidden
        log.warning("shift bar not measured: %s", exc)
        return {"shift_angstrom": PREREGISTERED_SHIFT_BAR, "measured_here": False,
                "why": f"the pocket-restricted superposition failed ({exc!r}); the pre-registered "
                       "bar of ESTADO 7.4 is used"}


def _record(path: Path, manifest: dict, pdb_id: str, pocket_key: str, run_id: str, block: dict):
    """Write the block inside the docking run it measures. The manifest is re-read first, as
    dock_ligand and build_md_package do."""
    try:
        current = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return f"The RMSD was measured but the manifest could not be re-read ({exc!r})."
    if current.get("created_at") != manifest.get("created_at"):
        return ("prepare_structures rewrote the manifest while this ran, so the RMSD was not "
                "recorded. Run cross_dock_rmsd again.")
    try:
        s = next(x for x in current["structures"] if (x.get("pdb_id") or "").upper() == pdb_id)
        run = s["docking"][pocket_key]["runs"][run_id]
        run["cross_docking"] = block
        # the pose list carries a null RMSD per pose by design (no reference at docking time); it
        # is filled here, with the block above saying where the reference came from
        for i, p in enumerate(run.get("poses") or []):
            if i < len(block["rmsd_angstrom_by_pose"]):
                p["rmsd_to_reference_pose_angstrom"] = block["rmsd_angstrom_by_pose"][i]
                p["rmsd_reference"] = block["reference_pose"]["from"]
        choice.write_manifest(path, current)
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        log.exception("Could not write %s", path)
        return f"The RMSD was measured but not recorded in the manifest ({exc!r})."
    return None


def compact(result: dict) -> dict:
    """Chat view: the numbers and their error bar, without the per-chain superposition table."""
    if not result.get("ok"):
        return result
    b = result["cross_docking"]
    return {"ok": True, "uniprot_id": result["uniprot_id"], "run_id": result["run_id"],
            "compound": f"{b.get('name') or b['compound_id']} ({b['compound_id']}), level T1",
            "docked_into": b["docked_into"],
            "reference_pose": b["reference_pose"]["from"],
            "superposition_rmsd_angstrom": b["superposition"]["rmsd_angstrom"],
            "superposition_pairs": b["superposition"]["pairs"],
            "top_pose_rmsd_angstrom": b["top_pose_rmsd_angstrom"],
            "shift_bar_angstrom": b["reference_pose_uncertainty"]["shift_angstrom"],
            "criterion_angstrom": b["criterion_angstrom"],
            "within_criterion": b["within_criterion"],
            "decides": b["within_criterion_with_bar"]["decides"],
            "best_pose": b["best_pose"],
            "negative_control": b["rmsd_negative_control"],
            "notes": [b["criterion_note"], b["reference_pose_uncertainty"]["note"], b["scope"]]}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m protein_mcp.files.crossdock",
                                description=__doc__.splitlines()[0])
    p.add_argument("uniprot_id")
    p.add_argument("--compound", required=True)
    p.add_argument("--structure", default="")
    p.add_argument("--pocket", default="")
    p.add_argument("--reference", default="", help="PDB id of the crystal that holds the reference "
                                                   "pose, when more than one does")
    a = p.parse_args(argv)
    res = cross_dock_rmsd(a.uniprot_id, a.compound, a.structure, a.pocket, a.reference)
    sys.stdout.buffer.write((json.dumps(res, ensure_ascii=False, indent=1) + "\n").encode("utf-8"))
    return 0 if res.get("ok") else 1


if __name__ == "__main__":
    raise SystemExit(main())
