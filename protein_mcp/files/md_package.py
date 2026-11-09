"""build_md_package + check_charmm_gui_output: hand a docked complex to CHARMM-GUI and verify
what comes back (Phase 3, step 4).

CHARMM-GUI has no public API for protein-ligand solution systems, so the delivery is manual: this
module writes the folder the user uploads, and reads the .tgz the wizard returns. Everything here
is pure computation, no network.

Measured on four real CHARMM-GUI jobs (ESTADO 8.4) and enforced here:
- the PDB must carry NO ANISOU records. CHARMM-GUI's missing-residue path reads them as ATOM
  lines (job 9418371365: 3,948 ATOM lines for 1,974 atoms, the extra ones at the ANISOU U values,
  e.g. 9855.000 9321.000 236.000), and Modeller then returns an expanded protein (Rg 18.6 ->
  22.8 A). The wizard reports no error, so only the output shows it;
- the PDB must carry the REMARK 465 block, or the wizard offers no "Model missing residues" at
  all (job 9417550068). The block is rebuilt here from _pdbx_unobs_or_zero_occ_residues of the
  mmCIF already downloaded (77 rows for 3OG7 = the 84 lines RCSB writes), never fetched again;
- the file must be named <PDB_ID>.pdb: CHARMM-GUI takes the entry id from the file name and asks
  RCSB for that sequence ("Fail to get the chain sequence from RCSB" with 3OG7_complex.pdb);
- the ligand goes in HETATM with CONECT and a SDF in the SAME atom order (job 9418739369 kept the
  33 heavy atoms to 0.0000 A and renamed only S1->S, CL1->CL plus equivalent hydrogens).

check_charmm_gui_output exists because two of those four jobs "finished" with a broken system (one
without the ligand, one with the protein expanded), so reaching the end of the wizard proves
nothing. Its thresholds come from the two measured jobs, good (9418739369) and bad (9418371365).

Pure function + CLI: python -m protein_mcp.files.md_package build P15056 --structure 3OG7
                     python -m protein_mcp.files.md_package check P15056 --archive charmm-gui.tgz
Imports common, files.choice and files.dock (its manifest resolution is reused so the package and
the docking run cannot disagree) (+ gemmi, rdkit, numpy).
"""
import argparse
import hashlib
import json
import math
import re
import shutil
import sys
import tarfile
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice, dock

LIGAND_CHAIN = "L"          # the cleaned receptor never uses it (checked before writing)
LIGAND_RESNAME = "LIG"
CALCULATED = ("Coordinates and parameters here are calculated or chosen by CHARMM-GUI, not "
              "experimental measurements.")
# the delivery folder holds no CHARMM-GUI output yet, and mixing the two provenances in one
# sentence is the project's provenance rule failing in the other direction
BUILD_NOTE = ("The receptor coordinates are experimental (the cleaned mmCIF, unchanged); the "
              "ligand pose is calculated by AutoDock Vina. Nothing here comes from CHARMM-GUI yet.")
# a heavy atom of the protein this close to a heavy atom of the ligand is a pocket contact
CONTACT_CUTOFF = 4.5
MIN_POCKET_OVERLAP = 0.5    # fraction of the delivered contact residues that must still touch it

# ---- thresholds of check_charmm_gui_output -------------------------------------------------
# Every one is set between the good job (9418739369) and the broken one (9418371365), both
# measured; the numbers of each are in the comment so a future change can see the margin.
RG_TOLERANCE = 0.5          # A. good: 0.02 / 0.06; broken: 4.15 / 4.93
MEAN_DD_TOLERANCE = 0.5     # A, mean |d(i,j) change| over CA pairs. good: 0.128 / 0.119; broken: 6.03 / 6.75
LIGAND_DD_TOLERANCE = 0.05  # A, same measure over the ligand's heavy atoms. good: 0.0000
MAX_PEPTIDE_BOND = 2.0      # A, C(i)-N(i+1). An unfilled gap shows up as 8-24 A
MIN_HEAVY_CONTACT = 2.0     # A, protein-ligand. good: 2.79 (LYS483:CE-C5)
# Distances are compared WITHOUT superposition (CA-CA distance matrix and radius of gyration):
# CHARMM-GUI reorients the system, and a least-squares fit on a distorted structure hides the
# distortion (the first measurement of job 9418371365 did exactly that).

_REMARK465_HEADER = (
    "REMARK 465",
    "REMARK 465 MISSING RESIDUES",
    "REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE",
    "REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
    "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
    "REMARK 465",
    "REMARK 465   M RES C SSSEQI",
)
# files of the .tgz that check_charmm_gui_output reads (suffix match inside the job folder)
_WANTED = ("step1_pdbreader.pdb", "step1_pdbreader.inp", "step1_pdbreader.out",
           "step3_pbcsetup.str", "step2.2_ions.str", "lig/lig.rtf", "lig/LIG.sdf", "lig/LIG.mol2",
           "namd/step3_input.pdb", "namd/step5_production.inp", "namd/sysinfo.dat")
_MAX_MEMBER_BYTES = 200 * 1024 * 1024


class MdPackageError(Exception):
    def __init__(self, message: str, **extra):
        super().__init__(message)
        self.extra = extra


# ------------------------------------------------------------------ pure pieces

def remark465_lines(cif_path: Path) -> tuple:
    """(lines, source) for the REMARK 465 block of an mmCIF.

    Rebuilt from _pdbx_unobs_or_zero_occ_residues, byte-identical to the block RCSB writes in its
    own PDB file (checked on 3OG7, 5JSM, 8C7Y, 2W2N, 3WD5, 6M0J and 1A2C; empty where RCSB's is).

    Only `occupancy_flag 1` rows are kept. That category also holds the residues that WERE
    modelled with zero occupancy, which RCSB writes in REMARK 475, not 465: they have
    coordinates, they survive cleaning, and calling them missing would tell the analyst to model
    a residue that is already there (1M17 GLU A 725, 1IGY and 12E8 entirely).

    `source` is "mmCIF" when the category was read, or a reason string when it was not, so the
    caller never reports "no missing residues" for an entry it could not read.
    """
    from gemmi import cif
    try:
        block = cif.read(str(cif_path)).sole_block()
    except Exception as exc:                                  # gemmi raises plain RuntimeError
        raise MdPackageError(f"Could not read {cif_path.name} as mmCIF: {exc!r}")
    tags = ["PDB_model_num", "polymer_flag", "auth_asym_id", "auth_comp_id", "auth_seq_id",
            "occupancy_flag", "PDB_ins_code"]
    table = block.find("_pdbx_unobs_or_zero_occ_residues.", tags)
    if not len(table):
        # absent category and "every residue observed" are different facts and must not be mixed:
        # block.find also returns nothing when one of the tags is missing from the file
        present = [t for t in tags if block.find_loop(f"_pdbx_unobs_or_zero_occ_residues.{t}")]
        if not present:
            return [], "no _pdbx_unobs_or_zero_occ_residues category in the mmCIF"
        missing = [t for t in tags if t not in present]
        if missing:
            return [], f"_pdbx_unobs_or_zero_occ_residues without {', '.join(missing)}"
    rows, zero_occupancy = [], 0
    for row in table:
        if cif.as_string(row[1]) != "Y":                      # polymer_flag: only polymer residues
            continue
        if cif.as_string(row[5]) != "1":                      # occupancy_flag: 0 = REMARK 475
            zero_occupancy += 1
            continue
        model = cif.as_string(row[0])
        model = "" if model in ("1", ".", "?", "") else model
        ins = cif.as_string(row[6])
        ins = "" if ins in (".", "?") else ins
        rows.append((model, cif.as_string(row[3]), cif.as_string(row[2]), cif.as_string(row[4]), ins))
    if zero_occupancy:
        # those residues are IN the file, with coordinates: saying nothing about them would turn
        # "no residue is missing" into "every residue is experimental", which is false
        note = (f"mmCIF; {zero_occupancy} residue(s) are in the file with zero occupancy "
                "(REMARK 475), modelled rather than located by the experiment")
    else:
        note = "mmCIF"
    if not rows:
        return [], note
    out = [line.ljust(80) for line in _REMARK465_HEADER]
    for model, resname, chain, seq, ins in rows:
        cols = [" "] * 80
        cols[0:10] = list("REMARK 465")
        if model:
            cols[14] = model[0]                               # column 15: model number
        cols[15:18] = list(f"{resname:>3}")                   # columns 16-18: residue name
        cols[19] = chain[0] if chain else " "                 # column 20: chain
        cols[21:26] = list(f"{seq:>5}")                       # columns 22-26: sequence number
        if ins:
            cols[26] = ins[0]                                 # column 27: insertion code
        out.append("".join(cols))
    return out, note


def missing_ranges(lines: list) -> list:
    """Contiguous runs of missing residues as {chain, first, last, length} from a REMARK 465 block."""
    runs = []
    for line in lines[len(_REMARK465_HEADER):]:
        chain, seq = line[19], line[21:26].strip()
        if not seq.lstrip("-").isdigit():
            continue
        num = int(seq)
        if runs and runs[-1]["chain"] == chain and num == runs[-1]["last"] + 1:
            runs[-1]["last"] = num
        else:
            runs.append({"chain": chain, "first": num, "last": num})
    for run in runs:
        run["length"] = run["last"] - run["first"] + 1
    return runs


def classify_ranges(runs: list, observed: dict, present=frozenset()) -> list:
    """Mark each run terminal or internal and say whether it is really absent from the receptor.

    `observed` maps chain -> (first, last) residue present; `present` is the set of
    (chain, residue number) the delivered receptor actually contains.

    CHARMM-GUI leaves the N-terminal runs unticked by default ("This is N-terminal missing
    residues"); only the internal ones are the gaps that would otherwise be bridged. A run is
    only offered for modelling when EVERY residue in it is absent from the file that is uploaded:
    the mmCIF says what the experiment did not locate, the uploaded file is what CHARMM-GUI reads,
    and the two can disagree.
    """
    out = []
    for run in runs:
        span = observed.get(run["chain"])
        if not span:
            kind = "chain not in the receptor"
        elif run["last"] < span[0]:
            kind = "N-terminal"
        elif run["first"] > span[1]:
            kind = "C-terminal"
        else:
            kind = "internal"
        still_there = [n for n in range(run["first"], run["last"] + 1) if (run["chain"], n) in present]
        entry = {**run, "kind": kind, "absent_from_the_receptor": not still_there}
        if still_there:
            entry["in_the_receptor"] = still_there
        out.append(entry)
    return out


def observed_spans(atom_lines: list) -> dict:
    """chain -> (first, last) residue number over ATOM records."""
    spans = {}
    for line in atom_lines:
        if not line.startswith("ATOM"):
            continue
        chain, num = line[21], line[22:26].strip()
        if not num.lstrip("-").isdigit():
            continue
        num = int(num)
        low, high = spans.get(chain, (num, num))
        spans[chain] = (min(low, num), max(high, num))
    return spans


def ligand_atom_names(mol) -> list:
    """Element + running index (C1, C2, ..., H1): unique, at most 4 characters, PDB-safe."""
    count, names = Counter(), []
    for atom in mol.GetAtoms():
        element = atom.GetSymbol().upper()
        count[element] += 1
        name = f"{element}{count[element]}"
        if len(name) > 4:
            raise MdPackageError(f"Atom name {name} does not fit the PDB column (4 characters).")
        names.append(name)
    if len(set(names)) != len(names):
        raise MdPackageError("Atom names are not unique.")
    return names


def complex_pdb_text(receptor_text: str, pose, names: list, remark: list) -> str:
    """Cleaned receptor + pose as HETATM LIG chain L, with CONECT, without ANISOU.

    The ligand's first serial continues after the receptor's highest one, TER records included
    (the prototype first reused the serial of chain B's TER).
    """
    keep = [line.rstrip("\n") for line in receptor_text.splitlines()
            if not line.startswith(("END", "CONECT", "MASTER", "ANISOU"))]
    if any(line.startswith(("ATOM", "HETATM")) and line[21] == LIGAND_CHAIN for line in keep):
        raise MdPackageError(f"The cleaned receptor already uses chain {LIGAND_CHAIN}; "
                             "the ligand would collide with it.")
    if remark:
        last = max((i for i, line in enumerate(keep) if line.startswith("REMARK")), default=-1)
        if last < 0:
            last = max((i for i, line in enumerate(keep)
                        if line.startswith(("HEADER", "TITLE", "COMPND", "EXPDTA", "KEYWDS"))), default=-1)
        keep = keep[:last + 1] + list(remark) + keep[last + 1:]
    serials = [int(line[6:11]) for line in keep
               if line.startswith(("ATOM", "HETATM", "TER")) and line[6:11].strip()]
    if not serials:
        raise MdPackageError("The cleaned receptor has no atom.")
    start = max(serials) + 1

    conf, hetatm, conect = pose.GetConformer(), [], []
    for i, atom in enumerate(pose.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        name = names[i]
        # PDB convention: one-letter element symbols start in column 14
        field = f" {name:<3}" if len(name) < 4 else name
        hetatm.append(f"HETATM{start + i:5d} {field:<4} {LIGAND_RESNAME} {LIGAND_CHAIN}   1    "
                      f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00          {atom.GetSymbol():>2}")
    for atom in pose.GetAtoms():
        neighbours = sorted(n.GetIdx() for n in atom.GetNeighbors())
        for k in range(0, len(neighbours), 4):
            conect.append("CONECT" + f"{start + atom.GetIdx():5d}"
                          + "".join(f"{start + j:5d}" for j in neighbours[k:k + 4]))
    ter = f"TER   {start + pose.GetNumAtoms():5d}      {LIGAND_RESNAME} {LIGAND_CHAIN}   1"
    return "\n".join(keep + hetatm + [ter] + conect + ["END"]) + "\n"


def ligand_sdf_text(pose, name: str = LIGAND_RESNAME) -> str:
    """The pose as a one-molecule SDF in the same atom order as the PDB, without the run's props."""
    from rdkit import Chem
    mol = Chem.Mol(pose)
    for prop in list(mol.GetPropNames()):
        mol.ClearProp(prop)
    mol.SetProp("_Name", name)
    return Chem.MolToMolBlock(mol) + "$$$$\n"


def readme_text(pdb_id: str, ranges: list, internal: list, ligand: dict, pose_line: str,
                remark_lines: int, source: str) -> str:
    terminal = [r for r in ranges if r["kind"] != "internal"]
    partial = [r for r in ranges if r["kind"] == "internal" and r not in internal]
    fmt = lambda runs: ", ".join(f"{r['chain']} {r['first']}-{r['last']}" for r in runs) or "none"
    if remark_lines:
        block = f"REMARK 465 ({remark_lines} lines, from {source})"
    else:
        block = f"**no REMARK 465 block** ({source}), so CHARMM-GUI will offer no loop modelling"
    step3 = (f'3. "Model missing residues": tick the internal runs ({fmt(internal)}).\n'
             f"   Leave the terminal runs unticked ({fmt(terminal)}): they are disordered tails, not gaps."
             if internal or terminal else
             '3. "Model missing residues": this entry has no run to tick.')
    if partial:
        # neither a gap to model nor a disordered tail: part of the run is already in the file
        step3 += (f"\n   NOT offered, and not tails: {fmt(partial)}. Some of their residues are "
                  "already in the uploaded file, so CHARMM-GUI cannot model the run. Whatever is "
                  "still absent will be bridged by a long bond, and check_charmm_gui_output will "
                  "say so. Look at the entry at RCSB before you simulate this system.")
    return f"""# CHARMM-GUI delivery: {pdb_id} with {ligand.get('comp_id') or LIGAND_RESNAME}

Generated by build_md_package ({utc_now()}). Solution Builder, one job.
{BUILD_NOTE}

## Files
- `{pdb_id}.pdb` - cleaned receptor + the docked pose as `{LIGAND_RESNAME}`, chain `{LIGAND_CHAIN}`,
  with CONECT, **no ANISOU**, and {block}. Keep the file name: CHARMM-GUI reads the entry id
  from it and asks RCSB for that sequence.
- `{pdb_id}_ligand.sdf` - the same pose, same atom order. {pose_line}

## Steps
1. Solution Builder -> upload `{pdb_id}.pdb` (PDB format).
2. "Reading Hetero Chain Residues": tick `{LIGAND_RESNAME}`, choose "Use CHARMM General Force
   Field ... (using ParamChem service)" and "the SDF file uploaded from" -> `{pdb_id}_ligand.sdf`.
{step3}
4. Terminal group patching, waterbox, ions and the rest: defaults.
5. Pick NAMD in the Input Generator and download the .tgz.
6. Check it before using it:
   `python -m protein_mcp.files.md_package check <UNIPROT> --archive <the .tgz>`

## What the simulation will contain that the crystal did not
The residues of the runs ticked in step 3 ({fmt(internal)}) are **modelled, not experimental**.
Any report must say so wherever they are discussed.
"""


# ------------------------------------------------------------------ geometry helpers

def read_atoms(path: Path, keep=None) -> list:
    """ATOM/HETATM records as dicts, by fixed PDB columns."""
    out = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            # a coordinate that does not fit its columns: not an atom this code can use. It
            # happens in the files CHARMM-GUI builds from ANISOU records read as ATOM lines
            # (job 9418371365), which malformed_atom_lines counts as the symptom it is.
            continue
        atom = {"name": line[12:16].strip(), "resname": line[17:21].strip(), "chain": line[21],
                "resid": line[22:26].strip(), "xyz": xyz,
                "segid": line[72:76].strip() if len(line) >= 76 else "",
                "element": (line[76:78].strip() if len(line) >= 78 else "") or line[12:16].strip()[:1]}
        if not atom["resid"].lstrip("-").isdigit():
            continue
        atom["resid"] = int(atom["resid"])
        if keep is None or keep(atom):
            out.append(atom)
    return out


def malformed_atom_lines(path: Path) -> int:
    """ATOM/HETATM lines whose coordinates do not parse or lie beyond 1000 A.

    This is the fingerprint of CHARMM-GUI reading ANISOU records as atoms: the U values land in
    the coordinate columns, giving either 9855.000 or a shifted field like '0 648.00'.
    """
    bad = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        try:
            if max(abs(float(line[30:38])), abs(float(line[38:46])), abs(float(line[46:54]))) > 1000:
                bad += 1
        except ValueError:
            bad += 1
    return bad


def radius_of_gyration(coords) -> float:
    xyz = np.asarray(coords, dtype=float)
    return float(np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(1).mean()))


def internal_distance_change(before, after) -> dict:
    """Mean and max |distance change| over every pair. No superposition, so a rotation, a
    translation or a re-centring of the whole system does not show up and a distortion does."""
    before, after = np.asarray(before, dtype=float), np.asarray(after, dtype=float)
    if len(before) != len(after) or len(before) < 2:
        raise MdPackageError(f"Cannot compare {len(before)} and {len(after)} positions.")
    d_before = np.linalg.norm(before[:, None, :] - before[None, :, :], axis=2)
    d_after = np.linalg.norm(after[:, None, :] - after[None, :, :], axis=2)
    diff = np.abs(d_before - d_after)[np.triu_indices(len(before), 1)]
    return {"pairs": int(diff.size), "mean_angstrom": round(float(diff.mean()), 4),
            "max_angstrom": round(float(diff.max()), 4)}


def min_distance(a, b) -> float:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2).min())


def long_peptide_bonds(atoms: list, limit: float = MAX_PEPTIDE_BOND) -> list:
    """C(i)-N(i+1) distances above `limit`: an unfilled gap bridged by an artificial bond."""
    by_residue = {}
    for atom in atoms:
        by_residue.setdefault((atom["segid"] or atom["chain"], atom["resid"]), {})[atom["name"]] = atom["xyz"]
    out = []
    for segment in sorted({key[0] for key in by_residue}):
        residues = sorted(num for seg, num in by_residue if seg == segment)
        for first, second in zip(residues, residues[1:]):
            tail, head = by_residue[(segment, first)].get("C"), by_residue[(segment, second)].get("N")
            if not tail or not head:
                continue
            distance = math.dist(tail, head)
            if distance > limit:
                out.append({"segment": segment, "from": first, "to": second,
                            "distance_angstrom": round(distance, 2)})
    return out


# ------------------------------------------------------------------ build_md_package

def contact_residues(protein: list, ligand: list, cutoff: float = CONTACT_CUTOFF) -> list:
    """Sorted "<chain or segment>:<number>" of the protein residues touching the ligand.

    This is what says the ligand is still IN its pocket. Comparing the shape of the pose is not
    enough: a pose translated 40 A away, or mirrored, keeps every internal distance.
    """
    heavy_p = [a for a in protein if not a["element"].startswith("H")]
    heavy_l = [a for a in ligand if not a["element"].startswith("H") and not a["name"].startswith("LP")]
    if not heavy_p or not heavy_l:
        return []
    coords_p = np.array([a["xyz"] for a in heavy_p], dtype=float)
    coords_l = np.array([a["xyz"] for a in heavy_l], dtype=float)
    close = (np.linalg.norm(coords_p[:, None, :] - coords_l[None, :, :], axis=2) < cutoff).any(axis=1)
    out = {f"{heavy_p[i]['segid'] or heavy_p[i]['chain']}:{heavy_p[i]['resid']}"
           for i in range(len(heavy_p)) if close[i]}
    return sorted(out)


def _contact_residues(pdb_path: Path) -> list:
    atoms = read_atoms(pdb_path)
    return contact_residues([a for a in atoms if a["resname"] != LIGAND_RESNAME],
                            [a for a in atoms if a["resname"] == LIGAND_RESNAME])


def _record_in_manifest(path: Path, manifest: dict, pdb_id: str, pocket_key: str, package: dict,
                        out_dir: Path, prep: Path):
    """Write an `md_package` block next to the docking block (decision of 2026-11-06).

    The manifest is re-read first: if prepare_structures rebuilt it meanwhile, the files stay and
    this says so, exactly as dock_ligand does.
    """
    try:
        current = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return f"The package is in {out_dir} but the manifest could not be re-read ({exc!r})."
    if current.get("created_at") != manifest.get("created_at"):
        return (f"The package is in {out_dir} but prepare_structures rewrote the manifest while it "
                "was being written, so it was not recorded there. Run build_md_package again.")
    try:
        s = next(x for x in current["structures"] if x.get("pdb_id") == pdb_id)
        block = (s.setdefault("md_package", {}))
        block[pocket_key] = {k: package[k] for k in ("generated_by", "generated_at", "run_id", "pose",
                                                     "missing_residues_source", "to_model", "files",
                                                     "checks", "note")}
        block[pocket_key]["dir"] = out_dir.relative_to(prep).as_posix()
        choice.write_manifest(path, current)
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        log.exception("Could not write %s", path)
        return f"The package is in {out_dir} but the manifest could not be written ({exc!r})."
    return None


def _pick_run(s: dict, pocket_key: str, conformer: str):
    block = (s.get("docking") or {}).get(pocket_key) or {}
    runs = block.get("runs") or {}
    if not runs:
        raise MdPackageError(f"No docking run for {s.get('pdb_id')} pocket {pocket_key}: run dock_ligand first.")
    if conformer:
        run_id = dock.run_id_of(pocket_key.split()[0], conformer.strip().lower())
        if run_id not in runs:
            raise MdPackageError(f"No run '{run_id}'. Available: {', '.join(sorted(runs))}.")
    elif len(runs) == 1:
        run_id = next(iter(runs))
    else:
        raise MdPackageError(f"{len(runs)} runs for this pocket; name one with conformer "
                             f"({', '.join(sorted(r.rsplit('_', 1)[-1] for r in runs))}).")
    return run_id, runs[run_id], block


def _inside(prep: Path, rel) -> Path:
    path = (prep / str(rel)).resolve()
    if prep.resolve() not in path.parents or not path.is_file():
        raise MdPackageError(f"{rel} is missing or outside prep/.")
    return path


def build_md_package(uniprot_id: str, structure: str = "", pocket: str = "", conformer: str = "",
                     pose: int = 1) -> dict:
    """Write the folder the user uploads to CHARMM-GUI. See the module docstring for the rules."""
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
    try:
        s, snap, origin = dock._resolve(manifest, structure, pocket)
        pocket_key = f"{snap['pocket']['comp_id']} {snap['pocket']['label_asym_id']}"
        run_id, run, block = _pick_run(s, pocket_key, conformer)

        poses_rel = next((v["path"] for k, v in (run.get("files") or {}).items() if k.endswith(".sdf")), None)
        if not poses_rel:
            raise MdPackageError(f"Run {run_id} records no poses SDF.")
        prep = path.parent
        poses_path = _inside(prep, poses_rel)
        recorded = (run["files"][Path(poses_rel).name] or {}).get("sha256")
        if recorded and dock.sha256(poses_path) != recorded:
            raise MdPackageError(f"{Path(poses_rel).name} changed since the manifest (SHA-256): dock again.")
        with open(poses_path, "rb") as fh:
            poses = [m for m in Chem.ForwardSDMolSupplier(fh, removeHs=False)]
        if any(m is None for m in poses):
            raise MdPackageError(f"{Path(poses_rel).name}: {sum(m is None for m in poses)} poses unreadable.")
        if len(poses) != len(run.get("poses") or []):
            raise MdPackageError(f"{Path(poses_rel).name} has {len(poses)} poses, the manifest records "
                                 f"{len(run.get('poses') or [])}.")
        if not 1 <= pose <= len(poses):
            raise MdPackageError(f"pose {pose} is outside 1..{len(poses)}.")
        chosen, record = Chem.Mol(poses[pose - 1]), run["poses"][pose - 1]

        receptor_path = _inside(prep, snap["receptor"])
        receptor_text = receptor_path.read_text(encoding="utf-8", errors="replace")
        raw_path = _inside(prep, s["file"]["path"])
        if dock.sha256(raw_path) != (s.get("file") or {}).get("sha256"):
            raise MdPackageError(f"{raw_path.name} changed since the manifest (SHA-256): run prepare_structures again.")
        if raw_path.suffix.lower() not in (".cif", ".mmcif"):
            raise MdPackageError(f"The raw file is {raw_path.name}; REMARK 465 is rebuilt from an mmCIF.")

        remark, remark_source = remark465_lines(raw_path)
        receptor_lines = receptor_text.splitlines()
        in_receptor = {(line[21], int(line[22:26]))
                       for line in receptor_lines
                       if line.startswith("ATOM") and line[22:26].strip().lstrip("-").isdigit()}
        ranges = classify_ranges(missing_ranges(remark), observed_spans(receptor_lines), in_receptor)
        to_model = [r for r in ranges if r["kind"] == "internal" and r["absent_from_the_receptor"]]
        names = ligand_atom_names(chosen)
        pdb_id = s["pdb_id"].upper()
        complex_text = complex_pdb_text(receptor_text, chosen, names, remark)
        sdf_text = ligand_sdf_text(chosen)
        pose_line = (f"Vina mode {record['mode']}, {record['vina_score_kcal_mol']:.2f} kcal/mol (calculated), "
                     f"RMSD to the crystal pose {record['rmsd_to_crystal_pose_angstrom']:.2f} A.")

        # written to a temporary folder and verified there: a package that fails verification must
        # not be left half-written next to the package.json of a previous one
        staging = Path(tempfile.mkdtemp(prefix="md_package_"))
        try:
            pdb_path, sdf_path = staging / f"{pdb_id}.pdb", staging / f"{pdb_id}_ligand.sdf"
            pdb_path.write_text(complex_text, encoding="ascii")
            sdf_path.write_text(sdf_text, encoding="ascii")
            checks = _verify_package(pdb_path, sdf_path, names, chosen)
            (staging / "README.md").write_text(
                readme_text(pdb_id, ranges, to_model, {"comp_id": snap["pocket"]["comp_id"]},
                            pose_line, len(remark), remark_source),
                encoding="utf-8")
            out_dir = prep / "md" / dock.dock_dir_name(pdb_id, snap["pocket"]["comp_id"],
                                                       snap["pocket"]["label_asym_id"])
            out_dir.mkdir(parents=True, exist_ok=True)
            written = []
            for name in (f"{pdb_id}.pdb", f"{pdb_id}_ligand.sdf", "README.md"):
                shutil.copyfile(staging / name, out_dir / name)
                written.append(out_dir / name)
        finally:
            shutil.rmtree(staging, ignore_errors=True)

        package = {
            "generated_by": "build_md_package", "generated_at": utc_now(), "uniprot_id": accession,
            "pdb_id": pdb_id, "pocket": pocket_key, "run_id": run_id, "selected_by": origin,
            "pose": {"index": pose, **record},
            "receptor": {"path": snap["receptor"], "chains_auth": snap.get("receptor_chains_auth"),
                         "coordinates": "experimental (the cleaned receptor of the mmCIF), unchanged"},
            "ligand": {"comp_id": snap["pocket"]["comp_id"], "name": snap["pocket"].get("name"),
                       "atom_names": names, "heavy_atoms": chosen.GetNumHeavyAtoms(),
                       "atoms": chosen.GetNumAtoms(),
                       "coordinates": [[round(c, 3) for c in chosen.GetConformer().GetAtomPosition(i)]
                                       for i in range(chosen.GetNumAtoms())]},
            "missing_residues_source": remark_source,
            "missing_residue_runs": ranges,
            "to_model": to_model,
            "pocket_residues_delivered": _contact_residues(out_dir / f"{pdb_id}.pdb"),
            "files": {p.name: {"path": p.relative_to(prep).as_posix(), "sha256": dock.sha256(p)}
                      for p in written},
            "checks": checks,
            "note": BUILD_NOTE,
        }
        (out_dir / "package.json").write_text(json.dumps(package, ensure_ascii=False, indent=1), encoding="utf-8")
        manifest_note = _record_in_manifest(path, manifest, pdb_id, pocket_key, package, out_dir, prep)
    except (MdPackageError, dock.DockError) as exc:
        return {"ok": False, "uniprot_id": accession, "error": str(exc), **getattr(exc, "extra", {})}
    except (OSError, ValueError, KeyError) as exc:
        log.exception("build_md_package failed")
        return {"ok": False, "uniprot_id": accession, "error": f"Package not written: {exc!r}"}

    label = lambda runs: [f"{r['chain']} {r['first']}-{r['last']}" for r in runs]
    fmt = label(to_model)
    terminal = label([r for r in ranges if r["kind"] != "internal"])
    partial = label([r for r in ranges if r["kind"] == "internal" and r not in to_model])
    if remark_source == "mmCIF":
        no_gap = "No internal gap: every residue of the system is experimental."
    elif remark_source.startswith("mmCIF;"):
        no_gap = ("No residue is unobserved, but " + remark_source.split("; ", 1)[1]
                  + ": those coordinates are not experimental and CHARMM-GUI will not touch them.")
    else:                      # absence of data is not the fact that there is nothing to model
        no_gap = (f"No record: {remark_source}, so this tool cannot say which residues are "
                  "missing. Check the entry at RCSB before trusting the system.")
    notes = [BUILD_NOTE,
             (f"The PDB has no ANISOU record and a REMARK 465 block of {len(remark)} lines, which "
              "is what makes CHARMM-GUI keep the protein undistorted and offer loop modelling "
              "(ESTADO 8.4)." if remark else
              f"The PDB has no ANISOU record, but also NO REMARK 465 block ({remark_source}), so "
              "CHARMM-GUI will not offer loop modelling and will bridge any gap with a long bond."),
             f"Keep the file name {pdb_id}.pdb: CHARMM-GUI takes the entry id from it."]
    if partial:
        notes.insert(0, f"Runs {', '.join(partial)} are internal gaps of the mmCIF whose residues "
                        "are partly in the uploaded file, so they are neither offered for "
                        "modelling nor disordered tails. Check them at RCSB.")
    if manifest_note:
        notes.insert(0, manifest_note)
    return {
        "ok": True, "uniprot_id": accession, "package_dir": str(out_dir),
        "built": f"{pdb_id} pocket {pocket_key}, {run_id}, pose {pose} of {len(poses)}",
        "upload": [f"{pdb_id}.pdb", f"{pdb_id}_ligand.sdf"],
        "recorded_in_manifest": manifest_note is None,
        "missing_residues_source": remark_source,
        "tick_model_missing_residues": fmt or no_gap,
        "leave_unticked": terminal or "none",
        **({"internal_runs_partly_in_the_file": partial} if partial else {}),
        "modelled_not_experimental": (
            f"After CHARMM-GUI fills them, {', '.join(fmt)} are modelled, not experimental; "
            "say so wherever the report discusses them." if fmt else no_gap),
        "checks": checks,
        "next": ("Follow README.md, then: python -m protein_mcp.files.md_package check "
                 f"{accession} --archive <the .tgz CHARMM-GUI returns>"),
        "notes": notes,
    }


def _verify_package(pdb_path: Path, sdf_path: Path, names: list, pose) -> dict:
    """Re-read what was just written; MdPackageError if it is not what CHARMM-GUI needs."""
    from rdkit import Chem

    text = pdb_path.read_text(encoding="ascii")
    # counted, not assumed: a value reported as a measurement must be one
    anisou = sum(1 for line in text.splitlines() if line.startswith("ANISOU"))
    if anisou:
        raise MdPackageError(f"The written PDB still has {anisou} ANISOU records.")
    ligand = [a for a in read_atoms(pdb_path) if a["resname"] == LIGAND_RESNAME]
    with open(sdf_path, "rb") as fh:
        back = next(Chem.ForwardSDMolSupplier(fh, removeHs=False), None)
    if back is None:
        raise MdPackageError("The written SDF does not read back.")
    if len(ligand) != back.GetNumAtoms() or back.GetNumAtoms() != pose.GetNumAtoms():
        raise MdPackageError(f"Atom counts differ: PDB {len(ligand)}, SDF {back.GetNumAtoms()}, "
                             f"pose {pose.GetNumAtoms()}.")
    conf = back.GetConformer()
    agree = sum(1 for i, atom in enumerate(ligand)
                if atom["name"] == names[i]
                and math.dist(atom["xyz"], tuple(conf.GetAtomPosition(i))) <= 1e-3)
    if agree != len(ligand):
        bad = next(i for i, atom in enumerate(ligand)
                   if atom["name"] != names[i]
                   or math.dist(atom["xyz"], tuple(conf.GetAtomPosition(i))) > 1e-3)
        raise MdPackageError(f"PDB and SDF disagree at atom {bad + 1} ({ligand[bad]['name']} vs "
                             f"{names[bad]}): CHARMM-GUI needs the same order and the same coordinates.")
    receptor = [a for a in read_atoms(pdb_path) if a["resname"] != LIGAND_RESNAME]
    if not receptor:
        raise MdPackageError("The written PDB has no receptor atom.")
    contact = min_distance([a["xyz"] for a in receptor if not a["element"].startswith("H")],
                           [a["xyz"] for a in ligand if not a["element"].startswith("H")])
    return {
        "anisou_records": anisou,
        "remark_465_lines": sum(1 for line in text.splitlines() if line.startswith("REMARK 465")),
        "ligand_atoms": len(ligand),
        "ligand_heavy_atoms": back.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge(back),
        "pdb_and_sdf_atoms_in_agreement": f"{agree} of {len(ligand)}",
        "min_receptor_ligand_heavy_distance_angstrom": round(contact, 2),
        "conect_records": sum(1 for line in text.splitlines() if line.startswith("CONECT")),
    }


# ------------------------------------------------------------------ check_charmm_gui_output

def _extract(archive: Path, dest: Path) -> Path:
    """Extract only the files this check reads, with tarfile's data filter. Returns the job folder."""
    try:
        with tarfile.open(archive, "r:gz") as tar:
            members, roots = [], set()
            for member in tar.getmembers():
                name = member.name.replace("\\", "/")
                parts = name.split("/")
                if len(parts) > 1:
                    roots.add(parts[0])
                if not member.isfile() or member.size > _MAX_MEMBER_BYTES:
                    continue
                if any(name.endswith("/" + w) or name == w for w in _WANTED) or name.endswith(".pdb"):
                    members.append(member)
            if len(roots) != 1:
                raise MdPackageError(f"The archive has {len(roots)} top folders; one CHARMM-GUI job expected.")
            tar.extractall(dest, members=members, filter="data")
    except tarfile.TarError as exc:
        raise MdPackageError(f"Could not read {archive.name} as a .tar.gz: {exc!r}")
    return dest / roots.pop()


def _parse_value(path: Path, pattern: str, cast=str):
    if not path.is_file():
        return None
    match = re.search(pattern, path.read_text(encoding="utf-8", errors="replace"), re.MULTILINE)
    if not match:
        return None
    try:
        return cast(match.group(1))
    except ValueError:
        return None


def _system_values(job: Path) -> dict:
    namd = job / "namd"
    production = namd / "step5_production.inp"
    ions = [a for a in read_atoms(namd / "step3_input.pdb") if a["segid"] == "IONS"] \
        if (namd / "step3_input.pdb").is_file() else []
    waters = sum(1 for a in read_atoms(namd / "step3_input.pdb") if a["segid"] == "SOLV" and a["name"] == "OH2") \
        if (namd / "step3_input.pdb").is_file() else 0
    return {
        "box_side_angstrom": _parse_value(job / "step3_pbcsetup.str", r"^\s*SET A = ([\d.]+)", float),
        "ions": dict(Counter(a["resname"] for a in ions)) or None,
        "waters": waters or None,
        "temperature_kelvin": _parse_value(production, r"^set temp\s+([\d.]+)", float),
        "timestep_fs": _parse_value(production, r"^timestep\s+([\d.]+)", float),
        "cutoff_angstrom": _parse_value(production, r"^cutoff\s+([\d.]+)", float),
        "production_steps": _parse_value(production, r"^numsteps\s+(\d+)", int),
        "lone_pairs": (json.loads((namd / "sysinfo.dat").read_text(encoding="utf-8")).get("lonepair")
                       if (namd / "sysinfo.dat").is_file() else None),
        "cgenff_penalties": _parse_value(
            job / "lig" / "lig.rtf",
            r"^RESI\s+\S+\s+[-\d.]+\s*!\s*param penalty=\s*([\d.]+\s*;\s*charge penalty=\s*[\d.]+)",
            lambda v: re.sub(r"\s+", " ", v.strip())),
    }


def check_charmm_gui_output(uniprot_id: str, archive: str, structure: str = "", pocket: str = "") -> dict:
    """Verify the .tgz CHARMM-GUI returns against the package build_md_package handed over.

    No `conformer`: the package folder is one per structure and pocket, so a second
    build_md_package for another run replaces it and there is nothing to choose between.
    `package.json` records which run the package came from.
    """
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    archive_path = Path(archive).expanduser()
    if not archive_path.is_file():
        return {"ok": False, "uniprot_id": accession, "error": f"No file at {archive_path}."}
    path = choice.manifest_path(accession)
    prep = path.parent
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        s, snap, _ = dock._resolve(manifest, structure, pocket)
        pdb_id = s["pdb_id"].upper()
        pkg_dir = prep / "md" / dock.dock_dir_name(pdb_id, snap["pocket"]["comp_id"], snap["pocket"]["label_asym_id"])
        pkg_file = pkg_dir / "package.json"
        if not pkg_file.is_file():
            raise MdPackageError(f"No package at {pkg_file}: run build_md_package first.")
        package = json.loads(pkg_file.read_text(encoding="utf-8"))
        delivered = _inside(prep, package["files"][f"{pdb_id}.pdb"]["path"])
        with tempfile.TemporaryDirectory(prefix="cgui_") as tmp:
            result = _check_job(_extract(archive_path, Path(tmp)), delivered, package)
    except (MdPackageError, dock.DockError) as exc:
        return {"ok": False, "uniprot_id": accession, "error": str(exc), **getattr(exc, "extra", {})}
    except (OSError, ValueError, KeyError) as exc:
        log.exception("check_charmm_gui_output failed")
        return {"ok": False, "uniprot_id": accession, "error": f"Archive not checked: {exc!r}"}
    return {"ok": True, "uniprot_id": accession, "archive": str(archive_path),
            "package_dir": str(pkg_dir), **result}


def _check_job(job: Path, delivered: Path, package: dict) -> dict:
    """Every measurement of one job folder, plus the verdict. Pure, so the tests can call it."""
    failures, checks, not_checked = [], {}, []
    step1 = job / "step1_pdbreader.pdb"
    final = job / "namd" / "step3_input.pdb"
    for needed in (step1, final):
        if not needed.is_file():
            raise MdPackageError(f"The archive has no {needed.relative_to(job).as_posix()}; "
                                 "a finished Solution Builder job with NAMD inputs is expected.")
    checks["job_id"] = job.name
    # Everything below is measured on the FINAL system (namd/step3_input.pdb): CHARMM-GUI
    # re-orients and re-centres between steps, so mixing it with step1_pdbreader.pdb compares
    # two different coordinate frames (it put the closest modelled loop at 7.28 A instead of 4.76).

    # 1. is this the job of this package?
    uploaded = next((p for p in job.glob("*.pdb") if p.name.lower() == delivered.name.lower()), None)
    digest = hashlib.sha256(uploaded.read_bytes()).hexdigest() if uploaded else None
    checks["uploaded_file_is_the_delivered_one"] = bool(digest and digest == package["files"][delivered.name]["sha256"])
    if not checks["uploaded_file_is_the_delivered_one"]:
        failures.append(f"the {delivered.name} inside the archive is not the file build_md_package wrote "
                        + ("(different SHA-256)" if digest else "(it is not there)")
                        + ": this .tgz is from another job or another input")

    # 2. the ANISOU trap: what CHARMM-GUI handed to Modeller
    intermediates = sorted(job.glob("*_pro?.orig.pdb"))
    if not intermediates:
        not_checked.append("the input CHARMM-GUI gives Modeller: the archive has no *_pro?.orig.pdb, "
                           "so the ANISOU trap could not be looked for")
    for intermediate in intermediates:
        bad = malformed_atom_lines(intermediate)
        checks.setdefault("modeller_input_impossible_coordinates", {})[intermediate.name] = bad
        if bad:
            failures.append(f"{intermediate.name} has {bad} atom lines with impossible or unreadable "
                            "coordinates: CHARMM-GUI read ANISOU records as atoms, so the modelled "
                            "structure cannot be trusted")

    # 3. ligand present and its pose untouched
    names = package["ligand"]["atom_names"]
    renamed = {"S1": "S", "CL1": "CL"}       # CGenFF strips the index of a lone sulfur/chlorine
    out_ligand = {a["name"]: a["xyz"] for a in read_atoms(final) if a["resname"] == LIGAND_RESNAME}
    checks["ligand_atoms_in_namd_system"] = len(out_ligand)
    checks["ligand_atoms_delivered"] = package["ligand"]["atoms"]
    if not out_ligand:
        failures.append("the system has no ligand: the hetero chain was left out (it is out by default, "
                        "without a warning)")
    else:
        heavy = [(i, n) for i, n in enumerate(names) if not n[:1].isdigit() and not n.startswith("H")]
        pairs = [(package["ligand"]["coordinates"][i], out_ligand[renamed.get(n, n)])
                 for i, n in heavy if renamed.get(n, n) in out_ligand]
        checks["ligand_heavy_atoms_matched"] = f"{len(pairs)} of {len(heavy)}"
        if len(pairs) < len(heavy):
            failures.append(f"only {len(pairs)} of {len(heavy)} ligand heavy atoms were found by name in the system")
        if len(pairs) >= 2:
            moved = internal_distance_change([p[0] for p in pairs], [p[1] for p in pairs])
            checks["ligand_pose_internal_change"] = moved
            if moved["max_angstrom"] > LIGAND_DD_TOLERANCE:
                failures.append(f"the docked pose changed shape (max {moved['max_angstrom']} A > "
                                f"{LIGAND_DD_TOLERANCE} A): the system is not simulating the pose that was docked")

    # 3b. the pose is still IN the pocket. Shape alone does not say so: a ligand translated away,
    # or mirrored, keeps every internal distance (both reproduced by the reviewer, 2026-11-08).
    delivered_pocket = package.get("pocket_residues_delivered") or []
    if out_ligand and delivered_pocket:
        out_atoms = read_atoms(final)
        returned = contact_residues([a for a in out_atoms if a["segid"].startswith("PRO")],
                                    [a for a in out_atoms if a["resname"] == LIGAND_RESNAME])
        protein_heavy = [a["xyz"] for a in out_atoms
                         if a["segid"].startswith("PRO") and not a["element"].startswith("H")]
        ligand_heavy = [a["xyz"] for a in out_atoms if a["resname"] == LIGAND_RESNAME
                        and not a["element"].startswith("H") and not a["name"].startswith("LP")]
        if protein_heavy and ligand_heavy:
            closest = round(min_distance(protein_heavy, ligand_heavy), 2)
            checks["min_protein_ligand_heavy_distance_angstrom"] = closest
            if closest < MIN_HEAVY_CONTACT:
                failures.append(f"the ligand overlaps the protein ({closest} A between heavy atoms, "
                                f"below {MIN_HEAVY_CONTACT} A): the system would blow up")
        as_chain = {f"{r.split(':')[0].removeprefix('PRO')}:{r.split(':')[1]}" for r in returned}
        kept = sorted(set(delivered_pocket) & as_chain)
        overlap = len(kept) / len(delivered_pocket)
        checks["pocket_residues"] = {"delivered": len(delivered_pocket), "returned": len(as_chain),
                                     "still_touching": len(kept), "overlap": round(overlap, 3)}
        if overlap < MIN_POCKET_OVERLAP:
            failures.append(f"the ligand is not in its pocket any more: only {len(kept)} of "
                            f"{len(delivered_pocket)} residues that touched it still do "
                            f"(overlap {overlap:.2f} < {MIN_POCKET_OVERLAP})")
    elif out_ligand:
        not_checked.append("whether the ligand is still in its pocket: the package records no "
                           "contact residues (it was built by an older version)")

    # 3c. the ligand file CHARMM-GUI actually used
    used_sdf = job / "lig" / "LIG.sdf"
    if not used_sdf.is_file():
        not_checked.append("the SDF CHARMM-GUI used: the archive has no lig/LIG.sdf, so the ligand "
                           "may have been read from RCSB, Mol2 or Antechamber instead of the "
                           "delivered file")
    elif not delivered.with_name(f"{package['pdb_id']}_ligand.sdf").is_file():
        # the archive may be perfectly fine: what is missing is our side of the comparison, and
        # blaming CHARMM-GUI for it would be a diagnosis nothing was measured for
        not_checked.append(f"whether lig/LIG.sdf is the delivered one: {package['pdb_id']}"
                           "_ligand.sdf is no longer in the package folder")
    else:
        delivered_sdf = delivered.with_name(f"{package['pdb_id']}_ligand.sdf")
        same = (used_sdf.read_text(encoding="utf-8", errors="replace").split()
                == delivered_sdf.read_text(encoding="utf-8", errors="replace").split())
        checks["ligand_sdf_is_the_delivered_one"] = same
        if not same:
            failures.append("lig/LIG.sdf is not the SDF build_md_package wrote: the bond orders and "
                            "protonation CGenFF parameterised are not the ones that were docked")

    # 4. protein geometry, without superposition
    before = {(a["chain"], a["resid"]): a["xyz"] for a in read_atoms(delivered)
              if a["name"] == "CA" and a["resname"] != LIGAND_RESNAME}
    after = {(a["segid"], a["resid"]): a["xyz"] for a in read_atoms(final) if a["name"] == "CA"}
    geometry = {}
    for chain in sorted({c for c, _ in before}):
        segid = f"PRO{chain}"
        keys = [k for k in sorted(before) if k[0] == chain and (segid, k[1]) in after]
        if len(keys) < 2:
            failures.append(f"chain {chain} could not be matched to segment {segid} in the output")
            continue
        old = [before[k] for k in keys]
        new = [after[(segid, k[1])] for k in keys]
        change = internal_distance_change(old, new)
        rg_before, rg_after = radius_of_gyration(old), radius_of_gyration(new)
        geometry[segid] = {"conserved_ca": len(keys), **change,
                           "radius_of_gyration_angstrom": {"delivered": round(rg_before, 2),
                                                           "returned": round(rg_after, 2)}}
        if change["mean_angstrom"] > MEAN_DD_TOLERANCE or abs(rg_after - rg_before) > RG_TOLERANCE:
            failures.append(f"{segid} came back distorted (mean CA-CA change {change['mean_angstrom']} A, "
                            f"radius of gyration {rg_before:.2f} -> {rg_after:.2f} A)")
    checks["protein_geometry"] = geometry

    # 5. gaps filled and no artificial bond left
    protein = [a for a in read_atoms(final) if a["segid"].startswith("PRO")]
    bonds = long_peptide_bonds(protein)
    checks["peptide_bonds_over_limit"] = bonds
    if bonds:
        failures.append(f"{len(bonds)} peptide bond(s) longer than {MAX_PEPTIDE_BOND} A "
                        f"(worst {max(b['distance_angstrom'] for b in bonds)} A): the gaps were bridged, not filled")
    present = {(a["segid"], a["resid"]) for a in protein}
    was_there = {(a["chain"], a["resid"]) for a in read_atoms(delivered) if a["resname"] != LIGAND_RESNAME}
    modelled, missing, never_missing = [], [], []
    for run in package.get("to_model") or []:
        segid, label = f"PRO{run['chain']}", f"{run['chain']} {run['first']}-{run['last']}"
        numbers = range(run["first"], run["last"] + 1)
        # "present in the output" only means modelled if it was absent from what was uploaded,
        # otherwise the check passes by construction for a run that was never a gap
        if any((run["chain"], n) in was_there for n in numbers):
            never_missing.append(label)
        elif [n for n in numbers if (segid, n) not in present]:
            missing.append(label)
        else:
            modelled.append(label)
    checks["modelled_runs_present"] = modelled
    if missing:
        failures.append(f"the runs {', '.join(missing)} are still absent: 'Model missing residues' was not applied")
    if never_missing:
        failures.append(f"the runs {', '.join(never_missing)} were already in the uploaded file, so "
                        "they cannot have been modelled: the package was built from a wrong list of gaps")

    # 6. the modelled loops must not sit on the ligand
    if out_ligand and modelled:
        ligand_xyz = [v for k, v in out_ligand.items() if not k.startswith(("H", "LP"))]
        worst = None
        for run in package.get("to_model") or []:
            segid = f"PRO{run['chain']}"
            loop = [a["xyz"] for a in protein
                    if a["segid"] == segid and run["first"] <= a["resid"] <= run["last"]
                    and not a["element"].startswith("H")]
            if loop and ligand_xyz:
                distance = min_distance(loop, ligand_xyz)
                entry = {"run": f"{run['chain']} {run['first']}-{run['last']}",
                         "min_distance_angstrom": round(distance, 2)}
                if worst is None or distance < worst["min_distance_angstrom"]:
                    worst = entry
        checks["closest_modelled_run_to_ligand"] = worst
        if worst and worst["min_distance_angstrom"] < MIN_HEAVY_CONTACT:
            failures.append(f"the modelled run {worst['run']} clashes with the ligand "
                            f"({worst['min_distance_angstrom']} A)")

    checks["system"] = _system_values(job)
    if not failures:
        verdict = "Every check that could run passed; the system can be submitted."
        if not_checked:
            verdict += (f" {len(not_checked)} check(s) could not run, listed in not_checked: read "
                        "them before submitting.")
    else:
        verdict = f"{len(failures)} check(s) failed: do NOT submit this system."
    return {
        "passes": not failures,
        "failures": failures,
        "not_checked": not_checked,
        "checks": checks,
        "verdict": verdict,
        "modelled_not_experimental": modelled or "none",
        "notes": [CALCULATED,
                  "The values under `system` (box, ions, waters, temperature, timestep, steps, "
                  "CGenFF penalties) are reported, not judged: no threshold is applied to them.",
                  "Geometry is compared without superposition (CA-CA distance matrix and radius of "
                  "gyration): CHARMM-GUI reorients the system and a fit on a distorted structure hides "
                  "the distortion.",
                  "Reaching the end of the wizard proves nothing: of the four jobs of the validation, "
                  "two finished with a broken system (ESTADO 8.4)."],
    }


# ------------------------------------------------------------------ CLI

def main(argv=None):
    """CLI: the same functions as the MCP tools, result as JSON (separate process, stdout is free)."""
    ap = argparse.ArgumentParser(prog="python -m protein_mcp.files.md_package",
                                 description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build", help="write the CHARMM-GUI delivery folder")
    check = sub.add_parser("check", help="verify the .tgz CHARMM-GUI returns")
    for parser in (build, check):
        parser.add_argument("uniprot_id")
        parser.add_argument("--structure", default="")
        parser.add_argument("--pocket", default="")
    build.add_argument("--conformer", default="", help="names the docking run when the pocket has "
                                                       "more than one")
    build.add_argument("--pose", type=int, default=1, help="1 = best scored (default)")
    check.add_argument("--archive", required=True)
    a = ap.parse_args(argv)
    if a.command == "build":
        result = build_md_package(a.uniprot_id, a.structure, a.pocket, a.conformer, a.pose)
    else:
        result = check_charmm_gui_output(a.uniprot_id, a.archive, a.structure, a.pocket)
    sys.stdout.buffer.write((json.dumps(result, ensure_ascii=False, indent=1) + "\n").encode("utf-8"))
    return 0 if result.get("ok") and result.get("passes", True) else 1


if __name__ == "__main__":
    sys.exit(main())
