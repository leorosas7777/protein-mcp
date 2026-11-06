"""dock_ligand: re-dock a crystal ligand into its own pocket with AutoDock Vina (Phase 3, step 2).

Receptor, box and ligand come from prep/manifest.json (the pocket is checked with
files.choice.resolve, the same validation as select_structure), never from the chat. The ligand
is the crystal pose of that pocket (prep/ligands/<PDB>_<comp>_<label>.sdf): either a conformer
regenerated from its own atoms and bonds (default; the pre-registered control, ESTADO 7.1) or
its crystal conformation (biased in favour). One ligand per call.

Measured in Phase 3 step 1 (ESTADO 8.3) and enforced here, before Vina is started:
- a pocket with large_box or spans_target_chains is refused (5JSM 6NB J did not finish in 26.7 min);
- a receptor that kept a nucleotide, cofactor or metal is refused naming it (Meeko has no
  template for GDP/ATP; never deleted here);
- the receptor is prepared by Meeko with --delete_bad_res_from_box_radius (it names every
  residue it deletes or refuses); -x / --delete_bad_res is never used;
- own guard: every ligand atom and the centroid inside the box, receptor atoms near its centre
  (Vina docks a ligand in another frame without complaint; only --score_only refuses it);
- RDKit reads files only through Python file objects (its path readers fail on non-ASCII paths
  in Windows) and mk_prepare_ligand is not used (it reads by path): Meeko's Python API is;
- the multi-pose output is parsed whole and the number of poses is checked against MODEL
  records (splitting on $$$$ lost 8 of 9 poses without an error);
- RMSD with rdMolAlign.CalcRMS (in place, symmetry-corrected), never GetBestRMS, plus a
  permanent negative control (the crystal pose moved 5 Å must give 5 Å).

Results go to prep/dock/<PDB>_<comp>_<label>/ and to a `docking` block of the structure in the
manifest (re-read before writing; refused if prepare_structures rebuilt it meanwhile).
prepare_structures keeps a block only if the raw file and the box are the same (carry_over).
Pure function + CLI: python -m protein_mcp.files.dock P15056 --structure 3OG7 --pocket "032 C".
Imports common, files.choice (+ rdkit, meeko, numpy).
"""
import argparse
import asyncio
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice

VINA_ENV = "PROTEIN_MCP_VINA"
# fixed by the pre-registration (ESTADO 7.1); exhaustiveness and seed can be changed per call
VINA_FIXED = {"num_modes": 9, "energy_range": 3, "min_rmsd": 1, "scoring": "vina"}
DEFAULT_EXHAUSTIVENESS = 32
MAX_EXHAUSTIVENESS = 64
DEFAULT_SEED = 42
CONFORMERS = ("regenerated", "crystal")
# Meeko deletes an unmatched residue only if ALL its atoms are more than this outside the box
# (meeko/polymer.py: minimum distance of any atom to the box surface; its --help says "center of
# mass", the code does not); nearer it refuses and names it. 8 Å = Vina's interaction cutoff, so
# no atom of a deleted residue was within reach of any point of the box. Chosen by Claude Code
# 2026-11-05 (step 1 measured 0 and 5).
DELETE_BAD_RES_RADIUS = 8.0
RECEPTOR_NEAR = 8.0            # Å from the box centre: the guard counts receptor atoms within it
RMSD_CRITERION = 2.0           # Å, pre-registered (ESTADO 7.1)
NEGATIVE_SHIFT = 5.0           # Å, permanent negative control
RECEPTOR_TIMEOUT = 300         # s
VINA_TIMEOUT = 900             # s; step 1: 70-95 s, 213-219 s with a loaded machine (2026-11-05)
PROGRESS_EVERY = 10.0          # s
RECEPTOR_LOG = "receptor_meeko.log"
WORKER_LOG = "worker.log"
# step 2.1: a run is a detached process, so the chat is not waiting for it. Measured in Desktop
# on 2026-11-05: the transport cuts the call at 60 s (remote bridge) or 240 s (normal chat), and
# the request is cancelled, which used to kill Vina with it (ESTADO 4, tool 12).
HEARTBEAT_EVERY = 10.0         # s: the worker rewrites its state file this often
STALE_AFTER = 60.0             # s without a heartbeat: the run counts as vanished (shutdown, kill)
STARTUP_GRACE = 180.0          # s allowed for the worker to report itself before that verdict:
# a process that has not started yet is not a process that died. Importing numpy, RDKit and Meeko
# in a cold process on a loaded 8 GB machine took over a minute once while the tests ran.
LAUNCH_CONFIRM_SECONDS = 6.0   # s waited for the worker to report "running"
STOP_WAIT_SECONDS = 20.0       # s waited for the worker to react to the stop marker
LOCK_NAME = ".dock_lock.json"  # one Vina at a time on this machine (12 cores, 8 GB)
RUN_STATES = ("queued", "running", "finished", "failed", "stopped", "vanished")
CALCULATED = ("Calculated by AutoDock Vina (scoring function 'vina'): a docking score is not a measured "
              "affinity and is not comparable with ChEMBL IC50/Ki/Kd values of the same dossier.")

_BAD_RES_RE = re.compile(r"Bad res (\S+?:\S+) is within radius of box edge \(([\d.]+) <=")
_IGNORED_RE = re.compile(r"Ignored (\S+?:\S+): ([\d.]+) A outside box")
_NO_MATCH_RE = re.compile(r"No template matched for residue_key='([^']+)'")
_VINA_VERSION_RE = re.compile(r"AutoDock Vina v?(\S+)")


class DockError(Exception):
    """A refusal or failure with a message for the chat; `extra` goes into the reply."""
    def __init__(self, message: str, **extra):
        super().__init__(message)
        self.extra = extra


class DockStopped(DockError):
    """The run was stopped on request (step 2.1): a state of its own, not a failure."""


# ---------------------------------------------------------------- environment

def find_vina() -> dict:
    """The Vina executable: PROTEIN_MCP_VINA, else `vina` on PATH."""
    env = os.environ.get(VINA_ENV, "").strip().strip('"')
    if env:
        path = Path(env)
        if not path.is_file():
            return {"found": False, "source": VINA_ENV, "error": f"{VINA_ENV} points to {env}, which is not a file."}
        return {"found": True, "source": VINA_ENV, "path": str(path)}
    found = shutil.which("vina")
    if found:
        return {"found": True, "source": "PATH", "path": found}
    return {"found": False, "source": None,
            "error": f"AutoDock Vina not found: {VINA_ENV} is not set and there is no `vina` on PATH."}


def vina_identity(path: str) -> dict:
    """Version as printed by the program, and the SHA-256 of the executable."""
    try:
        out = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=30,
                             stdin=subprocess.DEVNULL)
        m = _VINA_VERSION_RE.search(out.stdout + out.stderr)
        version = m.group(1) if m else None
    except (OSError, subprocess.SubprocessError) as exc:
        raise DockError(f"Vina could not be started ({path}): {exc!r}")
    if not version:
        raise DockError(f"{path} did not print an AutoDock Vina version.")
    return {"path": path, "version": version, "sha256": sha256(Path(path))}


def ascii_work_root(vina_path: str = "") -> Path:
    """A folder with an ASCII-only path for the intermediate files (ESTADO 8.3: RDKit and
    mk_prepare_* fail on non-ASCII paths in Windows). %TEMP% is under the user profile, so with an
    accented user name it is not ASCII: then <folder of the Vina executable>/tmp."""
    candidates = [Path(tempfile.gettempdir())]
    if vina_path:
        candidates.append(Path(vina_path).resolve().parent / "tmp")
    for c in candidates:
        if str(c).isascii():
            c.mkdir(parents=True, exist_ok=True)
            return c
    raise DockError("No folder with an ASCII-only path for temporary files: the temporary folder ("
                    f"{candidates[0]}) and the Vina folder are not ASCII. Install Vina in a folder like C:\\vina.")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------- pure steps

def read_sdf(path: Path):
    """The single molecule of an SDF, through a file object (never a path: non-ASCII paths)."""
    from rdkit import Chem
    with open(path, "rb") as fh:
        mols = list(Chem.ForwardSDMolSupplier(fh, removeHs=False))
    if len(mols) != 1 or mols[0] is None:
        raise DockError(f"{path.name}: expected one readable molecule, got {len(mols)} "
                        f"({sum(m is None for m in mols)} unreadable).")
    mol = mols[0]
    if mol.GetNumHeavyAtoms() == 0:
        raise DockError(f"{path.name}: no heavy atom.")
    return mol


def build_ligand(crystal, conformer: str, box: dict, seed: int):
    """(molecule with H and 3D coordinates, record). 'crystal' keeps the crystal coordinates;
    'regenerated' keeps atoms and bonds only, embeds a new conformer (ETKDGv3 + MMFF) and moves
    its centroid to the box centre (the translation is recorded)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    record = {"conformer": conformer, "smiles": Chem.MolToSmiles(Chem.RemoveHs(crystal)),
              "heavy_atoms": crystal.GetNumHeavyAtoms()}
    if conformer == "crystal":
        mol = Chem.AddHs(crystal, addCoords=True)
        record["note"] = ("Crystal conformation with hydrogens added: biased in favour of reproducing the "
                          "crystal pose; not the pre-registered control.")
        return mol, record
    mol = Chem.AddHs(Chem.RemoveHs(crystal))
    mol.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise DockError("RDKit could not generate a 3D conformer from the crystal ligand's atoms and bonds (ETKDGv3).")
    mmff = AllChem.MMFFOptimizeMolecule(mol) if AllChem.MMFFHasAllMoleculeParams(mol) else None
    conf = mol.GetConformer()
    pos = conf.GetPositions()
    shift = np.array(box["center"], dtype=float) - pos.mean(axis=0)
    for i, p in enumerate(pos + shift):
        conf.SetAtomPosition(i, p.tolist())
    record.update({
        "method": f"RDKit ETKDGv3 (randomSeed {int(seed)}) from the crystal ligand's atoms and bonds, then "
                  + ("MMFF94 minimisation" + (" (converged)" if mmff == 0 else " (not converged in the default iterations)")
                     if mmff is not None else "no MMFF minimisation (MMFF parameters missing)"),
        "translation_angstrom": [round(float(x), 3) for x in shift],
        "note": "Crystal coordinates discarded; the conformer was moved to the box centre and the translation recorded.",
    })
    return mol, record


def _in_box(pos, box):
    lo = np.array(box["center"], dtype=float) - np.array(box["size"], dtype=float) / 2
    hi = lo + np.array(box["size"], dtype=float)
    return ((pos >= lo) & (pos <= hi)).all(axis=1)


def receptor_coords(pdbqt_text: str):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in pdbqt_text.splitlines() if l.startswith(("ATOM", "HETATM"))])


def guard(ligand, box: dict, receptor_xyz, conformer: str) -> dict:
    """Own check before Vina: the ligand and the receptor are in the frame of the box."""
    pos = ligand.GetConformer().GetPositions()
    inside = _in_box(pos, box)
    centroid_in = bool(_in_box(pos.mean(axis=0)[None, :], box)[0])
    d = np.linalg.norm(receptor_xyz - np.array(box["center"], dtype=float), axis=1) if len(receptor_xyz) else np.array([])
    out = {
        "ligand_atoms_in_box": f"{int(inside.sum())} of {len(pos)}",
        "ligand_centroid_in_box": centroid_in,
        "receptor_atoms_within_8A_of_box_centre": int((d <= RECEPTOR_NEAR).sum()),
        "closest_receptor_atom_to_box_centre_angstrom": round(float(d.min()), 2) if len(d) else None,
        "note": ("With the crystal ligand both ligand checks pass by construction (the box was derived from it): "
                 "they only confirm that the SDF and the manifest share a frame. The receptor count is the "
                 "check that can fail." if conformer == "crystal" else
                 "The conformer was moved to the box centre, so the ligand checks confirm the move; "
                 "the receptor count is the check that can fail."),
    }
    if not inside.all() or not centroid_in:
        why = ("the regenerated conformer, centred on the box, does not fit inside it" if conformer == "regenerated"
               else "the SDF and the manifest are not in the same frame")
        raise DockError(f"Ligand outside the box ({out['ligand_atoms_in_box']} atoms inside): {why}. "
                        "Vina would dock it anyway, so it is not started.", guard=out)
    if out["receptor_atoms_within_8A_of_box_centre"] == 0:
        raise DockError(f"No receptor atom within {RECEPTOR_NEAR:g} Å of the box centre: receptor and box are "
                        "not in the same frame. Vina is not started.", guard=out)
    return out


def parse_meeko(stdout: str, stderr: str) -> dict:
    """What Meeko did with residues that match no template (it names them in both directions)."""
    text = stdout + "\n" + stderr
    return {
        "residues_without_template": sorted(set(_NO_MATCH_RE.findall(text)), key=_res_key),
        "deleted": [{"residue": r, "outside_box_angstrom": float(d)} for r, d in _IGNORED_RE.findall(text)],
        "refused": [{"residue": r, "distance_to_box_edge_angstrom": float(d)} for r, d in _BAD_RES_RE.findall(text)],
    }


def _res_key(r):
    chain, _, num = r.partition(":")
    m = re.match(r"-?\d+", num)
    return (chain, int(m.group()) if m else 0, num)


def prepare_receptor(receptor_pdb: Path, box: dict, work: Path) -> dict:
    """Meeko's mk_prepare_receptor (CLI, in this Python) with --delete_bad_res_from_box_radius."""
    args = [sys.executable, "-I", "-c", "from meeko.cli.mk_prepare_receptor import main; main()",
            "--read_pdb", receptor_pdb.name, "-o", "receptor", "-p",
            "--box_center", *map(str, box["center"]), "--box_size", *map(str, box["size"]),
            "--delete_bad_res_from_box_radius", str(DELETE_BAD_RES_RADIUS)]
    t0 = time.monotonic()
    try:
        r = subprocess.run(args, cwd=work, capture_output=True, text=True, timeout=RECEPTOR_TIMEOUT,
                           stdin=subprocess.DEVNULL, encoding="utf-8", errors="replace")
    except subprocess.TimeoutExpired:
        raise DockError(f"Meeko receptor preparation did not finish in {RECEPTOR_TIMEOUT} s.")
    (work / "receptor_meeko.log").write_text(f"$ mk_prepare_receptor {' '.join(args[4:])}\n--- stdout\n{r.stdout}\n--- stderr\n{r.stderr}",
                                             encoding="utf-8")
    meeko = parse_meeko(r.stdout, r.stderr)
    out = {"options": f"--read_pdb --delete_bad_res_from_box_radius {DELETE_BAD_RES_RADIUS:g}",
           "exit_code": r.returncode, "seconds": round(time.monotonic() - t0, 1), **meeko}
    if r.returncode != 0 or not (work / "receptor.pdbqt").is_file():
        if meeko["refused"]:
            msg = ("Meeko refused the receptor: residues without a matching template (incomplete side chains) "
                   f"lie within {DELETE_BAD_RES_RADIUS:g} Å of the box or inside it: "
                   + ", ".join(f"{x['residue']} ({x['distance_to_box_edge_angstrom']:g} Å)" for x in meeko["refused"])
                   + ". They are not deleted here (that would change the pocket).")
        else:
            tail = [l for l in (r.stdout + "\n" + r.stderr).splitlines() if l.strip()][-6:]
            msg = "Meeko could not prepare the receptor: " + " | ".join(tail)
        raise DockError(msg, receptor_preparation=out)
    warning = unaccounted_warning(meeko)
    if warning:
        out["warning"] = warning
    return out


def unaccounted_warning(meeko: dict):
    """The messages are parsed from text: if they do not add up, say so instead of 'nothing deleted'."""
    unaccounted = sorted(set(meeko["residues_without_template"]) - {d["residue"] for d in meeko["deleted"]}, key=_res_key)
    if unaccounted:
        return ("Meeko reported residues without a template that are not in its list of deleted residues: "
                + ", ".join(unaccounted) + f". Check {RECEPTOR_LOG} (its message format may have changed).")
    return None


def prepare_ligand_pdbqt(mol) -> tuple[str, int]:
    """Meeko's Python API (no file is read). Returns (pdbqt text, rotatable torsions)."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    with redirect_stdout(sys.stderr):
        setups = MoleculePreparation().prepare(mol)
        if len(setups) != 1:
            raise DockError(f"Meeko produced {len(setups)} ligand setups, expected 1.")
        text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise DockError(f"Meeko could not write the ligand PDBQT: {err}")
    return text, sum(1 for line in text.splitlines() if line.startswith("BRANCH"))


def vina_args(vina: str, box: dict, exhaustiveness: int, seed: int) -> list:
    c, s = box["center"], box["size"]
    return [vina, "--receptor", "receptor.pdbqt", "--ligand", "ligand.pdbqt",
            "--center_x", str(c[0]), "--center_y", str(c[1]), "--center_z", str(c[2]),
            "--size_x", str(s[0]), "--size_y", str(s[1]), "--size_z", str(s[2]),
            "--exhaustiveness", str(exhaustiveness), "--num_modes", str(VINA_FIXED["num_modes"]),
            "--energy_range", str(VINA_FIXED["energy_range"]), "--min_rmsd", str(VINA_FIXED["min_rmsd"]),
            "--seed", str(seed), "--scoring", VINA_FIXED["scoring"], "--out", "poses.pdbqt"]


def _cpu_seconds(handle) -> float:
    """kernel + user CPU seconds of a child process, from its handle (still readable after it
    exits). Windows only, None elsewhere. Validated against the real API on 2026-11-05: correct
    for a real executable (0.89 s of CPU for 0.87 s of wall time), but a trampoline reads ~0,
    because it spawns the real process (the .venv python.exe stub, and the .cmd that stands in
    for Vina in the offline tests, are trampolines)."""
    if os.name != "nt":
        return None
    try:
        import ctypes
        from ctypes import wintypes

        class _FileTime(ctypes.Structure):
            _fields_ = [("low", wintypes.DWORD), ("high", wintypes.DWORD)]

        times = [_FileTime() for _ in range(4)]       # creation, exit, kernel, user
        if not ctypes.windll.kernel32.GetProcessTimes(wintypes.HANDLE(int(handle)),
                                                      *(ctypes.byref(t) for t in times)):
            return None
        return round(sum(((t.high << 32) | t.low) / 1e7 for t in times[2:]), 1)
    except Exception:                                  # a measurement is never worth a failure
        log.debug("GetProcessTimes failed", exc_info=True)
        return None


def _kill_tree(proc) -> None:
    """Kill Vina and anything it started. Vina 1.2.7 starts nothing, but a wrapper does, and on
    Windows a surviving grandchild keeps the temporary folder locked (its working directory), so
    _remove_work would give up: measured in the offline test of a stopped run."""
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(proc.pid)], timeout=15,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.SubprocessError):
            log.debug("taskkill failed; falling back to kill()", exc_info=True)
    if proc.poll() is None:
        proc.kill()
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        log.warning("Vina did not exit after being killed")


async def run_vina(args: list, work: Path, progress=None, hooks=None) -> dict:
    """Vina in a child process, output to vina.log (never to the MCP stdio); polled, so the
    event loop keeps running, progress can be reported and a stop request is noticed.
    Returns its wall time and, on Windows, its own CPU seconds (step 2.1: this is what answers
    whether Vina is starved when it is launched from Claude Desktop)."""
    t0 = time.monotonic()
    with open(work / "vina.log", "wb") as out:
        try:
            proc = subprocess.Popen(args, cwd=work, stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT)
        except OSError as exc:
            raise DockError(f"Vina could not be started: {exc!r}")
        if hooks is not None:
            hooks.vina_started(proc.pid)
        last = 0.0
        try:
            while proc.poll() is None:
                await asyncio.sleep(1.0)
                elapsed = time.monotonic() - t0
                if hooks is not None and hooks.stop_requested():
                    _kill_tree(proc)
                    out.flush()
                    raise DockStopped(f"Stopped on request after {elapsed:.0f} s; Vina was killed and no "
                                      "pose was kept. Nothing was written to the dossier.")
                if elapsed > VINA_TIMEOUT:
                    _kill_tree(proc)
                    out.flush()
                    tail = [l for l in (work / "vina.log").read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()][-4:]
                    raise DockError(f"Vina did not finish in {VINA_TIMEOUT} s and was stopped; no pose was kept. "
                                    "Last lines of its log: " + " | ".join(tail))
                if progress and elapsed - last >= PROGRESS_EVERY:
                    last = elapsed
                    try:
                        await progress(elapsed)
                    except Exception:       # progress is a courtesy to the client, never a failure
                        log.debug("progress report failed", exc_info=True)
        finally:
            # the client cancelled the request (asyncio.CancelledError) or anything else went wrong:
            # never leave Vina running on all cores (reviewer 2026-11-05: it survived a cancellation)
            if proc.poll() is None:
                _kill_tree(proc)
                log.warning("Vina stopped: the docking request ended before Vina finished")
    seconds = round(time.monotonic() - t0, 1)
    if proc.returncode != 0 or not (work / "poses.pdbqt").is_file():
        tail = [l for l in (work / "vina.log").read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()][-6:]
        raise DockError(f"Vina failed (exit code {proc.returncode}): " + " | ".join(tail))
    cpu = _cpu_seconds(proc._handle) if os.name == "nt" else None
    timing = {"wall_seconds": seconds, "logical_cpus": os.cpu_count()}
    # Docking always burns CPU, so a zero is never a measurement of Vina: it means this handle
    # belongs to a wrapper that started the real process and exited. Reporting 0 would read as a
    # starved Vina, which is the very question these numbers exist to answer (reviewer H12/F1).
    # The wall time is NOT used to decide this: run_vina polls once a second, so a fast wrapper
    # lands on exactly 1.0 s and any threshold on it is a coin toss (reviewer F1).
    if cpu is None:
        timing["cpu_seconds"] = None
        timing["cpu_seconds_note"] = ("Not measured: the CPU time of the process could not be read on this "
                                      "system. Say nothing about speed from this run.")
    elif cpu <= 0:
        timing["cpu_seconds"] = None
        timing["cpu_seconds_note"] = ("Not measurable: the executable in PROTEIN_MCP_VINA used no CPU itself, "
                                      "so it started the real process and exited (a wrapper such as a .cmd or "
                                      ".bat). Say nothing about speed from this run.")
    else:
        timing["cpu_seconds"] = cpu
        if seconds:
            timing["cpu_per_wall_second"] = round(cpu / seconds, 2)
    return timing


def read_poses(pdbqt_text: str):
    """(molecule with one conformer per pose, scores). The whole multi-pose file is parsed and
    the count is checked against its MODEL records."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    expected = sum(1 for line in pdbqt_text.splitlines() if line.startswith("MODEL"))
    with redirect_stdout(sys.stderr):
        pmol = PDBQTMolecule(pdbqt_text, skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pmol)
    scores = [float(p.score) for p in pmol]
    if len(mols) != 1 or mols[0] is None:
        raise DockError(f"Vina's poses could not be rebuilt as one molecule ({len(mols)} returned).")
    mol = mols[0]
    if not expected or mol.GetNumConformers() != expected or len(scores) != expected:
        raise DockError(f"Pose count mismatch: {expected} MODEL records, {mol.GetNumConformers()} conformers, "
                        f"{len(scores)} scores.")
    return mol, scores


def poses_sdf(mol, scores) -> str:
    """One record per pose with its calculated score."""
    from rdkit import Chem
    blocks = []
    for i, (conf, score) in enumerate(zip(mol.GetConformers(), scores), 1):
        m = Chem.Mol(mol)
        m.SetProp("_Name", f"pose {i}")
        block = Chem.MolToMolBlock(m, confId=conf.GetId())
        blocks.append(block + f">  <vina_score_kcal_mol_calculated>\n{score}\n\n>  <pose>\n{i}\n\n$$$$\n")
    return "".join(blocks)


def count_sdf_records(path: Path) -> int:
    from rdkit import Chem
    with open(path, "rb") as fh:
        return sum(1 for m in Chem.ForwardSDMolSupplier(fh, removeHs=False) if m is not None)


def rmsd_to_crystal(poses, crystal) -> list:
    """In-place, symmetry-corrected heavy-atom RMSD (rdMolAlign.CalcRMS: no superposition)."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    ref, prb = Chem.RemoveHs(crystal), Chem.RemoveHs(poses)
    return [round(float(rdMolAlign.CalcRMS(prb, ref, prbId=c.GetId())), 4) for c in prb.GetConformers()]


def negative_control(crystal, measure=None) -> dict:
    """The crystal pose against itself (0) and moved NEGATIVE_SHIFT Å (must be NEGATIVE_SHIFT),
    measured with THE SAME function as the poses (rmsd_to_crystal), so a superposing RMSD fails
    here instead of passing (reviewer 2026-11-05)."""
    from rdkit import Chem
    measure = measure or rmsd_to_crystal
    probe = Chem.Mol(Chem.RemoveHs(crystal))
    moved = Chem.Conformer(probe.GetConformer())
    for i, p in enumerate(moved.GetPositions()):
        moved.SetAtomPosition(i, (p + np.array([NEGATIVE_SHIFT, 0.0, 0.0])).tolist())
    moved.SetId(1)
    probe.AddConformer(moved, assignId=False)
    same, shifted = measure(probe, crystal)
    ok = abs(same) < 0.01 and abs(shifted - NEGATIVE_SHIFT) < 0.01
    out = {"crystal_vs_itself_angstrom": round(same, 4),
           f"crystal_moved_{NEGATIVE_SHIFT:g}A_angstrom": round(shifted, 4), "passed": ok}
    if not ok:
        raise DockError("The RMSD negative control failed: the RMSD measurement cannot be trusted.", rmsd_negative_control=out)
    return out


# ---------------------------------------------------------------- the tool

def _resolve(manifest: dict, structure: str, pocket: str):
    """(structure record, snapshot from choice.resolve, origin) or DockError."""
    structure, pocket = (structure or "").strip(), (pocket or "").strip()
    origin = "arguments"
    if not structure:
        uc = manifest.get("user_choice")
        if not uc:
            raise DockError("No structure given and no user_choice in the manifest: name a structure and a pocket, "
                            "or record the user's choice with select_structure.", options=choice.options(manifest))
        if not uc.get("pocket"):
            raise DockError(f"The user's choice ({choice.label(uc)}) has no pocket, so there is no ligand and no box.",
                            options=choice.options(manifest))
        structure = uc["structure"]
        pocket = f"{uc['pocket']['comp_id']} {uc['pocket']['label_asym_id']}"
        origin = "user_choice"
    if not pocket:
        raise DockError("pocket is required (e.g. '032 C'): the ligand and the box come from it.",
                        options=choice.options(manifest))
    snap, err = choice.resolve(manifest, structure, pocket)
    if err:
        raise DockError(f"Not docked: {err}.", options=choice.options(manifest))
    if snap.get("kind") != "experimental" or not snap.get("pocket"):
        raise DockError("Docking needs an experimental structure with a pocket.", options=choice.options(manifest))
    s = next(x for x in manifest["structures"] if x["pdb_id"] == snap["structure"])
    if origin == "user_choice":
        snap["user_chains"] = [c["auth_asym_id"] for c in manifest["user_choice"].get("chains") or []]
    return s, snap, origin


def _receptor_note(snap: dict):
    """The receptor is the whole cleaned receptor (all target chains); say so when the user chose fewer."""
    rec = list(snap.get("receptor_chains_auth") or [])
    chosen = snap.get("user_chains") or []
    if chosen and set(rec) - set(chosen):
        return (f"The receptor docked is the whole cleaned receptor, chains {', '.join(rec)}; the user's choice names "
                f"chain(s) {', '.join(chosen)}. The other chains were not removed.")
    return None


def _refusals(s: dict, snap: dict) -> None:
    pk = snap["pocket"]
    name = f"{snap['structure']} {pk['comp_id']} {pk['label_asym_id']}"
    if pk.get("large_box") or pk.get("spans_target_chains"):
        flags = [f for f, v in (("large_box", pk.get("large_box")), ("spans_target_chains", pk.get("spans_target_chains"))) if v]
        raise DockError(f"Not docked: pocket {name} is marked {' and '.join(flags)} (box "
                        + " x ".join(f"{x:g}" for x in pk["box"]["size"]) + " Å; chains "
                        + ", ".join(pk.get("residues_by_chain") or []) + "). Measured in step 1: such a pocket "
                        "(5JSM 6NB J) did not finish docking in 26.7 min. No reduced box is invented here.")
    kept = snap.get("receptor_kept_nonpolymer") or []
    if kept:
        comps = {i["label_asym_id"]: i.get("comp_id") for i in s.get("nonpolymer_instances") or []}
        raise DockError("Not docked: the cleaned receptor keeps " + ", ".join(f"{comps.get(k, '?')} ({k})" for k in kept)
                        + ". Meeko has no template for nucleotides and cofactors (GDP, ATP) and this path is not "
                        "validated; they are not deleted here.")
    if (snap.get("receptor_format") or "pdb") != "pdb":
        raise DockError(f"Not docked: the cleaned receptor of {snap['structure']} is {snap['receptor_format']}, "
                        "not PDB; only the PDB path of Meeko (--read_pdb) is validated.")
    if not pk.get("sdf"):
        raise DockError(f"Not docked: pocket {name} has no crystal ligand SDF in the manifest.")


def _files(prep: Path, s: dict, snap: dict) -> dict:
    """Receptor, raw file and SDF, inside prep/ and with the SHA-256 the manifest records."""
    def inside(rel):
        p = (prep / str(rel)).resolve()
        if prep.resolve() not in p.parents or not p.is_file():
            raise DockError(f"Not docked: {rel} is missing or outside prep/.")
        return p
    raw = inside(s["file"]["path"])
    if sha256(raw) != s["file"].get("sha256"):
        raise DockError(f"Not docked: raw/{raw.name} changed since the manifest (SHA-256): run prepare_structures again.")
    inst = next((i for i in s.get("nonpolymer_instances") or [] if i.get("label_asym_id") == snap["pocket"]["label_asym_id"]), {})
    sdf = inst.get("sdf") or {}
    sdf_path = inside(sdf.get("path"))
    if sdf.get("sha256") and sha256(sdf_path) != sdf["sha256"]:
        raise DockError(f"Not docked: {sdf_path.name} changed since the manifest (SHA-256): run prepare_structures again.")
    return {"receptor": inside(snap["receptor"]), "raw": raw, "sdf": sdf_path, "sdf_record": sdf}


def _remove_work(work: Path, attempts: int = 10) -> None:
    """Delete the temporary folder. Right after Vina is killed, Windows can still hold its files
    for a moment (Desktop test 2026-11-05: the folder of a cancelled run stayed in %TEMP%), so retry."""
    for _ in range(attempts):
        shutil.rmtree(work, ignore_errors=True)
        if not work.exists():
            return
        time.sleep(0.3)
    log.warning("Temporary docking folder not removed: %s", work)


def _check_manifest(path: Path, manifest: dict, out_dir=None) -> dict:
    """The manifest as it is now; DockError if prepare_structures rebuilt it since it was read."""
    try:
        current = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise DockError(f"Could not re-read {path}: {exc}")
    if current.get("created_at") != manifest.get("created_at"):
        raise DockError("prepare_structures rewrote the manifest during docking; "
                        + (f"the files are in {out_dir} but the manifest was not updated. " if out_dir
                           else "nothing was written to the dossier. ")
                        + "Run dock_ligand again.")
    return current


def dock_dir_name(pdb: str, comp: str, label: str) -> str:
    return f"{pdb}_{comp}_{label}"


async def dock_ligand(uniprot_id: str, structure: str = "", pocket: str = "", conformer: str = "regenerated",
                      exhaustiveness: int = DEFAULT_EXHAUSTIVENESS, seed: int = DEFAULT_SEED, progress=None,
                      hooks=None) -> dict:
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    try:
        return await _dock(accession, structure, pocket, conformer, exhaustiveness, seed, progress, hooks)
    except DockError as exc:
        out = {"ok": False, "uniprot_id": accession, "error": str(exc), **exc.extra}
        if isinstance(exc, DockStopped):
            out["stopped"] = True
        return out


async def _dock(accession, structure, pocket, conformer, exhaustiveness, seed, progress, hooks=None) -> dict:
    # the same validation and refusals the launching call (start_dock) already applied
    inp = _inputs(accession, structure, pocket, conformer, exhaustiveness, seed)
    path, manifest, prep = inp["manifest_path"], inp["manifest"], inp["prep"]
    s, snap, origin, files = inp["structure_record"], inp["snapshot"], inp["origin"], inp["files"]
    conformer = inp["conformer"]
    ident = vina_identity(inp["vina"]["path"])
    pk = inp["pocket"]
    box = pk["box"]
    pdb, comp, label = inp["pdb"], inp["comp"], inp["label"]
    crystal = read_sdf(files["sdf"])

    root = ascii_work_root(ident["path"])
    work = Path(tempfile.mkdtemp(prefix="protein_mcp_dock_", dir=root))
    try:
        shutil.copyfile(files["receptor"], work / "receptor_in.pdb")
        if hooks is not None:
            hooks.phase("preparing the ligand and the receptor (Meeko)")
        ligand, lig_record = build_ligand(crystal, conformer, box, seed)
        receptor = await asyncio.to_thread(prepare_receptor, work / "receptor_in.pdb", box, work)
        rec_text = (work / "receptor.pdbqt").read_text(encoding="utf-8", errors="replace")
        checks = guard(ligand, box, receptor_coords(rec_text), conformer)
        lig_text, torsions = prepare_ligand_pdbqt(ligand)
        (work / "ligand.pdbqt").write_text(lig_text, encoding="utf-8")
        lig_record["rotatable_torsions"] = torsions
        rmsd_check = negative_control(crystal)
        args = vina_args(ident["path"], box, exhaustiveness, seed)
        if hooks is not None:
            hooks.phase(f"Vina running (exhaustiveness {exhaustiveness})")
        timing = await run_vina(args, work, progress, hooks)
        if hooks is not None:
            hooks.phase("reading the poses and writing the dossier")
        poses_text = (work / "poses.pdbqt").read_text(encoding="utf-8", errors="replace")
        poses, scores = read_poses(poses_text)
        rmsds = rmsd_to_crystal(poses, crystal)
        (work / "poses.sdf").write_text(poses_sdf(poses, scores), encoding="utf-8")
        if count_sdf_records(work / "poses.sdf") != len(scores):
            raise DockError("The poses SDF does not read back with every pose.")

        # refuse BEFORE touching prep/dock/ if prepare_structures rebuilt the manifest meanwhile:
        # otherwise a kept block of the new manifest could point to overwritten files (reviewer 2026-11-05)
        _check_manifest(path, manifest)
        # copy to the dossier (Python file APIs handle non-ASCII paths)
        run_id = f"{comp}_{conformer}"
        out_dir = prep / "dock" / dock_dir_name(pdb, comp, label)
        out_dir.mkdir(parents=True, exist_ok=True)
        names = {"receptor.pdbqt": "receptor.pdbqt", "receptor_meeko.log": "receptor_meeko.log",
                 "ligand.pdbqt": f"ligand_{run_id}.pdbqt", "poses.pdbqt": f"poses_{run_id}.pdbqt",
                 "poses.sdf": f"poses_{run_id}.sdf", "vina.log": f"vina_{run_id}.log"}
        written = {}
        for src, dst in names.items():
            shutil.copyfile(work / src, out_dir / dst)
            written[dst] = {"path": (out_dir / dst).relative_to(prep).as_posix(), "sha256": sha256(out_dir / dst)}
    finally:
        _remove_work(work)

    top = rmsds[0]
    # the pre-registration (ESTADO 7.1) fixes the conformer AND exhaustiveness 32, seed 42
    control_valid = conformer == "regenerated" and exhaustiveness == DEFAULT_EXHAUSTIVENESS and seed == DEFAULT_SEED
    now = utc_now()
    run = {
        "generated_by": "dock_ligand", "computed_at": now, "selected_by": origin,
        "ligand": {"source": "crystal pose of this pocket", "comp_id": comp, "label_asym_id": label,
                   "auth_asym_id": pk.get("auth_asym_id"), "name": pk.get("name"), "synonyms": pk.get("synonyms") or [],
                   "sdf": files["sdf_record"].get("path"), "sdf_sha256": sha256(files["sdf"]), **lig_record},
        "receptor": {"source": snap["receptor"], "chains_auth": snap.get("receptor_chains_auth"),
                     "user_choice_chains": snap.get("user_chains"), "chains_note": _receptor_note(snap), **receptor,
                     "residue_ids": "chain:number as in the cleaned receptor PDB (auth numbering)"},
        "guard": checks,
        "vina": {"version": ident["version"], "executable_sha256": ident["sha256"],
                 "exhaustiveness": exhaustiveness, "seed": seed, **VINA_FIXED, "cpu": "all (Vina default)",
                 **timing, "box": {"center": box["center"], "size": box["size"]}},
        "poses": [{"mode": i, "vina_score_kcal_mol": sc, "rmsd_to_crystal_pose_angstrom": r}
                  for i, (sc, r) in enumerate(zip(scores, rmsds), 1)],
        "redocking": {
            "top_pose_rmsd_angstrom": top, "criterion_angstrom": RMSD_CRITERION,
            "within_criterion": top <= RMSD_CRITERION,
            "is_preregistered_control": control_valid,
            "rmsd_method": "rdMolAlign.CalcRMS, heavy atoms, symmetry-corrected, in place (no superposition)",
            "rmsd_negative_control": rmsd_check,
            "scope": f"One re-docking of one ligand into one structure ({pdb} {comp} {label}); not the accuracy of the method.",
        },
        "files": written,
        "values": "calculated",
        "note": CALCULATED,
    }
    if conformer == "crystal":
        run["redocking"]["note"] = ("Started from the crystal conformation: biased in favour, so it does not "
                                    "validate the set-up; use conformer 'regenerated' for the control.")
    elif not control_valid:
        run["redocking"]["note"] = (f"exhaustiveness {exhaustiveness} / seed {seed} differ from the pre-registered "
                                    f"{DEFAULT_EXHAUSTIVENESS} / {DEFAULT_SEED}: the RMSD is reported against the same "
                                    "2.0 Å threshold, but this run is not the pre-registered control.")

    pocket_key = f"{comp} {label}"
    try:
        current = _check_manifest(path, manifest, out_dir)
        cs = next(x for x in current["structures"] if x.get("pdb_id") == pdb)
        block = (cs.setdefault("docking", {})).setdefault(pocket_key, {"box": box, "raw_sha256": s["file"]["sha256"], "runs": {}})
        block.update(box=box, raw_sha256=s["file"]["sha256"], dir=out_dir.relative_to(prep).as_posix())
        previous = block["runs"].get(run_id)
        if previous:
            run["replaced_run_computed_at"] = previous.get("computed_at")
        block["runs"][run_id] = run
        choice.write_manifest(path, current)
    except DockError:
        raise
    except (OSError, ValueError, StopIteration) as exc:
        log.exception("Could not write %s", path)
        raise DockError(f"Docking finished but the manifest could not be written: {exc!r}")
    return compact(accession, path, pdb, pocket_key, run)


def compact(accession: str, path: Path, pdb: str, pocket_key: str, run: dict) -> dict:
    rd = run["redocking"]
    rec = run["receptor"]
    lig = run["ligand"]
    return {
        "ok": True, "uniprot_id": accession, "manifest_path": str(path),
        "docked": f"{pdb} pocket {pocket_key}" + (f" ({lig['name']})" if lig.get("name") else "")
                  + f", conformer {lig['conformer']}" + (" (from user_choice)" if run["selected_by"] == "user_choice" else ""),
        "redocking": {
            "top_pose_rmsd_angstrom": rd["top_pose_rmsd_angstrom"],
            "criterion": f"<= {rd['criterion_angstrom']:g} Å (threshold pre-registered for the control "
                         f"run: regenerated conformer, exhaustiveness {DEFAULT_EXHAUSTIVENESS}, seed {DEFAULT_SEED})",
            "within_criterion": rd["within_criterion"],
            "is_preregistered_control": rd["is_preregistered_control"],
            **({"note": rd["note"]} if rd.get("note") else {}),
            "scope": rd["scope"],
        },
        "poses": [f"mode {p['mode']}: {p['vina_score_kcal_mol']:.2f} kcal/mol (calculated), "
                  f"RMSD to crystal pose {p['rmsd_to_crystal_pose_angstrom']:.2f} Å" for p in run["poses"]],
        "poses_returned": f"{len(run['poses'])} of num_modes {run['vina']['num_modes']} (energy_range {run['vina']['energy_range']} kcal/mol)",
        "receptor_preparation": {
            "receptor_chains": rec.get("chains_auth"),
            **({"chains_note": rec["chains_note"]} if rec.get("chains_note") else {}),
            **({"warning": rec["warning"]} if rec.get("warning") else {}),
            "deleted_residues": [f"{d['residue']} ({d['outside_box_angstrom']:g} Å outside the box)" for d in rec["deleted"]],
            "residues_without_template": rec["residues_without_template"],
        },
        "ligand": {"name": lig.get("name"), "synonyms": lig.get("synonyms") or [],
                   "smiles": lig["smiles"], "heavy_atoms": lig["heavy_atoms"], "rotatable_torsions": lig["rotatable_torsions"],
                   **({"translation_angstrom": lig["translation_angstrom"]} if "translation_angstrom" in lig else {})},
        "guard": {k: v for k, v in run["guard"].items() if k != "note"},
        "vina": {k: run["vina"][k] for k in ("version", "exhaustiveness", "seed", "wall_seconds",
                                             "cpu_seconds", "cpu_seconds_note", "cpu_per_wall_second",
                                             "logical_cpus")
                 if k in run["vina"]},
        "files": {k: v["path"] for k, v in run["files"].items()},
        "replaced_run_computed_at": run.get("replaced_run_computed_at"),
        "notes": [CALCULATED,
                  "RMSD is calculated in place against the crystal pose (no superposition); the negative control "
                  f"(crystal pose moved {NEGATIVE_SHIFT:g} Å) gave {rd['rmsd_negative_control'][f'crystal_moved_{NEGATIVE_SHIFT:g}A_angstrom']} Å.",
                  "The receptor is rigid, protonated by Meeko templates; missing loops stay missing. "
                  "Residue ids in receptor_preparation are chain:number of the cleaned receptor (auth numbering)."],
    }


def carry_over(old: dict, new: dict) -> dict:
    """Called by prepare_structures: keep a pocket's docking block if the raw file and the box are the same."""
    kept, dropped = [], []
    before = {s.get("pdb_id"): s for s in (old or {}).get("structures") or [] if s.get("docking")}
    for s in new.get("structures") or []:
        o = before.pop(s.get("pdb_id"), None)
        if not o:
            continue
        boxes = {f"{p['comp_id']} {p['label_asym_id']}": p.get("box")
                 for p in ((s.get("prepared") or {}).get("pockets") or []) if not p.get("error")}
        same_file = (s.get("file") or {}).get("ok") and (s.get("file") or {}).get("sha256")
        for key, block in o["docking"].items():
            if same_file and block.get("raw_sha256") == s["file"]["sha256"] and boxes.get(key) == block.get("box"):
                s.setdefault("docking", {})[key] = block
                kept.append(f"{s['pdb_id']} {key}")
            else:
                dropped.append(f"{s['pdb_id']} {key}: structure file or box changed, or pocket gone")
    for pdb, o in before.items():
        dropped += [f"{pdb} {key}: not in the new manifest" for key in o["docking"]]
    return {"kept": kept, "dropped": dropped}


# ---------------------------------------------------------------- background runs (step 2.1)
# Why: inside Claude Desktop the transport cancels the call at 60 s or 240 s (measured
# 2026-11-05) and Vina never finishes, so dock_ligand launches this module as a detached process
# and dock_status reads the result. The state of a run lives in a file next to its output; the
# manifest stays the durable record and is written by the worker exactly as before.

def _json_atomic(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
    os.replace(tmp, path)


def _read_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def run_id_of(comp: str, conformer: str) -> str:
    return f"{comp}_{conformer}"


def state_path(prep: Path, pdb: str, comp: str, label: str, conformer: str) -> Path:
    return prep / "dock" / dock_dir_name(pdb, comp, label) / f"run_{conformer}.state.json"


def stop_marker(state: Path) -> Path:
    return state.with_name(state.name.replace(".state.json", ".stop"))


def lock_path(accession: str) -> Path:
    """One lock for the whole installation: the limit is this machine's CPU and memory, not the
    protein, so a run of another accession must wait too."""
    return choice.manifest_path(accession).parents[2] / LOCK_NAME


def _age_seconds(stamp: str):
    """Seconds since an utc_now() stamp, or None if it cannot be read."""
    from datetime import datetime, timezone
    try:
        then = datetime.strptime(stamp, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    except (TypeError, ValueError):
        return None
    return max(0.0, (datetime.now(timezone.utc) - then).total_seconds())


def effective_state(state: dict) -> str:
    """What a recorded state means now: a run whose heartbeat stopped is "vanished" (the computer
    was shut down or the process was killed), never "running"."""
    recorded = (state or {}).get("state")
    if recorded in ("finished", "failed", "stopped"):
        return recorded
    if recorded not in ("queued", "running"):
        return "vanished"
    age = _age_seconds(state.get("heartbeat") or state.get("launched_at"))
    if age is None:
        return "vanished"
    # a run still waiting for its worker to report gets a longer grace than a run that was
    # reporting and stopped: "not started yet" and "died in the middle" are different facts
    limit = STARTUP_GRACE if (recorded == "queued" and not state.get("pid")) else STALE_AFTER
    if age <= limit:
        return recorded
    # A late heartbeat is not proof of death: on this 8 GB machine a worker importing RDKit and
    # Meeko while another run finishes has gone over a minute without writing (seen in the tests).
    # Declaring a live run vanished would also hand its lock to a second Vina, so the process is
    # checked before that verdict: its own pid, or the launcher's while it has not reported yet
    # (reviewer H8). The image name is checked too, because Windows reuses process ids.
    # its own pid once it has reported one; the launcher's only while it has not (reviewer F10),
    # because the launcher's process is the MCP server and would answer for any run
    pid = state.get("pid") or (state.get("launcher_pid") if recorded == "queued" else None)
    return recorded if _pid_alive(pid, expect="python") else "vanished"


class RunState:
    """The state file of one run, written by the worker (and created by the launcher).
    The worker writes its OWN pid here: the pid returned by Popen is the launcher's, and in this
    .venv python.exe is a trampoline that spawns the real process (measured 2026-11-05), so the
    launcher's pid identifies nothing."""

    def __init__(self, path: Path, data: dict = None):
        self.path = Path(path)
        self.data = dict(data or {})
        self._lock = threading.Lock()
        self._beat = None
        self._stop_seen = False

    def update(self, **fields) -> dict:
        with self._lock:
            self.data.update(fields)
            self.data["heartbeat"] = utc_now()
            started = self.data.get("started_at")
            if started and self.data.get("state") in ("queued", "running"):
                age = _age_seconds(started)
                if age is not None:
                    self.data["elapsed_seconds"] = round(age, 1)
            try:
                _json_atomic(self.path, self.data)
            except OSError:
                log.exception("Could not write the run state %s", self.path)
            return dict(self.data)

    # hooks used by dock_ligand / run_vina
    def phase(self, text: str) -> None:
        self.update(phase=text)

    def vina_started(self, pid: int) -> None:
        self.update(vina_pid=pid, phase="Vina running")

    def stop_requested(self) -> bool:
        if self._stop_seen:
            return True
        self._stop_seen = stop_marker(self.path).exists()
        return self._stop_seen

    def start_heartbeat(self) -> None:
        def beat():
            while not self._beat.stopped.wait(HEARTBEAT_EVERY):
                self.update()
        self._beat = threading.Thread(target=beat, name="dock-heartbeat", daemon=True)
        self._beat.stopped = threading.Event()
        self._beat.start()

    def stop_heartbeat(self) -> None:
        if self._beat is not None:
            self._beat.stopped.set()


def _pid_alive(pid, expect: str = "") -> bool:
    """Whether a process id is still running. `expect` is a fragment of the image name it must
    have (Windows: tasklist prints it), because Windows reuses process ids: without it, a run
    that died hours ago reads as alive as soon as something else gets its id, which would make
    dock_status report a phantom and keep its lock for ever (reviewer H7)."""
    if not pid:
        return False
    try:
        if os.name == "nt":
            out = subprocess.run(["tasklist", "/FI", f"PID eq {int(pid)}", "/NH"], timeout=20,
                                 capture_output=True, text=True).stdout
            if str(int(pid)) not in out:
                return False
            return expect.lower() in out.lower() if expect else True
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError, subprocess.SubprocessError):
        return False


def _lock_holder_live(held: dict) -> bool:
    """Whether the run that holds the lock is still going.
    A background run is judged by its state file: if the lock names one and it cannot be read
    (deleted with prep/dock/, or never written), the lock is expired, NOT alive. Judging it by
    the recorded pid would be judging the MCP server, which lives as long as Claude Desktop, and
    would leave the lock held for ever (reviewer H3). A run started by hand has no state file and
    is judged by its own process."""
    if held.get("state_file"):
        other = _read_json(Path(held["state_file"]))
        if other and effective_state(other) in ("queued", "running"):
            return True
        # the lock is taken a moment before the new state file is written, so a state file that is
        # OLDER than the lock belongs to the previous run and proves nothing about this holder.
        # Treat a just-taken lock as live for the startup grace (reviewer F11a).
        # the window is the launch itself (the state file is written right after the lock), not
        # the whole startup grace: a lock older than that with no readable state is expired (H3)
        taken = _age_seconds(held.get("locked_at"))
        if taken is not None and taken <= LAUNCH_CONFIRM_SECONDS:
            stamp = (other or {}).get("heartbeat") or (other or {}).get("launched_at") or ""
            if not other or stamp < (held.get("locked_at") or ""):
                return True
        return False
    return _pid_alive(held.get("pid"), expect="python")


def acquire_lock(accession: str, info: dict) -> str:
    """Create the lock file atomically (O_EXCL). Returns a token, or raises DockError naming the
    run that holds it. A lock whose run is over, or whose state cannot be read any more, is
    expired and taken over."""
    path = lock_path(accession)
    token = hashlib.sha256(f"{os.getpid()}{utc_now()}{info}".encode()).hexdigest()[:16]
    record = {"token": token, "locked_at": utc_now(), "pid": os.getpid(), **info}
    for attempt in (1, 2):
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            fd = os.open(str(path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            held = _read_json(path) or {}
            other = _read_json(Path(held["state_file"])) if held.get("state_file") else None
            live = _lock_holder_live(held)
            # also on the second attempt: another process may have taken the lock in between, and
            # deleting a live holder's lock would start a second Vina (reviewer H9)
            if live:
                raise DockError(
                    "Not started: a docking run is already in progress on this machine "
                    f"({held.get('label', 'unknown run')}, started {held.get('locked_at')}"
                    + (f", {other.get('elapsed_seconds')} s so far" if (other or {}).get("elapsed_seconds") else "")
                    + "). Only one Vina runs at a time on this machine. Check it with dock_status, or stop "
                    "it with dock_status(stop=True); a run started from the command line has to end on its own.",
                    running_run=held.get("label"), running_since=held.get("locked_at"))
            # the holder is over, or its state cannot be read any more. Take the lock over.
            log.warning("Taking over an expired docking lock held by %s", held.get("label"))
            try:
                os.unlink(str(path))
            except OSError:
                pass
            continue
        except OSError as exc:
            raise DockError(f"Not launched: the docking lock {path} could not be created: {exc}")
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(record, fh, ensure_ascii=False, indent=1)
        return token
    raise DockError(f"Not launched: the docking lock {path} is held and could not be taken over.")


def release_lock(accession: str, token: str) -> None:
    """Remove the lock only if it is still ours (never steal someone else's)."""
    path = lock_path(accession)
    held = _read_json(path)
    if held and held.get("token") == token:
        try:
            os.unlink(str(path))
        except OSError:
            log.warning("Could not remove the docking lock %s", path)


def _inputs(accession: str, structure: str, pocket: str, conformer: str, exhaustiveness, seed) -> dict:
    """Everything that can be checked without RDKit, Meeko or Vina: the same validation and the
    same refusals as a synchronous run, so the launching call answers in seconds."""
    conformer = (conformer or "").strip().lower()
    if conformer not in CONFORMERS:
        raise DockError(f"conformer must be one of {', '.join(CONFORMERS)}.")
    if isinstance(exhaustiveness, bool) or not isinstance(exhaustiveness, int) or not 1 <= exhaustiveness <= MAX_EXHAUSTIVENESS:
        raise DockError(f"exhaustiveness must be an integer from 1 to {MAX_EXHAUSTIVENESS}.")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 2**31 - 1:
        raise DockError("seed must be an integer from 0 to 2147483647 (it is always explicit: Vina's default is random).")
    path = choice.manifest_path(accession)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise DockError(f"No preparation manifest for {accession}: run prepare_structures first.")
    except (OSError, ValueError) as exc:
        raise DockError(f"Could not read {path}: {exc}")
    s, snap, origin = _resolve(manifest, structure, pocket)
    _refusals(s, snap)
    prep = path.parent
    files = _files(prep, s, snap)
    vina = find_vina()
    if not vina["found"]:
        raise DockError(vina["error"])
    pk = snap["pocket"]
    return {"manifest_path": path, "manifest": manifest, "prep": prep, "structure_record": s,
            "snapshot": snap, "origin": origin, "files": files, "vina": vina, "conformer": conformer,
            "pdb": snap["structure"], "comp": pk["comp_id"], "label": pk["label_asym_id"], "pocket": pk}


async def start_dock(uniprot_id: str, structure: str = "", pocket: str = "", conformer: str = "regenerated",
                     exhaustiveness: int = DEFAULT_EXHAUSTIVENESS, seed: int = DEFAULT_SEED) -> dict:
    """Validate, take the lock and launch this module as a detached process; answer at once."""
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    token, launched_ok = None, False
    try:
        inp = _inputs(accession, structure, pocket, conformer, exhaustiveness, seed)
        prep, pdb, comp, label = inp["prep"], inp["pdb"], inp["comp"], inp["label"]
        conformer = inp["conformer"]
        state_file = state_path(prep, pdb, comp, label, conformer)
        label_text = f"{accession} {pdb} {comp} {label} ({conformer})"
        running = effective_state(_read_json(state_file) or {})
        if running in ("queued", "running"):
            raise DockError(f"Not launched: this run ({label_text}) is already in progress. "
                            "Check it with dock_status.")
        out_dir = state_file.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        result_file = out_dir / f"run_{conformer}.result.json"
        # nothing of the previous run is destroyed until this one is really going: a request
        # refused by the lock used to delete the previous result file (reviewer F5)
        token = acquire_lock(accession, {"label": label_text, "state_file": str(state_file)})
        stop_marker(state_file).unlink(missing_ok=True)
        result_file.unlink(missing_ok=True)

        run_state = {
            "state": "queued", "uniprot_id": accession, "structure": pdb, "pocket": f"{comp} {label}",
            "conformer": conformer, "exhaustiveness": exhaustiveness, "seed": seed,
            "run_id": run_id_of(comp, conformer), "pocket_key": f"{comp} {label}",
            "selected_by": inp["origin"], "launched_at": utc_now(), "launched_by": "dock_ligand",
            "manifest_created_at": inp["manifest"].get("created_at"),
            "box": {"center": inp["pocket"]["box"]["center"], "size": inp["pocket"]["box"]["size"]},
            "result_file": result_file.name, "worker_log": WORKER_LOG, "lock_token": token,
            "dir": out_dir.relative_to(prep).as_posix(),
        }
        state = RunState(state_file, run_state)
        state.update()

        cmd = [sys.executable, "-m", "protein_mcp.files.dock", accession, "--structure", pdb,
               "--pocket", f"{comp} {label}", "--conformer", conformer,
               "--exhaustiveness", str(exhaustiveness), "--seed", str(seed),
               "--state", str(state_file), "--json-out", str(result_file)]
        launch = await asyncio.to_thread(_launch_worker, cmd, out_dir, lock_path(accession).parent)
        if not launch["ok"]:
            state.update(state="failed", error=launch["error"], finished_at=utc_now())
            release_lock(accession, token)
            raise DockError(launch["error"])
        # from here the lock belongs to the worker, which releases it when it ends: a cancellation
        # of this request must not take it away from a run that is already going
        launched_ok = True
        state.update(launch_flags=launch["flags"], launcher_pid=launch["pid"])

        confirmed = await _confirm_started(state_file)
        now = _read_json(state_file) or {}
        reply = {
            "ok": True, "uniprot_id": accession, "launched": True,
            "run": {"structure": pdb, "pocket": f"{comp} {label}", "conformer": conformer,
                    "exhaustiveness": exhaustiveness, "seed": seed, "selected_by": inp["origin"],
                    "state": effective_state(now), "worker_pid": now.get("pid"),
                    "launcher_pid": launch["pid"], "launch_flags": launch["flags"]},
            "box": run_state["box"],
            "receptor": {"file": inp["snapshot"]["receptor"], "chains_auth": inp["snapshot"].get("receptor_chains_auth"),
                         "chains_note": _receptor_note(inp["snapshot"])},
            "ligand": {"comp_id": comp, "label_asym_id": label, "name": inp["pocket"].get("name"),
                       "synonyms": inp["pocket"].get("synonyms") or []},
            "vina_executable": inp["vina"]["path"],
            "state_file": state_file.relative_to(inp["prep"]).as_posix(),
            "next": f"Call dock_status('{accession}') to see whether it is still running or to get the result.",
            "notes": [
                "Launched in a detached process, so it is not cut when this request or this chat ends, and "
                "nothing is written to the dossier until it finishes. Whether it also survives Claude "
                "Desktop closing has not been measured: if Desktop's job object kills its children, the run "
                "dies with it and dock_status will say so.",
                f"Measured on this machine: Vina alone took 60.9-124.4 s for 3OG7 032 C with "
                f"exhaustiveness {DEFAULT_EXHAUSTIVENESS}; wait at least a minute before asking.",
                "Only one docking run at a time on this machine; a second request is refused while "
                "this one lives.",
            ],
        }
        if not confirmed:
            reply["run"]["state"] = effective_state(now)
            reply["notes"].append(
                f"The worker had not reported itself {LAUNCH_CONFIRM_SECONDS:g} s after launch. It may still "
                "be starting; dock_status will say 'vanished' if it never did.")
        return reply
    except DockError as exc:
        if token:
            release_lock(accession, token)
        return {"ok": False, "uniprot_id": accession, "launched": False, "error": str(exc), **exc.extra}
    except BaseException:
        # the request was cancelled (the chat gave up) or anything else: never leave the lock
        # taken by a run that was never launched (reviewer F11b)
        if token and not launched_ok:
            release_lock(accession, token)
        raise


def _launch_worker(cmd: list, out_dir: Path, dossier_root: Path) -> dict:
    """Detached, like open_in_pymol launches PyMOL, and out of Claude Desktop's job object when
    the job allows it: a child of the job can be killed when the job closes (Desktop restarted
    its server five times on 2026-11-05). CREATE_BREAKAWAY_FROM_JOB fails with ERROR_ACCESS_DENIED
    when the job forbids it, so the launch is retried without it and `flags` records which one
    ran: that is how the first real Desktop run tells us whether breaking away is needed."""
    log_file = out_dir / WORKER_LOG
    # the worker must read the very dossier tree this launcher validated, not work it out again
    env = {**os.environ, "PYTHONUNBUFFERED": "1", "PROTEIN_MCP_DOSSIER_DIR": str(dossier_root)}
    attempts = [("detached+breakaway", 0x01000000), ("detached", 0)]
    if os.name != "nt":
        attempts = [("new_session", 0)]
    last = ""
    for name, extra in attempts:
        kwargs = {"start_new_session": True} if os.name != "nt" else {
            "creationflags": subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS | extra}
        try:
            with open(log_file, "ab") as out:
                out.write(f"\n--- {utc_now()} launching {name}: {' '.join(cmd)}\n".encode())
                out.flush()
                proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT,
                                        cwd=str(out_dir), close_fds=True, env=env, **kwargs)
            return {"ok": True, "pid": proc.pid, "flags": name}
        except OSError as exc:
            last = f"{name}: {exc}"
            log.warning("Detached launch failed (%s), trying the next flags", last)
    return {"ok": False, "error": f"Not launched: the docking process could not be started ({last})."}


async def _confirm_started(state_file: Path) -> bool:
    """Wait briefly for the worker to write its own pid and "running". The launcher's pid proves
    nothing (trampoline), so the state file is the only evidence."""
    deadline = time.monotonic() + LAUNCH_CONFIRM_SECONDS
    while time.monotonic() < deadline:
        await asyncio.sleep(0.25)
        data = _read_json(state_file) or {}
        if data.get("state") != "queued" or data.get("pid"):
            return True
    return False


def _manifest_run(prep: Path, accession: str, pdb: str, pocket_key: str, run_id: str, not_before: str = ""):
    """(manifest path, run block) of a recorded run, or (path, None). The manifest is the durable
    record: a state file is only what a process was reporting.
    `not_before` makes the block belong to THIS run: a run of the same pocket and conformer
    overwrites the same key, so without comparing `computed_at` a run that died would be handed
    the RMSD of the previous one as its own (reviewer F2)."""
    path = choice.manifest_path(accession) if accession else prep / "manifest.json"
    manifest = _read_json(path) or {}
    try:
        s = next(x for x in manifest.get("structures") or [] if x.get("pdb_id") == pdb)
        run = s["docking"][pocket_key]["runs"][run_id]
    except (StopIteration, KeyError, TypeError):
        return path, None
    if not_before and (run.get("computed_at") or "") < not_before:
        return path, None                   # an older run's block, not this one's
    return path, run


def _status_of(prep: Path, state_file: Path) -> dict:
    """One run, as the chat sees it: what it is doing, or its finished result from the manifest."""
    rel = state_file.relative_to(prep).as_posix()
    data = _read_json(state_file)
    if not data:
        # an unreadable or empty state file is not evidence of anything: saying "the computer was
        # shut down and nothing was kept" would be a story (reviewer F3)
        return {"ok": False, "state": "unknown", "state_file": rel,
                "error": f"The state file {rel} could not be read, so what that run did is unknown. "
                         "Its output may be in the manifest and in its folder; this tool cannot tell. "
                         f"Look at {rel.rsplit('/', 1)[0]}/{WORKER_LOG}."}
    state = effective_state(data)
    out = {"run": {k: data.get(k) for k in ("structure", "pocket", "conformer", "exhaustiveness", "seed",
                                            "selected_by", "launched_at", "launch_flags")},
           "state": state, "state_file": rel}
    raw_age = _age_seconds(data.get("heartbeat"))
    # a missing heartbeat is unknown, not 0 s ago (reviewer F7)
    age = round(raw_age, 1) if raw_age is not None else None
    said = f"{age:.1f} s ago" if age is not None else "at a time this tool cannot read"
    if state in ("queued", "running"):
        out.update(phase=data.get("phase"), elapsed_seconds=data.get("elapsed_seconds"),
                   worker_pid=data.get("pid"), vina_pid=data.get("vina_pid"),
                   heartbeat_age_seconds=age,
                   note="Still running in its own process; nothing is written to the dossier until it finishes.")
        if age is None or age > STALE_AFTER:
            # its process is alive (effective_state checked), but it has not reported for a while,
            # so elapsed_seconds and phase are that old. Windows also reuses process ids, so this
            # can be a dead run whose id went to another program (reviewer H7).
            out["warning"] = (f"Its last report was {said}, over the {STALE_AFTER:g} s limit, so `phase` and "
                              "`elapsed_seconds` are from then, not from now. A process with its id is "
                              "running, but this may be a run that died and an unrelated process that was "
                              "given the same id. Say that the run has not reported recently.")
        return out
    # the identifiers the manifest is looked up with; without them nothing can be checked
    ids = (data.get("structure"), data.get("pocket_key"), data.get("run_id"))
    started = data.get("started_at") or data.get("launched_at") or ""
    if state == "vanished":
        # before claiming nothing was written: the worker may have died between writing the
        # manifest and recording "finished" (reviewer H5). The block must be THIS run's (F2).
        path, run = _manifest_run(prep, data.get("uniprot_id"), *ids, not_before=started) if all(ids) else (None, None)
        if run:
            out.update(compact(data.get("uniprot_id"), path, ids[0], ids[1], run))
            out["state"] = "finished"
            out["note_on_state"] = (f"Its process stopped reporting {said} without recording the end of the "
                                    f"run, but the result below was written by this run itself (the manifest "
                                    f"block was computed at {run.get('computed_at')}, after this run started "
                                    f"at {started}). Nothing is missing except the final state.")
            return out
        if not all(ids):
            out["error"] = (f"The run stopped reporting {said} and its state file does not say which pocket or "
                            "run it was, so this tool cannot tell whether it wrote anything. Look in the "
                            f"manifest and in {rel.rsplit('/', 1)[0]}/ before launching it again.")
        elif data.get("state") == "queued" and not data.get("pid"):
            out["error"] = (f"The run never started: it was launched {said} and its process had still not "
                            f"reported itself (limit {STARTUP_GRACE:g} s), and no process of it is alive. The "
                            f"manifest has no block from it, so nothing was kept. Its output is in "
                            f"{rel.rsplit('/', 1)[0]}/{WORKER_LOG}; launch it again with dock_ligand.")
        else:
            out["error"] = (f"The run stopped without a result: its last report was {said}, over the "
                            f"{STALE_AFTER:g} s limit, and no process of it is alive. The computer may have been "
                            "shut down or the process killed. The manifest has no block from this run (a block "
                            "from an earlier run of the same pocket is not counted), so nothing was kept; "
                            "launch it again with dock_ligand.")
        out["last_phase"] = data.get("phase")
        return out
    if state in ("failed", "stopped"):
        out["error"] = data.get("error") or "No message was recorded."
        return out
    # finished: the result is in the manifest, read it from there (never from the state file)
    accession = data.get("uniprot_id")
    path, run = _manifest_run(prep, accession, *ids) if all(ids) else (None, None)
    if not run:
        out["error"] = ("The run reported success but its block is not in the manifest any more "
                        "(prepare_structures may have replaced it), or its state file does not say "
                        "which run it was. Launch it again.")
        return out
    out.update(compact(accession, path, ids[0], ids[1], run))
    out["state"] = "finished"
    return out


def _runs_in_manifest(prep: Path, accession: str) -> list:
    """Runs the manifest records, as (pdb, pocket_key, run_id, run). These are the durable record:
    a run started from the command line, or before the two-step tool existed, has no state file,
    and dock_status must not answer "nothing was ever launched" while its result sits in the
    manifest (reviewer H2)."""
    manifest = _read_json(choice.manifest_path(accession)) or {}
    found = []
    for s in manifest.get("structures") or []:
        for pocket_key, block in (s.get("docking") or {}).items():
            for run_id, run in (block.get("runs") or {}).items():
                found.append((s.get("pdb_id"), pocket_key, run_id, run))
    return found


async def dock_status(uniprot_id: str, structure: str = "", pocket: str = "", conformer: str = "",
                      stop: bool = False) -> dict:
    """Read the state of the background docking runs of one protein, or stop the live one."""
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    prep = choice.manifest_path(accession).parent
    folder = prep / "dock"
    wanted = (structure or "").strip().upper(), (pocket or "").strip(), (conformer or "").strip().lower()

    def matches(pdb, pocket_key, run_id):
        return not ((wanted[0] and (pdb or "").upper() != wanted[0])
                    or (wanted[1] and wanted[1].lower() not in (pocket_key or "").lower())
                    or (wanted[2] and not (run_id or "").endswith("_" + wanted[2])))

    all_states = sorted(folder.glob("*/run_*.state.json")) if folder.is_dir() else []
    runs, seen = [], set()
    for state_file in all_states:
        data = _read_json(state_file) or {}
        # every state file is "seen", match or not, or a filtered-out one would be counted twice:
        # once as a state file and again as a run only the manifest knows (reviewer F4)
        seen.add((data.get("structure"), data.get("pocket_key"), data.get("run_id")))
        if not matches(data.get("structure"), data.get("pocket"), data.get("run_id")):
            continue
        runs.append((state_file, data))
    # runs that only the manifest knows: started from the command line, or before this tool
    # existed. Their result is real and must not be answered with "nothing was launched" (H2).
    recorded = [(pdb, pk, rid, run) for pdb, pk, rid, run in _runs_in_manifest(prep, accession)
                if matches(pdb, pk, rid) and (pdb, pk, rid) not in seen]
    if not runs and not recorded:
        asked = ", ".join(x for x in wanted if x)
        return {"ok": False, "uniprot_id": accession,
                "error": (f"No docking run matches {asked}." if asked else
                          f"No docking run is on record for {accession}, in {folder.name}/ or in the manifest. "
                          "Launch one with dock_ligand."),
                "runs_on_record": sorted({f"{pdb} {pk} {rid}" for pdb, pk, rid, _ in
                                          _runs_in_manifest(prep, accession)}
                                         | {p.parent.name for p in all_states})}
    live = [(p, d) for p, d in runs if effective_state(d) in ("queued", "running")]
    if stop:
        if not live:
            return {"ok": True, "uniprot_id": accession, "stopped": False,
                    "note": "Nothing to stop: no matching run is live.",
                    "runs": [_status_of(prep, p) for p, _ in runs]}
        if len(live) > 1:                      # cannot happen with the lock, but never guess which one
            return {"ok": False, "uniprot_id": accession, "stopped": False,
                    "error": "More than one run looks live; name the structure and pocket to stop one.",
                    "runs": [_status_of(prep, p) for p, _ in live]}
        state_file, data = live[0]
        stop_marker(state_file).touch()
        deadline = time.monotonic() + STOP_WAIT_SECONDS
        while time.monotonic() < deadline:
            await asyncio.sleep(0.5)
            if effective_state(_read_json(state_file) or {}) not in ("queued", "running"):
                break
        final = _status_of(prep, state_file)
        # the three outcomes are different facts and must not share a note (reviewer H4)
        if final["state"] == "stopped":
            note = ("The run was stopped: Vina was killed, the temporary folder removed and nothing "
                    "written to the dossier.")
        elif final["state"] == "finished":
            note = ("It was NOT stopped: the run had already finished when the stop request was written, so "
                    "the result below is complete and is in the manifest. Report it as a finished run.")
        elif final["state"] == "failed":
            note = "It was not stopped by the request: the run had already failed. Its error is below."
        elif final["state"] in ("vanished", "unknown"):
            note = ("It was not stopped by the request: while the request was being written the run stopped "
                    "reporting on its own, so nothing of it is running any more. See its state below; do not "
                    "say it was stopped on request.")
        else:
            note = ("The stop request was written. The run had not reacted after "
                    f"{STOP_WAIT_SECONDS:g} s: ask dock_status again; if it stays live, its process "
                    f"(worker pid {data.get('pid')}, Vina pid {data.get('vina_pid')}) has to be ended "
                    "from the Task Manager. Windows reuses process ids, so check the image name first.")
        return {"ok": True, "uniprot_id": accession, "stopped": final["state"] == "stopped",
                "run": final, "note": note}
    total = len(all_states) + len({(r[0], r[1], r[2]) for r in _runs_in_manifest(prep, accession)
                                   if (r[0], r[1], r[2]) not in seen})
    return {"ok": True, "uniprot_id": accession,
            "selection": f"{len(runs) + len(recorded)} of {total} docking runs on record for {accession}",
            "live_runs": len(live),
            "runs": [_status_of(prep, p) for p, _ in runs]
                    + [_recorded_status(prep, accession, pdb, pk, rid, run) for pdb, pk, rid, run in recorded]}


def _recorded_status(prep: Path, accession: str, pdb: str, pocket_key: str, run_id: str, run: dict) -> dict:
    """A run the manifest records but no state file covers (command line, or an older version)."""
    return {"run": {"structure": pdb, "pocket": pocket_key, "conformer": run_id.rsplit("_", 1)[-1],
                    "exhaustiveness": (run.get("vina") or {}).get("exhaustiveness"),
                    "seed": (run.get("vina") or {}).get("seed"),
                    "selected_by": run.get("selected_by"), "computed_at": run.get("computed_at")},
            "state": "finished", "state_file": None,
            "note_on_state": "Recorded in the manifest with no background-run state file: it was not "
                             "started by dock_ligand (the command line, or a version before the two-step "
                             "tool). The result itself is the run's own.",
            **compact(accession, choice.manifest_path(accession), pdb, pocket_key, run)}


def _worker(a) -> dict:
    """The detached process: run the docking, keep the state file alive, release the lock."""
    state = RunState(Path(a.state), _read_json(Path(a.state)) or {})
    state.update(state="running", pid=os.getpid(), started_at=utc_now(), phase="starting")
    state.start_heartbeat()
    token = state.data.get("lock_token")
    accession = state.data.get("uniprot_id") or a.uniprot_id
    try:
        # the box and receptor dock_ligand already told the chat came from the manifest as it was
        # at launch; if prepare_structures replaced it in between, that reply is stale (H15)
        at_launch = state.data.get("manifest_created_at")
        now = (_read_json(choice.manifest_path(accession)) or {}).get("created_at")
        if at_launch and now != at_launch:
            raise DockError("prepare_structures replaced the manifest between the launch and the start of "
                            "this run, so the box reported when it was launched is not the current one. "
                            "Nothing was written to the dossier; launch it again with dock_ligand.")
        result = asyncio.run(dock_ligand(a.uniprot_id, a.structure, a.pocket, a.conformer,
                                         a.exhaustiveness, a.seed, hooks=state))
    except DockError as exc:
        result = {"ok": False, "error": str(exc), **exc.extra}
    except BaseException as exc:                       # never leave the state file saying "running"
        log.exception("The docking worker crashed")
        result = {"ok": False, "error": f"The docking process crashed: {exc!r}"}
    finally:
        state.stop_heartbeat()
    if a.json_out:
        try:
            _json_atomic(Path(a.json_out), result)
        except OSError:
            log.exception("Could not write %s", a.json_out)
    if result.get("ok"):
        final = "finished"
    elif result.get("stopped"):
        final = "stopped"
    else:
        final = "failed"
    state.update(state=final, finished_at=utc_now(), phase=None,
                 error=None if result.get("ok") else result.get("error"))
    stop_marker(Path(a.state)).unlink(missing_ok=True)
    if token:
        release_lock(accession, token)
    return result


def _cli(a) -> dict:
    """A run started by hand, in the foreground. It takes the same machine-wide lock as a
    background run, so a run from the command line and one from the chat cannot both have Vina
    going (one machine; measured: Vina alone uses every core and 382 MB)."""
    accession, err = parse_accession(a.uniprot_id)
    if err:
        return err
    label = f"{accession} {a.structure or 'user_choice'} {a.pocket} ({a.conformer}) [command line]"
    try:
        token = acquire_lock(accession, {"label": label})
    except DockError as exc:
        return {"ok": False, "uniprot_id": accession, "error": str(exc), **exc.extra}
    try:
        return asyncio.run(dock_ligand(a.uniprot_id, a.structure, a.pocket, a.conformer,
                                       a.exhaustiveness, a.seed))
    finally:
        release_lock(accession, token)


def main(argv=None):
    """CLI: the same function as the MCP tool, result as JSON (this is a separate process, not the
    MCP server, so stdout is free here)."""
    ap = argparse.ArgumentParser(prog="python -m protein_mcp.files.dock", description=__doc__.splitlines()[0])
    ap.add_argument("uniprot_id")
    ap.add_argument("--structure", default="")
    ap.add_argument("--pocket", default="")
    ap.add_argument("--conformer", default="regenerated", choices=CONFORMERS)
    ap.add_argument("--exhaustiveness", type=int, default=DEFAULT_EXHAUSTIVENESS)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--state", default="", help="state file of a background run: dock_ligand sets it when it "
                                                "launches this module as a detached process (step 2.1)")
    ap.add_argument("--json-out", dest="json_out", default="", help="write the result as JSON to this file")
    a = ap.parse_args(argv)
    if a.state:                     # background run: the state file and the lock are this process's job
        result = _worker(a)
    else:
        result = _cli(a)
        if a.json_out:
            _json_atomic(Path(a.json_out), result)
    sys.stdout.buffer.write((json.dumps(result, ensure_ascii=False, indent=1) + "\n").encode("utf-8"))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
