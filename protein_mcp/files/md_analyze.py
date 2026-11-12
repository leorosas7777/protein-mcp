"""Measure a finished NAMD trajectory: what the simulation did, never what it means.

This module reads the files a `build_md_job` run leaves behind (a CHARMM-GUI PSF and the DCD
trajectories) and reports five metrics. It is the counterpart of `md_run.py`: that module writes
the job, this one reads its output. Nothing here contacts a cluster or a database.

What it reports, and why each one:

1. `protein_rmsd`   - the backbone against the first analysed frame. A run whose protein has not
                      reached a plateau is not equilibrated, so every other number below is
                      premature. This is the control, not a result.
2. `ligand_rmsd`    - the ligand's heavy atoms, superposed on the POCKET backbone and never on
                      the ligand itself. Superposing on the ligand hides the very motion the
                      metric is for: the ligand would sit still by construction while leaving
                      the site.
3. `com_distance`   - between the centre of mass of the ligand and that of the pocket residues,
                      with periodic images undone by continuity between consecutive frames. This
                      separates "the ligand turned inside the site" from "the ligand left", which
                      RMSD alone cannot do.
4. `hbonds`         - per donor-H-acceptor triplet, the fraction of analysed frames in which the
                      geometric criterion holds. Hydrogen bonds form and break on a picosecond
                      scale while a DCD stores a frame every tens of picoseconds, so this is a
                      SAMPLE of the occupancy and is reported as such, never as a fraction of
                      elapsed time.
5. `rmsf`           - per pocket residue (CA), how much it moves about its own mean position.

Rules this module keeps, all of them learned from the real files of job 20070766:

- HEAVY ATOMS COME FROM MASS, NEVER FROM NAMES. The real ligand carries `LP1`, a CHARMM lone-pair
  particle of mass 0 on the chlorine. `not name H*` counts it and reports 34 heavy atoms where
  the docking record says 33.
- MINIMIZATION FRAMES ARE NOT DYNAMICS. `step4_equilibration.inp` holds `minimize 10000`, and with
  `dcdfreq 5000` that is 2 frames written before any dynamics. They carry a step number, so a
  reader takes them for the first 20 ps of simulation. The count is read from the input file and
  those frames are dropped.
- TRAJECTORIES WITH DIFFERENT FRAME SPACING ARE NOT CONCATENATED. Changing `dcdfreq` between
  chunks is a normal thing to do, and a chained read would silently stamp every chunk with the
  first one's spacing. Mixed spacing is refused by name.
- THE POCKET IS MEASURED, NOT NAMED. Its residues come from a declared distance cutoff around the
  ligand in the reference frame, and the full list is returned so it can be audited. A dossier
  label such as `3OG7_032_C` names a POCKET, not a segment of the built system.
- NO VERDICT. No function here decides whether a pose is stable. `reading()` offers one, marked as
  offered, and says the researcher's own judgement decides.
"""

from __future__ import annotations

import logging
import math
import re
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------------------------
# Criteria. Every one of these is a CHOSEN convention, not a measurement, and each is repeated in
# the output next to the number it produced so that a reader never has to come here to find it.
# ---------------------------------------------------------------------------------------------

HEAVY_MASS_MIN = 1.5
"""Atomic mass (u) above which an atom counts as heavy. Excludes hydrogen (~1.008) and the
massless virtual particles CHARMM uses for lone pairs."""

POCKET_CUTOFF = 5.0
"""Angstrom. A protein residue belongs to the pocket when any of its atoms is within this of a
ligand heavy atom in the reference frame. 5 A is the usual first-shell contact distance; it is a
convention, and `pocket_cutoff` changes it."""

HB_DA_CUTOFF = 3.5
"""Angstrom, donor to acceptor. The common geometric definition of a hydrogen bond."""

HB_ANGLE_CUTOFF = 150.0
"""Degrees, donor-hydrogen-acceptor. Below this the geometry is not counted as a bond."""

HB_NO_BY_MASS = "(mass 14.0 to 14.1 or mass 15.9 to 16.1)"
"""Nitrogen and oxygen, picked by mass for the same reason heavy atoms are: the real ligand's
atom types are CHARMM types (`OG2P1`, `NG2R51`), not elements, and this Universe carries no
element attribute at all, so a type- or element-based selection silently matches nothing."""

HB_H_BY_MASS = "mass 1.0 to 1.1"
"""Hydrogen by mass."""

PLATEAU_WINDOW = 0.25
"""The last quarter of the analysed frames is what `reading()` compares against the whole run when
it talks about a plateau. A quarter is a choice, nothing more."""

PLATEAU_MIN_FRAMES = 8
"""Fewest analysed frames for which the drift between the first and the last quarter is reported as
meaning anything. With `window = round(n/4)`, 8 frames is where each quarter first holds more than
one, so the figure is a difference between two AVERAGES rather than between two single frames. It
is a choice made here and is reported as one; below it the series carries `drift_warning`."""

CRITERIA = {
    "heavy_atoms": f"mass > {HEAVY_MASS_MIN} u (excludes hydrogen and massless lone pairs)",
    "pocket": f"protein residues with any atom within {POCKET_CUTOFF} A of a ligand heavy atom, "
              "in the reference frame",
    "hbond": f"acceptor = N or O by mass; donor = the N or O that the hydrogen is BONDED to in "
             f"the topology; geometry = donor-acceptor distance <= {HB_DA_CUTOFF} A and "
             f"donor-H-acceptor angle >= {HB_ANGLE_CUTOFF} deg, counted per frame. No partial "
             "charge takes part in deciding which atoms are eligible",
    "superposition": "the pocket backbone (N, CA, C, O of the pocket residues), mass-weighted; "
                     "never the ligand itself",
    "com_distance": "the ligand is first placed in the periodic image nearest its position in the "
                    "PREVIOUS analysed frame, using the frame's own box; the reference frame "
                    "itself is anchored on the pocket",
}


class MdAnalyzeError(Exception):
    """Raised with a message meant for the analyst, not a stack trace."""

    def __init__(self, message: str, **extra):
        super().__init__(message)
        self.extra = extra


# ---------------------------------------------------------------------------------------------
# 1. Finding the files
# ---------------------------------------------------------------------------------------------

def find_inputs(folder, topology=None, stage=None) -> dict:
    """The PSF, the DCDs and the NAMD inputs under `folder`, wherever build_md_job left them.

    `topology` overrides the search for a PSF, for the case where the topology lives somewhere
    else or is in another format MDAnalysis reads. It is the caller's choice and is reported as
    such: nothing is inferred from it.

    Looks in the folder itself and in `system/namd/`, which is where a build_md_job package keeps
    them, so either path works. DCDs are ordered by the chunk number in their name
    (`step5_2.dcd` after `step5_1.dcd`), because a plain alphabetical sort puts 10 before 2 and
    the time axis would come out shuffled.

    `stage` picks which CHARMM-GUI stage to analyse (`step4`, `step5`, ...). Left out, the LATEST
    stage present is used, and the others are named in `left_out`. Stages are never mixed: they
    were simulated at different points of the protocol, usually at different `dcdfreq`.
    """
    root = Path(folder)
    if not root.is_dir():
        raise MdAnalyzeError(f"{root} is not a folder.")
    places = [root, root / "system" / "namd", root / "run" / "system" / "namd", root / "namd"]
    seen, psfs, dcds, inps = set(), [], [], []
    files_seen = set()
    for place in places:
        if not place.is_dir() or place.resolve() in seen:
            continue
        seen.add(place.resolve())
        for target, pattern in ((psfs, "*.psf"), (dcds, "*.dcd"), (inps, "*.inp")):
            for path in sorted(place.glob(pattern)):
                # The same file can be reached through two of the places above: a Windows
                # junction, a symlink, a HARD link, or a `folder` that already IS the system/namd
                # directory. A DCD chained twice doubles the run and every mean, plateau window
                # and occupancy computed from it, and nothing in the output would show it.
                real = _file_identity(path)
                if real in files_seen:
                    continue
                files_seen.add(real)
                target.append(path)
    given = Path(topology) if topology else None
    if given is not None and not given.is_file():
        raise MdAnalyzeError(f"The topology you named is not a file: {given}")
    if not psfs and given is None:
        raise MdAnalyzeError(
            "No .psf found. The topology is what says which atom is which, so nothing can be "
            f"measured without it. Looked in: {', '.join(str(p) for p in places)}")
    if not dcds:
        raise MdAnalyzeError(
            "No .dcd found: there is no trajectory to measure yet. If the job is still running, "
            "the first chunk has not finished writing. "
            f"Looked in: {', '.join(str(p) for p in places)}")
    # step3_input.psf is the one CHARMM-GUI hands to NAMD; prefer it over the earlier stages
    # (step1_pdbreader.psf, step2.2_ions.psf), which describe a system that was not simulated.
    psf = given or next((p for p in psfs if "step3" in p.name), psfs[0])
    # Trajectories are grouped by the STAGE they belong to, and one stage is analysed.
    #
    # The earlier rule asked whether a name ended in `_<number>`, which is what build_md_job's own
    # chunks look like (step5_1.dcd ... step5_10.dcd), and called everything else "other". That is
    # a rule about a FILE NAME pretending to be a fact about a simulation, and it was wrong the
    # first time it met a trajectory this project did not generate: a 10 ns production run written
    # in one piece as `step5_production.dcd` was reported with `stage: other` and a note saying
    # "these are not production", about production data.
    #
    # Now: the stage is the text before the first `_` or `.` (`step5`), which is CHARMM-GUI's own
    # numbering, and the stage is NAMED BY ITS OWN NAMD INPUT. The words "production" and
    # "equilibration" never come from a rule here - they come from the input file's name, which
    # the build wrote.
    by_stage = {}
    for d in dcds:
        by_stage.setdefault(_stage_of(d), []).append(d)
    if stage:
        chosen = str(stage)
        if chosen not in by_stage:
            raise MdAnalyzeError(
                f"No trajectory belongs to stage {chosen!r}. The stages present are: "
                + "; ".join(f"{s} ({', '.join(d.name for d in v)})"
                            for s, v in sorted(by_stage.items())),
                stages_present={s: [d.name for d in v] for s, v in by_stage.items()})
    else:
        # The LATEST stage, because CHARMM-GUI numbers them in protocol order and the last one is
        # the one an analyst means. A declared rule, not an inference about the science.
        chosen = max(by_stage, key=lambda s: (_stage_rank(s), s))
    analysed = sorted(by_stage[chosen], key=_chunk_order)
    others = [d for d in dcds if d not in analysed]
    stage_inputs = [p for p in inps if _stage_of(p) == chosen]
    return {
        "psf": psf,
        # The names are NOT enough: a build_md_job package holds `step3_input.psf` both at the root
        # and under `run/system/namd/`, so "1 of 2 .psf found (step3_input.psf; the others:
        # step3_input.psf)" named neither. Regenerate CHARMM-GUI, update one copy only, and that
        # line would hide which system was measured. Paths, relative to the folder asked for.
        "psf_selection": (f"{_relative(psf, root)}, which you named; the search found "
                          f"{len(psfs)} .psf in the folder" if given is not None else
                          f"1 of {len(psfs)} .psf found ({_relative(psf, root)}; the others: "
                          + (", ".join(_relative(p, root) for p in psfs if p != psf) or "none")
                          + ")"),
        "dcds": analysed,
        "dcds_all": sorted(dcds, key=_chunk_order),
        "stage": chosen,
        "stage_inputs": [p.name for p in stage_inputs],
        "stage_named_by": ([p.stem for p in stage_inputs] or None),
        "stage_selection": (f"stage {chosen!r}, which you named" if stage else
                            f"stage {chosen!r}, the latest of "
                            f"{len(by_stage)} present ({', '.join(sorted(by_stage))})"
                            if len(by_stage) > 1 else
                            f"stage {chosen!r}, the only one present"),
        "stages_present": {s: [d.name for d in sorted(v, key=_chunk_order)]
                           for s, v in sorted(by_stage.items())},
        "stage_note": (
            f"{len(analysed)} of {len(dcds)} .dcd analysed: every trajectory of stage {chosen!r} "
            f"({', '.join(d.name for d in analysed)}). "
            + (f"Left out, as other stages of the same system: "
               f"{', '.join(f'{d.name} (stage {_stage_of(d)})' for d in others)}. They were "
               "simulated at another point of the protocol, usually at another dcdfreq, so "
               "chaining them here would run the time axis backwards. Pass `stage` to analyse one "
               "of them instead. " if others else "")
            + (f"The NAMD input of this stage is {', '.join(p.name for p in stage_inputs)}: any "
               "word like 'production' or 'equilibration' comes from THAT FILE'S NAME, written by "
               "whoever built the job, and not from any rule of this tool. Read "
               "`restraints` and `minimization` for what the input actually declares"
               if stage_inputs else
               f"NO .inp of stage {chosen!r} was found, so nothing here says which part of the "
               "protocol these trajectories are. Do not call them production or equilibration "
               "without the input file")),
        "left_out": [d.name for d in others],
        "psfs_all": psfs,
        "inputs": inps,
        "looked_in": [str(p) for p in places if p.is_dir()],
    }


def _file_identity(path: Path):
    """What makes two paths the same file: the device and inode, not the spelling.

    `resolve()` alone follows symlinks and junctions but NOT hard links, and `mklink /H` needs no
    privilege on Windows, so a hard-linked DCD would have been chained onto itself. The device and
    file index identify it in either case. Falls back to the resolved path when stat fails.
    """
    try:
        info = path.stat()
        if info.st_ino:
            return (info.st_dev, info.st_ino)
    except OSError:
        pass
    return path.resolve()


def _relative(path: Path, root: Path) -> str:
    """`path` as seen from `root`, so two files of the same name are told apart."""
    try:
        return str(path.relative_to(root)).replace("\\", "/")
    except ValueError:
        return str(path)


def _chunk_order(path: Path):
    """Sort key that reads the trailing number of `step5_10.dcd` as 10, not as the text "10".

    A file with no number sorts FIRST, not last: the one that occurs in practice is
    `step4_equilibration.dcd`, which was simulated before every chunk. Sorting it last put the
    time axis in reverse.
    """
    match = re.search(r"(\d+)(?=\.dcd$)", path.name)
    return (0, int(match.group(1))) if match else (-1, 0), path.name


# ---------------------------------------------------------------------------------------------
# 2. Frames that are not dynamics
# ---------------------------------------------------------------------------------------------

def _stage_of(path) -> str:
    """The CHARMM-GUI stage a file belongs to: the text before its first `_` or `.`.

    `step5_production.dcd`, `step5_1.dcd` and `step5_production.inp` are all stage `step5`. This
    is the one thing about these names that is a fact rather than a convention of this project.
    """
    return re.split(r"[_.]", Path(path).stem, maxsplit=1)[0]


def _stage_rank(stage: str) -> int:
    """The number inside a stage name (`step5` -> 5), or -1 when it has none.

    CHARMM-GUI numbers its stages in protocol order, so this orders them in time.
    """
    match = re.search(r"(\d+)", stage)
    return int(match.group(1)) if match else -1


def _stage_inputs(inputs, dcds=()):
    """The NAMD inputs that belong to the same stage as the trajectories being analysed.

    `step5_1.dcd` -> every `step5*.inp`. With no trajectories named, every input is a candidate.
    Returns `(candidates, stages)`.
    """
    stages = {s for s in (_stage_of(d) for d in map(str, dcds)) if s}
    candidates = [p for p in inputs if not stages or _stage_of(p) in stages]
    return candidates, stages


# The NAMD keywords that hold atoms in place, and the file and column each one reads. `constraints`
# ties atoms to their starting coordinates with a spring; `fixedAtoms` freezes them outright. Both
# were being described in the output as though both had been looked for, while only the first was.
_HOLD_KEYWORDS = {
    "constraints": ("consref", "conskfile", "conskcol", "spring (harmonic restraint)"),
    "fixedatoms": ("fixedatomsfile", "fixedatomsfile", "fixedatomscol", "frozen (fixed atoms)"),
}

# Where a force constant lives in a PDB line, by the letter `conskcol`/`fixedAtomsCol` names.
# NAMD reads the constant from one column of the reference PDB, and which one is a choice the
# input file makes: CHARMM-GUI writes B.
_PDB_COLUMNS = {"X": (30, 38), "Y": (38, 46), "Z": (46, 54), "O": (54, 60), "B": (60, 66)}


def _held_atom_flags(path: Path, column: str):
    """Per-atom flags and LABELS from a NAMD reference PDB, in file order.

    NAMD requires the reference file to carry the same atoms in the same order as the structure, so
    the position in this file IS the atom index. That is what makes it possible to say which of the
    atoms this module measures were held, instead of repeating the file name and hoping.

    Only `ATOM`/`HETATM` lines are counted, and that filter is not cosmetic: the real
    `prot_posres.ref` carries five `TER` lines, four of them INTERLEAVED, one immediately before
    the ligand's first atom. Counting them shifts every index after each one, and the report then
    says "1659 of 2174 backbone atoms held" - a credible number, with its "N of total", invented by
    an offset of five lines.

    The labels come back too, because a file of the right LENGTH can still describe another system:
    a posres from a regenerated CHARMM-GUI build has the same atom count and a different order, and
    without a check the answer is a confident "none of the 33 ligand atoms is held" about a ligand
    on a spring. `held_selection` compares them.

    Returns `(flags, labels, note)`; `flags` is None when nothing usable was read.
    """
    span = _PDB_COLUMNS.get((column or "B").strip().upper()[:1])
    if span is None:
        return None, None, (f"the force-constant column {column!r} is not one of X, Y, Z, O or B, "
                            "so the file was not read")
    try:
        text = path.read_text(errors="replace")
    except OSError as exc:
        return None, None, f"{path.name} could not be read ({type(exc).__name__})"
    flags, labels = [], []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chunk = line[span[0]:span[1]].strip()
        try:
            # abs() because a negative force constant is not a thing NAMD means by "free"; it is a
            # malformed file, and treating it as unheld would be the dangerous reading of the two.
            flags.append(abs(float(chunk)) > 0.0 if chunk else False)
        except ValueError:
            flags.append(False)
        labels.append((line[12:16].strip(), line[17:21].strip()))
    if not flags:
        return None, None, f"{path.name} holds no ATOM or HETATM line"
    return flags, labels, (f"{sum(flags)} of {len(flags)} atoms in {path.name} carry a non-zero "
                           f"force constant in column {span[0] + 1}-{span[1]} ({column})")


def restraints(inputs, dcds=()) -> dict:
    """Whether the NAMD input of this stage holds atoms in place, and WHICH atoms.

    This matters more than it looks, and more than naming the file would show. `step4_equilibration
    .inp` as CHARMM-GUI writes it carries `constraints on` with `consref restraints/prot_posres
    .ref`, and that file turned out to hold down not only the protein but **the ligand's 33 heavy
    atoms** at k = 1 kcal/mol/A^2 - exactly the 33 atoms this module measures. On such a trajectory
    neither the backbone RMSD nor the ligand RMSD nor the centre-of-mass range can move much, so
    all three controls pass by construction. An earlier version of this function named the file and
    marked only the protein line, which left the two ligand lines offered without reserve on a
    trajectory where the ligand is on a spring.

    So the reference file is READ, and the flags come back per atom in file order (NAMD requires
    the same atoms in the same order as the structure), which is what lets `_analyze` say how many
    of the atoms it measured were held rather than guessing from the file's name.

    `fixedAtoms` is checked too. The negative note used to claim it had been looked for while only
    `constraints` was, which is a sourced-looking statement about a check that never ran.
    """
    candidates, stages = _stage_inputs(inputs, dcds)
    evidence, found = [], []
    for path in candidates:
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        for keyword, (ref_key, k_key, col_key, description) in _HOLD_KEYWORDS.items():
            if not re.search(rf"^\s*{keyword}\s+on\b", text, re.MULTILINE | re.IGNORECASE):
                continue
            # `conskfile` wins over `consref` when both are given: NAMD takes the force constants
            # from conskfile and only the coordinates from consref, and it is the constants that
            # decide which atom is held.
            got_ref = re.search(rf"^\s*{k_key}\s+(\S+)", text, re.MULTILINE | re.IGNORECASE) \
                or re.search(rf"^\s*{ref_key}\s+(\S+)", text, re.MULTILINE | re.IGNORECASE)
            got_col = re.search(rf"^\s*{col_key}\s+(\S+)", text, re.MULTILINE | re.IGNORECASE)
            found.append({
                "input": path.name,
                "keyword": keyword,
                "kind": description,
                # the trailing `;` is CHARMM-GUI's own style on many lines (`rigidBonds all;`), so
                # a value read from one has to be stripped or the file name carries a semicolon
                "reference_file": got_ref.group(1).strip(";") if got_ref else None,
                "column": got_col.group(1).strip(";") if got_col else None,
                "column_source": (f"read from `{col_key}`" if got_col else
                                  "ASSUMED to be B, which is what CHARMM-GUI writes. The input "
                                  "declares no column, and NAMD's own default was NOT verified "
                                  "here, so the atom counts below are unverified"),
                "path": (path.parent / got_ref.group(1).strip(";")) if got_ref else None,
            })
        if found:
            for line in text.splitlines():
                if re.match(r"^\s*(constraints|consref|conskfile|conskcol|constraintScaling"
                            r"|fixedAtoms|fixedAtomsFile|fixedAtomsCol|consexp)\b",
                            line, re.IGNORECASE):
                    evidence.append(f"{path.name}: {line.strip()}")
            break

    on = bool(found)
    # When an input declares BOTH a spring and fixed atoms, the FROZEN one is reported as the
    # headline: freezing is the stronger statement, and reporting only the spring described the
    # wrong set of atoms while `evidence` quietly carried the other.
    primary = next((f for f in found if f["keyword"] == "fixedatoms"), found[0] if found else None)
    flags, labels, flags_note = (None, None, None)
    if primary is not None and primary["path"] is not None:
        flags, labels, flags_note = _held_atom_flags(primary["path"],
                                                     primary["column"] or "B")
    out = {
        "atoms_held": on,
        "positional_restraints": on,
        "kind": primary["kind"] if primary else None,
        "read_from": primary["input"] if primary else None,
        "reference_file": primary["reference_file"] if primary else None,
        "force_constant_column": (primary["column"] or "B") if primary else None,
        "force_constant_column_source": primary["column_source"] if primary else None,
        "declarations": [{k: v for k, v in f.items() if k != "path"} for f in found],
        "declarations_note": (
            f"{len(found)} kinds of hold are declared ({', '.join(f['kind'] for f in found)}). The "
            "atom counts in `held` come from the FROZEN one, which is the stronger statement; the "
            "other one's reference file was not read" if len(found) > 1 else
            f"one kind of hold is declared ({primary['kind']})" if primary else
            "no hold is declared"),
        "evidence": evidence,
        "inputs_read": [p.name for p in candidates],
        "stages_of_the_analysed_trajectories": sorted(stages),
        "reference_file_note": (
            flags_note if flags_note else
            f"the reference file {primary['reference_file']!r} named by {primary['input']} could "
            "not be found next to it, so WHICH atoms were held is unknown"
            if primary is not None and primary["path"] is not None else
            f"{primary['input']} turns the restraint on but names no reference file"
            if primary else None),
        "note": (
            f"{primary['input']} declares `{primary['keyword']}` on: the atoms its reference file "
            "marks are held in place. On such a trajectory an RMSD of those atoms CANNOT move "
            "much, so any control built on it passes by construction and is not evidence of "
            "anything. Which atoms were held is in `held`, read from the file itself"
            if on else
            "no NAMD input of this stage declares `constraints on` or `fixedAtoms on` (both were "
            "looked for), so nothing read here says the atoms were held down"
            if candidates else
            "NO NAMD input of this stage was found, so this is unknown rather than absent: if the "
            "trajectory came from a restrained equilibration, nothing here would show it"),
    }
    out["_flags"] = (flags, labels)
    return out


def held_selection(flags, labels, group, label: str, n_atoms: int) -> dict:
    """How many atoms of `group` the restraint file holds, by atom index. No guessing.

    Before answering, the file is checked against the system TWICE, because a wrong answer here is
    worse than no answer: `held: 0` reads as "the ligand was free" and removes the warning from
    the two lines that most need it.

    1. The file must describe as many atoms as the system has. NAMD aborts otherwise, so a
       mismatch means the file is not the one that ran.
    2. The atom and residue names at this group's own indices must agree with the topology. A
       posres file from a regenerated build has the right length and another order, and only the
       names catch that.

    Either check failing makes the answer UNKNOWN, loudly, instead of a confident zero.
    """
    if flags is None:
        return {"known": False, "note": f"whether the {label} was held is UNKNOWN: no restraint "
                                        "reference file was read"}
    if len(flags) != n_atoms:
        return {"known": False, "atoms_in_the_reference_file": len(flags), "atoms_in_the_system":
                int(n_atoms),
                "note": f"whether the {label} was held is UNKNOWN: the reference file describes "
                        f"{len(flags)} atoms and this system has {n_atoms}, so it is NOT the file "
                        "this trajectory ran with and no index in it means anything here"}
    mismatched = []
    if labels is not None:
        for index in group.indices:
            atom = group.universe.atoms[int(index)]
            name, resname = labels[int(index)]
            if name != atom.name or resname != atom.resname:
                mismatched.append(f"index {int(index)}: the file says {resname} {name}, the "
                                  f"topology says {atom.resname} {atom.name}")
            if len(mismatched) >= 3:
                break
    if mismatched:
        return {"known": False, "mismatches": mismatched,
                "note": f"whether the {label} was held is UNKNOWN: the reference file has the "
                        f"right number of atoms but they are NOT in this system's order "
                        f"({'; '.join(mismatched)}), so it describes another build"}
    held = [int(i) for i in group.indices if flags[int(i)]]
    return {
        "known": True,
        "held": len(held),
        "of": int(group.n_atoms),
        "note": (f"{len(held)} of {group.n_atoms} {label} atoms are held in place by the restraint "
                 f"file, so their motion is damped and NOT free. Any figure built on those atoms "
                 f"passes by construction" if held else
                 f"none of the {group.n_atoms} {label} atoms is held by the restraint file, whose "
                 f"atom and residue names were checked against this topology at each of their "
                 f"indices"),
        "alignment_checked": labels is not None,
    }


def minimization_frames(inputs, dcds=()) -> dict:
    """How many leading DCD frames a NAMD input writes from `minimize`, and the evidence for it.

    `dcds` are the trajectories actually being analysed, and only the inputs that belong to the
    SAME stage are read. Without that, a folder holding both `step4_equilibration.inp` (which has
    `minimize 10000`) and `step5_production.inp` (which has none) made every production chunk lose
    its first real frames to a minimization that chunk never ran - and the report then said so in
    writing. The real package folder is exactly that folder.

    NAMD counts minimization steps with the same counter as dynamics and writes them to the same
    DCD, so they arrive with a step number and a reader takes them for simulated time. They are
    not: minimization has no time axis. The count is `minimize / dcdfreq`, both read from the
    input file, and the lines are returned so the number can be checked.

    Returns `{"frames": 0, ...}` when no input declares a minimization, which is the normal case
    for a production chunk.
    """
    candidates, stages = _stage_inputs(inputs, dcds)
    out = {"frames": 0, "evidence": [], "read_from": None,
           "stages_of_the_analysed_trajectories": sorted(stages),
           "inputs_read": [p.name for p in candidates],
           "inputs_ignored": [p.name for p in inputs if p not in candidates],
           "matching_rule": "only the NAMD inputs whose name starts with the same stage as the "
                            "analysed trajectories are read (step5_1.dcd -> step5*.inp)"}
    for path in candidates:
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        # NAMD keywords are CASE-INSENSITIVE. Without IGNORECASE a file writing `DCDfreq` or
        # `Minimize` - both legal, and NAMD obeys them - read as "this input declares no
        # minimization", and the report then said so in writing while two minimization frames
        # entered the time axis as 20 ps of dynamics.
        mini = re.search(r"^\s*minimize\s+(\d+)", text, re.MULTILINE | re.IGNORECASE)
        freq = re.search(r"^\s*dcdfreq\s+(\d+)", text, re.MULTILINE | re.IGNORECASE)
        if not mini or not freq or int(freq.group(1)) == 0:
            continue
        frames = int(mini.group(1)) // int(freq.group(1))
        if frames > out["frames"]:
            out.update(frames=frames, read_from=path.name, evidence=[
                f"{path.name}: minimize {mini.group(1)}",
                f"{path.name}: dcdfreq {freq.group(1)}",
                f"{mini.group(1)} / {freq.group(1)} = {frames} frame(s) written before any dynamics",
            ])
    return out


# ---------------------------------------------------------------------------------------------
# 3. The universe, and the refusal to mix frame spacings
# ---------------------------------------------------------------------------------------------

def trajectory_plan(psf, dcds) -> dict:
    """Read each DCD on its own to learn its frame count and spacing, and refuse a mixed set.

    MDAnalysis will happily chain DCDs of different `dcdfreq` and stamp them all with the first
    file's spacing, which turns a time axis into fiction. Changing `dcdfreq` between chunks is a
    reasonable thing for an analyst to do, so this is checked rather than assumed.
    """
    import MDAnalysis as mda

    per_file = []
    for dcd in dcds:
        try:
            u = mda.Universe(str(psf), str(dcd))
        except Exception as exc:                                  # noqa: BLE001 - reported, not raised
            per_file.append({"file": dcd.name, "ok": False, "error": f"{type(exc).__name__}: {exc}"})
            continue
        try:
            per_file.append({
                "file": dcd.name,
                "ok": True,
                "frames": int(len(u.trajectory)),
                "dt_ps": round(float(u.trajectory.dt), 6),
                "atoms": int(len(u.atoms)),
            })
        finally:
            # On Windows an open DCD keeps a lock on the file, so every reader is closed even
            # though this function only wanted to look at the header.
            u.trajectory.close()
    good = [f for f in per_file if f.get("ok")]
    bad = [f for f in per_file if not f.get("ok")]
    if not good:
        raise MdAnalyzeError("No .dcd could be read.", per_file=per_file)
    if bad:
        # Reading the rest and carrying on would be worse than refusing: a chunk that cannot be
        # read is a HOLE IN THE TIME AXIS, and every series, plateau window and occupancy after it
        # would be computed over a trajectory with a silent gap. The earlier code noted the
        # failure here and then chained the broken file anyway, which crashed with a raw OSError.
        raise MdAnalyzeError(
            "These trajectories could not be read, so the time axis would have a gap in it: "
            + "; ".join(f"{f['file']} ({f['error']})" for f in bad)
            + ". Move them aside to measure the rest, and say in the report which chunks were "
              "left out.",
            per_file=per_file, unreadable=[f["file"] for f in bad])
    spacings = sorted({f["dt_ps"] for f in good})
    if len(spacings) > 1:
        groups = {s: [f["file"] for f in good if f["dt_ps"] == s] for s in spacings}
        raise MdAnalyzeError(
            "These trajectories do not share one frame spacing, so they cannot be read as one "
            "time axis: " + "; ".join(f"{s} ps -> {', '.join(n)}" for s, n in groups.items()) +
            ". This is what happens when dcdfreq changes between chunks. Analyse each spacing on "
            "its own (move the odd files aside, or pass the folder that holds only one set).",
            per_file=per_file, spacings=spacings)
    return {"per_file": per_file, "dt_ps": spacings[0],
            "frames_total": sum(f["frames"] for f in good),
            "selection": f"{len(good)} of {len(dcds)} .dcd read"}


# ---------------------------------------------------------------------------------------------
# 4. Selections
# ---------------------------------------------------------------------------------------------

def select_ligand(universe, resname: str, segid: str = "") -> dict:
    """The ligand's heavy atoms by MASS, with the count of what was dropped and why.

    Refuses when the system holds MORE THAN ONE copy of the residue. Two copies selected together
    are not a molecule: their joint centre of mass sits between them, in a place where there is no
    ligand at all, the RMSD mixes two molecules and the pocket comes out as the union of two
    sites - and `N heavy of M atoms` reads as one ligand the whole time. Several structures in
    this project's own dossier have two copies, so this is not hypothetical.
    """
    whole = universe.select_atoms(f"resname {resname}")
    if not len(whole):
        all_odd = sorted({r for r in set(universe.atoms.resnames) if r not in _STANDARD_RESNAMES})
        shown = all_odd[:25]
        raise MdAnalyzeError(
            f"No atom has resname {resname!r}. Non-standard residue names in this system "
            f"({len(shown)} of {len(all_odd)} shown): {', '.join(shown) or 'none'}. "
            "Pass `ligand_resname` with one of them.",
            resnames_available=all_odd)
    copies = [f"{r.segid}:{r.resid}" for r in whole.residues]
    if segid:
        whole = whole.select_atoms(f"segid {segid}")
        if not len(whole):
            raise MdAnalyzeError(
                f"No atom has resname {resname!r} in segment {segid!r}. The copies are: "
                f"{', '.join(copies)}.", copies=copies)
        copies_used = [f"{r.segid}:{r.resid}" for r in whole.residues]
        if len(copies_used) > 1:
            raise MdAnalyzeError(
                f"Segment {segid!r} still holds {len(copies_used)} copies of {resname!r}: "
                f"{', '.join(copies_used)}. This tool measures one ligand.", copies=copies_used)
    elif len(copies) > 1:
        raise MdAnalyzeError(
            f"This system holds {len(copies)} copies of resname {resname!r} ({', '.join(copies)}) "
            "and they are NOT one molecule: measured together, their centre of mass falls between "
            "them, where no ligand is. Pass `ligand_segid` with the segment of the copy to "
            "measure, and run the others separately.",
            copies=copies)
    heavy = whole.select_atoms(f"mass {HEAVY_MASS_MIN} to 100000")
    dropped = sorted(set(whole.names) - set(heavy.names))
    if not len(heavy):
        raise MdAnalyzeError(f"resname {resname!r} has no atom heavier than {HEAVY_MASS_MIN} u.")
    return {
        "atoms": heavy,
        "resname": resname,
        "selection": f"{len(heavy)} heavy of {len(whole)} atoms in resname {resname}"
                     + (f", segment {segid}" if segid else ""),
        "criterion": CRITERIA["heavy_atoms"],
        "copies_in_system": copies,
        "copies_note": (f"1 of {len(copies)} copies measured ({', '.join(copies)}); the others "
                        "were not and nothing here describes them"
                        if len(copies) > 1 else
                        f"one copy in the system ({copies[0]})"),
        "dropped": dropped,
        "dropped_note": (f"dropped by mass: {', '.join(dropped)}" if dropped else "nothing dropped"),
        "total_mass_u": round(float(whole.masses.sum()), 3),
        "segids": sorted(set(whole.segids)),
    }


_STANDARD_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HSD", "HSE", "HSP", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "TIP3", "TIP4", "WAT", "HOH", "SOD", "CLA", "POT", "CAL", "MG", "ZN",
}


def select_pocket(universe, ligand_heavy, cutoff: float, ref_frame: int) -> dict:
    """The protein residues around the ligand in the reference frame, with the list returned.

    Measured, not named: the caller gets every residue back as `SEGID:RESNAME RESID` so the
    selection can be checked against what is known about the site. The segments the pocket spans
    are reported too, because a ligand at a chain interface is a different situation from one
    buried in a single chain, and nothing here decides which it is.
    """
    universe.trajectory[ref_frame]
    pocket = universe.select_atoms(f"protein and byres around {cutoff} group L",
                                   L=ligand_heavy, updating=False)
    if not pocket.n_atoms:
        raise MdAnalyzeError(
            f"No protein residue lies within {cutoff} A of the ligand in frame {ref_frame}. "
            "Either the ligand is not in a site, or the reference frame is wrong.")
    residues = pocket.residues
    protein = universe.select_atoms("protein")
    spread = {}
    for r in residues:
        spread[r.segid] = spread.get(r.segid, 0) + 1
    backbone = pocket.select_atoms("name N CA C O")
    return {
        "atoms": pocket,
        "backbone": backbone,
        "residues": [f"{r.segid}:{r.resname}{r.resid}" for r in residues],
        "selection": f"{residues.n_residues} of {protein.residues.n_residues} protein residues",
        "criterion": CRITERIA["pocket"],
        "cutoff_angstrom": cutoff,
        "reference_frame": int(ref_frame),
        "segments": spread,
        "segments_note": (f"the pocket lies in one segment ({next(iter(spread))})" if len(spread) == 1
                          else f"the pocket spans {len(spread)} segments: " +
                               ", ".join(f"{k} ({v} residues)" for k, v in sorted(spread.items()))),
        "backbone_atoms": int(backbone.n_atoms),
    }


# ---------------------------------------------------------------------------------------------
# 4b. Which residues were BUILT rather than observed
# ---------------------------------------------------------------------------------------------

# CHARMM names histidine by protonation state; the deposited file says HIS. Nothing else is
# renamed: any other disagreement means the file does not describe this segment's numbering.
_RESNAME_ALIASES = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS"}

MODELLED_RULE = (
    "a residue of the topology that is ABSENT from the `<pdb>_<segid>.orig.pdb` CHARMM-GUI writes "
    "for that segment, or present there without all four backbone atoms N, CA, C, O, is counted "
    "as built. That file is the chain as CHARMM-GUI read it, before it filled gaps; the rule is "
    "read from the file names CHARMM-GUI uses, it is not a statement CHARMM-GUI makes. 'Observed' "
    "means present in the structure that was UPLOADED to CHARMM-GUI: it is experimental only if "
    "that structure was a deposited one, and this tool does not check what was uploaded. Only the "
    "folder of the topology and the one above it are searched. The residue NAMES of every residue "
    "the two share are compared (HSD/HSE/HSP read as HIS), and any disagreement makes the answer "
    "UNKNOWN, never 'none'")

_BACKBONE = {"N", "CA", "C", "O"}


def _ranges(resids) -> str:
    """[545, 546, 547, 597, 598] -> '545-547, 597-598'."""
    out, start, prev = [], None, None
    for r in sorted(resids):
        if start is None:
            start = prev = r
        elif r == prev + 1:
            prev = r
        else:
            out.append(f"{start}-{prev}" if prev != start else f"{start}")
            start = prev = r
    if start is not None:
        out.append(f"{start}-{prev}" if prev != start else f"{start}")
    return ", ".join(out)


def _orig_residues(path: Path):
    """({resid: resname}, {resid: atom names}) of an `.orig.pdb`, or a reason it cannot be trusted."""
    found, atoms = {}, {}
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 27:
            return None, f"{path.name} has an ATOM line too short to carry a residue number"
        if line[26].strip():
            return None, (f"{path.name} uses insertion codes ({line[17:27].strip()}), so its "
                          "residue numbers cannot be matched one to one with the topology")
        try:
            resid = int(line[22:26])
        except ValueError:
            return None, f"{path.name} has a residue number that is not an integer: {line[22:26]!r}"
        name = line[17:21].strip()
        if found.get(resid, name) != name:
            return None, f"{path.name} gives residue {resid} two names ({found[resid]}, {name})"
        found[resid] = name
        atoms.setdefault(resid, set()).add(line[12:16].strip())
    if not found:
        return None, f"{path.name} holds no ATOM or HETATM record"
    return (found, atoms), None


def modelled_residues(universe, search_dirs, root=None) -> dict:
    """Per protein segment, the residues CHARMM-GUI built instead of reading from the structure.

    Why this exists: on the first production run of this project the full-backbone RMSD of the
    pocket chain drifted 1.35 A, and 1.25 A of it came from 25 residues MODELLER had filled in
    (the activation loop 597-614 among them), which moved ~10 A while the 247 observed residues
    stayed under 0.5 A. A built loop is a guess with no experimental coordinates, and an RMSD over
    it says how much the guess relaxed, not whether the protein is equilibrated.

    Every segment answers `known: true` with the list (possibly empty), or `known: false` with
    the reason. An unknown is never reported as "nothing was built".
    """
    root = Path(root) if root is not None else None
    dirs, seen = [], set()
    for d in search_dirs:
        d = Path(d)
        if d.is_dir() and d.resolve() not in seen:
            seen.add(d.resolve())
            dirs.append(d)
    backbone = universe.select_atoms("protein and name N CA C O")
    segments = {}
    for seg in sorted(set(backbone.segids)):
        residues = universe.select_atoms(f"protein and segid {seg}").residues
        total = int(residues.n_residues)
        candidates, ids = [], set()
        for d in dirs:
            for path in sorted(d.glob(f"*_{seg.lower()}.orig.pdb")):
                ident = _file_identity(path)
                if ident not in ids:
                    ids.add(ident)
                    candidates.append(path)
        names = [_relative(p, root) if root else p.name for p in candidates]
        if not candidates:
            segments[seg] = {"known": False, "residues_total": total,
                             "reason": f"no *_{seg.lower()}.orig.pdb was found next to the "
                                       "topology, so which residues were built is UNKNOWN"}
            continue
        if len(candidates) > 1:
            segments[seg] = {"known": False, "residues_total": total,
                             "reason": f"{len(candidates)} different files could describe {seg} "
                                       f"({', '.join(names)}); none was chosen"}
            continue
        parsed, problem = _orig_residues(candidates[0])
        if problem:
            segments[seg] = {"known": False, "residues_total": total, "source_file": names[0],
                             "reason": problem}
            continue
        orig, orig_atoms = parsed
        topo = {}
        for r in residues:
            if int(r.resid) in topo:
                problem = f"the topology numbers two residues of {seg} as {r.resid}"
                break
            topo[int(r.resid)] = r.resname
        extra = sorted(set(orig) - set(topo))
        clash = [f"{k}: {orig[k]} in the file, {topo[k]} in the topology"
                 for k in sorted(set(orig) & set(topo))
                 if _RESNAME_ALIASES.get(topo[k], topo[k]) != _RESNAME_ALIASES.get(orig[k], orig[k])]
        if problem or extra or clash:
            reason = problem or (
                f"{len(extra)} residue number(s) of {names[0]} are not in segment {seg} "
                f"({_ranges(extra)})" if extra else
                f"{len(clash)} residue(s) carry a different name in {names[0]} than in the "
                f"topology ({'; '.join(clash[:5])}{'; ...' if len(clash) > 5 else ''})")
            segments[seg] = {"known": False, "residues_total": total, "source_file": names[0],
                             "reason": reason + ". The file does not describe this segment's "
                                                "numbering, so nothing is inferred from it"}
            continue
        absent = set(topo) - set(orig)
        # Present in the file but without a full backbone: MODELLER built the missing atoms, and
        # a backbone RMSD over them is partly a guess. They go with the built residues.
        partial = sorted(k for k in set(orig) & set(topo) if not _BACKBONE <= orig_atoms[k])
        built = sorted(absent | set(partial))
        segments[seg] = {
            "known": True,
            "source_file": names[0],
            "residues": built,
            "ranges": _ranges(built) or "none",
            "count": len(built),
            "residues_total": total,
            "partial_backbone": partial,
            "selection": f"{len(built)} of {total} residues of {seg} were built: {len(absent)} "
                         f"absent from {names[0]}"
                         + (f" and {len(partial)} present there without all of N, CA, C, O "
                            f"({_ranges(partial)})" if partial else ""),
        }
    return {"rule": MODELLED_RULE, "segments": segments}


# ---------------------------------------------------------------------------------------------
# 5. The single pass that produces four of the five metrics
# ---------------------------------------------------------------------------------------------

def _series(values, times, *, unit: str) -> dict:
    """A measured series with its own summary. No verdict, no smoothing, no fitting."""
    if not values:
        return {"values": [], "times_ps": [], "unit": unit, "n": 0}
    n = len(values)
    window = max(1, round(n * PLATEAU_WINDOW))
    tail = values[max(0, n - window):]
    head = values[:window]
    # Whether the drift means anything at all. An earlier version tested only whether the two
    # quarters OVERLAP, which with `window = round(n/4)` can happen at n = 1 and nowhere else, so
    # it waved through n = 2, where "the drift" is the difference between two consecutive frames
    # and "a plateau" is not a thing those two frames can show. The real condition is that each
    # window hold more than one frame, which needs at least PLATEAU_MIN_FRAMES.
    overlap = max(0, 2 * window - n)
    too_few = n < PLATEAU_MIN_FRAMES
    ordered = sorted(values)
    mid = n // 2
    return {
        "values": [round(float(v), 4) for v in values],
        "times_ps": [round(float(t), 3) for t in times],
        "unit": unit,
        "times_note": "picoseconds elapsed FROM THE FIRST ANALYSED FRAME, which is 0; not the "
                      "frame's own timestamp in the file",
        "n": n,
        "min": round(float(ordered[0]), 4),
        "max": round(float(ordered[-1]), 4),
        "range": round(float(ordered[-1] - ordered[0]), 4),
        "mean": round(float(sum(values) / n), 4),
        "median": round(float(ordered[mid] if n % 2 else (ordered[mid - 1] + ordered[mid]) / 2), 4),
        "first_quarter_mean": round(float(sum(head) / len(head)), 4),
        "last_quarter_mean": round(float(sum(tail) / len(tail)), 4),
        "drift": round(float(sum(tail) / len(tail) - sum(head) / len(head)), 4),
        "drift_is_meaningful": not too_few,
        "drift_warning": (
            f"THE DRIFT BELOW CANNOT SUPPORT A STATEMENT ABOUT A PLATEAU: with {n} frame(s) each "
            f"quarter holds {window} frame(s)"
            + (", and the first and the last quarter are THE SAME FRAME, so the drift is 0 by "
               "arithmetic and no test on it can fail. " if overlap else
               ", so the figure is a difference between single frames rather than between two "
               "averages. ")
            + f"At least {PLATEAU_MIN_FRAMES} analysed frames are needed for each quarter to hold "
              "more than one, and that figure is a choice made here, not a measurement"
            if too_few else None),
        "drift_is_meaningful_rule": f"at least {PLATEAU_MIN_FRAMES} analysed frames, so that each "
                                    f"quarter window holds more than one frame. A CHOSEN minimum",
        "drift_note": "last quarter mean minus FIRST quarter mean. Comparing the last quarter "
                      "with the mean of the whole run instead looks like a plateau test and is "
                      "not one: for a straight ramp of height h it only ever reaches 0.375 h, so "
                      "a backbone climbing steadily by 1.2 A would clear a 0.5 A test",
        "last_quarter_range": round(float(max(tail) - min(tail)), 4),
        "last_quarter_frames": len(tail),
        "last_quarter_note": f"the last {PLATEAU_WINDOW:.0%} of the analysed frames; "
                             "the window is a choice, not a measurement",
    }


def measure(universe, ligand, pocket, *, first_frame: int, ref_frame: int = None,
            modelled=None) -> dict:
    """One pass over the trajectory: backbone RMSD, ligand RMSD, COM distance and RMSF inputs.

    Every frame is superposed on the POCKET BACKBONE. The rotation that does it is then applied to
    the ligand, so the ligand RMSD answers "did it move relative to the site" instead of "did its
    internal geometry change", which is what fitting on the ligand itself would measure.

    `ref_frame` is the frame every RMSD is measured AGAINST; it defaults to `first_frame`, the
    first analysed frame. It is the same frame the pocket was selected in, so the reference pose
    and the reference site always come from one instant.

    Periodic images are undone BY CONTINUITY: each frame's ligand is placed in the image nearest
    to where the ligand was in the PREVIOUS frame, not nearest to the pocket. Anchoring on the
    pocket looks equivalent and is not: it cannot tell a real displacement of more than half a box
    from a wrapping jump, and it folds the real one back. Measured on a ligand genuinely walking
    0 -> 20 -> 40 -> 58 A out of a 60 A box, the pocket-anchored rule reported 0, 20, 20, 2 A and
    the reading said the pose had held. Continuity has no such ambiguity at these frame spacings:
    a ligand does not diffuse half a box in the tens of picoseconds between two stored frames,
    while a wrap moves it a whole box in one step.

    Only the small selections are kept in memory (33 ligand atoms, ~28 pocket CAs), so the cost
    does not grow with the 140k atoms of the box.

    `modelled` maps a segid to the residue numbers that were BUILT, not observed (see
    `modelled_residues`). For each segment with any, two more series are measured: the backbone
    of the OBSERVED residues superposed on itself, and the built residues alone, carried by that
    same superposition. The full-segment series stays the headline: it is the pre-registered
    control, and the split is an addition, not a replacement.
    """
    import numpy as np
    from MDAnalysis.analysis import align

    lig = ligand["atoms"]
    fit = pocket["backbone"]
    pocket_all = pocket["atoms"]
    # One CA per pocket residue, and the LABELS come from these atoms themselves. Zipping the RMSF
    # values against the pocket's residue list assumed every residue contributes exactly one CA;
    # a residue without one shifted every label by one from there on, and the "N pocket residues"
    # line stayed self-consistent while each number named the wrong residue.
    pocket_ca = pocket["atoms"].select_atoms("name CA")
    ca_labels = [f"{a.segid}:{a.resname}{a.resid}" for a in pocket_ca]
    missing_ca = [r for r in pocket["residues"] if r not in ca_labels]
    if len(ca_labels) != len(set(ca_labels)):
        raise MdAnalyzeError(
            "Some pocket residue has more than one atom named CA, so an RMSF value could not be "
            "tied to one residue.", pocket_residues=pocket["residues"])
    backbone_all = universe.select_atoms("protein and name N CA C O")
    if ref_frame is None:
        ref_frame = first_frame

    # ONE SERIES PER PROTEIN SEGMENT, never one for all of them at once. With `wrapAll on` NAMD
    # re-images whole molecules, and two chains are two molecules: when one is wrapped and the
    # other is not, a joint superposition reads tens of angstroms of RMSD while nothing has moved.
    # A segment superposed on ITSELF is immune to that, because subtracting its own centre of mass
    # removes any whole-molecule translation, including a box-length one.
    segids = sorted(set(backbone_all.segids))
    seg_sel = {s: backbone_all.select_atoms("segid " + s) for s in segids}
    seg_sel = {s: sel for s, sel in seg_sel.items() if sel.n_atoms}

    # Observed / built split, per segment with at least one built residue AND one observed one.
    split = {}
    for s, built in (modelled or {}).items():
        if s not in seg_sel or not built:
            continue
        sel = seg_sel[s]
        mask = np.isin(sel.resids, list(built))
        if mask.all() or not mask.any():
            continue
        split[s] = (sel[~mask], sel[mask])

    universe.trajectory[ref_frame]
    fit_masses = fit.masses
    ref_fit_com = fit.center_of_mass()
    ref_fit = fit.positions - ref_fit_com
    seg_ref = {s: (sel.positions - sel.center_of_mass(), sel.masses) for s, sel in seg_sel.items()}
    split_ref = {}
    for s, (obs, built) in split.items():
        com = obs.center_of_mass()
        split_ref[s] = (obs.positions - com, built.positions - com)
    obs_rmsd = {s: [] for s in split}
    built_rmsd = {s: [] for s in split}

    frames_without_box = 0
    non_orthorhombic = 0
    pbc_corrected = 0
    pocket_split = 0

    def _box_of_frame():
        """The frame's orthorhombic box lengths, or None with the reason counted."""
        nonlocal frames_without_box, non_orthorhombic
        dims = universe.trajectory.ts.dimensions
        if dims is None or not all(b > 0 for b in dims[:3]):
            frames_without_box += 1
            return None
        if not np.allclose(dims[3:], 90.0, atol=1e-3):
            non_orthorhombic += 1
            return None
        return np.asarray(dims[:3], dtype=float)

    # The reference frame is wrapped like any other, so its ligand is placed in the image nearest
    # the pocket. This is the ONE frame anchored on the pocket: it is the origin of the chain of
    # continuity that every later frame follows, and the ligand is in its site by definition here
    # (the pocket was selected around it in this very frame).
    ref_lig_pos = lig.positions
    _box = _box_of_frame()
    if _box is not None:
        ref_lig_pos = ref_lig_pos - _box * np.round(
            (ref_lig_pos.mean(axis=0) - pocket_all.center_of_mass()) / _box)
    ref_lig = ref_lig_pos - ref_fit_com
    ref_lig_mean = ref_lig_pos.mean(axis=0)
    # `_box_of_frame` has now counted the reference frame once; the per-frame loop counts it again
    # when the reference is inside the analysed range, so start the counters over.
    frames_without_box = 0
    non_orthorhombic = 0

    times, lig_rmsd, com_dist = [], [], []
    seg_rmsd = {s: [] for s in seg_sel}
    ca_frames = []
    prev_mean = None

    def unwrap(positions, pocket_com):
        """Place the ligand in the image nearest where it was in the PREVIOUS analysed frame.

        CHARMM-GUI writes `wrapAll on`, so NAMD moves whole molecules back into the central cell.
        A ligand sitting against a face of the box jumps a full box length between one frame and
        the next WITHOUT MOVING, and a raw RMSD reads that as the ligand leaving. The correction
        has to be applied to the COORDINATES, before the RMSD, so that the same numbers feed the
        RMSD and the centre-of-mass distance and the two can never contradict each other.

        The anchor is the previous frame, not the pocket. See the docstring of `measure`: the
        pocket anchor silently folds a real dissociation back into the site.
        """
        nonlocal pbc_corrected, prev_mean
        anchor = pocket_com if prev_mean is None else prev_mean
        box = _box_of_frame()
        if box is None:
            prev_mean = positions.mean(axis=0)
            return positions
        shift = box * np.round((positions.mean(axis=0) - anchor) / box)
        if np.any(shift):
            pbc_corrected += 1
        moved = positions - shift
        prev_mean = moved.mean(axis=0)
        return moved

    t0 = None
    for ts in universe.trajectory[first_frame:]:
        # Time is counted FROM THE FIRST ANALYSED FRAME, not from the file's own origin.
        # MDAnalysis spreads dt across every frame in the file, minimization frames included, so
        # `ts.time` put the last point of a 250 ps equilibration at 270 ps - an instant that never
        # existed, the exact length of the two minimization frames.
        t0 = float(ts.time) if t0 is None else t0
        times.append(round(float(ts.time) - t0, 6))

        # protein backbone against the reference, one segment at a time: the fit residual IS the
        # RMSD.
        for s, sel in seg_sel.items():
            ref_pos, masses = seg_ref[s]
            _, rmsd_seg = align.rotation_matrix(sel.positions - sel.center_of_mass(),
                                                ref_pos, weights=masses)
            seg_rmsd[s].append(float(rmsd_seg))
        # the same, on the observed residues only; the built ones ride on that superposition, so
        # their RMSD says how far the guess moved RELATIVE TO THE OBSERVED FRAME of the chain.
        for s, (obs, built) in split.items():
            obs_ref, built_ref = split_ref[s]
            obs_com = obs.center_of_mass()
            rot, rmsd_obs = align.rotation_matrix(obs.positions - obs_com, obs_ref,
                                                  weights=obs.masses)
            obs_rmsd[s].append(float(rmsd_obs))
            moved = (built.positions - obs_com) @ np.asarray(rot).T
            sq = ((moved - built_ref) ** 2).sum(axis=1)
            built_rmsd[s].append(float(np.sqrt((sq * built.masses).sum() / built.masses.sum())))

        # superpose on the pocket backbone, then carry the ligand along.
        pocket_com = pocket_all.center_of_mass()
        # A pocket spanning two chains can itself be torn apart by `wrapAll`, and then its centre
        # of mass and its superposition are both nonsense. Detected, never silently corrected.
        _pbox = universe.trajectory.ts.dimensions
        if _pbox is not None and all(b > 0 for b in _pbox[:3]):
            extent = pocket_all.positions.max(axis=0) - pocket_all.positions.min(axis=0)
            if np.any(extent > np.asarray(_pbox[:3], dtype=float) / 2.0):
                pocket_split += 1
        lig_pos = unwrap(lig.positions, pocket_com)
        fit_com = fit.center_of_mass()
        rot, _ = align.rotation_matrix(fit.positions - fit_com, ref_fit, weights=fit_masses)
        moved_lig = (lig_pos - fit_com) @ rot.T
        diff = moved_lig - ref_lig
        lig_rmsd.append(float(np.sqrt((diff ** 2).sum(axis=1).mean())))
        ca_frames.append(((pocket_ca.positions - fit_com) @ rot.T).copy())

        weighted = (lig_pos * lig.masses[:, None]).sum(axis=0) / lig.masses.sum()
        com_dist.append(float(np.linalg.norm(weighted - pocket_com)))

    ca = np.asarray(ca_frames)
    rmsf = np.sqrt(((ca - ca.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0))

    per_segment = {s: {**_series(vals, times, unit="angstrom"),
                       "atoms": int(seg_sel[s].n_atoms)}
                   for s, vals in seg_rmsd.items()}
    for s, (obs, built) in split.items():
        n_built = len(set(built.resids))
        per_segment[s]["observed_only"] = {
            **_series(obs_rmsd[s], times, unit="angstrom"),
            "atoms": int(obs.n_atoms),
            "selection": f"{obs.n_atoms} backbone atoms of the {len(set(obs.resids))} OBSERVED "
                         f"residues of {s}, superposed on themselves; the {n_built} built "
                         "residues are left out",
        }
        per_segment[s]["modelled_only"] = {
            **_series(built_rmsd[s], times, unit="angstrom"),
            "atoms": int(built.n_atoms),
            "selection": f"{built.n_atoms} backbone atoms of the {n_built} BUILT residues of {s}, "
                         "carried by the superposition on the observed ones (mass-weighted)",
        }
    # The headline series is the segment that moved MOST, named. Taking the mildest, or an average
    # over chains, would hide exactly the chain an analyst needs to see.
    worst = max(per_segment, key=lambda s: abs(per_segment[s]["drift"])) if per_segment else None
    ref_note = (f"frame {ref_frame} of this trajectory, the first analysed frame"
                if ref_frame == first_frame else
                f"frame {ref_frame} of this trajectory, which you named with `reference_frame`; "
                f"the first ANALYSED frame is {first_frame}")

    return {
        "frames_analysed": len(times),
        "protein_rmsd": {
            **(per_segment[worst] if worst else _series([], [], unit="angstrom")),
            "segid": worst,
            "selection": (f"{seg_sel[worst].n_atoms} backbone atoms (N, CA, C, O) of segment "
                          f"{worst}, the one of {len(per_segment)} protein segment(s) with the "
                          f"largest DRIFT - which is not the same as the one that moved most: a "
                          f"chain oscillating widely about one position has a drift near zero"
                      if worst else "no protein segment found"),
        "headline_rule": "largest absolute drift (last quarter mean minus first quarter mean), "
                         "NOT largest range or maximum. Compare `per_segment` for those",
            "per_segment": per_segment,
            "reference": ref_note,
            "criterion": "mass-weighted superposition of each segment on ITSELF; this is the "
                         "equilibration control, not a result about the ligand",
            "per_segment_note": (
                f"one series per protein segment ({', '.join(per_segment)}). They are NOT "
                "superposed together: with `wrapAll on` two chains are two molecules, and one of "
                "them re-imaged while the other is not would read as tens of angstroms of RMSD "
                "with nothing having moved. Motion of one chain RELATIVE to another is therefore "
                "not measured here at all"),
        },
        "ligand_rmsd": {
            **_series(lig_rmsd, times, unit="angstrom"),
            "selection": ligand["selection"],
            "reference": ref_note,
            "criterion": CRITERIA["superposition"],
            "superposed_on": f"{fit.n_atoms} pocket backbone atoms",
        },
        "com_distance": {
            **_series(com_dist, times, unit="angstrom"),
            "selection": f"centre of mass of {lig.n_atoms} ligand heavy atoms against "
                         f"{pocket_all.n_atoms} pocket atoms",
            "criterion": CRITERIA["com_distance"],
            "pbc_corrected_frames": pbc_corrected,
            "pbc_rule": "each frame's ligand is placed in the periodic image nearest to its "
                        "position in the PREVIOUS analysed frame (continuity), not nearest to the "
                        "pocket. A pocket anchor cannot tell a real displacement of more than half "
                        "a box from a wrapping jump and folds the real one back",
            "pbc_note": (
                f"{pbc_corrected} of {len(times)} frames had the ligand re-imaged before "
                "measuring, because it had jumped more than half a box since the previous frame. "
                "CHARMM-GUI writes `wrapAll on`, so NAMD does that re-imaging and it is not "
                "motion. The SAME corrected coordinates feed the ligand RMSD"
                if pbc_corrected else
                "no frame needed re-imaging: the ligand never moved more than half a box between "
                "two consecutive analysed frames. That is what was measured, and it is weaker "
                "than saying the ligand never crossed a boundary - a molecule can be re-imaged "
                "and come back within one frame without this test noticing"),
            "non_orthorhombic_frames": non_orthorhombic,
            "non_orthorhombic_note": (
                f"{non_orthorhombic} of {len(times)} frames have a box that is not rectangular, "
                "and NO re-imaging was applied to them: the rule used here assumes 90 deg angles. "
                "Treat their values as unverified" if non_orthorhombic else
                "every frame with a box has a rectangular one, which is what the re-imaging rule "
                "assumes"),
            "frames_without_box": frames_without_box,
            "frames_without_box_note": (
                f"{frames_without_box} of {len(times)} frames declare NO box at all, so no "
                "re-imaging could be applied to them and nothing here says whether the ligand was "
                "wrapped in them. A DCD written without unit-cell information does this"
                if frames_without_box else
                "every frame declares a box"),
            "pocket_split_frames": pocket_split,
            "pocket_split_note": (
                f"in {pocket_split} of {len(times)} frames the pocket's own atoms span more than "
                "half a box in some direction, which is what a pocket torn across a periodic "
                "boundary looks like. In those frames the pocket centre of mass and the "
                "superposition are both unreliable, and NOTHING was corrected for it. This is the "
                "situation of a site at a chain interface under `wrapAll on`"
                if pocket_split else
                "in no frame do the pocket atoms span more than half a box, so the pocket was not "
                "torn across a periodic boundary"),
        },
        "rmsf": {
            "per_residue": [{"residue": name, "rmsf_angstrom": round(float(v), 4)}
                            for name, v in zip(ca_labels, rmsf, strict=True)],
            "unit": "angstrom",
            "selection": f"the CA of {pocket_ca.n_atoms} of "
                         f"{len(pocket['residues'])} pocket residues",
            "criterion": "deviation about each atom's own mean position after the same pocket "
                         "backbone superposition; each value is labelled from the CA atom it was "
                         "computed from, not from a parallel list",
            "residues_without_ca": missing_ca,
            "residues_without_ca_note": (
                f"{len(missing_ca)} pocket residue(s) have no atom named CA and carry NO rmsf "
                f"value: {', '.join(missing_ca)}" if missing_ca else
                "every pocket residue contributed exactly one CA"),
            "max": {"residue": ca_labels[int(rmsf.argmax())],
                    "rmsf_angstrom": round(float(rmsf.max()), 4)} if len(rmsf) else None,
        },
    }


# ---------------------------------------------------------------------------------------------
# 6. Hydrogen bonds
# ---------------------------------------------------------------------------------------------

_GEOMETRIC_SELECTION = re.compile(
    r"\b(around|sphzone|sphlayer|cyzone|cylayer|point|prop|isolayer|byres)\b")


def hydrogen_bonds(universe, resname: str, *, first_frame: int, frames_analysed: int,
                   segid: str = "") -> dict:
    """Ligand-protein hydrogen bonds, as a sampled occupancy per donor-H-acceptor triplet.

    `segid` restricts the ligand side to ONE copy, and it has to: without it a system with two
    copies of the residue had its pocket, RMSD and centre of mass measured on the copy the caller
    named while this list quietly reported bonds of the other copy as well, under a heading that
    said one copy of two was measured.

    Donors and acceptors are found through the bonds the PSF already carries, so connectivity is
    read rather than guessed from distances.

    The occupancy is deliberately NOT called a fraction of time. A hydrogen bond breaks and
    re-forms on a picosecond scale and a DCD stores a frame every tens of picoseconds, so what is
    counted is the frames in which the geometry held, out of the frames that exist.
    """
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA

    # MDAnalysis's own guess_acceptors()/guess_hydrogens() are NOT used, and this is the single
    # most consequential decision in this function. They filter by PARTIAL CHARGE
    # (acceptors need charge < -0.5), which is not a hydrogen-bond criterion and cannot be
    # declared honestly next to a geometric one. On the real ligand it excluded BOTH sulfonamide
    # oxygens (-0.370) and the carbonyl (-0.482), hiding a bond to the DFG glycine that holds in
    # 36% of frames, while admitting an aromatic CARBON of charge -0.555 as an acceptor.
    # Here: acceptors are N and O by mass, and a donor is whatever N or O the hydrogen is bonded
    # to in the topology. Both rules are in CRITERIA["hbond"], and neither reads a charge.
    lig_sel = f"resname {resname}" + (f" and segid {segid}" if segid else "")
    scope = f"(({lig_sel}) or protein)"
    acceptors_sel = f"{scope} and {HB_NO_BY_MASS}"
    hydrogens_sel = f"{scope} and {HB_H_BY_MASS} and bonded {HB_NO_BY_MASS}"
    if not universe.select_atoms(hydrogens_sel).n_atoms:
        raise MdAnalyzeError(
            "No hydrogen in the protein or the ligand is bonded to a nitrogen or an oxygen, so no "
            "donor can be identified. A PSF carries bonds; a topology without them (a bare PDB, "
            "for instance) cannot support hydrogen-bond analysis.")
    # `update_selections=False`: the selections are TOPOLOGICAL (resname, protein, mass, bonded),
    # so re-evaluating them on every frame - MDAnalysis's default - returns the same atoms each
    # time and cost 119.5 s of a 135 s analysis on the real 140k-atom, 100-frame production,
    # against 1.6 s without it, with the same 191 rows (compared element by element, 2026-11-12).
    # That pushed the tool past the 60 s Claude Desktop's remote bridge waits, and the call died.
    # A geometric keyword in either selection would make per-frame updates necessary again, so
    # one is refused here rather than silently frozen at frame 0.
    for sel in (lig_sel, acceptors_sel, hydrogens_sel):
        if _GEOMETRIC_SELECTION.search(sel):
            raise MdAnalyzeError(
                f"The hydrogen-bond selection {sel!r} has a geometric term, and these selections "
                "are evaluated once, not per frame.")
    hba = HBA(universe=universe, between=[lig_sel, "protein"],
              donors_sel=None, hydrogens_sel=hydrogens_sel, acceptors_sel=acceptors_sel,
              d_a_cutoff=HB_DA_CUTOFF, d_h_a_angle_cutoff=HB_ANGLE_CUTOFF,
              update_selections=False)
    hba.run(start=first_frame, verbose=False)

    # `count_by_ids()` is NOT used, and deliberately. It converts to atom IDS
    # (`np.array([d.ids, h.ids, a.ids]).T`), which are not positional indices: indexing
    # `universe.atoms` with one happens to agree for a PSF, whose parser stores
    # `id = line number - 1`, and silently names the WRONG residue for a topology numbered any
    # other way - a hallucinated residue name wearing a measurement's face. It also raises
    # NoDataError on a topology that carries no ids at all.
    # `results.hbonds` stores INDICES (hbond_analysis.py lines 827-829), so the frames are
    # counted here, from the indices, and the labels come from the atoms those indices point at.
    rows = hba.results.hbonds
    counts = {}
    for row in rows:
        key = (int(row[1]), int(row[2]), int(row[3]))
        counts[key] = counts.get(key, 0) + 1

    bonds = []
    for (donor_ix, hydrogen_ix, acceptor_ix), count in counts.items():
        donor = universe.atoms[donor_ix]
        acceptor = universe.atoms[acceptor_ix]
        bonds.append({
            "donor": f"{donor.segid}:{donor.resname}{donor.resid}:{donor.name}",
            "hydrogen": universe.atoms[hydrogen_ix].name,
            "acceptor": f"{acceptor.segid}:{acceptor.resname}{acceptor.resid}:{acceptor.name}",
            "frames_present": int(count),
            "frames_analysed": frames_analysed,
            "occupancy_sampled": (round(int(count) / frames_analysed, 4)
                                  if frames_analysed else None),
            "ligand_side": "donor" if donor.resname == resname else "acceptor",
        })
    # most present first, then by name, so the order never depends on dictionary insertion
    bonds.sort(key=lambda b: (-b["frames_present"], b["donor"], b["acceptor"]))
    per_frame = [0] * frames_analysed
    for row in rows:
        idx = int(row[0]) - first_frame
        if 0 <= idx < frames_analysed:
            per_frame[idx] += 1
    return {
        "pairs": bonds,
        "pairs_found": len(bonds),
        "criterion": CRITERIA["hbond"],
        "occupancy_note": "occupancy_sampled is the share of ANALYSED FRAMES in which the "
                          "geometry held, not the share of elapsed time: hydrogen bonds break and "
                          "re-form faster than a DCD stores frames",
        # `_series` stamps every series with a times_note that promises picoseconds. Here the axis
        # is frame INDICES, so that note would be a false guarantee on the one field whose axis is
        # not time: a reader quoting it turned frame 24 into 24 ps of a 240 ps window. It is
        # replaced, not merely contradicted one line below.
        "count_per_frame": {**{k: v for k, v in
                               _series([float(c) for c in per_frame],
                                       list(range(frames_analysed)), unit="bonds").items()
                               if k not in ("times_ps", "times_note")},
                            # The KEY is renamed, not just annotated. Leaving `times_ps` in place
                            # with a note beside it still hands a consumer a field whose name
                            # promises picoseconds, and frame 24 of a 240 ps window became 24 ps.
                            "frame_indices": list(range(frames_analysed)),
                            "note": "this series has NO time axis: `frame_indices` are positions "
                                    "in the analysed range, not picoseconds. There is deliberately "
                                    "no `times_ps` field here"},
        "connectivity": "donors are the N or O each hydrogen is BONDED to in the topology, which "
                        "a PSF carries; acceptors are N and O by mass. No partial charge is used, "
                        "and nothing is guessed from distances",
        "eligible": {
            "acceptors": int(universe.select_atoms(acceptors_sel).n_atoms),
            "hydrogens": int(universe.select_atoms(hydrogens_sel).n_atoms),
            "acceptors_on_the_ligand": int(
                universe.select_atoms(f"({lig_sel}) and {HB_NO_BY_MASS}").n_atoms),
            "hydrogens_on_the_ligand": int(universe.select_atoms(
                f"({lig_sel}) and {HB_H_BY_MASS} and bonded {HB_NO_BY_MASS}").n_atoms),
            "note": "the atoms that WERE eligible before any geometry was tested. The two totals "
                    "count the protein AND the ligand together; only LIGAND-PROTEIN pairs were "
                    "ever reported, so most of those atoms could never produce a listed bond. The "
                    "two `_on_the_ligand` counts are the side that bounds what this list can hold",
        },
        "ligand_selection": lig_sel,
    }


# ---------------------------------------------------------------------------------------------
# 7. The reading, kept apart from the measurements
# ---------------------------------------------------------------------------------------------

def reading(result: dict) -> dict:
    """An offered interpretation, built only from the numbers above and from stated thresholds.

    This is the one place in the module that says anything resembling a conclusion, and it says it
    as an offer. Every line names the number it came from and the threshold it was compared with,
    so the reader can disagree with the threshold without having to recompute anything.

    The four thresholds are COMMON CONVENTIONS and nothing more: not measurements of this system,
    and not derived from any source this project retrieved. No literature is cited because none
    was consulted; an earlier version of this docstring attributed them to "the literature", which
    was a claim without evidence. They were written down in advance, with the reason for each, in
    ESTADO_PROYECTO.md section 7.2 on 2026-11-11, while the first production chunk of the first
    simulation was still running and no production trajectory existed yet.

    That pre-registration is weaker than it sounds and is not claimed to be more. A RESTRAINED
    equilibration of the same system had already been measured when the thresholds were written,
    so they were chosen with those four numbers in view. What the date rules out is tuning against
    a production result; it does not rule out tuning against the equilibration.
    """
    notes, thresholds = [], {
        "protein_plateau_angstrom": 0.5,
        "ligand_stable_angstrom": 2.0,
        "com_stable_range_angstrom": 2.0,
        "hbond_persistent_occupancy": 0.5,
    }
    prot, lig = result["protein_rmsd"], result["ligand_rmsd"]
    com, hb = result["com_distance"], result["hbonds"]

    restraint = result.get("trajectory", {}).get("restraints", {})
    by_atom = restraint.get("held", {})
    # "held" per group, from the restraint file read atom by atom. When the file was not read, the
    # answer is the keyword alone, which is weaker: it says SOMETHING was held, not what.
    def _held(group):
        entry = by_atom.get(group, {})
        if entry.get("known"):
            return entry.get("held", 0) > 0
        return bool(restraint.get("atoms_held", False))

    held = _held("protein_backbone")
    lig_held = _held("ligand_heavy")
    held_warning = (
        " THIS LINE CARRIES NO EVIDENCE: the NAMD input of this stage holds these atoms in place "
        "({}), so the figure cannot move and the comparison passes by construction. Read it again "
        "on an unrestrained production trajectory.")
    drift = abs(prot["drift"])
    notes.append(
        f"Protein backbone RMSD of segment {prot.get('segid')} (the one of "
        f"{len(prot.get('per_segment', {})) or 1} with the largest drift) averages "
        f"{prot['last_quarter_mean']} A over the last quarter against "
        f"{prot['first_quarter_mean']} A over the first (a drift of "
        f"{round(drift, 4)} A, compared with a chosen {thresholds['protein_plateau_angstrom']} A). "
        + ("That is consistent with a plateau." if drift <= thresholds["protein_plateau_angstrom"]
           else "That is a drift rather than a plateau, so the protein may not be equilibrated and "
                "everything below is premature.")
        + (held_warning.format(by_atom.get("protein_backbone", {}).get("note", "restraints on"))
           if held else "")
        + (" " + prot["drift_warning"] + "." if prot.get("drift_warning") else ""))
    # The observed/built split. Offered as its own line, never folded into the one above: the
    # full-backbone comparison is the pre-registered control and keeps its own verdict.
    limit = thresholds["protein_plateau_angstrom"]
    per_seg = prot.get("per_segment", {})
    split = {s: v for s, v in per_seg.items() if "observed_only" in v}
    built_info = result.get("modelled_residues", {}).get("segments", {})
    head = prot.get("segid")
    if split:
        parts, warnings_ = [], []
        for s, v in split.items():
            o, m = v["observed_only"], v["modelled_only"]
            info = built_info.get(s, {})
            o_drift = round(abs(o["drift"]), 4)
            parts.append(
                f"{s}: the observed residues drift {o_drift} A over {o['atoms']} atoms "
                f"({'within' if o_drift <= limit else 'ABOVE'} the chosen {limit} A; full "
                f"segment {round(abs(v['drift']), 4)} A), with {info.get('count', '?')} of "
                f"{info.get('residues_total', '?')} residues built and left out "
                f"({info.get('ranges', '?')}); the built residues alone drift "
                f"{round(m['drift'], 4)} A and average {m['last_quarter_mean']} A over the last "
                "quarter")
            for series in (o, m):
                if series.get("drift_warning"):
                    warnings_.append(f"{s}: {series['drift_warning']}")
        if warnings_:
            verdict = ("NO STATEMENT ABOUT A PLATEAU CAN BE MADE FROM THESE FIGURES: "
                       + "; ".join(dict.fromkeys(warnings_)) + ".")
        elif head not in split:
            verdict = (f"The segment of the line above, {head}, is not among them, so this split "
                       "says nothing about its drift.")
        else:
            h_full = abs(per_seg[head]["drift"])
            h_obs = abs(per_seg[head]["observed_only"]["drift"])
            verdict = (
                f"For {head}, the segment of the line above, the full backbone is within the "
                "chosen figure but its observed residues alone are NOT, so the built residues "
                "offset part of the observed drift and the full-backbone figure hides it."
                if h_full <= limit < h_obs else
                f"For {head}, the segment of the line above, the full backbone is within the "
                "chosen figure already, so the split changes nothing there." if h_full <= limit else
                f"For {head}, the segment of the line above, the observed residues are within "
                "the chosen figure, so what exceeds it comes from the residues CHARMM-GUI built."
                if h_obs <= limit else
                f"For {head}, the segment of the line above, the observed residues still drift "
                "beyond the chosen figure, so the built residues do not account for the whole "
                "drift.")
        notes.append(
            "Split into observed and BUILT residues (built = absent from, or without a full "
            "backbone in, CHARMM-GUI's .orig.pdb; observed = present in the structure uploaded to "
            "CHARMM-GUI, which is experimental only if that was a deposited structure - not "
            "checked here): " + "; ".join(parts) + ". " + verdict
            + " THIS SPLIT WAS NOT PRE-REGISTERED: it was added on 2026-11-12, after the "
              "full-backbone control above failed on this project's first production run. Its "
              "reason - built residues are a guess, not a measurement - does not depend on the "
              "result, but the decision to look came after seeing it. The line above remains the "
              "pre-registered control."
            + (held_warning.format(by_atom.get("protein_backbone", {}).get("note", "restraints on"))
               if held else ""))
    # Independent of the split above: every segment that drifts beyond the figure and whose built
    # residues are UNKNOWN says so, whatever happened in the other segments.
    for s, v in per_seg.items():
        info = built_info.get(s, {})
        if abs(v["drift"]) > limit and info and not info.get("known"):
            notes.append(
                f"Segment {s} drifts {round(abs(v['drift']), 4)} A over its full backbone, and "
                f"whether any of its residues were built rather than observed is UNKNOWN "
                f"({info.get('reason', 'not looked up')}), so that drift cannot be split into "
                "observed and built residues.")
    in_pocket = result.get("pocket", {}).get("modelled_residues_in_pocket") or []
    if result.get("pocket", {}).get("modelled_residues_known") is False:
        notes.append(
            "Whether any pocket residue was built rather than observed is UNKNOWN ("
            + result["pocket"].get("modelled_residues_in_pocket_note", "") + "). If some were, "
            "the ligand figures are superposed on a guessed backbone.")
    if in_pocket:
        notes.append(
            f"{len(in_pocket)} pocket residue(s) were built, not observed ({', '.join(in_pocket)}). "
            "The ligand RMSD is superposed on the pocket backbone, theirs included, and the "
            "pocket itself was selected around a ligand placed next to a guessed geometry.")
    notes.append(
        f"Ligand RMSD reaches {lig['max']} A at most and averages {lig['last_quarter_mean']} A over "
        f"the last quarter (chosen comparison: {thresholds['ligand_stable_angstrom']} A). "
        + ("The pose held within that figure." if lig["last_quarter_mean"] <=
           thresholds["ligand_stable_angstrom"] else
           "The pose moved beyond that figure, which is what a pose that does not survive in "
           "solution looks like.")
        + (held_warning.format(by_atom.get("ligand_heavy", {}).get("note", "restraints on"))
           if lig_held else ""))
    notes.append(
        f"Centre-of-mass distance stays in a {com['range']} A range (chosen comparison: "
        f"{thresholds['com_stable_range_angstrom']} A), mean {com['mean']} A. "
        + ("No dissociation is visible in this window." if com["range"] <=
           thresholds["com_stable_range_angstrom"] else
           "That spread is wide enough that the ligand may be leaving the site; the per-frame "
           "values say when.")
        + (held_warning.format(by_atom.get("ligand_heavy", {}).get("note", "restraints on"))
           if lig_held else ""))
    persistent = [b for b in hb["pairs"]
                  if (b["occupancy_sampled"] or 0) >= thresholds["hbond_persistent_occupancy"]]
    shown = persistent[:5]
    notes.append(
        f"{len(persistent)} of {hb['pairs_found']} hydrogen-bond pairs held in at least "
        f"{thresholds['hbond_persistent_occupancy']:.0%} of the analysed frames"
        + (": " + ", ".join(f"{b['donor']} -> {b['acceptor']} ({b['occupancy_sampled']:.0%})"
                            for b in shown)
           + (f" ({len(shown)} of {len(persistent)} listed; the full list is in hbonds.pairs)"
              if len(persistent) > len(shown) else "") if persistent else "") + "."
        # The three lines above carry a warning when the atoms are held; this one did not, which
        # left it as the only figure an analyst could quote believing it clean. A hydrogen bond
        # between two atoms that are both on springs is held in place along with them.
        + (" This figure is not marked as carrying no evidence, because an occupancy is not a "
           "displacement - but the atoms on both sides of these bonds are held in place "
           f"({by_atom.get('ligand_heavy', {}).get('note', 'see trajectory.restraints')}), so the "
           "geometry was kept near its starting value and the occupancy is an upper bound rather "
           "than a result." if lig_held or held else ""))
    if com.get("pbc_corrected_frames"):
        notes.append(
            f"{com['pbc_corrected_frames']} of {com['n']} frames had the ligand re-imaged before "
            "measuring, because it jumped more than half a box since the previous frame. That is "
            "CHARMM-GUI's `wrapAll on` at work and not motion; the numbers above already use the "
            "corrected coordinates, but a run where this happens often is one to look at frame by "
            "frame rather than through a summary.")
    if com.get("non_orthorhombic_frames"):
        notes.append(
            f"{com['non_orthorhombic_frames']} of {com['n']} frames have a box that is not "
            "rectangular, and the re-imaging rule used here assumes 90 deg angles. Those frames "
            "were left uncorrected, so treat both distance and RMSD as unverified.")
    if com.get("frames_without_box"):
        notes.append(
            f"{com['frames_without_box']} of {com['n']} frames declare no box at all, so no "
            "re-imaging could be applied to them. If the run used periodic boundaries, every "
            "number above is unverified for those frames.")
    if com.get("pocket_split_frames"):
        notes.append(
            f"In {com['pocket_split_frames']} of {com['n']} frames the pocket's own atoms span "
            "more than half a box, which is what a site torn across a periodic boundary looks "
            "like. The superposition and the centre of mass are both unreliable there, and "
            "nothing was corrected for it: the lines above should not be read for this run until "
            "that is understood.")

    return {
        "offered": notes,
        "thresholds_used": thresholds,
        "status": "READING OFFERED, NOT MEASURED",
        "thresholds_source": "common conventions, pre-registered with a reason for each in "
                             "ESTADO_PROYECTO.md section 7.2 on 2026-11-11, before any production "
                             "trajectory existed. No literature is cited because none was "
                             "consulted. The 2.0 A for the ligand RMSD matches the docking "
                             "validation threshold of section 7.1 BY CONVENTION, not by "
                             "derivation: that one compares a calculated pose with a crystal "
                             "pose, this one asks whether a pose persists over time",
        "disclaimer": "These lines are an interpretation built from the numbers in this same "
                      "result and from the thresholds listed here, which are conventions and not "
                      "measurements of this system. The measurements above stand on their own; "
                      "this reading does not. Apply your own judgement.",
    }


# ---------------------------------------------------------------------------------------------
# 8. The entry point
# ---------------------------------------------------------------------------------------------

def analyze_md(folder, *, topology=None, ligand_resname: str = "LIG", ligand_segid: str = "",
               pocket_cutoff: float = POCKET_CUTOFF, skip_frames: int = -1,
               reference_frame: int = -1, include_series: bool = True, stage=None) -> dict:
    """Measure the trajectory in `folder` and return every metric with its selection and criterion.

    `skip_frames` defaults to the minimization frames declared by the NAMD inputs, which is the
    honest default; pass 0 to keep them, or a number to drop more (to skip a settling period, for
    instance).

    `reference_frame` is the frame everything is measured AGAINST, and it moves three things at
    once: the pocket is selected around the ligand in that frame, both RMSDs use that frame's
    coordinates as their zero, and the ligand's periodic image in that frame is the origin of the
    continuity chain. It is an index into the WHOLE chained trajectory, counted from 0, not into
    the analysed part, and it defaults to the first analysed frame. Naming a frame outside the
    analysed range is allowed and sometimes right (a crystal pose kept as frame 0), but then the
    RMSD series does not start at 0 and the reported `reference` says so.

    Every refusal comes back as `{"ok": False, "error": ...}` with whatever the refusal knows
    attached, the same contract as the rest of the project's tools: a trajectory that cannot be
    measured is a result the analyst needs to read, not an exception for Claude Desktop to
    swallow.
    """
    try:
        return _analyze(folder, topology=topology, ligand_resname=ligand_resname,
                        ligand_segid=ligand_segid,
                        pocket_cutoff=pocket_cutoff, skip_frames=skip_frames,
                        reference_frame=reference_frame, include_series=include_series,
                        stage=stage)
    except MdAnalyzeError as exc:
        return {"ok": False, "folder": str(folder), "error": str(exc),
                **getattr(exc, "extra", {})}
    except (OSError, ValueError, KeyError, IndexError) as exc:
        logger.exception("analyze_md failed")
        return {"ok": False, "folder": str(folder),
                "error": f"{type(exc).__name__}: {exc}",
                "note": "this is a failure of the tool or of the files, not a measurement"}


def _span(measured) -> float:
    """Elapsed picoseconds between the first and the last analysed frame, from the frames."""
    times = measured["ligand_rmsd"]["times_ps"]
    return round(float(times[-1] - times[0]), 3) if times else 0.0


def _analyze(folder, *, topology, ligand_resname, ligand_segid, pocket_cutoff, skip_frames,
             reference_frame, include_series, stage=None) -> dict:
    """The body of analyze_md. Raises; analyze_md turns a refusal into `ok: false`."""
    import MDAnalysis as mda

    found = find_inputs(folder, topology, stage)
    plan = trajectory_plan(found["psf"], found["dcds"])
    mini = minimization_frames(found["inputs"], found["dcds"])
    held = restraints(found["inputs"], found["dcds"])

    first = mini["frames"] if skip_frames < 0 else int(skip_frames)
    if first >= plan["frames_total"]:
        raise MdAnalyzeError(
            f"Dropping {first} frame(s) leaves nothing of the {plan['frames_total']} in the "
            "trajectory.", frames_total=plan["frames_total"])

    universe = mda.Universe(str(found["psf"]), [str(d) for d in found["dcds"]])
    try:
        if len(universe.trajectory) != plan["frames_total"]:
            logger.warning("chained frame count %d differs from the per-file sum %d",
                           len(universe.trajectory), plan["frames_total"])

        ligand = select_ligand(universe, ligand_resname, ligand_segid)
        # A reference frame the trajectory does not have is the USER's mistake, not the tool's, so
        # it is refused with the frame count instead of escaping as "IndexError: Index 100 exceeds
        # length of trajectory" under a note calling it a failure of the tool or of the files.
        # A negative value other than the -1 default is refused too: Python would read it as
        # "counting from the end", and this code silently meant "the first analysed frame".
        ref = first if reference_frame == -1 else int(reference_frame)
        if ref < 0 or ref >= len(universe.trajectory):
            raise MdAnalyzeError(
                f"reference_frame {reference_frame} is not a frame of this trajectory, which has "
                f"{len(universe.trajectory)} (0 to {len(universe.trajectory) - 1}). Leave it out "
                "to use the first analysed frame; negative values are not counted from the end.",
                frames_total=int(len(universe.trajectory)))
        pocket = select_pocket(universe, ligand["atoms"], pocket_cutoff, ref)
        # CHARMM-GUI keeps the `.orig.pdb` files one level above `namd/` (in `system/`); the
        # folder the user named and the topology's own folder are searched too.
        # Only next to the topology: a `system/` elsewhere under `folder` can belong to another
        # build with the same numbering, and would be accepted as "known" for this one. "The
        # topology" includes every other .psf found that is BYTE-IDENTICAL to the one measured:
        # a build_md_job folder holds `step3_input.psf` at its root AND under run/system/namd/,
        # the search picks the root copy, and the `.orig.pdb` files sit in run/system/ - next to
        # the other copy of the same file, not of another build.
        import filecmp
        root = Path(folder)
        psf = Path(found["psf"])
        copies = [psf] + [p for p in found.get("psfs_all", [])
                          if _file_identity(p) != _file_identity(psf)
                          and filecmp.cmp(p, psf, shallow=False)]
        search = [d for p in copies for d in (p.parent, p.parent.parent)]
        built = modelled_residues(universe, search, root=root)
        built["searched"] = list(dict.fromkeys(_relative(d, root) for d in search))
        built["searched_note"] = (
            "the folder of the measured topology and the one above it"
            + (f", and the same for {len(copies) - 1} other .psf byte-identical to it"
               if len(copies) > 1 else ""))
        known_built = {s: set(v["residues"]) for s, v in built["segments"].items()
                       if v.get("known") and v.get("residues")}
        in_pocket = [label for label, r in zip(pocket["residues"], pocket["atoms"].residues)
                     if int(r.resid) in known_built.get(r.segid, ())]
        unknown_pocket_segs = sorted(s for s in pocket["segments"]
                                     if not built["segments"].get(s, {}).get("known"))
        # An empty list beside an UNKNOWN reads as "none". When any pocket segment is unknown, the
        # list is None unless something was found, and the note leads with the unknown.
        pocket["modelled_residues_known"] = not unknown_pocket_segs
        pocket["modelled_residues_in_pocket"] = (in_pocket or None) if unknown_pocket_segs \
            else in_pocket
        n_known = sum(1 for r in pocket["atoms"].residues if r.segid not in unknown_pocket_segs)
        found_part = (f"{len(in_pocket)} of {n_known} pocket residues were built, "
                      f"not observed: {', '.join(in_pocket)}" if in_pocket else
                      f"0 of {n_known} pocket residues were built")
        pocket["modelled_residues_in_pocket_note"] = (
            f"UNKNOWN whether any pocket residue of segment(s) {', '.join(unknown_pocket_segs)} "
            "was built (see modelled_residues.segments); "
            + ("in the other segment(s), " + found_part
               if len(unknown_pocket_segs) < len(pocket["segments"]) else
               "no other segment holds pocket residues")
            if unknown_pocket_segs else found_part)
        measured = measure(universe, ligand, pocket, first_frame=first, ref_frame=ref,
                           modelled=known_built)
        measured["hbonds"] = hydrogen_bonds(universe, ligand_resname, first_frame=first,
                                            frames_analysed=measured["frames_analysed"],
                                            segid=ligand_segid)
        # WHICH atoms the restraint file holds, matched by atom index against the atoms actually
        # measured. This is what turns "the input mentions a restraint file" into a statement
        # about these numbers.
        flags, labels = held.pop("_flags", (None, None))
        total = int(len(universe.atoms))
        held["held"] = {
            "ligand_heavy": held_selection(flags, labels, ligand["atoms"], "ligand heavy", total),
            "protein_backbone": held_selection(
                flags, labels, universe.select_atoms("protein and name N CA C O"),
                "protein backbone", total),
        }
    finally:
        universe.trajectory.close()
    held.pop("_flags", None)       # the per-atom arrays never leave this function

    out = {
        "ok": True,
        "folder": str(Path(folder)),
        "topology": {"file": found["psf"].name, "selection": found["psf_selection"],
                     "atoms": int(len(universe.atoms)),
                     "segments": sorted(set(universe.atoms.segids))},
        "trajectory": {
            **plan,
            "stage": found["stage"],
            "stage_inputs": found["stage_inputs"],
            "stage_selection": found["stage_selection"],
            "stages_present": found["stages_present"],
            "stage_note": found["stage_note"],
            "left_out": found["left_out"],
            "files": [d.name for d in found["dcds"]],
            "frames_dropped": first,
            "frames_analysed": measured["frames_analysed"],
            "frames_note": f"{measured['frames_analysed']} of {plan['frames_total']} frames "
                           f"analysed; {first} dropped from the start",
            "minimization": mini,
            "minimization_note": (
                f"{mini['frames']} leading frame(s) come from `minimize` in {mini['read_from']} "
                "and are NOT dynamics: minimization has no time axis. "
                + ("They were dropped." if first >= mini["frames"] else
                   "They were KEPT because skip_frames said so; the time axis below is fiction "
                   "for those frames.")
                if mini["frames"] else
                f"{len(mini['inputs_read'])} NAMD input(s) of this stage were read "
                f"({', '.join(mini['inputs_read'])}) and none declares a minimization, so no "
                "leading frame was treated as one"
                if mini["inputs_read"] else
                "NO NAMD input of this stage was found, so whether the trajectory opens with "
                "minimization frames is UNKNOWN, not absent. Nothing was dropped. Pass "
                "`skip_frames` yourself if you know the answer"),
            "restraints": held,
            "span_ps": _span(measured),
            "span_ns": round(_span(measured) / 1000, 5),
            "span_note": "the elapsed time between the FIRST and the LAST analysed frame, read "
                         "from the frames themselves. It is one frame spacing shorter than "
                         "frames x dt, which counts an interval that is not in the file",
        },
        "ligand": {k: v for k, v in ligand.items() if k != "atoms"},
        "pocket": {k: v for k, v in pocket.items() if k not in ("atoms", "backbone")},
        "modelled_residues": built,
        "protein_rmsd": measured["protein_rmsd"],
        "ligand_rmsd": measured["ligand_rmsd"],
        "com_distance": measured["com_distance"],
        "rmsf": measured["rmsf"],
        "hbonds": measured["hbonds"],
        "criteria": CRITERIA,
        "not_measured": [
            "Binding affinity. Nothing here computes a free energy, and none of these metrics is "
            "one: a stable pose is not a potent ligand.",
            "Bridging waters. Left out of this version on purpose.",
            "Whether the pose is correct. With no apo control, a change in the protein cannot be "
            "attributed to the ligand.",
        ],
    }
    out["reading"] = reading(out)
    if not include_series:
        for key in ("protein_rmsd", "ligand_rmsd", "com_distance"):
            out[key].pop("values", None)
            out[key].pop("times_ps", None)
        # per_segment carries a FULL series per protein chain, so dropping only the three series
        # above left the "summary only" reply with one series per chain in it: ten chains over a
        # thousand frames meant ten thousand numbers in a reply that said it had none.
        for series in out["protein_rmsd"].get("per_segment", {}).values():
            series.pop("values", None)
            series.pop("times_ps", None)
            for part in ("observed_only", "modelled_only"):
                series.get(part, {}).pop("values", None)
                series.get(part, {}).pop("times_ps", None)
        out["hbonds"]["count_per_frame"].pop("values", None)
        out["hbonds"]["count_per_frame"].pop("times_ps", None)
    return out
