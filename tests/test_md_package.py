"""Offline test for build_md_package and check_charmm_gui_output (protein_mcp/files/md_package.py,
Phase 3 step 4).

The manifest, the real 3OG7 mmCIF and the real cleaned receptor come from test_choice/test_dock;
the poses are the real Vina 1.2.7 output of 2026-11-05 (tests/fixtures/dock). The REMARK 465 block
built from the mmCIF is compared against the block RCSB writes in its own PDB file.

CHARMM-GUI's answer is a 42 MB .tgz, too big to keep here, so the archives are built by this test
from the delivered files, copying the shapes and the defects of the four REAL jobs of ESTADO 8.4:
segments PROA/PROB/HETA, S1->S and CL1->CL, the lone pair LP1 on the chlorine, a reoriented frame,
and, for the controls that MUST fail, the real defects measured on jobs 9401967361 (no ligand),
9417550068 (gaps bridged by bonds of 8.7-24.0 A) and 9418371365 (ANISOU read as atoms, protein
expanded from Rg 18.6 to 22.8 A).

The two real archives were also checked by hand on 2026-11-08: the good one (9418739369) passes
with mean CA-CA change 0.1285/0.1186 A and the closest modelled loop at 4.76 A from the ligand;
the broken one (9418371365) is rejected by six independent checks. Those numbers are asserted in
test_thresholds_separate_the_real_jobs so a change of threshold cannot silently accept the broken
job; the archives themselves are not in the repository.

Run: uv run python tests/test_md_package.py
"""
import gzip
import hashlib
import json
import math
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

from protein_mcp import dossier
from protein_mcp.files import choice, dock, md_package

sys.path.insert(0, str(Path(__file__).parent))
from test_choice import ACC, FIX  # noqa: E402
from test_dock import DFIX, make_prep  # noqa: E402

POCKET = "032 C"
# measured on the real jobs 9418739369 (good) and 9418371365 (broken), 2026-11-08
REAL_GOOD = {"mean_dd": (0.1285, 0.1186), "rg_change": (0.02, 0.06), "loop_to_ligand": 4.76,
             "ligand_dd_max": 0.0011, "min_contact": 2.79}
REAL_BROKEN = {"mean_dd": (6.0251, 6.7511), "rg_change": (4.15, 4.93), "garbage_atoms": (1974, 2028)}
REAL_BRIDGES = [8.72, 16.15, 10.21, 24.05, 9.54]     # job 9401967361, gaps not filled


# ----------------------------------------------------------------- set-up

def prep_with_a_docking_run():
    """The temp dossier of test_dock plus a docking block built from the real Vina output."""
    prep, manifest = make_prep("md_package")
    poses, scores = dock.read_poses((DFIX / "vina_3OG7_032_C_regenerated.poses.pdbqt").read_text(encoding="utf-8"))
    out_dir = prep / "dock" / "3OG7_032_C"
    out_dir.mkdir(parents=True)
    sdf = out_dir / "poses_032_regenerated.sdf"
    sdf.write_text(dock.poses_sdf(poses, scores), encoding="utf-8")
    crystal = dock.read_sdf(prep / "ligands" / "3OG7_032_C.sdf")
    rmsds = dock.rmsd_to_crystal(poses, crystal)
    s = next(x for x in manifest["structures"] if x["pdb_id"] == "3OG7")
    s["docking"] = {POCKET: {"raw_sha256": s["file"]["sha256"], "runs": {"032_regenerated": {
        "computed_at": "2026-11-05T00:00:00Z", "selected_by": "arguments",
        "poses": [{"mode": i, "vina_score_kcal_mol": sc, "rmsd_to_crystal_pose_angstrom": r}
                  for i, (sc, r) in enumerate(zip(scores, rmsds), 1)],
        "files": {"poses_032_regenerated.sdf": {"path": sdf.relative_to(prep).as_posix(),
                                                "sha256": dock.sha256(sdf)}},
    }}}}
    choice.write_manifest(prep / "manifest.json", manifest)
    return prep


def rotate(xyz, shift=(37.0, -12.5, 8.25)):
    """A rigid rotation plus a translation: CHARMM-GUI re-orients and re-centres the system, and
    every geometry check here must be blind to that."""
    x, y, z = xyz
    c, s = math.cos(0.7), math.sin(0.7)
    return (c * x - s * y + shift[0], s * x + c * y + shift[1], z + shift[2])


def pdb_line(serial, name, resname, segid, resid, xyz, element):
    field = f" {name:<3}" if len(name) < 4 else name
    # columns as CHARMM-GUI writes them: name 13-16, residue 18-20, number 23-26, segid 73-76
    return (f"ATOM  {serial:5d} {field:<4} {resname:<3}  {resid:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00      {segid:<4}{element:>2}")


def charmm_gui_system(delivered: Path, *, with_ligand=True, fill_gaps=True, expand=1.0,
                      clash_run=None, to_model=(), ligand_shift=(0.0, 0.0, 0.0), mirror=False,
                      ligand_scale=1.0):
    """The lines of a CHARMM-GUI output system built from the delivered complex.

    expand > 1 scales the protein about its centroid, which is how the broken job shows up
    (Rg 18.6 -> 22.8 A). fill_gaps=False leaves the gaps unmodelled, so the chain keeps the long
    C-N bridges the first two real jobs had.
    """
    atoms = md_package.read_atoms(delivered)
    protein = [a for a in atoms if a["resname"] != md_package.LIGAND_RESNAME]
    ligand = [a for a in atoms if a["resname"] == md_package.LIGAND_RESNAME]
    centre = [sum(a["xyz"][i] for a in protein) / len(protein) for i in range(3)]
    lines, serial = [], 1
    by_chain = {}
    for a in protein:
        by_chain.setdefault(a["chain"], []).append(a)
    for chain, chain_atoms in sorted(by_chain.items()):
        segid = f"PRO{chain}"
        built = list(chain_atoms)
        if fill_gaps:
            backbone = {}
            for a in chain_atoms:
                backbone.setdefault(a["resid"], {})[a["name"]] = a["xyz"]
            for run in to_model:
                if run["chain"] != chain:
                    continue
                before, after = backbone.get(run["first"] - 1), backbone.get(run["last"] + 1)
                if not before or not after:
                    continue
                start, end = before.get("C"), after.get("N")
                count = run["last"] - run["first"] + 1
                # every backbone atom of the run evenly spaced between the flanking C and N, so
                # each C(i)-N(i+1) is gap/(3*count+1): short, as a filled gap must be
                steps = 3 * count + 1
                for k, resid in enumerate(range(run["first"], run["last"] + 1)):
                    for j, (name, element) in enumerate((("N", "N"), ("CA", "C"), ("C", "C"))):
                        frac = (3 * k + j + 1) / steps
                        pos = tuple(start[i] + frac * (end[i] - start[i]) for i in range(3))
                        if clash_run == (chain, run["first"]) and ligand:
                            pos = ligand[0]["xyz"]      # the modelled loop sitting on the ligand
                        built.append({"name": name, "resname": "GLY", "chain": chain, "resid": resid,
                                      "xyz": pos, "element": element})
        built.sort(key=lambda a: (a["resid"], {"N": 0, "CA": 1, "C": 2}.get(a["name"], 3)))
        for a in built:
            xyz = tuple(centre[i] + expand * (a["xyz"][i] - centre[i]) for i in range(3))
            resname = "HSD" if a["resname"] == "HIS" else a["resname"]
            lines.append(pdb_line(serial, a["name"], resname, segid, a["resid"], rotate(xyz), a["element"]))
            serial += 1
    if with_ligand:
        renamed = {"S1": "S", "CL1": "CL"}
        lig_centre = [sum(a["xyz"][i] for a in ligand) / len(ligand) for i in range(3)]

        def place(xyz):
            out = tuple(lig_centre[i] + ligand_scale * (xyz[i] - lig_centre[i]) for i in range(3))
            if mirror:                                 # a reflection keeps every internal distance
                out = (2 * lig_centre[0] - out[0], out[1], out[2])
            return rotate(tuple(out[i] + ligand_shift[i] for i in range(3)))

        for a in ligand:
            lines.append(pdb_line(serial, renamed.get(a["name"], a["name"]), md_package.LIGAND_RESNAME,
                                  "HETA", 1, place(a["xyz"]), a["element"]))
            serial += 1
        chlorine = next((a for a in ligand if a["name"] == "CL1"), None)
        if chlorine:                                  # the lone pair CGenFF adds on the halogen
            lines.append(pdb_line(serial, "LP1", md_package.LIGAND_RESNAME, "HETA", 1,
                                  place(chlorine["xyz"]), ""))
            serial += 1
    return "\n".join(lines) + "\nEND\n"


def make_archive(tmp: Path, delivered: Path, name="job.tgz", *, job="charmm-gui-9418739369",
                 uploaded=None, anisou_garbage=False, ligand_sdf=None, **kw):
    """A .tgz with the files check_charmm_gui_output reads, in the real layout."""
    root = tmp / "build" / job
    (root / "namd").mkdir(parents=True, exist_ok=True)
    (root / "lig").mkdir(exist_ok=True)
    system = charmm_gui_system(delivered, **kw)
    (root / "step1_pdbreader.pdb").write_text(system, encoding="ascii")
    (root / "namd" / "step3_input.pdb").write_text(system, encoding="ascii")
    shutil.copy(uploaded or delivered, root / delivered.name)
    (root / "step1_pdbreader.inp").write_text(
        "generate PROA setup warn first NTER last CTER\ngenerate HETA setup warn first none last none\n",
        encoding="ascii")
    (root / "step1_pdbreader.out").write_text(" Total charge =   10.00000\n", encoding="ascii")
    (root / "step3_pbcsetup.str").write_text(" SET A = 114\n SET B = 114\n SET C = 114\n", encoding="ascii")
    (root / "step2.2_ions.str").write_text(" SET NIONTOT = 256\n", encoding="ascii")
    # lig/LIG.sdf is the ligand file CHARMM-GUI really used; in the real job it is byte-identical
    # (ignoring line endings) to the one build_md_package delivered
    shutil.copy(ligand_sdf or delivered.with_name(delivered.stem + "_ligand.sdf"),
                root / "lig" / "LIG.sdf")
    (root / "lig" / "lig.rtf").write_text(
        "RESI LIG            0.000 ! param penalty=   7.000 ; charge penalty=   6.086\n", encoding="ascii")
    (root / "namd" / "step5_production.inp").write_text(
        "set temp                303.15;\ncutoff                  12.0;\ntimestep                2.0;\n"
        "numsteps                500000;\n", encoding="ascii")
    (root / "namd" / "sysinfo.dat").write_text('{"lonepair": true}', encoding="ascii")
    # what CHARMM-GUI hands to Modeller. The good job has one clean line per atom; the broken one
    # (9418371365) has the ANISOU U values in the coordinate columns, the third field shifted
    clean = ("ATOM      1  N   ASP   449     -15.047   9.355 -11.670  0.00  0.00           N  \n"
             "ATOM      3  CA  ASP   449     -14.538   9.164 -13.021  0.00  0.00           C  \n")
    garbage = ("ATOM      2  N   ASP   449    9855.0009321.000 236.000  0.00  0.00           N  \n"
               "ATOM      4  CA  ASP   449    10196.0009808.000 648.000  0.00  0.00           C  \n")
    (root / "3og7_proa.orig.pdb").write_text(clean + (garbage if anisou_garbage else ""),
                                             encoding="ascii")
    archive = tmp / name
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(root, arcname=job)
    return archive


# ----------------------------------------------------------------- pure pieces

MINIMAL_CIF = """data_TEST
loop_
_pdbx_unobs_or_zero_occ_residues.id
_pdbx_unobs_or_zero_occ_residues.PDB_model_num
_pdbx_unobs_or_zero_occ_residues.polymer_flag
_pdbx_unobs_or_zero_occ_residues.occupancy_flag
_pdbx_unobs_or_zero_occ_residues.auth_asym_id
_pdbx_unobs_or_zero_occ_residues.auth_comp_id
_pdbx_unobs_or_zero_occ_residues.auth_seq_id
_pdbx_unobs_or_zero_occ_residues.PDB_ins_code
_pdbx_unobs_or_zero_occ_residues.label_asym_id
_pdbx_unobs_or_zero_occ_residues.label_comp_id
_pdbx_unobs_or_zero_occ_residues.label_seq_id
1 1 Y 1 A MET 100 ? A MET 1
2 1 Y 0 A GLU 725 ? A GLU 2
3 1 N 1 A HOH 900 ? B HOH 3
#
"""


def write_cif(text: str) -> Path:
    path = Path(tempfile.mkdtemp()) / "x.cif"
    path.write_text(text, encoding="ascii")
    return path


def test_zero_occupancy_residues_are_not_called_missing():
    """occupancy_flag 0 is REMARK 475, not 465: those residues HAVE coordinates.

    Calling them missing told the analyst to model a residue that is already in the file
    (1M17 GLU A 725; 1IGY and 12E8 entirely). Found by the reviewer, 2026-11-08.
    """
    lines, source = md_package.remark465_lines(write_cif(MINIMAL_CIF))
    # the filtered row is reported, not dropped in silence: it has coordinates but is not experimental
    assert source.startswith("mmCIF; 1 residue(s) are in the file with zero occupancy"), source
    residues = [line[15:18].strip() + " " + line[21:26].strip()
                for line in lines[len(md_package._REMARK465_HEADER):]]
    assert residues == ["MET 100"], residues      # GLU 725 (flag 0) and HOH (not polymer) are out


def test_a_missing_category_is_not_the_same_as_no_gaps():
    """Absence of data must never be reported as "every residue is experimental"."""
    empty, source = md_package.remark465_lines(write_cif("data_TEST\n_cell.length_a 1.0\n"))
    assert empty == [] and "no _pdbx_unobs_or_zero_occ_residues category" in source, source
    without_tag = MINIMAL_CIF.replace("_pdbx_unobs_or_zero_occ_residues.occupancy_flag\n", "") \
                             .replace(" 1 A MET", " A MET").replace(" 0 A GLU", " A GLU") \
                             .replace(" 1 A HOH", " A HOH")
    lines, source = md_package.remark465_lines(write_cif(without_tag))
    assert lines == [] and "without occupancy_flag" in source, source


def test_remark465_matches_the_rcsb_block():
    cif = Path(tempfile.mkdtemp()) / "3OG7.cif"
    cif.write_bytes(gzip.decompress((FIX / "3OG7.cif.gz").read_bytes()))
    lines, source = md_package.remark465_lines(cif)
    assert source == "mmCIF"
    assert len(lines) == 84, len(lines)
    # columns of the PDB format, as RCSB writes them in files.rcsb.org/download/3OG7.pdb
    assert lines[7] == "REMARK 465     MET A   432".ljust(80), repr(lines[7])
    assert lines[-1] == "REMARK 465     SER B   630".ljust(80), repr(lines[-1])
    assert all(len(line) == 80 for line in lines)
    runs = md_package.missing_ranges(lines)
    assert [(r["chain"], r["first"], r["last"]) for r in runs] == [
        ("A", 432, 448), ("A", 545, 547), ("A", 597, 614), ("A", 627, 630),
        ("B", 432, 448), ("B", 601, 614), ("B", 627, 630)], runs
    assert runs[0]["length"] == 17
    kinds = md_package.classify_ranges(runs, {"A": (449, 720), "B": (449, 720)})
    assert [k["kind"] for k in kinds] == ["N-terminal", "internal", "internal", "internal",
                                          "N-terminal", "internal", "internal"], kinds
    assert all(k["absent_from_the_receptor"] for k in kinds)
    # a run whose residues ARE in the uploaded file is not a gap, whatever the mmCIF says
    with_one = md_package.classify_ranges(runs, {"A": (449, 720), "B": (449, 720)}, {("A", 546)})
    gap = next(k for k in with_one if (k["chain"], k["first"]) == ("A", 545))
    assert not gap["absent_from_the_receptor"] and gap["in_the_receptor"] == [546], gap


def test_an_entry_without_unobserved_residues_gets_no_block():
    cif = Path(tempfile.mkdtemp()) / "8TVK.cif"
    cif.write_bytes(gzip.decompress((FIX / "8TVK.cif.gz").read_bytes()))
    lines, source = md_package.remark465_lines(cif)
    # whatever this entry has, an empty block must never be written as a header with no residue
    assert lines == [] or len(lines) > len(md_package._REMARK465_HEADER), len(lines)
    assert source, "the source of the answer must always be stated"


def test_long_peptide_bonds_find_the_real_bridges():
    """The five bridges of job 9401967361, as C-N distances."""
    atoms = []
    for i, bridge in enumerate(REAL_BRIDGES):
        base = i * 1000.0
        atoms.append({"name": "C", "resname": "ALA", "chain": "A", "segid": "PROA",
                      "resid": 1 + i * 10, "xyz": (base, 0.0, 0.0), "element": "C"})
        atoms.append({"name": "N", "resname": "ALA", "chain": "A", "segid": "PROA",
                      "resid": 2 + i * 10, "xyz": (base + bridge, 0.0, 0.0), "element": "N"})
    found = md_package.long_peptide_bonds(atoms)
    assert [f["distance_angstrom"] for f in found] == REAL_BRIDGES, found


def test_internal_distance_change_ignores_rotation_but_not_distortion():
    points = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (0.0, 4.0, 0.0), (1.0, 2.0, 3.0)]
    moved = [rotate(p) for p in points]
    assert md_package.internal_distance_change(points, moved)["max_angstrom"] < 1e-6
    stretched = [(p[0] * 1.2, p[1] * 1.2, p[2] * 1.2) for p in points]
    assert md_package.internal_distance_change(points, stretched)["max_angstrom"] > 0.5


# ----------------------------------------------------------------- build_md_package

def test_build_writes_what_charmm_gui_needs():
    prep = prep_with_a_docking_run()
    result = md_package.build_md_package(ACC, "3OG7", POCKET)
    assert result["ok"], result
    out = Path(result["package_dir"])
    assert result["upload"] == ["3OG7.pdb", "3OG7_ligand.sdf"]
    pdb = (out / "3OG7.pdb").read_text(encoding="ascii")
    # the three rules the four real jobs established
    assert "ANISOU" not in pdb
    assert sum(1 for line in pdb.splitlines() if line.startswith("REMARK 465")) == 84
    assert (out / "3OG7.pdb").name == "3OG7.pdb"        # the entry id is read from the file name
    checks = result["checks"]
    assert checks["ligand_atoms"] == 51 and checks["ligand_heavy_atoms"] == 33
    assert checks["formal_charge"] == 0 and checks["conect_records"] == 51
    assert checks["pdb_and_sdf_atoms_in_agreement"] == "51 of 51"
    assert result["missing_residues_source"] == "mmCIF"
    assert result["recorded_in_manifest"], result["notes"]
    assert checks["min_receptor_ligand_heavy_distance_angstrom"] == 2.98, checks
    assert result["tick_model_missing_residues"] == ["A 545-547", "A 597-614", "A 627-630",
                                                     "B 601-614", "B 627-630"], result
    assert result["leave_unticked"] == ["A 432-448", "B 432-448"], result
    assert "modelled, not experimental" in result["modelled_not_experimental"]
    package = json.loads((out / "package.json").read_text(encoding="utf-8"))
    assert package["ligand"]["atom_names"][:4] == ["C1", "C2", "C3", "S1"], package["ligand"]["atom_names"][:4]
    assert len(package["ligand"]["coordinates"]) == 51
    assert len(package["pocket_residues_delivered"]) == 19, package["pocket_residues_delivered"]
    # the decision of 2026-11-06: the package is recorded in the manifest, not only on disk
    manifest = json.loads((dossier.DOSSIER_DIR / ACC / "prep" / "manifest.json").read_text(encoding="utf-8"))
    block = next(x for x in manifest["structures"] if x["pdb_id"] == "3OG7")["md_package"][POCKET]
    assert block["run_id"] == "032_regenerated" and block["dir"] == "md/3OG7_032_C"
    readme = (out / "README.md").read_text(encoding="utf-8")
    assert "A 597-614" in readme and "modelled, not experimental" in readme
    return prep, out / "3OG7.pdb"


def test_build_never_offers_a_run_that_is_in_the_file():
    """Control for the WIRING, not only for classify_ranges.

    3OG7 has no zero-occupancy residue, so every run is genuinely absent and a build that forgot
    to pass the residues it already has would look identical. This forces the question by making
    one residue of A 545-547 present in the receptor: the run must leave tick_model_missing_
    residues, must NOT be called a disordered tail, and must be reported on its own.
    """
    prep_with_a_docking_run()
    receptor = dossier.DOSSIER_DIR / ACC / "prep" / "clean" / "3OG7_receptor.pdb"
    lines = receptor.read_text(encoding="utf-8").splitlines()
    donor = next(i for i, line in enumerate(lines)
                 if line.startswith("ATOM") and line[21] == "A" and line[22:26].strip() == "549")
    lines.insert(donor, lines[donor][:22] + " 546" + lines[donor][26:])
    receptor.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = md_package.build_md_package(ACC, "3OG7", POCKET)
    assert result["ok"], result
    assert "A 545-547" not in result["tick_model_missing_residues"], result["tick_model_missing_residues"]
    assert "A 545-547" not in result["leave_unticked"], result["leave_unticked"]
    assert result["internal_runs_partly_in_the_file"] == ["A 545-547"], result
    readme = (Path(result["package_dir"]) / "README.md").read_text(encoding="utf-8")
    assert "NOT offered, and not tails: A 545-547" in readme, readme


def test_build_says_when_the_mmcif_has_no_record():
    """Absence of data must never reach the chat or the README as "there are no gaps"."""
    prep_with_a_docking_run()
    real = md_package.remark465_lines
    try:
        md_package.remark465_lines = lambda path: ([], "no _pdbx_unobs_or_zero_occ_residues category in the mmCIF")
        result = md_package.build_md_package(ACC, "3OG7", POCKET)
    finally:
        md_package.remark465_lines = real
    assert result["ok"], result
    assert "No record" in result["tick_model_missing_residues"], result["tick_model_missing_residues"]
    assert "No record" in result["modelled_not_experimental"]
    assert result["checks"]["remark_465_lines"] == 0
    assert any("NO REMARK 465 block" in n for n in result["notes"]), result["notes"]
    readme = (Path(result["package_dir"]) / "README.md").read_text(encoding="utf-8")
    assert "no REMARK 465 block" in readme and "no run to tick" in readme, readme


def test_zero_occupancy_only_is_not_all_experimental():
    """An entry whose only rows are occupancy_flag 0 has no gap AND no all-experimental claim."""
    lines, source = md_package.remark465_lines(write_cif(
        MINIMAL_CIF.replace("1 1 Y 1 A MET 100 ? A MET 1\n", "")))
    assert lines == [] and source.startswith("mmCIF;") and "zero occupancy" in source, source
    prep_with_a_docking_run()
    real = md_package.remark465_lines
    try:
        md_package.remark465_lines = lambda path: ([], source)
        result = md_package.build_md_package(ACC, "3OG7", POCKET)
    finally:
        md_package.remark465_lines = real
    answer = result["modelled_not_experimental"]
    assert "every residue of the system is experimental" not in answer, answer
    assert "not experimental" in answer and "zero occupancy" in answer, answer


def test_build_refuses_what_it_cannot_write():
    prep_with_a_docking_run()
    assert md_package.build_md_package(ACC, "3OG7", POCKET, pose=0)["error"].startswith("pose 0")
    assert "outside 1..9" in md_package.build_md_package(ACC, "3OG7", POCKET, pose=99)["error"]
    bad = md_package.build_md_package(ACC, "3OG7", POCKET, conformer="crystal")
    assert "No run" in bad["error"], bad
    # a receptor that already uses chain L would silently merge with the ligand
    receptor = next(p for p in (dossier.DOSSIER_DIR / ACC / "prep" / "clean").glob("3OG7_receptor.pdb"))
    text = receptor.read_text(encoding="utf-8")
    receptor.write_text("\n".join(line[:21] + "L" + line[22:] if line.startswith("ATOM") else line
                                  for line in text.splitlines()) + "\n", encoding="utf-8")
    clash = md_package.build_md_package(ACC, "3OG7", POCKET)
    assert "already uses chain L" in clash["error"], clash


# ----------------------------------------------------------------- check_charmm_gui_output

def test_check_accepts_a_good_job():
    prep, delivered = test_build_writes_what_charmm_gui_needs()
    package = json.loads((delivered.parent / "package.json").read_text(encoding="utf-8"))
    tmp = Path(tempfile.mkdtemp())
    archive = make_archive(tmp, delivered, to_model=package["to_model"])
    result = md_package.check_charmm_gui_output(ACC, str(archive), "3OG7", POCKET)
    assert result["ok"] and result["passes"], result.get("failures") or result
    checks = result["checks"]
    assert checks["uploaded_file_is_the_delivered_one"]
    assert checks["ligand_atoms_in_namd_system"] == 52            # 51 + the lone pair on the Cl
    assert checks["ligand_pose_internal_change"]["max_angstrom"] <= md_package.LIGAND_DD_TOLERANCE
    assert checks["peptide_bonds_over_limit"] == []
    assert checks["modelled_runs_present"] == ["A 545-547", "A 597-614", "A 627-630",
                                               "B 601-614", "B 627-630"], checks
    assert result["not_checked"] == [], result["not_checked"]
    assert checks["ligand_sdf_is_the_delivered_one"]
    assert checks["pocket_residues"]["overlap"] >= md_package.MIN_POCKET_OVERLAP, checks["pocket_residues"]
    for segid, geometry in checks["protein_geometry"].items():
        rg = geometry["radius_of_gyration_angstrom"]
        assert abs(rg["returned"] - rg["delivered"]) <= md_package.RG_TOLERANCE, (segid, rg)
    assert checks["system"] == {"box_side_angstrom": 114.0, "ions": None, "waters": None,
                                "temperature_kelvin": 303.15, "timestep_fs": 2.0,
                                "cutoff_angstrom": 12.0, "production_steps": 500000,
                                "lone_pairs": True,
                                "cgenff_penalties": "7.000 ; charge penalty= 6.086"}, checks["system"]
    return prep, delivered, package


def test_check_rejects_every_real_defect():
    """One control per defect measured on the real jobs. Each MUST fail, and for its own reason."""
    prep, delivered, package = test_check_accepts_a_good_job()
    tmp = Path(tempfile.mkdtemp())
    to_model = package["to_model"]

    def failures(**kw):
        archive = make_archive(tmp, delivered, name=f"{len(list(tmp.glob('*.tgz')))}.tgz", **kw)
        result = md_package.check_charmm_gui_output(ACC, str(archive), "3OG7", POCKET)
        assert result["ok"], result
        assert not result["passes"], f"accepted a broken job: {kw}"
        return " | ".join(result["failures"])

    # 1. job 9401967361: the hetero chain is out by default and nothing says so
    assert "no ligand" in failures(with_ligand=False, to_model=to_model)
    # 2. jobs 9401967361 and 9417550068: gaps bridged, not filled
    bridged = failures(fill_gaps=False, to_model=to_model)
    assert "longer than" in bridged and "still absent" in bridged, bridged
    # 3. job 9418371365: the protein comes back expanded
    assert "distorted" in failures(expand=1.22, to_model=to_model)
    # 4. job 9418371365: ANISOU read as atoms in what was handed to Modeller
    assert "ANISOU" in failures(anisou_garbage=True, to_model=to_model)
    # 5. an archive from another job or another input
    other = delivered.with_name("other.pdb")
    other.write_text(delivered.read_text(encoding="ascii").replace("END\n", "REMARK\nEND\n"), encoding="ascii")
    assert "not the file build_md_package wrote" in failures(uploaded=other, to_model=to_model)
    # 6. a modelled loop built on top of the ligand
    assert "clashes with the ligand" in failures(to_model=to_model, clash_run=("A", 597))
    # 7. the ligand ends up outside its pocket. Its internal distances are untouched, so the
    #    shape check alone accepted this (reviewer, 2026-11-08)
    assert "not in its pocket" in failures(to_model=to_model, ligand_shift=(40.0, 0.0, 0.0))
    # 8. a mirrored ligand: every internal distance is identical and it stays in the pocket, so
    #    only a real protein-ligand distance catches it
    mirrored = failures(to_model=to_model, mirror=True)
    assert "overlaps the protein" in mirrored or "not in its pocket" in mirrored, mirrored
    # 9. the pose distorted: the control LIGAND_DD_TOLERANCE did not have
    assert "changed shape" in failures(to_model=to_model, ligand_scale=1.2)
    # 10. CHARMM-GUI parameterised another ligand file (RCSB, Mol2, Antechamber or OpenFF)
    other_sdf = delivered.with_name("other_ligand.sdf")
    other_sdf.write_text(delivered.with_name("3OG7_ligand.sdf").read_text(encoding="ascii")
                         .replace("LIG", "XXX", 1), encoding="ascii")
    assert "not the SDF build_md_package wrote" in failures(to_model=to_model, ligand_sdf=other_sdf)
    # 11. a run that was never a gap cannot have been modelled (it would pass by construction)
    package["to_model"] = to_model + [{"chain": "A", "first": 500, "last": 502, "length": 3,
                                       "kind": "internal", "absent_from_the_receptor": True}]
    (delivered.parent / "package.json").write_text(json.dumps(package), encoding="utf-8")
    assert "cannot have been modelled" in failures(to_model=to_model)


def test_check_refuses_before_measuring_anything():
    prep_with_a_docking_run()
    missing = md_package.check_charmm_gui_output(ACC, "no-such-file.tgz", "3OG7", POCKET)
    assert not missing["ok"] and "No file at" in missing["error"]
    tmp = Path(tempfile.mkdtemp())
    not_a_tar = tmp / "x.tgz"
    not_a_tar.write_bytes(b"not a tar at all")
    assert md_package.build_md_package(ACC, "3OG7", POCKET)["ok"]
    bad = md_package.check_charmm_gui_output(ACC, str(not_a_tar), "3OG7", POCKET)
    assert not bad["ok"] and "tar.gz" in bad["error"], bad
    # the package has to exist first
    shutil.rmtree(dossier.DOSSIER_DIR / ACC / "prep" / "md", ignore_errors=True)
    none = md_package.check_charmm_gui_output(ACC, str(not_a_tar), "3OG7", POCKET)
    assert "run build_md_package first" in none["error"], none


def test_thresholds_separate_the_real_jobs():
    """The thresholds must sit between the two real archives of 2026-11-08 (not in the repo).

    Without this, a later change could widen a threshold until the broken job passes.
    """
    for measured in REAL_GOOD["mean_dd"]:
        assert measured < md_package.MEAN_DD_TOLERANCE
    for measured in REAL_BROKEN["mean_dd"]:
        assert measured > md_package.MEAN_DD_TOLERANCE
    for measured in REAL_GOOD["rg_change"]:
        assert measured < md_package.RG_TOLERANCE
    for measured in REAL_BROKEN["rg_change"]:
        assert measured > md_package.RG_TOLERANCE
    assert REAL_GOOD["ligand_dd_max"] < md_package.LIGAND_DD_TOLERANCE
    assert REAL_GOOD["loop_to_ligand"] > md_package.MIN_HEAVY_CONTACT
    assert REAL_GOOD["min_contact"] > md_package.MIN_HEAVY_CONTACT
    assert min(REAL_BRIDGES) > md_package.MAX_PEPTIDE_BOND
    assert all(count > 0 for count in REAL_BROKEN["garbage_atoms"])


def main():
    test_zero_occupancy_residues_are_not_called_missing()
    test_a_missing_category_is_not_the_same_as_no_gaps()
    test_remark465_matches_the_rcsb_block()
    test_an_entry_without_unobserved_residues_gets_no_block()
    test_long_peptide_bonds_find_the_real_bridges()
    test_internal_distance_change_ignores_rotation_but_not_distortion()
    test_build_writes_what_charmm_gui_needs()
    test_build_never_offers_a_run_that_is_in_the_file()
    test_build_says_when_the_mmcif_has_no_record()
    test_zero_occupancy_only_is_not_all_experimental()
    test_build_refuses_what_it_cannot_write()
    test_check_accepts_a_good_job()
    test_check_rejects_every_real_defect()
    test_check_refuses_before_measuring_anything()
    test_thresholds_separate_the_real_jobs()
    print("ALL OK")


if __name__ == "__main__":
    main()
