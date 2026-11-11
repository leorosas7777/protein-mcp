"""Offline test for analyze_md (protein_mcp/files/md_analyze.py, Phase 3 step 6).

The real trajectory of job 20070766 is 45 MB and is not in the repository, so the fixture is a
synthetic system written with MDAnalysis itself. It is synthetic on purpose: every metric here has
an ANALYTIC answer, so a wrong number is a failed assertion rather than a judgement call. What the
fixture does copy from the real files are their defects, because those are what broke the naive
version:

- a ligand carrying `LP1`, a lone pair of mass 0, so that "heavy atoms" by name and by mass give
  different counts (5 against 4 here; 34 against 33 in the real system);
- a NAMD input with `minimize` before `run`, so that the first DCD frames are not dynamics;
- two chunks written with different `dcdfreq`, which is what happens the moment an analyst raises
  the sampling rate between chunks.

The motion is a RIGID TRANSLATION of the ligand by 0, 1, 2 and 3 A with the protein held still.
That is what makes the controls below possible, because each metric has one right answer:

    ligand RMSD        0, 1, 2, 3 A          (fit on the pocket; a fit on the LIGAND gives 0)
    COM distance       rises by exactly 3 A
    protein RMSD       0 at every frame      (the protein never moves)
    pocket RMSF        0 at every residue    (same reason)

Controls that must fail:

- `test_a_fit_on_the_ligand_would_hide_the_motion`: superposing on the ligand instead of the
  pocket returns ~0 A for a ligand that moved 3 A. The real code must not agree with it.
- `test_heavy_atoms_by_name_would_count_the_lone_pair`: the name-based selection is computed here
  and asserted to be WRONG, so the mass-based rule cannot quietly be replaced by it.
- `test_keeping_the_minimization_frames_changes_the_answer`: with the frames kept, the numbers
  move. A fixture where dropping them made no difference would prove nothing.
- `test_mixed_frame_spacing_is_refused`: two chunks of different dcdfreq must raise, not average.
- `test_a_ligand_that_leaves_is_not_called_stable`: the offered reading must change its wording
  when the ligand dissociates, or it is decoration.

Run: uv run python tests/test_md_analyze.py
"""
import json
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

from protein_mcp.files import md_analyze  # noqa: E402
from protein_mcp.files.md_analyze import MdAnalyzeError, analyze_md  # noqa: E402

# ---------------------------------------------------------------------------------------------
# The fixture
# ---------------------------------------------------------------------------------------------

# Three ALA residues in one segment, each N, H, CA, C, O. Residue 3 is deliberately parked far
# away so that the pocket comes out as 2 of 3 and the cutoff is doing real work.
PROTEIN = [
    ("PROA", 1, "ALA", [("N", "N", 14.007, (0.0, 0.0, 0.0)),
                        ("HN", "H", 1.008, (1.0, 0.0, 0.0)),
                        ("CA", "C", 12.011, (-1.0, 1.0, 0.0)),
                        ("C", "C", 12.011, (-1.0, 2.0, 0.0)),
                        ("O", "O", 15.999, (-2.0, 2.5, 0.0))]),
    ("PROA", 2, "ALA", [("N", "N", 14.007, (0.0, 4.0, 0.0)),
                        ("HN", "H", 1.008, (1.0, 4.0, 0.0)),
                        ("CA", "C", 12.011, (-1.0, 5.0, 0.0)),
                        ("C", "C", 12.011, (-1.0, 6.0, 0.0)),
                        ("O", "O", 15.999, (-2.0, 6.5, 0.0))]),
    ("PROA", 3, "ALA", [("N", "N", 14.007, (0.0, 40.0, 0.0)),
                        ("HN", "H", 1.008, (1.0, 40.0, 0.0)),
                        ("CA", "C", 12.011, (-1.0, 41.0, 0.0)),
                        ("C", "C", 12.011, (-1.0, 42.0, 0.0)),
                        ("O", "O", 15.999, (-2.0, 42.5, 0.0))]),
]

# The ligand. O1 sits 2.8 A from the HN of residue 1 along +x and collinear with N-HN, so the
# geometry is a hydrogen bond at frame 0 (2.8 A, 180 deg) and is broken by the translation.
# LP1 has mass 0: it is the lone pair that the name-based selection wrongly counts.
# H2 is a hydrogen ON A CARBON. It exists so that "a donor is the N or O the hydrogen is bonded
# to" is a claim the fixture can falsify: drop the `bonded` clause from the selection and H2
# becomes an eligible donor, which the real ligand would repeat 16 times over.
LIGAND = [("N1", "N", 14.007, (3.8, 0.0, 0.0)),
          ("H1", "H", 1.008, (4.3, 0.9, 0.0)),
          ("C1", "C", 12.011, (4.8, -0.8, 0.0)),
          ("CL", "CL", 35.453, (6.2, -0.4, 0.0)),
          ("LP1", "LP", 0.0, (6.9, -0.2, 0.0)),
          ("O1", "O", 15.999, (3.8, 0.0, 0.0)),
          ("H2", "H", 1.008, (5.0, -1.9, 0.0))]
LIGAND[5] = ("O1", "O", 15.999, (2.8, 0.0, 0.0))

# bonds inside the ligand, as indices into LIGAND: N1-H1, N1-C1, C1-CL, CL-LP1, C1-O1, C1-H2.
LIGAND_BONDS = [(0, 1), (0, 2), (2, 3), (3, 4), (2, 5), (2, 6)]

SHIFTS_DYNAMICS = [0.0, 1.0, 2.0, 3.0]
SHIFTS_MINIMIZATION = [-15.0, -8.0]       # absurd on purpose: keeping them must change the numbers
BOX = [60.0, 60.0, 60.0, 90.0, 90.0, 90.0]

NAMD_INPUT = """\
set temp                303.15;
dcdfreq                 1000;
timestep                2.0;
minimize                2000
numsteps                90000000
run                     4000
"""


def _universe(extra_massless=None):
    """A Universe with the topology attributes the analyser needs, and nothing else.

    `extra_massless` adds one more zero-mass particle to the ligand under a name of the caller's
    choosing, which is how the "heavy atoms by mass, never by a name rule" claim is tested against
    a name that no `not name LP*` pattern would catch.
    """
    import MDAnalysis as mda

    ligand = list(LIGAND) + ([extra_massless] if extra_massless else [])
    names, types, masses, resids, resnames, segids, positions = [], [], [], [], [], [], []
    residue_of_atom = []
    res_index = 0
    for segid, resid, resname, atoms in PROTEIN:
        for name, typ, mass, xyz in atoms:
            names.append(name); types.append(typ); masses.append(mass)
            positions.append(xyz); residue_of_atom.append(res_index)
        resids.append(resid); resnames.append(resname); segids.append(segid)
        res_index += 1
    for name, typ, mass, xyz in ligand:
        names.append(name); types.append(typ); masses.append(mass)
        positions.append(xyz); residue_of_atom.append(res_index)
    resids.append(1); resnames.append("LIG"); segids.append("HETA")

    n_atoms = len(names)
    u = mda.Universe.empty(n_atoms, n_residues=len(resids), n_segments=len(set(segids)),
                           atom_resindex=np.array(residue_of_atom),
                           residue_segindex=np.array([0, 0, 0, 1]),
                           trajectory=True)
    u.add_TopologyAttr("name", names)
    u.add_TopologyAttr("type", types)
    u.add_TopologyAttr("element", [t if t not in ("LP", "DUM") else "" for t in types])
    u.add_TopologyAttr("mass", masses)
    u.add_TopologyAttr("resid", resids)
    u.add_TopologyAttr("resname", resnames)
    u.add_TopologyAttr("segid", sorted(set(segids), key=segids.index))
    u.atoms.positions = np.array(positions, dtype=np.float32)
    # bonds: the backbone of each residue plus the ligand, so donors are read and not guessed
    bonds = []
    for r in range(3):
        base = r * 5
        bonds += [(base, base + 1), (base, base + 2), (base + 2, base + 3), (base + 3, base + 4)]
    lig0 = len(PROTEIN) * 5
    bonds += [(lig0 + a, lig0 + b) for a, b in LIGAND_BONDS]
    u.add_TopologyAttr("bonds", bonds)
    return u


def write_fixture(folder: Path, *, shifts=None, dt=2.0, name="step5_1.dcd",
                  with_minimization=True, inp=True, box=None) -> Path:
    """Write a PSF-equivalent topology, a DCD of the given shifts, and the NAMD input.

    The topology goes out as a PDB rather than a PSF: MDAnalysis cannot write a PSF, and what the
    analyser needs from the topology (names, masses, bonds, segids) a PDB carries. `find_inputs`
    looks for `*.psf`, so the test points the analyser at the files directly where it must, and
    the PSF path is exercised against the real file in the live check instead.
    """
    import MDAnalysis as mda

    folder.mkdir(parents=True, exist_ok=True)
    u = _universe()
    lig = u.select_atoms("resname LIG")
    base = lig.positions.copy()
    frames = (list(SHIFTS_MINIMIZATION) if with_minimization else []) + list(
        SHIFTS_DYNAMICS if shifts is None else shifts)
    path = folder / name
    with mda.Writer(str(path), n_atoms=len(u.atoms), dt=dt) as writer:
        for shift in frames:
            lig.positions = base + np.array([shift, 0.0, 0.0], dtype=np.float32)
            # `box=None` keeps the default rectangular cell; a caller can hand a triclinic one, or
            # all zeros, which is how a DCD written without unit-cell information arrives.
            u.trajectory.ts.dimensions = np.array(BOX if box is None else box, dtype=np.float32)
            writer.write(u.atoms)
    if inp:
        (folder / "step5_production.inp").write_text(
            NAMD_INPUT if with_minimization else NAMD_INPUT.replace("minimize                2000\n", ""))
    return path


# ---------------------------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------------------------

def test_the_fixture_is_what_it_claims():
    """The synthetic system has the defects the real one has, or the rest proves nothing."""
    u = _universe()
    lig = u.select_atoms("resname LIG")
    assert len(lig) == len(LIGAND), len(lig)
    # a hydrogen on carbon, which must never become a donor
    assert "H2" in set(lig.names)
    h2 = lig.select_atoms("name H2")[0]
    assert all(b.mass > 11.0 for a in (h2,) for b in a.bonded_atoms), "H2 must hang off a carbon"
    by_mass = lig.select_atoms(f"mass {md_analyze.HEAVY_MASS_MIN} to 100000")
    by_name = lig.select_atoms("not name H*")
    assert len(by_mass) == 4, len(by_mass)
    assert len(by_name) == 5, len(by_name)
    assert len(by_mass) != len(by_name), "the fixture must reproduce the lone-pair trap"
    assert "LP1" in set(by_name.names) and "LP1" not in set(by_mass.names)
    assert u.select_atoms("protein").residues.n_residues == 3
    print("ok  the fixture carries the lone pair and the name/mass disagreement")


def test_minimization_frames_are_read_from_the_input():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        mini = md_analyze.minimization_frames(list(folder.glob("*.inp")))
    assert mini["frames"] == 2, mini
    assert "minimize 2000" in " ".join(mini["evidence"])
    assert "dcdfreq 1000" in " ".join(mini["evidence"])
    assert mini["read_from"] == "step5_production.inp"
    print("ok  minimize 2000 / dcdfreq 1000 is read as 2 frames, with its evidence")


def test_an_input_without_minimization_claims_none():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        mini = md_analyze.minimization_frames(list(folder.glob("*.inp")))
    assert mini["frames"] == 0 and mini["read_from"] is None, mini
    print("ok  a production input declares no minimization frames")


def test_heavy_atoms_by_name_would_count_the_lone_pair():
    """The control for the first real defect: the name rule is asserted to be WRONG here."""
    u = _universe()
    chosen = md_analyze.select_ligand(u, "LIG")
    assert chosen["selection"] == f"4 heavy of {len(LIGAND)} atoms in resname LIG", chosen["selection"]
    assert "LP1" in chosen["dropped"], chosen["dropped"]
    wrong = len(u.select_atoms("resname LIG and not name H*"))
    assert wrong == 5 and wrong != chosen["atoms"].n_atoms, (wrong, chosen["atoms"].n_atoms)
    print(f"ok  heavy atoms by mass = 4, by name = {wrong} (the lone pair); the code uses mass")


def test_the_pocket_is_measured_not_named():
    u = _universe()
    lig = md_analyze.select_ligand(u, "LIG")
    pocket = md_analyze.select_pocket(u, lig["atoms"], md_analyze.POCKET_CUTOFF, 0)
    assert pocket["selection"] == "2 of 3 protein residues", pocket["selection"]
    assert pocket["residues"] == ["PROA:ALA1", "PROA:ALA2"], pocket["residues"]
    assert "PROA:ALA3" not in pocket["residues"], "the far residue must stay out"
    assert pocket["segments"] == {"PROA": 2}, pocket["segments"]
    assert md_analyze.POCKET_CUTOFF == 5.0, md_analyze.POCKET_CUTOFF   # a literal, on purpose
    assert "5.0 A" in pocket["criterion"], pocket["criterion"]
    print("ok  the pocket is 2 of 3 residues by distance, and the list comes back")


def test_a_ligand_that_is_not_there_is_refused():
    u = _universe()
    try:
        md_analyze.select_ligand(u, "ZZZ")
    except MdAnalyzeError as exc:
        assert "ZZZ" in str(exc) and "LIG" in str(exc), str(exc)
        print("ok  an absent resname is refused and the real ones are listed")
        return
    raise AssertionError("a ligand that is not in the system must be refused")


def test_the_metrics_match_the_analytic_answer():
    """A rigid 0-1-2-3 A translation of the ligand, with the protein held still."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        top = _run_topology(folder)
        r = _run(folder, top)

    assert r["trajectory"]["frames_analysed"] == 4, r["trajectory"]
    assert r["trajectory"]["frames_dropped"] == 2, r["trajectory"]

    lig = r["ligand_rmsd"]["values"]
    assert np.allclose(lig, [0.0, 1.0, 2.0, 3.0], atol=2e-3), lig

    # The COM distance is the NORM of a vector, so a 3 A push along x does not add 3 A to it
    # unless the vector points along x. The expected series is computed from the fixture's own
    # geometry rather than guessed, which is the whole point of a synthetic system.
    expected = _expected_com_distances()
    com = r["com_distance"]["values"]
    assert np.allclose(com, expected, atol=2e-3), (com, expected)
    assert abs(r["com_distance"]["range"] - (max(expected) - min(expected))) < 2e-3
    assert max(r["protein_rmsd"]["values"]) < 1e-3, r["protein_rmsd"]["values"]
    assert all(e["rmsf_angstrom"] < 1e-3 for e in r["rmsf"]["per_residue"]), r["rmsf"]
    print("ok  ligand RMSD 0/1/2/3 A, COM matches the analytic norm, protein RMSD 0, RMSF 0: the analytic answer")


def test_a_fit_on_the_ligand_would_hide_the_motion():
    """The control that must fail: superposing on the ligand returns ~0 for a 3 A move."""
    from MDAnalysis.analysis import align
    import MDAnalysis as mda

    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        dcd = write_fixture(folder)
        top = _run_topology(folder)
        u = mda.Universe(str(top), str(dcd))
        try:
            lig = u.select_atoms(f"resname LIG and mass {md_analyze.HEAVY_MASS_MIN} to 100000")
            u.trajectory[2]
            ref = lig.positions - lig.center_of_mass()
            wrong = []
            for _ in u.trajectory[2:]:
                mob = lig.positions - lig.center_of_mass()
                _, rms = align.rotation_matrix(mob, ref, weights=lig.masses)
                wrong.append(float(rms))
        finally:
            u.trajectory.close()      # Windows keeps a lock on an open DCD
        r = _run(folder, top)

    assert max(wrong) < 1e-3, f"a fit on the ligand should hide a rigid move, got {wrong}"
    right = r["ligand_rmsd"]["values"]
    assert max(right) > 2.9, right
    assert "never the ligand itself" in r["ligand_rmsd"]["criterion"]
    print(f"ok  fitting on the ligand gives {max(wrong):.4f} A for a 3 A move; the code gives "
          f"{max(right):.3f} A")


def test_keeping_the_minimization_frames_changes_the_answer():
    """If dropping them made no difference, the fixture would not be testing anything."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        top = _run_topology(folder)
        dropped = _run(folder, top)
        # reference_frame=2 keeps the pocket defined from a dynamics frame: in this fixture the
        # minimization frames park the ligand 15 A away on purpose, and defining a pocket from
        # one of them is refused (see test_a_reference_frame_without_a_pocket_is_refused).
        kept = _run(folder, top, skip_frames=0, reference_frame=2)

    assert dropped["trajectory"]["frames_analysed"] == 4
    assert kept["trajectory"]["frames_analysed"] == 6
    # The two frames must move the numbers, or this control is empty. The figure to beat is not
    # arbitrary: it is how far the fixture parks the ligand in the first minimization frame,
    # measured from the reference frame, so keeping them MUST reach it and dropping them cannot
    # come anywhere near it.
    # frame 2 of the whole trajectory is the FIRST dynamics frame: the file holds
    # len(SHIFTS_MINIMIZATION) minimization frames before it.
    parked = abs(SHIFTS_MINIMIZATION[0] - SHIFTS_DYNAMICS[2 - len(SHIFTS_MINIMIZATION)])
    assert abs(kept["ligand_rmsd"]["max"] - parked) < 1e-6, (
        kept["ligand_rmsd"]["max"], parked)
    assert dropped["ligand_rmsd"]["max"] < parked / 2, (
        dropped["ligand_rmsd"]["max"], parked)
    assert kept["com_distance"]["range"] > 2 * dropped["com_distance"]["range"], (
        kept["com_distance"]["range"], dropped["com_distance"]["range"])
    assert "They were dropped." in dropped["trajectory"]["minimization_note"]
    assert "KEPT" in kept["trajectory"]["minimization_note"]
    assert "fiction" in kept["trajectory"]["minimization_note"]
    print(f"ok  dropping the 2 minimization frames changes the COM range from "
          f"{kept['com_distance']['range']} to {dropped['com_distance']['range']} A")


def test_a_reference_frame_without_a_pocket_is_refused():
    """A reference frame where the ligand is nowhere near the protein must fail loudly.

    Inventing a pocket from a frame that has none would put a residue list in the report that no
    measurement supports, which is the failure this whole module is built to avoid.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        top = _run_topology(folder)
        r = _run(folder, top, skip_frames=0, reference_frame=0)
    assert r["ok"] is False, r
    assert "within 5.0 A" in r["error"] and "frame 0" in r["error"], r["error"]
    print("ok  a reference frame with no pocket is refused instead of guessed")


def test_mixed_frame_spacing_is_refused():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, name="step5_1.dcd", dt=2.0)
        write_fixture(folder, name="step5_2.dcd", dt=10.0, inp=False)
        top = _run_topology(folder)
        r = _run(folder, top)
    assert r["ok"] is False, r
    assert "step5_1.dcd" in r["error"] and "step5_2.dcd" in r["error"], r["error"]
    assert "dcdfreq" in r["error"], r["error"]
    assert r["spacings"] == [2.0, 10.0], r.get("spacings")
    print("ok  two chunks of different frame spacing are refused by name, as ok: false")


def test_the_chunks_are_ordered_by_number_not_by_text():
    paths = [Path(f"step5_{i}.dcd") for i in (1, 2, 10, 11, 3)]
    ordered = [p.name for p in sorted(paths, key=md_analyze._chunk_order)]
    assert ordered == ["step5_1.dcd", "step5_2.dcd", "step5_3.dcd",
                       "step5_10.dcd", "step5_11.dcd"], ordered
    print("ok  step5_10.dcd comes after step5_2.dcd, not after step5_1.dcd")


def test_a_folder_without_a_trajectory_says_so():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        folder.mkdir(exist_ok=True)
        (folder / "step3_input.psf").write_text("not a psf")
        try:
            md_analyze.find_inputs(folder)
        except MdAnalyzeError as exc:
            assert "No .dcd" in str(exc) and "still running" in str(exc), str(exc)
            print("ok  a folder with no .dcd is told so, with the running job as the likely cause")
            return
    raise AssertionError("a folder without a trajectory must be refused")


def test_the_hydrogen_bond_occupancy_is_a_sample_of_frames():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        top = _run_topology(folder)
        r = _run(folder, top)
    hb = r["hbonds"]
    assert "not the share of elapsed time" in hb["occupancy_note"], hb["occupancy_note"]
    # LITERALS, not str(CONSTANT): comparing a constant with a string the module built from that
    # same constant passes whatever the constant is. Changing 150 deg to 90 deg left this test
    # green. These two numbers are the geometric definition the module declares, so they are
    # written out here and either side moving must break this.
    assert md_analyze.HB_DA_CUTOFF == 3.5, md_analyze.HB_DA_CUTOFF
    assert md_analyze.HB_ANGLE_CUTOFF == 150.0, md_analyze.HB_ANGLE_CUTOFF
    assert "3.5 A" in hb["criterion"], hb["criterion"]
    assert "150.0 deg" in hb["criterion"], hb["criterion"]
    assert "nothing is guessed from distances" in hb["connectivity"]
    for pair in hb["pairs"]:
        assert pair["frames_analysed"] == 4, pair
        assert 0 < pair["occupancy_sampled"] <= 1, pair
        assert pair["frames_present"] <= 4, pair
    print(f"ok  {hb['pairs_found']} hydrogen-bond pair(s), occupancy declared as a frame sample")


def test_a_ligand_that_leaves_is_not_called_stable():
    """The reading must change when the ligand dissociates, or it is decoration."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, shifts=[0.0, 6.0, 12.0, 20.0])
        top = _run_topology(folder)
        gone = _run(folder, top)
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        # Gentler than the analytic fixture on purpose: 0-1-2-3 A already crosses the 2 A the
        # reading compares against, so a ligand that HOLDS has to move less than that.
        write_fixture(folder, shifts=[0.0, 0.3, 0.5, 0.4])
        top = _run_topology(folder)
        held = _run(folder, top)

    assert any("held within" in n for n in held["reading"]["offered"]), held["reading"]["offered"]
    assert any("moved beyond" in n for n in gone["reading"]["offered"]), gone["reading"]["offered"]
    assert any("may be leaving the site" in n for n in gone["reading"]["offered"])
    assert gone["ligand_rmsd"]["max"] > 10 > held["ligand_rmsd"]["max"]
    print("ok  a ligand that walks 20 A is not described with the same words as one that holds")


def test_the_reading_is_marked_as_offered_and_names_its_thresholds():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        top = _run_topology(folder)
        r = _run(folder, top)
    reading = r["reading"]
    assert reading["status"] == "READING OFFERED, NOT MEASURED"
    assert "Apply your own judgement." in reading["disclaimer"]
    assert "conventions and not measurements" in reading["disclaimer"]
    # the thresholds must point at their pre-registration and must NOT claim a source that was
    # never consulted: an earlier version attributed them to "the literature" with no citation
    assert "section 7.2" in reading["thresholds_source"], reading["thresholds_source"]
    assert "No literature is cited" in reading["thresholds_source"]
    assert "literature" not in reading["disclaimer"]
    for name, value in reading["thresholds_used"].items():
        assert any(str(value) in n for n in reading["offered"]), f"{name} is not shown to the reader"
    # no metric block may carry a verdict of its own
    for key in ("protein_rmsd", "ligand_rmsd", "com_distance", "rmsf", "hbonds"):
        text = json.dumps(r[key]).lower()
        for word in ("stable", "unstable", "good", "favourable", "poor", "pass", "fail"):
            assert word not in text, f"{key} carries a verdict: {word!r}"
    print("ok  the reading is marked as offered, shows every threshold, and no metric judges")


def test_nothing_claims_an_affinity():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        top = _run_topology(folder)
        r = _run(folder, top)
    joined = " ".join(r["not_measured"]).lower()
    assert "affinity" in joined and "free energy" in joined, r["not_measured"]
    assert "apo control" in joined, r["not_measured"]
    body = json.dumps({k: v for k, v in r.items() if k != "not_measured"}).lower()
    for word in ("kcal", "kd", "ic50", "binding energy", "affinity"):
        assert word not in body, f"the result mentions {word!r} without measuring it"
    print("ok  no affinity, no free energy, and the gaps are listed")


def test_the_result_survives_json():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        top = _run_topology(folder)
        r = _run(folder, top)
        text = json.dumps(r)
        light = _run(folder, top, include_series=False)
    assert len(text) > 500
    assert "values" not in json.dumps(light["ligand_rmsd"])
    assert light["ligand_rmsd"]["mean"] == r["ligand_rmsd"]["mean"]
    print("ok  the result is JSON, and include_series=False drops the series but not the summary")


# ---------------------------------------------------------------------------------------------
# The mutations that the first round of tests let through
#
# Each test below exists because a specific sabotage of md_analyze.py passed all 18 of the
# original tests. The fixture there moved only the ligand, by translation, in a rectangular box
# with no wrapping, so everything the module does about ROTATION, about PERIODIC IMAGES and about
# the SHAPE of a curve was unmeasured, and the one occupancy figure was never pinned to a number.
# ---------------------------------------------------------------------------------------------

def _rotation(degrees: float):
    a = np.radians(degrees)
    return np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])


def write_rotating_fixture(folder: Path, degrees=(0.0, 10.0, 25.0, 40.0)) -> Path:
    """The WHOLE system turns rigidly: protein and ligand together, nothing moving relative.

    A superposition that fits translation only cannot undo a rotation, so the ligand RMSD against
    a reference frame comes out large. A correct fit on the pocket backbone returns it to ~0.
    This is what `rot = np.eye(3)` quietly broke while every earlier test still passed.
    """
    import MDAnalysis as mda

    folder.mkdir(parents=True, exist_ok=True)
    u = _universe()
    base = u.atoms.positions.copy()
    centre = base.mean(axis=0)
    path = folder / "step5_1.dcd"
    with mda.Writer(str(path), n_atoms=len(u.atoms), dt=2.0) as writer:
        for deg in degrees:
            u.atoms.positions = (base - centre) @ _rotation(deg).T + centre
            u.trajectory.ts.dimensions = np.array(BOX, dtype=np.float32)
            writer.write(u.atoms)
    (folder / "step5_production.inp").write_text(
        NAMD_INPUT.replace("minimize                2000\n", ""))
    return path


def write_wrapped_fixture(folder: Path) -> Path:
    """The ligand jumps a FULL BOX LENGTH in two frames and moves 0.2 A in the others.

    This is `wrapAll on`, which CHARMM-GUI writes into both real inputs: NAMD puts a molecule
    that drifts past a face of the box back into the central cell. Nothing physical happens, and
    a raw RMSD reads it as the ligand travelling 60 A.
    """
    import MDAnalysis as mda

    folder.mkdir(parents=True, exist_ok=True)
    u = _universe()
    lig = u.select_atoms("resname LIG")
    base = lig.positions.copy()
    path = folder / "step5_1.dcd"
    with mda.Writer(str(path), n_atoms=len(u.atoms), dt=2.0) as writer:
        for shift in (0.0, 0.2, BOX[0], BOX[0] + 0.2):
            lig.positions = base + np.array([shift, 0.0, 0.0], dtype=np.float32)
            u.trajectory.ts.dimensions = np.array(BOX, dtype=np.float32)
            writer.write(u.atoms)
    (folder / "step5_production.inp").write_text(
        NAMD_INPUT.replace("minimize                2000\n", ""))
    return path


def test_a_rotating_protein_needs_the_rotation():
    """M1: `rot = np.eye(3)` survived all 18 original tests. It must not survive this one."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_rotating_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None)
    assert r["ok"], r
    lig = r["ligand_rmsd"]["values"]
    assert max(lig) < 0.05, (
        "the system turns rigidly, so a fit on the pocket must return the ligand to its "
        f"reference: got {lig}")
    prot = r["protein_rmsd"]["values"]
    assert max(prot) < 0.05, prot
    print(f"ok  a 40 deg rigid rotation leaves the ligand RMSD at {max(lig):.4f} A "
          "(translation-only fitting cannot do this)")


def test_a_wrapped_ligand_is_not_called_motion():
    """M2: deleting the minimum-image correction survived all 18 original tests."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_wrapped_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None)
    assert r["ok"], r
    lig, com = r["ligand_rmsd"], r["com_distance"]
    assert com["pbc_corrected_frames"] == 2, com["pbc_corrected_frames"]
    assert max(lig["values"]) < 0.5, (
        "a ligand re-imaged by wrapAll has not moved, so the RMSD must stay small: "
        f"got {lig['values']}")
    assert com["range"] < 0.5, com["values"]
    # and the reading must not contradict itself about the same frames
    offered = " ".join(r["reading"]["offered"])
    assert "moved beyond" not in offered, offered
    assert "re-imaged" in offered, offered
    print(f"ok  a full-box jump gives RMSD {max(lig['values']):.4f} A and is reported as "
          "re-imaging, not as motion")


def test_a_monotonic_drift_is_not_a_plateau():
    """M3: `reading()` always saying "consistent with a plateau" survived all 18 tests.

    The old statistic compared the last quarter with the mean of the WHOLE run, which for a
    straight ramp of height h only ever reaches 0.375 h: a backbone climbing steadily by 1.2 A
    passed a 0.5 A test. The statistic is now last quarter against FIRST quarter.
    """
    n = 24
    for height, expected in ((0.1, "plateau"), (0.8, "drift"), (1.2, "drift"), (2.0, "drift")):
        series = md_analyze._series([height * i / (n - 1) for i in range(n)],
                                    list(range(n)), unit="angstrom")
        fake = {"protein_rmsd": series, "ligand_rmsd": series,
                "com_distance": {**series, "pbc_corrected_frames": 0,
                                 "non_orthorhombic_frames": 0},
                "hbonds": {"pairs": [], "pairs_found": 0}}
        first = md_analyze.reading(fake)["offered"][0]
        said_plateau = "consistent with a plateau" in first
        assert said_plateau == (expected == "plateau"), (height, expected, first)
    print("ok  a monotonic ramp of 0.8, 1.2 or 2.0 A is called a drift; 0.1 A is a plateau")


def test_heavy_atoms_are_chosen_by_mass_not_by_any_name_rule():
    """M4: `not name H* and not name LP*` survived, because the fixture's only massless particle
    was called LP1. A massless particle under any other name must still be excluded."""
    u = _universe(extra_massless=("DUM", "DUM", 0.0, (6.5, 0.6, 0.0)))
    chosen = md_analyze.select_ligand(u, "LIG")
    assert chosen["atoms"].n_atoms == 4, chosen["selection"]
    assert "DUM" in chosen["dropped"] and "LP1" in chosen["dropped"], chosen["dropped"]
    by_name = len(u.select_atoms("resname LIG and not name H* and not name LP*"))
    assert by_name == 5 and by_name != chosen["atoms"].n_atoms, by_name
    print(f"ok  a massless particle named DUM is dropped by mass; any name rule keeps {by_name}")


def test_mixed_frame_spacing_is_caught_in_a_middle_chunk():
    """M5: comparing only the first and the last DCD survived all 18 tests."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, name="step5_1.dcd", dt=2.0, with_minimization=False)
        write_fixture(folder, name="step5_2.dcd", dt=10.0, with_minimization=False, inp=False)
        write_fixture(folder, name="step5_3.dcd", dt=2.0, with_minimization=False, inp=False)
        _run_topology(folder)
        r = _run(folder, None)
    assert r["ok"] is False, r
    assert "step5_2.dcd" in r["error"], r["error"]
    print("ok  a middle chunk with a different spacing is caught, not skipped over")


def test_minimization_frames_round_down():
    """M7: rounding up survived, because 2000/1000 is exact. A partial frame is not a frame."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        folder.mkdir(exist_ok=True)
        (folder / "step5_production.inp").write_text(
            NAMD_INPUT.replace("minimize                2000", "minimize                2500"))
        mini = md_analyze.minimization_frames(list(folder.glob("*.inp")))
    assert mini["frames"] == 2, mini
    print("ok  minimize 2500 / dcdfreq 1000 is 2 whole frames, not 3")


def test_the_occupancy_denominator_is_the_analysed_frames():
    """M10: dividing by `frames_analysed + 2` survived, because no test fixed an occupancy."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None)
    pairs = r["hbonds"]["pairs"]
    assert pairs, "the fixture must produce at least one hydrogen bond"
    pair = pairs[0]
    # the fixture's geometry holds in exactly one of the four analysed frames: 2.8 A at shift 0,
    # then 3.8, 4.8 and 5.8 A, all past the 3.5 A cutoff
    assert pair["frames_present"] == 1, pair
    assert pair["frames_analysed"] == 4, pair
    assert pair["occupancy_sampled"] == 0.25, pair
    print("ok  occupancy is exactly 1/4 = 0.25: the denominator is the analysed frames")


def test_the_equilibration_is_not_chained_onto_the_production_chunks():
    """The real folder holds both, at different dcdfreq, and the equilibration ran FIRST."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, name="step5_1.dcd", dt=2.0, with_minimization=False)
        write_fixture(folder, name="step5_2.dcd", dt=2.0, with_minimization=False, inp=False)
        write_fixture(folder, name="step4_equilibration.dcd", dt=10.0, inp=False)
        _run_topology(folder)
        found = md_analyze.find_inputs(folder)
        r = _run(folder, None)
    assert [p.name for p in found["dcds"]] == ["step5_1.dcd", "step5_2.dcd"], found["dcds"]
    assert found["left_out"] == ["step4_equilibration.dcd"], found["left_out"]
    assert r["ok"], r
    # the stage is `step5`, read from the file names, and NOT the word "production", which this
    # tool is no longer allowed to assert about a trajectory it did not generate
    assert r["trajectory"]["stage"] == "step5", r["trajectory"]["stage"]
    assert r["trajectory"]["stages_present"] == {
        "step4": ["step4_equilibration.dcd"],
        "step5": ["step5_1.dcd", "step5_2.dcd"]}, r["trajectory"]["stages_present"]
    assert "step4_equilibration.dcd (stage step4)" in r["trajectory"]["stage_note"]
    assert "the latest of 2 present" in r["trajectory"]["stage_selection"], \
        r["trajectory"]["stage_selection"]
    print("ok  the equilibration is left out as stage step4 and named in the reply")


def test_a_production_run_written_in_one_piece_is_not_called_another_stage():
    """THE DEFECT THIS RULE WAS CHANGED FOR, and it is the file the user's own run produces.

    The old rule asked whether a name ended in `_<number>`, which only build_md_job's own chunks
    do. A 10 ns production run written in one piece as `step5_production.dcd` therefore came back
    with `stage: "other"` and a note reading "these are not production: say which stage they are
    before quoting any number from them" - about production data. The numbers were right and the
    label was false, which is the worst of the two ways to be wrong here.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, name="step5_production.dcd", with_minimization=False, inp=False)
        (folder / "step5_production.inp").write_text(
            NAMD_INPUT.replace("minimize                2000\n", ""))
        _run_topology(folder)
        found = md_analyze.find_inputs(folder)
        r = _run(folder, None)
    assert [p.name for p in found["dcds"]] == ["step5_production.dcd"], found["dcds"]
    assert r["ok"], r
    assert r["trajectory"]["stage"] == "step5", r["trajectory"]["stage"]
    assert r["trajectory"]["stage_inputs"] == ["step5_production.inp"], r["trajectory"]
    assert r["trajectory"]["left_out"] == [], r["trajectory"]["left_out"]
    # the wording the old rule produced must be gone
    note = r["trajectory"]["stage_note"]
    assert "are not production" not in note, note
    assert "_<number>" not in note, note
    # and the tool must say where the word "production" comes from, instead of claiming it
    assert "comes from THAT FILE'S NAME" in note, note
    assert "step5_production.inp" in note, note
    assert r["trajectory"]["frames_analysed"] == 4, r["trajectory"]["frames_note"]
    print("ok  a one-piece step5_production.dcd is stage step5, named by its own input, with no "
          "claim that it is not production")


def test_a_named_stage_is_analysed_and_an_absent_one_is_refused():
    """`stage` lets the analyst pick, instead of the tool guessing from a file name."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, name="step5_production.dcd", with_minimization=False, inp=False)
        write_fixture(folder, name="step4_equilibration.dcd", dt=10.0, inp=False)
        (folder / "step5_production.inp").write_text(
            NAMD_INPUT.replace("minimize                2000\n", ""))
        (folder / "step4_equilibration.inp").write_text(NAMD_INPUT)
        _run_topology(folder)
        latest = _run(folder, None)
        picked = _run(folder, None, stage="step4")
        absent = _run(folder, None, stage="step9")
    # by default the LATEST stage, which is the declared rule
    assert latest["trajectory"]["stage"] == "step5", latest["trajectory"]["stage"]
    # and naming step4 analyses the equilibration instead, with its own minimization frames read
    assert picked["trajectory"]["stage"] == "step4", picked["trajectory"]["stage"]
    assert picked["trajectory"]["stage_inputs"] == ["step4_equilibration.inp"], picked["trajectory"]
    assert picked["trajectory"]["minimization"]["frames"] == 2, picked["trajectory"]["minimization"]
    assert "which you named" in picked["trajectory"]["stage_selection"], picked["trajectory"]
    assert picked["trajectory"]["left_out"] == ["step5_production.dcd"], picked["trajectory"]
    # a stage that is not there is a refusal that lists what is
    assert absent["ok"] is False, absent
    assert "No trajectory belongs to stage 'step9'" in absent["error"], absent["error"]
    assert set(absent["stages_present"]) == {"step4", "step5"}, absent["stages_present"]
    print("ok  the latest stage is the default, `stage` picks another, and an absent one is "
          "refused with the list")


def test_an_unnumbered_trajectory_sorts_first():
    paths = [Path(n) for n in ("step5_2.dcd", "step4_equilibration.dcd", "step5_1.dcd")]
    assert [p.name for p in sorted(paths, key=md_analyze._chunk_order)] == [
        "step4_equilibration.dcd", "step5_1.dcd", "step5_2.dcd"]
    print("ok  an unnumbered trajectory sorts before the chunks, not after them")


def test_a_refusal_comes_back_as_ok_false():
    """Every other tool in the project answers a refusal with ok: false, and so must this one."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None, skip_frames=0, reference_frame=0)
    assert r["ok"] is False, r
    assert "within 5.0 A" in r["error"], r["error"]
    assert "folder" in r
    print("ok  a refusal is ok: false with its reason, not an exception")


def test_the_time_axis_starts_at_the_first_analysed_frame():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None)
    times = r["ligand_rmsd"]["times_ps"]
    assert times[0] == 0.0, times
    assert times == [0.0, 2.0, 4.0, 6.0], times
    assert "FROM THE FIRST ANALYSED FRAME" in r["ligand_rmsd"]["times_note"]
    print("ok  the time axis starts at 0 at the first analysed frame, not at the file's origin")


def test_the_hbond_criterion_does_not_mention_a_charge_and_says_what_was_eligible():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None)
    hb = r["hbonds"]
    assert "No partial charge" in hb["criterion"], hb["criterion"]
    assert "No partial charge is used" in hb["connectivity"], hb["connectivity"]
    assert hb["eligible"]["acceptors"] > 0 and hb["eligible"]["hydrogens"] > 0, hb["eligible"]
    # The eligible set is mass-based, so it is every N and O of the protein and the ligand and
    # nothing else: 3 ALA x (N, O) = 6, plus the ligand's N1 and O1 = 8. A charge-based rule drops
    # some of those, which on the real ligand hid both sulfonamide oxygens and the carbonyl.
    assert hb["eligible"]["acceptors"] == 8, hb["eligible"]
    # hydrogens: the 3 backbone HN plus the ligand's H1, all bonded to an N or an O. The ligand
    # also carries hydrogens on carbon in the real system, and those must never become donors.
    assert hb["eligible"]["hydrogens"] == 4, hb["eligible"]
    print(f"ok  the hydrogen-bond criterion names no charge; {hb['eligible']['acceptors']} "
          "acceptors were eligible before any geometry")


def test_the_thresholds_match_their_pre_registration():
    """The four thresholds in the code must be the four written down in ESTADO_PROYECTO.md 7.2.

    Pre-registering a criterion and then quietly changing it in the code is worse than not
    pre-registering it at all, because the document would still vouch for the old number. This
    test fails if either side drifts.
    """
    doc = Path(__file__).resolve().parents[1] / "ESTADO_PROYECTO.md"
    if not doc.is_file():                      # the test must not depend on the document existing
        print("ok  (skipped: ESTADO_PROYECTO.md is not here)")
        return
    text = doc.read_text(encoding="utf-8")
    section = text.split("### 7.2")[1].split("## 8.")[0] if "### 7.2" in text else ""
    assert section, "ESTADO_PROYECTO.md has no section 7.2: the thresholds are not pre-registered"
    fake = {"protein_rmsd": md_analyze._series([0.0], [0.0], unit="angstrom"),
            "ligand_rmsd": md_analyze._series([0.0], [0.0], unit="angstrom"),
            "com_distance": {**md_analyze._series([0.0], [0.0], unit="angstrom"),
                             "pbc_corrected_frames": 0, "non_orthorhombic_frames": 0},
            "hbonds": {"pairs": [], "pairs_found": 0}}
    used = md_analyze.reading(fake)["thresholds_used"]
    assert used == {"protein_plateau_angstrom": 0.5, "ligand_stable_angstrom": 2.0,
                    "com_stable_range_angstrom": 2.0, "hbond_persistent_occupancy": 0.5}, used
    for number in ("0,5 Å", "2,0 Å", "≥ 0,5"):
        assert number in section, f"{number} is not in the pre-registration"
    assert "2026-11-11" in section
    assert "no se consultó ninguna" in section, "the section must say no literature was consulted"
    print("ok  the four thresholds in the code are the four pre-registered in section 7.2")


def _two_copy_universe():
    """The same system with a SECOND copy of the ligand in another segment, 15 A away.

    `8C7Y TXV U/X` and `5JSM 6NB J/R` in this project's own dossier are real cases of this.
    """
    import MDAnalysis as mda

    u = _universe()
    second = _universe()
    second.atoms.positions = second.atoms.positions + np.array([15.0, 0.0, 0.0])
    lig2 = second.select_atoms("resname LIG")
    merged = mda.Merge(u.atoms, lig2)
    merged.add_TopologyAttr("segid", ["PROA", "HETA", "HETB"])
    merged.residues[-1].segment = merged.segments[2]
    return merged


def test_two_copies_of_the_ligand_are_refused_not_merged():
    """Guard 1: two copies measured together are not a molecule.

    Their joint centre of mass falls between them, where there is no ligand, and the earlier code
    reported "8 heavy of 12 atoms in resname LIG" as though it were one.
    """
    u = _two_copy_universe()
    assert u.select_atoms("resname LIG").residues.n_residues == 2
    try:
        md_analyze.select_ligand(u, "LIG")
    except MdAnalyzeError as exc:
        assert "2 copies" in str(exc), str(exc)
        assert "centre of mass falls between them" in str(exc), str(exc)
        assert "ligand_segid" in str(exc), str(exc)
        assert len(exc.extra["copies"]) == 2, exc.extra
    else:
        raise AssertionError("two copies of the ligand must be refused")

    # and naming one segment measures exactly that copy
    one = md_analyze.select_ligand(u, "LIG", "HETB")
    assert one["atoms"].n_atoms == 4, one["selection"]
    assert "segment HETB" in one["selection"], one["selection"]
    assert "1 of 2 copies measured" in one["copies_note"], one["copies_note"]
    centre = one["atoms"].center_of_mass()
    both = u.select_atoms("resname LIG and mass 1.5 to 100000").center_of_mass()
    assert abs(centre[0] - both[0]) > 5, (centre, both)
    print("ok  two ligand copies are refused; ligand_segid measures one and says 1 of 2")


def test_an_unreadable_chunk_is_refused_instead_of_leaving_a_gap():
    """Guard 3: the earlier code noted the failure and then chained the broken file anyway.

    This test prints an "Exception ignored in: ... AttributeError: 'DCDReader' object has no
    attribute '_file'" on stderr. That is MDAnalysis's own `__del__` running on a reader whose
    constructor failed, which is the point of the test; it is noise from the garbage collector,
    not a failure, and the assertions below are what decides.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, name="step5_1.dcd", with_minimization=False)
        (folder / "step5_2.dcd").write_bytes(b"this is not a DCD" * 40)
        _run_topology(folder)
        r = _run(folder, None)
    assert r["ok"] is False, r
    assert "step5_2.dcd" in r["error"], r["error"]
    assert "gap in it" in r["error"], r["error"]
    assert r["unreadable"] == ["step5_2.dcd"], r.get("unreadable")
    print("ok  an unreadable chunk is refused by name: a missing chunk is a hole in the time axis")


def test_every_rmsf_value_is_labelled_from_its_own_atom():
    """Guard 4: the labels used to come from a parallel list, so one missing CA shifted them all."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None)
    rmsf = r["rmsf"]
    labels = [e["residue"] for e in rmsf["per_residue"]]
    assert labels == r["pocket"]["residues"], (labels, r["pocket"]["residues"])
    assert rmsf["residues_without_ca"] == [], rmsf["residues_without_ca"]
    assert "every pocket residue contributed exactly one CA" in rmsf["residues_without_ca_note"]
    assert "not from a parallel list" in rmsf["criterion"]
    assert rmsf["max"]["residue"] in labels
    assert rmsf["selection"] == "the CA of 2 of 2 pocket residues", rmsf["selection"]
    print("ok  each RMSF value is labelled from the CA it was computed from, and gaps are named")


# ---------------------------------------------------------------------------------------------
# Third review round: the defects the second round found, each with the control that catches it
# ---------------------------------------------------------------------------------------------

def test_a_ligand_that_really_leaves_is_not_folded_back_into_the_site():
    """THE REGRESSION. An earlier periodic-image rule anchored each frame on the POCKET.

    That rule cannot tell a real displacement of more than half a box from a wrapping jump, and it
    chose wrong: a ligand walking 0 -> 20 -> 40 -> 58 A out of a 60 A box came back as
    0, 20, 20, 2 A, and the offered reading said "The pose held within that figure" about a ligand
    58 A from its site. The rule now anchors on the PREVIOUS FRAME, where the two cases are
    distinguishable, and this test is the one that fails if anyone anchors on the pocket again.
    """
    walk = [0.0, 20.0, 40.0, 58.0]            # steps of 20, 20, 18 A; half the box is 30 A
    assert max(walk) > BOX[0] / 2, "the walk must pass half a box, or it tests nothing"
    assert max(np.diff(walk)) < BOX[0] / 2, "no single step may pass half a box"
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, shifts=walk, with_minimization=False)
        _run_topology(folder)
        r = _run(folder, None)
    values = r["ligand_rmsd"]["values"]
    # a rigid translation of the ligand with the protein still: the RMSD IS the displacement
    assert np.allclose(values, walk, atol=1e-3), (values, walk)
    assert r["com_distance"]["pbc_corrected_frames"] == 0, r["com_distance"]["pbc_note"]
    offered = " ".join(r["reading"]["offered"])
    assert "The pose held" not in offered, offered
    assert "beyond that figure" in offered, offered
    print(f"ok  a ligand walking {walk[-1]:.0f} A out of a {BOX[0]:.0f} A box reads as "
          f"{values[-1]:.0f} A, not folded back")


def test_a_wrapped_ligand_and_a_walking_one_are_told_apart():
    """The same two trajectories through the same rule, side by side.

    This is the control the pocket anchor could never pass: under it, a wrap and a 58 A walk both
    came out small. One must be ~0 A and the other ~58 A, measured by the same code.
    """
    with tempfile.TemporaryDirectory() as tmp:
        wrapped_dir, walking_dir = Path(tmp) / "w", Path(tmp) / "k"
        write_wrapped_fixture(wrapped_dir)
        _run_topology(wrapped_dir)
        wrapped = _run(wrapped_dir, None)
        write_fixture(walking_dir, shifts=[0.0, 20.0, 40.0, 58.0], with_minimization=False)
        _run_topology(walking_dir)
        walking = _run(walking_dir, None)
    assert wrapped["ligand_rmsd"]["max"] < 1.0, wrapped["ligand_rmsd"]["values"]
    assert walking["ligand_rmsd"]["max"] > 50.0, walking["ligand_rmsd"]["values"]
    assert wrapped["com_distance"]["pbc_corrected_frames"] > 0
    assert walking["com_distance"]["pbc_corrected_frames"] == 0
    print(f"ok  a wrap reads {wrapped['ligand_rmsd']['max']:.2f} A and a real walk "
          f"{walking['ligand_rmsd']['max']:.1f} A, through the same rule")


def _two_chain_universe():
    """The fixture's protein in TWO segments, PROA and PROB, 20 A apart. PROB has no ligand."""
    import MDAnalysis as mda

    u = _universe()
    other = _universe()
    other.atoms.positions = other.atoms.positions + np.array([0.0, 0.0, 20.0])
    merged = mda.Merge(u.atoms, other.select_atoms("protein"))
    merged.add_TopologyAttr("segid", ["PROA", "HETA", "PROB"])
    for residue in merged.residues[4:]:
        residue.segment = merged.segments[2]
    return merged


def test_two_protein_chains_are_never_superposed_together():
    """With `wrapAll on` two chains are two molecules, and one can be re-imaged alone.

    A single joint backbone RMSD then reads a box length while nothing has moved. Each segment is
    superposed on ITSELF, which removes any whole-molecule translation, so the right answer for
    this trajectory is 0 A per segment. The joint fit is computed here and asserted to be WRONG,
    so it cannot quietly come back.
    """
    import MDAnalysis as mda
    from MDAnalysis.analysis import align

    u = _two_chain_universe()
    assert sorted(set(u.select_atoms("protein").segids)) == ["PROA", "PROB"], \
        sorted(set(u.select_atoms("protein").segids))
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        folder.mkdir(parents=True, exist_ok=True)
        prob = u.select_atoms("protein and segid PROB")
        base = prob.positions.copy()
        path = folder / "step5_1.dcd"
        # frame 1 re-images PROB alone by a whole box in z: NAMD does exactly this, per molecule
        with mda.Writer(str(path), n_atoms=len(u.atoms), dt=2.0) as writer:
            for jump in (0.0, BOX[2]):
                prob.positions = base + np.array([0.0, 0.0, jump], dtype=np.float32)
                u.trajectory.ts.dimensions = np.array(BOX, dtype=np.float32)
                writer.write(u.atoms)
        (folder / "step5_production.inp").write_text(
            NAMD_INPUT.replace("minimize                2000\n", ""))
        traj = mda.Universe(u._topology, str(path))
        try:
            ligand = md_analyze.select_ligand(traj, "LIG")
            pocket = md_analyze.select_pocket(traj, ligand["atoms"], md_analyze.POCKET_CUTOFF, 0)
            measured = md_analyze.measure(traj, ligand, pocket, first_frame=0)
            # what a JOINT fit would have said about the same two frames
            both = traj.select_atoms("protein and name N CA C O")
            traj.trajectory[0]
            ref = both.positions - both.center_of_mass()
            masses = both.masses
            traj.trajectory[1]
            _, joint = align.rotation_matrix(both.positions - both.center_of_mass(), ref,
                                             weights=masses)
        finally:
            traj.trajectory.close()

    per_segment = measured["protein_rmsd"]["per_segment"]
    assert set(per_segment) == {"PROA", "PROB"}, set(per_segment)
    for segid, series in per_segment.items():
        assert max(series["values"]) < 1e-3, (segid, series["values"])
    assert joint > 5.0, f"the joint fit must be visibly wrong, got {joint}"
    assert "two chains are two molecules" in measured["protein_rmsd"]["per_segment_note"]
    print(f"ok  a chain re-imaged alone reads 0.00 A per segment; a joint fit would say "
          f"{joint:.1f} A")


def test_the_hydrogen_bonds_follow_the_named_ligand_copy():
    """`ligand_segid` chose the copy for the pocket and the RMSD and was ignored here.

    The result then said "1 of 2 copies measured" above a list of bonds that included the other
    copy's, which is a sourced-looking claim about a molecule nothing in the report measured.
    """
    u = _two_copy_universe()
    both = md_analyze.hydrogen_bonds(u, "LIG", first_frame=0, frames_analysed=1)
    one = md_analyze.hydrogen_bonds(u, "LIG", first_frame=0, frames_analysed=1, segid="HETA")
    assert one["ligand_selection"] == "resname LIG and segid HETA", one["ligand_selection"]
    assert both["eligible"]["acceptors_on_the_ligand"] == \
        2 * one["eligible"]["acceptors_on_the_ligand"], (both["eligible"], one["eligible"])
    assert one["eligible"]["hydrogens_on_the_ligand"] == 1, one["eligible"]
    for bond in one["pairs"]:
        assert "HETB" not in bond["donor"] and "HETB" not in bond["acceptor"], bond
    print(f"ok  ligand_segid reaches the hydrogen bonds: {one['eligible']['acceptors_on_the_ligand']}"
          f" of {both['eligible']['acceptors_on_the_ligand']} ligand acceptors were eligible")


def test_hydrogen_bond_labels_survive_ids_that_are_not_indices():
    """`count_by_ids()` returns ATOM IDS; the code indexed `universe.atoms` with them.

    It agreed with the truth only because a PSF numbers its atoms from 1 and MDAnalysis's PSF
    parser stores `id - 1`. Shift the ids and every donor, acceptor and residue name in the report
    moves by one atom; on a topology with no ids at all it raises NoDataError. That is a
    hallucinated residue name with a measurement's face on it.

    `count_by_ids()` is gone: `results.hbonds` already stores INDICES, so the frames are counted
    from those. This test holds the labels still while the ids move underneath them.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        psf = _run_topology(folder)
        import MDAnalysis as mda
        u = mda.Universe(str(psf), str(folder / "step5_1.dcd"))
        try:
            straight = md_analyze.hydrogen_bonds(u, "LIG", first_frame=0, frames_analysed=4)
            # the same system, renumbered: ids no longer coincide with positions
            u.atoms.ids = u.atoms.ids + 1000
            shifted = md_analyze.hydrogen_bonds(u, "LIG", first_frame=0, frames_analysed=4)
        finally:
            u.trajectory.close()
    assert straight["pairs_found"] > 0, straight
    naming = [(b["donor"], b["hydrogen"], b["acceptor"]) for b in straight["pairs"]]
    assert [(b["donor"], b["hydrogen"], b["acceptor"]) for b in shifted["pairs"]] == naming, \
        (naming, shifted["pairs"])
    # and a topology with NO ids at all must still work, which count_by_ids() cannot do
    no_ids = md_analyze.hydrogen_bonds(_universe(), "LIG", first_frame=0, frames_analysed=1)
    assert no_ids["pairs_found"] > 0, no_ids
    print(f"ok  the {straight['pairs_found']} hydrogen-bond label(s) survive renumbering the "
          "atoms, and a topology with no ids at all is measured instead of raising")


def test_a_hydrogen_on_a_carbon_is_never_a_donor():
    """The `bonded` clause is what keeps C-H out, and the fixture now carries one to prove it."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        _run_topology(folder)
        r = _run(folder, None)
    eligible = r["hbonds"]["eligible"]
    u = _universe()
    all_h = u.select_atoms(f"(resname LIG or protein) and {md_analyze.HB_H_BY_MASS}")
    polar_h = u.select_atoms(f"(resname LIG or protein) and {md_analyze.HB_H_BY_MASS} "
                             f"and bonded {md_analyze.HB_NO_BY_MASS}")
    # 3 backbone HN + the ligand's H1 and H2 = 5 hydrogens; only 4 hang off an N or an O
    assert len(all_h) == 5, len(all_h)
    assert len(polar_h) == 4, len(polar_h)
    assert len(all_h) != len(polar_h), "the fixture must carry a hydrogen on carbon"
    assert eligible["hydrogens"] == len(polar_h), (eligible, len(polar_h))
    assert "H2" not in set(polar_h.names), set(polar_h.names)
    print(f"ok  {len(all_h) - len(polar_h)} hydrogen on carbon is excluded: "
          f"{eligible['hydrogens']} eligible of {len(all_h)}")


def test_a_pocket_residue_without_a_ca_is_named_and_not_mislabelled():
    """A residue with no CA must lose its RMSF value and be listed, not shift every label after it."""
    import MDAnalysis as mda

    u = _universe()
    # drop the CA of pocket residue 2, which is what a truncated or non-standard residue looks like
    keep = u.select_atoms("not (resid 2 and segid PROA and name CA)")
    trimmed = mda.Merge(keep)
    trimmed.add_TopologyAttr("segid", ["PROA", "HETA"])
    trimmed.residues[-1].segment = trimmed.segments[1]
    trimmed.load_new(keep.positions[None, :, :], order="fac")
    ligand = md_analyze.select_ligand(trimmed, "LIG")
    pocket = md_analyze.select_pocket(trimmed, ligand["atoms"], md_analyze.POCKET_CUTOFF, 0)
    measured = md_analyze.measure(trimmed, ligand, pocket, first_frame=0)
    rmsf = measured["rmsf"]
    labels = [e["residue"] for e in rmsf["per_residue"]]
    assert len(pocket["residues"]) == 2, pocket["residues"]
    assert labels == ["PROA:ALA1"], labels
    assert rmsf["residues_without_ca"] == ["PROA:ALA2"], rmsf["residues_without_ca"]
    assert "carry NO rmsf value" in rmsf["residues_without_ca_note"]
    assert rmsf["selection"] == "the CA of 1 of 2 pocket residues", rmsf["selection"]
    print("ok  a pocket residue without a CA is named and dropped, and no label shifts onto it")


def test_a_triclinic_box_is_not_treated_as_rectangular():
    """The re-imaging rule assumes 90 deg angles. A box that has none must be counted and declared."""
    triclinic = [60.0, 60.0, 60.0, 90.0, 60.0, 90.0]
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, box=triclinic)
        _run_topology(folder)
        r = _run(folder, None)
    com = r["com_distance"]
    assert com["non_orthorhombic_frames"] == 4, com
    assert com["pbc_corrected_frames"] == 0, com
    assert "not rectangular" in com["non_orthorhombic_note"]
    assert any("not\nrectangular" in n or "not rectangular" in n for n in r["reading"]["offered"]), \
        r["reading"]["offered"]
    print(f"ok  {com['non_orthorhombic_frames']} triclinic frames are counted, left uncorrected "
          "and said so in the reading")


def test_a_frame_without_a_box_is_not_called_rectangular():
    """A DCD written with no unit cell must be reported as unknown, not as a 90 deg box."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, box=[0.0] * 6)
        _run_topology(folder)
        r = _run(folder, None)
    com = r["com_distance"]
    assert com["frames_without_box"] == 4, com
    assert com["non_orthorhombic_frames"] == 0, com
    assert com["pbc_corrected_frames"] == 0, com
    assert "declare NO box at all" in com["frames_without_box_note"]
    assert any("no box at all" in n for n in r["reading"]["offered"]), r["reading"]["offered"]
    print("ok  4 frames without a unit cell are reported as unknown, not as rectangular")


def test_the_span_is_the_elapsed_time_between_the_first_and_the_last_frame():
    """4 frames 2 ps apart span 6 ps, not 8: `frames x dt` counts an interval that is not there."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, dt=2.0)
        _run_topology(folder)
        r = _run(folder, None)
    traj = r["trajectory"]
    assert traj["frames_analysed"] == 4, traj["frames_analysed"]
    assert traj["dt_ps"] == 2.0, traj["dt_ps"]
    assert traj["span_ps"] == 6.0, traj["span_ps"]
    assert traj["span_ps"] == r["ligand_rmsd"]["times_ps"][-1], traj["span_ps"]
    assert traj["span_ns"] == 0.006, traj["span_ns"]
    print(f"ok  4 frames 2 ps apart span {traj['span_ps']} ps, not {4 * 2.0}")


def test_the_hydrogen_bonds_start_at_the_same_frame_as_everything_else():
    """If HBA ignored `start`, the minimization frames would feed the occupancy denominator."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)                      # 2 minimization frames + 4 dynamics
        _run_topology(folder)
        dropped = _run(folder, None)
        kept = _run(folder, None, skip_frames=0, reference_frame=2)
    assert dropped["trajectory"]["frames_analysed"] == 4
    assert kept["trajectory"]["frames_analysed"] == 6
    # the minimization frames park the ligand 15 A away, where no hydrogen bond can hold, so
    # counting them can only ADD frames to the denominator and never to the numerator
    pairs_dropped = {(b["donor"], b["acceptor"]): b for b in dropped["hbonds"]["pairs"]}
    pairs_kept = {(b["donor"], b["acceptor"]): b for b in kept["hbonds"]["pairs"]}
    assert pairs_dropped, dropped["hbonds"]
    for key, bond in pairs_dropped.items():
        assert bond["frames_analysed"] == 4, bond
        assert pairs_kept[key]["frames_present"] == bond["frames_present"], (bond, pairs_kept[key])
        assert pairs_kept[key]["occupancy_sampled"] < bond["occupancy_sampled"], (
            bond, pairs_kept[key])
    print("ok  the hydrogen bonds start at the first analysed frame: the two skipped frames "
          "change the denominator and not the count")


RESTRAINED_INPUT = """\
dcdfreq                 1000;
constraints             on
consexp                 2
consref                 restraints/prot_posres.ref
conskfile               restraints/prot_posres.ref
conskcol                B
constraintScaling       1.0
run                     4000
"""


def test_positional_restraints_in_the_input_are_declared():
    """H7. `step4_equilibration.inp` holds the backbone with springs; the plateau cannot fail.

    Reporting a 0.04 A drift from such a trajectory as evidence of equilibration is the kind of
    claim this project exists to prevent, so the restraint is read from the input and the reading
    says in writing that the line carries no evidence.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, inp=False)
        (folder / "step5_production.inp").write_text(RESTRAINED_INPUT)
        _run_topology(folder)
        held = _run(folder, None)
        (folder / "step5_production.inp").write_text(
            NAMD_INPUT.replace("minimize                2000\n", ""))
        free = _run(folder, None)
    assert held["trajectory"]["restraints"]["positional_restraints"] is True, \
        held["trajectory"]["restraints"]
    evidence = " ".join(held["trajectory"]["restraints"]["evidence"])
    assert "constraints             on" in evidence, evidence
    assert "prot_posres.ref" in evidence, evidence
    assert "CARRIES NO EVIDENCE" in held["reading"]["offered"][0], held["reading"]["offered"][0]

    assert free["trajectory"]["restraints"]["positional_restraints"] is False, \
        free["trajectory"]["restraints"]
    assert "CARRIES NO EVIDENCE" not in free["reading"]["offered"][0], free["reading"]["offered"][0]
    print("ok  `constraints on` is read from the input and the plateau line is marked as "
          "carrying no evidence")


def test_no_input_for_this_stage_is_unknown_and_not_absent():
    """H4. "No input declares a minimization" and "no input was read" are different facts."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, inp=False)
        _run_topology(folder)
        silent = _run(folder, None)
        (folder / "step5_production.inp").write_text(
            NAMD_INPUT.replace("minimize                2000\n", ""))
        read = _run(folder, None)
    assert silent["trajectory"]["minimization"]["inputs_read"] == []
    assert "UNKNOWN, not absent" in silent["trajectory"]["minimization_note"], \
        silent["trajectory"]["minimization_note"]
    assert "unknown rather than absent" in silent["trajectory"]["restraints"]["note"], \
        silent["trajectory"]["restraints"]["note"]
    assert read["trajectory"]["minimization"]["inputs_read"] == ["step5_production.inp"]
    assert "none declares a minimization" in read["trajectory"]["minimization_note"], \
        read["trajectory"]["minimization_note"]
    print("ok  a missing NAMD input is reported as unknown, not as an input that declared nothing")


def test_the_same_trajectory_reached_twice_is_read_once():
    """H6. `folder` and `folder/system/namd` can be the same place, and a chunk chained twice
    doubles the run: every mean, plateau window and occupancy would be computed over it."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        _run_topology(folder)
        # a directory junction, which `mklink /J` makes on Windows without any privilege and
        # which `Path.resolve()` follows: `folder/namd` IS `folder`, so every file is reachable
        # twice. On a POSIX box the same thing is a directory symlink.
        link = folder / "namd"
        try:
            link.symlink_to(folder, target_is_directory=True)
        except (OSError, NotImplementedError):
            import subprocess
            done = subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(folder)],
                                  capture_output=True)
            if done.returncode != 0 or not link.exists():
                print("ok  (skipped: this machine makes neither a symlink nor a junction)")
                return
        found = md_analyze.find_inputs(folder)
        r = _run(folder, None)
    assert [d.name for d in found["dcds"]] == ["step5_1.dcd"], [str(d) for d in found["dcds"]]
    assert r["trajectory"]["frames_analysed"] == 4, r["trajectory"]["frames_note"]
    print("ok  the same DCD reached through two folders is read once, not chained onto itself")


def test_the_reference_frame_moves_the_rmsd_zero_as_well_as_the_pocket():
    """H9. `reference_frame` used to move only the pocket selection, which no docstring said.

    The RMSD then measured against the first ANALYSED frame while the pocket came from another
    one, so the reported `reference` named a frame the numbers were not taken against.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        _run_topology(folder)
        first = _run(folder, None)
        third = _run(folder, None, reference_frame=2)
    # the fixture translates the ligand by 0, 1, 2, 3 A: against frame 0 that is 0,1,2,3
    assert np.allclose(first["ligand_rmsd"]["values"], SHIFTS_DYNAMICS, atol=1e-3), \
        first["ligand_rmsd"]["values"]
    # against frame 2 it is |shift - 2| = 2, 1, 0, 1
    expected = [abs(s - SHIFTS_DYNAMICS[2]) for s in SHIFTS_DYNAMICS]
    assert np.allclose(third["ligand_rmsd"]["values"], expected, atol=1e-3), \
        (third["ligand_rmsd"]["values"], expected)
    assert "first analysed frame" in first["ligand_rmsd"]["reference"]
    assert "you named with `reference_frame`" in third["ligand_rmsd"]["reference"], \
        third["ligand_rmsd"]["reference"]
    print(f"ok  reference_frame=2 moves the RMSD zero: {third['ligand_rmsd']['values']} "
          f"instead of {first['ligand_rmsd']['values']}")


# ---------------------------------------------------------------------------------------------
# Fourth review round: the restraint file, and the 11 mutations the third round let through
# ---------------------------------------------------------------------------------------------

def write_ref(folder: Path, *, hold_ligand=True, column="B", keyword="constraints",
              declare_column=True, shift=0, extra_atoms=0, semicolons=False,
              consref_only=False, atoms=True, name="prot_posres.ref") -> Path:
    """A NAMD restraint reference PDB, in the layout CHARMM-GUI's `prot_posres.ref` uses.

    The force constant lives in ONE column of a PDB line and the input file says which. The real
    file of job 20070766 marks 4389 of 139989 atoms, and among them **the ligand's 33 heavy
    atoms**, which is the fact that made this whole function necessary: the tool used to name the
    file and mark only the protein line.

    INTERLEAVED `TER` LINES ARE WRITTEN, because the real file has five of them and four are in
    the middle - one immediately before the ligand's first atom (line 8718, with `C1 LIG` on
    8719). A fixture without them is cleaner than the real file on the one axis that decides
    whether every index after each TER is off by one, and the lesson of this project is that a
    fixture must not be cleaner than the file it stands for.

    The keyword arguments are the ways the file can be wrong, each one a case the tool has to
    survive: `shift` reorders the atoms (a posres from a regenerated build), `extra_atoms` makes it
    describe a bigger system, `declare_column=False` leaves the column unsaid, `semicolons` adds
    CHARMM-GUI's trailing `;`, `consref_only` omits `conskfile`, and `atoms=False` writes a file
    with no ATOM line at all.
    """
    u = _universe()
    folder.mkdir(parents=True, exist_ok=True)
    ref = folder / "restraints"
    ref.mkdir(exist_ok=True)
    path = ref / name
    order = list(u.atoms)
    if shift:
        order = order[shift:] + order[:shift]
    lines = []
    for position, atom in enumerate(order):
        ligand = atom.resname == "LIG"
        heavy = atom.mass > md_analyze.HEAVY_MASS_MIN
        held = 1.0 if (heavy and (not ligand or hold_ligand)) else 0.0
        x, y, z = atom.position
        # the columns NAMD reads: 31-38 X, 39-46 Y, 47-54 Z, 55-60 occupancy, 61-66 B
        values = {"X": x, "Y": y, "Z": z, "O": 1.0, "B": 0.0}
        values[column.upper()] = held
        # a TER before the ligand, exactly where the real file has one
        if ligand and position and order[position - 1].resname != "LIG":
            lines.append("TER")
        lines.append(
            f"ATOM  {position + 1:5d} {atom.name:<4s} {atom.resname:<3s} A{atom.resid:4d}    "
            f"{values['X']:8.3f}{values['Y']:8.3f}{values['Z']:8.3f}"
            f"{values['O']:6.2f}{values['B']:6.2f}")
    for extra in range(extra_atoms):
        lines.append(f"ATOM  {len(order) + extra + 1:5d} OH2  TIP3 B{extra + 1:4d}    "
                     f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}")
    lines.insert(0, "REMARK  written by test_md_analyze")
    path.write_text(("\n".join(lines) if atoms else "REMARK nothing here") + "\nTER\nEND\n")

    col_key = "conskcol" if keyword == "constraints" else "fixedAtomsCol"
    file_key = "conskfile" if keyword == "constraints" else "fixedAtomsFile"
    end = ";" if semicolons else ""
    text = f"dcdfreq                 1000;\n{keyword}             on\n"
    text += f"consref                 restraints/{name}{end}\n"
    if not consref_only:
        text += f"{file_key}               restraints/{name}{end}\n"
    if declare_column:
        text += f"{col_key}                {column}{end}\n"
    text += "constraintScaling       1.0\nrun                     4000\n"
    (folder / "step5_production.inp").write_text(text)
    return path


def test_the_restraint_file_is_read_and_says_which_atoms_were_held():
    """THE SECOND REGRESSION, and the same shape as the first: a claim the files contradict.

    `step4_equilibration.inp` of job 20070766 names `restraints/prot_posres.ref`, and that file
    holds down **33 of the ligand's 52 atoms** - exactly the 33 heavy atoms this module measures,
    at k = 1 kcal/mol/A^2. The earlier version marked only the protein line as carrying no
    evidence, so the two ligand lines were offered without reserve about a ligand on a spring, and
    ESTADO_PROYECTO.md said in writing that the restraints "do not hold the ligand".

    So the reference file is read, atom by atom, and the three lines are marked from what it says.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, inp=False)
        write_ref(folder, hold_ligand=True)
        _run_topology(folder)
        both = _run(folder, None)
        write_ref(folder, hold_ligand=False)
        protein_only = _run(folder, None)

    held = both["trajectory"]["restraints"]
    assert held["atoms_held"] is True, held
    assert held["reference_file"] == "restraints/prot_posres.ref", held
    assert held["force_constant_column"] == "B", held
    lig_held = held["held"]["ligand_heavy"]
    assert lig_held["known"] is True and lig_held["held"] == lig_held["of"] == 4, lig_held
    assert held["held"]["protein_backbone"]["held"] > 0, held["held"]
    marked = [n for n in both["reading"]["offered"] if "CARRIES NO EVIDENCE" in n]
    assert len(marked) == 3, [n[:80] for n in both["reading"]["offered"]]

    # and when the file holds the protein only, the LIGAND lines must come back unmarked: a blanket
    # warning on every line would be as useless as none at all
    lig_free = protein_only["trajectory"]["restraints"]["held"]["ligand_heavy"]
    assert lig_free["held"] == 0, lig_free
    marked_free = [n for n in protein_only["reading"]["offered"] if "CARRIES NO EVIDENCE" in n]
    assert len(marked_free) == 1, [n[:80] for n in protein_only["reading"]["offered"]]
    assert "Protein backbone" in marked_free[0], marked_free[0]
    print(f"ok  the restraint file is read: {lig_held['held']} of {lig_held['of']} ligand heavy "
          "atoms held marks 3 lines, and holding the protein alone marks 1")


def test_fixed_atoms_are_looked_for_and_not_merely_claimed():
    """The negative note used to say `fixedAtoms` had been looked for while only `constraints` was.

    `fixedAtoms on` freezes atoms outright, so a plateau over such a trajectory is even emptier
    than one over a restrained trajectory, and the tool was saying it had checked.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, inp=False)
        write_ref(folder, keyword="fixedAtoms")
        _run_topology(folder)
        frozen = _run(folder, None)
    held = frozen["trajectory"]["restraints"]
    assert held["atoms_held"] is True, held
    assert "frozen" in held["kind"], held["kind"]
    assert held["held"]["ligand_heavy"]["held"] == 4, held["held"]
    assert any("CARRIES NO EVIDENCE" in n for n in frozen["reading"]["offered"])
    print(f"ok  `fixedAtoms on` is detected as {held['kind']}, not reported as nothing held")


def test_namd_keywords_are_read_case_insensitively():
    """NAMD does not care about case, and `DCDfreq` is legal. The code did care.

    With `DCDfreq`, `minimization_frames` found nothing and the report then said "an input was
    read and declares no minimization" while two minimization frames entered the time axis as
    20 ps of dynamics.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, inp=False)
        _run_topology(folder)
        results = {}
        for label, text in (
                ("lower", NAMD_INPUT),
                ("DCDfreq", NAMD_INPUT.replace("dcdfreq", "DCDfreq")),
                ("Minimize", NAMD_INPUT.replace("minimize", "Minimize")),
                ("MINIMIZE", NAMD_INPUT.replace("minimize", "MINIMIZE").replace("dcdfreq",
                                                                                "DcdFreq"))):
            (folder / "step5_production.inp").write_text(text)
            results[label] = _run(folder, None)
    for label, r in results.items():
        assert r["trajectory"]["minimization"]["frames"] == 2, (label, r["trajectory"]["minimization"])
        assert r["trajectory"]["frames_analysed"] == 4, (label, r["trajectory"]["frames_note"])
        assert "They were dropped." in r["trajectory"]["minimization_note"], label
    print(f"ok  {len(results)} spellings of `minimize`/`dcdfreq` all read as 2 frames")


def test_too_few_frames_cannot_support_a_plateau_and_says_so():
    """A quarter window of ONE frame makes "the drift" a difference between two single frames.

    The earlier guard tested whether the two quarters OVERLAP, which with `window = round(n/4)`
    can only happen at n = 1. So n = 2 - where the drift is literally frame 1 minus frame 0 - came
    back as `drift_is_meaningful: True` and the reading offered "That is consistent with a
    plateau". The condition is now that each quarter hold more than one frame.
    """
    assert md_analyze.PLATEAU_MIN_FRAMES == 8, md_analyze.PLATEAU_MIN_FRAMES   # a literal
    # the arithmetic the old guard relied on: the windows overlap at n = 1 and NOWHERE else, which
    # is why it could not catch n = 2. Computed here so the claim is not taken on trust.
    overlaps = {n: max(0, 2 * max(1, round(n * md_analyze.PLATEAU_WINDOW)) - n)
                for n in (1, 2, 3, 4, 8, 25, 27)}
    assert [n for n, o in overlaps.items() if o] == [1], overlaps

    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, shifts=[float(i) * 0.1 for i in range(8)], with_minimization=False)
        _run_topology(folder)
        runs = {n: _run(folder, None, skip_frames=8 - n) for n in (1, 2, 4, 8)}
    for n in (1, 2, 4):
        prot = runs[n]["protein_rmsd"]
        assert runs[n]["trajectory"]["frames_analysed"] == n, runs[n]["trajectory"]["frames_note"]
        assert prot["drift_is_meaningful"] is False, (n, prot["drift"], prot["drift_warning"])
        assert "CANNOT SUPPORT A STATEMENT ABOUT A PLATEAU" in prot["drift_warning"], (n, prot)
        assert str(md_analyze.PLATEAU_MIN_FRAMES) in prot["drift_warning"], (n, prot)
        assert any("CANNOT SUPPORT A STATEMENT" in line for line in runs[n]["reading"]["offered"]), \
            (n, runs[n]["reading"]["offered"][0])
    # and at the chosen minimum the warning must be gone, or it would fire for every run
    enough = runs[8]["protein_rmsd"]
    assert runs[8]["trajectory"]["frames_analysed"] == 8
    assert enough["drift_is_meaningful"] is True, enough
    assert enough["drift_warning"] is None, enough["drift_warning"]
    assert not any("CANNOT SUPPORT" in line for line in runs[8]["reading"]["offered"])
    print(f"ok  1, 2 and 4 analysed frames cannot support a plateau and say so; "
          f"{md_analyze.PLATEAU_MIN_FRAMES} can. The old overlap test only fired at n=1")


def write_flexing_fixture(folder: Path) -> Path:
    """One pocket CA oscillates; everything else is still. Gives RMSF an analytic shape.

    This is what the earlier fixture could not do: with RMSF 0 everywhere, `argmax` and `argmin`
    name the same residue and a deviation about frame 0 equals one about the mean.
    """
    import MDAnalysis as mda

    folder.mkdir(parents=True, exist_ok=True)
    u = _universe()
    ca2 = u.select_atoms("protein and resid 2 and name CA")
    base = ca2.positions.copy()
    path = folder / "step5_1.dcd"
    with mda.Writer(str(path), n_atoms=len(u.atoms), dt=2.0) as writer:
        for wobble in (-1.0, 1.0, -1.0, 1.0):
            ca2.positions = base + np.array([0.0, 0.0, wobble], dtype=np.float32)
            u.trajectory.ts.dimensions = np.array(BOX, dtype=np.float32)
            writer.write(u.atoms)
    (folder / "step5_production.inp").write_text(
        NAMD_INPUT.replace("minimize                2000\n", ""))
    return path


def test_the_rmsf_names_the_most_flexible_residue_and_is_taken_about_the_mean():
    """Two mutations at once: `argmax` -> `argmin`, and the mean -> the first frame.

    Both survived while the fixture held every CA still, which made every RMSF 0 and the two
    answers identical. One residue now wobbles, so the right answer is one residue and the wrong
    answers are visibly different numbers, computed here.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_flexing_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None)
    rows = {e["residue"]: e["rmsf_angstrom"] for e in r["rmsf"]["per_residue"]}
    assert set(rows) == {"PROA:ALA1", "PROA:ALA2"}, rows
    assert rows["PROA:ALA2"] > 3 * max(rows["PROA:ALA1"], 1e-6), rows
    assert r["rmsf"]["max"]["residue"] == "PROA:ALA2", (r["rmsf"]["max"], rows)
    # the wobble is +-1 A about its own mean, so the answer about the MEAN is ~1.0 A; about the
    # FIRST frame the deviations are 0, 2, 0, 2, which is sqrt(2) ~ 1.414 A. The code must match
    # the first and not the second.
    assert abs(rows["PROA:ALA2"] - 1.0) < 0.2, rows
    assert abs(rows["PROA:ALA2"] - np.sqrt(2.0)) > 0.2, rows
    print(f"ok  the wobbling residue is named as the most flexible at "
          f"{rows['PROA:ALA2']:.3f} A (about its mean, not about frame 0)")


def write_marching_fixture(folder: Path) -> Path:
    """The WHOLE system translates rigidly, protein and ligand together.

    Nothing moves relative to anything, so every metric must read 0 and the centre-of-mass
    distance must stay put. A pocket centre of mass frozen at the reference frame instead of
    recomputed per frame turns this into a ligand marching away from its site.
    """
    import MDAnalysis as mda

    folder.mkdir(parents=True, exist_ok=True)
    u = _universe()
    base = u.atoms.positions.copy()
    path = folder / "step5_1.dcd"
    with mda.Writer(str(path), n_atoms=len(u.atoms), dt=2.0) as writer:
        for march in (0.0, 4.0, 8.0, 12.0):
            u.atoms.positions = base + np.array([0.0, march, 0.0], dtype=np.float32)
            u.trajectory.ts.dimensions = np.array(BOX, dtype=np.float32)
            writer.write(u.atoms)
    (folder / "step5_production.inp").write_text(
        NAMD_INPUT.replace("minimize                2000\n", ""))
    return path


def test_a_rigid_translation_of_everything_moves_no_metric():
    """The pocket centre of mass is taken per frame, not once at the reference frame."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_marching_fixture(folder)
        _run_topology(folder)
        r = _run(folder, None)
    assert max(r["ligand_rmsd"]["values"]) < 1e-3, r["ligand_rmsd"]["values"]
    assert r["com_distance"]["range"] < 1e-3, r["com_distance"]["values"]
    assert max(r["protein_rmsd"]["values"]) < 1e-3, r["protein_rmsd"]["values"]
    # and the control: a pocket centre frozen at frame 0 would have read the march itself
    assert max(r["com_distance"]["values"]) < 12.0, r["com_distance"]["values"]
    print(f"ok  the whole system marching 12 A leaves the COM distance in a "
          f"{r['com_distance']['range']} A range")


def test_the_angle_cutoff_decides_whether_a_bond_is_counted():
    """`d_h_a_angle_cutoff` was never tied to behaviour: setting it to 0 kept every test green.

    The criterion string would then declare ">= 150 deg" about a run that tested no angle at all,
    which is a declared criterion that is simply false. The geometry here is within the distance
    cutoff and BENT, so the angle is the only thing that can decide it.
    """
    import MDAnalysis as mda

    u = _universe()
    # bend the N1-H1 bond away from O of residue 1 while keeping N1 within 3.5 A of it
    acceptor = u.select_atoms("protein and resid 1 and name O")[0]
    donor = u.select_atoms("resname LIG and name N1")[0]
    hydrogen = u.select_atoms("resname LIG and name H1")[0]
    donor.position = acceptor.position + np.array([3.0, 0.0, 0.0])
    straight = acceptor.position - donor.position
    straight = straight / np.linalg.norm(straight)
    hydrogen.position = donor.position + straight          # pointing AT the acceptor: ~180 deg
    aligned = md_analyze.hydrogen_bonds(u, "LIG", first_frame=0, frames_analysed=1)
    hydrogen.position = donor.position - straight          # pointing AWAY: ~0 deg
    bent = md_analyze.hydrogen_bonds(u, "LIG", first_frame=0, frames_analysed=1)

    pair = ("HETA:LIG1:N1", "PROA:ALA1:O")
    found = {(b["donor"], b["acceptor"]) for b in aligned["pairs"]}
    assert pair in found, found
    assert pair not in {(b["donor"], b["acceptor"]) for b in bent["pairs"]}, bent["pairs"]
    # distance is identical in both, so only the angle can have decided it
    assert float(np.linalg.norm(acceptor.position - donor.position)) < md_analyze.HB_DA_CUTOFF
    print("ok  the same donor-acceptor distance is counted when aligned and not when bent: "
          "the angle criterion is doing work")


def test_the_ligand_side_of_each_bond_is_named_correctly():
    """`ligand_side` was asserted nowhere, so inverting it was invisible."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        _run_topology(folder)
        r = _run(folder, None)
    pairs = r["hbonds"]["pairs"]
    assert pairs, r["hbonds"]
    for bond in pairs:
        on_ligand = "LIG" in bond["donor"]
        assert bond["ligand_side"] == ("donor" if on_ligand else "acceptor"), bond
        assert ("LIG" in bond["acceptor"]) != on_ligand, bond
    print(f"ok  the ligand side of all {len(pairs)} bond(s) matches which end carries the ligand")


def test_the_per_frame_bond_count_is_aligned_with_the_analysed_frames():
    """The series was indexed by absolute frame number, so skipping frames shifted it.

    Nothing asserted the series, so leaving out `- first_frame` left zeros at the start and
    dropped the last frames off the end, silently.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)                     # 2 minimization frames, ligand parked 15 A away
        _run_topology(folder)
        r = _run(folder, None, include_series=True)
    series = r["hbonds"]["count_per_frame"]
    total = sum(b["frames_present"] for b in r["hbonds"]["pairs"])
    assert len(series["values"]) == 4, series["values"]
    assert sum(series["values"]) == total, (series["values"], total)
    # the bond exists at the first analysed frame, so the series cannot start at 0
    assert series["values"][0] > 0, series["values"]
    # The KEY must be gone, not merely annotated: a consumer reads the number, not the note, and
    # `times_ps` beside a frame index turned frame 24 of a 240 ps window into 24 ps.
    assert "times_ps" not in series, series.keys()
    assert "times_note" not in series, series.keys()
    assert series["frame_indices"] == [0, 1, 2, 3], series["frame_indices"]
    assert "NO time axis" in series["note"], series["note"]
    print(f"ok  the per-frame bond count sums to the {total} frame-presences and starts where the "
          "analysis starts")


def test_the_headline_segment_is_the_largest_drift_and_says_so():
    """A chain swinging widely about one position has a drift near zero.

    The reading used to call the headline segment "the one that moved most", which the code never
    computed. Here PROB has 25x the range of PROA and PROA has all the drift, so the two readings
    name different chains and the wording has to match the rule.
    """
    import MDAnalysis as mda

    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        folder.mkdir(parents=True, exist_ok=True)
        u = _two_chain_universe()
        prob = u.select_atoms("protein and segid PROB and name CA")
        proa = u.select_atoms("protein and segid PROA and resid 1 and name CA")
        base_b, base_a = prob.positions.copy(), proa.positions.copy()
        path = folder / "step5_1.dcd"
        with mda.Writer(str(path), n_atoms=len(u.atoms), dt=2.0) as writer:
            # PROB swings out and comes BACK, so its first and last quarter coincide and its drift
            # is 0 while its range is large. PROA ramps steadily, so all the drift is its own.
            for swing, ramp in ((0.0, 0.0), (8.0, 0.1), (8.0, 0.2), (0.0, 0.3)):
                prob.positions = base_b + np.array([0.0, 0.0, swing], dtype=np.float32)
                proa.positions = base_a + np.array([0.0, 0.0, ramp], dtype=np.float32)
                u.trajectory.ts.dimensions = np.array(BOX, dtype=np.float32)
                writer.write(u.atoms)
        traj = mda.Universe(u._topology, str(path))
        try:
            ligand = md_analyze.select_ligand(traj, "LIG")
            pocket = md_analyze.select_pocket(traj, ligand["atoms"], md_analyze.POCKET_CUTOFF, 0)
            measured = md_analyze.measure(traj, ligand, pocket, first_frame=0)
        finally:
            traj.trajectory.close()
    prot = measured["protein_rmsd"]
    per = prot["per_segment"]
    assert per["PROB"]["range"] > 10 * per["PROA"]["range"], (per["PROA"], per["PROB"])
    assert abs(per["PROA"]["drift"]) > abs(per["PROB"]["drift"]), (per["PROA"], per["PROB"])
    assert prot["segid"] == "PROA", (prot["segid"], per)
    assert "largest DRIFT" in prot["selection"], prot["selection"]
    assert "largest absolute drift" in prot["headline_rule"], prot["headline_rule"]
    fake = {"protein_rmsd": prot, "ligand_rmsd": measured["ligand_rmsd"],
            "com_distance": measured["com_distance"], "hbonds": {"pairs": [], "pairs_found": 0}}
    line = md_analyze.reading(fake)["offered"][0]
    assert "largest drift" in line, line
    assert "moved most" not in line, line
    print(f"ok  PROB swings {per['PROB']['range']} A and PROA drifts "
          f"{per['PROA']['drift']} A; the headline is PROA and the wording says drift")


def test_a_pocket_torn_across_a_periodic_boundary_is_detected():
    """`pocket_split_frames` was in the server docstring as a confidence signal and had no test."""
    import MDAnalysis as mda

    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        folder.mkdir(parents=True, exist_ok=True)
        u = _universe()
        res2 = u.select_atoms("protein and resid 2")
        base = res2.positions.copy()
        path = folder / "step5_1.dcd"
        with mda.Writer(str(path), n_atoms=len(u.atoms), dt=2.0) as writer:
            # frame 1 puts pocket residue 2 a whole box away: the pocket now spans the cell
            for jump in (0.0, 0.0, BOX[0] - 5.0, 0.0):
                res2.positions = base + np.array([jump, 0.0, 0.0], dtype=np.float32)
                u.trajectory.ts.dimensions = np.array(BOX, dtype=np.float32)
                writer.write(u.atoms)
        (folder / "step5_production.inp").write_text(
            NAMD_INPUT.replace("minimize                2000\n", ""))
        _run_topology(folder)
        r = _run(folder, None)
    com = r["com_distance"]
    assert com["pocket_split_frames"] == 1, com["pocket_split_note"]
    assert "torn across a periodic boundary" in com["pocket_split_note"]
    assert any("torn across a periodic boundary" in n for n in r["reading"]["offered"]), \
        r["reading"]["offered"]
    print("ok  1 frame with the pocket spanning more than half a box is counted and said so")


def test_the_same_file_in_two_different_folders_is_read_once():
    """The folder-level dedup cannot catch this: the folders really are different."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        _run_topology(folder)
        nested = folder / "system" / "namd"
        nested.mkdir(parents=True)
        target, link = folder / "step5_1.dcd", nested / "step5_1.dcd"
        # A HARD link, which `os.link` and `mklink /H` make on Windows without any privilege and
        # which `Path.resolve()` does NOT see through: the two paths resolve differently and name
        # one file. That is why the identity is the device and inode, not the spelling.
        import os
        try:
            os.link(target, link)
        except (OSError, NotImplementedError, AttributeError):
            print("ok  (skipped: this filesystem cannot create a hard link)")
            return
        assert link.resolve() != target.resolve(), "the two paths must resolve differently"
        found = md_analyze.find_inputs(folder)
        r = _run(folder, None)
    assert len(found["dcds"]) == 1, [str(d) for d in found["dcds"]]
    assert r["trajectory"]["frames_analysed"] == 4, r["trajectory"]["frames_note"]
    print("ok  a hard-linked .dcd in a second folder is read once: 4 frames, not 8, although the "
          "two paths resolve differently")


def test_a_label_cannot_shift_onto_the_wrong_residue():
    """The CA that goes missing is the FIRST pocket residue's, not the last.

    Truncating a parallel list (`pocket["residues"][:n]`) gives the right answer when the gap is at
    the end, which is what the earlier test happened to exercise. With the gap at the front, a
    parallel list names residue 1 while the value belongs to residue 2.
    """
    import MDAnalysis as mda

    u = _universe()
    keep = u.select_atoms("not (resid 1 and segid PROA and name CA)")
    trimmed = mda.Merge(keep)
    trimmed.add_TopologyAttr("segid", ["PROA", "HETA"])
    trimmed.residues[-1].segment = trimmed.segments[1]
    trimmed.load_new(keep.positions[None, :, :], order="fac")
    ligand = md_analyze.select_ligand(trimmed, "LIG")
    pocket = md_analyze.select_pocket(trimmed, ligand["atoms"], md_analyze.POCKET_CUTOFF, 0)
    measured = md_analyze.measure(trimmed, ligand, pocket, first_frame=0)
    rmsf = measured["rmsf"]
    labels = [e["residue"] for e in rmsf["per_residue"]]
    assert pocket["residues"] == ["PROA:ALA1", "PROA:ALA2"], pocket["residues"]
    # the ONE value belongs to residue 2. A parallel list truncated to one entry says residue 1.
    assert labels == ["PROA:ALA2"], labels
    assert labels != pocket["residues"][:len(labels)], (labels, pocket["residues"])
    assert rmsf["residues_without_ca"] == ["PROA:ALA1"], rmsf["residues_without_ca"]
    print("ok  with the gap at the FRONT of the pocket, the label still names its own atom")


def test_the_hydrogen_bonds_of_the_other_copy_are_really_excluded():
    """`ligand_segid` had to change a BOND, not just a string and a count.

    The earlier two-copy fixture parked the second ligand 15 A from any protein, so it could never
    produce a bond: dropping `segid` from the selection changed nothing measurable. Here each copy
    sits in its own protein, so each one really does make a bond.
    """
    import MDAnalysis as mda

    first = _universe()
    second = _universe()
    second.atoms.positions = second.atoms.positions + np.array([0.0, 0.0, 25.0])
    merged = mda.Merge(first.atoms, second.atoms)
    merged.add_TopologyAttr("segid", ["PROA", "HETA", "PROB", "HETB"])
    for residue, index in zip(merged.residues, (0, 0, 0, 1, 2, 2, 2, 3), strict=True):
        residue.segment = merged.segments[index]

    def ligand_segids(result):
        """The segment of whichever end of each bond carries the ligand."""
        out = set()
        for bond in result["pairs"]:
            end = bond["donor"] if "LIG" in bond["donor"] else bond["acceptor"]
            out.add(end.split(":")[0])
        return out

    both = md_analyze.hydrogen_bonds(merged, "LIG", first_frame=0, frames_analysed=1)
    one = md_analyze.hydrogen_bonds(merged, "LIG", first_frame=0, frames_analysed=1, segid="HETA")
    assert both["pairs_found"] == 2, both["pairs"]
    assert ligand_segids(both) == {"HETA", "HETB"}, both["pairs"]
    assert one["pairs_found"] == 1, one["pairs"]
    assert ligand_segids(one) == {"HETA"}, one["pairs"]
    print(f"ok  both copies make a bond ({both['pairs_found']}); ligand_segid reports "
          f"{one['pairs_found']}, and it is the named copy's")


def test_a_summary_only_reply_carries_no_series_anywhere():
    """`include_series=False` left a full series per protein chain inside `per_segment`."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        _run_topology(folder)
        summary = _run(folder, None, include_series=False)
        full = _run(folder, None, include_series=True)
    def series_keys(result):
        found = []
        def walk(node, path=""):
            if isinstance(node, dict):
                for key, value in node.items():
                    if key in ("values", "times_ps"):
                        found.append(f"{path}.{key}")
                    walk(value, f"{path}.{key}")
            elif isinstance(node, list):
                for item in node:
                    walk(item, path)
        walk(result)
        return found
    assert series_keys(summary) == [], series_keys(summary)
    assert series_keys(full), "include_series=True must still return the series"
    assert summary["protein_rmsd"]["per_segment"], "the per-segment summaries must survive"
    assert summary["protein_rmsd"]["per_segment"]["PROA"]["mean"] is not None
    print(f"ok  include_series=False leaves no series at all; True returns "
          f"{len(series_keys(full))} of them")


def test_a_reference_frame_the_trajectory_does_not_have_is_refused():
    """It is the user's choice, so it comes back as a refusal, not as an IndexError."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        _run_topology(folder)
        far = _run(folder, None, reference_frame=100)
        negative = _run(folder, None, reference_frame=-3)
    for label, r in (("100", far), ("-3", negative)):
        assert r["ok"] is False, (label, r)
        assert "is not a frame of this trajectory" in r["error"], (label, r["error"])
        assert "which has 4" in r["error"], (label, r["error"])
        assert r.get("frames_total") == 4, (label, r)
        assert "IndexError" not in r["error"], (label, r["error"])
        assert "failure of the tool" not in r.get("note", ""), (label, r)
    print("ok  reference_frame 100 and -3 are refused with the frame count, not raised")


def test_two_topologies_of_the_same_name_are_told_apart():
    """A build_md_job package holds step3_input.psf twice; the reply named neither."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False)
        _run_topology(folder)
        nested = folder / "system" / "namd"
        nested.mkdir(parents=True)
        write_psf(nested / "step3_input.psf")
        found = md_analyze.find_inputs(folder)
        r = _run(folder, None)
    selection = found["psf_selection"]
    assert selection == r["topology"]["selection"]
    assert "step3_input.psf; the others: system/namd/step3_input.psf" in selection, selection
    print(f"ok  the two topologies are named by path: {selection}")


# ---------------------------------------------------------------------------------------------
# Fifth review round: the restraint reader's own failure modes
# ---------------------------------------------------------------------------------------------

def _restrained(folder: Path, **kwargs) -> dict:
    """analyze_md on the fixture with a restraint file written by `write_ref(**kwargs)`."""
    write_fixture(folder, with_minimization=False, inp=False)
    write_ref(folder, **kwargs)
    _run_topology(folder)
    return md_analyze.analyze_md(folder)


def test_interleaved_ter_lines_do_not_shift_every_index():
    """The real prot_posres.ref has 5 TER lines, one immediately before the ligand's first atom.

    Counting them as atoms shifts every index after each one. On the real file that produced
    "1659 of 2174 backbone atoms held" instead of 2174 of 2174 - a credible number, with its "N of
    total", invented by an offset of five lines. The fixture now carries a TER in the same place.
    """
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        r = _restrained(folder, hold_ligand=True)
        ref = folder / "restraints" / "prot_posres.ref"
        text = ref.read_text()
        # the fixture must really contain the hazard, or this test proves nothing
        assert "\nTER\nATOM" in text, "the fixture needs a TER in the middle, before the ligand"
        assert text.count("TER") >= 2, text.count("TER")
        # what a reader without the ATOM/HETATM filter would have counted
        unfiltered = sum(1 for line in text.splitlines() if line.strip())
        filtered = sum(1 for line in text.splitlines() if line.startswith(("ATOM", "HETATM")))
        assert unfiltered > filtered, (unfiltered, filtered)
    lig = r["trajectory"]["restraints"]["held"]["ligand_heavy"]
    assert lig["known"] is True and lig["held"] == lig["of"] == 4, lig
    assert r["trajectory"]["restraints"]["held"]["protein_backbone"]["held"] == 12, \
        r["trajectory"]["restraints"]["held"]["protein_backbone"]
    assert f"of {filtered} atoms" in r["trajectory"]["restraints"]["reference_file_note"], \
        r["trajectory"]["restraints"]["reference_file_note"]
    print(f"ok  {unfiltered - filtered} non-ATOM lines are skipped, so the ligand still reads "
          f"{lig['held']} of {lig['of']} and not a shifted count")


def test_a_reference_file_for_another_system_is_unknown_not_zero():
    """`held: 0` reads as "the ligand was free" and removes the warning. It must not be guessable.

    Two ways a file can be the wrong one: the wrong number of atoms, and the right number in
    another order (a posres from a regenerated CHARMM-GUI build). The second is the dangerous one,
    because the length check passes and every index lands on another atom.
    """
    with tempfile.TemporaryDirectory() as tmp:
        bigger, reordered, right = Path(tmp) / "b", Path(tmp) / "r", Path(tmp) / "ok"
        too_big = _restrained(bigger, extra_atoms=7)
        shifted = _restrained(reordered, shift=3)
        good = _restrained(right)

    big = too_big["trajectory"]["restraints"]["held"]["ligand_heavy"]
    assert big["known"] is False, big
    assert big["atoms_in_the_reference_file"] == big["atoms_in_the_system"] + 7, big
    assert "NOT the file this trajectory ran with" in big["note"], big["note"]

    mix = shifted["trajectory"]["restraints"]["held"]["ligand_heavy"]
    assert mix["known"] is False, mix
    assert "NOT in this system's order" in mix["note"], mix["note"]
    assert mix["mismatches"], mix
    # the trap: the atoms it WOULD have reported are not held, so the old code said exactly this
    assert "held" not in mix, mix

    # and when the file is right, the answer is a number again
    assert good["trajectory"]["restraints"]["held"]["ligand_heavy"]["known"] is True
    assert good["trajectory"]["restraints"]["held"]["ligand_heavy"]["alignment_checked"] is True

    # with the answer unknown, the ligand lines must NOT be left clean either: the fallback is the
    # keyword, which still says something was held
    for label, result in (("bigger", too_big), ("reordered", shifted)):
        marked = [n for n in result["reading"]["offered"] if "CARRIES NO EVIDENCE" in n]
        assert len(marked) == 3, (label, [n[:70] for n in result["reading"]["offered"]])
    print("ok  a longer file and a reordered file are both UNKNOWN, not 0, and the lines stay "
          "marked from the keyword")


def test_an_assumed_force_constant_column_is_not_presented_as_read():
    """With no `conskcol`, B is assumed because CHARMM-GUI writes it. That is a guess, and says so.

    It matters on the real system: the occupancy column is 1.00 in all 139989 atoms, so reading O
    instead of B would mark the entire box. The counts under an assumed column are unverified.
    """
    with tempfile.TemporaryDirectory() as tmp:
        declared, assumed, occupancy = Path(tmp) / "d", Path(tmp) / "a", Path(tmp) / "o"
        said = _restrained(declared, column="B", declare_column=True)
        guessed = _restrained(assumed, column="B", declare_column=False)
        # the same force constants written in the OCCUPANCY column, and declared as such
        in_o = _restrained(occupancy, column="O", declare_column=True)

    assert "read from `conskcol`" in said["trajectory"]["restraints"][
        "force_constant_column_source"], said["trajectory"]["restraints"]
    source = guessed["trajectory"]["restraints"]["force_constant_column_source"]
    assert "ASSUMED" in source and "unverified" in source, source
    assert guessed["trajectory"]["restraints"]["force_constant_column"] == "B"
    # a column other than B must actually be read and used, or _PDB_COLUMNS is decoration
    assert in_o["trajectory"]["restraints"]["force_constant_column"] == "O", \
        in_o["trajectory"]["restraints"]
    assert in_o["trajectory"]["restraints"]["held"]["ligand_heavy"]["held"] == 4, \
        in_o["trajectory"]["restraints"]["held"]
    assert "column 55-60 (O)" in in_o["trajectory"]["restraints"]["reference_file_note"], \
        in_o["trajectory"]["restraints"]["reference_file_note"]
    print("ok  a declared column is reported as read, an undeclared one as ASSUMED and "
          "unverified, and column O is really read")


def test_both_a_spring_and_fixed_atoms_are_reported_and_the_frozen_one_leads():
    """Freezing is the stronger statement, so its file is the one the atom counts come from."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, inp=False)
        write_ref(folder, hold_ligand=False, name="spring.ref")          # spring: protein only
        write_ref(folder, hold_ligand=True, name="frozen.ref", keyword="fixedAtoms")
        (folder / "step5_production.inp").write_text(
            "dcdfreq                 1000;\n"
            "constraints             on\n"
            "conskfile               restraints/spring.ref\n"
            "conskcol                B\n"
            "fixedAtoms              on\n"
            "fixedAtomsFile          restraints/frozen.ref\n"
            "fixedAtomsCol           B\n"
            "run                     4000\n")
        _run_topology(folder)
        r = md_analyze.analyze_md(folder)
    held = r["trajectory"]["restraints"]
    kinds = {d["kind"] for d in held["declarations"]}
    assert len(held["declarations"]) == 2, held["declarations"]
    assert kinds == {"spring (harmonic restraint)", "frozen (fixed atoms)"}, kinds
    assert "frozen" in held["kind"], held["kind"]
    assert held["reference_file"] == "restraints/frozen.ref", held["reference_file"]
    assert "stronger statement" in held["declarations_note"], held["declarations_note"]
    # the frozen file holds the ligand, the spring file does not: taking the spring would have
    # reported the ligand as free while it was frozen
    assert held["held"]["ligand_heavy"]["held"] == 4, held["held"]["ligand_heavy"]
    assert any("fixedAtoms              on" in line for line in held["evidence"]), held["evidence"]
    assert any("constraints             on" in line for line in held["evidence"]), held["evidence"]
    print(f"ok  both holds are listed ({len(held['declarations'])}), the frozen one leads, and the "
          "ligand reads 4 held instead of 0")


def test_the_force_constant_file_wins_over_the_coordinate_file():
    """NAMD takes the constants from `conskfile` and only the coordinates from `consref`."""
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder, with_minimization=False, inp=False)
        write_ref(folder, hold_ligand=False, name="coords.ref")      # consref: ligand NOT held
        write_ref(folder, hold_ligand=True, name="constants.ref")    # conskfile: ligand held
        (folder / "step5_production.inp").write_text(
            "dcdfreq                 1000;\n"
            "constraints             on\n"
            "consref                 restraints/coords.ref;\n"
            "conskfile               restraints/constants.ref;\n"
            "conskcol                B;\n"
            "run                     4000\n")
        _run_topology(folder)
        r = md_analyze.analyze_md(folder)
    held = r["trajectory"]["restraints"]
    # the trailing `;` is CHARMM-GUI's style, and an unstripped one makes the path unopenable
    assert held["reference_file"] == "restraints/constants.ref", held["reference_file"]
    assert ";" not in held["reference_file"], held["reference_file"]
    assert held["force_constant_column"] == "B", held["force_constant_column"]
    assert held["held"]["ligand_heavy"]["held"] == 4, held["held"]["ligand_heavy"]
    print("ok  conskfile beats consref and the trailing semicolons are stripped: the ligand reads "
          "4 held, not 0")


def test_a_reference_file_with_nothing_in_it_is_unknown():
    """An absent, empty or oddly-named column file must make the answer unknown, not zero."""
    with tempfile.TemporaryDirectory() as tmp:
        empty, missing, odd = Path(tmp) / "e", Path(tmp) / "m", Path(tmp) / "c"
        no_atoms = _restrained(empty, atoms=False)

        write_fixture(missing, with_minimization=False, inp=False)
        (missing / "step5_production.inp").write_text(
            "dcdfreq 1000;\nconstraints on\nconskfile restraints/gone.ref\nconskcol B\n"
            "run 4000\n")
        _run_topology(missing)
        absent = md_analyze.analyze_md(missing)

        write_fixture(odd, with_minimization=False, inp=False)
        write_ref(odd, column="B")
        _run_topology(odd)
        (odd / "step5_production.inp").write_text(
            (odd / "step5_production.inp").read_text().replace("conskcol                B",
                                                               "conskcol                Q"))
        bad_column = md_analyze.analyze_md(odd)

    for label, result, expected in (("empty", no_atoms, "holds no ATOM or HETATM line"),
                                    ("missing", absent, "could not be read (FileNotFoundError)"),
                                    ("column Q", bad_column, "not one of X, Y, Z, O or B")):
        held = result["trajectory"]["restraints"]
        assert held["atoms_held"] is True, (label, held)
        assert expected in held["reference_file_note"], (label, held["reference_file_note"])
        assert held["held"]["ligand_heavy"]["known"] is False, (label, held["held"])
        assert "UNKNOWN" in held["held"]["ligand_heavy"]["note"], (label, held["held"])
        # the keyword alone still says something was held, so the lines stay marked
        assert sum("CARRIES NO EVIDENCE" in n for n in result["reading"]["offered"]) == 3, \
            (label, [n[:70] for n in result["reading"]["offered"]])
    print("ok  an empty file, a missing file and an unknown column are all UNKNOWN with the lines "
          "still marked")


def test_the_hydrogen_bond_line_is_not_left_clean_when_the_atoms_are_held():
    """It was the only figure an analyst could quote believing it unqualified."""
    with tempfile.TemporaryDirectory() as tmp:
        restrained, free = Path(tmp) / "h", Path(tmp) / "f"
        held = _restrained(restrained, hold_ligand=True)
        write_fixture(free, with_minimization=False)
        _run_topology(free)
        loose = md_analyze.analyze_md(free)
    hbond_line = [n for n in held["reading"]["offered"] if "hydrogen-bond pairs" in n]
    assert len(hbond_line) == 1, held["reading"]["offered"]
    assert "upper bound rather than a result" in hbond_line[0], hbond_line[0]
    assert "an occupancy is not a displacement" in hbond_line[0], hbond_line[0]
    clean = [n for n in loose["reading"]["offered"] if "hydrogen-bond pairs" in n]
    assert "upper bound" not in clean[0], clean[0]
    print("ok  the hydrogen-bond line is qualified when the atoms are held and left alone when "
          "they are not")


# ---------------------------------------------------------------------------------------------
# The topology, written as a real PSF
# ---------------------------------------------------------------------------------------------

# MDAnalysis cannot write a PSF, and the formats it can write lose exactly what this test needs:
# a PDB keeps no charges and no bonds, and hydrogen-bond detection reads both. So the fixture
# writes the PSF itself, in the XPLOR layout CHARMM-GUI produces
# (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8), which also means `find_inputs` is exercised on
# a real `.psf` instead of being handed a path.
CHARGES = {"N": -0.47, "H": 0.31, "C": 0.07, "O": -0.51, "CL": -0.15, "LP": 0.0}


def write_psf(path: Path) -> Path:
    atoms, bonds = [], []
    for segid, resid, resname, ats in PROTEIN:
        for name, typ, mass, _ in ats:
            atoms.append((segid, str(resid), resname, name, typ, CHARGES[typ], mass))
    for name, typ, mass, _ in LIGAND:
        atoms.append(("HETA", "1", "LIG", name, typ, CHARGES[typ], mass))
    for r in range(len(PROTEIN)):
        b = r * 5
        bonds += [(b, b + 1), (b, b + 2), (b + 2, b + 3), (b + 3, b + 4)]
    lig0 = len(PROTEIN) * 5
    bonds += [(lig0 + a, lig0 + b) for a, b in LIGAND_BONDS]

    lines = ["PSF", "", "       1 !NTITLE",
             " REMARKS fixture for test_md_analyze", "", f"{len(atoms):8d} !NATOM"]
    for i, (segid, resid, resname, name, typ, charge, mass) in enumerate(atoms, 1):
        lines.append(f"{i:8d} {segid:<4s} {resid:<4s} {resname:<4s} {name:<4s} {typ:<4s} "
                     f"{charge:14.6f}{mass:14.6f}{0:8d}")
    lines += ["", f"{len(bonds):8d} !NBOND: bonds"]
    row = []
    for a, b in bonds:
        row += [f"{a + 1:8d}", f"{b + 1:8d}"]
        if len(row) == 8:
            lines.append("".join(row))
            row = []
    if row:
        lines.append("".join(row))
    lines.append("")
    path.write_text("\n".join(lines))
    return path


def _run_topology(folder: Path) -> Path:
    path = folder / "step3_input.psf"
    return path if path.exists() else write_psf(path)


def _expected_com_distances():
    """|d0 + shift * x| for each dynamics shift, straight from the fixture's coordinates.

    d0 is the vector from the pocket's centre of mass to the ligand's, both taken exactly as the
    module takes them: the ligand's HEAVY atoms against ALL atoms of the pocket residues.
    """
    u = _universe()
    lig = md_analyze.select_ligand(u, "LIG")
    pocket = md_analyze.select_pocket(u, lig["atoms"], md_analyze.POCKET_CUTOFF, 0)
    d0 = lig["atoms"].center_of_mass() - pocket["atoms"].center_of_mass()
    return [float(np.linalg.norm(d0 + np.array([s, 0.0, 0.0]))) for s in SHIFTS_DYNAMICS]


def _run(folder: Path, topology: Path, **kwargs) -> dict:
    """analyze_md on the fixture folder: the PSF is found by the module, not handed to it."""
    return md_analyze.analyze_md(folder, **kwargs)


def main():
    test_the_fixture_is_what_it_claims()
    test_minimization_frames_are_read_from_the_input()
    test_an_input_without_minimization_claims_none()
    test_heavy_atoms_by_name_would_count_the_lone_pair()
    test_the_pocket_is_measured_not_named()
    test_a_ligand_that_is_not_there_is_refused()
    test_the_metrics_match_the_analytic_answer()
    test_a_fit_on_the_ligand_would_hide_the_motion()
    test_keeping_the_minimization_frames_changes_the_answer()
    test_a_reference_frame_without_a_pocket_is_refused()
    test_mixed_frame_spacing_is_refused()
    test_the_chunks_are_ordered_by_number_not_by_text()
    test_a_folder_without_a_trajectory_says_so()
    test_the_hydrogen_bond_occupancy_is_a_sample_of_frames()
    test_a_ligand_that_leaves_is_not_called_stable()
    test_the_reading_is_marked_as_offered_and_names_its_thresholds()
    test_nothing_claims_an_affinity()
    test_the_result_survives_json()
    # the mutations the first round let through
    test_a_rotating_protein_needs_the_rotation()
    test_a_wrapped_ligand_is_not_called_motion()
    test_a_monotonic_drift_is_not_a_plateau()
    test_heavy_atoms_are_chosen_by_mass_not_by_any_name_rule()
    test_mixed_frame_spacing_is_caught_in_a_middle_chunk()
    test_minimization_frames_round_down()
    test_the_occupancy_denominator_is_the_analysed_frames()
    test_the_equilibration_is_not_chained_onto_the_production_chunks()
    test_a_production_run_written_in_one_piece_is_not_called_another_stage()
    test_a_named_stage_is_analysed_and_an_absent_one_is_refused()
    test_an_unnumbered_trajectory_sorts_first()
    test_a_refusal_comes_back_as_ok_false()
    test_the_time_axis_starts_at_the_first_analysed_frame()
    test_the_hbond_criterion_does_not_mention_a_charge_and_says_what_was_eligible()
    # the guards from the second review round
    test_the_thresholds_match_their_pre_registration()
    test_two_copies_of_the_ligand_are_refused_not_merged()
    test_an_unreadable_chunk_is_refused_instead_of_leaving_a_gap()
    test_every_rmsf_value_is_labelled_from_its_own_atom()
    # the third review round: the regression above all
    test_a_ligand_that_really_leaves_is_not_folded_back_into_the_site()
    test_a_wrapped_ligand_and_a_walking_one_are_told_apart()
    test_two_protein_chains_are_never_superposed_together()
    test_the_hydrogen_bonds_follow_the_named_ligand_copy()
    test_hydrogen_bond_labels_survive_ids_that_are_not_indices()
    test_a_hydrogen_on_a_carbon_is_never_a_donor()
    test_a_pocket_residue_without_a_ca_is_named_and_not_mislabelled()
    test_a_triclinic_box_is_not_treated_as_rectangular()
    test_a_frame_without_a_box_is_not_called_rectangular()
    test_the_span_is_the_elapsed_time_between_the_first_and_the_last_frame()
    test_the_hydrogen_bonds_start_at_the_same_frame_as_everything_else()
    test_positional_restraints_in_the_input_are_declared()
    test_no_input_for_this_stage_is_unknown_and_not_absent()
    test_the_same_trajectory_reached_twice_is_read_once()
    test_the_reference_frame_moves_the_rmsd_zero_as_well_as_the_pocket()
    # the fourth round: the restraint file, and the mutations the third round let through
    test_the_restraint_file_is_read_and_says_which_atoms_were_held()
    test_fixed_atoms_are_looked_for_and_not_merely_claimed()
    test_namd_keywords_are_read_case_insensitively()
    test_too_few_frames_cannot_support_a_plateau_and_says_so()
    test_the_rmsf_names_the_most_flexible_residue_and_is_taken_about_the_mean()
    test_a_rigid_translation_of_everything_moves_no_metric()
    test_the_angle_cutoff_decides_whether_a_bond_is_counted()
    test_the_ligand_side_of_each_bond_is_named_correctly()
    test_the_per_frame_bond_count_is_aligned_with_the_analysed_frames()
    test_the_headline_segment_is_the_largest_drift_and_says_so()
    test_a_pocket_torn_across_a_periodic_boundary_is_detected()
    test_the_same_file_in_two_different_folders_is_read_once()
    test_a_label_cannot_shift_onto_the_wrong_residue()
    test_the_hydrogen_bonds_of_the_other_copy_are_really_excluded()
    test_a_summary_only_reply_carries_no_series_anywhere()
    test_a_reference_frame_the_trajectory_does_not_have_is_refused()
    test_two_topologies_of_the_same_name_are_told_apart()
    # the fifth round: the restraint reader's own failure modes
    test_interleaved_ter_lines_do_not_shift_every_index()
    test_a_reference_file_for_another_system_is_unknown_not_zero()
    test_an_assumed_force_constant_column_is_not_presented_as_read()
    test_both_a_spring_and_fixed_atoms_are_reported_and_the_frozen_one_leads()
    test_the_force_constant_file_wins_over_the_coordinate_file()
    test_a_reference_file_with_nothing_in_it_is_unknown()
    test_the_hydrogen_bond_line_is_not_left_clean_when_the_atoms_are_held()
    print("ALL OK")


if __name__ == "__main__":
    main()
