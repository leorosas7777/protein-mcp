"""Offline test for build_console_page (protein_mcp/console.py, Phase 3 step 10).

The page is the one artefact of this project a reader can click through instead of reading, so
the only thing it may never do is add something the dossier does not hold. Everything here is
built around that: the embedded payload is parsed back out of the HTML and compared field by
field against the fixture, and every path that reads a file from disk is checked with the file
changed under it.

The fixture copies the shapes of the real BRAF V600E dossier (2026-11-18): 8C7Y with the TXV U
pocket, one T1 compound run with its cross_docking block and one crystallographic re-docking,
an ADMET block, and an md_analysis stage with a series file. The numbers are made up so that a
value printed from the wrong field is visible.

Controls that must fail (each is run with the condition inverted inside the test):
- a pose SDF or a series file whose SHA-256 no longer matches the manifest is refused, and the
  page says so instead of drawing the old coordinates;
- a file whose manifest path leaves prep/ is not embedded;
- only the FIRST record of a multi-pose SDF reaches the page;
- the compound of a pocket key is looked up in ligand_choice, never taken from a name;
- per-variant counts stay null: the page never fills in a number nobody retrieved.

Run: uv run python tests/test_console.py
"""
import hashlib
import json
import re
import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from protein_mcp import console, dossier  # noqa: E402

ACC = "P15056"
POCKET = "TXV U"
RUN_T1 = "TXV_regenerated_chembl1229517"
RUN_REDOCK = "TXV_regenerated"

# A hand-written SDF: the coordinates are literals so the shift to the box centre is checked
# against arithmetic done here, not against the same code that produced it. Two records, so the
# "only the first pose is embedded" rule has something to fail on.
def sdf_record(title: str, x0: float) -> str:
    atoms = [("C", x0, 20.0, 30.0), ("N", x0 + 1.0, 20.0, 30.0), ("O", x0 + 2.0, 20.0, 30.0)]
    lines = [title, "  protein-mcp test", "",
             f"{len(atoms):3d}{2:3d}  0  0  0  0  0  0  0  0999 V2000"]
    for sym, x, y, z in atoms:
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {sym:<3s} 0  0  0  0  0  0  0  0  0  0  0  0")
    lines += ["  1  2  1  0", "  2  3  1  0", "M  END", "$$$$"]
    return "\n".join(lines) + "\n"


POSES_SDF = sdf_record("pose 1 of this run", 10.0) + sdf_record("SECOND_RECORD_MARKER", 90.0)
BOX_CENTRE = [0.02, 110.789, 15.175]

SERIES = {"tool": "analyze_md", "result": {
    "protein_rmsd": {"per_segment": {
        "PROA": {"values": [0.1, 1.35], "times_ps": [0.0, 10000.0], "min": 0.1, "max": 1.35,
                 "drift": 1.3534,
                 "observed_only": {"values": [0.1, 0.42], "times_ps": [0.0, 10000.0],
                                   "min": 0.1, "max": 0.42, "drift": 0.425}}}},
    "ligand_rmsd": {"values": [0.2, 1.9], "times_ps": [0.0, 10000.0], "min": 0.2, "max": 1.9},
    "com_distance": {"values": [3.0, 3.2], "times_ps": [0.0, 10000.0], "min": 3.0, "max": 3.2},
    # the real defect of the data: the H-bond counts carry no time axis of their own
    "hbonds": {"count_per_frame": {"values": [4, 2]},
               "pairs": [{"donor": "LIG N2", "acceptor": "GLN530 O", "frames_present": 22,
                          "frames_analysed": 25, "occupancy_sampled": 0.88}]},
    "rmsf": {"per_residue": [{"residue": "GLU501", "rmsf_angstrom": 0.91}]},
    "reading": {"equilibration_control": {"status": "NOT MET", "drift_angstrom": 1.3534,
                                          "threshold_angstrom": 0.5, "segid": "PROA"},
                "offered": ["plateau", "ligand_rmsd"]},
    "trajectory": {"span_ns": 10.0, "frames_note": "100 of 100 frames analysed",
                   "restraints": {"note": "the production input declares no restraints"}},
}}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def manifest(prep: Path) -> dict:
    poses = {"mode": 1, "vina_score_kcal_mol": -11.191, "rmsd_to_crystal_pose_angstrom": None,
             "rmsd_to_reference_pose_angstrom": 1.9152, "rmsd_reference": "3OG7 032 C"}
    files = {"poses_TXV_regenerated_chembl1229517.sdf": {
        "path": f"dock/8C7Y_TXV_U/poses_{RUN_T1}.sdf",
        "sha256": sha(POSES_SDF.encode())}}
    t1 = {
        "ligand": {"compound_id": "CHEMBL1229517", "name": "VEMURAFENIB", "level": "T1",
                   "slug": "chembl1229517", "inchikey": "GPXBXXGIAQBQNI-UHFFFAOYSA-N",
                   "smiles_as_recorded": "CCCS(=O)(=O)Nc1ccc(F)c1F"},
        "compound_docking": {"compound_id": "CHEMBL1229517", "name": "VEMURAFENIB", "level": "T1",
                             "top_pose_score_kcal_mol": -11.191, "heavy_atoms": 33,
                             "ligand_efficiency_kcal_mol_per_heavy_atom": -0.3391,
                             "box_limited": False,
                             "box_fit": {"min_clearance_angstrom": 4.2,
                                         "limiting_box_axis_angstrom": 22.0,
                                         "limiting_ligand_extent_angstrom": 17.8},
                             "scope": "one compound in one pocket of one receptor"},
        "cross_docking": {"top_pose_rmsd_angstrom": 1.9152, "criterion_angstrom": 2.0,
                          "reference_pose": {"from": "3OG7 032 C",
                                             "matched_by": "InChIKey of the crystallographic SDF"},
                          "reference_pose_uncertainty": {"shift_angstrom": 0.767},
                          "within_criterion_with_bar": {"range_angstrom": [1.148, 2.682],
                                                        "decides": False},
                          "superposition": {"rmsd_angstrom": 1.671, "pairs": 241}},
        "poses": [poses], "files": files, "vina": {"wall_seconds": 234.8},
    }
    redock = {
        "ligand": {"comp_id": "TXV", "name": "TXV"},
        "redocking": {"top_pose_rmsd_angstrom": 1.053, "criterion_angstrom": 2.0,
                      "within_criterion": True, "scope": "re-docking of the crystal ligand"},
        "poses": [{"mode": 1, "vina_score_kcal_mol": -10.4,
                   "rmsd_to_crystal_pose_angstrom": 1.053}],
        "files": {"poses_TXV_regenerated.sdf": {"path": f"dock/8C7Y_TXV_U/poses_{RUN_REDOCK}.sdf",
                                                "sha256": sha(POSES_SDF.encode())}},
        "vina": {"wall_seconds": 170.0},
    }
    structure = {
        "pdb_id": "8C7Y", "title": "BRAF V600E with TXV", "resolution_angstrom": 2.1,
        "prepared": {"ok": True, "receptor": {"path": "clean/8C7Y_receptor.pdb"},
                     "pockets": [{"comp_id": "TXV", "label_asym_id": "U", "auth_asym_id": "A",
                                  "heavy_atoms": 30, "large_box": False,
                                  "spans_target_chains": False,
                                  "box": {"center": BOX_CENTRE, "size": [21.8, 16.27, 22.0]},
                                  "residues": [
                                      {"resname": "GLU", "auth_seq_id": 501, "min_distance": 3.1},
                                      {"resname": "TRP", "auth_seq_id": 531, "min_distance": 4.4}]}]},
        "docking": {POCKET: {"box": {"center": BOX_CENTRE}, "runs": {RUN_T1: t1,
                                                                    RUN_REDOCK: redock}}},
        "admet": {f"{POCKET} CHEMBL1229517": {
            "compound_identity": {"compound_id": "CHEMBL1229517", "name": "VEMURAFENIB"},
            "chembl_record": {"chembl_id": "CHEMBL1229517", "pref_name": "VEMURAFENIB",
                              "link": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL1229517/",
                              "status": "approved drug",
                              "fields": {"max_phase": 4.0, "first_approval": 2011,
                                         "withdrawn_flag": False, "black_box_warning": False},
                              "atc_classifications": ["L01EC01"]},
            "calculated_properties": {"rdkit": {"molecular_weight": 489.9, "clogp": 5.544,
                                                "tpsa": 92.4, "hbd": 2, "hba": 5,
                                                "rotatable_bonds": 6, "qed": 0.33,
                                                "heavy_atoms": 33}},
            "measured_adme": {"retrieved": True, "rows_read": 41,
                              "adme_rows": [{"adme_class": "clearance", "standard_type": "CL",
                                             "value": 6.61, "relation": "=",
                                             "units": "ml kg-1", "assay_description": "dog"}]},
            "rule_reading": {"results": {"Lipinski": {"violations": 1, "violated": ["cLogP"]},
                                         "Veber": {"violations": 0, "violated": []}}},
        }},
        "md_package": {f"{POCKET} CHEMBL1229517": {
            "run_id": RUN_T1, "pose": {"mode": 1, "vina_score_kcal_mol": -11.191,
                                       "rmsd_to_crystal_pose_angstrom": None}}},
        "md_job": {f"{POCKET} CHEMBL1229517": {
            "verified": {"checks": {"system": {"box_side_angstrom": 114.0, "waters": 41000,
                                               "temperature_kelvin": 303.15, "timestep_fs": 2.0}}},
            "slurm": {"nanoseconds_total": 10.0}}},
        "md_analysis": {f"{POCKET} CHEMBL1229517": {"production": {
            "result": SERIES["result"],
            "series_file": {"path": "md/8C7Y_TXV_U_chembl1229517/analysis_production.json",
                            "sha256": sha(json.dumps(SERIES).encode())}}}},
    }
    return {
        "generated_by": "prepare_structures", "uniprot_id": ACC,
        "created_at": "2026-11-18T17:00:00Z", "mutation_filter": "V600E",
        "user_choice": {"structure": "8C7Y", "resolution_angstrom": 2.1,
                        "reason": "the most recent inhibitor complex"},
        "ligand_choice": {"reason": "five candidates and the anchor",
                          "compounds": [{"id": "CHEMBL1229517", "name": "VEMURAFENIB",
                                         "level": "T1"},
                                        {"id": "CHEMBL3577124", "name": "LY-3009120",
                                         "level": "T3"}]},
        "structures": [structure],
    }


PROFILE = {"tool": "get_protein_profile", "retrieved_at": "2026-11-01T10:00:00Z", "result": {
    "uniprot": {"protein_name": "Serine/threonine-protein kinase B-raf", "gene": "BRAF",
                "organism": "Homo sapiens", "length_aa": 766, "entry_name": "BRAF_HUMAN",
                "source_url": "https://www.uniprot.org/uniprotkb/P15056",
                "function": [{"text": "Protein kinase involved in the MAPK pathway."}],
                "natural_variants": [{"change": "V -> E", "start": 600,
                                      "description": "in a melanoma sample", "dbsnp": ["rs113488022"]},
                                     {"change": "G -> R", "start": 466, "description": ""}]},
    "open_targets": {"ensembl_gene_id": "ENSG00000157764", "data_version": "25.06",
                     "associations": {"total": 1317, "examined": 500,
                                      "selection": "10 of 1317",
                                      "rows": [{"disease_id": "MONDO_0005105", "name": "melanoma",
                                                "score": 0.81, "is_oncology": True,
                                                "datatype_scores": {"genetic_association": 0.7}}]},
                     "tractability": {"SM": True}, "modality_names": {"SM": "Small molecule"}}}}

STRUCTURES = {"tool": "get_structures", "retrieved_at": "2026-11-01T10:01:00Z", "result": {
    "experimental": {"total_entries": 162, "examined": 100, "confirmed_in_examined": 44,
                     "selection": "7 of 162",
                     "structures": [{"pdb_id": "8C7Y", "title": "BRAF V600E with TXV",
                                     "resolution_angstrom": 2.1, "release_date": "2023-02-01",
                                     "method": "X-RAY DIFFRACTION",
                                     "ligands": [{"id": "TXV"}],
                                     "links": {"rcsb": "https://www.rcsb.org/structure/8C7Y"}}],
                     "discarded_not_confirmed": [{"pdb_id": "9GLU", "title": "KRasG12D-GDP"}]},
    "alphafold": {"entry_id": "AF-P15056-F1", "mean_plddt": 79.2, "confidence_band": "confident",
                  "fraction_of_residues": {"very_high": 0.42}}}}

LIGANDS = {"tool": "get_ligands", "retrieved_at": "2026-11-01T10:02:00Z", "result": {
    "selection": "15 of 1740", "ligands_with_qualifying_data": 348,
    "ligands_with_qualifying_data_is_lower_bound": True,
    "target": {"chembl_id": "CHEMBL5145"},
    "activities": {"mechanism_compounds_download": "complete"},
    "ligands": [{"chembl_id": "CHEMBL1229517", "name": "VEMURAFENIB", "max_phase": 4.0,
                 "first_approval": 2011, "best_pchembl": 8.5, "n_measurements": 63,
                 "link": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL1229517/",
                 "mechanism": {"mechanism_of_action": "Serine/threonine-protein kinase B-raf inhibitor"},
                 "identity": {"molecular_weight": 489.9, "inchi_key": "GPXBXXGIAQBQNI-UHFFFAOYSA-N"},
                 "properties": {"alogp": 5.54, "psa": 92.4, "ro5_violations": 1},
                 "potency": [{"mutation": "V600E", "standard_type": "IC50", "assay_class": "biochemical",
                              "n_measurements": 12, "best_nM": 13.0, "median_nM": 31.0,
                              "worst_nM": 59.0, "assay_examples": ["Inhibition of BRAF V600E"]}]}]}}

LITERATURE = {"tool": "search_literature", "retrieved_at": "2026-11-01T10:03:00Z", "result": {
    "query_used": "BRAF AND inhibitor", "total_hits": 4210,
    "papers": [{"title": "Clinical efficacy of a RAF inhibitor", "journal": "Nature", "year": 2010,
                "pmid": "20823850", "doi": "10.1038/nature09454", "open_access": False,
                "abstract": "A phase 1 trial of PLX4032 in metastatic melanoma."}]}}


def build_dossier(root: Path) -> Path:
    """dossiers/P15056 with the four dossier files, the manifest, and the files the page reads."""
    dossier.DOSSIER_DIR = root / "dossiers"
    folder = dossier.DOSSIER_DIR / ACC
    prep = folder / "prep"
    (prep / "clean").mkdir(parents=True, exist_ok=True)
    (prep / "dock" / "8C7Y_TXV_U").mkdir(parents=True, exist_ok=True)
    (prep / "md" / "8C7Y_TXV_U_chembl1229517").mkdir(parents=True, exist_ok=True)
    for name, doc in (("profile.json", PROFILE), ("structures_v600e.json", STRUCTURES),
                      ("ligands_v600e.json", LIGANDS),
                      ("literature_v600e-inhibitor.json", LITERATURE)):
        (folder / name).write_text(json.dumps(doc), encoding="utf-8")
    (prep / "clean" / "8C7Y_receptor.pdb").write_text("ATOM  RECEPTOR_TEXT_8C7Y\nEND\n",
                                                      encoding="utf-8")
    # written as bytes, and hashed as the same bytes: write_text would translate the newlines
    # on Windows and the SHA-256 in the manifest would never match the file on disk
    for run in (RUN_T1, RUN_REDOCK):
        (prep / "dock" / "8C7Y_TXV_U" / f"poses_{run}.sdf").write_bytes(POSES_SDF.encode())
    (prep / "md" / "8C7Y_TXV_U_chembl1229517" / "analysis_production.json").write_bytes(
        json.dumps(SERIES).encode())
    (prep / "manifest.json").write_text(json.dumps(manifest(prep)), encoding="utf-8")
    return folder


def payload(html: str) -> dict:
    """The data the page carries, parsed back out of the HTML."""
    line = next(l for l in html.splitlines() if l.startswith("const DATA = "))
    return json.loads(line[len("const DATA = "):].rstrip(";").replace("<\\/", "</"))


def run(folder: Path):
    res = console.build_console_page(ACC)
    assert res["ok"], res
    html = Path(res["file"]).read_text(encoding="utf-8")
    return res, html, payload(html)


# --------------------------------------------------------------------------- the page itself
def test_every_value_in_the_page_comes_from_the_dossier(folder: Path):
    res, html, data = run(folder)
    assert res["file"] == str(folder / f"console_{ACC}.html")
    assert res["holds"] == {"associations": 1, "variants": 2, "ligands_retrieved": 1,
                            "docking_runs": 2, "admet_results": 1, "md_pockets": 1,
                            "md_stages_measured": 1}, res["holds"]

    ctx = data["context"]
    assert ctx["identity"]["gene"] == "BRAF" and ctx["identity"]["length"] == 766
    assert ctx["assoc"]["total"] == 1317 and ctx["assoc"]["selection"] == "10 of 1317"
    assert ctx["structures"]["total"] == 162 and ctx["structures"]["rows"][0]["pdb"] == "8C7Y"
    assert ctx["ligands"]["rows"][0]["potency"][0]["median"] == 31.0
    assert ctx["literature"]["rows"][0]["pmid"] == "20823850"
    assert ctx["retrieved"]["ligands"] == "2026-11-01T10:02:00Z"

    exp = data["experiments"]
    runs = {r["run"]: r for r in exp["runs"]}
    assert set(runs) == {RUN_T1, RUN_REDOCK}
    t1 = runs[RUN_T1]
    assert t1["kind"] == "compound" and t1["level"] == "T1" and t1["score"] == -11.191
    assert t1["le"] == -0.3391 and t1["heavy"] == 33 and t1["clearance"] == 4.2
    assert t1["anchor"]["rmsd"] == 1.9152 and t1["anchor"]["bar"] == 0.767
    assert t1["anchor"]["decides"] is False and t1["anchor"]["range"] == [1.148, 2.682]
    assert t1["poses"][0]["rmsd"] == 1.9152, "a T1 run must print the RMSD to the reference pose"
    # the re-docking keeps the crystal RMSD, and carries no anchor block at all
    assert runs[RUN_REDOCK]["anchor"] is None and runs[RUN_REDOCK]["redock_rmsd"] == 1.053
    assert runs[RUN_REDOCK]["poses"][0]["rmsd"] == 1.053
    assert runs[RUN_REDOCK]["score"] is None, "a re-docking has no compound_docking score"

    assert exp["admet"][0]["props"]["clogp"] == 5.544
    assert exp["admet"][0]["rules"]["Lipinski"] == {"violations": 1, "violated": ["cLogP"]}
    assert exp["admet"][0]["adme_rows"][0]["units"] == "ml kg-1"
    md = exp["md"][0]
    assert md["compound"] == "VEMURAFENIB" and md["level"] == "T1" and md["ns_total"] == 10.0
    assert md["box_side"] == 114.0 and md["pose_score"] == -11.191
    stage = md["stages"][0]
    assert stage["label"] == "production" and stage["status"] == "NOT MET"
    assert stage["drift"] == 1.3534 and stage["threshold"] == 0.5
    assert stage["hbonds"][0]["share"] == 0.88 and stage["rmsf"][0]["rmsf"] == 0.91
    assert exp["choice"]["structure"] == "8C7Y" and exp["compound_choice"]["compounds"][1]["id"] \
           == "CHEMBL3577124"

    # the page is self-contained: nothing is fetched when it is opened. The page's own code is
    # the template (3Dmol.js does carry a fetch for the remote loading this page never calls),
    # and no element in the file may point at a URL to load.
    from protein_mcp.console_assets import TEMPLATE
    for forbidden in ("fetch(", "XMLHttpRequest", "cdn.", "import("):
        assert forbidden not in TEMPLATE, forbidden
    assert not re.search(r"<(script|link|img|iframe)[^>]*(src|href)\s*=\s*[\"']?(https?:)?//",
                         html), "the page loads something from the network"
    assert "</script" not in json.dumps(data), "a closing tag could end the data block early"
    print("every value comes from the dossier OK")


def test_the_default_dossier_files_are_the_ones_read(folder: Path):
    """get_ligands("P15056") writes ligands.json and get_ligands(..., mutation="V600E") writes
    ligands_v600e.json. A page that only ever looked for the second one showed an empty panel
    for a dossier built with default calls, and silently preferred the filtered query when both
    were there."""
    res, _html, data = run(folder)
    assert data["context"]["files"] == {
        "profile": "profile.json", "structures": "structures_v600e.json",
        "ligands": "ligands_v600e.json", "literature": "literature_v600e-inhibitor.json"}
    assert res["dossier_files_read"]["ligands"] == "ligands_v600e.json"

    # the default call's file is written too: it must win, and the page must say which it read
    default = json.loads(json.dumps(LIGANDS))
    default["result"]["selection"] = "15 of 2104"
    default["result"]["ligands"] = [dict(LIGANDS["result"]["ligands"][0], name="DEFAULT CALL")]
    (folder / "ligands.json").write_bytes(json.dumps(default).encode())
    try:
        res, _html, data = run(folder)
        assert data["context"]["files"]["ligands"] == "ligands.json"
        assert data["context"]["ligands"]["selection"] == "15 of 2104"
        assert data["context"]["ligands"]["rows"][0]["name"] == "DEFAULT CALL"
        assert res["dossier_files_read"]["ligands"] == "ligands.json"
    finally:
        (folder / "ligands.json").unlink()
    # control: with neither file the panel is empty and says nothing it did not read
    moved = folder / "ligands_v600e.json"
    saved = moved.read_bytes()
    moved.unlink()
    try:
        _res, _html, data = run(folder)
        assert data["context"]["files"]["ligands"] is None
        assert data["context"]["ligands"]["rows"] == []
    finally:
        moved.write_bytes(saved)
    print("the default dossier files are the ones read OK")


def test_a_number_nobody_retrieved_stays_empty(folder: Path):
    """The variants table has per-variant columns on purpose, and they must stay null: counting
    PDB entries or papers per variant is a question this project never asked any source."""
    _res, _html, data = run(folder)
    variants = data["context"]["variants"]
    assert [v["change"] for v in variants] == ["V -> E", "G -> R"]
    assert all(v["pdb"] is None and v["chembl"] is None and v["papers"] is None for v in variants)
    # control: the fields exist, so a page that filled them in would differ here
    assert set(variants[0]) >= {"pdb", "chembl", "papers"}
    print("unretrieved counts stay empty OK")


# --------------------------------------------------------------------------- files on disk
def test_a_pose_file_that_changed_is_not_drawn(folder: Path):
    prep = folder / "prep"
    _res, _html, data = run(folder)
    geo = next(r["geometry"] for r in data["experiments"]["runs"] if r["run"] == RUN_T1)
    # the control that must hold first: with the recorded hash the atoms are read, and they are
    # the literal coordinates of the first record minus the box centre
    assert [a["e"] for a in geo["atoms"]] == ["C", "N", "O"], geo
    assert geo["atoms"][0] == {"e": "C", "x": round(10.0 - BOX_CENTRE[0], 3),
                               "y": round(20.0 - BOX_CENTRE[1], 3),
                               "z": round(30.0 - BOX_CENTRE[2], 3)}
    assert geo["bonds"] == [[0, 1], [1, 2]]

    path = prep / "dock" / "8C7Y_TXV_U" / f"poses_{RUN_T1}.sdf"
    saved = path.read_bytes()
    path.write_bytes(saved.replace(b"   10.0000", b"  777.0000", 1))
    try:
        res, html, data = run(folder)
        geo = next(r["geometry"] for r in data["experiments"]["runs"] if r["run"] == RUN_T1)
        assert "atoms" not in geo and "no longer matches the SHA-256" in geo["error"], geo
        # and the 3D panel, where a wrong molecule would look most convincing, refuses it too
        assert [p["run"] for p in data["three_d"]["poses"]] == [RUN_REDOCK]
        assert "777.0000" not in html, "the tampered coordinates reached the page"
        assert res["embedded_3d"]["poses_not_embedded"] == [
            {"run": RUN_T1, "file": f"dock/8C7Y_TXV_U/poses_{RUN_T1}.sdf",
             "reason": "no longer matches the SHA-256 in the manifest"}]
    finally:
        path.write_bytes(saved)
    print("a changed pose file is refused OK")


def test_only_the_first_record_of_a_pose_file_is_embedded(folder: Path):
    """The report is about the top-scoring pose; embedding the whole file would put eight more
    poses in the page with nothing saying which one every figure is about."""
    _res, html, data = run(folder)
    assert "SECOND_RECORD_MARKER" not in html
    poses = {p["run"]: p for p in data["three_d"]["poses"]}
    assert set(poses) == {RUN_T1, RUN_REDOCK}
    assert poses[RUN_T1]["sdf"].count("$$$$") == 1
    assert "pose 1 of this run" in poses[RUN_T1]["sdf"]
    assert poses[RUN_T1]["level"] == "T1" and poses[RUN_T1]["kind"] == "compound"
    assert poses[RUN_REDOCK]["kind"] == "redock"
    # control: the file on disk really does hold the second record
    assert "SECOND_RECORD_MARKER" in POSES_SDF
    print("only the first pose record is embedded OK")


def test_a_series_file_that_changed_is_not_plotted(folder: Path):
    prep = folder / "prep"
    _res, _html, data = run(folder)
    series = data["experiments"]["md"][0]["stages"][0]["series"]
    assert series["backbone PROA"]["v"] == [0.1, 1.35]
    assert series["backbone PROA"]["t"] == [0.0, 10.0], "times are printed in ns"
    assert series["backbone PROA observed"]["drift"] == 0.425
    # the H-bond counts have no time axis of their own: the ligand axis is reused only because
    # the two series have the same length, which is checked
    assert series["H-bonds"]["v"] == [4.0, 2.0] and series["H-bonds"]["t"] == [0.0, 10.0]

    path = prep / "md" / "8C7Y_TXV_U_chembl1229517" / "analysis_production.json"
    saved = path.read_bytes()
    path.write_bytes(saved.replace(b"1.35", b"9.99"))
    try:
        _res, _html, data = run(folder)
        series = data["experiments"]["md"][0]["stages"][0]["series"]
        assert "no longer matches the SHA-256" in series["error"], series
        assert "backbone PROA" not in series
    finally:
        path.write_bytes(saved)
    print("a changed series file is refused OK")


def test_a_series_without_a_matching_axis_is_dropped(folder: Path):
    """A value list and a time list of different lengths are not a series: plotting them would
    pair a frame with the wrong time."""
    prep = folder / "prep"
    path = prep / "md" / "8C7Y_TXV_U_chembl1229517" / "analysis_production.json"
    broken = json.loads(json.dumps(SERIES))
    broken["result"]["ligand_rmsd"]["times_ps"] = [0.0]
    broken["result"]["hbonds"]["count_per_frame"]["values"] = [4, 2, 1]
    text = json.dumps(broken)
    path.write_bytes(text.encode())
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    block = man["structures"][0]["md_analysis"][f"{POCKET} CHEMBL1229517"]["production"]
    saved_sha = block["series_file"]["sha256"]
    block["series_file"]["sha256"] = sha(text.encode())
    (prep / "manifest.json").write_text(json.dumps(man), encoding="utf-8")
    try:
        _res, _html, data = run(folder)
        series = data["experiments"]["md"][0]["stages"][0]["series"]
        assert "ligand RMSD" not in series and "H-bonds" not in series, series
        assert series["COM distance"]["v"] == [3.0, 3.2], "the sound series must survive"
    finally:
        path.write_bytes(json.dumps(SERIES).encode())
        block["series_file"]["sha256"] = saved_sha
        (prep / "manifest.json").write_text(json.dumps(man), encoding="utf-8")
    print("a series without a matching axis is dropped OK")


def test_a_file_outside_prep_is_not_embedded(folder: Path):
    prep = folder / "prep"
    (folder / "outside_receptor.pdb").write_text("ATOM  OUTSIDE_RECEPTOR\nEND\n", encoding="utf-8")
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    saved = json.dumps(man)
    # the SAME bytes as the recorded pose, so its SHA-256 matches the manifest: then only the
    # path guard can refuse it, and the control cannot pass thanks to the hash check instead
    (folder / "outside_poses.sdf").write_bytes(POSES_SDF.encode())
    man["structures"][0]["prepared"]["receptor"]["path"] = "../outside_receptor.pdb"
    runs = man["structures"][0]["docking"][POCKET]["runs"]
    runs[RUN_T1]["files"]["poses_TXV_regenerated_chembl1229517.sdf"]["path"] = "../outside_poses.sdf"
    (prep / "manifest.json").write_text(json.dumps(man), encoding="utf-8")
    try:
        res, html, data = run(folder)
        assert data["three_d"]["receptors"] == {}
        assert "OUTSIDE_RECEPTOR" not in html
        assert res["embedded_3d"]["receptors_not_embedded"] == [
            {"pdb": "8C7Y", "reason": "no cleaned receptor in prep/"}]
        # the same rule for a pose file: not embedded in the 3D panel and not drawn in the
        # compound panel, both of which read a path straight out of the manifest
        assert [p["run"] for p in data["three_d"]["poses"]] == [RUN_REDOCK]
        assert res["embedded_3d"]["poses_not_embedded"] == [
            {"run": RUN_T1, "file": "../outside_poses.sdf", "reason": "missing from prep/"}]
        geo = next(r["geometry"] for r in data["experiments"]["runs"] if r["run"] == RUN_T1)
        assert "atoms" not in geo and "is missing" in geo["error"], geo
    finally:
        (prep / "manifest.json").write_text(saved, encoding="utf-8")
    print("a file outside prep/ is not embedded OK")


def test_the_receptor_cap_is_reported_not_silent(folder: Path):
    prep = folder / "prep"
    man = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    saved = json.dumps(man)
    second = json.loads(json.dumps(man["structures"][0]))
    second["pdb_id"] = "3OG7"
    second["prepared"]["receptor"]["path"] = "clean/3OG7_receptor.pdb"
    (prep / "clean" / "3OG7_receptor.pdb").write_text("ATOM  RECEPTOR_TEXT_3OG7\nEND\n",
                                                      encoding="utf-8")
    man["structures"].append(second)
    (prep / "manifest.json").write_text(json.dumps(man), encoding="utf-8")
    cap = console.MAX_EMBEDDED_RECEPTORS
    console.MAX_EMBEDDED_RECEPTORS = 1
    try:
        res, html, data = run(folder)
        assert sorted(data["three_d"]["receptors"]) == ["8C7Y"]
        assert "RECEPTOR_TEXT_8C7Y" in html and "RECEPTOR_TEXT_3OG7" not in html
        left = res["embedded_3d"]["receptors_not_embedded"]
        assert left and left[0]["pdb"] == "3OG7" and "would be embedded" in left[0]["reason"]
        # the runs of the structure left out are still listed: only its receptor is missing
        assert len([r for r in data["experiments"]["runs"] if r["pdb"] == "3OG7"]) == 2
    finally:
        console.MAX_EMBEDDED_RECEPTORS = cap
        (prep / "manifest.json").write_text(saved, encoding="utf-8")
    print("the receptor cap is reported OK")


# --------------------------------------------------------------------------- lookups
def test_the_compound_of_a_pocket_is_looked_up_not_guessed():
    m = {"ligand_choice": {"compounds": [{"id": "CHEMBL1229517", "name": "VEMURAFENIB"},
                                         {"id": "CHEMBL3577124", "name": "LY-3009120"}]}}
    assert console._compound_of(m, "TXV U CHEMBL3577124")["name"] == "LY-3009120"
    assert console._compound_of(m, "TXV U chembl1229517")["name"] == "VEMURAFENIB"
    # a key that names no recorded compound gets nothing, never the first one in the list
    assert console._compound_of(m, "TXV U") == {}
    assert console._compound_of(m, "TXV U VEMURAFENIB") == {}, "matched by name, not by id"
    assert console._compound_of({}, "TXV U CHEMBL1229517") == {}
    print("the compound of a pocket is looked up OK")


# --------------------------------------------------------------------------- missing inputs
def test_without_a_manifest_the_page_says_so(folder: Path):
    prep = folder / "prep"
    saved = (prep / "manifest.json").read_text(encoding="utf-8")
    (prep / "manifest.json").unlink()
    try:
        res, _html, data = run(folder)
        assert data["experiments"] is None and data["three_d"] is None
        assert res["holds"]["docking_runs"] == 0 and res["holds"]["md_stages_measured"] == 0
        assert res["embedded_3d"] == {"receptors": [], "poses": 0,
                                      "receptors_not_embedded": [], "poses_not_embedded": []}
        # the phase 1 panels are still there: a missing manifest is not a missing dossier
        assert data["context"]["identity"]["gene"] == "BRAF"
    finally:
        (prep / "manifest.json").write_text(saved, encoding="utf-8")

    (prep / "manifest.json").write_text("{not json", encoding="utf-8")
    try:
        bad = console.build_console_page(ACC)
        assert not bad["ok"] and "could not be read" in bad["error"], bad
    finally:
        (prep / "manifest.json").write_text(saved, encoding="utf-8")
    print("a missing manifest is reported OK")


def test_without_a_dossier_nothing_is_written(folder: Path):
    res = console.build_console_page("Q99999")
    assert not res["ok"] and "No dossier for Q99999" in res["error"]
    assert not (dossier.DOSSIER_DIR / "Q99999").exists()
    bad = console.build_console_page("not-an-accession")
    assert not bad["ok"], bad
    print("a missing dossier is refused OK")


def test_the_documents_it_links_are_the_ones_that_exist(folder: Path):
    res, html, _data = run(folder)
    assert res["documents_linked"] == [] and res["documents_missing"] == [
        f"report_{ACC}.pdf", f"experiments_{ACC}.pdf"]
    (folder / f"experiments_{ACC}.pdf").write_bytes(b"%PDF-1.4 stub")
    try:
        res, html, data = run(folder)
        assert res["documents_linked"] == [f"experiments_{ACC}.pdf"]
        assert res["documents_missing"] == [f"report_{ACC}.pdf"]
        assert [d["exists"] for d in data["documents"]] == [False, True]
    finally:
        (folder / f"experiments_{ACC}.pdf").unlink()
    print("the linked documents are the ones that exist OK")


def test_the_viewer_library_is_embedded_or_named(folder: Path):
    res, html, data = run(folder)
    assert "3Dmol" in res["viewer"] and "KB)" in res["viewer"], res["viewer"]
    assert console.VIEWER_ASSET.is_file(), "the 3Dmol asset is missing from the package"
    assert "$3Dmol" in html or "3Dmol" in html
    # the page redistributes the bundle, so the BSD-3-Clause notice has to be in the page and
    # the licence files have to be in the package
    assert "Copyright (c) 2014, University of Pittsburgh" in html
    assert "BSD-3-Clause" in html and console.VIEWER_LICENSE in res["viewer"]
    assert data["viewer_credit"] == f"3D viewer: {console.VIEWER_NAME}, {console.VIEWER_LICENSE}"
    assets = console.VIEWER_ASSET.parent
    for name in ("3Dmol-LICENSE.txt", "3Dmol-min.js.LICENSE.txt", "PROVENANCE.md"):
        assert (assets / name).is_file(), name
    # the file the page carries is the build PROVENANCE.md names, by hash and not by name
    recorded = (assets / "PROVENANCE.md").read_text(encoding="utf-8")
    assert hashlib.sha256(console.VIEWER_ASSET.read_bytes()).hexdigest() in recorded
    assert console.VIEWER_NAME.split()[-1] in recorded
    saved = console.VIEWER_ASSET
    console.VIEWER_ASSET = folder / "no_such_library.js"
    try:
        res, html, _data = run(folder)
        assert "not embedded" in res["viewer"], res["viewer"]
        # and it does not fall back to a CDN
        assert "src=\"http" not in html
    finally:
        console.VIEWER_ASSET = saved
    print("the viewer library is embedded or named OK")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        folder = build_dossier(Path(tmp))
        test_every_value_in_the_page_comes_from_the_dossier(folder)
        test_the_default_dossier_files_are_the_ones_read(folder)
        test_a_number_nobody_retrieved_stays_empty(folder)
        test_a_pose_file_that_changed_is_not_drawn(folder)
        test_only_the_first_record_of_a_pose_file_is_embedded(folder)
        test_a_series_file_that_changed_is_not_plotted(folder)
        test_a_series_without_a_matching_axis_is_dropped(folder)
        test_a_file_outside_prep_is_not_embedded(folder)
        test_the_receptor_cap_is_reported_not_silent(folder)
        test_the_compound_of_a_pocket_is_looked_up_not_guessed()
        test_without_a_manifest_the_page_says_so(folder)
        test_without_a_dossier_nothing_is_written(folder)
        test_the_documents_it_links_are_the_ones_that_exist(folder)
        test_the_viewer_library_is_embedded_or_named(folder)
    print("ALL OK")


main()
