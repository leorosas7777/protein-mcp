"""Offline test for analyze_admet (protein_mcp/files/admet.py + the ADME part of
protein_mcp/sources/chembl.py, Phase 3 step 7).

No network. The ChEMBL pages are fixtures, and every shape and every defect in them was COPIED
FROM THE REAL API on 2026-11-11 (ESTADO 8.5), because those are what a naive version gets wrong:

- `assay_type=A` is NOT ADMET. Vemurafenib's 17 "A" rows are 9 `Activity` (HepG2 cytotoxicity,
  % viability), 3 `IC50` (cytotoxicity), 2 `Ratio` and only 3 real ADME rows. The fixture keeps
  that mix, so a filter that trusts ChEMBL's own label reports cytotoxicity as ADME;
- the only aqueous solubility is filed under `assay_type=P`, not A;
- the vocabulary is case-inconsistent: `Papp` and `permeability` in the same molecule;
- units are incompatible inside one standard_type: the fixture's two `CL` rows are
  `uL.min-1.(10^6cells)-1` and `ml kg-1`, so anything that averages them is wrong by construction;
- values arrive CENSORED (`standard_relation` '>'), and absent (`standard_value` None with
  `activity_comment` 'Not Determined');
- `max_phase` is the STRING '4.0' and every `molecule_properties` number is a STRING ('5.54');
- `pref_name` is None for a research compound.

Controls that must fail (each asserts the WRONG answer explicitly, so the right code cannot be
replaced by the naive one without a red test):

- `test_without_the_allowlist_cytotoxicity_would_be_reported_as_adme`
- `test_querying_only_assay_type_a_would_lose_the_solubility`
- `test_averaging_the_two_clearance_rows_would_invent_a_number`
- `test_a_connectivity_block_key_is_refused`
- `test_a_mismatched_pose_and_crystal_never_reach_chembl`

ADDED after the external review of 2026-11-11, which found that 30 passing tests touched neither
`analyze_admet` nor the report section, so several mutations survived them. The end-to-end group
below builds a real manifest in a temporary dossier, runs the TOOL with `chembl._get` stubbed,
and renders section 3 of the experiment report from what it recorded:

- `test_the_tool_runs_end_to_end_and_records_the_block`
- `test_the_report_section_three_is_rendered_from_the_block`
- `test_the_report_never_calls_an_excluded_row_not_adme`
- `test_a_truncated_lookup_says_so_in_the_selection_and_the_report`
- `test_the_links_are_built_from_chembl_ids`
- `test_the_applied_cutoffs_and_the_published_cutoffs_cannot_drift`
"""
import asyncio
import json
import os
import re
import sys
import tempfile
from pathlib import Path

import httpx

from protein_mcp.files import admet
from protein_mcp.sources import chembl

# ---------------------------------------------------------------- fixtures (real shapes)
VEM_KEY = "GPXBXXGIAQBQNI-UHFFFAOYSA-N"
VEM_SMILES = "CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(-c4ccc(Cl)cc4)cc23)c1F"
TXV_SMILES = "CC(C)(C#N)c1cccc(C(=O)Nc2ccc(F)c(C(=O)c3c[nH]c4ncc(Cl)cc34)c2F)c1"
# a sentence that mentions prediction is acceptable only if it DENIES it
_NEGATED = re.compile(r"\b(not|no|never|nothing|neither|none|reserved)\b", re.I)

MOLECULE_PAGE = {
    "page_meta": {"total_count": 1},
    "molecules": [{
        "molecule_chembl_id": "CHEMBL1229517",
        "pref_name": "VEMURAFENIB",
        "max_phase": "4.0",                      # a STRING in the real API, not a number
        "first_approval": 2011,
        "withdrawn_flag": False, "black_box_warning": 0,
        "oral": True, "parenteral": False, "topical": False,
        "molecule_type": "Small molecule", "natural_product": 0, "prodrug": 0,
        "therapeutic_flag": True, "availability_type": 1, "chirality": 2,
        "dosed_ingredient": True, "indication_class": None, "usan_year": 2010,
        "usan_stem": "-rafenib", "first_in_class": 1, "structure_type": "MOL",
        "inorganic_flag": 0, "polymer_flag": 0, "orphan": 0, "chemical_probe": 0,
        "veterinary": 0,
        "atc_classifications": ["L01EC01"],
        "cross_references": [{"xref_id": "vemurafenib", "xref_name": "vemurafenib",
                              "xref_src": "DailyMed"}],
        "molecule_structures": {"standard_inchi_key": VEM_KEY, "canonical_smiles": VEM_SMILES},
        # the real API carries these; the fixture used to omit molecule_synonyms, which is
        # EXACTLY the field a name-inventing mutation would reach for, so the test that claims
        # to cover "no invented name" could not have failed (reviewer finding, 2026-11-11)
        "molecule_synonyms": [
            {"molecule_synonym": "PLX-4032", "syn_type": "RESEARCH_CODE"},
            {"molecule_synonym": "ZELBORAF", "syn_type": "TRADE_NAME"},
            {"molecule_synonym": "RO 5185426", "syn_type": "RESEARCH_CODE"},
        ],
        "usan_stem_definition": "protein kinase inhibitors", "usan_substem": "-rafenib",
        "biotherapeutic": None, "helm_notation": None,
        "molecule_properties": {                 # every number a STRING, as ChEMBL returns them
            "alogp": "5.54", "aromatic_rings": 4, "full_mwt": "489.93", "hba": 4, "hbd": 2,
            "heavy_atoms": 33, "num_ro5_violations": 1, "psa": "91.92",
            "qed_weighted": "0.33", "ro3_pass": "N", "rtb": 7,
        },
        "molecule_hierarchy": {"molecule_chembl_id": "CHEMBL1229517",
                               "parent_chembl_id": "CHEMBL1229517"},
    }],
}

ADME_ROWS_A = [
    # the three real ADME rows, with the case inconsistency and the incompatible units
    {"standard_type": "CL", "standard_relation": None, "standard_value": "6.61",
     "standard_units": "uL.min-1.(10^6cells)-1", "assay_chembl_id": "CHEMBL3301372",
     "assay_description": "ASTRAZENECA: Intrinsic clearance measured in human hepatocytes",
     "document_chembl_id": "CHEMBL3301361", "bao_label": "cell-based format"},
    {"standard_type": "CL", "standard_relation": None, "standard_value": "12.0",
     "standard_units": "ml kg-1", "assay_chembl_id": "CHEMBL9999001",
     "assay_description": "Clearance in rat after intravenous dose",
     "document_chembl_id": "CHEMBL9999002", "target_organism": "Rattus norvegicus"},
    {"standard_type": "permeability", "standard_relation": "=", "standard_value": "2.9",
     "standard_units": "10^-6 cm/s", "assay_chembl_id": "CHEMBL3301373",
     "assay_description": "Permeability in human Caco-2 cells",
     "document_chembl_id": "CHEMBL3301361"},
    {"standard_type": "Papp", "standard_relation": "=", "standard_value": "3.0",
     "standard_units": "10^-6 cm/s", "assay_chembl_id": "CHEMBL3301374",
     "assay_description": "Apparent permeability apical to basolateral in LLC-PK1 cells",
     "document_chembl_id": "CHEMBL3301361"},
    # a censored row that must keep its relation
    {"standard_type": "T1/2", "standard_relation": ">", "standard_value": "120.0",
     "standard_units": "hr", "assay_chembl_id": "CHEMBL9999010",
     "assay_description": "Half-life in human liver microsomes",
     "document_chembl_id": "CHEMBL9999011"},
]
NON_ADME_ROWS_A = (
    [{"standard_type": "Activity", "standard_relation": ">", "standard_value": "75.0",
      "standard_units": "%", "assay_chembl_id": f"CHEMBL336551{i}",
      "assay_description": "Cytotoxicity against human HepG2 cells assessed as cell viability",
      "document_chembl_id": "CHEMBL3352089"} for i in range(9)]
    + [{"standard_type": "IC50", "standard_relation": None, "standard_value": None,
        "standard_units": None, "assay_chembl_id": f"CHEMBL434724{i}",
        "activity_comment": "Not Determined",          # absent value, real shape
        "assay_description": "Cytotoxicity against human HELF cells",
        "document_chembl_id": "CHEMBL4346662"} for i in range(3)]
    + [{"standard_type": "Ratio", "standard_relation": None, "standard_value": "2.0",
        "standard_units": None, "assay_chembl_id": f"CHEMBL999902{i}",
        "assay_description": "Ratio of something that is not ADME",
        "document_chembl_id": "CHEMBL9999021"} for i in range(2)]
)
ADME_ROWS_P = [
    {"standard_type": "Solubility", "standard_relation": "=", "standard_value": "2.0",
     "standard_units": "ug.mL-1", "assay_chembl_id": "CHEMBL9999030",
     "assay_description": "Aqueous solubility of the compound in FaSSIF at pH 6.5",
     "document_chembl_id": "CHEMBL9999031"},
]
ROWS = {"A": ADME_ROWS_A + NON_ADME_ROWS_A, "P": ADME_ROWS_P}
assert len(ROWS["A"]) == 19 and len(ROWS["P"]) == 1


def fake_get(pages=None, molecule=None):
    """Stand-in for chembl._get that pages exactly as the real endpoint does."""
    calls = []

    async def _get(client, endpoint, params, deadline=None):
        # the real _get takes the shared budget as its fourth argument (F1); a fixture that did
        # not accept it would hide every call that passes one
        calls.append((endpoint, dict(params), deadline))
        if endpoint == "molecule":
            return molecule if molecule is not None else MOLECULE_PAGE
        rows = (pages if pages is not None else ROWS).get(params.get("assay_type"), [])
        off, lim = int(params.get("offset", 0)), int(params.get("limit", 20))
        return {"page_meta": {"total_count": len(rows)}, "activities": rows[off:off + lim]}

    _get.calls = calls
    return _get


def run(coro):
    return asyncio.run(coro)


def adme(**kw):
    original = chembl._get
    chembl._get = fake_get(**kw)
    try:
        return run(chembl.fetch_adme_activities("CHEMBL1229517")), chembl._get.calls
    finally:
        chembl._get = original


def mol_from_smiles(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=42)
    return m


def write_sdf(path: Path, smiles):
    from rdkit import Chem
    w = Chem.SDWriter(str(path))
    w.write(mol_from_smiles(smiles))
    w.close()


# ---------------------------------------------------------------- the fixture is what it claims
def test_the_fixture_copies_the_real_defects():
    types = [r["standard_type"] for r in ROWS["A"]]
    assert types.count("Activity") == 9 and types.count("IC50") == 3
    assert "Papp" in types and "permeability" in types, "case inconsistency must be in the fixture"
    cl_units = {r["standard_units"] for r in ROWS["A"] if r["standard_type"] == "CL"}
    assert len(cl_units) == 2, "the two CL rows must have incompatible units"
    assert ROWS["P"][0]["standard_type"] == "Solubility", "solubility must sit under P, not A"
    assert MOLECULE_PAGE["molecules"][0]["max_phase"] == "4.0", "max_phase is a string"
    assert all(isinstance(MOLECULE_PAGE["molecules"][0]["molecule_properties"][k], str)
               for k in ("alogp", "psa", "full_mwt"))


# ---------------------------------------------------------------- the allowlist
def test_only_the_allowlisted_types_survive():
    res, _ = adme()
    kept = sorted(r["standard_type"] for r in res["adme_rows"])
    assert kept == ["CL", "CL", "Papp", "Solubility", "T1/2", "permeability"], kept
    assert res["rows_excluded_by_allowlist_by_standard_type"] == {"Activity": 9, "IC50": 3, "Ratio": 2}
    assert res["rows_read"] == 20


def test_without_the_allowlist_cytotoxicity_would_be_reported_as_adme():
    """CONTROL THAT MUST FAIL: trusting ChEMBL's assay_type=A label reports HepG2 cytotoxicity
    and 'Ratio' as ADME. The real code must not agree with this count."""
    naive = [r for r in ROWS["A"] + ROWS["P"]]            # what ChEMBL labels A and P, unfiltered
    assert len(naive) == 20
    assert any(r["standard_type"] == "Activity" for r in naive)
    res, _ = adme()
    assert len(res["adme_rows"]) == 6 < len(naive)
    assert not any(r["standard_type"] in ("Activity", "Ratio") for r in res["adme_rows"])


def test_a_type_that_is_not_adme_has_no_class():
    for t in ("Hammett constant", "Saluretic potency", "Tail vein toxicity", "Activity",
              "No. of animals protected/animals tested", "", None):
        assert chembl.adme_class(t) is None, t
    assert chembl.adme_class("PAPP") == chembl.adme_class("Papp") == "permeability"
    assert chembl.adme_class("  CL  ") == "clearance", "whitespace must not change the class"


# ---------------------------------------------------------------- assay types queried
def test_querying_only_assay_type_a_would_lose_the_solubility():
    """CONTROL THAT MUST FAIL: with A alone the only aqueous solubility disappears, which is the
    bug measured on the real API on 2026-11-11."""
    a_only, _ = adme(pages={"A": ROWS["A"], "P": []})
    assert not any(r["adme_class"] == "solubility" for r in a_only["adme_rows"])
    both, calls = adme()
    assert any(r["adme_class"] == "solubility" for r in both["adme_rows"])
    assert {c[1]["assay_type"] for c in calls if c[0] == "activity"} == {"A", "P"}


def test_what_was_not_queried_is_named():
    res, _ = adme()
    assert set(res["assay_types_queried"]) == {"A", "P"}
    assert set(res["assay_types_not_queried"]) == {"B", "F", "T", "U"}
    assert "not queried" in res["note"]
    assert res["rows_total_by_assay_type"] == {"A": 19, "P": 1}


def test_the_page_size_is_small_enough_for_the_api():
    """limit=1000 and limit=200 returned HTTP 500 three times out of three on the real API."""
    assert chembl.ADME_PAGE <= 50, "a bigger page makes ChEMBL return 500 on this query"
    res, calls = adme()
    assert all(c[1]["limit"] == chembl.ADME_PAGE for c in calls if c[0] == "activity")
    assert len(res["adme_rows"]) == 6, "paging must not drop or duplicate rows"


def test_paging_reads_every_row():
    many = {"A": [dict(ADME_ROWS_A[0], assay_chembl_id=f"CHEMBL{i}") for i in range(47)], "P": []}
    res, calls = adme(pages=many)
    assert res["rows_read"] == 47 and len(res["adme_rows"]) == 47
    assert len([c for c in calls if c[0] == "activity"]) == 4, "3 pages of A + 1 of P"


# ---------------------------------------------------------------- units and censoring
def test_every_row_keeps_its_own_units_and_relation():
    res, _ = adme()
    cl = [r for r in res["adme_rows"] if r["standard_type"] == "CL"]
    assert {r["units"] for r in cl} == {"uL.min-1.(10^6cells)-1", "ml kg-1"}
    half = next(r for r in res["adme_rows"] if r["standard_type"] == "T1/2")
    assert half["relation"] == ">" and half["value"] == "120.0"
    for r in res["adme_rows"]:
        assert "assay_description" in r and r["assay_description"]
        assert r["assay_type"] in ("A", "P") and r["assay_type_meaning"]


def test_averaging_the_two_clearance_rows_would_invent_a_number():
    """CONTROL THAT MUST FAIL: the mean of 6.61 uL/min/10^6cells and 12.0 ml/kg is 9.3 of
    nothing. No aggregate of any kind may appear in the output."""
    res, _ = adme()
    assert (6.61 + 12.0) / 2 == 9.305                      # the number that must never appear
    text = json.dumps(res)
    assert "9.3" not in text, "an average of incompatible units appeared in the output"
    # no field may BE an aggregate ('mean residence time' is an ADME class, not an aggregate,
    # and 'assay_type_meaning' is a label, so the check is on the key names, not the text)
    def keys(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                yield k
                yield from keys(v)
        elif isinstance(obj, list):
            for v in obj:
                yield from keys(v)
    banned = {"mean", "average", "median", "aggregate", "aggregated", "pooled", "summary",
              "mean_value", "value_mean", "best", "worst", "n_measurements"}
    for key in keys(res):
        assert key.lower() not in banned, key
        assert not key.lower().startswith(("mean_", "avg_", "median_")), key
    # every value stays exactly the string ChEMBL returned: no arithmetic happened to it
    assert all(isinstance(r["value"], (str, type(None))) for r in res["adme_rows"])
    assert {r["value"] for r in res["adme_rows"] if r["standard_type"] == "CL"} == {"6.61", "12.0"}


# ---------------------------------------------------------------- identity
def test_a_connectivity_block_key_is_refused():
    """CONTROL THAT MUST FAIL: the 14-character connectivity block also matches another salt,
    tautomer or stereoisomer, so it must never be sent as an exact key."""
    for bad in (VEM_KEY[:14], "", None, "not-a-key", VEM_KEY.lower() + "X",
                "GPXBXXGIAQBQNI-UHFFFAOYSA", "GPXBXXGIAQBQNI-UHFFFAOYSA-1"):
        res = run(chembl.fetch_molecule_by_inchikey(bad))
        assert res["found"] is False and "error" in res, bad
        assert "not queried" in res["error"], "a malformed key must not reach the network"


def test_a_lowercase_key_is_accepted_and_normalised():
    original = chembl._get
    chembl._get = fake_get()
    try:
        res = run(chembl.fetch_molecule_by_inchikey(VEM_KEY.lower()))
        assert res["found"] and res["inchi_key"] == VEM_KEY
        assert chembl._get.calls[0][1]["molecule_structures__standard_inchi_key"] == VEM_KEY
    finally:
        chembl._get = original


def test_two_molecules_with_the_same_key_choose_neither():
    page = {"page_meta": {"total_count": 2},
            "molecules": [{"molecule_chembl_id": "CHEMBL1"}, {"molecule_chembl_id": "CHEMBL2"}]}
    original = chembl._get
    chembl._get = fake_get(molecule=page)
    try:
        res = run(chembl.fetch_molecule_by_inchikey(VEM_KEY))
        assert res["found"] is False and res["matches"] == 2
        assert res["candidates"] == ["CHEMBL1", "CHEMBL2"]
    finally:
        chembl._get = original


def test_an_unnamed_research_compound_is_not_given_a_name():
    """The fixture DOES carry molecule_synonyms (PLX-4032, ZELBORAF), so a fallback from
    pref_name to a synonym is the mutation this test has to catch, and it can."""
    page = json.loads(json.dumps(MOLECULE_PAGE))
    page["molecules"][0]["pref_name"] = None
    page["molecules"][0]["max_phase"] = None
    assert page["molecules"][0]["molecule_synonyms"], "the name-inventing path must be reachable"
    original = chembl._get
    chembl._get = fake_get(molecule=page)
    try:
        res = run(chembl.fetch_molecule_by_inchikey(VEM_KEY))
        assert res["found"] and res["pref_name"] is None
        text = json.dumps(res)
        for name in ("VEMURAFENIB", "PLX-4032", "ZELBORAF", "RO 5185426"):
            assert name not in text, f"{name} was carried into the result as a name"
    finally:
        chembl._get = original


# ---------------------------------------------------------------- the local calculation
def test_a_minus_one_flag_is_named_and_not_interpreted():
    """MEASURED 2026-11-11: TXV comes back with prodrug, chirality, availability_type,
    first_in_class and inorganic_flag all -1. Reporting -1 as a number reads as a fact, and
    reporting it as 'no' would be an interpretation this project cannot support."""
    page = json.loads(json.dumps(MOLECULE_PAGE))
    mol = page["molecules"][0]
    for f in ("prodrug", "chirality", "availability_type", "first_in_class", "inorganic_flag"):
        mol[f] = -1
    original = chembl._get
    chembl._get = fake_get(molecule=page)
    try:
        res = run(chembl.fetch_molecule_by_inchikey(VEM_KEY))
        assert res["record_fields_minus_one"] == [
            "availability_type", "chirality", "first_in_class", "inorganic_flag", "prodrug"]
        assert res["record"]["prodrug"] == -1, "the raw value is kept, not rewritten"
    finally:
        chembl._get = original
    # a molecule with real 0/1 values must report an EMPTY list, or the marker means nothing
    chembl._get = fake_get()
    try:
        assert run(chembl.fetch_molecule_by_inchikey(VEM_KEY))["record_fields_minus_one"] == []
    finally:
        chembl._get = original


def test_the_descriptors_name_the_function_that_produced_them():
    desc = admet.descriptors(mol_from_smiles(VEM_SMILES))
    assert desc["heavy_atoms"] == 33 and desc["molecular_formula"] == "C23H18ClF2N3O3S"
    assert round(desc["molecular_weight"], 1) == 489.9
    assert 5.0 < desc["clogp"] < 6.0 and 85 < desc["tpsa"] < 95
    assert desc["engine"].startswith("RDKit")
    for field in ("clogp", "tpsa", "hbd", "hba", "rotatable_bonds"):
        assert desc["functions"][field].startswith("rdkit"), field


def test_the_two_engines_are_reported_apart_and_never_merged():
    desc = admet.descriptors(mol_from_smiles(VEM_SMILES))
    block = {"rdkit": desc, "chembl": MOLECULE_PAGE["molecules"][0]["molecule_properties"],
             "note": admet.CALCULATED_NOTE}
    assert block["rdkit"]["clogp"] != float(block["chembl"]["alogp"]), \
        "the two engines differ; a fixture where they matched would hide a merge"
    assert "never averaged" in admet.CALCULATED_NOTE
    assert "confirm" in admet.CALCULATED_NOTE


# ---------------------------------------------------------------- the rules
def test_vemurafenib_breaks_lipinski_on_logp_only():
    reading = admet.rule_reading(admet.descriptors(mol_from_smiles(VEM_SMILES)))
    lip = reading["results"]["Lipinski"]
    assert lip["violations"] == 1, lip
    assert len(lip["violated"]) == 1 and lip["violated"][0].startswith("clogp "), lip["violated"]
    assert "> 5" in lip["violated"][0], "the violated line must name the cut-off it exceeded"
    assert reading["results"]["Veber"]["violations"] == 0
    assert reading["results"]["Egan"]["violations"] == 0


def test_the_rules_are_offered_not_measured():
    reading = admet.rule_reading(admet.descriptors(mol_from_smiles(TXV_SMILES)))
    assert reading["status"] == "READING OFFERED, NOT MEASURED"
    assert set(reading["rules_used"]) == {"Lipinski", "Veber", "Egan"}
    for name, rule in reading["rules_used"].items():
        assert rule["cutoffs"], name
    text = json.dumps(reading).lower()
    for forbidden in ("verdict", "passes", "drug-like", "druglike", "good", "bad", "safe"):
        assert forbidden not in text, forbidden
    low = reading["disclaimer"].lower()
    assert "not measurements" in low and "not predictions" in low
    assert "not citations" in low, "the rule names must be declared as labels, not citations"
    assert "conventions" in low and "judgement decides" in low


def test_a_missing_descriptor_makes_the_count_unknown_not_zero():
    reading = admet.rule_reading({"molecular_weight": 400.0, "clogp": None, "hbd": 1, "hba": 2,
                                  "rotatable_bonds": 3, "tpsa": 50.0})
    assert reading["results"]["Lipinski"]["violations"] is None
    assert reading["results"]["Lipinski"]["not_evaluated"] == ["clogp"]
    assert reading["results"]["Veber"]["violations"] == 0, "Veber does not need clogp"
    assert reading["results"]["Egan"]["violations"] is None


def test_a_molecule_over_every_cutoff_counts_every_violation():
    reading = admet.rule_reading({"molecular_weight": 900.0, "clogp": 9.0, "hbd": 9, "hba": 20,
                                  "rotatable_bonds": 25, "tpsa": 300.0})
    assert reading["results"]["Lipinski"]["violations"] == 4
    assert reading["results"]["Veber"]["violations"] == 2
    assert reading["results"]["Egan"]["violations"] == 2


def test_the_cutoffs_are_the_pre_registered_ones():
    """ESTADO 7.3, written before this code. A silent change here is a changed criterion."""
    assert admet.RULES["Lipinski"]["cutoffs"] == ["molecular_weight <= 500", "clogp <= 5",
                                                  "hbd <= 5", "hba <= 10"]
    assert admet.RULES["Veber"]["cutoffs"] == ["rotatable_bonds <= 10", "tpsa <= 140"]
    assert admet.RULES["Egan"]["cutoffs"] == ["tpsa <= 131.6", "clogp <= 5.88"]


# ---------------------------------------------------------------- the molecule that is described
def test_the_pose_and_the_crystal_sdf_give_the_same_key():
    with tempfile.TemporaryDirectory() as tmp:
        a, b = Path(tmp) / "pose.sdf", Path(tmp) / "crystal.sdf"
        write_sdf(a, VEM_SMILES)
        write_sdf(b, VEM_SMILES)
        assert admet._inchikey(admet._first_mol(a)) == admet._inchikey(admet._first_mol(b))
        assert admet._inchikey(admet._first_mol(a)) == VEM_KEY


def test_a_mismatched_pose_and_crystal_never_reach_chembl():
    """CONTROL THAT MUST FAIL: if the two files are different compounds the run and the pocket
    have drifted apart, and attributing ChEMBL data to either would be a fabricated link."""
    with tempfile.TemporaryDirectory() as tmp:
        prep = Path(tmp)
        (prep / "ligands").mkdir()
        (prep / "dock").mkdir()
        pose, crystal = prep / "dock" / "poses_X_regenerated.sdf", prep / "ligands" / "X_C.sdf"
        write_sdf(pose, VEM_SMILES)
        write_sdf(crystal, TXV_SMILES)       # a DIFFERENT compound on purpose
        assert admet._inchikey(admet._first_mol(pose)) != \
            admet._inchikey(admet._first_mol(crystal))
        s = {"pdb_id": "XXXX",
             "docking": {"X C": {"runs": {"X_regenerated": {
                 "files": {"poses_X_regenerated.sdf":
                           {"path": "dock/poses_X_regenerated.sdf"}}}}}},
             "nonpolymer_instances": [{"comp_id": "X", "label_asym_id": "C",
                                       "sdf": {"ok": True, "path": "ligands/X_C.sdf"}}],
             "prepared": {"pockets": [{"comp_id": "X", "label_asym_id": "C",
                                       "heavy_atoms": 33}]}}
        snap = {"structure": "XXXX", "pocket": {"comp_id": "X", "label_asym_id": "C"}}
        try:
            admet.identity(prep, s, snap, "X C", "", 1)
            raise AssertionError("a mismatch must be refused, not reported")
        except admet.AdmetError as exc:
            assert "different InChIKeys" in str(exc)
            assert exc.extra["compound_identity"]["identity_check"] == "MISMATCH"


def test_without_a_docking_run_the_crystal_sdf_is_used_and_marked():
    with tempfile.TemporaryDirectory() as tmp:
        prep = Path(tmp)
        (prep / "ligands").mkdir()
        write_sdf(prep / "ligands" / "X_C.sdf", VEM_SMILES)
        s = {"pdb_id": "XXXX",
             "nonpolymer_instances": [{"comp_id": "X", "label_asym_id": "C",
                                       "sdf": {"ok": True, "path": "ligands/X_C.sdf"}}],
             "prepared": {"pockets": []}}
        snap = {"structure": "XXXX", "pocket": {"comp_id": "X", "label_asym_id": "C"}}
        ident, mol = admet.identity(prep, s, snap, "X C", "", 1)
        assert ident["source"].startswith("crystallographic")
        assert "no docking run" in ident["source"]
        assert ident["run_id"] is None and ident["pose"] is None
        assert ident["inchikey"] == VEM_KEY and ident["source_sha256"]
        assert ident["identity_check"] == "only one SDF was available, so there was nothing to compare"


def test_a_pose_file_that_changed_since_the_manifest_is_refused():
    with tempfile.TemporaryDirectory() as tmp:
        prep = Path(tmp)
        (prep / "dock").mkdir()
        write_sdf(prep / "dock" / "p.sdf", VEM_SMILES)
        s = {"pdb_id": "XXXX",
             "docking": {"X C": {"runs": {"X_regenerated": {"files": {
                 "p.sdf": {"path": "dock/p.sdf", "sha256": "0" * 64}}}}}},
             "nonpolymer_instances": [], "prepared": {"pockets": []}}
        snap = {"structure": "XXXX", "pocket": {"comp_id": "X", "label_asym_id": "C"}}
        try:
            admet.identity(prep, s, snap, "X C", "", 1)
            raise AssertionError("a changed SDF must be refused")
        except admet.AdmetError as exc:
            assert "SHA-256" in str(exc)


def test_a_file_outside_prep_is_refused():
    with tempfile.TemporaryDirectory() as tmp:
        prep = Path(tmp) / "prep"
        prep.mkdir()
        outside = Path(tmp) / "elsewhere.sdf"
        write_sdf(outside, VEM_SMILES)
        s = {"pdb_id": "XXXX", "nonpolymer_instances": [
            {"comp_id": "X", "label_asym_id": "C",
             "sdf": {"ok": True, "path": "../elsewhere.sdf"}}], "prepared": {"pockets": []}}
        snap = {"structure": "XXXX", "pocket": {"comp_id": "X", "label_asym_id": "C"}}
        try:
            admet.identity(prep, s, snap, "X C", "", 1)
            raise AssertionError("a path outside prep/ must be refused")
        except admet.AdmetError as exc:
            assert "outside prep" in str(exc)


# ---------------------------------------------------------------- what the output must not say
def test_nothing_in_the_output_is_called_predicted():
    res, _ = adme()
    desc = admet.descriptors(mol_from_smiles(VEM_SMILES))
    reading = admet.rule_reading(desc)
    blob = {"measured": res, "calculated": desc, "reading": reading}

    def walk(obj, path=""):
        if isinstance(obj, dict):
            for k, v in obj.items():
                assert "predict" not in k.lower(), f"{path}.{k}"
                yield from walk(v, f"{path}.{k}")
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                yield from walk(v, f"{path}[{i}]")
        elif isinstance(obj, str):
            yield path, obj

    # 'predict' may appear ONLY inside a denial; a bare claim of prediction is the defect
    for path, value in walk(blob):
        for sentence in value.replace(";", ".").split("."):
            if "predict" in sentence.lower():
                assert "not " in sentence.lower() or "no " in sentence.lower(), f"{path}: {sentence}"
    assert "PREDICTED" not in res.get("status", "") + reading["status"]
    assert "not predictions" in reading["disclaimer"].lower()


def test_a_failed_lookup_is_unknown_not_no_data():
    async def boom(client, endpoint, params, deadline=None):
        raise RuntimeError("ChEMBL activity failed: HTTP 500")
    original = chembl._get
    chembl._get = boom
    try:
        try:
            run(chembl.fetch_adme_activities("CHEMBL1229517"))
            raise AssertionError("the failure must propagate, not return an empty list")
        except RuntimeError as exc:
            assert "500" in str(exc)
    finally:
        chembl._get = original


def test_a_molecule_with_no_adme_rows_says_how_many_were_read():
    res, _ = adme(pages={"A": NON_ADME_ROWS_A, "P": []})
    assert res["adme_rows"] == []
    assert res["rows_read"] == 14
    assert sum(res["rows_excluded_by_allowlist_by_standard_type"].values()) == 14
    assert "14" in res["note"]


# ------------------------------------------------- end to end: the TOOL and the REPORT
# The review of 2026-11-11 found that nothing called analyze_admet, compact, _selection or
# _record, and that the report section was never rendered with content. These build a real
# manifest in a temporary dossier so the tool and the document are exercised, not imitated.
def temp_dossier(tmp: Path, *, runs=True):
    """A minimal but REAL manifest: the shape dock_ligand and prepare_structures write."""
    import importlib

    from protein_mcp import dossier as dossier_mod
    os.environ["PROTEIN_MCP_DOSSIER_DIR"] = str(tmp / "dossiers")
    importlib.reload(dossier_mod)
    from protein_mcp.files import choice as choice_mod
    importlib.reload(choice_mod)

    prep = tmp / "dossiers" / "P15056" / "prep"
    (prep / "ligands").mkdir(parents=True)
    (prep / "dock" / "3OG7_032_C").mkdir(parents=True)
    pose = prep / "dock" / "3OG7_032_C" / "poses_032_regenerated.sdf"
    crystal = prep / "ligands" / "3OG7_032_C.sdf"
    write_sdf(pose, VEM_SMILES)
    write_sdf(crystal, VEM_SMILES)
    (prep / "clean").mkdir()
    (prep / "clean" / "3OG7_receptor.pdb").write_text("REMARK fixture\nEND\n", encoding="utf-8")

    structure = {
        "pdb_id": "3OG7",
        "resolution_angstrom": 2.45,
        "file": {"path": "raw/3OG7.cif", "sha256": "f" * 64},
        "construct": {"mutations": [], "mutations_from_title": ["V600E"],
                      "mutations_unparsed": []},
        "nonpolymer_instances": [{"comp_id": "032", "label_asym_id": "C", "category": "ligand",
                                  "auth_asym_id": "A", "name": "Vemurafenib",
                                  "synonyms": ["PLX4032"], "formula_weight": 489.92,
                                  "sdf": {"ok": True, "path": "ligands/3OG7_032_C.sdf"}}],
        "prepared": {"ok": True,
                     "receptor": {"path": "clean/3OG7_receptor.pdb", "format": "pdb",
                                  "chains_auth": ["A"], "atoms": 2100,
                                  "kept_nonpolymer_label_asym_ids": []},
                     "metrics": {"r_free": 0.2575,
                                 "chains": [{"auth_asym_id": "A", "label_asym_id": "A",
                                             "modelled_uniprot_range": "444-721",
                                             "missing_uniprot_ranges": ["597-614"]}]},
                     "pockets": [{"comp_id": "032", "label_asym_id": "C", "auth_asym_id": "A",
                                  "heavy_atoms": 33, "n_residues": 22, "large_box": False,
                                  "spans_target_chains": False, "residues_by_chain": {"A": 22},
                                  "cutoff": 5.0, "residues": [], "partner_residues": [],
                                  "adjacent_missing_ranges": [], "alphafold_plddt": None,
                                  "box": {"center": [2.6, -2.3, -19.4], "size": [26.3, 15.1, 16.4]},
                                  "sdf": "ligands/3OG7_032_C.sdf"}]},
    }
    if runs:
        structure["docking"] = {"032 C": {"runs": {"032_regenerated": {
            "computed_at": "2026-11-11T00:00:00Z", "generated_by": "dock_ligand",
            "poses": [{"mode": 1, "score": -11.0}],
            "files": {"poses_032_regenerated.sdf":
                      {"path": "dock/3OG7_032_C/poses_032_regenerated.sdf"}}}}}}
    manifest = {"version": 1, "created_at": "2026-11-11T00:00:00Z",
                "generated_by": "prepare_structures", "mutation_filter": "V600E",
                "sequence": {"header": "sp|P15056|BRAF_HUMAN"},
                "selection": "1 of 1", "structures": [structure],
                "user_choice": None, "user_choice_history": [], "warnings": [],
                "decisions": [], "pending": []}
    choice_mod.write_manifest(choice_mod.manifest_path("P15056"), manifest)
    return prep, choice_mod


def run_tool(tmp, *, pages=None, molecule=None, **kw):
    import importlib

    from protein_mcp.files import admet as admet_mod
    prep, _ = temp_dossier(tmp, runs=kw.pop("runs", True))
    importlib.reload(admet_mod)
    original = chembl._get
    chembl._get = fake_get(pages=pages, molecule=molecule)
    try:
        return admet_mod, prep, run(admet_mod.analyze_admet("P15056", "3OG7", "032 C", **kw))
    finally:
        chembl._get = original
        os.environ.pop("PROTEIN_MCP_DOSSIER_DIR", None)


def test_the_tool_runs_end_to_end_and_records_the_block():
    with tempfile.TemporaryDirectory() as tmp:
        mod, prep, res = run_tool(Path(tmp))
        assert res["ok"] is True, res
        ident = res["compound_identity"]
        assert ident["inchikey"] == VEM_KEY and ident["source"] == "docked pose"
        assert ident["identity_check"].startswith("the docked pose and the crystallographic")
        assert ident["heavy_atoms_agree"] is True and ident["heavy_atoms_in_manifest_pocket"] == 33
        assert res["chembl_record"]["status"].startswith("RECORD")
        assert res["measured_adme"]["status"] == "MEASURED (retrieved from ChEMBL)"
        assert len(res["measured_adme"]["adme_rows"]) == 6
        assert res["rule_reading"]["results"]["Lipinski"]["violations"] == 1
        assert "4 of 18" not in res["selection"]        # this fixture is 6 of 20
        assert "6 of 20 rows read" in res["selection"], res["selection"]
        assert "not on the allowlist" in res["selection"]
        # recorded in the manifest, and readable back
        m = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
        block = m["structures"][0]["admet"]["032 C"]
        assert block["compound_identity"]["inchikey"] == VEM_KEY
        assert block["generated_by"] == "analyze_admet"
        assert res["recorded_in"].endswith("manifest.json")
        # compact() must keep every block and only drop bulk
        small = mod.compact(res)
        assert "functions" not in small["calculated_properties"]["rdkit"]
        assert small["measured_adme"]["adme_rows"] == res["measured_adme"]["adme_rows"]
        assert small["chembl_record"] == res["chembl_record"]


def test_the_tool_output_never_calls_anything_predicted():
    """Now walks the REAL output of the tool, not a dict assembled by the test."""
    with tempfile.TemporaryDirectory() as tmp:
        _mod, _prep, res = run_tool(Path(tmp))
    for path, value in _strings(res):
        for sentence in value.replace(";", ".").split("."):
            if "predict" in sentence.lower():
                assert _NEGATED.search(sentence), f"{path}: {sentence}"
    assert "PREDICTED" not in res["measured_adme"]["status"]
    assert "PREDICTED" not in res["calculated_properties"]["status"]


def test_the_tool_never_merges_the_two_engines():
    """Also on the real output: a mean of the two engines would appear here, and the test of the
    hand-built block could not have seen it."""
    with tempfile.TemporaryDirectory() as tmp:
        _mod, _prep, res = run_tool(Path(tmp))
    calc = res["calculated_properties"]
    assert set(calc["rdkit"]) & {"alogp", "psa", "full_mwt"} == set(), "ChEMBL keys inside rdkit"
    assert "alogp" in calc["chembl"] and "clogp" in calc["rdkit"]
    for path, _v in _strings(res):
        assert "agree" not in path.lower()
    banned = {"mean", "average", "median", "aggregate", "pooled", "consensus", "combined"}
    for key in _keys(res):
        assert key.lower() not in banned, key


def test_the_links_are_built_from_chembl_ids():
    with tempfile.TemporaryDirectory() as tmp:
        _mod, _prep, res = run_tool(Path(tmp))
    assert res["chembl_record"]["link"] == \
        "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL1229517/"
    assert VEM_KEY not in res["chembl_record"]["link"], "the link must not be built from the key"
    for r in res["measured_adme"]["adme_rows"]:
        assert r["assay_link"].endswith(f"/{r['assay_chembl_id']}/")
        assert r["document_link"].endswith(f"/{r['document_chembl_id']}/")
        assert "assay_report_card" in r["assay_link"]
    ids = {r["assay_chembl_id"] for r in res["measured_adme"]["adme_rows"]}
    for row in ADME_ROWS_A + ADME_ROWS_P:
        if chembl.adme_class(row["standard_type"]):
            assert row["assay_chembl_id"] in ids


def test_a_truncated_lookup_says_so_in_the_selection_and_the_report():
    """CONTROL: with more rows than the cap, the counts are a lower bound and both the tool and
    the document have to say it. They did not, and a 1143-row molecule is real (CHEMBL53463)."""
    many = {"A": [dict(ADME_ROWS_A[0], assay_chembl_id=f"CHEMBL{i:07d}") for i in range(60)],
            "P": list(ADME_ROWS_P)}
    old_cap = chembl.MAX_ADME_ROWS
    chembl.MAX_ADME_ROWS = 25
    try:
        with tempfile.TemporaryDirectory() as tmp:
            _mod, prep, res = run_tool(Path(tmp), pages=many)
            meas = res["measured_adme"]
            assert meas["truncated"] is True
            assert meas["rows_read"] == 25, meas["rows_read"]   # a cap on ROWS, not on offset
            # the cap ran out inside type A, so type P was NEVER asked about: its total must be
            # None ("unknown"), never 0, and the overall total is itself a lower bound
            assert meas["rows_total_by_assay_type"] == {"A": 60, "P": None}
            assert meas["assay_types_skipped_at_the_row_cap"] == ["P"]
            assert meas["rows_total"] == 60 and meas["rows_not_read"] == 35
            assert meas["rows_total_is_lower_bound"] is True
            assert "P" not in meas["assay_types_queried"]
            assert "P" in meas["assay_types_not_queried"]
            assert "UNKNOWN, not zero" in meas["note"]
            assert "TRUNCATED" in res["selection"] and "LOWER BOUND" in res["selection"]
            assert "TRUNCATED" in meas["note"]
            md = _render_section_three(prep)
            assert "TRUNCATED" in md and "LOWER BOUND" in md
    finally:
        chembl.MAX_ADME_ROWS = old_cap


def _render_section_three(prep: Path) -> str:
    import importlib

    from protein_mcp import experiment_report as er
    importlib.reload(er)
    res = er.build_experiment_report("P15056")
    assert res.get("ok"), res
    md = Path(res["markdown_path"]).read_text(encoding="utf-8")
    return md[md.index("## 3. ADMET"):md.index("## 4.")]


def test_the_report_section_three_is_rendered_from_the_block():
    with tempfile.TemporaryDirectory() as tmp:
        _mod, prep, res = run_tool(Path(tmp))
        sec = _render_section_three(prep)
    assert "Not run" not in sec, "the section must have content, not the declared gap"
    assert VEM_KEY in sec and "CHEMBL1229517" in sec
    for heading in ("Compound identity", "ChEMBL record", "Measured ADME",
                    "Calculated properties", "Property rules"):
        assert heading in sec, heading
    # the two clearance rows, with their incompatible units, each printed once and unconverted
    assert "uL.min-1.(10^6cells)-1" in sec and "ml kg-1" in sec
    assert "6.61" in sec and "12.0" in sec and "9.3" not in sec
    assert "> 120.0" in sec, "a censored value must keep its relation in the document"
    assert "READING OFFERED, NOT MEASURED" in sec
    assert "molecular_weight <= 500" in sec
    # 'predict' may appear in the document ONLY inside a denial
    for sentence in sec.replace("\n", " ").replace(";", ".").split("."):
        if "predict" in sentence.lower():
            assert _NEGATED.search(sentence), sentence


def test_the_report_never_calls_an_excluded_row_not_adme():
    """CONTROL THAT MUST FAIL: the excluded `Ratio` rows of vemurafenib are P-gp and BCRP efflux
    ratios and one excluded `Activity` row is CYP450 induction -- all of them ADME. The document
    used to tell the analyst they were "NOT ADME"."""
    with tempfile.TemporaryDirectory() as tmp:
        _mod, prep, _res = run_tool(Path(tmp))
        sec = _render_section_three(prep)
    assert "NOT ADME" not in sec, "the report claims the source's data is not ADME"
    assert "that are not ADME" not in sec
    assert "NOT on this project's allowlist" in sec
    assert "does not mean they are not ADME" in sec
    assert "Activity (9)" in sec and "Ratio (2)" in sec, "the excluded rows must still be counted"
    # F3: and the rows the paragraph warns about are now in the document the analyst reads,
    # under a heading that says they are unclassified and not part of the ADME count
    assert "The excluded rows themselves" in sec and "unclassified" in sec, sec[-1500:]
    assert "NOT part of the count above" in sec
    for fixture_row in ADME_ROWS_A:
        if not chembl.adme_class(fixture_row["standard_type"]):
            head = (fixture_row["assay_description"] or "")[:40]
            assert head and head in sec, f"excluded row missing from the report: {head}"


def test_the_record_defaults_are_marked_in_the_tool_and_the_report():
    """A False flag on an unannotated entry is not a finding, and the -1 fields are not either."""
    page = json.loads(json.dumps(MOLECULE_PAGE))
    mol = page["molecules"][0]
    mol["pref_name"] = None
    mol["max_phase"] = None
    mol["molecule_type"] = None
    mol["atc_classifications"] = []
    for f in ("prodrug", "chirality", "availability_type", "first_in_class", "inorganic_flag",
              "orphan", "veterinary"):
        mol[f] = -1
    for f in ("oral", "parenteral", "topical", "therapeutic_flag", "dosed_ingredient",
              "withdrawn_flag"):
        mol[f] = False
    with tempfile.TemporaryDirectory() as tmp:
        _mod, prep, res = run_tool(Path(tmp), molecule=page)
        rec = res["chembl_record"]
        assert len(rec["fields_minus_one"]) == 7, rec["fields_minus_one"]
        assert "oral" in rec["boolean_fields_false"] and "withdrawn_flag" in rec["boolean_fields_false"]
        assert rec["fields"]["oral"] is False, "the raw value must be kept"
        sec = _render_section_three(prep)
    assert "(default, see below)" in sec and "(-1, see below)" in sec
    assert "none of them is evidence" in sec
    # and no invented name reached the document
    for name in ("VEMURAFENIB", "PLX-4032", "ZELBORAF"):
        assert name not in sec, name


def test_a_failed_molecule_lookup_is_unknown_in_the_tool_not_an_exception():
    """It used to escape as a RuntimeError through the MCP boundary while ESTADO claimed the
    block went to UNKNOWN."""
    async def boom(client, endpoint, params, deadline=None):
        raise RuntimeError("ChEMBL molecule failed: HTTP 500")
    with tempfile.TemporaryDirectory() as tmp:
        import importlib

        from protein_mcp.files import admet as admet_mod
        temp_dossier(Path(tmp))
        importlib.reload(admet_mod)
        original = chembl._get
        chembl._get = boom
        try:
            res = run(admet_mod.analyze_admet("P15056", "3OG7", "032 C"))
        finally:
            chembl._get = original
            os.environ.pop("PROTEIN_MCP_DOSSIER_DIR", None)
    assert res["ok"] is True
    assert res["chembl_record"]["status"] == "UNKNOWN (the ChEMBL lookup failed)"
    assert res["measured_adme"]["status"] == "UNKNOWN (nothing was retrieved)"
    assert res["measured_adme"]["retrieved"] is False
    assert "UNKNOWN" in res["selection"] and "FAILED" in res["selection"]
    assert "no record" not in res["selection"].lower()
    # the calculation does not depend on ChEMBL and must still be there
    assert res["calculated_properties"]["rdkit"]["heavy_atoms"] == 33
    assert res["rule_reading"]["results"]["Lipinski"]["violations"] == 1


def test_a_compound_absent_from_chembl_is_not_reported_as_a_failed_lookup():
    empty = {"page_meta": {"total_count": 0}, "molecules": []}
    with tempfile.TemporaryDirectory() as tmp:
        _mod, prep, res = run_tool(Path(tmp), molecule=empty)
        assert res["chembl_record"]["status"] == "NO RECORD"
        assert res["measured_adme"]["no_chembl_entry"] is True
        assert "has no ChEMBL entry" in res["selection"]
        assert "FAILED" not in res["selection"]
        sec = _render_section_three(prep)
    assert "no ChEMBL molecule has this exact" in sec
    assert "the ChEMBL lookup did not complete" not in sec


def test_a_calculated_value_in_the_measured_block_is_flagged():
    """Real rows exist whose own description says the value was calculated: 'Acid dissociation
    constant was calculated based on carbon atom CH2'."""
    pages = {"A": [{"standard_type": "pKa", "standard_relation": "=", "standard_value": "7.4",
                    "standard_units": None, "assay_chembl_id": "CHEMBL9999040",
                    "assay_description": "Acid dissociation constant was calculated based on "
                                         "carbon atom CH2",
                    "document_chembl_id": "CHEMBL9999041"}] + ADME_ROWS_A[:1], "P": []}
    with tempfile.TemporaryDirectory() as tmp:
        _mod, prep, res = run_tool(Path(tmp), pages=pages)
        rows = res["measured_adme"]["adme_rows"]
        pka = next(r for r in rows if r["standard_type"] == "pKa")
        cl = next(r for r in rows if r["standard_type"] == "CL")
        assert pka["description_says_value_was_calculated"] is True
        assert pka["calculated_wording_found"] == ["calculated"]
        assert cl["description_says_value_was_calculated"] is False
        assert "say in their own assay description that the value was CALCULATED" in res["selection"]
        sec = _render_section_three(prep)
    assert "[depositor says calculated]" in sec


def test_the_applied_cutoffs_and_the_published_cutoffs_cannot_drift():
    """CONTROL: the numbers APPLIED and the strings PUBLISHED used to be two separate literals,
    so the pre-registration test only checked the ones the report prints."""
    for name, items in admet.RULE_LIMITS.items():
        published = admet.RULES[name]["cutoffs"]
        assert published == [f"{f} <= {lim}" for f, lim in items], name
    # and the applied numbers are the pre-registered ones (ESTADO 7.3)
    assert admet.RULE_LIMITS["Lipinski"] == [("molecular_weight", 500), ("clogp", 5),
                                             ("hbd", 5), ("hba", 10)]
    assert admet.RULE_LIMITS["Veber"] == [("rotatable_bonds", 10), ("tpsa", 140)]
    assert admet.RULE_LIMITS["Egan"] == [("tpsa", 131.6), ("clogp", 5.88)]


def test_a_value_between_two_rules_cutoffs_separates_them():
    """A TPSA of 135 is over Egan's 131.6 and under Veber's 140. Mutating one cut-off to the
    other's value is invisible without a case in between."""
    base = {"molecular_weight": 400.0, "clogp": 3.0, "hbd": 1, "hba": 2, "rotatable_bonds": 3}
    r = admet.rule_reading({**base, "tpsa": 135.0})["results"]
    assert r["Egan"]["violations"] == 1 and r["Egan"]["violated"] == ["tpsa 135.0 > 131.6"]
    assert r["Veber"]["violations"] == 0
    r = admet.rule_reading({**base, "tpsa": 145.0})["results"]
    assert r["Egan"]["violations"] == 1 and r["Veber"]["violations"] == 1
    # cLogP between Lipinski's 5 and Egan's 5.88
    r = admet.rule_reading({**base, "tpsa": 50.0, "clogp": 5.5})["results"]
    assert r["Lipinski"]["violations"] == 1 and r["Egan"]["violations"] == 0


def test_a_value_exactly_on_a_cutoff_is_not_a_violation():
    """'<= 500' means 500 passes. A '>' mutated to '>=' is invisible without this case."""
    r = admet.rule_reading({"molecular_weight": 500.0, "clogp": 5.0, "hbd": 5, "hba": 10,
                            "rotatable_bonds": 10, "tpsa": 140.0})["results"]
    assert r["Lipinski"]["violations"] == 0, r["Lipinski"]
    assert r["Veber"]["violations"] == 0, r["Veber"]
    r = admet.rule_reading({"molecular_weight": 500.01, "clogp": 5.0, "hbd": 5, "hba": 10,
                            "rotatable_bonds": 10, "tpsa": 140.0})["results"]
    assert r["Lipinski"]["violations"] == 1


def test_a_missing_file_and_an_escaping_path_are_told_apart():
    with tempfile.TemporaryDirectory() as tmp:
        prep = Path(tmp) / "prep"
        prep.mkdir()
        write_sdf(Path(tmp) / "elsewhere.sdf", VEM_SMILES)
        try:
            admet._inside(prep, "../elsewhere.sdf")
            raise AssertionError("an escaping path must be refused")
        except admet.AdmetError as exc:
            assert "outside prep/" in str(exc) and "missing" not in str(exc)
        try:
            admet._inside(prep, "absent.sdf")
            raise AssertionError("a missing file must be refused")
        except admet.AdmetError as exc:
            assert "missing" in str(exc) and "outside" not in str(exc)


def _keys(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from _keys(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _keys(v)


def _strings(obj, path=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _strings(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _strings(v, f"{path}[{i}]")
    elif isinstance(obj, str):
        yield path, obj


# ---- F1: the time budget (Desktop test of 2026-11-12) --------------------------------------
class _Resp:
    def __init__(self, code):
        self.status_code = code

    def raise_for_status(self):
        raise AssertionError("not reached: the fixture only returns codes _get retries on")

    def json(self):
        return {}


class _SlowClient:
    """A client that always answers 500, like ChEMBL did three times of three at limit=200."""

    def __init__(self):
        self.attempts = []

    async def get(self, url, params=None, **kw):
        self.attempts.append(kw.get("timeout"))
        return _Resp(500)


def test_the_budget_stops_the_retries_instead_of_waiting_three_times():
    """CONTROL: with a spent budget the first attempt is the last, and the error says so.

    Without the budget _get retried three times with 1.5 + 3 + 4.5 s of back-off on top of a
    60 s read timeout each: 189 s for one endpoint, while the caller waits. The second half of
    this test is the control that tells the two apart -- it asserts the OLD behaviour is still
    what happens when no budget is passed, so this test cannot pass by the budget being ignored.
    """
    client = _SlowClient()
    spent = chembl.Deadline(0)                      # already expired
    try:
        run(chembl._get(client, "activity", {"x": 1}, spent))
        raise AssertionError("a spent budget must not be queried at all")
    except RuntimeError as exc:
        assert "gave up" in str(exc) and "budget" in str(exc), exc
    assert client.attempts == [], "nothing may be requested after the budget is gone"

    # control: no budget -> the old three attempts (sleeps neutralised so the suite stays fast)
    naps = []

    async def _no_sleep(sec):
        naps.append(sec)

    original = asyncio.sleep
    asyncio.sleep = _no_sleep
    try:
        client2 = _SlowClient()
        try:
            run(chembl._get(client2, "activity", {"x": 1}))
            raise AssertionError("three 500s must still fail")
        except RuntimeError as exc:
            assert "HTTP 500" in str(exc)
        assert len(client2.attempts) == 3, client2.attempts
        assert naps == [1.5, 3.0, 4.5], naps
        assert client2.attempts == [None, None, None], "no budget means no per-request override"
    finally:
        asyncio.sleep = original


def test_a_request_never_waits_past_the_budget():
    """The per-request timeout is capped by what is left, so one slow call cannot eat the rest."""
    client = _SlowClient()
    naps = []

    async def _no_sleep(sec):
        naps.append(sec)

    original = asyncio.sleep
    asyncio.sleep = _no_sleep
    try:
        run(chembl._get(client, "activity", {"x": 1}, chembl.Deadline(0.5)))
    except RuntimeError:
        pass
    finally:
        asyncio.sleep = original
    assert client.attempts, "the first attempt must happen"
    for t in client.attempts:
        # t.read <= READ_TIMEOUT compares the threshold against itself and cannot fail
        # (reviewer finding 2026-11-12): what matters is that it follows what is LEFT
        assert t is not None and t.read <= 0.5 + 1.0, t
    assert all(n <= 0.5 for n in naps), f"back-off outlived the budget: {naps}"


def test_the_measured_numbers_behind_the_budget_are_pinned():
    """CONTROL: READ_TIMEOUT=600 and ADME_BUDGET=4000 passed the whole suite (reviewer finding
    2026-11-12), so the two numbers that make F1 work could drift in silence. They are the
    reason the Desktop call stopped hanging, and each has a measurement behind it: a healthy
    ChEMBL answers in under 2 s, a failing one took 63.2 s, and the four fixed calls measured
    26.7-40.1 s. If either number is changed, change this test and say why in ESTADO."""
    assert chembl.READ_TIMEOUT == 20.0, chembl.READ_TIMEOUT
    assert chembl.CONNECT_TIMEOUT == 10.0, chembl.CONNECT_TIMEOUT
    assert chembl.ADME_BUDGET == 40.0, chembl.ADME_BUDGET
    # the whole point: the budget plus the connect slack has to stay under a chat client's window
    assert chembl.ADME_BUDGET + chembl.CONNECT_TIMEOUT <= 55.0
    assert chembl.READ_TIMEOUT < chembl.ADME_BUDGET, "one request must not be able to eat it all"
    # get_ligands keeps the long timeout: pages that END WELL above 20 s were measured
    assert chembl.TIMEOUT == 60.0, "get_ligands must not inherit the ADME read cap"
    assert isinstance(chembl.ADME_TIMEOUT, httpx.Timeout)


def test_the_whole_admet_call_shares_one_budget():
    """Every ChEMBL request of one analyze_admet call carries the SAME budget object: two
    independent budgets would add up to twice the wait the caller was promised."""
    import importlib

    from protein_mcp.files import admet as admet_mod
    with tempfile.TemporaryDirectory() as tmp:
        temp_dossier(Path(tmp))
        importlib.reload(admet_mod)
        fake = fake_get()
        original = chembl._get
        chembl._get = fake
        try:
            res = run(admet_mod.analyze_admet("P15056", "3OG7", "032 C"))
        finally:
            chembl._get = original
            os.environ.pop("PROTEIN_MCP_DOSSIER_DIR", None)
    assert res["ok"] is True, res
    seen = [c[2] for c in fake.calls]
    assert len(fake.calls) >= 3, fake.calls          # molecule + one page per assay type
    assert seen, "no call recorded"
    assert all(isinstance(d, chembl.Deadline) for d in seen), seen
    assert len({id(d) for d in seen}) == 1, "the molecule lookup and the activities must share it"


def test_a_type_left_unread_by_the_clock_is_not_called_empty():
    """CONTROL, and the honest half of F1: a type the clock never reached is UNKNOWN for a
    DIFFERENT reason than the row cap, and the two must not be reported as one. Blaming the cap
    for a slow source would be a false statement about the data, which is the failure mode this
    whole step exists to avoid."""

    class ExpiresAfter(chembl.Deadline):
        def __init__(self, n):
            super().__init__(10_000)
            self.left = n

        def expired(self):
            if self.left <= 0:
                return True
            self.left -= 1
            return False

    original = chembl._get
    chembl._get = fake_get()
    try:
        # one clean pass (assay type A), then the clock is gone before P is asked about at all
        res = run(chembl.fetch_adme_activities("CHEMBL1229517", deadline=ExpiresAfter(1)))
    finally:
        chembl._get = original
    assert res["rows_total_by_assay_type"]["P"] is None, res["rows_total_by_assay_type"]
    assert res["assay_types_skipped_at_the_time_budget"] == ["P"], res
    assert res["assay_types_skipped_at_the_row_cap"] is None, "the cap was not the reason"
    assert res["rows_total_is_lower_bound"] is True and res["truncated"] is True
    assert "time budget ran out" in res["note"] and "UNKNOWN, not zero" in res["note"]
    assert "P" not in res["assay_types_queried"], res["assay_types_queried"]
    assert res["assay_types_not_queried"].get("P") == "physicochemical"


# ---- F2: no name from any retrieved source -------------------------------------------------
def test_a_budget_spent_before_any_query_is_unknown_not_zero_rows():
    """CONTROL THAT MUST FAIL, and it did: F1 itself created a new way to say zero. With the
    budget gone before any assay type could be asked about, this returned
    retrieved=True / rows_total=0 / "MEASURED (retrieved from ChEMBL)" and the selection read
    "ChEMBL holds 0 rows in the assay types queried" -- indistinguishable from TXV's REAL zero,
    which is the regression this whole step exists to prevent (reviewer finding 2026-11-12)."""
    original = chembl._get
    chembl._get = fake_get()
    try:
        spent = chembl.Deadline(chembl.ADME_BUDGET)
        spent._start -= chembl.ADME_BUDGET + 1          # the molecule lookup ate all of it
        res = run(chembl.fetch_adme_activities("CHEMBL1229517", deadline=spent))
    finally:
        chembl._get = original
    assert res["retrieved"] is False, res
    assert "rows_total" not in res, "a total of 0 must not be reported when nothing was asked"
    assert res["assay_types_skipped_at_the_time_budget"] == ["A", "P"], res
    low = res["note"].lower()
    assert "unknown" in low and "not 'no data'" in low and "not zero" in low, res["note"]
    assert "time budget ran out" in low
    # and the tool's own status must follow it to UNKNOWN, not MEASURED
    with tempfile.TemporaryDirectory() as tmp:
        import importlib

        from protein_mcp.files import admet as admet_mod
        temp_dossier(Path(tmp))
        importlib.reload(admet_mod)
        chembl._get = fake_get()
        real_fetch = chembl.fetch_adme_activities

        async def spent_fetch(cid, client=None, deadline=None):
            d = chembl.Deadline(chembl.ADME_BUDGET)
            d._start -= chembl.ADME_BUDGET + 1
            return await real_fetch(cid, client=client, deadline=d)

        chembl.fetch_adme_activities = spent_fetch
        try:
            res2 = run(admet_mod.analyze_admet("P15056", "3OG7", "032 C"))
            sec = _render_section_three(Path(tmp) / "dossiers" / "P15056" / "prep")
        finally:
            chembl.fetch_adme_activities = real_fetch
            chembl._get = original
            os.environ.pop("PROTEIN_MCP_DOSSIER_DIR", None)
    assert res2["measured_adme"]["status"] == "UNKNOWN (nothing was retrieved)", res2
    # the document is what the analyst reads, and it printed none of this
    assert "UNKNOWN, not 'no data' and not zero rows" in sec, sec[-1500:]
    assert "time budget ran out" in sec
    assert "No row in the assay types queried is on this project's ADME allowlist" not in sec


def test_a_clock_that_stops_mid_type_does_not_blame_the_row_cap():
    """CONTROL: with the clock dying after a type's first page, `total` is known, so the type is
    not 'skipped' -- and the truncation sentence used to say "(cap 300)" with the cap nowhere
    near. Blaming our own cap for a slow source is a false statement about the data."""

    class ExpiresAfter(chembl.Deadline):
        def __init__(self, n):
            super().__init__(10_000)
            self.left = n

        def expired(self):
            if self.left <= 0:
                return True
            self.left -= 1
            return False

    many = {"A": [dict(ADME_ROWS_A[0], assay_chembl_id=f"CHEMBL{i:07d}") for i in range(60)],
            "P": list(ADME_ROWS_P)}
    original = chembl._get
    chembl._get = fake_get(pages=many)
    try:
        # two clean page reads of type A, then the clock dies with its total already known
        res = run(chembl.fetch_adme_activities("CHEMBL1229517", deadline=ExpiresAfter(2)))
    finally:
        chembl._get = original
    assert res["truncated"] is True and res["rows_read"] < res["rows_total"]
    assert res["assay_types_cut_off_at_the_time_budget"] == ["A"], res
    assert res["rows_read"] < chembl.MAX_ADME_ROWS, "the cap was never near"
    note = res["note"]
    assert "the time budget ran out, NOT the row cap" in note, note
    assert f"(cap {chembl.MAX_ADME_ROWS})" not in note


def test_a_compound_no_source_names_is_not_named_and_not_described():
    """CONTROL for the hallucination of the Desktop test of 2026-11-12: with the lookup down and
    the pocket name None, a compound was reported as 'the chemically linked vemurafenib
    inhibitor (VEM-3-VEM)' -- a name and a description from no source at all."""
    from protein_mcp.files import admet as admet_mod
    failed = admet_mod._name_status({"lookup_failed": True})
    absent = admet_mod._name_status({"found": False})
    unnamed = admet_mod._name_status({"found": True, "chembl_id": "CHEMBL5440746",
                                      "pref_name": None})
    named = admet_mod._name_status({"found": True, "chembl_id": "CHEMBL1229517",
                                    "pref_name": "VEMURAFENIB"})
    assert failed.startswith("UNKNOWN") and "NO NAME IS AVAILABLE" in failed
    assert absent.startswith("NO NAME IN CHEMBL")
    assert unnamed.startswith("UNNAMED IN CHEMBL")
    for text in (failed, absent, unnamed):
        low = text.lower()
        assert "do not invent one" in low
        assert "do not give it a code of your own" in low
        # the exact shapes the invented description took
        for word in ("analogue", "dimer", "linked"):
            assert word in low, (word, text)
        # and the claim stays inside what this output can support: wwPDB's systematic name WAS
        # retrieved by get_structures and is in the manifest, so "no retrieved source names this
        # compound" was an over-claim of the kind this project punishes (reviewer, 2026-11-12)
        assert "this output holds no name" in low, text
        assert "wwpdb" in low and "systematic" in low, text
        assert "no retrieved source" not in low, "the over-claim is back"
    assert named.startswith("NAMED BY THE SOURCE") and "VEMURAFENIB" in named
    assert "do not invent one" not in named.lower()
    # and it must not send the reader to a synonym list that does not exist
    assert "synonyms are in" not in named


def test_the_unnamed_case_carries_the_rule_in_the_tool_output():
    """On the REAL output of the tool, in both branches that lack a name."""
    with tempfile.TemporaryDirectory() as tmp:
        _mod, _prep, res = run_tool(Path(tmp), molecule={"molecules": [],
                                                         "page_meta": {"total_count": 0}})
    assert res["chembl_record"]["status"] == "NO RECORD"
    assert res["chembl_record"]["name_status"].startswith("NO NAME")
    assert any("holds NO name for this compound" in n for n in res["notes"]), res["notes"]

    with tempfile.TemporaryDirectory() as tmp:
        _mod, _prep, res = run_tool(Path(tmp))     # vemurafenib: a real pref_name
    assert res["chembl_record"]["name_status"].startswith("NAMED BY THE SOURCE")
    assert not any("holds NO name for this compound" in n for n in res["notes"])


# ---- F3: the excluded rows travel with their own description -------------------------------
def test_an_excluded_row_travels_with_its_own_description():
    """F3: counts alone left the analyst unable to see WHAT was excluded, and some of those rows
    hold real ADME (P-gp and BCRP efflux ratios, CYP450 induction). They come back one by one,
    unclassified, and they must not be added to the allowlisted rows."""
    res, _calls = adme()
    counts = res["rows_excluded_by_allowlist_by_standard_type"]
    rows = res["rows_excluded_by_allowlist"]
    assert rows and len(rows) == sum(counts.values()), (len(rows), counts)
    kept_ids = {r["assay_chembl_id"] for r in res["adme_rows"]}
    for r in rows:
        assert r["on_adme_allowlist"] is False
        assert "adme_class" not in r, "an excluded row must not be classified"
        assert r["assay_description"], "the description is the whole point of F3"
        assert chembl.adme_class(r["standard_type"]) is None
        assert r["assay_chembl_id"] not in kept_ids, "an excluded row is not an ADME row"
        assert r["assay_link"].endswith(f"/{r['assay_chembl_id']}/")
    assert res["rows_excluded_shown"] == f"{len(rows)} of {sum(counts.values())}"
    low = res["rows_excluded_note"].lower()
    assert "not on this project's adme allowlist" in low
    assert "do not call them 'not adme'" in low
    assert "do not count them as adme measurements" in low
    # the selection is unchanged: showing them must not change what counts as ADME
    assert len(res["adme_rows"]) + len(rows) == res["rows_read"]


def test_the_excluded_rows_are_capped_and_say_so():
    """CONTROL: a 1143-row molecule is real, so the list is capped and the cap is reported."""
    many = {"A": [dict(ADME_ROWS_A[0], standard_type="Activity",
                       assay_chembl_id=f"CHEMBL{i:07d}") for i in range(60)],
            "P": []}
    old = chembl.MAX_EXCLUDED_SHOWN
    chembl.MAX_EXCLUDED_SHOWN = 5
    try:
        res, _ = adme(pages=many)
    finally:
        chembl.MAX_EXCLUDED_SHOWN = old
    assert len(res["rows_excluded_by_allowlist"]) == 5
    assert res["rows_excluded_shown"] == "5 of 60"
    assert res["rows_excluded_by_allowlist_by_standard_type"] == {"Activity": 60}


def main():
    test_the_fixture_copies_the_real_defects()
    test_only_the_allowlisted_types_survive()
    test_without_the_allowlist_cytotoxicity_would_be_reported_as_adme()
    test_a_type_that_is_not_adme_has_no_class()
    test_querying_only_assay_type_a_would_lose_the_solubility()
    test_what_was_not_queried_is_named()
    test_the_page_size_is_small_enough_for_the_api()
    test_paging_reads_every_row()
    test_every_row_keeps_its_own_units_and_relation()
    test_averaging_the_two_clearance_rows_would_invent_a_number()
    test_a_connectivity_block_key_is_refused()
    test_a_lowercase_key_is_accepted_and_normalised()
    test_two_molecules_with_the_same_key_choose_neither()
    test_an_unnamed_research_compound_is_not_given_a_name()
    test_a_minus_one_flag_is_named_and_not_interpreted()
    test_the_descriptors_name_the_function_that_produced_them()
    test_the_two_engines_are_reported_apart_and_never_merged()
    test_vemurafenib_breaks_lipinski_on_logp_only()
    test_the_rules_are_offered_not_measured()
    test_a_missing_descriptor_makes_the_count_unknown_not_zero()
    test_a_molecule_over_every_cutoff_counts_every_violation()
    test_the_cutoffs_are_the_pre_registered_ones()
    test_the_pose_and_the_crystal_sdf_give_the_same_key()
    test_a_mismatched_pose_and_crystal_never_reach_chembl()
    test_without_a_docking_run_the_crystal_sdf_is_used_and_marked()
    test_a_pose_file_that_changed_since_the_manifest_is_refused()
    test_a_file_outside_prep_is_refused()
    test_nothing_in_the_output_is_called_predicted()
    test_a_failed_lookup_is_unknown_not_no_data()
    test_a_molecule_with_no_adme_rows_says_how_many_were_read()
    # end to end (added after the external review of 2026-11-11)
    test_the_tool_runs_end_to_end_and_records_the_block()
    test_the_tool_output_never_calls_anything_predicted()
    test_the_tool_never_merges_the_two_engines()
    test_the_links_are_built_from_chembl_ids()
    test_a_truncated_lookup_says_so_in_the_selection_and_the_report()
    test_the_report_section_three_is_rendered_from_the_block()
    test_the_report_never_calls_an_excluded_row_not_adme()
    test_the_record_defaults_are_marked_in_the_tool_and_the_report()
    test_a_failed_molecule_lookup_is_unknown_in_the_tool_not_an_exception()
    test_a_compound_absent_from_chembl_is_not_reported_as_a_failed_lookup()
    test_a_calculated_value_in_the_measured_block_is_flagged()
    test_the_applied_cutoffs_and_the_published_cutoffs_cannot_drift()
    test_a_value_between_two_rules_cutoffs_separates_them()
    test_a_value_exactly_on_a_cutoff_is_not_a_violation()
    test_a_missing_file_and_an_escaping_path_are_told_apart()
    # F1, F2, F3 (Desktop test of 2026-11-12)
    test_the_budget_stops_the_retries_instead_of_waiting_three_times()
    test_a_request_never_waits_past_the_budget()
    test_the_measured_numbers_behind_the_budget_are_pinned()
    test_the_whole_admet_call_shares_one_budget()
    test_a_type_left_unread_by_the_clock_is_not_called_empty()
    test_a_budget_spent_before_any_query_is_unknown_not_zero_rows()
    test_a_clock_that_stops_mid_type_does_not_blame_the_row_cap()
    test_a_compound_no_source_names_is_not_named_and_not_described()
    test_the_unnamed_case_carries_the_rule_in_the_tool_output()
    test_an_excluded_row_travels_with_its_own_description()
    test_the_excluded_rows_are_capped_and_say_so()
    print("ALL OK")


if __name__ == "__main__":
    main()
