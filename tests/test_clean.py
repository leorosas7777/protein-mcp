"""Offline test for protein_mcp/files/clean.py (receptor, chain metrics, pocket, box, AlphaFold pLDDT).
Fixtures are REAL files (2026-10-22), gzipped: 3OG7 (BRAF V600E + vemurafenib; chain A has
597-614 unmodelled, chain B has GLU 600), 8TVK (KRAS G12C + GDP + Mg, no drug) and the
AlphaFold model AF-P15056-F1 v6, 8C7Y (two TXV, crystal-symmetry contacts), UniProt FASTA of
P15056 and P01116, plus RCSB chain replies for the manifest records.
Expected values were read from the real files with gemmi during the 2026-10-22 validation.
Run: uv run python tests/test_clean.py
"""
import json
import tempfile
from pathlib import Path

from protein_mcp.files import clean
from protein_mcp.sources import rcsb

FIX = Path(__file__).parent / "fixtures" / "prep"


def canonical(acc: str) -> str:
    lines = (FIX / f"uniprot_{acc}.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines[1:])


BRAF, KRAS = canonical("P15056"), canonical("P01116")


def records(fixture: str, accession: str) -> dict:
    data = json.loads((FIX / fixture).read_text(encoding="utf-8"))
    out = {}
    for e in data["data"]["entries"]:
        rec = rcsb._parse_chains(e, accession)
        rec["pdb_id"] = e["rcsb_id"]
        for inst in rec["nonpolymer_instances"]:        # same shape prepare_structures stores
            inst["neighbor_target_chains"] = {c: {"auth_asym_id": None, "closest_angstrom": d["closest_angstrom"]}
                                              for c, d in inst["neighbor_target_chains"].items()}
        out[e["rcsb_id"]] = rec
    return out


def het_names(path: Path):
    return sorted({line[17:20].strip() for line in path.read_text().splitlines() if line.startswith("HETATM")})


def test_3og7(out: Path):
    rec = records("rcsb_chains_3OG7_5JSM_7ZR0.json", "P15056")["3OG7"]
    r = clean.clean_structure(FIX / "3OG7.cif.gz", out, rec, "V600E", BRAF)
    assert r["ok"], r
    assert r["metrics"]["r_free"] == 0.2575 and r["metrics"]["models_in_file"] == 1
    a, b = r["metrics"]["chains"]
    assert (a["auth_asym_id"], b["auth_asym_id"]) == ("A", "B")
    assert a["variant_check"] == {"uniprot_position": 600, "uniprot_wild_type": "V", "observed": None,
                                  "status": "not modelled"}
    assert b["variant_check"]["status"] == "variant" and b["variant_check"]["observed"] == "GLU 600"
    assert [597, 614] in a["missing_uniprot_ranges"] and a["residues_modelled"] == 247
    (p,) = r["pockets"]
    assert (p["comp_id"], p["auth_asym_id"], p["heavy_atoms"]) == ("032", "A", 33)
    assert len(p["residues"]) == 22 and p["partner_residues"] == []
    got = {(x["resname"], x["uniprot"]) for x in p["residues"]}
    assert {("LYS", 483), ("THR", 529), ("CYS", 532), ("ASP", 594), ("PHE", 595), ("GLY", 596)} <= got
    assert all(str(x["uniprot"]) == x["auth_seq_id"] for x in p["residues"])     # BRAF: auth == UniProt here
    assert p["adjacent_missing_ranges"] == [{"auth_asym_id": "A", "missing_uniprot_range": [597, 614]}]
    assert p["box"]["center"] == [2.643, -2.28, -19.403] and p["box"]["size"] == [26.3, 15.13, 16.4]
    assert not p["large_box"]
    receptor = Path(r["receptor"]["path"])
    assert receptor.suffix == ".pdb" and r["receptor"]["chains_auth"] == ["A", "B"]
    assert het_names(receptor) == []                                             # no 032, no water

    # the wild-type letter comes from UniProt, not from the filter text (reviewer: "E600V" said wild type)
    r2 = clean.clean_structure(FIX / "3OG7.cif.gz", out, rec, "E600V", BRAF)
    vb = r2["metrics"]["chains"][1]["variant_check"]
    assert vb["status"] == "other residue" and "UniProt has V at 600" in vb["warning"]

    # RCSB chain ids absent from the file: an error, not an empty receptor with ok true
    wrong = json.loads(json.dumps(rec))
    for c in wrong["target_entities"][0]["chains"]:
        c["label_asym_id"] = "Z" + c["label_asym_id"]
    r3 = clean.clean_structure(FIX / "3OG7.cif.gz", out, wrong, "V600E", BRAF)
    assert not r3["ok"] and "no atom of the target chains" in r3["error"]
    return {"3OG7": r["pockets"]}


def test_8c7y_symmetry(out: Path):
    """NeighborSearch with the unit cell returns symmetry images; they must not enter the pocket.
    Before the fix TXV U (chain A) had B:S630 at 47.42 Å, TXV X had A:R603 at 45.69 Å."""
    rec = records("rcsb_chains_8C7Y.json", "P15056")["8C7Y"]
    r = clean.clean_structure(FIX / "8C7Y.cif.gz", out, rec, "V600E", BRAF)
    assert r["ok"], r
    by_label = {p["label_asym_id"]: p for p in r["pockets"]}
    assert {x["auth_asym_id"] for x in by_label["U"]["residues"]} == {"A"}
    assert {x["auth_asym_id"] for x in by_label["X"]["residues"]} == {"B"}
    assert all(x["min_distance"] <= clean.POCKET_CUTOFF for p in r["pockets"] for x in p["residues"])
    assert all(x.get("reason") for x in r["pocket_not_computed"])                 # every omission says why


def test_8tvk(out: Path):
    rec = records("rcsb_chains_8TVK.json", "P01116")["8TVK"]
    r = clean.clean_structure(FIX / "8TVK.cif.gz", out, rec, "G12C", KRAS)
    assert r["ok"], r
    vc = r["metrics"]["chains"][0]["variant_check"]
    assert vc["status"] == "variant" and vc["uniprot_wild_type"] == "G" and "warning" not in vc   # CYS 12 modelled
    assert r["pockets"] == [] and r["pocket_not_computed"] == []                 # GDP is a nucleotide, not a drug
    assert het_names(Path(r["receptor"]["path"])) == ["GDP", "MG"]                # kept: GDP and its Mg2+
    # wrong filter letters (reviewer: "D13G" called the real GLY 13 a variant): UniProt G13 -> wild type
    r = clean.clean_structure(FIX / "8TVK.cif.gz", out, rec, "D13G", KRAS)
    vc = r["metrics"]["chains"][0]["variant_check"]
    assert vc["status"] == "wild type" and "UniProt has G at 13" in vc["warning"]
    # a bare position / non-substitution filter gives no variant check
    r = clean.clean_structure(FIX / "8TVK.cif.gz", out, rec, "G12", KRAS)
    assert "variant_check" not in r["metrics"]["chains"][0]


def test_alphafold(pockets):
    m = clean.alphafold_metrics(FIX / "AF-P15056-F1-model_v6.cif.gz", pockets, "V600E", BRAF)
    assert m["ok"]
    assert m["variant_position"]["residue"] == "VAL" and m["variant_position"]["plddt"] == 49.06
    assert m["variant_position"]["note"].startswith("Wild-type residue")
    bad = clean.alphafold_metrics(FIX / "AF-P15056-F1-model_v6.cif.gz", pockets, "E600V", BRAF)
    assert bad["variant_position"]["note"].startswith("Wild-type")            # compared with UniProt, not the filter
    assert clean.alphafold_metrics(FIX / "AF-P15056-F1-model_v6.cif.gz", pockets, "V600E", "")[
        "variant_position"]["note"].startswith("Not compared")
    (p,) = m["pockets"]
    assert p["residues_compared"] == 22 and p["mean_plddt"] == 86.44


def test_failure(out: Path):
    bad = out / "broken.cif"
    bad.write_text("<!DOCTYPE html>")
    r = clean.clean_structure(bad, out, {"pdb_id": "XXXX", "target_entities": []}, "")
    assert not r["ok"] and "cleaning failed" in r["error"]
    assert not clean.alphafold_metrics(bad, {}, "")["ok"]


def main():
    out = Path(tempfile.mkdtemp())
    pockets = test_3og7(out)
    test_8c7y_symmetry(out)
    test_8tvk(out)
    test_alphafold(pockets)
    test_failure(out)
    print("ALL OK")


main()
