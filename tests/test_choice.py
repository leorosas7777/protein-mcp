"""Offline test for select_structure (protein_mcp/files/choice.py), its carry-over when
prepare_structures replaces the manifest, and report section 6 "Structure preparation".
The manifest is built like prepare_structures builds it, from REAL files: RCSB chain replies
for 3OG7/5JSM/7ZR0 and 8C7Y, the real mmCIFs cleaned with clean.clean_structure and the real
UniProt FASTA (tests/fixtures/prep/). One SDF failure has the real shape of 6V34 QOP C.
Run: uv run python tests/test_choice.py
"""
import copy
import json
import re
import tempfile
from pathlib import Path

from protein_mcp import dossier, report
from protein_mcp.files import choice, clean, prepare
from protein_mcp.sources import rcsb

FIX = Path(__file__).parent / "fixtures" / "prep"
ACC = "P15056"


def build_manifest(out: Path) -> dict:
    canonical = "".join(l.strip() for l in (FIX / "uniprot_P15056.fasta").read_text().splitlines()[1:])
    chains = {}
    for name in ("rcsb_chains_3OG7_5JSM_7ZR0.json", "rcsb_chains_8C7Y.json"):
        for e in json.loads((FIX / name).read_text(encoding="utf-8"))["data"]["entries"]:
            chains[e["rcsb_id"]] = rcsb._parse_chains(e, ACC)
    structures = []
    for pdb_id, res in (("8C7Y", 1.65), ("5JSM", 2.19), ("3OG7", 2.45), ("7ZR0", 3.4)):
        rec = {"pdb_id": pdb_id, "resolution_angstrom": res,
               "links": {"rcsb": f"https://www.rcsb.org/structure/{pdb_id}"},
               "file": {"ok": True, "path": f"raw/{pdb_id}.cif", "sha256": f"sha-{pdb_id}"},
               "construct": {"mutations": ["I544A", "F667E"] if pdb_id == "3OG7" else [],
                             # 8C7Y: variant named only in the title (real dossier); empty field != native
                             "mutations_from_title": ["V600E"] if pdb_id == "8C7Y" else [],
                             "mutations_unparsed": []}}
        rec.update(copy.deepcopy(chains[pdb_id]))
        for inst in rec["nonpolymer_instances"]:
            inst["sdf"] = None if inst["category"] == "additive" else {
                "ok": True, "path": f"ligands/{pdb_id}_{inst['comp_id']}_{inst['label_asym_id']}.sdf"}
        if pdb_id == "7ZR0":       # real 404 of 2026-10-22 test: not downloaded, not cleaned
            rec["file"] = {"ok": False, "url": rcsb.cif_url("7ZR0"), "error": "HTTP 404"}
            rec["prepared"] = None
        else:
            rec["prepared"] = clean.clean_structure(FIX / f"{pdb_id}.cif.gz", out, rec, "V600E", canonical)
            rec["prepared"]["receptor"]["path"] = f"clean/{pdb_id}_receptor.pdb"
        structures.append(rec)
    # AlphaFold pLDDT per pocket, from the real AF model, as prepare_structures attaches it
    pockets = {r["pdb_id"]: r["prepared"]["pockets"] for r in structures if (r["prepared"] or {}).get("ok")}
    af = clean.alphafold_metrics(FIX / "AF-P15056-F1-model_v6.cif.gz", pockets, "V600E", canonical)
    hits = {(x["pdb_id"], x["label_asym_id"]): x for x in af["pockets"]}
    for pdb_id, plist in pockets.items():
        for pk in plist:
            hit = hits.get((pdb_id, pk["label_asym_id"]))
            if hit:
                pk["alphafold_plddt"] = {k: hit[k] for k in ("mean_plddt", "residues_below_70", "residues_compared")}
    # real failure shape (6V34 QOP C, 2026-10-26), put on 8C7Y TXV X
    txv_x = next(i for i in structures[0]["nonpolymer_instances"] if i["label_asym_id"] == "X")
    txv_x["sdf"] = {"ok": False, "url": rcsb.ligand_sdf_url("8C7Y", "X"),
                    "error": "RemoteProtocolError('Server disconnected without sending a response.')"}
    return {
        "manifest_version": 1, "uniprot_id": ACC, "created_at": "2026-10-28T10:00:00Z",
        "generated_by": "prepare_structures",
        "source_dossier": {"file": "structures_v600e.json", "retrieved_at": "2026-10-21T21:18:49Z"},
        "mutation_filter": "V600E", "selection": "4 of 7 experimental structures in structures_v600e.json",
        "summary": prepare.summarize(structures), "structures": structures,
        "alphafold": {"entry_id": "AF-P15056-F1", "is_canonical_sequence": True, "mean_plddt": 66.38,
                      "file": {"ok": True, "path": "raw/AF-P15056-F1-model_v6.cif", "sha256": "sha-af",
                               "url": "https://alphafold.ebi.ac.uk/files/AF-P15056-F1-model_v6.cif"},
                      "variant_position": {"uniprot_position": 600, "residue": "VAL", "plddt": 49.06}},
        "failed_downloads": ["cif/7ZR0: HTTP 404"], "warnings": [], "decisions": prepare.DECISIONS,
        "pending": list(prepare.PENDING), "user_choice": None,
    }


def write(m):
    choice.write_manifest(choice.manifest_path(ACC), m)


def read():
    return json.loads(choice.manifest_path(ACC).read_text(encoding="utf-8"))


def test_select(base: dict):
    sel = choice.select_structure
    missing = sel(ACC, "3OG7", reason="x")
    assert not missing["ok"] and "run prepare_structures first" in missing["error"]
    write(base)
    before = choice.manifest_path(ACC).read_text(encoding="utf-8")

    # rejected: nothing written, options given
    bad = {
        ("1ABC", "", None): "not in the manifest",
        ("7ZR0", "", ("K",)): "no cleaned receptor",                       # 404 in the manifest
        ("8C7Y", "", None): "has 2 target chains",
        ("8C7Y", "TXV", None): "matches 2 pockets",
        ("8C7Y", "TXV U", ("B",)): "hold no residue of pocket TXV U",
        ("8C7Y", "TXV U", ("a",)): "not a target chain",                    # auth ids are case-sensitive
        ("3OG7", "BEN", None): "no pocket 'BEN'",
        ("5JSM", "BEN H", None): "no pocket",                              # too small: no pocket in the manifest
        ("AlphaFold", "032 C", None): "no ligand",
        ("3OG7", "032 C D", None): "is not '<comp_id> <label_asym_id>'",
    }
    for (structure, pocket, chains), text in bad.items():
        r = sel(ACC, structure, pocket, chains, reason="x")
        assert not r["ok"] and text in r["error"], (structure, pocket, r)
        assert set(r["options"]["structures"]) == {"8C7Y", "5JSM", "3OG7", "AlphaFold"}
        assert r["options"]["not_selectable"] == {"7ZR0": "HTTP 404"}
    # Desktop test 2026-10-28: the choice is checked before the reason, so a bad choice is never hidden
    assert "has 2 target chains" in sel(ACC, "8C7Y", reason="")["error"]
    nr = sel(ACC, "3OG7", "032", reason="  ")
    assert not nr["ok"] and "reason is required" in nr["error"] and nr["valid_choice"] == "3OG7 chain A, pocket 032 C"
    assert "longer than" in sel(ACC, "3OG7", "032", reason="x" * 1001)["error"]
    assert not sel(ACC, "", reason="x")["ok"]
    assert choice.manifest_path(ACC).read_text(encoding="utf-8") == before

    # accepted: 032 alone is unique in 3OG7; the summary format "<PDB> comp label" works too
    r = sel(ACC, "3og7", "032", reason="Pose cristalográfica de vemurafenib")
    assert r["ok"] and r["recorded"] == "3OG7 chain A, pocket 032 C" and r["replaced_choice"] is None
    uc = read()["user_choice"]
    assert uc == r["user_choice"] and uc["reason"] == "Pose cristalográfica de vemurafenib"
    assert uc["receptor"] == "clean/3OG7_receptor.pdb" and uc["raw_sha256"] == "sha-3OG7"
    assert uc["chains"] == [{"auth_asym_id": "A", "label_asym_id": "A", "modelled_uniprot_range": [449, 720],
                             "missing_uniprot_ranges": [[545, 547], [597, 614], [627, 630]],
                             "variant_status": "not modelled", "variant_observed": None}]
    pk = uc["pocket"]
    assert (pk["n_residues"], pk["residues_by_chain"], pk["synonyms"]) == (22, {"A": 22}, ["Vemurafenib", "PLX4032"])
    assert pk["sdf"] == "ligands/3OG7_032_C.sdf" and pk["adjacent_missing_ranges"] and uc["warnings"] == []
    assert uc["construct"] == {"mutations": ["I544A", "F667E"], "mutations_from_title": [], "mutations_unparsed": []}
    assert uc["receptor_atoms"] > 0 and uc["receptor_kept_nonpolymer"] == []
    assert choice.STEP4_PENDING not in read()["pending"]

    # choosing again replaces it; history keeps the old one
    r = sel(ACC, "8C7Y", "8C7Y TXV X", reason="cadena B")
    assert r["ok"] and r["replaced_choice"] == "3OG7 chain A, pocket 032 C" and r["previous_choices"] == 1
    assert any("no SDF" in w and "RemoteProtocolError" in w for w in r["user_choice"]["warnings"])
    assert r["user_choice"]["pocket"]["sdf"] is None
    # bivalent pocket with one of its chains: accepted, box not reduced (said)
    r = sel(ACC, "5JSM", "6NB J", ["A"], reason="sitio A")
    assert r["ok"] and r["user_choice"]["pocket"]["residues_by_chain"] == {"A": 24, "B": 23}
    assert any("NOT reduced" in w for w in r["user_choice"]["warnings"])
    assert r["recorded"] == "5JSM chain A, pocket 6NB J (pocket spans chains A, B; box not reduced)"
    assert choice.label({**r["user_choice"], "chains": [{"auth_asym_id": "A"}, {"auth_asym_id": "B"}]}) == "5JSM chains A, B, pocket 6NB J"
    # chain without pocket: allowed, said
    r = sel(ACC, "8C7Y", chains=["A"], reason="apo")
    assert r["ok"] and r["user_choice"]["pocket"] is None and "no pocket" in r["recorded"]
    assert any("No pocket chosen" in w for w in r["user_choice"]["warnings"])
    r = sel(ACC, "alphafold", reason="modelo completo")
    assert r["ok"] and r["user_choice"]["kind"] == "alphafold" and r["user_choice"]["file_sha256"] == "sha-af"
    assert len(read()["user_choice_history"]) == 4


def test_carry_over(base: dict):
    write(base)
    choice.select_structure(ACC, "3OG7", "032 C", reason="pose")
    old = read()

    same = copy.deepcopy(base)
    same["created_at"] = "2026-10-29T00:00:00Z"
    status = choice.carry_over(old, same)
    assert status == {"kept": "3OG7 chain A, pocket 032 C", "changed_fields": []}
    assert same["user_choice"]["reason"] == "pose" and same["user_choice"]["chosen_at"] == old["user_choice"]["chosen_at"]
    assert same["user_choice"]["revalidated_at"] == "2026-10-29T00:00:00Z"
    assert choice.STEP4_PENDING not in same["pending"]

    moved = copy.deepcopy(base)               # e.g. another cutoff: the pocket data change, the choice stays
    pocket = next(s for s in moved["structures"] if s["pdb_id"] == "3OG7")["prepared"]["pockets"][0]
    pocket["box"] = dict(pocket["box"], size=[1.0, 1.0, 1.0])
    assert choice.carry_over(old, moved)["changed_fields"] == ["pocket"]

    subset = copy.deepcopy(base)              # prepare_structures(pdb_ids=["8C7Y"])
    subset["structures"] = [s for s in subset["structures"] if s["pdb_id"] != "3OG7"]
    status = choice.carry_over(old, subset)
    assert status["dropped"] == "3OG7 chain A, pocket 032 C" and "not in the manifest" in status["reason"]
    assert subset["user_choice"] is None and subset["user_choice_history"][-1]["dropped_reason"] == status["reason"]
    assert subset["pending"][0] == choice.STEP4_PENDING

    newer = copy.deepcopy(base)               # RCSB re-released the entry
    next(s for s in newer["structures"] if s["pdb_id"] == "3OG7")["file"]["sha256"] = "other"
    assert "SHA-256" in choice.carry_over(old, newer)["reason"]

    fresh = copy.deepcopy(base)
    assert choice.carry_over(None, fresh) is None and fresh["user_choice"] is None and fresh["user_choice_history"] == []


def test_report(base: dict):
    folder = dossier.DOSSIER_DIR / ACC
    (folder / "structures_v600e.json").write_text(json.dumps({
        "tool": "get_structures", "params": {"mutation": "V600E"}, "retrieved_at": "2026-10-21T21:18:49Z",
        "result": {"ok": True, "experimental": {"ok": True, "structures": []}, "alphafold": {"ok": False}}}),
        encoding="utf-8")
    write(base)
    r = report.build_report(ACC)
    assert r["ok"] and r["preparation"] == "prep/manifest.json", r
    md = Path(r["markdown_path"]).read_text(encoding="utf-8")
    sec = md.split("## 6. Structure preparation (Phase 2)")[1].split("## 7. Methods")[0]
    assert "No structure chosen (select_structure not called)." in sec
    assert "4 of 7 experimental structures" in sec and "Note:" not in sec
    # 5 pockets of 10 ligand instances in 3 prepared entries; 7ZR0 ATP (nucleotide) is not counted
    assert "5 pockets computed of 10 ligand instances in the 3 of 4 entries that were prepared" in sec, sec
    assert "| 7ZR0 |" not in sec and "[7ZR0](https://www.rcsb.org/structure/7ZR0)" in sec and "not prepared: HTTP 404" in sec
    assert "A: modelled 448-721; no internal gap; variant position: variant (GLU 600)" in sec
    # reviewer 2026-10-28: title-only variant is shown (not "0"); pLDDT counts UniProt positions, not residues
    assert "| title: V600E |" in sec and "| field: I544A, F667E |" in sec and "none listed (not proof of wild type)" in sec
    assert "87.28 (at 24 UniProt positions)" in sec and " of 47 residues" not in sec
    assert report._chain_line({"auth_asym_id": "A", "label_asym_id": "A", "modelled_uniprot_range": None,
                               "missing_uniprot_ranges": []}) == "A: modelled range no record; gaps: no record"
    assert "spans several target chains" in sec and "missing 597-614 next to the pocket (chain A)" in sec
    assert "5JSM BEN H: fewer than 10 heavy atoms" in sec and "cif/7ZR0: HTTP 404" in sec
    assert "Structure preparation rules" in md and "## 6. Structure" in md

    choice.select_structure(ACC, "3OG7", "032", reason="Pose de vemurafenib")
    choice.select_structure(ACC, "8C7Y", "TXV U", reason="V600E en la cadena del bolsillo")
    md = Path(report.build_report(ACC)["markdown_path"]).read_text(encoding="utf-8")
    sec = md.split("### Structure chosen by the user")[1].split("## 7. Methods")[0]
    assert "- Choice: 8C7Y chain A, pocket TXV U" in sec and "- Construct mutations: title: V600E" in sec
    assert 'Reason recorded by select_structure: "V600E en la cadena del bolsillo"' in sec
    assert "1 earlier choice(s) in user_choice_history: 3OG7 chain A, pocket 032 C (replaced)" in sec
    assert "Ligand pose (SDF): ligands/8C7Y_TXV_U.sdf" in sec
    # anti-hallucination: every link in section 6 is in the manifest
    blob = choice.manifest_path(ACC).read_text(encoding="utf-8")
    sec6 = md.split("## 6. Structure preparation")[1].split("## 7. Methods")[0]
    for url in re.findall(r"\]\((https?://[^)\s]+)\)", sec6):
        assert url in blob, url

    # the dossier file was fetched again after the manifest: said
    data = json.loads((folder / "structures_v600e.json").read_text(encoding="utf-8"))
    data["retrieved_at"] = "2026-10-28T09:00:00Z"
    (folder / "structures_v600e.json").write_text(json.dumps(data), encoding="utf-8")
    md = Path(report.build_report(ACC)["markdown_path"]).read_text(encoding="utf-8")
    assert "was retrieved 2026-10-28T09:00:00Z; the manifest was built from the copy retrieved 2026-10-21T21:18:49Z" in md

    # no manifest: no section 6, methods keep their content
    choice.manifest_path(ACC).unlink()
    r = report.build_report(ACC)
    md = Path(r["markdown_path"]).read_text(encoding="utf-8")
    assert r["preparation"] is None and "## 6. Structure" not in md and "## 7. Methods" in md


def main():
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())
    base = build_manifest(Path(tempfile.mkdtemp()))
    test_select(base)
    test_carry_over(base)
    test_report(base)
    print("ALL OK")


main()
