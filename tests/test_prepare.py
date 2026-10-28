"""Offline test for prepare_structures (protein_mcp/files/): chain parsing, content checks,
manifest and failure handling. Fixtures in tests/fixtures/prep/ are real replies (2026-10-22):
RCSB chains GraphQL for 3OG7/5JSM/7ZR0, a ModelServer SDF and the HTTP-200 SDF with no molecule
that ModelServer returns for an unknown label_asym_id. The structures dossier is the real
BRAF V600E one, cut to 3OG7, 5JSM and 7ZR0.
Run: uv run python tests/test_prepare.py
"""
import asyncio
import gzip
import json
import tempfile
from pathlib import Path

import httpx

from protein_mcp import dossier
from protein_mcp.files import choice, clean, download, prepare, viewer
from protein_mcp.sources import rcsb

FIX = Path(__file__).parent / "fixtures" / "prep"
CHAINS = json.loads((FIX / "rcsb_chains_3OG7_5JSM_7ZR0.json").read_text(encoding="utf-8"))
SDF = (FIX / "modelserver_3OG7_C.sdf").read_bytes()
EMPTY_SDF = (FIX / "modelserver_3OG7_ZZ_empty.sdf").read_bytes()
REAL_CIF = {"3OG7": gzip.decompress((FIX / "3OG7.cif.gz").read_bytes())}     # 5JSM is a stub: cleaning fails
REAL_AF = gzip.decompress((FIX / "AF-P15056-F1-model_v6.cif.gz").read_bytes())
FASTA = (FIX / "uniprot_P15056.fasta").read_bytes()                            # real UniProt reply

# get_structures result for BRAF V600E (real values, cut to three entries)
STRUCTURES = {
    "tool": "get_structures", "params": {"max_models": 7, "mutation": "V600E"}, "retrieved_at": "2026-10-21T21:18:49Z",
    "result": {"ok": True, "uniprot_id": "P15056",
               "experimental": {"ok": True, "total_entries": 38, "mutation_filter": "V600E", "returned": 3, "structures": [
                   {"pdb_id": "5JSM", "resolution_angstrom": 2.19,
                    "title": "BRAFV600E Kinase Domain In Complex with Chemically Linked Vemurafenib Inhibitor VEM-3-VEM",
                    "mutations": ["A688R", "F667E", "I543A", "I544S", "I551K", "K630S", "K723G", "L588N", "L706S",
                                  "L716E", "P722S", "Q562R", "Q709R", "S713E", "S720E", "V600E", "Y673S"],
                    "mutations_from_title": [], "mutations_unparsed": [], "links": {"rcsb": "https://www.rcsb.org/structure/5JSM"}},
                   {"pdb_id": "3OG7", "resolution_angstrom": 2.45,
                    "mutations": ["A688R", "F667E", "I544A", "I551K", "K630S", "L588N", "L706S", "L716E", "Q562R",
                                  "Q709R", "S713E", "S720E", "Y673S"],                 # V600E only in the title
                    "mutations_from_title": ["V600E"], "mutations_unparsed": [], "links": {}},
                   {"pdb_id": "7ZR0", "resolution_angstrom": 3.4, "mutations": [], "mutations_from_title": ["V600E"],
                    "mutations_unparsed": [], "links": {}},
               ]},
               "alphafold": {"ok": True, "entry_id": "AF-P15056-F1", "mean_plddt": 66.38, "is_canonical_sequence": True,
                             "links": {"cif": "https://alphafold.ebi.ac.uk/files/AF-P15056-F1-model_v6.cif"}}},
}


download.RETRY_WAIT = 0.0
CALLS = {}      # url -> requests seen, for the retry checks


def handler(request: httpx.Request) -> httpx.Response:
    url = str(request.url)
    CALLS[url] = CALLS.get(url, 0) + 1
    if request.method == "POST" and url == rcsb.RCSB_GRAPHQL:
        return httpx.Response(200, json=CHAINS)
    if url.startswith("https://files.rcsb.org/download/"):
        pdb_id = url.rsplit("/", 1)[1].removesuffix(".cif")
        if pdb_id == "7ZR0":                        # simulated failure
            return httpx.Response(404, text="<!DOCTYPE html><html></html>")
        return httpx.Response(200, content=REAL_CIF.get(pdb_id) or f"data_{pdb_id}\n# \n_entry.id   {pdb_id}\n".encode())
    if url.startswith("https://alphafold.ebi.ac.uk/files/"):
        return httpx.Response(200, content=REAL_AF)
    if url.endswith("/P15056.fasta"):
        return httpx.Response(200, content=FASTA)
    if url.startswith("https://models.rcsb.org/v1/"):
        if "3OG7" in url and CALLS[url] == 1:       # real failure of 6V34 QOP C, 2026-10-26: retried
            raise httpx.RemoteProtocolError("Server disconnected without sending a response.", request=request)
        if "5JSM" in url and "label_asym_id=P" in url:
            return httpx.Response(503, text="Service Unavailable")      # still down after the retries
        # 5JSM BEN instance H: the real empty reply (HTTP 200, no molecule)
        return httpx.Response(200, content=EMPTY_SDF if "label_asym_id=H" in url else SDF)
    return httpx.Response(500, text=f"unexpected {url}")


def test_parse_chains():
    entries = {e["rcsb_id"]: rcsb._parse_chains(e, "P15056") for e in CHAINS["data"]["entries"]}
    z = entries["7ZR0"]
    assert z["target_entities"][0]["chains"] == [{"label_asym_id": "D", "auth_asym_id": "K"}]   # label != auth
    assert z["target_entities"][0]["uniprot_alignment"] == [{"entity_begin": 63, "uniprot_begin": 1, "length": 766}]
    assert {o["description"] for o in z["other_polymer_entities"]} == {"Hsp90 co-chaperone Cdc37", "Heat shock protein HSP 90-beta"}
    atp = z["nonpolymer_instances"]
    assert [i["category"] for i in atp] == ["nucleotide", "nucleotide"]
    assert all(i["neighbor_target_chains"] == {} and i["neighbor_other_chains"] for i in atp)    # ATP sits on HSP90 only
    o = entries["3OG7"]
    assert [c["auth_asym_id"] for c in o["target_entities"][0]["chains"]] == ["A", "B"]
    (lig,) = o["nonpolymer_instances"]                                                           # one 032, chain A only
    assert (lig["comp_id"], lig["label_asym_id"], lig["auth_asym_id"]) == ("032", "C", "A")
    assert list(lig["neighbor_target_chains"]) == ["A"]
    j = entries["5JSM"]
    six = [i for i in j["nonpolymer_instances"] if i["comp_id"] == "6NB"]
    assert [i["label_asym_id"] for i in six] == ["J", "R"]                                       # one SDF per instance
    assert set(six[0]["neighbor_target_chains"]) == {"A", "B"}
    assert {i["comp_id"] for i in j["nonpolymer_instances"] if i["category"] == "additive"} == {"CL", "DMS", "GOL"}
    assert rcsb._parse_chains({"rcsb_id": "X", "polymer_entities": [{"rcsb_polymer_entity_container_identifiers":
                                                                     {"entity_id": "1", "uniprot_ids": None}}]},
                              "P15056")["target_entities"] == []                                # Fab chains: uniprot_ids None


def test_checks():
    assert download.check_sdf(SDF) is None
    assert "without a molecule" in download.check_sdf(EMPTY_SDF)
    assert download.check_mmcif("3OG7")(b"data_3OG7\n# \n") is None
    assert "expected" in download.check_mmcif("3OG7")(b"data_4XV1\n")
    assert "not an mmCIF" in download.check_mmcif("3OG7")(b"<!DOCTYPE html>")
    assert download.check_fasta("P15056")(FASTA) is None
    assert download.check_fasta("P15056")(b">sp|Q5IBP5|X\nMA") is not None


def network_down(request: httpx.Request) -> httpx.Response:
    raise httpx.ConnectError("down", request=request)


def cif_down(request: httpx.Request) -> httpx.Response:
    if request.url.host in ("files.rcsb.org", "alphafold.ebi.ac.uk"):
        return httpx.Response(503, text="Service Unavailable")
    return handler(request)


async def run_prepare(transport=handler, *args):
    original = httpx.AsyncClient
    httpx.AsyncClient = lambda **kw: original(transport=httpx.MockTransport(transport), **kw)
    try:
        return await prepare.prepare_structures("P15056", "V600E", *args)
    finally:
        httpx.AsyncClient = original


def write_dossier(data):
    (dossier.DOSSIER_DIR / "P15056" / "structures_v600e.json").write_text(json.dumps(data), encoding="utf-8")


def test_prepare():
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())
    missing = asyncio.run(prepare.prepare_structures("P15056", "V600E"))
    assert not missing["ok"] and 'mutation="V600E"' in missing["error"]

    (dossier.DOSSIER_DIR / "P15056").mkdir()
    write_dossier(STRUCTURES)
    reply = asyncio.run(run_prepare())
    assert reply["ok"] and reply["selection"].startswith("3 of 3 experimental structures in structures_v600e.json + the AlphaFold")
    base = dossier.DOSSIER_DIR / "P15056" / "prep"
    m = json.loads((base / "manifest.json").read_text(encoding="utf-8"))
    assert m["source_dossier"]["file"] == "structures_v600e.json" and m["user_choice"] is None and m["pending"]
    assert m["sequence"]["header"].startswith(">sp|P15056|") and m["sequence"]["length"] == 766

    by_id = {s["pdb_id"]: s for s in m["structures"]}
    assert by_id["3OG7"]["file"]["ok"] and (base / by_id["3OG7"]["file"]["path"]).is_file()
    assert by_id["3OG7"]["construct"]["mutations_from_title"] == ["V600E"]
    z = by_id["7ZR0"]
    assert not z["file"]["ok"] and z["file"]["error"] == "HTTP 404" and not (base / "raw" / "7ZR0.cif").exists()
    assert z["target_entities"] is not None                    # chain data still recorded when the file fails
    ben = next(i for i in by_id["5JSM"]["nonpolymer_instances"] if i["label_asym_id"] == "H")
    assert not ben["sdf"]["ok"] and "without a molecule" in ben["sdf"]["error"]
    assert not (base / "ligands" / "5JSM_BEN_H.sdf").exists()
    assert all(i["sdf"] is None for i in by_id["5JSM"]["nonpolymer_instances"] if i["category"] == "additive")
    assert not list(base.rglob("*.part"))
    assert len(reply["failed_downloads"]) == 3 and any("7ZR0" in f for f in reply["failed_downloads"])
    assert m["failed_downloads"] == reply["failed_downloads"]             # also at the top of the manifest
    # retries: a cut connection is retried and recovers; a 503 gives up after 3 attempts; a rejected
    # content (HTTP 200 without a molecule) and a 404 are not retried
    o_sdf = next(i for i in by_id["3OG7"]["nonpolymer_instances"] if i["comp_id"] == "032")["sdf"]
    assert o_sdf["ok"] and o_sdf["attempts"] == 2
    ben_p = next(i for i in by_id["5JSM"]["nonpolymer_instances"] if i["label_asym_id"] == "P")["sdf"]
    assert not ben_p["ok"] and ben_p["error"] == "HTTP 503 (after 3 attempts)"
    assert CALLS[rcsb.ligand_sdf_url("5JSM", "P")] == 3 and CALLS[rcsb.ligand_sdf_url("5JSM", "H")] == 1
    assert CALLS[rcsb.cif_url("7ZR0")] == 1 and "attempts" not in by_id["3OG7"]["file"]

    cz = next(s for s in reply["structures"] if s["pdb_id"] == "7ZR0")
    assert cz["download_error"] == "HTTP 404" and cz["sample_uniprot_ranges"] == ["1-766"]
    c3 = next(s for s in reply["structures"] if s["pdb_id"] == "3OG7")
    assert len(c3["declared_mutations"]) == 13 and "V600E" not in c3["declared_mutations"]
    (vem,) = c3["ligand_instances"]                     # wwPDB synonyms only (no DrugBank ones)
    assert vem["synonyms"] == ["Vemurafenib", "PLX4032"] and vem["formula_weight"] == 489.922
    assert vem["name"].startswith("N-(3-{[5-(4-chlorophenyl)") and len(vem["name"]) <= 100
    c5 = next(s for s in reply["structures"] if s["pdb_id"] == "5JSM")
    assert c5["title"].endswith("VEM-3-VEM")              # the only hint that 6NB is a linked vemurafenib
    six = next(x for x in c5["ligand_instances"] if x["comp_id"] == "6NB")
    assert six["synonyms"] == [] and six["formula_weight"] == 1057.095
    assert cz["declared_mutations"] == [] and "receptor" not in cz                 # not downloaded, not cleaned
    # step 2 (cleaning) on the real 3OG7; the 5JSM stub fails and is reported, not hidden
    assert len(reply["warnings"]) == 1 and reply["warnings"][0].startswith("5JSM: cleaning failed: the file has no model"), reply["warnings"]
    assert next(s for s in reply["structures"] if s["pdb_id"] == "5JSM")["clean_error"].startswith("cleaning failed")
    assert c3["receptor"] == "clean/3OG7_receptor.pdb" and (base / c3["receptor"]).is_file()
    assert [(c["auth_asym_id"], c["variant"]) for c in c3["chain_checks"]] == [("A", "not modelled"), ("B", "variant")]
    (pk,) = c3["pockets"]
    assert pk["n_residues"] == 22 and "A:K483" in pk["residues"] and "A:G596" in pk["residues"]
    sm = reply["summary"]                                 # same summary in the reply and the manifest
    assert sm == m["summary"] and sm["structures_not_compared"] == ["5JSM", "7ZR0"]
    assert sm["pockets_with_variant_in_pocket_chain"] == [] and sm["pockets_without_variant_in_pocket_chain"] == ["3OG7 032 C: chain A not modelled (22 pocket residues); single-chain pocket"]
    assert sm["by_resolution"] == ["5JSM 2.19 (not compared)", "3OG7 2.45", "7ZR0 3.4 (not compared)"]
    assert pk["residues_by_chain"] == {"A": 22} and pk["heavy_atoms"] == 33 and "spans_target_chains" not in pk
    assert pk["adjacent_missing_ranges"] == [{"auth_asym_id": "A", "missing_uniprot_range": [597, 614]}]
    assert pk["alphafold_mean_plddt"] == 86.44 and reply["alphafold"]["variant_position"]["residue"] == "VAL"
    assert m["structures"][1]["prepared"]["pockets"][0]["residues"][0]["uniprot"] == 463   # full detail in the manifest
    assert not any("step 2" in x or "step 3" in x for x in m["pending"])
    assert [l["neighbor_target_auth_chains"] for l in cz["ligand_instances"]] == [[], []]
    assert reply["alphafold"]["ok"] and m["alphafold"]["file"]["path"] == "raw/AF-P15056-F1-model_v6.cif"

    # open_in_pymol on this manifest (step 3): only 3OG7 was cleaned (5JSM is a stub, 7ZR0 a 404)
    view = asyncio.run(viewer.open_in_pymol("P15056", start=False))
    assert view["ok"] and [s["pdb_id"] for s in view["structures"]] == ["3OG7"]
    assert {s["pdb_id"] for s in view["skipped"]} == {"5JSM", "7ZR0"}
    (v3,) = view["structures"]
    assert v3["variant_position"]["chains"] == ["B"] and [c["status"] for c in v3["variant_position"]["by_chain"]] == [
        "not modelled", "variant"]
    assert [(p["selection"], p["n_residues"], p["box_object"]) for p in v3["pockets"]] == [
        ("pocket_3OG7_032_C", 22, "box_3OG7_032_C")]
    pml = Path(view["pml_path"]).read_text(encoding="ascii")
    assert "chain A and resi 463+" in pml and "+596))" in pml and "AF-P15056-F1-model_v6.cif" in pml

    # a second run says what it replaced; IDs outside the dossier are never downloaded
    again = asyncio.run(prepare.prepare_structures("P15056", "V600E", ["1UWH"], include_alphafold=False))
    assert not again["ok"] and again["pdb_ids_not_in_dossier"] == ["1UWH"]
    assert reply["user_choice"] is None
    assert choice.select_structure("P15056", "3OG7", "032 C", reason="pose de vemurafenib")["ok"]
    reply2 = asyncio.run(run_prepare())
    assert reply2["replaced_manifest_from"]["file"] == "structures_v600e.json"
    # step 4: the choice survives a rebuild with the same files, and is dropped (kept in history) without 3OG7
    assert reply2["user_choice"] == {"kept": "3OG7 chain A, pocket 032 C", "changed_fields": []}, reply2["user_choice"]
    kept = json.loads((base / "manifest.json").read_text(encoding="utf-8"))
    assert kept["user_choice"]["reason"] == "pose de vemurafenib" and choice.STEP4_PENDING not in kept["pending"]
    sub = asyncio.run(run_prepare(handler, ["5JSM"]))
    assert sub["user_choice"]["dropped"] == "3OG7 chain A, pocket 032 C" and "not in the manifest" in sub["user_choice"]["reason"]
    dropped = json.loads((base / "manifest.json").read_text(encoding="utf-8"))
    assert dropped["user_choice"] is None and dropped["user_choice_history"][0]["reason"] == "pose de vemurafenib"
    reply2 = asyncio.run(run_prepare())
    assert reply2["user_choice"] is None and len(json.loads((base / "manifest.json").read_text(encoding="utf-8"))["user_choice_history"]) == 1
    assert choice.select_structure("P15056", "3OG7", "032 C", reason="otra vez")["ok"]

    # network down / every structure file failing: ok false, the good manifest is kept
    good = (base / "manifest.json").read_text(encoding="utf-8")
    down = asyncio.run(run_prepare(network_down))
    assert not down["ok"] and "kept" in down["error"]
    failed = asyncio.run(run_prepare(cif_down))
    assert not failed["ok"] and "Every structure download failed" in failed["error"] and failed["failed_downloads"]
    assert (base / "manifest.json").read_text(encoding="utf-8") == good
    assert json.loads(good)["user_choice"]["reason"] == "otra vez"       # a failed rerun keeps the choice too

    # RCSB failed when the dossier was made (get_structures still saved it because AlphaFold worked)
    broken = json.loads(json.dumps(STRUCTURES))
    broken["result"]["experimental"] = {"ok": False, "error": "RCSB search failed: timeout"}
    write_dossier(broken)
    r = asyncio.run(run_prepare())
    assert r["ok"] and "unknown, not zero" in r["warnings"][0] and "None" not in r["selection"]
    assert r["selection"].startswith("0 of 0") and "RCSB failed" in r["selection"]

    # AlphaFold requested but failed in the dossier: said, not silently dropped
    broken = json.loads(json.dumps(STRUCTURES))
    broken["result"]["alphafold"] = {"ok": False, "error": "AlphaFold DB request failed: 500"}
    write_dossier(broken)
    r = asyncio.run(run_prepare())
    assert r["ok"] and r["alphafold"] is None and any("AlphaFold not used" in w for w in r["warnings"])
    r = asyncio.run(run_prepare(handler, None, False))           # not requested: no warning
    assert r["ok"] and not any("AlphaFold" in w for w in r["warnings"])


def test_summary():
    """Real 8C7Y, 5JSM and 3OG7 files: the cases where Claude generalised in Desktop (2026-10-26)."""
    canonical = "".join(l.strip() for l in FASTA.decode().splitlines()[1:])
    chains = {e["rcsb_id"]: rcsb._parse_chains(e, "P15056") for e in CHAINS["data"]["entries"]}
    chains.update({e["rcsb_id"]: rcsb._parse_chains(e, "P15056") for e in
                   json.loads((FIX / "rcsb_chains_8C7Y.json").read_text(encoding="utf-8"))["data"]["entries"]})
    out, structures = Path(tempfile.mkdtemp()), []
    for pdb_id, res in (("3OG7", 2.45), ("5JSM", 2.19), ("8C7Y", 1.65), ("7ZR0", 3.4)):
        rec = dict(chains[pdb_id], pdb_id=pdb_id, resolution_angstrom=res)
        rec["prepared"] = (clean.clean_structure(FIX / f"{pdb_id}.cif.gz", out, rec, "V600E", canonical)
                           if pdb_id != "7ZR0" else {"ok": False, "error": "download failed"})
        structures.append(rec)
    s = prepare.summarize(structures)
    assert s["structures_compared"] == ["3OG7", "5JSM", "8C7Y"] and s["structures_not_compared"] == ["7ZR0"]
    assert s["by_resolution"] == ["8C7Y 1.65", "5JSM 2.19", "3OG7 2.45", "7ZR0 3.4 (not compared)"]
    assert s["scope"].startswith("3 of 4 structures compared")
    assert s["by_r_free"] == ["8C7Y 0.175", "3OG7 0.2575", "5JSM 0.2597"]
    assert s["chains_without_internal_gaps"] == ["8C7Y:A (modelled 448-721)"]   # 597-614 is NOT everywhere; 444-447 unmodelled
    assert s["ligands_without_pocket"] == ["5JSM BEN H: fewer than 10 heavy atoms", "5JSM EOH I: fewer than 10 heavy atoms",
                                           "5JSM EOH M: fewer than 10 heavy atoms", "5JSM EOH N: fewer than 10 heavy atoms",
                                           "5JSM BEN P: fewer than 10 heavy atoms"]   # too small: listed, not dropped
    st = s["variant_position_by_chain_status"]
    assert st["variant"] == ["3OG7:B", "5JSM:A", "5JSM:B", "5JSM:C", "5JSM:D", "8C7Y:A"]
    assert st["not modelled"] == ["3OG7:A", "8C7Y:B"]
    # one text format for all lists (Desktop test 3: Claude misread the former list of objects)
    assert s["pockets_with_variant_in_pocket_chain"] == [
        "5JSM 6NB J: chain A variant (24 pocket residues), chain B variant (23 pocket residues); spans several chains",
        "5JSM 6NB R: chain C variant (24 pocket residues), chain D variant (23 pocket residues); spans several chains",
        "8C7Y TXV U: chain A variant (26 pocket residues); single-chain pocket"]
    assert s["pockets_without_variant_in_pocket_chain"] == [
        "3OG7 032 C: chain A not modelled (22 pocket residues); single-chain pocket", "8C7Y TXV X: chain B not modelled (28 pocket residues); single-chain pocket"]
    assert s["pocket_counts"] == {"with_variant_in_pocket_chain": 3, "ambiguous": 0,
                                  "without_variant_in_pocket_chain": 2, "not_classified": 0}
    assert s["pockets_with_ambiguous_variant_position"] == [] and s["pockets_not_classified"] == []

    # reviewer edge cases (2026-10-26): nothing is dropped silently and "ambiguous" is not "without"
    edge = json.loads(json.dumps(structures[:3]))
    o, j, y = edge
    o["prepared"]["metrics"]["chains"][0]["variant_check"]["status"] = "ambiguous (microheterogeneity)"
    y["prepared"]["pockets"][0] = {"comp_id": "TXV", "label_asym_id": "U", "error": "not found in the file"}
    j["prepared"]["pockets"][0]["residues_by_chain"] = {"A": 24, "C": 1}         # variant chain + other
    j["prepared"]["metrics"]["chains"][0]["variant_check"]["status"] = "not modelled"
    e = prepare.summarize(edge)
    assert e["pockets_with_ambiguous_variant_position"] == [
        "3OG7 032 C: chain A ambiguous (microheterogeneity) (22 pocket residues); single-chain pocket"]
    assert e["pockets_not_classified"] == ["8C7Y TXV U: not found in the file"]
    assert "5JSM 6NB J: chain A not modelled (24 pocket residues), chain C variant (1 pocket residue); spans several chains"         in e["pockets_with_variant_in_pocket_chain"]             # the variant chain holds 1 residue: visible
    # no substitution filter: no variant fields, the rest is still there
    for r in structures[:3]:
        r["prepared"] = clean.clean_structure(FIX / f"{r['pdb_id']}.cif.gz", out, r, "", canonical)
    s = prepare.summarize(structures)
    assert "variant_position_by_chain_status" not in s and s["chains_without_internal_gaps"] == ["8C7Y:A (modelled 448-721)"]
    empty = prepare.summarize([])
    assert empty["scope"].startswith("0 of 0") and empty["by_resolution"] == [] and "variant_position_by_chain_status" not in empty


def main():
    test_parse_chains()
    test_checks()
    test_prepare()
    test_summary()
    print("ALL OK")


main()
