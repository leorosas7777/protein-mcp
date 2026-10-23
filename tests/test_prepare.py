"""Offline test for prepare_structures (protein_mcp/files/): chain parsing, content checks,
manifest and failure handling. Fixtures in tests/fixtures/prep/ are real replies (2026-10-22):
RCSB chains GraphQL for 3OG7/5JSM/7ZR0, a ModelServer SDF and the HTTP-200 SDF with no molecule
that ModelServer returns for an unknown label_asym_id. The structures dossier is the real
BRAF V600E one, cut to 3OG7, 5JSM and 7ZR0.
Run: uv run python tests/test_prepare.py
"""
import asyncio
import json
import tempfile
from pathlib import Path

import httpx

from protein_mcp import dossier
from protein_mcp.files import download, prepare
from protein_mcp.sources import rcsb

FIX = Path(__file__).parent / "fixtures" / "prep"
CHAINS = json.loads((FIX / "rcsb_chains_3OG7_5JSM_7ZR0.json").read_text(encoding="utf-8"))
SDF = (FIX / "modelserver_3OG7_C.sdf").read_bytes()
EMPTY_SDF = (FIX / "modelserver_3OG7_ZZ_empty.sdf").read_bytes()
FASTA = b">sp|P15056|BRAF_HUMAN Serine/threonine-protein kinase B-raf OS=Homo sapiens\nMAALSGGGGG\nAEPGQ\n"

# get_structures result for BRAF V600E (real values, cut to three entries)
STRUCTURES = {
    "tool": "get_structures", "params": {"max_models": 7, "mutation": "V600E"}, "retrieved_at": "2026-10-21T21:18:49Z",
    "result": {"ok": True, "uniprot_id": "P15056",
               "experimental": {"ok": True, "total_entries": 38, "mutation_filter": "V600E", "returned": 3, "structures": [
                   {"pdb_id": "5JSM", "resolution_angstrom": 2.19,
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
               "alphafold": {"ok": True, "entry_id": "AF-P15056-F1", "mean_plddt": 66.38,
                             "links": {"cif": "https://alphafold.ebi.ac.uk/files/AF-P15056-F1-model_v6.cif"}}},
}


def handler(request: httpx.Request) -> httpx.Response:
    url = str(request.url)
    if request.method == "POST" and url == rcsb.RCSB_GRAPHQL:
        return httpx.Response(200, json=CHAINS)
    if url.startswith("https://files.rcsb.org/download/"):
        pdb_id = url.rsplit("/", 1)[1].removesuffix(".cif")
        if pdb_id == "7ZR0":                        # simulated failure
            return httpx.Response(404, text="<!DOCTYPE html><html></html>")
        return httpx.Response(200, content=f"data_{pdb_id}\n# \n_entry.id   {pdb_id}\n".encode())
    if url.startswith("https://alphafold.ebi.ac.uk/files/"):
        return httpx.Response(200, content=b"data_AF-P15056-F1\n#\n")
    if url.endswith("/P15056.fasta"):
        return httpx.Response(200, content=FASTA)
    if url.startswith("https://models.rcsb.org/v1/"):
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
    assert m["sequence"]["header"].startswith(">sp|P15056|") and m["sequence"]["length"] == 15

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
    assert len(reply["failed_downloads"]) == 2 and any("7ZR0" in f for f in reply["failed_downloads"])

    cz = next(s for s in reply["structures"] if s["pdb_id"] == "7ZR0")
    assert cz["download_error"] == "HTTP 404" and cz["sample_uniprot_ranges"] == ["1-766"]
    c3 = next(s for s in reply["structures"] if s["pdb_id"] == "3OG7")
    assert len(c3["declared_mutations"]) == 13 and "V600E" not in c3["declared_mutations"]
    assert cz["declared_mutations"] == [] and reply["warnings"] == []
    assert [l["neighbor_target_auth_chains"] for l in cz["ligand_instances"]] == [[], []]
    assert reply["alphafold"]["ok"] and m["alphafold"]["file"]["path"] == "raw/AF-P15056-F1-model_v6.cif"

    # a second run says what it replaced; IDs outside the dossier are never downloaded
    again = asyncio.run(prepare.prepare_structures("P15056", "V600E", ["1UWH"], include_alphafold=False))
    assert not again["ok"] and again["pdb_ids_not_in_dossier"] == ["1UWH"]
    reply2 = asyncio.run(run_prepare())
    assert reply2["replaced_manifest_from"]["file"] == "structures_v600e.json"

    # network down / every structure file failing: ok false, the good manifest is kept
    good = (base / "manifest.json").read_text(encoding="utf-8")
    down = asyncio.run(run_prepare(network_down))
    assert not down["ok"] and "kept" in down["error"]
    failed = asyncio.run(run_prepare(cif_down))
    assert not failed["ok"] and "Every structure download failed" in failed["error"] and failed["failed_downloads"]
    assert (base / "manifest.json").read_text(encoding="utf-8") == good

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


def main():
    test_parse_chains()
    test_checks()
    test_prepare()
    print("ALL OK")


main()
