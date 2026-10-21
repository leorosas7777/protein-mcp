"""Offline test for target_profile.py (UniProt + Open Targets) using REAL responses
saved on 2026-10-20 in tests/fixtures/: KRAS entry (DISEASE comments without a
`disease` block, versioned Ensembl gene IDs, the G10 -> GG insertion), Open Targets
KRAS (germline syndromes first by score), the null answer to a versioned ID and an
inactive UniProt entry that answers HTTP 200.
Run: uv run python tests/test_profile.py
"""
import asyncio
import copy
import json
import sys
import tempfile
from pathlib import Path

import httpx

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))  # project root
import common  # noqa: E402
import target_profile  # noqa: E402

FIX = ROOT / "tests" / "fixtures"


def load(name):
    return json.loads((FIX / name).read_text(encoding="utf-8"))


KRAS = load("uniprot_P01116_profile.json")
OT_KRAS = load("opentargets_ENSG00000133703.json")
OT_NULL = load("opentargets_versioned_id_null.json")
INACTIVE = load("uniprot_inactive_Q0ZZZ9.json")
# real error text from a page size of 10000 (2026-10-20)
OT_ERROR = {"errors": [{"message": "Argument 'page' has invalid value: At path '': Error during input "
                                   "parameter check. Violations:\n\nThere was a pagination error."}]}


def make_handler(uniprot_entry, ot_body, calls):
    def handler(req: httpx.Request):
        calls.append(req.url.host)
        if req.url.host == "rest.uniprot.org":
            acc = req.url.path.rsplit("/", 1)[-1]
            if uniprot_entry is None:
                return httpx.Response(404, json={"url": f"http://rest.uniprot.org/uniprotkb/{acc}",
                                                 "messages": ["Resource not found"]})
            return httpx.Response(200, json=uniprot_entry)
        if req.url.host == "api.platform.opentargets.org":
            gid = json.loads(req.content)["variables"]["id"]
            assert "." not in gid, gid                     # never send a versioned ID
            return httpx.Response(200, json=ot_body)
        return httpx.Response(404)
    return handler


async def run(uniprot_entry=KRAS, ot_body=OT_KRAS, acc="P01116", **kw):
    calls = []
    transport = httpx.MockTransport(make_handler(uniprot_entry, ot_body, calls))
    async with httpx.AsyncClient(transport=transport) as c:
        r = await target_profile.get_protein_profile(acc, _client=c, **kw)
    return r, calls


async def main():
    tmp = tempfile.mkdtemp()
    common.DOSSIER_DIR = Path(tmp)

    r, calls = await run()
    assert r["ok"] and r["dossier_path"], r
    up, ot = r["uniprot"], r["open_targets"]
    assert up["gene"] == "KRAS" and up["ensembl_gene_ids"] == ["ENSG00000133703"]
    assert up["chembl_target_ids"] == ["CHEMBL2189121"]
    assert up["function"]["truncated"] and "by similarity" in up["function"]["evidence"]
    assert any("Isoform 2B" in s for s in up["subcellular_locations"])
    # 10 diseases, 2 without a name: name stays None and the note is quoted
    assert up["diseases"]["count"] == 10
    unnamed = [d for d in up["diseases"]["items"] if d["name"] is None]
    assert len(unnamed) == 2 and all(d["note"] for d in unnamed), unnamed
    assert "pylocytic astrocytoma" in unnamed[0]["note"]
    aml = up["diseases"]["items"][0]
    assert aml["name"] == "Leukemia, acute myelogenous" and aml["links"]["omim"].endswith("/601626")
    # variants: substitutions get a change, the G10 -> GG insertion keeps its feature id
    assert up["natural_variants"]["count"] == 35
    assert {"G12C", "G12D", "Q61H"} <= set(up["natural_variants"]["changes"])
    assert "VAR_034601" in up["natural_variants"]["changes"]
    assert up["mutagenesis_count"] == 8

    assert ot["found"] and ot["uniprot_crosscheck"] and ot["data_version"] == "26.09"
    a = ot["associations"]
    assert (a["total"], a["examined"]) == (2262, 10)
    assert a["top_overall"][0]["name"] == "Noonan syndrome"              # germline first by score
    assert a["top_oncology"][0]["name"] == "non-small cell lung carcinoma"
    assert a["oncology_in_examined"] == len(a["top_oncology"])
    sm = ot["tractability"]["small molecule"]
    assert "Approved Drug" in sm["yes"] and "Druggable Family" in sm["no"]
    assert a["top_oncology"][0]["link"].endswith("/evidence/ENSG00000133703/MONDO_0005233")
    assert calls == ["rest.uniprot.org", "api.platform.opentargets.org"]

    # the dossier holds the full result; the reply stays compact
    full = json.loads(Path(r["dossier_path"]).read_text(encoding="utf-8"))
    assert Path(r["dossier_path"]).parent.name == "P01116"
    assert len(full["open_targets"]["associations"]["rows"]) == 10
    assert full["uniprot"]["natural_variants"][0]["pubmed"] == ["17468812"]
    reply_kb = len(json.dumps(r)) / 1024
    assert reply_kb < 15, reply_kb

    r2, _ = await run(max_diseases=3)
    assert len(r2["open_targets"]["associations"]["top_overall"]) == 3

    # Open Targets: null target (versioned ID / missing) -> found False, not an error
    r3, _ = await run(ot_body=OT_NULL)
    assert r3["ok"] and r3["open_targets"]["found"] is False and "sin registro" in r3["open_targets"]["note"]

    # Open Targets GraphQL error -> isolated, UniProt part still returned
    r4, _ = await run(ot_body=OT_ERROR)
    assert r4["ok"] and r4["open_targets"]["ok"] is False and "pagination" in r4["open_targets"]["error"]
    assert r4["uniprot"]["gene"] == "KRAS"

    # non-human entry: Open Targets is never called
    mouse = copy.deepcopy(KRAS)
    mouse["organism"]["taxonId"] = 10090
    r5, calls5 = await run(uniprot_entry=mouse)
    assert r5["open_targets"]["found"] is False and calls5 == ["rest.uniprot.org"]

    # UniProt: inactive entry answers 200 -> error; 404 -> error; bad accession -> error
    r6, _ = await run(uniprot_entry=INACTIVE, acc="Q0ZZZ9")
    assert r6["ok"] is False and "inactive" in r6["error"], r6
    r7, _ = await run(uniprot_entry=None, acc="P0DZZ9")
    assert r7["ok"] is False and "not found" in r7["error"]
    r8 = await target_profile.get_protein_profile("not an id")
    assert r8["ok"] is False

    print(json.dumps({"reply_kb": round(reply_kb, 1), "associations": {k: a[k] for k in ("total", "examined", "oncology_in_examined")}}, indent=1))
    print("ALL OK")


asyncio.run(main())
