"""Offline test for the UniProt accession check (uniprot.check_accession) and how get_ligands,
get_structures, get_protein_profile and search_literature use it. The answers in
fixtures/uniprot_accession_checks.json are real UniProt responses (2026-11-04): existing,
nonexistent, deleted, merged (303 with an "Inactive" body) and malformed.
Run: uv run python tests/test_accession.py
"""
import asyncio
import json
import tempfile
from pathlib import Path

import httpx

from protein_mcp import dossier, server
from protein_mcp.sources import chembl, uniprot

ROOT = Path(__file__).resolve().parent.parent
REAL = json.loads((ROOT / "tests" / "fixtures" / "uniprot_accession_checks.json").read_text(encoding="utf-8"))


def real_uniprot(req):
    acc = req.url.path.rsplit("/", 1)[-1]
    answer = REAL[acc]
    headers = {"location": answer["location"]} if answer["location"] else {}
    return httpx.Response(answer["status"], json=answer["json"], headers=headers)


async def check(acc, handler=real_uniprot, **client_args):
    async with httpx.AsyncClient(transport=httpx.MockTransport(handler), **client_args) as c:
        return await uniprot.check_accession(acc, _client=c)


async def check_answers():
    assert await check("P01116") == {"ok": True}
    r = await check("P0DZZ9")
    assert r["ok"] is False and "not found" in r["error"]
    r = await check("Q0ZZZ9")
    assert r["ok"] is False and "inactive" in r["error"] and "Not part of a reference proteome" in r["error"]
    r = await check("P01118")   # merged: named, never replaced by the new accession
    assert r["ok"] is False and "MERGED" in r["error"] and "P01116" in r["error"]
    r = await check("ABCDEF")
    assert r["ok"] is False and "not a valid" in r["error"]

    # UniProt not answering proves nothing about the accession: warning, not an error
    def down(req):
        raise httpx.ConnectError("connection refused")

    for handler, text in ((down, "ConnectError"),
                          (lambda req: httpx.Response(503, text="<html>down</html>"), "HTTP 503"),
                          (lambda req: httpx.Response(200, text="<html>not json</html>"), "JSONDecodeError"),
                          (lambda req: httpx.Response(303, json={"entryType": "?"}), "HTTP 303")):
        r = await check("P01116", handler)
        assert r["ok"] is True and "could not confirm" in r["warning"] and text in r["warning"], (text, r)


async def fetch_entry_refuses_merged():
    """get_protein_profile and search_literature follow redirects: without this check they saved
    P01116's data under P01118. Real flow: 303 -> the P01116 entry."""
    def handler(req):
        if "P01118" in req.url.path:
            return real_uniprot(req)
        return httpx.Response(200, json=REAL["P01116"]["json"])

    async with httpx.AsyncClient(transport=httpx.MockTransport(handler), follow_redirects=True) as c:
        r = await uniprot.fetch_entry_json(c, "P01118")
        assert r["ok"] is False and "secondary accession of P01116" in r["error"], r
        assert (await uniprot.fetch_entry_json(c, "P01116"))["ok"] is True


async def tools_use_the_check():
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())
    calls = []
    verdict = {}

    async def fake_check(acc):
        return verdict

    async def fake_ligands(acc, mutation, max_ligands):
        calls.append("chembl")
        return {"ok": True, "uniprot_id": acc, "target": None, "message": "No single-protein ChEMBL target"}

    async def fake_pdb(client, acc, limit, mutation):
        calls.append("rcsb")
        return {"ok": True, "total_entries": 0, "models": []}

    async def fake_af(client, acc):
        calls.append("alphafold")
        return {"ok": False, "error": "no model"}

    uniprot.check_accession, chembl.get_ligands = fake_check, fake_ligands
    server.fetch_pdb, server.fetch_alphafold = fake_pdb, fake_af

    # rejected: the error reaches the chat, no source is called, nothing is saved
    verdict.update(ok=False, error="UniProt accession P0DZZ9 not found.")
    for tool in (server.get_ligands, server.get_structures):
        r = await tool("P0DZZ9")
        assert r["ok"] is False and "not found" in r["error"] and r["uniprot_id"] == "P0DZZ9", r
    assert not calls and not (dossier.DOSSIER_DIR / "P0DZZ9").exists()

    # a malformed id never reaches UniProt or the dossier
    assert (await server.get_ligands("not an id"))["ok"] is False

    # confirmed: an empty result is a real "no data", no warning
    verdict.clear()
    verdict.update(ok=True)
    r = await server.get_ligands("P01116")
    assert r["ok"] and r["target"] is None and "accession_warning" not in r and Path(r["dossier_path"]).exists()

    # UniProt down: the sources are still asked and the warning goes to the chat and to the
    # notes that build_report prints, in both tools
    verdict.update(warning="UniProt could not confirm that accession P01116 exists (HTTP 503)")
    for tool in (server.get_ligands, server.get_structures):
        r = await tool("P01116")
        assert "could not confirm" in r["accession_warning"] and r["accession_warning"] in r["notes"], r
        saved = json.loads(Path(r["dossier_path"]).read_text(encoding="utf-8"))
        assert saved["result"]["accession_warning"] and saved["result"]["accession_warning"] in saved["result"]["notes"]
    assert {"chembl", "rcsb", "alphafold"} <= set(calls)


def main():
    asyncio.run(check_answers())
    asyncio.run(fetch_entry_refuses_merged())
    asyncio.run(tools_use_the_check())
    print("ALL OK")


main()
