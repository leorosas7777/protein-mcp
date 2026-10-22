"""Offline test for the dossier helpers in protein_mcp/dossier.py (file names, saving, failures).
Run: uv run python tests/test_dossier.py
"""
import json
import tempfile
from pathlib import Path

from protein_mcp import dossier


def main():
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())

    # one file per parameter combination; the default call is just the base name
    assert dossier.dossier_name("ligands", "") == "ligands"
    assert dossier.dossier_name("ligands", "G12C") == "ligands_g12c"
    assert dossier.dossier_name("literature", "crystal structure", "recent", "") == "literature_crystal-structure_recent"
    assert dossier.dossier_name("literature", "../../etc") == "literature_etc"   # no path traversal
    assert dossier.dossier_name("structures", "  ") == "structures"

    result = {"ok": True, "uniprot_id": "P01116", "papers": [{"title": "Targeting KRAS for cancer therapy."}]}
    reply = dossier.attach_dossier(result, "P01116", "literature", "search_literature", {"topic": ""})
    path = Path(reply["dossier_path"])
    assert path == dossier.DOSSIER_DIR / "P01116" / "literature.json"
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["tool"] == "search_literature" and saved["result"] == result and saved["retrieved_at"]
    assert "dossier_path" not in result                       # the original is not modified

    # compact goes to the chat, the full result to disk
    reply = dossier.attach_dossier(result, "P01116", "x", "t", {}, compact=lambda r: {"ok": True, "n": len(r["papers"])})
    assert reply == {"ok": True, "n": 1, "dossier_path": str(dossier.DOSSIER_DIR / "P01116" / "x.json")}

    # failures are returned unchanged and never saved; neither is a call without accession
    err = {"ok": False, "error": "UniProt accession P0DZZ9 not found."}
    assert dossier.attach_dossier(err, "P0DZZ9", "literature", "t", {}) is err
    assert not (dossier.DOSSIER_DIR / "P0DZZ9").exists()
    assert dossier.attach_dossier(result, None, "literature", "t", {}) is result

    # a write failure keeps the reply and says so
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp()) / "file"
    dossier.DOSSIER_DIR.write_text("not a folder")
    reply = dossier.attach_dossier(result, "P01116", "literature", "t", {})
    assert reply["dossier_path"] is None and "could not be written" in reply["dossier_note"]
    print("ALL OK")


main()
