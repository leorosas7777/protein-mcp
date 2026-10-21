"""Offline test for the dossier helpers in common.py (file names, saving, failures).
Run: uv run python tests/test_dossier.py
"""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # project root
import common  # noqa: E402


def main():
    common.DOSSIER_DIR = Path(tempfile.mkdtemp())

    # one file per parameter combination; the default call is just the base name
    assert common.dossier_name("ligands", "") == "ligands"
    assert common.dossier_name("ligands", "G12C") == "ligands_g12c"
    assert common.dossier_name("literature", "crystal structure", "recent", "") == "literature_crystal-structure_recent"
    assert common.dossier_name("literature", "../../etc") == "literature_etc"   # no path traversal
    assert common.dossier_name("structures", "  ") == "structures"

    result = {"ok": True, "uniprot_id": "P01116", "papers": [{"title": "Targeting KRAS for cancer therapy."}]}
    reply = common.attach_dossier(result, "P01116", "literature", "search_literature", {"topic": ""})
    path = Path(reply["dossier_path"])
    assert path == common.DOSSIER_DIR / "P01116" / "literature.json"
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["tool"] == "search_literature" and saved["result"] == result and saved["retrieved_at"]
    assert "dossier_path" not in result                       # the original is not modified

    # compact goes to the chat, the full result to disk
    reply = common.attach_dossier(result, "P01116", "x", "t", {}, compact=lambda r: {"ok": True, "n": len(r["papers"])})
    assert reply == {"ok": True, "n": 1, "dossier_path": str(common.DOSSIER_DIR / "P01116" / "x.json")}

    # failures are returned unchanged and never saved; neither is a call without accession
    err = {"ok": False, "error": "UniProt accession P0DZZ9 not found."}
    assert common.attach_dossier(err, "P0DZZ9", "literature", "t", {}) is err
    assert not (common.DOSSIER_DIR / "P0DZZ9").exists()
    assert common.attach_dossier(result, None, "literature", "t", {}) is result

    # a write failure keeps the reply and says so
    common.DOSSIER_DIR = Path(tempfile.mkdtemp()) / "file"
    common.DOSSIER_DIR.write_text("not a folder")
    reply = common.attach_dossier(result, "P01116", "literature", "t", {})
    assert reply["dossier_path"] is None and "could not be written" in reply["dossier_note"]
    print("ALL OK")


main()
