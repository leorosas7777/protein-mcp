"""Dossier on disk: every tool saves its FULL result to dossiers/<accession>/<name>.json
(git-ignored); the chat only gets a compact summary plus the path, and build_report
reads from disk. Imports only common.
"""
import json
import os
import re
from pathlib import Path

from protein_mcp.common import log, utc_now

# Default: <repo>/dossiers (the folder that holds pyproject.toml). PROTEIN_MCP_DOSSIER_DIR
# overrides it, e.g. when the package is installed outside the repository.
DOSSIER_DIR = Path(os.environ.get("PROTEIN_MCP_DOSSIER_DIR") or Path(__file__).resolve().parents[1] / "dossiers")


def save_dossier(accession: str, name: str, data: dict):
    """Write `data` to dossiers/<accession>/<name>.json. Return the path, or None on failure.

    `accession` must already be validated by parse_accession (letters and digits only).
    """
    folder = DOSSIER_DIR / accession
    path = folder / f"{name}.json"
    try:
        folder.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
    except OSError:
        log.exception("Could not write dossier file %s", path)
        return None
    return str(path)


def read_dossier(accession: str, name: str):
    """The saved {tool, params, retrieved_at, result} of dossiers/<accession>/<name>.json, or None
    if it is missing or unreadable. `name` must come from dossier_name (no path separators)."""
    path = DOSSIER_DIR / accession / f"{name}.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except (OSError, ValueError):
        log.exception("Could not read dossier file %s", path)
        return None
    return data if isinstance(data, dict) and "result" in data else None


def dossier_name(base: str, *parts) -> str:
    """File name for one call, e.g. ("ligands", "G12C") -> "ligands_g12c". Empty parts are
    skipped, so the default call is just `base`; other parameters get their own file."""
    slugs = [re.sub(r"[^a-z0-9]+", "-", str(p).lower()).strip("-")[:40] for p in parts if p]
    return "_".join([base] + [s for s in slugs if s])


def attach_dossier(result: dict, accession, name: str, tool: str, params: dict, compact=None) -> dict:
    """Save the FULL `result` to the dossier and return what goes to the chat:
    `compact(result)` (or the result unchanged) plus `dossier_path`.
    Failed results are not saved and are returned as they are."""
    if not accession or not result.get("ok"):
        return result
    path = save_dossier(accession, name, {
        "tool": tool, "params": params, "retrieved_at": utc_now(), "result": result,
    })
    reply = compact(result) if compact else dict(result)
    reply["dossier_path"] = path
    if path is None:
        reply["dossier_note"] = "The dossier could not be written (see server log); only this reply exists."
    return reply
