"""Shared helpers: logging, accession validation, text cleaning, safe wrapper, dossier.

IMPORTANT: never use print() anywhere in this project (stdout is the MCP
protocol channel). Use `log` (goes to stderr).
"""
import html
import json
import logging
import re
from datetime import datetime, timezone
from pathlib import Path

log = logging.getLogger("protein-mcp")

# Full tool results are written here, one folder per accession (git-ignored).
# The chat only gets a compact summary plus the path; build_report reads from disk.
DOSSIER_DIR = Path(__file__).resolve().parent / "dossiers"

ACCESSION_RE = re.compile(r"[A-Z0-9]{6,10}")
_INLINE_TAG_RE = re.compile(r"</?(?:sup|sub|i|b|em|strong|u|span)\b[^>]*>", re.I)
_TAG_RE = re.compile(r"</?[a-zA-Z][^>]*>")


def parse_accession(uniprot_id: str):
    """Return (accession, None) if valid, else (None, error_dict)."""
    accession = uniprot_id.strip().upper()
    if not ACCESSION_RE.fullmatch(accession):
        return None, {
            "ok": False,
            "error": (
                f"'{uniprot_id}' does not look like a UniProt accession "
                "(e.g. P01116). Use resolve_protein first."
            ),
        }
    return accession, None


def clean_text(value) -> str:
    """Clean HTML that Europe PMC / RCSB leave in titles and abstracts.

    Entities are decoded FIRST: Europe PMC often sends escaped tags such as
    "&lt;sup&gt;", which would turn into real tags if decoded after stripping.
    Inline tags (sup/sub/i/b...) are removed without a space (KRAS<sup>G13C</sup>
    -> KRASG13C); block tags (h4, p...) become a space.
    """
    if not value:
        return ""
    text = html.unescape(str(value))
    text = _INLINE_TAG_RE.sub("", text)
    text = _TAG_RE.sub(" ", text)
    return re.sub(r"\s+", " ", text).strip()


def short(value, limit: int = 80):
    text = clean_text(value)
    if len(text) > limit:
        return text[: limit - 3] + "..."
    return text or None


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


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


async def safe(coro, label: str) -> dict:
    """The API shapes are partly unverified: never let one source crash the tool."""
    try:
        return await coro
    except Exception as exc:  # noqa: BLE001
        log.exception("%s failed", label)
        return {"ok": False, "error": f"{label} failed unexpectedly: {exc!r}"}
