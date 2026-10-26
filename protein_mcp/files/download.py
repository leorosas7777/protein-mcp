"""Download one file to disk with its provenance (URL, size, SHA-256, time) after checking
that the content is what was asked for. Imports only common.

Checks are needed because some servers answer HTTP 200 without the data: ModelServer
returns a property-only SDF with no molecule for an unknown label_asym_id (3OG7 "ZZ").
"""
import asyncio
import hashlib
from pathlib import Path

import httpx

from protein_mcp.common import log, utc_now

MAX_FILE_BYTES = 300 * 1024 * 1024  # large cryo-EM mmCIF files can exceed 100 MB
RETRIES = 2          # extra attempts after a transient failure
RETRY_WAIT = 1.0     # seconds before the first retry, doubled for the next
_TRANSIENT_HTTP = {429, 500, 502, 503, 504}
_HEAD_BYTES = 4096


def check_mmcif(entry_id: str = ""):
    """mmCIF starts with 'data_<ID>' (RCSB: 'data_3OG7'; AlphaFold: 'data_AF-P15056-F1')."""
    def check(head: bytes):
        first = head.split(b"\n", 1)[0].strip().decode("ascii", "replace")
        if not first.lower().startswith("data_"):
            return f"not an mmCIF file (first line: {first[:60]!r})"
        if entry_id and first[5:].upper() != entry_id.upper():
            return f"mmCIF is for '{first[5:]}', expected '{entry_id}'"
        return None
    return check


def check_fasta(accession: str):
    """UniProt FASTA header: '>sp|P15056|BRAF_HUMAN ...'."""
    def check(head: bytes):
        first = head.split(b"\n", 1)[0].decode("ascii", "replace")
        if not first.startswith(">") or f"|{accession}|" not in first:
            return f"not the UniProt FASTA of {accession} (header: {first[:60]!r})"
        return None
    return check


def check_sdf(head: bytes):
    """Molfile counts line (4th line) must declare at least one atom. ModelServer omits the
    'V2000' tag there, so only the atom count is read."""
    lines = head.decode("ascii", "replace").splitlines()
    try:
        atoms = int(lines[3][:3])
    except (IndexError, ValueError):
        return "SDF without a molecule (the server found no atoms for this selection)"
    return None if atoms > 0 else "SDF with 0 atoms"


async def download(client: httpx.AsyncClient, url: str, path: Path, check) -> dict:
    """Stream `url` to `path`. Return {ok, url, path, bytes, sha256, retrieved_at} or
    {ok: False, url, error}. A failed or rejected download never replaces an existing file.
    Transient failures (connection cut, timeout, HTTP 429/5xx) are retried RETRIES times;
    `attempts` is added when more than one was needed. A rejected content is not retried."""
    for attempt in range(1, RETRIES + 2):
        result = await _download_once(client, url, path, check)
        transient = result.pop("transient", False)
        if result["ok"] or not transient or attempt > RETRIES:
            if attempt > 1:
                result["attempts"] = attempt
                if not result["ok"]:
                    result["error"] += f" (after {attempt} attempts)"
            return result
        log.warning("Transient download failure %s (attempt %d): %s", url, attempt, result["error"])
        await asyncio.sleep(RETRY_WAIT * 2 ** (attempt - 1))


async def _download_once(client: httpx.AsyncClient, url: str, path: Path, check) -> dict:
    part = path.with_name(path.name + ".part")
    digest, size = hashlib.sha256(), 0
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        async with client.stream("GET", url) as resp:
            if resp.status_code != 200:
                return {"ok": False, "url": url, "error": f"HTTP {resp.status_code}",
                        "transient": resp.status_code in _TRANSIENT_HTTP}
            with part.open("wb") as fh:
                async for chunk in resp.aiter_bytes():
                    size += len(chunk)
                    if size > MAX_FILE_BYTES:
                        break
                    digest.update(chunk)
                    fh.write(chunk)
        if size > MAX_FILE_BYTES:
            part.unlink(missing_ok=True)
            return {"ok": False, "url": url, "error": f"larger than {MAX_FILE_BYTES // 2**20} MB, not downloaded"}
        with part.open("rb") as fh:
            problem = check(fh.read(_HEAD_BYTES)) if check else None
        if problem:
            part.unlink(missing_ok=True)
            return {"ok": False, "url": url, "error": problem}
        part.replace(path)
    except (httpx.HTTPError, OSError) as exc:
        log.error("Download failed %s: %r", url, exc)
        part.unlink(missing_ok=True)
        # 6V34 QOP C (2026-10-26): RemoteProtocolError "Server disconnected" on one SDF of 18
        return {"ok": False, "url": url, "error": f"download failed: {exc!r}",
                "transient": isinstance(exc, httpx.TransportError)}
    return {"ok": True, "url": url, "path": str(path), "bytes": size,
            "sha256": digest.hexdigest(), "retrieved_at": utc_now()}
