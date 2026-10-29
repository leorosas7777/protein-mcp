"""PDBe PISA: crystal interfaces and proposed multimers of a PDB entry (Phase 2, step 5).

Validated 2026-10-28 (ESTADO_PROYECTO.md 8.2): only the legacy XML services answer
(cgi-bin/interfaces.pisa and multimers.pisa); the JSON endpoints /pdbe/api/pisa/interfaces and
/interfacelist give 404. Per interface: the two molecules (auth chain, class, symmetry operator
as PISA writes it, e.g. "Z,X,Y"), interface area, css (PISA's complexation significance score)
and per residue the buried surface area. Multimers come in sets: each set divides ALL the chains
of the crystal into assemblies (3WD5: set 1 = A[3] + LH, set 2 = A[3] + L + H, set 3 = LH + A),
so the sets are PISA's alternatives and the assemblies of one set coexist. PISA's verdict is
reported as PISA's, never as a fact. The two services are fetched separately: one can fail
alone (cryo-EM 5A1A: interfaces Ok, multimers "No symmetry operations"). Imports only common.
"""
import xml.etree.ElementTree as ET

import httpx

from protein_mcp.common import log

PISA_INTERFACES = "https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?{pdb}"
PISA_MULTIMERS = "https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?{pdb}"
POLYMER_CLASSES = {"Protein", "DNA", "RNA"}


def _float(text):
    try:
        return float(text) + 0.0          # PISA writes "-0" for some css (4Y6O): show 0.0
    except (TypeError, ValueError):
        return None


def _int(text):
    try:
        return int(text)
    except (TypeError, ValueError):
        return None


def _text(node, tag):
    value = node.findtext(tag) if node is not None else None
    return " ".join(value.split()) if value else None


def _entry(root):
    entry = root.find("pdb_entry")
    status = _text(entry, "status") if entry is not None else _text(root, "status")
    return (entry, None) if status == "Ok" else (None, f"PISA status: {status or 'missing'}")


def parse_interfaces(xml: bytes) -> dict:
    """Polymer-polymer interfaces only (ligand-protein ones are PISA's too, but pockets cover ligands).
    An interface without an id or two molecules is skipped and counted."""
    entry, error = _entry(ET.fromstring(xml))
    if error:
        return {"ok": False, "error": error}
    out, unreadable = [], 0
    for itf in entry.iter("interface"):
        mols = itf.findall("molecule")
        if len(mols) != 2 or any(_text(m, "class") not in POLYMER_CLASSES for m in mols):
            continue
        iid = _int(itf.findtext("id"))
        if iid is None or any(not _text(m, "chain_id") or not _text(m, "symop") for m in mols):
            unreadable += 1
            continue
        out.append({
            "id": iid,
            "area": _float(itf.findtext("int_area")),
            "css": _float(itf.findtext("css")),
            "n_occ": _int(itf.findtext("n_occ")),
            "molecules": [{
                "chain": _text(m, "chain_id"), "class": _text(m, "class"), "symop": _text(m, "symop"),
                # residues with buried surface (PISA's own interface residue list), auth numbers + ins code
                "residues": [f"{_text(r, 'seq_num')}{(_text(r, 'ins_code') or '').strip()}" for r in m.iter("residue")
                             if (_float(r.findtext("bsa")) or 0) > 0],
            } for m in mols],
        })
    return {"ok": True, "interfaces": out, "unreadable_interfaces": unreadable}


def parse_multimers(xml: bytes) -> dict:
    """Every set with ALL its assemblies (an assembly can be in several sets: 3WD5 A[3] in 1 and 2)."""
    entry, error = _entry(ET.fromstring(xml))
    if error:
        return {"ok": False, "error": error}
    sets = []
    for s in entry.iter("asm_set"):
        sets.append({"set": _text(s, "ser_no"), "assemblies": [{
            "id": _text(a, "id"), "size": _text(a, "size"),
            "composition": _text(a, "composition") or _text(a, "formula"),
            "verdict": _text(a, "score"),
        } for a in s.findall("assembly")]})
    return {"ok": True, "sets": sets, "total_assemblies": _int(entry.findtext("total_asm"))}


async def fetch_pisa(client, pdb_id: str) -> dict:
    """{sources, interfaces | interfaces_error, sets + total_assemblies | multimers_error}. Each
    service is independent; any failure (HTTP, XML, unexpected content) is reported, never raised."""
    pdb = pdb_id.lower()
    urls = {"interfaces": PISA_INTERFACES.format(pdb=pdb), "multimers": PISA_MULTIMERS.format(pdb=pdb)}
    out = {"sources": list(urls.values())}
    for key, url in urls.items():
        try:
            resp = await client.get(url, timeout=120.0)
            resp.raise_for_status()
            parsed = (parse_interfaces if key == "interfaces" else parse_multimers)(resp.content)
        except Exception as exc:          # noqa: BLE001 - a PISA failure must never stop the tool
            log.error("PISA %s failed for %s: %s", key, pdb_id, exc)
            parsed = {"ok": False, "error": f"request or reading failed: {exc!r}"}
        if not parsed["ok"]:
            out[f"{key}_error"] = parsed["error"]
        elif key == "interfaces":
            out["interfaces"] = parsed["interfaces"]
            if parsed["unreadable_interfaces"]:
                out["interfaces_note"] = f"{parsed['unreadable_interfaces']} interface(s) without id/chain/operator skipped"
        else:
            out["sets"], out["total_assemblies"] = parsed["sets"], parsed["total_assemblies"]
    return out
