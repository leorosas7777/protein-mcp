"""pymol_view: change the view of an OPEN PyMOL session from the chat (Phase 2, step 4b).

The .pml written by open_in_pymol starts a receiver inside PyMOL (files.viewer._RECEIVER) that
writes prep/pymol/control.json (port on 127.0.0.1, random token, session id) and accepts only
viewer.LIVE_OPS with data arguments. Here each action is checked against the session's
scene_<session>.json and the manifest, turned into those operations and sent; the chat never
sends PyMOL text. UniProt numbers are translated to the author numbering PyMOL uses with gemmi
on the raw mmCIF, after checking its SHA-256 against the manifest. Nothing is computed about
the structure and nothing is saved: closing PyMOL or reopening the session loses the changes.
Imports common, files.{clean, viewer, prepare (prep_dir)}.
"""
import hashlib
import json
import re
import socket

from protein_mcp.common import log, parse_accession
from protein_mcp.files import clean, viewer
from protein_mcp.files.prepare import prep_dir

ACTIONS = ("status", "show_structure", "focus_pocket", "color_chain", "show_residues", "hide", "show", "reset")
TARGETS = ("ligands", "boxes", "alphafold")
# Plain names the chat may ask for -> PyMOL colour names (a closed list: no colour text from the chat)
COLORS = {"red": "red", "green": "green", "blue": "blue", "yellow": "yellow", "orange": "orange",
          "magenta": "magenta", "cyan": "cyan", "white": "white", "grey": "grey80", "gray": "grey80",
          "pink": "pink", "purple": "purple", "salmon": "salmon", "wheat": "wheat", "slate": "slate"}
# not used by the scene (grey receptor, yellow ligands, light blue pockets, magenta variant,
# AlphaFold pLDDT blue/light blue/yellow/orange): orange was lost on low-pLDDT AlphaFold (2026-10-28)
DEFAULT_RESIDUE_COLOR = "green"
MAX_RESIDUES = 100
TIMEOUT_SECONDS = 10.0
_RANGE_RE = re.compile(r"\s*(\d{1,5})\s*(?:-\s*(\d{1,5}))?\s*")
NOT_SEEN = "PyMOL accepted the commands; the server cannot see the window, so it cannot confirm what is drawn."


class _Reject(Exception):
    """A request that does not match the session or the manifest: nothing is sent."""
    def __init__(self, error, options=None):
        super().__init__(error)
        self.error, self.options = error, options


def _read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def send(control: dict, ops: list) -> dict:
    """Send [[function, args], ...] to the receiver. {ok, results} or {ok: False, error}."""
    try:
        s = socket.create_connection(("127.0.0.1", int(control["port"])), timeout=TIMEOUT_SECONDS)
    except OSError as exc:
        # a closed PyMOL leaves control.json behind; nothing listens on its port any more
        return {"ok": False, "not_open": True,
                "error": f"PyMOL is not open ({exc.__class__.__name__}). Open it again with open_in_pymol."}
    try:
        with s:
            s.sendall((json.dumps({"token": control["token"], "ops": ops}) + "\n").encode())
            line = s.makefile("rb").readline(10_000_000)
    except OSError as exc:
        # reviewer 2026-10-28: once sent, PyMOL may have applied the commands even if no answer came
        return {"ok": False, "error": f"The commands were sent but PyMOL gave no answer within {TIMEOUT_SECONDS:g} s "
                                      f"({exc.__class__.__name__}): they may or may not have been applied. Ask the "
                                      "user to look at the window; use status to check."}
    try:
        reply = json.loads(line)
    except ValueError:
        return {"ok": False, "error": "The process on the control port did not answer as the PyMOL receiver."}
    if reply.get("session") != control.get("session"):
        return {"ok": False, "error": "The process on the control port belongs to another session; reopen with open_in_pymol."}
    return reply


def _structure(scene: dict, name: str):
    """('AlphaFold', None) or ('experimental', scene entry)."""
    entries = {e["pdb_id"]: e for e in scene.get("structures") or []}
    name = (name or "").strip()
    af = scene.get("alphafold")
    if af and name.lower() in ("alphafold", str(af.get("entry_id") or "").lower()):
        return "alphafold", None
    if name.upper() in entries:
        return "experimental", entries[name.upper()]
    raise _Reject(f"structure '{name}' is not in this PyMOL session" if name else "structure is required for this action",
                  {"structures": list(entries) + (["AlphaFold"] if af else [])})


def _show_only(scene: dict, group=None, af=False) -> list:
    """Enable one group (or the AlphaFold object) and disable the others."""
    ops = [["disable", [e["group"]]] for e in scene.get("structures") or [] if e["group"] != group]
    if scene.get("alphafold") and not af:
        ops.append(["disable", [scene["alphafold"]["object"]]])
    ops.append(["enable", [group or scene["alphafold"]["object"]]])
    return ops


def _pocket(entry: dict, text: str) -> dict:
    """Pocket row of the scene from 'comp label', 'comp' (if unique) or '<PDB> comp label'."""
    rows = [p for p in entry.get("pockets") or [] if p.get("selection")]
    names = [f"{p['comp_id']} {p['label_asym_id']}" for p in rows]
    parts = (text or "").split()
    if len(parts) == 3 and parts[0].upper() == entry["pdb_id"]:
        parts = parts[1:]
    if not parts or len(parts) > 2:
        raise _Reject(f"pocket '{text}' is not '<comp_id> <label_asym_id>' (e.g. '032 C')", {"pockets": names})
    hits = [p for p in rows if p["comp_id"].upper() == parts[0].upper()
            and (len(parts) == 1 or p["label_asym_id"] == parts[1])]
    if len(hits) != 1:
        raise _Reject((f"no pocket '{text}' shown for {entry['pdb_id']} in this session" if not hits else
                       f"'{text}' matches {len(hits)} pockets in {entry['pdb_id']}: give the instance too"),
                      {"pockets": names})
    return hits[0]


def _manifest_record(manifest: dict, pdb_id: str) -> dict:
    rec = next((s for s in manifest.get("structures") or [] if s.get("pdb_id") == pdb_id), None)
    if not rec or not (rec.get("prepared") or {}).get("ok"):
        raise _Reject(f"{pdb_id} has no prepared record in the manifest")
    return rec


def _target_chains(rec: dict) -> list:
    return [c["auth_asym_id"] for c in rec["prepared"]["metrics"]["chains"]]


def _chain(rec: dict, chain: str) -> str:
    chains = _target_chains(rec)
    chain = (chain or "").strip()
    if not chain:
        if len(chains) == 1:
            return chains[0]
        raise _Reject(f"{rec['pdb_id']} has {len(chains)} target chains ({', '.join(chains)}): name the chain",
                      {"chains": chains})
    if chain not in chains:
        raise _Reject(f"chain {chain} is not a target chain of {rec['pdb_id']} (auth ids, case-sensitive)",
                      {"chains": chains})
    return chain


def _color(color: str, default=None) -> str:
    key = (color or "").strip().lower()
    if not key and default:
        return default
    if key not in COLORS:
        raise _Reject(f"colour '{color}' is not in the list" if key else "color is required", {"colors": sorted(COLORS)})
    return COLORS[key]


def _numbers(text: str) -> list:
    """'594-596, 600' -> [594, 595, 596, 600]."""
    out = []
    for part in str(text or "").replace(";", ",").split(","):
        if not part.strip():
            continue
        m = _RANGE_RE.fullmatch(part)
        if not m:
            raise _Reject(f"residues '{text}' must be UniProt numbers or ranges, e.g. '594-596, 600'")
        a, b = int(m.group(1)), int(m.group(2) or m.group(1))
        if b < a or b - a >= MAX_RESIDUES:
            raise _Reject(f"range '{part.strip()}' is not valid (at most {MAX_RESIDUES} residues)")
        out += [n for n in range(a, b + 1) if n not in out]
    if not out:
        raise _Reject("residues is required, e.g. '594-596, 600' (UniProt numbering)")
    if len(out) > MAX_RESIDUES:
        raise _Reject(f"at most {MAX_RESIDUES} residues per call")
    return out


def _raw_file(prep, record: dict, what: str):
    """The raw mmCIF of the manifest, only if its SHA-256 still matches the manifest."""
    f = record.get("file") or {}
    path = (prep / str(f.get("path") or "")).resolve()
    if not f.get("ok") or prep.resolve() not in path.parents or not path.is_file():
        raise _Reject(f"the raw file of {what} is not in prep/")
    if hashlib.sha256(path.read_bytes()).hexdigest() != f.get("sha256"):
        raise _Reject(f"the raw file of {what} changed since the manifest was written (SHA-256); run prepare_structures again")
    return path


def plan(action, scene, manifest, prep, structure="", pocket="", chain="", residues="", color="", target=""):
    """(ops, information for the reply). Raises _Reject."""
    info = {}
    if action == "status":
        return [["get_names", ["objects", 1]]], info
    if action == "reset":
        # objects start enabled when loaded and the scene only disables some: enable every group,
        # ligand and box first (hide/show may have disabled them), then replay the scene
        objects = [o for e in scene.get("structures") or [] for o in
                   [e["group"]] + [l["object"] for l in e["ligands"]] + [p["box_object"] for p in e["pockets"] if p.get("box_object")]]
        # the scene's "hide everything" also hides the boxes' cgo representation, which only
        # load_cgo (not live) set: show it again (reviewer 2026-10-28, checked with get_vis in PyMOL 3.2)
        boxes = [["show", ["cgo", p["box_object"]]] for e in scene.get("structures") or []
                 for p in e["pockets"] if p.get("box_object")]
        return ([["enable", [o]] for o in objects] + list(scene.get("reset_ops") or []) + boxes,
                {"done": "the scene of open_in_pymol was drawn again"})
    if action in ("hide", "show"):
        key = (target or "").strip().lower()
        if key not in TARGETS:
            raise _Reject(f"target must be one of {', '.join(TARGETS)}", {"targets": list(TARGETS)})
        if key == "alphafold":
            if not scene.get("alphafold"):
                raise _Reject("the AlphaFold model is not loaded in this session")
            objects = [scene["alphafold"]["object"]]
        else:
            entries = [_structure(scene, structure)[1]] if structure else scene.get("structures") or []
            if any(e is None for e in entries):
                raise _Reject("AlphaFold has no ligands or boxes")
            objects = ([l["object"] for e in entries for l in e["ligands"]] if key == "ligands"
                       else [p["box_object"] for e in entries for p in e["pockets"] if p.get("box_object")])
            if not objects:
                raise _Reject(f"no {key} loaded" + (f" for {structure}" if structure else ""))
        info["objects"] = objects
        ops = [["disable" if action == "hide" else "enable", [o]] for o in objects]
        if action == "show" and key == "boxes":          # in case a reset hid their cgo representation
            ops += [["show", ["cgo", o]] for o in objects]
        return ops, info

    kind, entry = _structure(scene, structure)
    stale = manifest.get("created_at") != scene.get("manifest_created_at")
    if stale and action in ("color_chain", "show_residues"):
        raise _Reject("prepare_structures rewrote the manifest after this session was opened: reopen it with "
                      "open_in_pymol so that chains and residues match")
    if action == "show_structure":
        if kind == "alphafold":
            obj = scene["alphafold"]["object"]
            return _show_only(scene, af=True) + [["orient", [obj]]], {"shown": "AlphaFold " + obj}
        return _show_only(scene, entry["group"]) + [["orient", [entry["receptor"]]]], {"shown": entry["pdb_id"]}
    if action == "focus_pocket":
        if kind == "alphafold":
            raise _Reject("the AlphaFold model has no ligand, so no pocket")
        row = _pocket(entry, pocket)
        ops = _show_only(scene, entry["group"])
        if row.get("box_object"):
            ops.append(["enable", [row["box_object"]]])
        info["pocket"] = {k: row.get(k) for k in ("comp_id", "label_asym_id", "auth_asym_id", "n_residues",
                                                  "n_residues_selected", "selection", "box_object")}
        return ops + [["zoom", [row["selection"], 5]]], info
    if action == "color_chain":
        if kind == "alphafold":
            raise _Reject("color_chain works on the experimental structures; AlphaFold is coloured by pLDDT")
        c = _chain(_manifest_record(manifest, entry["pdb_id"]), chain)
        col = _color(color)
        info["colored"] = f"{entry['pdb_id']} chain {c}: {col} (carbon atoms and cartoon; other elements keep their colours)"
        return [["color", [col, f"{entry['receptor']} and chain {c} and elem C"]]], info
    if action == "show_residues":
        wanted = _numbers(residues)
        col = _color(color, DEFAULT_RESIDUE_COLOR)
        if kind == "alphafold":
            af = manifest.get("alphafold") or {}
            if af.get("is_canonical_sequence") is not True:
                raise _Reject("the AlphaFold model is not the canonical sequence: its numbers are not UniProt numbers")
            # long proteins come as fragments (F2, F3...) whose numbering may not start at 1 (reviewer
            # 2026-10-28, not reproduced): only a model covering the whole canonical sequence is used
            length = (manifest.get("sequence") or {}).get("length")      # residues is text, e.g. "1-766"
            if not str(af.get("entry_id") or "").endswith("-F1") or not length or af.get("residues") != f"1-{length}":
                raise _Reject("residues are shown on AlphaFold only for a single model covering the whole canonical "
                              f"sequence (model residues {af.get('residues')}, UniProt length {length})")
            path, record, c, obj, name = _raw_file(prep, af, "the AlphaFold model"), None, "", scene["alphafold"]["object"], "AlphaFold"
        else:
            rec = _manifest_record(manifest, entry["pdb_id"])
            c, obj, name = _chain(rec, chain), entry["receptor"], entry["pdb_id"]
            path, record = _raw_file(prep, rec, entry["pdb_id"]), rec
        mapping = clean.residue_map(path, record, c)
        if not mapping["ok"]:
            raise _Reject(mapping["error"])
        found = [n for n in wanted if n in mapping["residues"]]
        info["residues"] = [f"{n}: {mapping['residues'][n][0]} (auth {mapping['residues'][n][1]})" for n in found]
        info["not_modelled"] = [n for n in wanted if n not in mapping["residues"]]
        if mapping["ambiguous"]:
            info["ambiguous"] = {str(n): v for n, v in mapping["ambiguous"].items() if n in found}
        if not found:
            raise _Reject(f"none of the residues {residues} is modelled in {name}" + (f" chain {c}" if c else ""),
                          {"not_modelled": info["not_modelled"]})
        sel, _, dropped = viewer._selection(obj, [(c or "A", mapping["residues"][n][1]) for n in found])
        if kind == "alphafold":       # AlphaFold DB models have one chain, A
            sel = sel.replace("(chain A and ", "(")
        if dropped or not sel:
            raise _Reject(f"PyMOL cannot select these residue ids: {dropped}")
        named = viewer._name("res", name)
        info["selection"] = named
        enable = ["enable", [entry["group"] if entry else obj]]
        return [enable, ["select", [named, sel]], ["show", ["sticks", f"{named} and not name N+C+O"]],
                ["color", [col, f"{named} and elem C"]], ["zoom", [named, 5]], ["deselect", []]], info
    raise _Reject(f"action '{action}' is not handled")


def pymol_view(uniprot_id: str, action: str, structure: str = "", pocket: str = "", chain: str = "",
               residues: str = "", color: str = "", target: str = "") -> dict:
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    action = (action or "").strip().lower()
    if action not in ACTIONS:
        return {"ok": False, "error": f"action must be one of: {', '.join(ACTIONS)}"}
    prep = prep_dir(accession)
    folder = prep / "pymol"
    try:
        control = _read_json(folder / viewer.CONTROL_FILE)
        scene = _read_json(folder / f"scene_{control['session']}.json")
    except FileNotFoundError:
        return {"ok": False, "uniprot_id": accession,
                "error": f"No PyMOL session with chat control for {accession}: open it with open_in_pymol first "
                         "(sessions opened before this feature have no control)."}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        log.exception("Could not read the control files in %s", folder)
        return {"ok": False, "uniprot_id": accession, "error": f"Could not read the session files in {folder}: {exc!r}"}
    try:
        manifest = _read_json(prep / "manifest.json")
    except (OSError, ValueError):
        manifest = {}
    try:
        ops, info = plan(action, scene, manifest, prep, structure, pocket, chain, residues, color, target)
    except _Reject as r:
        out = {"ok": False, "uniprot_id": accession, "action": action, "error": f"Nothing sent: {r.error}."}
        if r.options:
            out["options"] = r.options
        return out

    reply = send(control, ops)
    if not reply.get("ok"):
        out = {"ok": False, "uniprot_id": accession, "action": action, "error": reply.get("error")}
        if reply.get("not_open"):
            out["pymol_open"] = False
            # control.json is rewritten by the receiver once the new session has loaded (7 structures +
            # AlphaFold took more than 5 s on 2026-10-28): an older control.json may mean "still loading"
            try:
                if (folder / "session.pml").stat().st_mtime > (folder / viewer.CONTROL_FILE).stat().st_mtime:
                    out["error"] += (" A newer session was written after this control file: if PyMOL was just "
                                     "opened, it may still be loading; try again in a few seconds.")
            except OSError:
                pass
        return out
    results = reply.get("results") or []
    failed = [f"{op[0]}{tuple(op[1])[:2]}: {r.get('error')}" for op, r in zip(ops, results) if not r.get("ok")]
    out = {"ok": not failed, "uniprot_id": accession, "action": action, "session_id": control.get("session"),
           **info, "sent": len(ops), "applied": sum(1 for r in results if r.get("ok")), "failed": failed,
           "note": NOT_SEEN + " Changes are not saved: reopening the session restores the open_in_pymol scene."}
    if action == "status":
        # enabled flags are per object: a ligand stays "enabled" inside a hidden group, so the
        # state is given per structure (2026-10-28, real PyMOL), never as a raw object list
        enabled = set((results[0].get("result") or []) if results and results[0].get("ok") else [])
        out["structures"] = [{
            "pdb_id": e["pdb_id"], "shown": e["group"] in enabled,
            "pockets": [f"{p['comp_id']} {p['label_asym_id']}" for p in e["pockets"] if p.get("selection")],
            "ligands_hidden": [l["object"] for l in e["ligands"] if l["object"] not in enabled],
            "boxes_hidden": [p["box_object"] for p in e["pockets"] if p.get("box_object") and p["box_object"] not in enabled],
        } for e in scene.get("structures") or []]
        if scene.get("alphafold"):
            out["alphafold"] = {"object": scene["alphafold"]["object"], "shown": scene["alphafold"]["object"] in enabled}
        out["session_opened_at"] = scene.get("created_at")
        if manifest.get("created_at") != scene.get("manifest_created_at"):
            out["warning"] = "prepare_structures rewrote the manifest after this session was opened: reopen it."
    return out
