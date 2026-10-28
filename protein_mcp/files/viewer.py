"""open_in_pymol: a PyMOL session (.pml) built from prep/manifest.json, launched only if PyMOL
is found (Phase 2, step 3). Without PyMOL the tool returns the path of the .pml.

The scene is a list of operations (function name + arguments), each a `cmd.<function>` call
from OPS, rendered by render_pml. The .pml also starts a small receiver (step 4b) that lets
files.control (pymol_view) send LIVE_OPS to the open window: only listed functions with data
arguments, never PyMOL text from the chat. The scene is saved as pymol/scene_<session>.json.

PyMOL defects found on Windows (2026-10-26, ESTADO_PROYECTO.md section 8.2): it reads a .pml in
the local code page (cp1252), its `load` command cannot open paths with non-ASCII characters,
and a BOM breaks the first line. So the .pml is pure ASCII and every call goes through a
`python` block with arguments written by ascii(). PyMOL writes to stdout: its output goes to a
log file, never to the MCP stdio. Without Qt it exits within seconds ("compile with --glut"),
so the process is checked after LAUNCH_CHECK_SECONDS.
Imports common and dossier (via files.prepare.prep_dir). Never touches the network.
"""
import asyncio
import json
import os
import re
import secrets
import shutil
import subprocess
from pathlib import Path

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files.prepare import prep_dir

PYMOL_ENV = "PROTEIN_MCP_PYMOL"
LAUNCH_CHECK_SECONDS = 5.0
LOG_TAIL_LINES = 8

# Every function a session may call (cmd.<name>). render_pml rejects anything else.
OPS = frozenset({
    "load", "load_cgo", "hide", "show", "color", "set_color", "select", "group", "enable",
    "disable", "align", "orient", "deselect", "bg_color",
})

# Functions the chat control (files.control, step 4b) may call in a RUNNING session; the receiver
# started by the .pml rejects anything else. Validated in PyMOL 3.2 on 2026-10-28: cmd.load runs
# .pml/.py/.pym files and cmd.label evaluates its expression as Python, so neither is here
# (selections are not Python: "__import__('os')" is an invalid selection). set_color evaluates a
# text RGB with PyMOL's safe_list_eval (the reviewer reported a sandbox escape, not reproduced
# here) and align takes a matrix file: neither is needed live (reset keeps the session's colours
# and AlphaFold stays where it was aligned), so both are out too.
LIVE_OPS = frozenset((OPS - {"load", "load_cgo", "set_color", "align"})
                     | {"zoom", "center", "get_names", "count_atoms"})
CONTROL_FILE = "control.json"
_LISTENING = "protein-mcp: control listening"

# CGO codes (pymol.cgo BEGIN, LINES, VERTEX, END, COLOR), read from PyMOL 3.2 on 2026-10-26.
_CGO_BEGIN, _CGO_LINES, _CGO_VERTEX, _CGO_END, _CGO_COLOR = 2.0, 1.0, 4.0, 3.0, 6.0

# AlphaFold DB pLDDT bands and colours, as on the AlphaFold DB site. Painted in this order, each
# over the previous one, so a value on a boundary (exactly 70) keeps the lower band.
PLDDT_BANDS = [
    ("plddt_very_low", [1.0, 0.490, 0.271], "all", "<= 50 (very low)"),
    ("plddt_low", [1.0, 0.859, 0.075], "b > 50", "50-70 (low)"),
    ("plddt_confident", [0.396, 0.796, 0.953], "b > 70", "70-90 (confident)"),
    ("plddt_very_high", [0.0, 0.325, 0.839], "b > 90", "> 90 (very high)"),
]

LEGEND = {
    "receptor": "cartoon, grey (one colour for every target chain); kept nucleotides/cofactors as sticks, metals as spheres",
    "ligands": "sticks, yellow carbons (crystal pose from ligands/*.sdf). Additives (ions, buffers) are not "
               "loaded, but small compounds the manifest classes as ligands are (has_pocket false). A ligand with on_target_chain false sits on a partner chain that is not in the "
               "receptor, so it floats next to it",
    "pockets": "side chains as sticks, light blue carbons; the cartoon is light blue there too (residues of "
               "the manifest pocket, within 5 Å)",
    "variant_position": "sticks, magenta carbons: the residue modelled at the variant position "
                        "(whatever it is: variant, wild type or other, see the manifest)",
    "boxes": "white lines: the box of each pocket (ligand + padding), in the mmCIF frame",
    "alphafold": "cartoon coloured by pLDDT: " + ", ".join(f"{label} {name}" for name, _, _, label in PLDDT_BANDS),
}

_CHAIN_RE = re.compile(r"[A-Za-z0-9]{1,4}")
_RESI_RE = re.compile(r"-?\d+[A-Za-z]?")
_OBSERVED_RE = re.compile(r"[A-Z0-9]{1,3} (-?\d+[A-Za-z]?)")      # clean.py: "GLU 600" (auth number)


def _name(*parts) -> str:
    """PyMOL object name: letters, digits and underscores only."""
    return re.sub(r"[^A-Za-z0-9_]", "_", "_".join(str(p) for p in parts if p))


def _resi(value: str) -> str:
    """auth_seq_id for a selection; a negative number needs a backslash in PyMOL."""
    return "\\" + value if value.startswith("-") else value


def _selection(obj: str, residues):
    """(`obj and ((chain A and resi 1+2) or ...)` or None, residues selected, "chain:resi" pairs
    dropped because PyMOL cannot select those IDs)."""
    by_chain, dropped = {}, []
    for chain, seq in residues:
        if _CHAIN_RE.fullmatch(str(chain or "")) and _RESI_RE.fullmatch(str(seq or "")):
            by_chain.setdefault(chain, {})[_resi(str(seq))] = None
        else:
            dropped.append(f"{chain}:{seq}")
    if not by_chain:
        return None, 0, dropped
    parts = [f"(chain {c} and resi {'+'.join(r)})" for c, r in by_chain.items()]
    return f"{obj} and ({' or '.join(parts)})", sum(len(r) for r in by_chain.values()), dropped


def _box_cgo(box: dict):
    """The 12 edges of the box as CGO lines, or None if center/size are not three numbers."""
    center, size = box.get("center"), box.get("size")
    if not (isinstance(center, list) and isinstance(size, list) and len(center) == len(size) == 3
            and all(isinstance(v, (int, float)) for v in center + size)):
        return None
    lo = [center[i] - size[i] / 2 for i in range(3)]
    hi = [center[i] + size[i] / 2 for i in range(3)]
    corner = [[(lo, hi)[(k >> i) & 1][i] for i in range(3)] for k in range(8)]
    obj = [_CGO_COLOR, 1.0, 1.0, 1.0, _CGO_BEGIN, _CGO_LINES]
    for a in range(8):
        for i in range(3):
            b = a | (1 << i)
            if b != a:
                obj += [_CGO_VERTEX, *corner[a], _CGO_VERTEX, *corner[b]]
    return obj + [_CGO_END]


def _file(prep: Path, rel):
    """Absolute path of a manifest file (relative to prep/), or None if it is outside prep/ or missing."""
    if not isinstance(rel, str) or not rel:
        return None
    path = (prep / rel).resolve()
    if prep.resolve() not in path.parents or not path.is_file():
        return None
    return path


def build_scene(manifest: dict, prep: Path, pdb_ids=None, include_alphafold=True) -> dict:
    """Operations and the description of what they show. Only the first structure is enabled."""
    ops, structures, skipped, notes = [], [], [], []
    wanted = [p.strip().upper() for p in pdb_ids] if pdb_ids else None
    listed = [s.get("pdb_id") for s in manifest.get("structures") or []]
    ops.append(("bg_color", ["black"]))
    ops.append(("hide", ["everything"]))
    for name, rgb, _, _ in PLDDT_BANDS:
        ops.append(("set_color", [name, rgb]))

    for s in manifest.get("structures") or []:
        pdb = s.get("pdb_id")
        if wanted is not None and pdb not in wanted:
            continue
        prepared = s.get("prepared") or {}
        receptor = _file(prep, (prepared.get("receptor") or {}).get("path")) if prepared.get("ok") else None
        if receptor is None:
            reason = (prepared.get("error") or "no cleaned receptor in the manifest") if not prepared.get("ok") \
                else "receptor file missing from prep/"
            skipped.append({"pdb_id": pdb, "reason": reason})
            continue
        rec = _name("rec", pdb)
        members, entry = [rec], {"pdb_id": pdb, "group": _name("s", pdb), "receptor": rec, "ligands": [],
                                 "ligands_not_loaded": [], "pockets": [], "pockets_not_shown": [],
                                 "variant_position": None}
        with_pocket = {(p.get("comp_id"), p.get("label_asym_id")) for p in prepared.get("pockets") or []
                       if not p.get("error")}
        # label_asym_id -> RCSB entity description of the non-target polymer chains (7ZR0: A, B = HSP90-beta)
        partner_of = {c.get("label_asym_id"): e.get("description") for e in s.get("other_polymer_entities") or []
                      for c in e.get("chains") or []}
        ops += [("load", [receptor.as_posix(), rec]), ("show", ["cartoon", rec]),
                ("color", ["grey80", f"{rec} and elem C"]),
                ("show", ["sticks", f"{rec} and organic"]), ("show", ["spheres", f"{rec} and inorganic"])]

        for inst in s.get("nonpolymer_instances") or []:
            record = inst.get("sdf") or {}
            sdf = _file(prep, record.get("path"))
            if sdf is None:
                if inst.get("category") != "additive":          # additives get no SDF by design
                    entry["ligands_not_loaded"].append({
                        "comp_id": inst.get("comp_id"), "label_asym_id": inst.get("label_asym_id"),
                        "reason": "SDF download failed" if record and not record.get("ok")
                        else "SDF file not found in prep/" if record else "no SDF in the manifest"})
                continue
            lig = _name("lig", pdb, inst.get("comp_id"), inst.get("label_asym_id"))
            members.append(lig)
            entry["ligands"].append({"object": lig, "comp_id": inst.get("comp_id"), "name": inst.get("name"),
                                     "auth_asym_id": inst.get("auth_asym_id"),
                                     # False: RCSB lists no target-chain neighbour (7ZR0 ATP sits on HSP90)
                                     "on_target_chain": bool(inst.get("neighbor_target_chains")),
                                     "has_pocket": (inst.get("comp_id"), inst.get("label_asym_id")) in with_pocket,
                                     # partner (non-target) chains it touches, named only from the manifest
                                     "partner_chains": [
                                         f"{(v or {}).get('auth_asym_id') or k}: "
                                         f"{partner_of.get(k) or 'no description in the manifest'}"
                                         for k, v in sorted((inst.get("neighbor_other_chains") or {}).items())]})
            ops += [("load", [sdf.as_posix(), lig]), ("show", ["sticks", lig]),
                    ("color", ["yellow", f"{lig} and elem C"])]

        for pocket in prepared.get("pockets") or []:
            if pocket.get("error"):
                entry["pockets_not_shown"].append({"comp_id": pocket.get("comp_id"),
                                                   "label_asym_id": pocket.get("label_asym_id"),
                                                   "reason": f"pocket error in the manifest: {pocket['error']}"})
                continue
            tag = (pdb, pocket.get("comp_id"), pocket.get("label_asym_id"))
            sel, n_sel, dropped = _selection(
                rec, [(r.get("auth_asym_id"), r.get("auth_seq_id")) for r in pocket.get("residues") or []])
            row = {"comp_id": pocket.get("comp_id"), "label_asym_id": pocket.get("label_asym_id"),
                   "auth_asym_id": pocket.get("auth_asym_id"), "n_residues": len(pocket.get("residues") or []),
                   "n_residues_selected": n_sel, "selection": None, "box_object": None}
            if dropped:
                row["residues_not_selectable"] = dropped
            if sel:
                row["selection"] = _name("pocket", *tag)
                ops += [("select", [row["selection"], sel]),
                        ("show", ["sticks", f"{row['selection']} and (sidechain or name CA)"]),
                        ("color", ["lightblue", f"{row['selection']} and elem C"])]
            cgo = _box_cgo(pocket.get("box") or {})
            if cgo:
                row["box_object"] = _name("box", *tag)
                members.append(row["box_object"])
                ops.append(("load_cgo", [cgo, row["box_object"]]))
            entry["pockets"].append(row)

        positions, by_chain = [], []
        for chain in (prepared.get("metrics") or {}).get("chains") or []:
            check = chain.get("variant_check") or {}
            if not check:
                continue
            by_chain.append({"auth_asym_id": chain.get("auth_asym_id"), "observed": check.get("observed"),
                             "status": check.get("status")})
            m = _OBSERVED_RE.fullmatch(check.get("observed") or "")
            if m:
                positions.append((chain.get("auth_asym_id"), m.group(1)))
        sel, _, _ = _selection(rec, positions)
        if sel:
            entry["variant_position"] = {"selection": _name("variant_pos", pdb),
                                         "chains": [c for c, _ in positions], "by_chain": by_chain}
            ops += [("select", [entry["variant_position"]["selection"], sel]),
                    ("show", ["sticks", f"{entry['variant_position']['selection']} and not name N+C+O"]),
                    ("color", ["magenta", f"{entry['variant_position']['selection']} and elem C"])]

        ops.append(("group", [entry["group"], " ".join(members)]))
        structures.append(entry)

    first = structures[0] if structures else None
    for entry in structures[1:]:
        ops.append(("disable", [entry["group"]]))
    for i, entry in enumerate(structures):
        entry["enabled"] = i == 0

    alphafold = None
    af = manifest.get("alphafold") or {}
    af_path = _file(prep, (af.get("file") or {}).get("path")) if include_alphafold else None
    if af_path:
        obj = _name("af", manifest.get("uniprot_id"))
        ops += [("load", [af_path.as_posix(), obj]), ("show", ["cartoon", obj])]
        for name, _, sel, _ in PLDDT_BANDS:
            ops.append(("color", [name, f"{obj} and ({sel})"]))
        alphafold = {"object": obj, "entry_id": af.get("entry_id"), "enabled": False, "aligned_to": None}
        if first:
            ops.append(("align", [obj, first["receptor"]]))
            alphafold["aligned_to"] = first["receptor"]
            notes.append(f"{obj} is superposed on {first['receptor']} in the session only (PyMOL align); "
                         "the files are unchanged. Other structures and the boxes keep their mmCIF frames.")
        ops.append(("disable", [obj]))
    elif include_alphafold:
        notes.append("No AlphaFold file in the manifest (or it is missing from prep/): not loaded.")

    if first:
        focus = next((p["selection"] for p in first["pockets"] if p["selection"]), first["receptor"])
        ops.append(("orient", [focus]))
    ops.append(("deselect", []))

    off = [f"{e['pdb_id']} {l['comp_id']} {l['object'].rsplit('_', 1)[-1]}"
           + (f" (on {'; '.join(l['partner_chains'])})" if l["partner_chains"] else "")
           for e in structures for l in e["ligands"] if not l["on_target_chain"]]
    if off:
        notes.append("Loaded but not on a target chain (RCSB lists no target neighbour; the partner chain is "
                     "not in the receptor, so they float next to it): " + ", ".join(off) + ".")
    if len(structures) > 1:
        notes.append("Structures are not superposed: each keeps the frame of its mmCIF (the frame of its "
                     "boxes). Only the first one is shown; the others are enabled from the object panel.")
    return {
        "ops": ops, "structures": structures, "alphafold": alphafold, "skipped": skipped, "notes": notes,
        "pdb_ids_not_in_manifest": [p for p in wanted if p not in listed] if wanted else [],
    }


# One failed call (e.g. an unreadable file) must not stop the rest: a python block ends at the
# first exception. The error goes to the PyMOL console.
_RUNNER = [
    "from pymol import cmd",
    "def _do(fn, *args):",
    "    try:",
    "        getattr(cmd, fn)(*args)",
    "    except Exception as exc:",
    "        print(' protein-mcp: cmd.%s%r failed: %s' % (fn, args[:2], exc))",
]
_FAILED = "protein-mcp: cmd."
_LOADED = "protein-mcp: session loaded"

# Receiver for the chat control (step 4b), run inside PyMOL after the scene. Standard library only:
# one socket on 127.0.0.1 (ephemeral port), a random token, one JSON request per connection
# {"token", "ops": [[function, args]]}; only LIVE_OPS with text/number arguments (lists of numbers
# are also accepted, none of LIVE_OPS needs one). Port, token and session id go to pymol/control.json. Validated
# 2026-10-28 (ESTADO_PROYECTO.md 8.2): cmd calls from these threads work and the window keeps
# redrawing; the receiver dies with PyMOL. The token is created here, at run time: never in the .pml.
_RECEIVER = r'''
import json as _json, os as _os, secrets as _secrets, socket as _socket, threading as _threading
def _pmcp_receiver(control, session, live):
    token = _secrets.token_hex(16)
    srv = _socket.socket(_socket.AF_INET, _socket.SOCK_STREAM)
    srv.bind(("127.0.0.1", 0))
    srv.listen(8)
    def scalar(a):
        return isinstance(a, (str, int, float)) and not isinstance(a, bool)
    def handle(conn):
        with conn:
            f = conn.makefile("rwb")
            def reply(obj):
                obj["session"] = session
                f.write((_json.dumps(obj) + "\n").encode()); f.flush()
            try:
                req = _json.loads(f.readline(1000000))
            except Exception as exc:
                return reply({"ok": False, "error": "bad request: %s" % exc})
            if not isinstance(req, dict) or not _secrets.compare_digest(str(req.get("token")), token):
                return reply({"ok": False, "error": "bad token"})
            results = []
            for op in req.get("ops") or []:
                fn, args = (op + [None, None])[:2] if isinstance(op, list) else (None, None)
                if fn not in live or not isinstance(args, list) or not all(
                        scalar(a) or (isinstance(a, list) and a and all(scalar(x) and not isinstance(x, str) for x in a))
                        for a in args):
                    results.append({"fn": fn, "ok": False, "error": "not allowed"})
                    continue
                try:
                    r = getattr(cmd, fn)(*args)
                    results.append({"fn": fn, "ok": True,
                                    "result": r if isinstance(r, (str, int, float, list, type(None))) else repr(r)})
                except Exception as exc:
                    results.append({"fn": fn, "ok": False, "error": str(exc).strip()})
            reply({"ok": True, "results": results})
    def serve():
        while True:
            conn, _ = srv.accept()
            _threading.Thread(target=handle, args=(conn,), daemon=True).start()
    _threading.Thread(target=serve, daemon=True).start()
    with open(control + ".part", "w") as out:
        _json.dump({"port": srv.getsockname()[1], "token": token, "session": session, "pid": _os.getpid()}, out)
    _os.replace(control + ".part", control)
    print(" protein-mcp: control listening")
'''


def render_pml(ops, header, control=None) -> str:
    """ASCII .pml: comment header + one python block with a _do(<function>, ...) call per operation.
    control = (control.json path, session id) also starts the chat-control receiver."""
    lines = [f"# {line}".encode("ascii", "backslashreplace").decode() for line in header]
    lines += ["python", *_RUNNER]
    for fn, args in ops:
        if fn not in OPS:
            raise ValueError(f"PyMOL function not allowed: {fn!r}")
        lines.append(f"_do({', '.join(ascii(a) for a in [fn, *args])})")
    if control:
        lines += _RECEIVER.strip("\n").splitlines()
        lines.append("try:")
        lines.append(f"    _pmcp_receiver({ascii(str(control[0]))}, {ascii(control[1])}, {ascii(sorted(LIVE_OPS))})")
        lines.append("except Exception as exc:")
        lines.append("    print(' protein-mcp: control not started: %s' % exc)")
    lines += [f"print(' {_LOADED}')", "python end"]
    return "\n".join(lines) + "\n"


def find_pymol() -> dict:
    """The command that starts PyMOL: PROTEIN_MCP_PYMOL (a pymol executable, or a Python with
    pymol installed: then `-m pymol`), else `pymol` on PATH."""
    env = os.environ.get(PYMOL_ENV, "").strip().strip('"')
    if env:
        path = Path(env)
        if not path.is_file():
            return {"found": False, "source": PYMOL_ENV,
                    "error": f"{PYMOL_ENV} points to {env}, which is not a file."}
        if path.name.lower().startswith("python"):
            # python.exe (a venv launcher) opens an empty console window next to PyMOL; pythonw.exe of
            # the same environment does not, and its output still goes to launch.log (2026-10-28)
            windowless = path.with_name("pythonw.exe")
            exe = windowless if path.name.lower() == "python.exe" and windowless.is_file() else path
            command = [str(exe), "-m", "pymol"]
        else:
            command = [str(path)]
        return {"found": True, "source": PYMOL_ENV, "command": command}
    found = shutil.which("pymol")
    if found:
        return {"found": True, "source": "PATH", "command": [found]}
    return {"found": False, "source": None,
            "error": f"PyMOL not found: {PYMOL_ENV} is not set and there is no `pymol` on PATH."}


def _log_lines(path: Path) -> list[str]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    return [line for line in text.splitlines() if line.strip()]


async def launch(command, pml: Path, log_path: Path) -> dict:
    """Start PyMOL detached (outlives the server), output to log_path; report an early exit."""
    kwargs = {"start_new_session": True}
    if os.name == "nt":
        kwargs = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS}
    # unbuffered, so the log already holds the session's messages when it is read below
    kwargs["env"] = {**os.environ, "PYTHONUNBUFFERED": "1"}
    try:
        with open(log_path, "wb") as out:
            proc = subprocess.Popen(command + [str(pml)], stdin=subprocess.DEVNULL, stdout=out,
                                    stderr=subprocess.STDOUT, cwd=pml.parent, close_fds=True, **kwargs)
    except OSError as exc:
        return {"launched": False, "error": f"Could not start PyMOL: {exc}"}
    await asyncio.sleep(LAUNCH_CHECK_SECONDS)
    code = proc.poll()
    lines = _log_lines(log_path)
    # the session's own messages start with "protein-mcp:"; PyMOL's echo of the .pml source does not
    errors = [line.strip() for line in lines
              if line.strip().startswith((_FAILED, "protein-mcp: control not started"))]
    if code is None:
        loaded = any(line.strip() == _LOADED for line in lines)
        out = {"launched": True, "pid": proc.pid, "session_loaded": loaded, "session_errors": errors,
               "control_listening": any(line.strip() == _LISTENING for line in lines),
               "note": f"PyMOL is still running after {LAUNCH_CHECK_SECONDS:g} s; the server cannot see the "
                       "window itself."}
        if not loaded:
            out["note"] += (" The session had not finished loading at that point: later errors are only "
                            "in log_path.")
        return out
    tail = lines[-LOG_TAIL_LINES:]
    out = {"launched": False, "exit_code": code, "log_tail": tail,
           "error": f"PyMOL exited after less than {LAUNCH_CHECK_SECONDS:g} s (exit code {code})."}
    if any("compile with --glut" in t or "Qt not available" in t for t in tail):
        out["hint"] = ("This PyMOL has no Qt: install PyQt5 or PySide6 in the same Python environment "
                       "(pip/uv pip install PyQt5).")
    return out


async def open_in_pymol(uniprot_id: str, pdb_ids=None, include_alphafold=True, start=True) -> dict:
    accession, error = parse_accession(uniprot_id)
    if error:
        return error
    prep = prep_dir(accession)
    manifest_path = prep / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {"ok": False, "error": f"No preparation manifest for {accession}: run prepare_structures first."}
    except (OSError, ValueError) as exc:
        log.exception("Could not read %s", manifest_path)
        return {"ok": False, "error": f"Could not read {manifest_path}: {exc}"}

    scene = build_scene(manifest, prep, pdb_ids, include_alphafold)
    if not scene["structures"] and not scene["alphafold"]:
        return {"ok": False, "error": "Nothing to show: no requested structure has a cleaned receptor.",
                "skipped": scene["skipped"], "pdb_ids_not_in_manifest": scene["pdb_ids_not_in_manifest"]}

    folder = prep / "pymol"
    pml = folder / "session.pml"
    # what this session holds, for pymol_view (files.control): one file per session, so a later
    # call that only writes a new .pml does not change what the open window is matched against
    session = secrets.token_hex(8)
    scene_file = folder / f"scene_{session}.json"
    scene_record = {
        "session_id": session, "created_at": utc_now(), "manifest_created_at": manifest.get("created_at"),
        "structures": scene["structures"], "alphafold": scene["alphafold"],
        # "reset" replays the scene without loading anything
        "reset_ops": [[fn, args] for fn, args in scene["ops"] if fn in LIVE_OPS],
    }
    header = [
        f"PyMOL session for {accession}, written by protein-mcp open_in_pymol at {utc_now()}.",
        f"Built from prep/manifest.json created {manifest.get('created_at')}.",
        "Open it in PyMOL with File > Run Script..., or: pymol session.pml",
        "Every file is loaded with cmd.load in a python block: PyMOL cannot load non-ASCII paths otherwise.",
    ]
    try:
        folder.mkdir(parents=True, exist_ok=True)
        scene_file.write_text(json.dumps(scene_record, ensure_ascii=False, indent=1), encoding="utf-8")
        pml.write_text(render_pml(scene["ops"], header, (folder / CONTROL_FILE, session)),
                       encoding="ascii", newline="\n")
    except OSError as exc:
        log.exception("Could not write %s", pml)
        return {"ok": False, "error": f"Could not write {pml}: {exc}"}

    pymol = find_pymol()
    if start and pymol["found"]:
        pymol.update(await launch(pymol["command"], pml, folder / "launch.log"))
        pymol["log_path"] = str(folder / "launch.log")
    elif not start:
        pymol["launched"] = False
        pymol["note"] = "start=False: the session file was written, PyMOL was not started."
    else:
        pymol["launched"] = False

    shown = [e["pdb_id"] for e in scene["structures"]]
    return {
        "ok": True, "uniprot_id": accession, "pml_path": str(pml),
        "selection": f"{len(shown)} of {len(manifest.get('structures') or [])} structures in the manifest"
                     + (" + the AlphaFold model" if scene["alphafold"] else ""),
        "structures": scene["structures"], "alphafold": scene["alphafold"], "skipped": scene["skipped"],
        "pdb_ids_not_in_manifest": scene["pdb_ids_not_in_manifest"], "legend": LEGEND, "notes": scene["notes"],
        "pymol": pymol,
        "control": {"session_id": session,
                    "note": "While this PyMOL window is open, pymol_view can change its view from the chat "
                            "(started by the .pml, also when it is opened by hand with File > Run Script)."},
        "how_to_open": "In PyMOL: File > Run Script... and choose pml_path. To let the server start PyMOL, set "
                       f"{PYMOL_ENV} in claude_desktop_config.json (env) to pymol.exe or to a Python with "
                       "pymol and PyQt5 installed.",
    }
