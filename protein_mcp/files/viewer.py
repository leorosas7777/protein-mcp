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
import hashlib
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
    # trajectory sessions (step 10b). load_traj mounts a DCD on a topology already loaded;
    # create/delete build the small object and drop the solvent; mset/mplay/set drive the movie;
    # intra_fit removes the whole-system drift so the ligand's own motion is what is seen.
    "load_traj", "create", "delete", "mset", "mplay", "set", "intra_fit", "zoom",
    # a PSF carries no secondary structure, so a trajectory draws as loops until dss assigns it
    "dss",
})

# NOT derived from OPS any more. LIVE_OPS is the surface the chat control can reach in a RUNNING
# window, and it was security-probed as a fixed list on 2026-10-28; adding a trajectory function
# to OPS must not silently widen it. Anything new stays out until it is probed the same way.
_LIVE_EXCLUDED = {"load", "load_cgo", "set_color", "align",
                  "load_traj", "create", "delete", "mset", "mplay", "set", "intra_fit", "dss"}

# Functions the chat control (files.control, step 4b) may call in a RUNNING session; the receiver
# started by the .pml rejects anything else. Validated in PyMOL 3.2 on 2026-10-28: cmd.load runs
# .pml/.py/.pym files and cmd.label evaluates its expression as Python, so neither is here
# (selections are not Python: "__import__('os')" is an invalid selection). set_color evaluates a
# text RGB with PyMOL's safe_list_eval (the reviewer reported a sandbox escape, not reproduced
# here) and align takes a matrix file: neither is needed live (reset keeps the session's colours
# and AlphaFold stays where it was aligned), so both are out too.
LIVE_OPS = frozenset((OPS - _LIVE_EXCLUDED) | {"zoom", "center", "get_names", "count_atoms"})
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


def render_pml(ops, header, control=None, tail=None) -> str:
    """ASCII .pml: comment header + one python block with a _do(<function>, ...) call per operation.
    control = (control.json path, session id) also starts the chat-control receiver.
    tail = extra python lines appended after the scene, for a batch session that has to do
    something no single cmd call expresses (the frame loop of the renderer). They are code, not
    data, so a caller builds them from a constant template with ascii()-written arguments; this
    function does not and cannot check them."""
    # PyMOL splits a line on ";" BEFORE deciding it is a comment, so a semicolon in the header
    # is run as a command ("Error: ambiguous command", measured here 2026-11-19). It is a
    # comment, so the character is simply dropped rather than the sentence rewritten.
    lines = [f"# {line}".replace(";", ",").encode("ascii", "backslashreplace").decode()
             for line in header]
    lines += ["python", *_RUNNER]
    for fn, args in ops:
        if fn not in OPS:
            raise ValueError(f"PyMOL function not allowed: {fn!r}")
        lines.append(f"_do({', '.join(ascii(a) for a in [fn, *args])})")
    if tail:
        lines += list(tail)
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


# CHARMM-GUI names the ligand it builds `segi HETA` / `resn LIG`, measured on the real topology
# of this project (3OG7 032 C: 52 atoms, 33 of them heavy). It is a convention of the builder, not
# of the file format, so the session states the assumption and still shows the protein if it
# matches nothing.
LIGAND_SELECTION = "segi HETA or resn LIG"
MAX_TRAJ_STATES = 200          # states actually mounted; a cap, not a statement about the run


def _hash_problem(path: Path, recorded):
    """None when the file still has the SHA-256 recorded for it, else why it cannot be used.

    A file with NO recorded hash and a file whose hash changed are refused alike, but they are
    not the same fact and the reason must not say the second when it is the first (the reviewer
    found this same invention on 2026-11-12 with md_job.psf).
    """
    if not recorded:
        return "no SHA-256 was recorded for it, so there is nothing to check it against"
    if hashlib.sha256(path.read_bytes()).hexdigest() != recorded:
        return "no longer has the SHA-256 that was recorded for it"
    return None


def md_records(manifest: dict) -> list:
    """(structure, pocket key, stage label, record) for every analysis analyze_md recorded."""
    out = []
    for s in manifest.get("structures") or []:
        analyses = s.get("md_analysis")
        if not isinstance(analyses, dict):
            continue
        for key, stages in sorted(analyses.items()):
            if not isinstance(stages, dict):
                continue
            for label, rec in sorted(stages.items()):
                if isinstance(rec, dict) and rec.get("trajectories"):
                    out.append((s, key, label, rec))
    return out


def build_trajectory_scene(prep: Path, rec: dict, every: int, max_states: int,
                           play: bool = True) -> dict:
    """Operations for one analysed trajectory, and what they will show.

    The files come from the md_analysis block analyze_md wrote, and each one is checked against
    the SHA-256 recorded there BEFORE anything is loaded: a trajectory that changed on disk is
    not the trajectory that was measured, and opening it would show something no figure in the
    experiment record describes.

    The session keeps protein + ligand and deletes the full system, because the solvent is most
    of the file and nothing here looks at it: on this project's own run that is 8 768 atoms
    instead of 139 989.
    """
    topology = (rec.get("topology") or {}).get("file")
    psf = _file(prep, topology)
    if psf is None:
        return {"error": f"the topology {topology or 'recorded for this run'} is missing from prep/"}
    problem = _hash_problem(psf, (rec.get("topology") or {}).get("sha256"))
    if problem:
        return {"error": f"the topology {topology} {problem}, so it cannot be shown as the "
                         "topology that was measured"}
    checked = [{"file": topology, "role": "topology"}]
    problems = []

    ops = [("bg_color", ["black"]), ("hide", ["everything"]),
           ("load", [psf.as_posix(), "system"])]
    frames, mounted = 0, 0
    for entry in rec.get("trajectories") or []:
        path = _file(prep, entry.get("file"))
        if path is None:
            problems.append({"file": entry.get("file"), "reason": "missing from prep/"})
            continue
        problem = _hash_problem(path, entry.get("sha256"))
        if problem:
            problems.append({"file": entry.get("file"), "reason": problem})
            continue
        # the cap is what is LEFT, not the whole cap again: cmd.load_traj applies `max` to its
        # own call, so passing the full cap per file would mount up to max states per file
        room = int(max_states) - mounted
        if room <= 0:
            problems.append({"file": entry.get("file"),
                             "reason": f"the cap of {max_states} states was already reached"})
            continue
        # positional, in cmd.load_traj's own order: filename, object, state, format, interval,
        # average, start, stop, max (read from PyMOL 3.2.0a with inspect, not assumed).
        # state=0 APPENDS after the last state; state=1 overwrites from the first one, so with
        # a production written in several DCDs only the last file survived while the tool
        # reported the sum (measured in PyMOL 3.2 on 2026-11-19: two loads at state=1 left 3
        # states, at state=0 left 6).
        ops.append(("load_traj", [path.as_posix(), "system", 0, "", int(every), 1, 1, -1, room]))
        checked.append({"file": entry.get("file"), "role": "trajectory",
                        "frames_recorded": entry.get("frames")})
        frames += entry.get("frames") or 0
        # every `every`-th frame is mounted, so a file of N frames gives ceil(N / every) states,
        # never N: counting them as N would make the movie claim frames that were not loaded
        n = entry.get("frames")
        mounted += min(room, -(-int(n) // int(every)) if n else room)
    if not any(fn == "load_traj" for fn, _ in ops):
        return {"error": "no trajectory of this run can be opened: "
                         + "; ".join(f"{p['file']}: {p['reason']}" for p in problems)}

    keep = f"system and (polymer or ({LIGAND_SELECTION}))"
    states = max(1, min(int(max_states), mounted or int(max_states)))
    ops += [
        ("create", ["traj", keep, 0, 0]),
        ("delete", ["system"]),
        # the topology has no secondary structure (a PSF carries none), so without this the
        # protein draws as loops everywhere. It is PyMOL's own assignment from the coordinates
        # of the state shown, and the result says so: it is not read from any source.
        ("dss", ["traj"]),
        ("hide", ["everything"]),
        ("show", ["cartoon", "traj and polymer"]),
        ("color", ["grey70", "traj and polymer and elem C"]),
        ("show", ["sticks", f"traj and ({LIGAND_SELECTION})"]),
        ("color", ["yellow", f"traj and ({LIGAND_SELECTION}) and elem C"]),
        ("hide", ["everything", "traj and hydro"]),
        # every state is fitted on the protein backbone of the first one, so what moves on screen
        # is the ligand inside the pocket and not the whole box drifting through space
        ("intra_fit", ["traj and polymer and name CA", 1]),
        ("orient", [f"traj and ({LIGAND_SELECTION})"]),
        ("zoom", [f"traj and ({LIGAND_SELECTION})", 8.0]),
        ("set", ["movie_fps", 12]),
        # frame N is state N from here on, which is what the renderer steps through
        ("mset", [f"1 -{states}"]),
    ]
    if play:
        ops.append(("mplay", []))
    return {"ops": ops, "files_checked": checked, "files_refused": problems,
            "frames_recorded": frames, "states": states,
            "ligand_selection": LIGAND_SELECTION,
            "secondary_structure": "assigned by PyMOL (dss) from the coordinates: a PSF carries "
                                   "none, and nothing here reads it from a source",
            "note": "Protein and ligand only: the solvent and the ions come in with the topology "
                    "and are dropped straight after, because nothing here looks at them. Every "
                    "state is superposed on the backbone of the first, and the cartoon's "
                    "secondary structure is PyMOL's own assignment from the coordinates."}


def read_manifest(accession: str):
    """(manifest, None) or (None, the refusal to return). One reader for every session tool."""
    prep = prep_dir(accession)
    try:
        return json.loads((prep / "manifest.json").read_text(encoding="utf-8")), None
    except FileNotFoundError:
        return None, {"ok": False, "error": f"No preparation manifest for {accession}: run "
                                            "prepare_structures first."}
    except (OSError, ValueError) as exc:
        log.exception("Could not read the manifest of %s", accession)
        return None, {"ok": False, "error": f"Could not read the manifest: {exc}"}


def select_md_record(manifest: dict, structure: str, pocket: str, stage: str):
    """((structure, pocket key, stage label, record), every run listed) or (None, the refusal).

    Shared by the tool that opens a trajectory and the one that renders it: which run is chosen
    changes what the user sees, so an ambiguous request is refused with the list, never guessed.
    """
    records = md_records(manifest)
    listed = [f"{s.get('pdb_id')} | {key} | {label}" for s, key, label, _r in records]
    if not records:
        return None, {"ok": False, "available": [],
                      "error": "No analysed trajectory is recorded for this accession. The chain "
                               "is build_md_package, CHARMM-GUI, build_md_job, the cluster, and "
                               "then analyze_md with uniprot_id, which records the analysis this "
                               "tool reads."}
    wanted = [(s, key, label, rec) for s, key, label, rec in records
              if (not structure or (s.get("pdb_id") or "").upper() == structure.strip().upper())
              and (not pocket or key.upper() == pocket.strip().upper())
              and (not stage or label.lower() == stage.strip().lower())]
    if not wanted:
        return None, {"ok": False, "available": listed,
                      "error": "No recorded trajectory matches that structure, pocket and stage."}
    if len(wanted) > 1:
        return None, {"ok": False,
                      "available": [f"{s.get('pdb_id')} | {key} | {label}"
                                    for s, key, label, _r in wanted],
                      "error": f"{len(wanted)} recorded trajectories match. Name one with "
                               "structure, pocket and stage."}
    return (wanted[0], listed), None


async def open_md_in_pymol(uniprot_id: str, structure: str = "", pocket: str = "",
                           stage: str = "", every: int = 1, max_states: int = MAX_TRAJ_STATES,
                           start: bool = True) -> dict:
    """Open an analysed trajectory in PyMOL and play it. Only a run analyze_md already recorded."""
    accession, error = parse_accession(uniprot_id)
    if error:
        return error
    prep = prep_dir(accession)
    manifest, refusal = read_manifest(accession)
    if refusal:
        return refusal

    chosen, refusal = select_md_record(manifest, structure, pocket, stage)
    if refusal:
        return {**refusal, "uniprot_id": accession}
    (s, key, label, rec), listed = chosen

    scene = build_trajectory_scene(prep, rec, max(1, int(every)), max(1, int(max_states)))
    if scene.get("error"):
        return {"ok": False, "uniprot_id": accession, "error": scene["error"],
                "run": f"{s.get('pdb_id')} {key} {label}", "available": listed}

    folder = prep / "pymol"
    pml = folder / f"trajectory_{s.get('pdb_id')}_{key.replace(' ', '_')}_{label}.pml"
    header = [
        f"PyMOL trajectory session for {accession}, written by protein-mcp open_md_in_pymol "
        f"at {utc_now()}.",
        f"Run: {s.get('pdb_id')} pocket {key}, stage {label}, as recorded by analyze_md.",
        "Every file was checked against the SHA-256 analyze_md recorded before this was written.",
        f"Ligand selection: {scene['ligand_selection']} (CHARMM-GUI's own naming).",
        f"Mounted: every {every} frame(s), at most {max_states} states.",
        "The movie plays protein + ligand only; the solvent is dropped after loading.",
    ]
    try:
        folder.mkdir(parents=True, exist_ok=True)
        pml.write_text(render_pml(scene["ops"], header), encoding="ascii", newline="\n")
    except OSError as exc:
        log.exception("Could not write %s", pml)
        return {"ok": False, "error": f"Could not write {pml}: {exc}"}

    pymol = find_pymol()
    if start and pymol["found"]:
        pymol.update(await launch(pymol["command"], pml, folder / "launch.log"))
        pymol["log_path"] = str(folder / "launch.log")
    else:
        pymol["launched"] = False
        if not start:
            pymol["note"] = "start=False: the session file was written, PyMOL was not started."

    return {
        "ok": True, "uniprot_id": accession,
        "run": f"{s.get('pdb_id')} pocket {key}, stage {label}",
        "pml_path": str(pml), "pymol": pymol,
        "files_checked": scene["files_checked"], "files_refused": scene["files_refused"],
        "frames_recorded": scene["frames_recorded"], "states_mounted": scene["states"],
        "mounted": f"every {every} frame(s), at most {max_states} states",
        "shows": scene["note"], "ligand_selection": scene["ligand_selection"],
        "available": listed,
        "not_a_measurement": "This is the trajectory itself, not a measurement of it. Every figure "
                             "about this run is in the experiment record; watching a pose move is "
                             "not evidence that it is right, and the equilibration control of that "
                             "run still applies to everything seen here.",
    }


# One colour per docked compound. Picked to be told apart on the black background and NOT to
# carry any meaning: the colour says which compound, never how good it is.
POSE_COLORS = ["yellow", "cyan", "magenta", "salmon", "palegreen", "orange", "slate", "wheat",
               "deepteal", "lightpink"]


def docking_blocks(manifest: dict) -> list:
    """(structure, pocket key, block) for every pocket that has at least one docking run."""
    out = []
    for s in manifest.get("structures") or []:
        blocks = s.get("docking")
        if not isinstance(blocks, dict):
            continue
        for key, block in sorted(blocks.items()):
            if isinstance(block, dict) and (block.get("runs") or {}):
                out.append((s, key, block))
    return out


def build_poses_scene(manifest: dict, prep: Path, s: dict, key: str, block: dict,
                      compounds=None) -> dict:
    """The receptor of one pocket with the top pose of each of its docking runs, side by side.

    Every pose file is checked against the SHA-256 the docking run recorded before it is loaded:
    a file that changed is not the pose Vina wrote, and a picture of it would illustrate a score
    that no longer belongs to it.

    The whole SDF of a run is loaded, so the object holds every pose Vina kept and PyMOL shows
    state 1, which is the top-scoring one -- the pose the experiment record is about. The object
    panel's state slider walks the others, and the result says so rather than hiding them.
    """
    pdb = s.get("pdb_id")
    prepared = s.get("prepared") or {}
    receptor_rec = prepared.get("receptor") or {}
    receptor = _file(prep, receptor_rec.get("path"))
    if receptor is None:
        return {"error": f"the cleaned receptor of {pdb} is not in prep/: run prepare_structures"}

    rec = _name("rec", pdb)
    ops = [("bg_color", ["black"]), ("hide", ["everything"]),
           ("load", [receptor.as_posix(), rec]), ("show", ["cartoon", rec]),
           ("color", ["grey70", f"{rec} and elem C"]),
           ("show", ["spheres", f"{rec} and inorganic"])]
    members, runs, refused = [rec], [], []

    wanted = {c.strip().upper() for c in compounds} if compounds else None
    seen = []
    for i, (run_id, run) in enumerate(sorted((block.get("runs") or {}).items())):
        ligand = run.get("ligand") or {}
        cd = run.get("compound_docking") or {}
        rd = run.get("redocking") or {}
        cid = ligand.get("compound_id") or ligand.get("comp_id") or ""
        # the ChEMBL id when the record has no name (CHEMBL525191 has none), never the run id
        # dressed up as a compound
        name = ligand.get("name") or cid or run_id
        seen.append(run_id)
        if wanted is not None and not ({str(cid).upper(), str(name).upper(), run_id.upper()}
                                       & wanted):
            continue
        sdf_rec = next((f for n, f in (run.get("files") or {}).items() if n.endswith(".sdf")), None)
        if not sdf_rec or not sdf_rec.get("path"):
            refused.append({"run": run_id, "reason": "no pose SDF is recorded for this run"})
            continue
        path = _file(prep, sdf_rec["path"])
        if path is None:
            refused.append({"run": run_id, "file": sdf_rec["path"],
                            "reason": "missing from prep/"})
            continue
        problem = _hash_problem(path, sdf_rec.get("sha256"))
        if problem:
            refused.append({"run": run_id, "file": sdf_rec["path"], "reason": problem})
            continue
        obj = _name("pose", cid or run_id)
        colour = POSE_COLORS[len(runs) % len(POSE_COLORS)]
        members.append(obj)
        ops += [("load", [path.as_posix(), obj]), ("show", ["sticks", obj]),
                ("color", [colour, f"{obj} and elem C"]),
                ("hide", ["everything", f"{obj} and hydro"])]
        runs.append({
            "run": run_id, "object": obj, "color": colour,
            "compound": name, "compound_id": ligand.get("compound_id"),
            "level": ligand.get("level"),
            "kind": "compound" if cd else "redocking" if rd else "unknown",
            # a re-docking has no compound_docking block: its top-pose score is the score of
            # pose 1 and its heavy-atom count is on the ligand. Reading only compound_docking
            # printed "no record" for a score the manifest holds.
            "top_pose_score_kcal_mol": cd.get("top_pose_score_kcal_mol") if cd else
                                       ((run.get("poses") or [{}])[0].get("vina_score_kcal_mol")),
            "ligand_efficiency": cd.get("ligand_efficiency_kcal_mol_per_heavy_atom"),
            "heavy_atoms": cd.get("heavy_atoms") or ligand.get("heavy_atoms"),
            "box_limited": bool(cd.get("box_limited")),
            "rmsd_to_crystal_pose_angstrom": rd.get("top_pose_rmsd_angstrom"),
            "rmsd_to_reference_pose_angstrom": (run.get("cross_docking") or {})
                                               .get("top_pose_rmsd_angstrom"),
            "poses_in_file": len(run.get("poses") or []),
            "file": sdf_rec["path"],
        })
    if not runs:
        return {"error": f"no pose of {pdb} pocket {key} can be opened",
                "files_refused": refused, "runs_in_the_manifest": seen}

    pocket = next((p for p in prepared.get("pockets") or []
                   if " ".join(x for x in (p.get("comp_id"), p.get("label_asym_id")) if x) == key),
                  None)
    pocket_row = None
    if pocket:
        sel, n_sel, dropped = _selection(
            rec, [(r.get("auth_asym_id"), r.get("auth_seq_id")) for r in pocket.get("residues") or []])
        pocket_row = {"residues_in_the_manifest": len(pocket.get("residues") or []),
                      "residues_selected": n_sel, "residues_not_selectable": dropped,
                      "selection": None, "box_object": None}
        if sel:
            pocket_row["selection"] = _name("pocket", pdb, key.replace(" ", "_"))
            ops += [("select", [pocket_row["selection"], sel]),
                    ("show", ["sticks", f"{pocket_row['selection']} and (sidechain or name CA)"]),
                    ("color", ["grey50", f"{pocket_row['selection']} and elem C"])]
        cgo = _box_cgo(block.get("box") or pocket.get("box") or {})
        if cgo:
            pocket_row["box_object"] = _name("box", pdb, key.replace(" ", "_"))
            members.append(pocket_row["box_object"])
            ops.append(("load_cgo", [cgo, pocket_row["box_object"]]))

    focus = " or ".join(r["object"] for r in runs)
    ops += [("group", [_name("poses", pdb, key.replace(" ", "_")), " ".join(members)]),
            ("orient", [focus]), ("zoom", [focus, 4.0]), ("deselect", [])]
    return {
        "ops": ops, "runs": runs, "files_refused": refused, "receptor_object": rec,
        "receptor_file": receptor_rec.get("path"),
        "receptor_sha256_recorded": bool(receptor_rec.get("sha256")),
        "pocket": pocket_row, "runs_in_the_manifest": seen,
        "note": "One receptor and one object per docking run, all in the frame of this crystal. "
                "State 1 of each object is the top-scoring pose; the other poses Vina kept are "
                "in the same object and the state slider walks them.",
    }


async def open_poses_in_pymol(uniprot_id: str, structure: str = "", pocket: str = "",
                              compounds=None, start: bool = True) -> dict:
    """Open the docked poses of one pocket together in PyMOL, one colour per compound."""
    accession, error = parse_accession(uniprot_id)
    if error:
        return error
    prep = prep_dir(accession)
    manifest, refusal = read_manifest(accession)
    if refusal:
        return refusal

    blocks = docking_blocks(manifest)
    listed = [f"{s.get('pdb_id')} | {key} | {len(block.get('runs') or {})} run(s)"
              for s, key, block in blocks]
    if not blocks:
        return {"ok": False, "uniprot_id": accession, "available": [],
                "error": "No docking run is recorded for this accession: run dock_ligand first."}
    wanted = [(s, key, block) for s, key, block in blocks
              if (not structure or (s.get("pdb_id") or "").upper() == structure.strip().upper())
              and (not pocket or key.upper() == pocket.strip().upper())]
    if not wanted:
        return {"ok": False, "uniprot_id": accession, "available": listed,
                "error": "No docked pocket matches that structure and pocket."}
    if len(wanted) > 1:
        # the poses of two pockets are in different frames: showing them together would put two
        # receptors on screen as if they were superposed, which they are not
        return {"ok": False, "uniprot_id": accession,
                "available": [f"{s.get('pdb_id')} | {key}" for s, key, _b in wanted],
                "error": f"{len(wanted)} docked pockets match. Name one with structure and pocket."}
    s, key, block = wanted[0]

    scene = build_poses_scene(manifest, prep, s, key, block, compounds)
    if scene.get("error"):
        return {"ok": False, "uniprot_id": accession, "error": scene["error"],
                "files_refused": scene.get("files_refused", []),
                "runs_in_the_manifest": scene.get("runs_in_the_manifest", []), "available": listed}

    folder = prep / "pymol"
    pml = folder / f"poses_{s.get('pdb_id')}_{key.replace(' ', '_')}.pml"
    header = [
        f"PyMOL pose session for {accession}, written by protein-mcp open_poses_in_pymol "
        f"at {utc_now()}.",
        f"Pocket: {s.get('pdb_id')} {key}. {len(scene['runs'])} of "
        f"{len(scene['runs_in_the_manifest'])} docking run(s) of this pocket.",
        "Every pose file was checked against the SHA-256 its docking run recorded.",
        "State 1 of each object is the top-scoring pose of that run.",
        "The colour says which compound; it says nothing about the score.",
    ]
    try:
        folder.mkdir(parents=True, exist_ok=True)
        pml.write_text(render_pml(scene["ops"], header), encoding="ascii", newline="\n")
    except OSError as exc:
        log.exception("Could not write %s", pml)
        return {"ok": False, "error": f"Could not write {pml}: {exc}"}

    pymol = find_pymol()
    if start and pymol["found"]:
        pymol.update(await launch(pymol["command"], pml, folder / "launch.log"))
        pymol["log_path"] = str(folder / "launch.log")
    else:
        pymol["launched"] = False
        if not start:
            pymol["note"] = "start=False: the session file was written, PyMOL was not started."

    return {
        "ok": True, "uniprot_id": accession, "pocket": f"{s.get('pdb_id')} {key}",
        "selection": f"{len(scene['runs'])} of {len(scene['runs_in_the_manifest'])} docking runs "
                     f"of this pocket",
        "pml_path": str(pml), "pymol": pymol,
        "runs": scene["runs"], "files_refused": scene["files_refused"],
        "pocket_residues": scene["pocket"], "receptor": {
            "object": scene["receptor_object"], "file": scene["receptor_file"],
            "sha256_recorded": scene["receptor_sha256_recorded"],
            "note": "The cleaned receptor carries no SHA-256 in the manifest, so it is the only "
                    "file here that is not checked before it is loaded."
            if not scene["receptor_sha256_recorded"] else None},
        "shows": scene["note"], "available": listed,
        "colors_repeat": len(scene["runs"]) > len(POSE_COLORS),
        "not_a_measurement": "Seeing where a compound sits is not evidence that it binds there. "
                             "The scores beside each pose are calculated values of this project's "
                             "docking runs, they are not affinities, and the set has no negative "
                             "control, so nothing here orders the compounds.",
    }


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
