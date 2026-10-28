"""Offline test for the chat control of PyMOL (step 4b): the receiver written into session.pml by
open_in_pymol and pymol_view (protein_mcp/files/control.py).
The manifest is the one of test_choice (real 3OG7, 5JSM, 8C7Y, AlphaFold files and RCSB chain
replies from tests/fixtures/prep/), with the real mmCIFs in raw/ and their real SHA-256. The
.pml python block is run with a recording stand-in for pymol.cmd: the real receiver code then
listens on 127.0.0.1 and every pymol_view call goes through a real socket.
Real PyMOL (2026-10-28, ESTADO_PROYECTO.md 8.2): cmd.load runs .pml/.py files and cmd.label
evaluates Python, so the receiver must refuse both.
Run: uv run python tests/test_control.py
"""
import gzip
import hashlib
import json
import os
import shutil
import socket
import sys
import tempfile
import threading
import time
import types
from pathlib import Path

from protein_mcp import dossier
from protein_mcp.files import choice, control, viewer

sys.path.insert(0, str(Path(__file__).parent))
from test_choice import ACC, FIX, build_manifest  # noqa: E402


class FakeCmd:
    """Records every call; answers get_names like PyMOL (enabled objects)."""
    def __init__(self):
        self.calls, self.lock = [], threading.Lock()
        self.disabled = set()

    def __getattr__(self, name):
        def call(*args):
            with self.lock:
                self.calls.append((name, args))
            if name == "disable":
                self.disabled.add(args[0])
            elif name == "enable":
                self.disabled.discard(args[0])
            elif name == "get_names":
                return ["rec_3OG7", "s_8C7Y", "s_3OG7", "box_3OG7_032_C", "lig_3OG7_032_C", "af_P15056"]
            elif name == "color" and args[1] == "nonexistent":
                raise Exception(' Error: Invalid selection name "nonexistent".')
            return None
        return call


def make_prep():
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())
    prep = dossier.DOSSIER_DIR / ACC / "prep"
    m = build_manifest(prep / "clean")
    (prep / "raw").mkdir(parents=True)
    (prep / "ligands").mkdir()
    shutil.copy(FIX / "modelserver_3OG7_C.sdf", prep / "ligands" / "3OG7_032_C.sdf")
    for s in m["structures"]:
        if s["file"].get("ok"):
            data = gzip.decompress((FIX / f"{s['pdb_id']}.cif.gz").read_bytes())
            (prep / s["file"]["path"]).write_bytes(data)
            s["file"]["sha256"] = hashlib.sha256(data).hexdigest()
    af = gzip.decompress((FIX / "AF-P15056-F1-model_v6.cif.gz").read_bytes())
    (prep / m["alphafold"]["file"]["path"]).write_bytes(af)
    m["alphafold"]["file"]["sha256"] = hashlib.sha256(af).hexdigest()
    m["alphafold"]["residues"] = "1-766"                 # real forms: text range, FASTA length
    m["sequence"] = {"ok": True, "length": 766}
    choice.write_manifest(prep / "manifest.json", m)
    return prep, m


def run_pml(pml: Path) -> FakeCmd:
    """Run the python block of session.pml with a stand-in pymol.cmd (as PyMOL would)."""
    lines = pml.read_text(encoding="ascii").splitlines()
    block = "\n".join(lines[lines.index("python") + 1: lines.index("python end")])
    fake = FakeCmd()
    sys.modules["pymol"] = types.SimpleNamespace(cmd=fake)
    exec(compile(block, str(pml), "exec"), {"__name__": "pymol_session"})
    return fake


def test_control():
    prep, m = make_prep()
    view = control.pymol_view
    r = view(ACC, "status")
    assert not r["ok"] and "open it with open_in_pymol first" in r["error"]
    assert "action must be one of" in view(ACC, "label")["error"]

    import asyncio
    opened = asyncio.run(viewer.open_in_pymol(ACC, start=False))
    assert opened["ok"], opened
    session = opened["control"]["session_id"]
    pml = Path(opened["pml_path"])
    text = pml.read_text(encoding="ascii")
    assert "_pmcp_receiver(" in text and "token_hex(16)" in text
    scene = json.loads((pml.parent / f"scene_{session}.json").read_text(encoding="utf-8"))
    assert scene["manifest_created_at"] == m["created_at"]
    assert all(fn in viewer.LIVE_OPS for fn, _ in scene["reset_ops"]) and not any(fn == "load" for fn, _ in scene["reset_ops"])
    assert "load" not in viewer.LIVE_OPS and "load_cgo" not in viewer.LIVE_OPS and "label" not in viewer.LIVE_OPS

    fake = run_pml(pml)
    ctl = json.loads((pml.parent / viewer.CONTROL_FILE).read_text())
    assert ctl["session"] == session and ctl["token"] not in text           # the token is made at run time
    assert any(c[0] == "load" for c in fake.calls)                          # the scene itself still loads
    fake.calls.clear()

    # receiver: closed list, data arguments only, token, session
    raw = control.send(ctl, [["load", ["evil.pml"]], ["label", ["all", "__import__('os')"]],
                             ["color", ["red", {"a": 1}]], ["color", ["red", True]],
                             ["set_color", ["c1", [1.0, 0.0, 0.0]]], ["color", ["red", "nonexistent"]], "bad"])
    # set_color and align are out of LIVE_OPS too (reviewer 2026-10-28: text RGB is evaluated; matrix file)
    assert [x["ok"] for x in raw["results"]] == [False] * 7, raw
    assert raw["results"][5]["error"].startswith("Error: Invalid selection")
    assert [c[0] for c in fake.calls] == ["color"]
    assert "set_color" not in viewer.LIVE_OPS and "align" not in viewer.LIVE_OPS
    assert control.send(dict(ctl, token="x" * 32), [["get_names", []]])["error"] == "bad token"
    assert "another session" in control.send(dict(ctl, session="other"), [["get_names", []]])["error"]
    fake.calls.clear()

    # actions
    r = view(ACC, "show_structure", structure="8c7y")
    assert r["ok"] and r["applied"] == r["sent"] and r["shown"] == "8C7Y" and "cannot see the window" in r["note"]
    assert fake.calls[-2:] == [("enable", ("s_8C7Y",)), ("orient", ("rec_8C7Y",))]
    assert ("disable", ("s_3OG7",)) in fake.calls and ("disable", ("af_P15056",)) in fake.calls
    r = view(ACC, "focus_pocket", structure="8C7Y", pocket="TXV")
    assert not r["ok"] and "matches 2 pockets" in r["error"] and r["options"]["pockets"] == ["TXV U", "TXV X"]
    r = view(ACC, "focus_pocket", structure="8C7Y", pocket="8C7Y TXV U")
    assert r["ok"] and fake.calls[-1] == ("zoom", ("pocket_8C7Y_TXV_U", 5)) and ("enable", ("box_8C7Y_TXV_U",)) in fake.calls
    assert not view(ACC, "focus_pocket", structure="AlphaFold", pocket="032")["ok"]
    r = view(ACC, "color_chain", structure="3OG7", chain="B", color="Red")
    assert r["ok"] and fake.calls[-1] == ("color", ("red", "rec_3OG7 and chain B and elem C"))
    assert view(ACC, "color_chain", structure="3OG7", chain="b", color="red")["options"] == {"chains": ["A", "B"]}
    assert "chartreuse" in view(ACC, "color_chain", structure="3OG7", chain="B", color="chartreuse")["error"]
    assert "name the chain" in view(ACC, "color_chain", structure="3OG7", color="red")["error"]
    assert "not in this PyMOL session" in view(ACC, "show_structure", structure="7ZR0")["error"]    # not prepared

    # residues: UniProt -> auth from the real file; 600 is not modelled in 3OG7 chain A (597-614 missing)
    r = view(ACC, "show_residues", structure="3OG7", chain="A", residues="594-596, 600")
    assert r["ok"] and r["residues"] == ["594: ASP (auth 594)", "595: PHE (auth 595)", "596: GLY (auth 596)"]
    assert r["not_modelled"] == [600]
    assert ("select", ("res_3OG7", "rec_3OG7 and ((chain A and resi 594+595+596))")) in fake.calls
    assert ("color", ("green", "res_3OG7 and elem C")) in fake.calls and fake.calls[-1] == ("deselect", ())
    r = view(ACC, "show_residues", structure="3OG7", chain="B", residues="600")
    assert r["residues"] == ["600: GLU (auth 600)"]
    r = view(ACC, "show_residues", structure="AlphaFold", residues="600", color="cyan")
    assert r["ok"] and r["residues"] == ["600: VAL (auth 600)"]
    assert ("select", ("res_AlphaFold", "af_P15056 and ((resi 600))")) in fake.calls
    assert "none of the residues" in view(ACC, "show_residues", structure="3OG7", chain="A", residues="600-605")["error"]
    # an AlphaFold fragment (numbering may not start at 1) is refused, not guessed
    frag = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    frag["alphafold"]["residues"] = "1-700"
    choice.write_manifest(prep / "manifest.json", frag)
    assert "whole canonical" in view(ACC, "show_residues", structure="AlphaFold", residues="600")["error"]
    choice.write_manifest(prep / "manifest.json", m)
    for bad in ("abc", "600-500", "1-500", ""):
        assert not view(ACC, "show_residues", structure="3OG7", chain="A", residues=bad)["ok"], bad

    # hide/show/reset/status
    r = view(ACC, "hide", target="boxes", structure="3OG7")
    assert r["ok"] and r["objects"] == ["box_3OG7_032_C"]
    assert view(ACC, "hide", target="ligands")["objects"] == ["lig_3OG7_032_C"]       # only 3OG7 has an SDF here
    assert not view(ACC, "hide", target="surface")["ok"]
    fake.calls.clear()
    r = view(ACC, "reset")
    assert r["ok"] and ("enable", ("box_3OG7_032_C",)) in fake.calls and ("enable", ("lig_3OG7_032_C",)) in fake.calls
    # the scene's hide everything also hides the cgo of the boxes: shown again after it (reviewer 2026-10-28)
    hide_all = max(i for i, c in enumerate(fake.calls) if c == ("hide", ("everything",)))
    assert ("show", ("cgo", "box_3OG7_032_C")) in fake.calls[hide_all:]
    fake.calls.clear()
    assert view(ACC, "show", target="boxes", structure="3OG7")["ok"]
    assert fake.calls == [("enable", ("box_3OG7_032_C",)), ("show", ("cgo", "box_3OG7_032_C"))]
    assert not any(c[0] in ("load", "load_cgo") for c in fake.calls)
    r = view(ACC, "status")
    s = {x["pdb_id"]: x for x in r["structures"]}
    assert s["8C7Y"]["shown"] and s["3OG7"]["shown"] and not s["5JSM"]["shown"] and r["alphafold"]["shown"]
    assert s["8C7Y"]["pockets"] == ["TXV U", "TXV X"] and s["3OG7"]["boxes_hidden"] == []
    assert "enabled_objects" not in r

    # the raw file changed after the manifest: refused (no numbering from an unverified file)
    raw3 = prep / "raw" / "3OG7.cif"
    raw3.write_bytes(raw3.read_bytes() + b"\n")
    assert "SHA-256" in view(ACC, "show_residues", structure="3OG7", chain="A", residues="594")["error"]
    # prepare_structures rewrote the manifest: chains and residues refused, status warns
    m2 = json.loads((prep / "manifest.json").read_text(encoding="utf-8"))
    m2["created_at"] = "2026-10-29T00:00:00Z"
    choice.write_manifest(prep / "manifest.json", m2)
    assert "reopen it" in view(ACC, "color_chain", structure="3OG7", chain="B", color="red")["error"]
    assert "warning" in view(ACC, "status")

    # sent but no answer in time (reviewer 2026-10-28): NOT reported as closed
    slow = socket.socket()
    slow.bind(("127.0.0.1", 0))
    slow.listen(1)
    got = []

    def answer_late():
        conn, _ = slow.accept()
        got.append(conn.makefile("rb").readline())
        time.sleep(1.5)
        conn.close()
    threading.Thread(target=answer_late, daemon=True).start()
    old_timeout, control.TIMEOUT_SECONDS = control.TIMEOUT_SECONDS, 0.5
    r = control.send(dict(ctl, port=slow.getsockname()[1]), [["get_names", []]])
    control.TIMEOUT_SECONDS = old_timeout
    assert not r["ok"] and "not_open" not in r and "may or may not have been applied" in r["error"] and got, r

    # PyMOL closed: control.json stays, nothing listens on its port
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    closed_port = s.getsockname()[1]
    s.close()
    (pml.parent / viewer.CONTROL_FILE).write_text(json.dumps(dict(ctl, port=closed_port)))
    r = view(ACC, "status")
    assert not r["ok"] and r["pymol_open"] is False and "still be loading" not in r["error"]
    time.sleep(0.05)
    os.utime(pml)                                       # a newer session.pml: maybe still loading
    assert "may still be loading" in view(ACC, "status")["error"]


def test_find_pythonw():
    """python.exe (venv launcher) -> pythonw.exe of the same folder: no extra console window."""
    folder = Path(tempfile.mkdtemp())
    for name in ("python.exe", "pythonw.exe"):
        (folder / name).write_bytes(b"")
    old = os.environ.get(viewer.PYMOL_ENV)
    try:
        os.environ[viewer.PYMOL_ENV] = str(folder / "python.exe")
        assert viewer.find_pymol()["command"] == [str(folder / "pythonw.exe"), "-m", "pymol"]
        (folder / "pythonw.exe").unlink()
        assert viewer.find_pymol()["command"] == [str(folder / "python.exe"), "-m", "pymol"]
    finally:
        if old is None:
            os.environ.pop(viewer.PYMOL_ENV, None)
        else:
            os.environ[viewer.PYMOL_ENV] = old


def main():
    test_control()
    test_find_pythonw()
    print("ALL OK")


main()
