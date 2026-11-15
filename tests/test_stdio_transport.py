"""The server must answer a tool call that loads a native library without waiting for another message.

Found 2026-11-14 while diagnosing an analyze_md call that never answered in Claude Desktop: the
MCP stdio transport keeps a blocking read pending on stdin in a worker thread, and on Windows the
first import of scipy.linalg (pulled in by MDAnalysis) waits until that read ends, i.e. until the
client sends something else. `server.detach_stdin` gives the transport a duplicate of stdin and
points fd 0 at NUL.

The test drives the real server over stdio with the MCP client and calls analyze_md on the
synthetic trajectory of test_md_analyze.py (no uniprot_id: nothing is recorded). The control runs
the same server WITHOUT detach_stdin and must NOT answer within CONTROL_WAIT seconds; if it does,
this machine does not reproduce the hang and the first check proves nothing, so the test fails.

Windows only. Run with: uv run python tests/test_stdio_transport.py   (about 30-60 s)
"""
import asyncio
import json
import os
import sys
import tempfile
import time
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

sys.path.insert(0, str(Path(__file__).parent))
from test_md_analyze import _run_topology, write_fixture  # noqa: E402

FIXED = "from protein_mcp.server import main; main()"
UNFIXED = "from protein_mcp.server import mcp; mcp.run()"
ANSWER_WITHIN = 60.0   # generous: the call itself takes a few seconds on the fixture
CONTROL_WAIT = 20.0    # the unfixed server must still be silent after this long


async def _call(code: str, folder: Path, timeout: float, nudge: bool = False):
    """Seconds until analyze_md answers on a fresh server, or None if it did not within `timeout`.

    With `nudge`, after a timeout the client writes one more message (a ping) and the third value
    is how long the pending call then took to answer, or None: that is what shows the hang waits
    on stdin and not on something else.
    """
    params = StdioServerParameters(command=sys.executable, args=["-c", code], env=dict(os.environ))
    outcome = {}
    try:
        with open(os.devnull, "w") as errlog:
            async with stdio_client(params, errlog=errlog) as (r, w):
                async with ClientSession(r, w) as session:
                    await session.initialize()
                    t = time.perf_counter()
                    task = asyncio.ensure_future(
                        session.call_tool("analyze_md", {"folder": str(folder)}))
                    done, _ = await asyncio.wait({task}, timeout=timeout)
                    if done:
                        outcome["result"] = (time.perf_counter() - t, task.result(), None)
                    elif not nudge:
                        task.cancel()
                        outcome["result"] = (None, None, None)
                    else:
                        t2 = time.perf_counter()
                        ping = asyncio.ensure_future(session.send_ping())
                        done, _ = await asyncio.wait({task}, timeout=ANSWER_WITHIN)
                        after = time.perf_counter() - t2 if done else None
                        for f in (task, ping):
                            f.cancel()
                        outcome["result"] = (None, None, after)
    except BaseExceptionGroup:
        # Closing the client's end of stdin releases a hung server, whose late answer then
        # arrives on a closed session. Only tolerated after a timeout was already recorded.
        if outcome.get("result", (0,))[0] is not None:
            raise
    return outcome["result"]


def test_the_server_answers_without_another_message():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        write_fixture(folder)
        _run_topology(folder)

        took, res, _ = asyncio.run(_call(FIXED, folder, ANSWER_WITHIN))
        assert took is not None, f"the server did not answer within {ANSWER_WITHIN} s"
        body = json.loads(res.content[0].text)
        assert not res.isError and body.get("ok") is True, res.content[0].text[:500]

        silent, _, after_nudge = asyncio.run(_call(UNFIXED, folder, CONTROL_WAIT, nudge=True))
        assert silent is None, (
            f"control: the server WITHOUT detach_stdin answered in {silent:.1f} s, so this "
            "machine does not reproduce the hang and the check above proves nothing")
        assert after_nudge is not None, (
            "control: the unfixed server stayed silent even after another message arrived, so "
            "it hangs for some other reason and this test does not isolate the stdin cause")
    print(f"ok  analyze_md answered in {took:.1f} s; without detach_stdin no answer in "
          f"{CONTROL_WAIT:.0f} s, and the pending call answered {after_nudge:.1f} s after one "
          "more message reached stdin (the control fails as it must, for the stated reason)")


def main():
    if sys.platform != "win32":
        print("skipped: the hang is specific to Windows")
        return
    test_the_server_answers_without_another_message()
    print("ALL OK")


if __name__ == "__main__":
    main()
