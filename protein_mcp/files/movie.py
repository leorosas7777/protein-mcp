"""render_md_movie: an analysed trajectory as an .mp4 file (Phase 3, step 10c).

open_md_in_pymol plays a trajectory in a window on this machine. That window cannot be attached
to a report, sent to a reviewer or opened on a machine without PyMOL, and the PDF can only hold
a still. This module writes the same scene to a video file: PyMOL renders one image per state in
batch mode and the ffmpeg binary that ships with imageio-ffmpeg encodes them.

What it does NOT do, on purpose:
  * no network, and no model: the frames are the recorded trajectory drawn as it is;
  * no measuring. A video is not evidence. The result carries the equilibration control of the
    run exactly as analyze_md recorded it, because every frame is premature when it failed;
  * no file it has not checked. The topology and every DCD are matched against the SHA-256 in
    the manifest before PyMOL is started, by the same function the window session uses.

PyMOL has to run as a separate process in batch mode (`-cq`): rendering is a loop over states
that no single cmd call expresses, and the window session (files.viewer.launch) is detached and
never waited for. So this is the only place in the project that waits for PyMOL.

Imports common, files.prepare (prep_dir), files.viewer and imageio_ffmpeg. Never touches the
network.
"""
import asyncio
import json
import shutil
import time
from pathlib import Path

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import viewer
from protein_mcp.files.prepare import prep_dir

# Measured on this machine (2026-11-19, PyMOL 3.2 in batch, the real 3OG7 production trajectory,
# 8 768 atoms shown, ray-traced): 0.9 to 3.1 s per frame across four renders at 640x480 and
# 960x720. It is a range and not a figure, because the same settings gave 0.87 s and 2.53 s per
# frame depending on what else was running. A whole 100-frame trajectory is therefore 1.5 to 5
# minutes, which is why `every` defaults to auto: the default call mounts at most
# TARGET_FRAMES states so that it comes back in a time Claude Desktop is likely to wait for.
# Nothing is hidden by that: the states mounted and the sampling are in the result and in the
# .pml header, and every=1 renders all of them.
MAX_FRAMES = 240
TARGET_FRAMES = 60
SECONDS_PER_FRAME_MEASURED = 3.1
RENDER_TIMEOUT_SECONDS = 1800
DEFAULT_WIDTH, DEFAULT_HEIGHT, DEFAULT_FPS = 960, 720, 12
FRAME_GLOB = "frame_*.png"
_RENDERED = "protein-mcp: rendered"

# The frame loop, appended to the .pml as code (render_pml cannot check it, so it is a constant
# here and its arguments are written with ascii()). cmd.frame(n) selects state n because the
# scene ran mset 1 -N first; cmd.png with an explicit width and height does not need a window,
# which is why the render is ray-traced and not a screen grab.
_RENDER = r'''
import os as _os
def _pmcp_render(folder, first, last, width, height, ray):
    _os.makedirs(folder, exist_ok=True)
    done = 0
    for state in range(first, last + 1):
        cmd.frame(state)
        cmd.png(_os.path.join(folder, "frame_%05d.png" % state), width=width, height=height,
                ray=ray)
        done += 1
    print(" protein-mcp: rendered %d of %d frames" % (done, last - first + 1))
'''


def _tail(folder: Path, states: int, width: int, height: int, ray: int) -> list:
    lines = _RENDER.strip("\n").splitlines()
    lines.append("try:")
    lines.append(f"    _pmcp_render({ascii(str(folder))}, 1, {states}, {width}, {height}, {ray})")
    lines.append("except Exception as exc:")
    lines.append("    print(' protein-mcp: render failed: %s' % exc)")
    return lines


def ffmpeg_binary():
    """(path, version line) of the ffmpeg that ships with imageio-ffmpeg, or (None, the reason)."""
    try:
        import imageio_ffmpeg
    except ImportError as exc:
        return None, f"imageio-ffmpeg is not installed ({exc}): uv sync"
    try:
        return imageio_ffmpeg.get_ffmpeg_exe(), f"imageio-ffmpeg {imageio_ffmpeg.__version__}"
    except Exception as exc:                                   # it downloads on some platforms
        return None, f"imageio-ffmpeg has no ffmpeg binary available ({exc})"


async def _run(command, log_path: Path, timeout: float) -> dict:
    """Run a program to completion with its output in a file. Killed at the timeout."""
    started = time.monotonic()
    with open(log_path, "wb") as out:
        try:
            proc = await asyncio.create_subprocess_exec(
                *command, stdin=asyncio.subprocess.DEVNULL, stdout=out,
                stderr=asyncio.subprocess.STDOUT, cwd=str(log_path.parent))
        except OSError as exc:                 # the binary is not where it was said to be
            return {"ok": False, "error": f"could not start {command[0]}: {exc}",
                    "seconds": round(time.monotonic() - started, 1)}
        try:
            code = await asyncio.wait_for(proc.wait(), timeout)
        except asyncio.TimeoutError:
            proc.kill()
            await proc.wait()
            return {"ok": False, "timed_out": True, "seconds": round(time.monotonic() - started, 1),
                    "error": f"killed after {timeout:g} s"}
    return {"ok": code == 0, "exit_code": code, "seconds": round(time.monotonic() - started, 1)}


def _tail_lines(path: Path, n: int = 12) -> list:
    try:
        return [l.strip() for l in path.read_text(encoding="utf-8", errors="replace").splitlines()
                if l.strip()][-n:]
    except OSError:
        return []


# the window session and the video quote the same control from the same record: it lives in
# files.viewer so there is one reading of it and not two that could drift apart
_control = viewer.equilibration_control


async def render_md_movie(uniprot_id: str, structure: str = "", pocket: str = "", stage: str = "",
                          every: int = 0, width: int = DEFAULT_WIDTH, height: int = DEFAULT_HEIGHT,
                          fps: int = DEFAULT_FPS, keep_frames: bool = False) -> dict:
    """One analysed trajectory -> an .mp4 beside the reports of that accession."""
    accession, error = parse_accession(uniprot_id)
    if error:
        return error
    every, fps = max(0, int(every)), max(1, min(60, int(fps)))
    width, height = max(160, min(1920, int(width))), max(120, min(1440, int(height)))
    prep = prep_dir(accession)
    manifest, refusal = viewer.read_manifest(accession)
    if refusal:
        return refusal
    chosen, refusal = viewer.select_md_record(manifest, structure, pocket, stage)
    if refusal:
        return {**refusal, "uniprot_id": accession}
    (s, key, label, rec), listed = chosen
    run_name = f"{s.get('pdb_id')} {key} {label}"
    recorded = sum(e.get("frames") or 0 for e in rec.get("trajectories") or [])
    # every = 0 is "auto": one frame in N, so that the default call renders TARGET_FRAMES at
    # most and still spans the whole trajectory. Chosen here, reported below, never silent.
    sampling = "asked for"
    if not every:
        every = max(1, -(-recorded // TARGET_FRAMES)) if recorded else 1
        sampling = f"chosen here so that at most {TARGET_FRAMES} frames are rendered"

    ffmpeg, ffmpeg_note = ffmpeg_binary()
    if ffmpeg is None:
        return {"ok": False, "uniprot_id": accession, "run": run_name, "error": ffmpeg_note}
    pymol = viewer.find_pymol()
    if not pymol["found"]:
        return {"ok": False, "uniprot_id": accession, "run": run_name,
                "error": pymol["error"] + " A movie needs PyMOL on this machine: it is what draws "
                                          "the frames.",
                "how_to_fix": f"Set {viewer.PYMOL_ENV} in claude_desktop_config.json (env) to "
                              "pymol.exe or to a Python with pymol installed."}

    scene = viewer.build_trajectory_scene(prep, rec, every, MAX_FRAMES, play=False)
    if scene.get("error"):
        return {"ok": False, "uniprot_id": accession, "run": run_name, "error": scene["error"],
                "available": listed}
    states = scene["states"]
    if states > MAX_FRAMES:                       # build_trajectory_scene caps it; belt and braces
        return {"ok": False, "uniprot_id": accession, "run": run_name,
                "error": f"{states} frames is more than this tool renders in one call "
                         f"({MAX_FRAMES}, about {MAX_FRAMES * SECONDS_PER_FRAME_MEASURED:g} s at "
                         "the rate measured on this machine). Raise `every` to mount fewer."}

    folder = prep / "pymol"
    tag = f"{s.get('pdb_id')}_{key.replace(' ', '_')}_{label}"
    frames_dir = folder / f"frames_{tag}"
    pml = folder / f"movie_{tag}.pml"
    out = prep.parent / f"trajectory_{accession}_{tag}.mp4"
    header = [
        f"PyMOL batch render for {accession}, written by protein-mcp render_md_movie at "
        f"{utc_now()}.",
        f"Run: {run_name}, as recorded by analyze_md.",
        "Every file was checked against the SHA-256 analyze_md recorded before this was written.",
        f"Mounted: every {every} frame(s), {states} states; rendered at {width}x{height}.",
        "This script renders images and writes no video: ffmpeg does that afterwards.",
    ]
    try:
        frames_dir.mkdir(parents=True, exist_ok=True)
        for old in frames_dir.glob(FRAME_GLOB):   # a half-finished render must not reach the video
            old.unlink()
        # without this PyMOL ray-traces onto a transparent background and the frames carry an
        # alpha channel the encoder then flattens to whatever it likes (measured here on
        # 2026-11-19: the PNG opened white and the video came out black)
        ops = scene["ops"] + [("set", ["ray_opaque_background", 1])]
        pml.write_text(viewer.render_pml(ops, header,
                                         tail=_tail(frames_dir, states, width, height, 1)),
                       encoding="ascii", newline="\n")
    except OSError as exc:
        log.exception("Could not prepare the render of %s", run_name)
        return {"ok": False, "uniprot_id": accession, "run": run_name,
                "error": f"Could not write {pml}: {exc}"}

    render_log = folder / f"movie_{tag}.render.log"
    drawn = await _run(pymol["command"] + ["-cq", str(pml)], render_log,
                       RENDER_TIMEOUT_SECONDS)
    frames = sorted(frames_dir.glob(FRAME_GLOB))
    tail = _tail_lines(render_log)
    # PyMOL can stop half way and still exit 0, so three things are checked, not one: that it
    # ended well, that the frame loop reached its own end marker, and that there are as many
    # images as states. Without this a video of a third of the trajectory came back as ok with
    # the full state count beside it.
    finished = any(line.startswith(_RENDERED) for line in tail)
    if not frames or not drawn.get("ok") or not finished or len(frames) != states:
        if not frames:
            why = "PyMOL wrote no frame."
        elif not drawn.get("ok"):
            why = (f"PyMOL ended badly ({drawn.get('error') or 'exit code ' + str(drawn.get('exit_code'))})"
                   f" after {len(frames)} of {states} frames.")
        elif not finished:
            why = (f"PyMOL wrote {len(frames)} of {states} frames and never reached the end of "
                   "the frame loop.")
        else:
            why = f"PyMOL wrote {len(frames)} frames, not the {states} states it was asked for."
        return {"ok": False, "uniprot_id": accession, "run": run_name,
                "error": why + " No video was written: a video of part of a trajectory would "
                               "look like the whole of it.",
                "frames_written": len(frames), "states_asked_for": states,
                "frames_dir": str(frames_dir), "pymol": {**pymol, **drawn},
                "log_path": str(render_log), "log_tail": tail, "pml_path": str(pml)}

    encode_log = folder / f"movie_{tag}.ffmpeg.log"
    control = _control(rec)
    # the video leaves this machine; the chat message does not go with it, so the run, the
    # sampling and the equilibration control of that run are written into the file itself
    comment = (f"{run_name}, {accession}. {len(frames)} frames, one in {every} of the "
               f"{recorded or 'recorded'} frames analyze_md recorded, played at {fps} fps: "
               f"the seconds of this video are playback, not simulated time. "
               f"Equilibration control of this run: {control.get('status')}"
               + (f" (drift {control.get('drift_angstrom')} A against "
                  f"{control.get('threshold_angstrom')} A on {control.get('segid')})"
                  if control.get("drift_angstrom") is not None else "")
               + ". A video is not a measurement; the figures are in the experiment record.")
    encoded = await _run([ffmpeg, "-y", "-framerate", str(fps), "-start_number", "1",
                          "-i", str(frames_dir / "frame_%05d.png"),
                          # even dimensions, which yuv420p needs, without changing the view
                          "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                          "-c:v", "libx264", "-pix_fmt", "yuv420p",
                          "-metadata", f"title={run_name} ({accession})",
                          "-metadata", f"comment={comment}",
                          str(out)],
                         encode_log, RENDER_TIMEOUT_SECONDS)
    if not encoded["ok"] or not out.is_file():
        return {"ok": False, "uniprot_id": accession, "run": run_name,
                "error": "ffmpeg did not write the video.", "ffmpeg": {**encoded,
                                                                       "binary": ffmpeg_note},
                "frames_written": len(frames), "frames_dir": str(frames_dir),
                "log_path": str(encode_log), "log_tail": _tail_lines(encode_log)}

    kept = keep_frames
    if not keep_frames:
        try:
            shutil.rmtree(frames_dir)
        except OSError as exc:                    # not fatal: the video is already written
            log.warning("Could not remove %s (%s)", frames_dir, exc)
            kept = True

    return {
        "ok": True, "uniprot_id": accession, "run": f"{s.get('pdb_id')} pocket {key}, stage {label}",
        "file": str(out), "size_mb": round(out.stat().st_size / 1048576, 2),
        "frames": len(frames), "fps": fps, "seconds_long": round(len(frames) / fps, 1),
        "mounted": f"every {every} frame(s) of the recorded trajectory, {states} states",
        "sampling": sampling, "frames_in_the_record": recorded,
        "frames_recorded": scene["frames_recorded"],
        "resolution": f"{width}x{height}", "ray_traced": True,
        "render_seconds": drawn.get("seconds"), "encode_seconds": encoded.get("seconds"),
        "files_checked": scene["files_checked"], "files_refused": scene["files_refused"],
        "shows": scene["note"], "ligand_selection": scene["ligand_selection"],
        "equilibration_control": control,
        "embedded_in_the_file": comment,
        "frames_dir": str(frames_dir) if kept else None,
        "frames_removed": not kept,
        "pymol": {"command": pymol["command"], "source": pymol["source"],
                  "exit_code": drawn.get("exit_code"), "log_path": str(render_log)},
        "ffmpeg": {"binary": ffmpeg_note, "exit_code": encoded.get("exit_code"),
                   "log_path": str(encode_log)},
        "pml_path": str(pml),
        "available": listed,
        "not_a_measurement": "A video of a trajectory is not a measurement of it. Every figure "
                             "about this run is in the experiment record, and the equilibration "
                             "control above applies to everything in these frames: if it reads "
                             "NOT MET, what the video shows is premature, however convincing it "
                             "looks.",
    }


if __name__ == "__main__":  # python -m protein_mcp.files.movie P15056 3OG7 "032 C" production
    import sys
    args = sys.argv[1:]
    print(json.dumps(asyncio.run(render_md_movie(*args)), indent=1, ensure_ascii=False))
