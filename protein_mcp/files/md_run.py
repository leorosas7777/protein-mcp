"""build_md_job: turn a verified CHARMM-GUI archive into a SLURM job for any cluster
(Phase 3, step 5).

The server never reaches the cluster: this writes files. The user copies the folder over and
runs sbatch, exactly as the CHARMM-GUI upload is manual.

Two kinds of value meet here and are never mixed:

1. THE SIMULATION, read from the archive, never chosen here. The run order and the chunk scheme
   come from the archive's own namd/README (a csh script: equilibration once, then `cntmax`
   production chunks, each chained to the previous one by the same `sed` renaming); the steps,
   the timestep and the rest from the NAMD inputs. The chunk length comes from the `run` line,
   which is what NAMD executes, not from `numsteps`. If that README is not the shape this code
   knows, the job is refused rather than guessed.

2. THE CLUSTER, about which this code knows NOTHING. No partition, no account, no module and no
   time limit is built in. Called without `site_report`, the tool writes no script and returns
   the commands for the user to run on their own machine (`sinfo`, `sacctmgr`, `module avail`,
   `scontrol show config`); the pasted output is then parsed strictly, each field is proposed
   with the command it came from, the user may override any of it, and anything unreadable
   becomes a CONFIRM slot that stops the script. The pasted text is data: nothing executes it.

Nothing is predicted anywhere. `--time` is the limit the cluster reports, used as a CEILING, and
the way to replace it is to measure the first chunk with sacct. No wall time per nanosecond is
estimated for any machine.

The archive is re-checked here with check_charmm_gui_output's own measurements, and a job is
refused if they fail: two of the four jobs of the step-4 validation finished with a broken
system (ESTADO 8.4), so a submission script for an unchecked archive would hand the cluster a
broken system.

Pure function + CLI: python -m protein_mcp.files.md_run build P15056 --archive charmm-gui.tgz
Imports common, files.choice, files.dock (manifest resolution) and files.md_package (archive
reading and the output checks).
"""
import argparse
import json
import re
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice, dock, md_package

# extraction caps: the archive is untrusted input (45 MB / 212 members in the validation job)
MAX_MEMBERS = 5000
MAX_TOTAL_BYTES = 4 * 1024 ** 3
GPU_RESIDENT_LINE = "CUDASOAintegrate on"


class MdRunError(Exception):
    def __init__(self, message: str, **extra):
        super().__init__(message)
        self.extra = extra


# ---------------------------------------------------------------- the site
# Nothing about any cluster is known here. Every site value is either READ from the output the
# user pastes (and then carries the command it came from), CHOSEN by the user, or left as a
# CONFIRM slot that stops the script. The four marks used in the script:
#   REPORTED  read from the pasted output of a command, which is named
#   CHOSEN    the user's decision (a choice the cluster does not make for them)
#   CONFIRM   not known: the script refuses to run until it is filled in
#   EDITABLE  a knob whose default comes from the archive or from a stated rule
# The first two print `|`-delimited columns on purpose: a line of free text, a note or a module
# separator cannot then be mistaken for a row, and the order of the columns is fixed by the
# command instead of being guessed from what the output looks like.
SINFO = 'sinfo -o "%P|%l|%G|%D|%c|%m"'
SACCTMGR = 'sacctmgr -n -P show assoc where user=$USER format=account,partition,qos'
SITE_COMMANDS = [
    {"command": SINFO,
     "gives": "the partitions, their time limit, their GPUs (GRES), and the cores and memory of "
              "a node",
     "fills": "--partition, --time (its ceiling), --gres/--gpus, --cpus-per-task, --mem"},
    {"command": SACCTMGR,
     "gives": "the accounts you may charge, if the site uses them",
     "fills": "--account (left out when the site reports none)"},
    {"command": 'module avail 2>&1 | grep -i namd ; echo "---" ; which namd3 namd2 2>/dev/null',
     "gives": "the NAMD modules, and any NAMD already on your PATH",
     "fills": "the module to load and the binary to run; a NAMD you compiled yourself is given "
              "as a path instead"},
    {"command": 'scontrol show config | grep -iE "MaxArraySize|GresTypes"',
     "gives": "the largest job array the site allows and whether it knows about GPUs as a GRES",
     "fills": "--array (refused if it does not fit); GresTypes is reported next to the GPU "
              "request so you can see whether this site knows `gpu` at all"},
]
HOW_TO_MEASURE_TIME = ("sacct -j <jobid> --format=JobID,JobName,Elapsed,State,MaxRSS -X")

# a value that stops the script until it is filled in
CONFIRM_MODULE = "CONFIRM_NAMD_MODULE"
CONFIRM_BINARY = "CONFIRM_NAMD_BINARY"
CONFIRM_PARTITION = "CONFIRM_PARTITION"
CONFIRM_ACCOUNT = "CONFIRM_ACCOUNT"

# From the user's own working template (2026-11-09): a multicore-CUDA NAMD 3 build is run
# DIRECTLY, with no charmrun, taking its thread count from the allocation. The archive's README
# instead assumes `charmrun namd2 +p4`, a CPU build. Which one a given cluster provides is
# unknown, so this is a default the script says to check, not a fact about anywhere.
DEFAULT_NAMD_OPTS = "+p${SLURM_CPUS_PER_TASK:-1} +setcpuaffinity +idlepoll"
NAMD_OPTS_SOURCE = ("a multicore-CUDA NAMD 3 build, run directly with no charmrun; threads from "
                    "the allocation. Check it against the build your cluster provides")
# with no reported limit to go by, a chunked run does not need more: it restarts chunk by chunk
FALLBACK_TIME = "24:00:00"


def site_questions(reason: str) -> dict:
    """What to ask the user to paste, and why. Returned instead of guessing a cluster."""
    return {
        "site_unknown": reason,
        "paste_the_output_of": SITE_COMMANDS,
        "then": "call build_md_job again with site_report set to everything that came back "
                "(all four outputs in one string is fine; they are told apart by their shape).",
        "note": "Nothing is run on your machine or on the cluster by this tool. The output you "
                "paste is read as data: it is parsed, echoed back, and whatever cannot be parsed "
                "is reported instead of guessed.",
        "if_you_compiled_namd_yourself": "say the full path to the binary (namd_binary=...); no "
                                         "module is then loaded.",
    }


def _time_to_seconds(value: str):
    """SLURM time limits: 1-00:00:00, 24:00:00, 6:00:00, 30:00, infinite. None if unreadable."""
    value = (value or "").strip().lower()
    if value in ("infinite", "unlimited", "n/a", ""):
        return None
    days, dash, rest = value.partition("-")
    if not dash:
        days, rest = "0", value
    try:
        parts = [int(p) for p in rest.split(":")]
        day_count = int(days)
    except ValueError:
        return None
    # with a day part SLURM counts from the left (D-HH, D-HH:MM, D-HH:MM:SS); without one, from
    # the right (MM, MM:SS, HH:MM:SS). Reading "1-12" as 12 minutes was wrong by a day and a half.
    if dash:
        hours, minutes, seconds = (parts + [0, 0, 0])[:3]
    elif len(parts) == 3:
        hours, minutes, seconds = parts
    elif len(parts) == 2:
        hours, minutes, seconds = 0, parts[0], parts[1]
    elif len(parts) == 1:
        hours, minutes, seconds = 0, parts[0], 0
    else:
        return None
    return day_count * 86400 + hours * 3600 + minutes * 60 + seconds


def _gpus_in_gres(gres: str):
    """GPUs per node from a GRES string: gpu:a100:4, gpu:4, gpu:4(S:0-1), (null).

    0 when the string names no GPU, and None when it names one WITHOUT a count (`gpu:a100`):
    counting that as 1 would invent a number and attribute it to sinfo, and everything derived
    from it (the share of cores and memory, and the refusal of a larger --gpus) would be wrong.
    """
    total, unknown = 0, False
    for field in (gres or "").split(","):
        field = field.split("(")[0]
        parts = field.strip().split(":")
        if not parts or parts[0].lower() != "gpu":
            continue
        count = parts[-1]
        if count.isdigit():
            total += int(count)
        else:
            unknown = True
    if unknown and not total:
        return None
    return total


def _sinfo_row(fields: list, partitions: list) -> bool:
    """Append one sinfo row, or say it is not one. Every column is checked for its own shape."""
    name, limit, gres, nodes, cpus, memory = fields
    if not re.fullmatch(r"[\w.-]+\*?", name or ""):
        return False
    if _time_to_seconds(limit) is None and limit.lower() not in ("infinite", "unlimited"):
        return False
    # sinfo appends `+` when the nodes of the partition are not all the same: the number is then
    # the smallest of them, which is the honest one to size a request with, and it is reported
    # as a minimum rather than passed off as the node's count
    if not cpus.rstrip("+").isdigit() or not memory.rstrip("+").isdigit():
        return False
    partitions.append({
        "name": name.rstrip("*"), "is_default": name.endswith("*"),
        "time_limit": limit, "time_limit_seconds": _time_to_seconds(limit),
        "gres": gres, "gpus_per_node": _gpus_in_gres(gres), "nodes": nodes,
        "cpus_per_node": int(cpus.rstrip("+")), "cpus_is_minimum": cpus.endswith("+"),
        "memory_mb_per_node": int(memory.rstrip("+")),
        "memory_is_minimum": memory.endswith("+"),
    })
    return True


def parse_site_report(text: str) -> dict:
    """Read the pasted command output. Strict: a line that does not match a known shape is kept
    in `unparsed` and never turned into a value.

    This is text from another machine, pasted by hand. It is data: nothing here executes it, and
    a field that cannot be read becomes a CONFIRM slot rather than a guess.
    """
    partitions, accounts, modules, binaries, unparsed = [], [], [], [], []
    max_array, gres_types = None, None
    for raw in (text or "").splitlines():
        line = raw.strip()
        if not line or line == "---":
            continue
        match = re.search(r"MaxArraySize\s*=\s*(\d+)", line, re.IGNORECASE)
        if match:
            max_array = int(match.group(1))
            continue
        match = re.search(r"GresTypes\s*=\s*(\S+)", line, re.IGNORECASE)
        if match:
            gres_types = match.group(1)
            continue
        if "|" in line:
            # a delimited row: which command it came from is decided by the column count the
            # command was asked for, not by what the text looks like
            fields = [f.strip() for f in line.split("|")]
            if len(fields) == 6 and _sinfo_row(fields, partitions):
                continue
            if len(fields) == 3 and re.fullmatch(r"[A-Za-z][\w.-]*", fields[0]):
                accounts.append({"account": fields[0], "partition": fields[1], "qos": fields[2]})
                continue
            unparsed.append(line)
            continue
        if re.fullmatch(r"/[\w./+@:-]*namd[\w./+@:-]*", line, re.IGNORECASE):
            binaries.append(line)
            continue
        # a module line. `module avail` pads several names across one line and the names can be
        # hierarchical (`apps/namd/3.0.3`), with markers like (D) or (L) after them. Taking the
        # whole token and keeping only the ones that name NAMD is what matches the real output;
        # a pattern anchored on `namd/` missed `apps/namd/3.0.3` entirely.
        found = [token for token in re.split(r"\s+", line)
                 if "namd" in token.lower() and not token.startswith("/")
                 and re.fullmatch(r"[\w.+-]+(?:/[\w.+-]+)*", token)]
        if found:
            modules.extend(found)
            continue
        # anything else is noise, a header or a message: it is reported, never turned into a
        # value. A short line of free text used to become an account.
        unparsed.append(line)
    return {"partitions": partitions, "accounts": accounts,
            "namd_modules": sorted(set(modules)), "namd_binaries": sorted(set(binaries)),
            "max_array_size": max_array, "gres_types": gres_types, "unparsed": unparsed}


def gpu_types(gres: str) -> list:
    """The GPU model names a GRES string lists: `gpu:a100:12(S:1,3)` -> ["a100"]."""
    out = []
    for field in (gres or "").split(","):
        parts = field.split("(")[0].strip().split(":")
        if len(parts) >= 3 and parts[0].lower() == "gpu":
            out.append(parts[1])
    return out


def group_partitions(partitions: list) -> list:
    """One entry per partition NAME, holding every row sinfo printed for it.

    A real cluster prints one row per kind of node: on the paste that exposed this, `gpu` had
    six rows (rtx_2080, rtx_3090, a100, V100...). Treating each row as its own partition made
    the code pick whichever came first — the weakest card on the list — and count the same
    partition several times in "N of total".
    """
    grouped = {}
    for row in partitions:
        entry = grouped.setdefault(row["name"], {"name": row["name"], "variants": []})
        entry["variants"].append(row)
        entry["is_default"] = entry.get("is_default") or row["is_default"]
    for entry in grouped.values():
        seconds = [v["time_limit_seconds"] for v in entry["variants"]]
        entry["time_limit"] = entry["variants"][0]["time_limit"]
        entry["time_limit_seconds"] = None if None in seconds else max(seconds)
        entry["gpu_types"] = sorted({t for v in entry["variants"] for t in gpu_types(v["gres"])})
        entry["has_gpu"] = any(v["gpus_per_node"] or v["gpus_per_node"] is None
                               for v in entry["variants"])
    return list(grouped.values())


def _node_count(entry: dict) -> int:
    total = 0
    for v in entry.get("variants") or []:
        digits = (v.get("nodes") or "").rstrip("+")
        total += int(digits) if digits.isdigit() else 0
    return total


def _pick_partition(partitions: list):
    """The proposal: a GPU partition, the default one if it has GPUs, else the largest.

    The rule is stated on the line it writes, because it IS a choice this code makes and the
    user may well want another one. Ordering by time limit alone was worse than it looked: on a
    real cluster it picked a 5-node partition with a 15-day limit over the general GPU one with
    18 nodes and 14 days, which is very likely a group's own queue. Size is a weak signal too,
    but it is a stated one, and every other GPU partition is listed beside it.
    """
    grouped = group_partitions(partitions)
    with_gpus = [p for p in grouped if p["has_gpu"]]
    pool = with_gpus or grouped
    if not pool:
        return None, "no partition could be read from the pasted output", len(grouped), []
    def key(p):
        return (p["is_default"], _node_count(p), p["time_limit_seconds"] or 0)
    best = sorted(pool, key=key, reverse=True)[0]
    if not with_gpus:
        why = "no partition reports a GPU at all, so this is simply the largest; check it"
    elif best["is_default"]:
        why = "the cluster's default partition, and it has GPUs"
    else:
        why = (f"the GPU partition with the most nodes ({_node_count(best)}), which is usually the "
               "general one rather than a group's own queue")
    rest = [f"{p['name']} ({p['time_limit']}, {_node_count(p)} nodes"
            + (f", {', '.join(p['gpu_types'])}" if p["gpu_types"] else ", no GPU") + ")"
            for p in grouped if p["name"] != best["name"]]
    return best, why, len(grouped), rest


def _pick_variant(entry: dict, gpu_type: str):
    """Which kind of node of that partition the request is sized for, and what is left open.

    With several kinds and no `gpu_type`, nothing is sized: taking the share from whichever row
    came first is how `gpu` on a real cluster ended up proposing its weakest card.
    """
    variants = entry.get("variants") or []
    if gpu_type.strip():
        wanted = gpu_type.strip().lower()
        matching = [v for v in variants if wanted in [t.lower() for t in gpu_types(v["gres"])]]
        if not matching:
            raise MdRunError(
                f"No node type of partition {entry['name']} reports a GPU called "
                f"{gpu_type.strip()!r}.", gpu_types_reported=entry.get("gpu_types"))
        if len(matching) > 1:
            return None, (f"Partition {entry['name']} reports {len(matching)} kinds of node with "
                          f"a {gpu_type.strip()} GPU and they differ in cores or memory, so "
                          "neither was sized: set cpus_per_task yourself.")
        return matching[0], None
    if len(variants) == 1:
        return variants[0], None
    if not variants:
        return None, None
    types = entry.get("gpu_types") or []
    return None, (f"Partition {entry['name']} is made of {len(variants)} kinds of node"
                  + (f" with GPUs of different models ({', '.join(types)})" if len(types) > 1 else "")
                  + ". Nothing was sized from them, and the request does not name a model, so the "
                    "job lands on whichever node the scheduler gives it. Name one with gpu_type "
                    "to size the request and pin the model.")


def resolve_site(parsed: dict, systems: int, *, partition: str = "", account: str = "",
                 time_limit: str = "", cpus_per_task: int = 0, gpus: int = 1, gpu_type: str = "",
                 namd_module: str = "", namd_binary: str = "") -> dict:
    """Turn the pasted output plus the user's overrides into the header of the script.

    Every field comes back as {value, mark, source}: REPORTED (read from a named command),
    CHOSEN (the user said so), CONFIRM (unknown, blocks the script) or EDITABLE (a rule, stated).
    """
    out, notes = {}, []
    sinfo = SINFO
    grouped = group_partitions(parsed["partitions"])
    if grouped:
        chosen, why, total, others = _pick_partition(parsed["partitions"])
    else:
        chosen, why, total, others = None, "no partition could be read", 0, []
    named = next((p for p in grouped if p["name"] == partition.strip()),
                 None) if partition.strip() else None
    if partition.strip() and not named:
        if grouped:
            raise MdRunError(
                f"The pasted output does not list a partition called {partition.strip()!r}.",
                partitions=[p["name"] for p in grouped])
        named = {"name": partition.strip(), "time_limit": None, "time_limit_seconds": None,
                 "variants": [], "gpu_types": [], "has_gpu": False, "is_default": False}
        notes.append(f"Partition {named['name']} is used because you named it; nothing was "
                     "parsed from the pasted output to check it against.")
    picked = named or chosen
    if named:
        others = [p["name"] for p in grouped if p["name"] != picked["name"]]
    variant = None
    if not picked:
        out["partition"] = {"value": CONFIRM_PARTITION, "mark": "CONFIRM",
                            "source": "no partition could be read from the pasted output"}
    else:
        # a partition can be several kinds of node, and which one the job lands on is what
        # decides the GPU model, the cores and the memory: it is named or left open
        variant, variant_note = _pick_variant(picked, gpu_type)
        out["partition"] = {
            "value": picked["name"],
            "mark": "CHOSEN" if named else "REPORTED",
            "source": (("you named it" if named else f"{sinfo}: {why}")
                       + f" (1 of {total} partitions read"
                       + (f"; the others: {', '.join(others)})" if others else ")")),
            "node_types_reported": [
                {"gres": v["gres"], "cpus_per_node": v["cpus_per_node"],
                 "cpus_is_minimum": v.get("cpus_is_minimum", False),
                 "memory_mb_per_node": v["memory_mb_per_node"],
                 "gpus_per_node": v["gpus_per_node"]} for v in picked["variants"]],
            "gpu_types_reported": picked["gpu_types"],
        }
        if variant_note:
            notes.append(variant_note)

    # --time: the partition reports a LIMIT. How long the run takes is not known by anyone yet,
    # so the value is a ceiling and the way to fix it is to measure the first chunk.
    limit = (picked or {}).get("time_limit")
    if time_limit.strip():
        asked = _time_to_seconds(time_limit.strip())
        reported = (picked or {}).get("time_limit_seconds")
        if asked is None:
            raise MdRunError(f"{time_limit.strip()!r} is not a SLURM time (D-HH:MM:SS or HH:MM:SS).")
        if reported is not None and asked > reported:
            raise MdRunError(f"{time_limit.strip()} is over the limit the cluster reports for "
                             f"partition {picked['name']} ({limit}).")
        out["time"] = {"value": time_limit.strip(), "mark": "CHOSEN", "source": "you asked for it",
                       "ceiling_reported": limit}
    elif limit and _time_to_seconds(limit) is not None:
        out["time"] = {"value": limit, "mark": "REPORTED",
                       "source": f"{sinfo}: the limit of partition {picked['name']}, used as a "
                                 "CEILING, not an estimate"}
    else:
        out["time"] = {"value": FALLBACK_TIME, "mark": "EDITABLE",
                       "source": "no usable limit was reported (unlimited or unreadable), so a "
                                 "chunked run is asked for 24 h; it restarts chunk by chunk"}
    notes.append("--time is a ceiling, never an estimate: no wall time per nanosecond has been "
                 f"measured for this system anywhere. Measure the first chunk with `{HOW_TO_MEASURE_TIME}` "
                 "and lower it.")

    # GPUs: --gres=gpu:N is the widest supported form and is what the user's own template used;
    # --gpus=N goes in commented, because some sites want that instead
    reported_gpus = (variant or {}).get("gpus_per_node")
    if gpus < 1:
        raise MdRunError("gpus must be at least 1.")
    if reported_gpus and gpus > reported_gpus:
        raise MdRunError(f"partition {picked['name']} reports {reported_gpus} GPU(s) per node and "
                         f"{gpus} were asked for.")
    if reported_gpus:
        gpu_source = f"{sinfo}: the node type reports {reported_gpus} GPU(s) per node"
    elif variant is None and picked and (picked.get("variants") or []):
        gpu_source = (f"{sinfo}: {len(picked['variants'])} kinds of node in this partition and "
                      "none was settled on, so how many GPUs one has was not taken from the "
                      "output")
    elif reported_gpus is None:
        gpu_source = (f"{sinfo}: the GRES of this node type names a GPU but not how many "
                      f"({(variant or {}).get('gres')}), so the count is unknown")
    else:
        gpu_source = "no GPU was reported for this partition; check that it has one"
    known = parsed["gres_types"]
    if known and "gpu" not in known.lower():
        gpu_source += (f". The site's GresTypes is {known!r}, which does not include `gpu`: "
                       "--gres=gpu may be rejected here")
    # sinfo prints the model one way and the site may have configured it under a longer name
    # (`a100` against `a100-sxm4-40gb`): both strings are here, so the mismatch is reported
    # rather than left for sbatch to reject
    if gpu_type.strip() and known:
        configured = [k.strip() for k in known.split(",") if k.strip()]
        wanted = gpu_type.strip().lower()
        if not any(k.lower() == wanted for k in configured):
            close = [k for k in configured if wanted in k.lower() or k.lower() in wanted]
            gpu_source += (f". NOTE: sinfo calls this model {gpu_type.strip()!r} but the site's "
                           f"GresTypes does not list that exact name"
                           + (f"; the closest it lists is {', '.join(close)}, so if sbatch "
                              "rejects the request, try that spelling" if close else
                              ", so sbatch may reject the request"))
    requested_type = ""
    if gpu_type.strip():
        requested_type = next((t for t in (picked or {}).get("gpu_types") or []
                               if t.lower() == gpu_type.strip().lower()), gpu_type.strip())
        gpu_source += f". Model pinned to {requested_type} because you named it"
    elif len(((picked or {}).get("gpu_types") or [])) > 1:
        gpu_source += (". The request does not name a model, so the job lands on whichever of "
                       f"{', '.join(picked['gpu_types'])} the scheduler gives it")
    out["gpus"] = {"value": gpus, "mark": "CHOSEN" if gpus != 1 or requested_type else "EDITABLE",
                   "source": gpu_source, "gres_known_to_the_site": known,
                   "gres": f"gpu:{requested_type}:{gpus}" if requested_type else f"gpu:{gpus}"}

    # cores and memory: the share of a node that goes with the GPUs asked for
    # the share of a node that goes with ONE GPU; `share_over` is None when the GRES did not
    # say how many a node has, and then nothing is divided and nothing is proposed
    share_over = reported_gpus if reported_gpus else (None if reported_gpus is None else 1)
    per_node_cpus = (variant or {}).get("cpus_per_node")
    if cpus_per_task:
        if per_node_cpus and cpus_per_task > per_node_cpus:
            raise MdRunError(f"partition {picked['name']} reports {per_node_cpus} cores per node "
                             f"and {cpus_per_task} were asked for.")
        out["cpus_per_task"] = {"value": cpus_per_task, "mark": "CHOSEN",
                                "source": "you asked for it"}
    elif per_node_cpus and share_over:
        least = " at least" if (variant or {}).get("cpus_is_minimum") else ""
        out["cpus_per_task"] = {"value": max(1, per_node_cpus // share_over), "mark": "EDITABLE",
                                "source": f"{sinfo}:{least} {per_node_cpus} cores per node divided "
                                          f"by {share_over} GPU(s) per node (the share of ONE GPU; "
                                          "raise it yourself if you asked for more)"}
    else:
        out["cpus_per_task"] = {"value": None, "mark": "EDITABLE",
                                "source": ("the GRES does not say how many GPUs a node has, so no "
                                           "share can be worked out" if per_node_cpus else
                                           "no single kind of node was settled on, so no core "
                                           "count was taken from the output")
                                          + "; the line is commented out and NAMD then runs on "
                                            "one thread unless you set it"}
    per_node_mb = (variant or {}).get("memory_mb_per_node")
    if per_node_mb and share_over:
        out["mem"] = {"value": f"{max(1, (per_node_mb // share_over) // 1024)}G", "mark": "EDITABLE",
                      "source": f"{sinfo}: {per_node_mb} MB per node divided by {share_over} "
                                "GPU(s) per node (the share of ONE GPU)"}
    else:
        if per_node_mb:
            why_mem = "the GRES does not say how many GPUs a node has, so no share can be worked out"
        elif variant is None and picked and (picked.get("variants") or []):
            why_mem = "no single kind of node was settled on, so no memory was taken from the output"
        else:
            why_mem = "no memory was reported"
        out["mem"] = {"value": None, "mark": "EDITABLE",
                      "source": why_mem + "; the line is commented out and the site's default "
                                          "applies"}

    # account: only if the site reports one
    if account.strip():
        out["account"] = {"value": account.strip(), "mark": "CHOSEN", "source": "you named it"}
    else:
        names = sorted({a["account"] for a in parsed["accounts"]})
        if len(names) == 1:
            out["account"] = {"value": names[0], "mark": "REPORTED",
                              "source": f"{SACCTMGR}: the only account reported for you"}
        elif names:
            out["account"] = {"value": CONFIRM_ACCOUNT, "mark": "CONFIRM",
                              "source": f"{SACCTMGR} reports {len(names)} accounts; name the one to "
                                        f"charge: {', '.join(names)}"}
        else:
            out["account"] = {"value": None, "mark": "EDITABLE",
                              "source": f"{SACCTMGR} reported no account, so the line is left "
                                        "commented out; uncomment it if sbatch asks for one"}

    # NAMD: a module, or a binary the user compiled; never both invented
    if namd_binary.strip():
        out["namd_binary"] = {"value": namd_binary.strip(), "mark": "CHOSEN",
                              "source": "you gave the path"}
    elif len(parsed["namd_binaries"]) == 1:
        out["namd_binary"] = {"value": parsed["namd_binaries"][0], "mark": "REPORTED",
                              "source": "which: the only NAMD on your PATH (1 of 1 read)"}
    else:
        if parsed["namd_binaries"]:
            why_bin = (f"which found {len(parsed['namd_binaries'])} NAMD binaries: "
                       f"{', '.join(parsed['namd_binaries'])}")
        elif parsed["namd_modules"]:
            why_bin = ("no NAMD was on your PATH, and the module that provides it was not loaded "
                       "when you ran `which`. Load it and run `which namd3 namd2` again, or give "
                       "the path: the name of the executable a module installs is not something "
                       "this tool can know")
        else:
            why_bin = "no NAMD was found on your PATH"
        out["namd_binary"] = {"value": CONFIRM_BINARY, "mark": "CONFIRM", "source": why_bin}
    if namd_module.strip().lower() in ("none", "no", "-"):
        out["namd_module"] = {"value": None, "mark": "CHOSEN",
                              "source": "you said no module is to be loaded"}
    elif namd_module.strip():
        out["namd_module"] = {"value": namd_module.strip(), "mark": "CHOSEN",
                              "source": "you named it"}
    elif namd_binary.strip() and len(parsed["namd_modules"]) == 1:
        # giving a path used to drop the module silently, and that is how a job dies in the queue:
        # a binary installed by a module usually needs it for its CUDA libraries, so the module
        # stays loaded and the reason says why
        out["namd_module"] = {"value": parsed["namd_modules"][0], "mark": "REPORTED",
                              "source": "module avail: the only NAMD module. It is loaded even "
                                        "though you gave a binary path, because a binary a module "
                                        "installed normally needs it on LD_LIBRARY_PATH; pass "
                                        "namd_module='none' if it must not be loaded"}
    elif namd_binary.strip():
        out["namd_module"] = {"value": None, "mark": "CHOSEN",
                              "source": "you gave a binary path and no NAMD module was listed, so "
                                        "none is loaded"}
    elif len(parsed["namd_modules"]) == 1:
        out["namd_module"] = {"value": parsed["namd_modules"][0], "mark": "REPORTED",
                              "source": "module avail: the only NAMD module (1 of 1 read)"}
    elif parsed["namd_modules"]:
        out["namd_module"] = {"value": CONFIRM_MODULE, "mark": "CONFIRM",
                              "source": f"module avail lists {len(parsed['namd_modules'])}: "
                                        f"{', '.join(parsed['namd_modules'])}"}
    elif out["namd_binary"]["mark"] in ("REPORTED", "CHOSEN"):
        out["namd_module"] = {"value": None, "mark": "REPORTED",
                              "source": "no NAMD module was listed, and a binary was found on "
                                        "your PATH, so no module is loaded"}
    else:
        out["namd_module"] = {"value": CONFIRM_MODULE, "mark": "CONFIRM",
                              "source": "module avail listed no NAMD module"}

    # the name of the build often says which flags are right, and that is evidence from the
    # pasted output rather than an assumption: a multicore build is launched directly, while an
    # mpi or ofi one needs charmrun or srun and these flags would be wrong
    named_build = f"{out['namd_binary']['value'] or ''} {out['namd_module']['value'] or ''}".lower()
    opts_source = NAMD_OPTS_SOURCE
    if "multicore" in named_build:
        opts_source += (". The build you gave is named `multicore"
                        + ("-CUDA`" if "cuda" in named_build else "`")
                        + ", which is the kind these flags are for")
    elif any(kind in named_build for kind in ("mpi", "ofi", "verbs", "ucx")):
        opts_source += (". WARNING: the build you gave looks like an MPI or network one, which is "
                        "normally started with charmrun or srun instead; these flags are for a "
                        "multicore build")
    out["namd_opts"] = {"value": DEFAULT_NAMD_OPTS, "mark": "EDITABLE", "source": opts_source}
    # the script runs on the cluster with the user's own credentials, so a value that reaches it
    # must be a plain path or module name. Nothing here executes the pasted text, but writing it
    # into a bash assignment unchecked would hand the shell whatever it contained.
    for key in ("namd_binary", "namd_module"):
        value = out[key]["value"]
        if value and not re.fullmatch(r"[\w./+@:-]+", value):
            raise MdRunError(f"{value!r} is not a plain path or module name; it would be written "
                             "into the script and run on the cluster.")

    # the array has to fit what the site allows
    if parsed["max_array_size"] is not None and systems > parsed["max_array_size"]:
        raise MdRunError(f"{systems} systems need an array of that size, and the site reports "
                         f"MaxArraySize={parsed['max_array_size']}. Split them into several "
                         "submissions.", max_array_size=parsed["max_array_size"])
    out["max_array_size"] = (
        {"value": parsed["max_array_size"], "mark": "REPORTED",
         "source": f"scontrol show config; the array of {systems} system(s) fits in it"}
        if parsed["max_array_size"] is not None else
        {"value": None, "mark": "EDITABLE",
         "source": "scontrol show config reported no MaxArraySize, so the size of the array was "
                   "NOT checked against the site's limit"})
    to_confirm = {k: v["source"] for k, v in out.items() if v["mark"] == "CONFIRM"}
    if parsed["unparsed"]:
        notes.append(f"{len(parsed['unparsed'])} pasted line(s) were not understood and were "
                     "ignored rather than guessed; they are in `unparsed`.")
    notes.append("What the cluster reports is not a promise that the job runs: that NAMD works "
                 "there, and that --gres=gpu is the form this site wants, only the first sbatch "
                 "will tell you.")
    return {"fields": out, "to_confirm": to_confirm, "notes": notes,
            "unparsed": parsed["unparsed"]}


# ------------------------------------------------------------------ archive

def extract_archive(archive: Path, dest: Path) -> Path:
    """Extract the whole archive under `dest` and return its single top folder.

    tarfile's "data" filter plus an explicit check that every member stays inside one top
    folder: this is a file the user downloaded from a web service, so it is untrusted input.
    """
    try:
        with tarfile.open(archive, "r:gz") as tar:
            members, roots, total = [], set(), 0
            for member in tar.getmembers():
                name = member.name.replace("\\", "/")
                if name.startswith("/") or ".." in name.split("/"):
                    raise MdRunError(f"The archive has an unsafe path ({member.name}).")
                roots.add(name.split("/")[0])
                if not (member.isfile() or member.isdir()):
                    continue
                if member.size > md_package._MAX_MEMBER_BYTES:
                    raise MdRunError(f"{member.name} is larger than the {md_package._MAX_MEMBER_BYTES} "
                                     "byte limit of one member.")
                total += member.size
                members.append(member)
            if len(roots) != 1:
                raise MdRunError(f"The archive has {len(roots)} top folders; one CHARMM-GUI job expected.")
            if len(members) > MAX_MEMBERS:
                raise MdRunError(f"The archive has {len(members)} members, over the {MAX_MEMBERS} limit.")
            if total > MAX_TOTAL_BYTES:
                raise MdRunError(f"The archive expands to {total} bytes, over the {MAX_TOTAL_BYTES} limit.")
            tar.extractall(dest, members=members, filter="data")
    except tarfile.TarError as exc:
        raise MdRunError(f"Could not read {archive.name} as a .tar.gz: {exc!r}")
    job = dest / roots.pop()
    if not job.is_dir():
        raise MdRunError("The archive's top entry is not a folder.")
    return job


def _job_id(folder):
    """The CHARMM-GUI job number inside the folder name, or None.

    The folder is `charmm-gui-9418739369`; the job is 9418739369. Calling the folder name a job
    id is wrong in a record that an analyst may quote, and guessing one when the name does not
    carry it would be worse.
    """
    match = re.fullmatch(r"charmm-gui-(\d+)", (folder or "").strip())
    return match.group(1) if match else None


def readme_scheme(readme: Path) -> dict:
    """The run order, read from the archive's namd/README (CHARMM-GUI writes it as csh).

    Refused rather than guessed when the file is not that shape: the names of the inputs and
    the chaining by `sed` are what the script reproduces, so reading them wrong would silently
    run the wrong thing.
    """
    if not readme.is_file():
        raise MdRunError("The archive has no namd/README, which is where the run order is written. "
                         "A Solution Builder job with NAMD inputs is expected.")
    text = readme.read_text(encoding="utf-8", errors="replace")
    out = {}
    for key in ("equi_prefix", "prod_prefix", "prod_step"):
        match = re.search(rf"^\s*set\s+{key}\s*=\s*(\S+)", text, re.MULTILINE)
        if not match:
            raise MdRunError(f"namd/README does not set {key}: this is not the CHARMM-GUI README "
                             "this tool knows how to translate.")
        out[key] = match.group(1)
    match = re.search(r"^\s*set\s+cntmax\s*=\s*(\d+)", text, re.MULTILINE)
    if not match:
        raise MdRunError("namd/README does not set cntmax, which is the number of production "
                         "chunks; this tool will not invent one.")
    out["chunks"] = int(match.group(1))
    out["source"] = "namd/README of the archive"
    return out


def _line_of(path: Path, pattern: str) -> dict:
    """{line, text, value} of the first line matching `pattern`, or None.

    The caller shows this to the analyst so an edit is made in a located line with a known
    current value, instead of being described from memory.
    """
    if not path.is_file():
        return None
    for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        match = re.match(pattern, line)
        if match:
            return {"line": number, "text": line.strip(), "value": match.group(1)}
    return None


def input_values(namd: Path, scheme: dict) -> dict:
    """What the NAMD inputs of the archive actually say. Every value carries its line.

    The real production input ends with BOTH `numsteps 500000;` and `run 500000;` (lines 147 and
    148 of job 9418739369). `run` is the command that executes the steps, so it is what the
    length of a chunk is read from, and `numsteps` alone would be the wrong line to tell an
    analyst to edit: raising it while `run` stays would keep every chunk at its old length
    while the report claimed otherwise. Both lines are reported, and a disagreement between
    them is carried out as a warning instead of being resolved here.
    """
    production = namd / f"{scheme['prod_prefix']}.inp"
    equilibration = namd / f"{scheme['equi_prefix']}.inp"
    for needed in (production, equilibration):
        if not needed.is_file():
            raise MdRunError(f"The archive has no namd/{needed.name}, named by its own README.")
    numsteps = _line_of(production, r"\s*numsteps\s+(\d+)")
    run = _line_of(production, r"\s*run\s+(\d+)")
    timestep = _line_of(production, r"\s*timestep\s+([\d.]+)")
    steps = run or numsteps
    if not steps or not timestep:
        raise MdRunError(f"namd/{production.name} has no `run` or `numsteps` line, or no "
                         "`timestep` line, so the length of a chunk cannot be read.")
    warning = None
    if run and numsteps and run["value"] != numsteps["value"]:
        warning = (f"namd/{production.name} says `run {run['value']}` on line {run['line']} and "
                   f"`numsteps {numsteps['value']}` on line {numsteps['line']}. NAMD runs the "
                   "steps of `run`, so that is what is reported here; the two should be made to "
                   "agree before this is simulated.")
    return {
        "production_input": production.name,
        "equilibration_input": equilibration.name,
        "steps_per_chunk": steps,
        "steps_per_chunk_from": "run" if run else "numsteps",
        "numsteps_line": numsteps,
        "run_line": run,
        "timestep_fs": timestep,
        "nanoseconds_per_chunk": round(int(steps["value"]) * float(timestep["value"]) / 1_000_000, 4),
        "length_warning": warning,
        "equilibration_minimize": _line_of(equilibration, r"\s*minimize\s+(\d+)"),
        "equilibration_run": _line_of(equilibration, r"\s*run\s+(\d+)"),
        "equilibration_outputname": _line_of(equilibration, r"\s*set\s+outputname\s+(\S+?);?\s*$"),
        "restartfreq": _line_of(production, r"\s*restartfreq\s+(\d+)"),
        "dcdfreq": _line_of(production, r"\s*dcdfreq\s+(\d+)"),
    }


def enable_gpu_resident(production: Path) -> dict:
    """Insert `CUDASOAintegrate on` into the production input, as an edit that is declared.

    It goes BEFORE the first `minimize` / `run` / `numsteps` line, never at the end of the file:
    NAMD executes the simulation when it reaches `run`, so a keyword written after it would
    either be ignored while this tool reported the mode as enabled, or abort the job.

    Only the production input. CHARMM-GUI's equilibration carries positional restraints, and
    whether GPU-resident NAMD supports them has NOT been verified here, so that file is left
    exactly as it came.
    """
    before = dock.sha256(production)
    lines = production.read_text(encoding="utf-8", errors="replace").splitlines()
    if any(re.match(r"\s*CUDASOAintegrate\s+on", line, re.IGNORECASE) for line in lines):
        return {"file": f"namd/{production.name}", "change": "already present",
                "sha256_before": before, "sha256_after": before, "line": None}
    where = next((i for i, line in enumerate(lines)
                  if re.match(r"\s*(minimize|run|numsteps)\s+\d", line)), None)
    if where is None:
        raise MdRunError(f"namd/{production.name} has no minimize/run/numsteps line, so there is "
                         "no safe place to put CUDASOAintegrate: it was not edited.")
    block = ["# Added by build_md_job, NOT written by CHARMM-GUI: GPU-resident integration",
             GPU_RESIDENT_LINE, ""]
    production.write_text("\n".join(lines[:where] + block + lines[where:]) + "\n", encoding="utf-8")
    # the line numbers reported are the ones IN THE FILE AS WRITTEN: the insertion pushed
    # everything below it down, and the caller re-reads the inputs for the same reason
    return {"file": f"namd/{production.name}",
            "change": f"inserted '{GPU_RESIDENT_LINE}' at line {where + 2} of the file as written, "
                      f"before '{lines[where].strip()}' (now line {where + len(block) + 1}), "
                      "which is the line that executes the run",
            "line": where + 2, "sha256_before": before, "sha256_after": dock.sha256(production)}


# ------------------------------------------------------------------ the script

def system_entry(rel: str, label: str, scheme: dict, values: dict) -> dict:
    """One system of the array: where it is and the names its own archive uses."""
    return {
        "dir": rel, "label": label,
        "equi": scheme["equi_prefix"], "prod": scheme["prod_prefix"], "step": scheme["prod_step"],
        # the file that says the equilibration finished is named by the input's own `set
        # outputname` when it has one, not by the README's prefix; they agree in the real
        # archive and nothing here assumes they must
        "done": (values["equilibration_outputname"] or {}).get("value") or scheme["equi_prefix"],
        "chunks": scheme["chunks"],
        "nanoseconds_per_chunk": values["nanoseconds_per_chunk"],
        "production_input": values["production_input"],
        "run_line": (values["run_line"] or {}).get("line"),
        "numsteps_line": (values["numsteps_line"] or {}).get("line"),
    }


def _column(systems: list, key: str, quote: bool = True) -> str:
    return "\n".join(f'  "{s[key]}"' if quote else f"  {s[key]}" for s in systems)


def _sbatch(flag, field, extra: str = "") -> str:
    """One #SBATCH line with its mark and where the value came from.

    A value of None comments the line out: a flag this code cannot fill honestly is left for the
    site's default rather than guessed.
    """
    value, mark, source = field["value"], field["mark"], field["source"]
    body = f"{flag}={value}"
    if value is None:
        return f"##SBATCH {flag}=<unset>   # {mark}: {source}{extra}"
    return f"#SBATCH {body}   # {mark}: {source}{extra}"


def slurm_text(job_name: str, systems: list, site: dict) -> str:
    """The batch script, for one system or for several as a SLURM array.

    One array task = one system, run independently, which is the shape of the user's own working
    template (2026-11-09). A single system is the same script with --array=0-0, so there is one
    generator and one body.

    Nothing about the cluster is known to this code: every #SBATCH line comes from `site`, which
    resolve_site built out of the output the user pasted, and carries its own mark:
      REPORTED  read from the pasted output of a command, which is named on the line
      CHOSEN    the user's decision
      CONFIRM   not known: the script refuses to run until it is filled in
      EDITABLE  a knob whose default comes from the archive or from a stated rule
    The per-system names (inputs, prefixes, chunk count) are read from each archive, so systems
    built by different CHARMM-GUI jobs can sit in the same array without sharing names.
    """
    last = len(systems) - 1
    f = site["fields"]
    gres = f["gpus"]
    header = "\n".join([
        f"#SBATCH --job-name={job_name}   # EDITABLE",
        _sbatch("--partition", f["partition"]),
        _sbatch("--time", f["time"]),
        f"#SBATCH --gres={gres['gres']}   # {gres['mark']}: {gres['source']}",
        f"##SBATCH --gpus={gres['value']}   # the same request in the newer form; some sites want "
        "this one instead of --gres, and only sbatch will say",
        "#SBATCH --nodes=1   # EDITABLE: one node per system",
        "#SBATCH --ntasks=1   # EDITABLE: a multicore NAMD build is one task",
        _sbatch("--cpus-per-task", f["cpus_per_task"]),
        _sbatch("--mem", f["mem"]),
        f"#SBATCH --array=0-{last}   # EDITABLE: one task per system, see SYSTEMS below",
        "#SBATCH --output=%x_%A_%a.out",
        "#SBATCH --error=%x_%A_%a.err",
        _sbatch("--account", f["account"]),
    ])
    listing = "\n".join(f"#   task {i}: {s['label']} -> {s['dir']} "
                        f"({s['chunks']} x {s['nanoseconds_per_chunk']} ns)"
                        for i, s in enumerate(systems))
    length_note = "\n".join(
        f"#   {s['label']}: BOTH `run` (line {s['run_line']}) and `numsteps` (line "
        f"{s['numsteps_line']}) of {s['dir']}/{s['production_input']}"
        for s in systems)
    return f"""#!/bin/bash
{header}
#
# Written by build_md_job (protein-research MCP). Nothing here was run or verified on any
# cluster: the lines above come from the output you pasted, and the ones marked CONFIRM stop
# this script until you fill them in. --time is a CEILING, not an estimate; measure the first
# chunk with `{HOW_TO_MEASURE_TIME}` and lower it.
#
# Systems in this array:
{listing}
#
# Submit from the folder that holds this script and those system folders:  sbatch <this file>
set -euo pipefail

# ---- NAMD ------------------------------------------------------------------------------
# {f["namd_module"]["mark"]}: {f["namd_module"]["source"]}
NAMD_MODULE="{f["namd_module"]["value"] or ""}"
# {f["namd_binary"]["mark"]}: {f["namd_binary"]["source"]}
NAMD_BIN="{f["namd_binary"]["value"]}"
# EDITABLE: {f["namd_opts"]["source"]}
NAMD_OPTS="{f["namd_opts"]["value"]}"

# ---- EDITABLE: the systems, one per array task -----------------------------------------
# Paths are relative to the submission directory, so the tree can be copied anywhere. Only
# systems whose archive passed check_charmm_gui_output are listed here; a folder added by hand
# is one that nothing verified.
SYSTEMS=(
{_column(systems, "dir")}
)
LABELS=(
{_column(systems, "label")}
)
EQUI=(
{_column(systems, "equi")}
)
PROD=(
{_column(systems, "prod")}
)
STEP=(
{_column(systems, "step")}
)
EQUI_DONE=(
{_column(systems, "done")}
)
# EDITABLE: how many nanoseconds. Each chunk restarts from the previous one. These counts are
# `cntmax` from each archive's own namd/README; they were not chosen by this tool. To make the
# chunks themselves longer, change in the input of that system:
{length_note}
# NAMD executes the steps of `run`, so changing only `numsteps` leaves the chunk as it is. Those
# are STEPS, not seconds. Editing an input means that system is no longer the one CHARMM-GUI
# returned and check_charmm_gui_output verified.
CHUNKS=(
{_column(systems, "chunks", quote=False)}
)

# ---- refuse to run while a hole is left ------------------------------------------------
# The commented ##SBATCH --account line is excluded on purpose: it is inert, and if the site
# needs no account it is meant to stay as it is.
placeholder="CONFIRM""_"
if grep -n "$placeholder" "$0" | grep -qv '^[0-9]*:##SBATCH'; then
    echo "This script still has CONFIRM placeholders:" >&2
    grep -n "$placeholder" "$0" | grep -v '^[0-9]*:##SBATCH' >&2
    echo "Fill in the NAMD module and binary and submit again. The commented --account line" >&2
    echo "does not have to be touched." >&2
    exit 2
fi

# ---- one array task = one system -------------------------------------------------------
task="${{SLURM_ARRAY_TASK_ID:-0}}"
if [ "$task" -ge "${{#SYSTEMS[@]}}" ]; then
    echo "Array task $task has no system: SYSTEMS has ${{#SYSTEMS[@]}} entries." >&2
    exit 2
fi
# Slurm runs a COPY of this script from its spool directory, so $0 is not this folder
here="${{SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}}"
sysdir="$here/${{SYSTEMS[$task]}}"
label="${{LABELS[$task]}}"
equi="${{EQUI[$task]}}"; prod="${{PROD[$task]}}"; step="${{STEP[$task]}}"
equi_done="${{EQUI_DONE[$task]}}"; chunks="${{CHUNKS[$task]}}"

echo "=== task $task : $label : $sysdir ==="
if [ ! -d "$sysdir" ]; then
    echo "ERROR: $sysdir is not there. Submit from the folder that holds the system folders." >&2
    exit 2
fi
cd "$sysdir"
if [ -n "$NAMD_MODULE" ]; then
    module load "$NAMD_MODULE"
fi
echo "[build_md_job] $(date -Is) host $(hostname) job ${{SLURM_JOB_ID:-none}} task $task"
echo "[build_md_job] NAMD: $NAMD_BIN $NAMD_OPTS"

# Seconds the last PRODUCTION chunk took, kept in a file so that a re-submitted job starts with
# a measurement instead of running its first chunk unguarded. The equilibration never feeds it:
# it is a different amount of work, so it would under-estimate a chunk.
timing_file=".build_md_job_chunk_seconds"
last_seconds=0
stopped_early=0
[ -f "$timing_file" ] && last_seconds=$(cat "$timing_file")

run_namd () {{   # $1 = input, $2 = output log; $3 = "chunk" when the timing is a chunk's
    echo "[build_md_job] $(date -Is) start $1"
    local t0=$SECONDS
    # shellcheck disable=SC2086
    if ! "$NAMD_BIN" $NAMD_OPTS "$1" > "$2"; then
        echo "ERROR: $1 failed for $label ($sysdir) - stopping this task only." >&2
        exit 1
    fi
    local took=$(( SECONDS - t0 ))
    echo "[build_md_job] $(date -Is) finished $1 in ${{took}} s -> $2"
    if [ "${{3:-}}" = "chunk" ]; then
        last_seconds=$took
        echo "$took" > "$timing_file"
    fi
}}

fits () {{        # $1 = seconds allowed for the next chunk; $2 = what the last one measured
    if [ -z "${{SLURM_JOB_END_TIME:-}}" ]; then
        echo "[build_md_job] SLURM_JOB_END_TIME is not set: no wall-time guard, the job may be cut mid-chunk"
        return 0
    fi
    local left=$(( SLURM_JOB_END_TIME - $(date +%s) ))
    if [ "$left" -lt "$1" ]; then
        echo "[build_md_job] ${{left}} s of wall time left; the last chunk took $2 s and this one is"
        echo "[build_md_job] allowed $1 s (that 20 % margin is chosen, not measured): stopping here."
        echo "[build_md_job] Nothing is lost: submit this script again and it continues at the next chunk."
        return 1
    fi
    return 0
}}

# ---- 1. equilibration, as CHARMM-GUI wrote it ------------------------------------------
if [ -f "${{equi_done}}.coor" ]; then
    echo "[build_md_job] ${{equi_done}}.coor is already there: equilibration skipped"
else
    run_namd "${{equi}}.inp" "${{equi}}.out"
fi

# ---- 2. production, chunk by chunk -----------------------------------------------------
# The scheme of the archive's README, in bash: chunk 1 restarts from the equilibration and
# every later chunk from the chunk before it, by renaming inside a copy of the input.
for (( cnt=1; cnt<=chunks; cnt++ )); do
    out="${{step}}_${{cnt}}"
    if [ -f "${{out}}.coor" ]; then
        echo "[build_md_job] chunk ${{cnt}} (${{out}}) is already done: skipped"
        continue
    fi
    # 20 % over the last measured chunk: a margin chosen here, not a measurement of anything
    if [ "$last_seconds" -gt 0 ] && ! fits $(( last_seconds * 12 / 10 )) "$last_seconds"; then
        stopped_early=1
        break
    fi
    if [ "$cnt" -eq 1 ]; then
        sed "s/${{prod}}/${{out}}/" "${{prod}}.inp" > "${{step}}_run.inp"
    else
        prev="${{step}}_$(( cnt - 1 ))"
        sed "s/${{equi}}/${{prev}}/" "${{prod}}.inp" | sed "s/${{prod}}/${{out}}/" > "${{step}}_run.inp"
    fi
    run_namd "${{step}}_run.inp" "${{out}}.out" chunk
done

done_chunks=0
for (( cnt=1; cnt<=chunks; cnt++ )); do
    [ -f "${{step}}_${{cnt}}.coor" ] && done_chunks=$(( done_chunks + 1 ))
done
echo "=== task $task : $label : ${{done_chunks}} of ${{chunks}} chunks done ==="
ls -l "${{step}}"_*.dcd 2>/dev/null || echo "  (no .dcd written)"
if [ "$stopped_early" -eq 1 ] || [ "$done_chunks" -lt "$chunks" ]; then
    echo "[build_md_job] NOT finished: submit this script again to continue." >&2
    exit 3
fi
echo "=== task $task : $label complete ==="
"""


def readme_text(package: dict, scheme: dict, values: dict, verified: dict, job_name: str,
                run_dir: Path, edits: list, site: dict) -> str:
    modelled = package.get("to_model") or []
    runs = ", ".join(f"{r['chain']} {r['first']}-{r['last']}" for r in modelled) or "none"
    holes = ("\n".join(f"- `{v}` - {why}" for v, why in site["to_confirm"].items())
             if site["to_confirm"] else
             "- none: every field could be filled in from what you pasted or from what you chose.")
    site_table = "\n".join(
        ["| field | value | mark | where it comes from |", "|---|---|---|---|"]
        + [f"| `{name}` | `{field['value']}` | {field['mark']} | {field['source']} |"
           for name, field in site["fields"].items()])
    commands = "\n".join(f"- `{c['command']}` - {c['gives']}" for c in SITE_COMMANDS)
    edited = ("\n".join(f"- `{e['file']}`: {e['change']}" for e in edits)
              if edits else "- none: the archive's inputs are exactly as CHARMM-GUI returned them.")
    total = round(scheme["chunks"] * values["nanoseconds_per_chunk"], 4)
    not_checked = "\n".join(f"- {n}" for n in verified.get("not_checked") or []) or "- (none)"
    return f"""# Molecular dynamics job: {package['pdb_id']} with {package['ligand'].get('comp_id') or 'LIG'}

Written by build_md_job ({utc_now()}). Nothing has been simulated and nothing was sent anywhere:
this folder is what you copy to the cluster, and the header was built from the output you pasted.

CHARMM-GUI job `{_job_id(verified.get('checks', {}).get('job_id')) or '(not in the folder name)'}`
(folder `{verified.get('checks', {}).get('job_id')}`), checked here before this script was written
({verified.get('verdict')}).

## Files
- `system/` - the CHARMM-GUI archive, extracted. The NAMD inputs live in `system/namd/`.
- `run_md.slurm` - the batch script. It runs `{values['equilibration_input']}` once and then
  `CHUNKS` production chunks, which is the scheme of `system/namd/README` translated to bash.
- `job.json` - every value used, where it was read, and what was left open.

## Where each line of the header comes from
Nothing about your cluster is known to this tool. Every `#SBATCH` line carries its own mark:

| mark | meaning |
|---|---|
| `REPORTED` | read from the output you pasted. The line names the command it came from. |
| `CHOSEN` | you said so. |
| `CONFIRM` | **not known.** The script refuses to run until it is filled in. |
| `EDITABLE` | a knob. Its default comes from the archive, or from a rule stated on the line. |

{site_table}

A line written as `##SBATCH` is commented out on purpose: that flag could not be filled in
honestly, so the site's own default applies. What your cluster reports is not a promise that the
job runs: whether NAMD works there, and whether this site wants `--gres=gpu` or `--gpus`, only
the first `sbatch` will tell you.

## Running several systems at once
One array task = one system. This script lists its own, with `--array=0-0`. `../../run_all.slurm`
(in `prep/md/`) is rebuilt on every `build_md_job` call with **every** system whose archive passed
the checks, one task each, and with this same header: its site values are not worked out again. Adding a folder to `SYSTEMS` by hand adds one that nothing verified.

## Before you submit
`run_md.slurm` aborts with exit code 2 while any `CONFIRM` placeholder is left in it, except a
commented `##SBATCH` line, which is inert.

{holes}

If the job log says `module: command not found`, the batch shell has no module system loaded and
the site's init line has to be sourced before `module load`. That is site configuration, not
something this tool can know.

To fill in or re-check any of the header, these are the commands the values came from:

{commands}

## Transfer and submit
This folder is at `{run_dir}`. Copy it whole, with whatever your machine has; neither line below
was tried from here:

```
scp -r "{run_dir}" <user>@<cluster>:<path>/{job_name}
rsync -av "{run_dir.as_posix()}/" <user>@<cluster>:<path>/{job_name}/   # from WSL: /mnt/c/...
```
Then:
```
ssh <user>@<cluster>
cd <path>/{job_name} && sbatch run_md.slurm
```
Re-submitting the same script continues where it stopped: a chunk whose `.coor` exists is skipped.

## How long it runs
{scheme['chunks']} chunks x {values['nanoseconds_per_chunk']} ns = **{total} ns**. The {scheme['chunks']} is `cntmax` in
`system/namd/README`, written by CHARMM-GUI; it was not chosen here.

- **more nanoseconds**: raise `CHUNKS` in `run_md.slurm`;
- **longer chunks**: in `system/namd/{values['production_input']}`, change BOTH of these lines
  together, because NAMD executes the steps of `run` and leaving it behind would keep the chunk
  at its current length:
  - line {values['run_line']['line'] if values['run_line'] else '-'}: `{values['run_line']['text'] if values['run_line'] else 'no run line'}`
  - line {values['numsteps_line']['line'] if values['numsteps_line'] else '-'}: `{values['numsteps_line']['text'] if values['numsteps_line'] else 'no numsteps line'}`

  {values['steps_per_chunk']['value']} steps x {values['timestep_fs']['value']} fs = {values['nanoseconds_per_chunk']} ns per chunk. These are STEPS, not seconds.

Editing that input means the system is no longer the one CHARMM-GUI returned and
check_charmm_gui_output verified; say so in any report.

`--time` is `{site["fields"]["time"]["value"]}` ({site["fields"]["time"]["source"]}). It is a
ceiling, not an estimate: how long THIS system takes has not been measured anywhere. Measure the
first chunk and lower it:

```
{HOW_TO_MEASURE_TIME}
```

The script times each production chunk
and keeps the last measurement in `system/namd/.build_md_job_chunk_seconds`, so that a
re-submitted job is guarded from its first chunk. The measurement is printed in the job's own
output (`<job-name>_<jobid>_<task>.out` in the folder you submit from) as
`[build_md_job] ... finished <input> in N s`; NAMD's own log of each step is written next to the
inputs, as `system/namd/<prefix>.out`, never `.log`. The guard allows the next chunk 20 % more than
the last one took; that margin was chosen here and is not a measurement. The FIRST chunk of the
first submission runs unguarded, because nothing has been timed yet.

## Inputs edited by this tool
{edited}

## What the simulation contains that the crystal did not
Residues {runs} were modelled by CHARMM-GUI, not located by the experiment. Any report must say so
wherever they are discussed.

{md_package.CALCULATED}

The receptor coordinates are experimental and the ligand pose was calculated by AutoDock Vina;
everything CHARMM-GUI added (the modelled residues, the water box, the ions, the force-field
parameters) is calculated or chosen, not measured. The values of the system are reported without
any threshold applied to them.

Checks that could not run on this archive:

{not_checked}
"""


# ------------------------------------------------------------------ build_md_job

def verified_systems(manifest: dict, prep: Path) -> list:
    """Every system in the manifest whose archive passed the checks and whose files are there.

    This is what the shared array script is built from: a folder is listed only when its
    `job.json` records `verified.passes` and its `system/namd` is on disk, so submitting
    run_all.slurm cannot start a system that nothing verified.
    """
    out = []
    for s in manifest.get("structures") or []:
        for pocket_key, block in sorted((s.get("md_job") or {}).items()):
            rel = block.get("dir")
            if not rel:
                continue
            run_dir = prep / rel
            record = _read_json(run_dir / "job.json")
            if not record or not ((record.get("verified") or {}).get("passes")):
                continue
            if not (run_dir / "system" / "namd").is_dir():
                continue
            # a job.json written by another version, or half written, must not take down a call
            # that has already written everything: it is simply not listed
            if not (record.get("scheme") and record.get("inputs")):
                continue
            try:
                out.append(system_entry(
                    f"{Path(rel).parent.name}/{run_dir.name}/system/namd",
                    f"{s.get('pdb_id')} {pocket_key}", record["scheme"], record["inputs"]))
            except (KeyError, TypeError):
                log.warning("md_job record at %s is not readable; left out of run_all.slurm", rel)
    return out


def _read_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def write_array_script(prep: Path, manifest: dict, site: dict) -> dict:
    """(Re)write prep/md/run_all.slurm with every verified system, one array task each.

    Regenerated on every call instead of appended to, so a system whose files are gone or whose
    archive no longer passes simply stops being listed.
    """
    systems = verified_systems(manifest, prep)
    if not systems:
        return {"written": False, "reason": "no verified system is recorded in the manifest"}
    path = prep / "md" / "run_all.slurm"
    # the SAME resolved site as the single-system script: re-deriving it turned a value the user
    # chose, or a fallback this code invented, into one marked REPORTED with a command named
    # beside it, which is the exact falsehood this design exists to prevent
    limit = site["fields"]["max_array_size"]["value"]
    if limit is not None and len(systems) > limit:
        return {"written": False,
                "reason": f"{len(systems)} systems need an array of that size and the site "
                          f"reports MaxArraySize={limit}; submit them in several batches."}
    path.write_text(slurm_text("md_all", systems, site), encoding="utf-8", newline="\n")
    return {"written": True, "path": str(path), "systems": [s["label"] for s in systems],
            "array": f"0-{len(systems) - 1}",
            "header_is": "the same one as run_md.slurm: the site values are not worked out again",
            "submit_from": str(prep / "md"),
            "note": "One array task per system. Copy the whole folder so that the relative paths "
                    "in SYSTEMS still resolve, and submit from the folder that holds it."}


def _record_in_manifest(path: Path, manifest: dict, pdb_id: str, pocket_key: str, job: dict,
                        run_dir: Path, prep: Path):
    """Write an `md_job` block beside `md_package`, re-reading the manifest first (as dock and
    build_md_package do: if prepare_structures rebuilt it meanwhile, the files stay and this
    says so)."""
    try:
        current = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return f"The job is in {run_dir} but the manifest could not be re-read ({exc!r})."
    if current.get("created_at") != manifest.get("created_at"):
        return (f"The job is in {run_dir} but prepare_structures rewrote the manifest while it was "
                "being written, so it was not recorded there. Run build_md_job again.")
    try:
        s = next(x for x in current["structures"] if x.get("pdb_id") == pdb_id)
        block = s.setdefault("md_job", {})
        block[pocket_key] = {k: job[k] for k in ("generated_by", "generated_at", "archive",
                                                 "job_id", "job_folder", "site",
                                                 "verified", "scheme", "inputs", "inputs_edited",
                                                 "slurm", "to_confirm", "files")}
        block[pocket_key]["dir"] = run_dir.relative_to(prep).as_posix()
        choice.write_manifest(path, current)
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        log.exception("Could not write %s", path)
        return f"The job is in {run_dir} but the manifest could not be written ({exc!r})."
    return None


def build_md_job(uniprot_id: str, archive: str, structure: str = "", pocket: str = "",
                 gpu_resident: bool = False, site_report: str = "", partition: str = "",
                 account: str = "", time_limit: str = "", cpus_per_task: int = 0, gpus: int = 1,
                 gpu_type: str = "", namd_module: str = "", namd_binary: str = "") -> dict:
    """Write the SLURM job folder for a checked CHARMM-GUI archive. See the module docstring.

    Without `site_report` nothing about the cluster is known, so no script is written and the
    commands to paste come back instead.
    """
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    archive_path = Path(archive).expanduser()
    if not archive_path.is_file():
        return {"ok": False, "uniprot_id": accession, "error": f"No file at {archive_path}."}
    path = choice.manifest_path(accession)
    prep = path.parent
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"ok": False, "uniprot_id": accession,
                "error": f"No readable manifest at {path}: run prepare_structures first."}
    staging, pkg_dir = None, None
    try:
        s, snap, origin = dock._resolve(manifest, structure, pocket)
        pdb_id = s["pdb_id"].upper()
        pocket_key = f"{snap['pocket']['comp_id']} {snap['pocket']['label_asym_id']}"
        pkg_dir = prep / "md" / dock.dock_dir_name(pdb_id, snap["pocket"]["comp_id"],
                                                   snap["pocket"]["label_asym_id"])
        pkg_file = pkg_dir / "package.json"
        if not pkg_file.is_file():
            raise MdRunError(f"No package at {pkg_file}: run build_md_package first, then CHARMM-GUI.")
        package = json.loads(pkg_file.read_text(encoding="utf-8"))
        delivered = md_package._inside(prep, package["files"][f"{pdb_id}.pdb"]["path"])

        # the archive is extracted once and the step-4 checks run on that same tree, so the job
        # is written from the files that were verified, not from a second reading
        staging = Path(tempfile.mkdtemp(prefix="md_run_"))
        job = extract_archive(archive_path, staging)
        verified = md_package._check_job(job, delivered, package)
        if not verified["passes"]:
            return {"ok": False, "uniprot_id": accession, "archive": str(archive_path),
                    "error": "No job was written: this archive does not pass the checks of "
                             "check_charmm_gui_output, so the system must not be simulated.",
                    "failures": verified["failures"], "not_checked": verified["not_checked"],
                    "verdict": verified["verdict"], "checks": verified["checks"]}

        if not (site_report or "").strip():
            return {"ok": False, "uniprot_id": accession, "archive": str(archive_path),
                    "error": "No job was written: nothing is known about the cluster this will "
                             "run on, and no value of its configuration will be guessed.",
                    "archive_passed_the_checks": True,
                    "verified": {"passes": True, "verdict": verified["verdict"],
                                 "not_checked": verified["not_checked"]},
                    **site_questions("build_md_job was called without site_report")}

        namd = job / "namd"
        scheme = readme_scheme(namd / "README")
        values = input_values(namd, scheme)
        edits = []
        if gpu_resident:
            edits = [enable_gpu_resident(namd / values["production_input"])]
            # the insertion moved every line below it: re-read, or the tool would name line 14
            # for `run` when the keyword it just wrote is there
            values = input_values(namd, scheme)

        job_name = f"{pdb_id}_{snap['pocket']['comp_id']}_{snap['pocket']['label_asym_id']}_md"
        run_dir = pkg_dir / "run"
        # built beside the final folder and swapped in only when it is complete: a call that
        # fails half way must not destroy a verified job that was already there (build_md_package
        # stages for the same reason)
        new_dir = pkg_dir / "run.new"
        if new_dir.is_dir():
            shutil.rmtree(new_dir)
        new_dir.mkdir(parents=True)
        shutil.move(str(job), str(new_dir / "system"))
        mine = system_entry("system/namd", f"{pdb_id} {pocket_key}", scheme, values)
        parsed = parse_site_report(site_report)
        site = resolve_site(parsed, 1, partition=partition, account=account,
                            time_limit=time_limit, cpus_per_task=cpus_per_task, gpus=gpus,
                            gpu_type=gpu_type, namd_module=namd_module, namd_binary=namd_binary)
        (new_dir / "run_md.slurm").write_text(slurm_text(job_name, [mine], site),
                                              encoding="utf-8", newline="\n")
        record = {
            "generated_by": "build_md_job", "generated_at": utc_now(), "uniprot_id": accession,
            "pdb_id": pdb_id, "pocket": pocket_key, "selected_by": origin,
            "archive": {"path": str(archive_path), "sha256": dock.sha256(archive_path)},
            # the folder name as it comes, and the numeric id only when the name really carries
            # one: "charmm-gui-9418739369" is a folder, 9418739369 is the job
            "job_folder": verified["checks"].get("job_id"),
            "job_id": _job_id(verified["checks"].get("job_id")),
            "package": {"dir": pkg_dir.relative_to(prep).as_posix(),
                        "run_id": package.get("run_id"), "pose": package.get("pose"),
                        "to_model": package.get("to_model")},
            "verified": {"passes": verified["passes"], "verdict": verified["verdict"],
                         "not_checked": verified["not_checked"], "checks": verified["checks"]},
            "scheme": scheme,
            "inputs": values,
            # the simulation and the cluster are kept apart on purpose: different owner,
            # different provenance, different way of being wrong
            "site": site["fields"],
            "site_report_unparsed": site["unparsed"],
            "inputs_edited": edits,
            "slurm": {"job_name": job_name, "chunks": scheme["chunks"],
                      "nanoseconds_total": round(scheme["chunks"] * values["nanoseconds_per_chunk"], 4),
                      "array": "0-0",
                      "site_values_from": "the output the user pasted; every field in `site` "
                                          "carries its own mark and source"},
            "to_confirm": site["to_confirm"],
            "files": {},
        }
        (new_dir / "README.md").write_text(
            readme_text(package, scheme, values, verified, job_name, run_dir, edits, site),
            encoding="utf-8")
        record["files"] = {p.name: {"path": (run_dir / p.name).relative_to(prep).as_posix(),
                                    "sha256": dock.sha256(p)}
                           for p in (new_dir / "run_md.slurm", new_dir / "README.md")}
        (new_dir / "job.json").write_text(json.dumps(record, ensure_ascii=False, indent=1),
                                          encoding="utf-8")
        # the swap is the last thing that happens, so everything above can fail harmlessly
        replaced = run_dir.is_dir()
        old_dir = pkg_dir / "run.previous"
        if replaced:
            if old_dir.is_dir():
                shutil.rmtree(old_dir)
            run_dir.rename(old_dir)
        new_dir.rename(run_dir)
        if replaced:
            shutil.rmtree(old_dir, ignore_errors=True)
        manifest_note = _record_in_manifest(path, manifest, pdb_id, pocket_key, record, run_dir, prep)
        # the shared array script is rebuilt from the manifest as it stands AFTER this run was
        # recorded, so it lists exactly the verified systems and nothing else
        array = write_array_script(prep, _read_json(path) or manifest, site)
    except (MdRunError, dock.DockError, md_package.MdPackageError) as exc:
        # MdPackageError comes out of _check_job when the archive is not a finished job: it is a
        # refusal with a readable reason, not a crash, and the docstring promises ok: false
        if pkg_dir:
            shutil.rmtree(pkg_dir / "run.new", ignore_errors=True)
        return {"ok": False, "uniprot_id": accession, "error": str(exc), **getattr(exc, "extra", {})}
    except (OSError, ValueError, KeyError) as exc:
        log.exception("build_md_job failed")
        if pkg_dir:
            shutil.rmtree(pkg_dir / "run.new", ignore_errors=True)
        return {"ok": False, "uniprot_id": accession,
                "error": f"Job not written: {exc!r}",
                "note": "Nothing was replaced: a job folder that was already there is untouched."}
    finally:
        if staging:
            shutil.rmtree(staging, ignore_errors=True)

    notes = [
        "Nothing was simulated and nothing was sent: the transfer and sbatch are manual, like the "
        "CHARMM-GUI upload.",
        f"run_md.slurm aborts with exit code 2 while any CONFIRM placeholder is left in it "
        "a commented ##SBATCH line is inert and does not block it.",
f"--time is {site['fields']['time']['value']} ({site['fields']['time']['source']}): a "
        "ceiling, not an estimate. How long THIS system takes has not been measured anywhere. "
        f"Measure the first chunk with `{HOW_TO_MEASURE_TIME}` and lower it, or read the line "
        "`[build_md_job] ... finished ... in N s` the script prints in the job output. NAMD's own "
        "logs are system/namd/<prefix>.out, not .log. The script times each production chunk and "
        "keeps the last measurement, so a re-submitted job is guarded from its first chunk; the "
        "very first chunk of the first submission runs unguarded.",
        f"The chunk length is read from the `run` line of {values['production_input']}, not from "
        "`numsteps`: NAMD executes the steps of `run`, and the real inputs carry both.",
        f"{scheme['chunks']} chunks x {values['nanoseconds_per_chunk']} ns comes from cntmax in the "
        "archive's own namd/README; it was not chosen by this tool.",
    ]
    if edits:
        notes.insert(0, "GPU-resident integration was added to " + values["production_input"]
                     + f" ('{GPU_RESIDENT_LINE}'). That input is no longer byte-identical to what "
                       "CHARMM-GUI returned, and the equilibration was left untouched because "
                       "whether GPU-resident NAMD supports its positional restraints is not verified.")
    if verified["not_checked"]:
        notes.insert(0, f"{len(verified['not_checked'])} check(s) could not run on this archive; "
                        "they are in not_checked and are listed in README.md too.")
    if values["length_warning"]:
        notes.insert(0, values["length_warning"])
    if array["written"]:
        notes.append(f"To run several systems at once there is now {array['path']}: a SLURM array "
                     f"with one task per system ({', '.join(array['systems'])}), rebuilt on every "
                     "call from the manifest, so it lists only systems whose archive passed the "
                     "checks.")
    notes.append("Every #SBATCH line is marked REPORTED (read from the output you pasted, with "
                 "the command named), CHOSEN (your decision), CONFIRM (not known; the script will "
                 "not run until it is filled in) or EDITABLE (a knob, with its default's origin "
                 "stated). `site` lists them field by field.")
    notes.append("What your cluster reports is not a promise that the job runs: whether NAMD works "
                 "there, and whether the site wants --gres=gpu or --gpus, only the first sbatch "
                 "will tell you.")
    notes.extend(site["notes"])
    if replaced:
        notes.append("A previous run folder for this pocket was replaced.")
    if manifest_note:
        notes.insert(0, manifest_note)
    modelled = package.get("to_model") or []
    return {
        "ok": True, "uniprot_id": accession, "run_dir": str(run_dir),
        "job_id": record["job_id"], "job_folder": record["job_folder"],
        "built": f"{pdb_id} pocket {pocket_key}, CHARMM-GUI job "
                 f"{record['job_id'] or record['job_folder']}",
        # from the measurement, never a literal: the verdict and this flag cannot disagree
        "verified": {"passes": verified["passes"], "verdict": verified["verdict"],
                     "not_checked": verified["not_checked"]},
        "to_confirm": site["to_confirm"],
        # the cluster and the simulation are reported apart, always: different owner, different
        # provenance, different way of being wrong
        "site": site["fields"],
        "site_report_unparsed": site["unparsed"],
        "simulation": {
            "chunks": scheme["chunks"], "nanoseconds_per_chunk": values["nanoseconds_per_chunk"],
            "nanoseconds_total": round(scheme["chunks"] * values["nanoseconds_per_chunk"], 4),
            "chunks_from": scheme["source"],
            "equilibration": f"{values['equilibration_input']} as CHARMM-GUI wrote it"
                             + (f" (minimize {values['equilibration_minimize']['value']}, run "
                                f"{values['equilibration_run']['value']} steps)"
                                if values["equilibration_minimize"] and values["equilibration_run"] else ""),
        },
        "edit_to_change_the_length": {
            "more_nanoseconds": f"CHUNKS in {run_dir / 'run_md.slurm'} (now {scheme['chunks']})",
            "longer_chunks": f"in {run_dir / 'system' / 'namd' / values['production_input']}, change "
                             "BOTH lines together, because NAMD executes the steps of `run`: "
                             + " and ".join(f"line {line['line']} ({line['text']})"
                                            for line in (values["run_line"], values["numsteps_line"])
                                            if line)
                             + f". {values['steps_per_chunk']['value']} steps x "
                               f"{values['timestep_fs']['value']} fs = "
                               f"{values['nanoseconds_per_chunk']} ns per chunk; these are steps, "
                               "not seconds",
            "warning": "Editing the NAMD input means the system is no longer the one CHARMM-GUI "
                       "returned and check_charmm_gui_output verified; say so in any report.",
        },
        # neither line was tried from this machine: they are what to copy, not a verified command
        "transfer": {"folder": str(run_dir),
                     "scp": f"scp -r \"{run_dir}\" <user>@<cluster>:<path>/{job_name}",
                     "rsync_from_wsl_or_linux":
                         f"rsync -av \"{run_dir.as_posix()}/\" <user>@<cluster>:<path>/{job_name}/",
                     "note": "Neither command was run from here; on Windows rsync usually needs WSL "
                             "and a /mnt/c path."},
        "submit": f"cd <path>/{job_name} && sbatch run_md.slurm",
        "inputs_edited": edits,
        "run_several_at_once": array,
        "recorded_in_manifest": manifest_note is None,
        "modelled_not_experimental": (
            ", ".join(f"{r['chain']} {r['first']}-{r['last']}" for r in modelled)
            + " were modelled by CHARMM-GUI, not located by the experiment; say so wherever a "
              "report discusses them." if modelled else
            "No run was modelled: every residue of the system is as it was delivered."),
        "next": ("Fill in what is still CONFIRM, copy the folder to the cluster and submit it. "
                 "Nothing here reaches it." if site["to_confirm"] else
                 "Copy the folder to the cluster and submit it. Nothing here reaches it."),
        "notes": notes,
    }


# ------------------------------------------------------------------ CLI

def main(argv=None):
    """CLI: the same function as the MCP tool, result as JSON (separate process, stdout is free)."""
    ap = argparse.ArgumentParser(prog="python -m protein_mcp.files.md_run",
                                 description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build", help="write the SLURM job folder for a checked archive")
    build.add_argument("uniprot_id")
    build.add_argument("--archive", required=True)
    build.add_argument("--structure", default="")
    build.add_argument("--pocket", default="")
    build.add_argument("--gpu-resident", dest="gpu_resident", action="store_true",
                       help="insert CUDASOAintegrate on into the production input (declared edit)")
    build.add_argument("--site-report", default="", help="path to a file holding the output of the "
                                                         "four commands, or '-' to read stdin")
    build.add_argument("--partition", default="")
    build.add_argument("--account", default="")
    build.add_argument("--time-limit", dest="time_limit", default="")
    build.add_argument("--cpus-per-task", dest="cpus_per_task", type=int, default=0)
    build.add_argument("--gpus", type=int, default=1)
    build.add_argument("--gpu-type", dest="gpu_type", default="", help="the GPU model to pin, when "
                                                                      "the partition has several")
    build.add_argument("--namd-module", dest="namd_module", default="")
    build.add_argument("--namd-binary", dest="namd_binary", default="")
    a = ap.parse_args(argv)
    report = ""
    if a.site_report == "-":
        report = sys.stdin.read()
    elif a.site_report:
        report = Path(a.site_report).expanduser().read_text(encoding="utf-8", errors="replace")
    result = build_md_job(a.uniprot_id, a.archive, a.structure, a.pocket, a.gpu_resident, report,
                          a.partition, a.account, a.time_limit, a.cpus_per_task, a.gpus,
                          a.gpu_type, a.namd_module, a.namd_binary)
    sys.stdout.buffer.write((json.dumps(result, ensure_ascii=False, indent=1) + "\n").encode("utf-8"))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
