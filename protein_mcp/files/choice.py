"""select_structure: the user's structure choice, written to prep/manifest.json (Phase 2, step 4).

The choice is checked against the manifest (structure, chains and pocket must exist there) and
stored with a copy of the manifest data Phase 3 needs (receptor, chains, pocket, box, SDF). Nothing
is computed here. prepare_structures calls carry_over() when it replaces the manifest: the choice
is kept if it is still valid in the new manifest, otherwise it moves to user_choice_history with
the reason. Imports common and dossier.
"""
import json

from protein_mcp import dossier
from protein_mcp.common import log, parse_accession, utc_now

STEP4_PENDING = "step 4: user's structure choice (select_structure)"
MAX_REASON = 1000
# fields that record the choice itself; everything else is copied from the manifest
_META = ("reason", "chosen_at", "recorded_by", "manifest_created_at", "revalidated_at", "changed_fields", "warnings")


def manifest_path(accession: str):
    return dossier.DOSSIER_DIR / accession / "prep" / "manifest.json"


def write_manifest(path, manifest: dict):
    """Atomic write (a .part file replaced in one step): a reader never sees half a manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
    part.replace(path)


def pending(items, has_choice: bool) -> list:
    rest = [p for p in items or [] if p != STEP4_PENDING]
    return rest if has_choice else [STEP4_PENDING] + rest


def label(choice: dict) -> str:
    """'3OG7 chain A, pocket 032 C' (for replies and the report)."""
    if not choice:
        return "none"
    chosen = [c["auth_asym_id"] for c in choice.get("chains") or []]
    pk = choice.get("pocket")
    text = (choice["structure"] + (f" chain{'s' if len(chosen) > 1 else ''} {', '.join(chosen)}" if chosen else "")
            + (f", pocket {pk['comp_id']} {pk['label_asym_id']}" if pk else ", no pocket"))
    # Desktop test 2026-10-28: "5JSM chain A, pocket 6NB J" read as a chain-A-only pocket; its box covers A and B
    spans = list((pk or {}).get("residues_by_chain") or {})
    if len(spans) > 1 and set(spans) - set(chosen):
        text += f" (pocket spans chains {', '.join(spans)}; box not reduced)"
    return text


def _pocket_name(pk: dict) -> str:
    return f"{pk['comp_id']} {pk['label_asym_id']}"


def options(manifest: dict) -> dict:
    """What can be chosen: structures with a cleaned receptor, their target chains (auth) and
    pockets ('comp label: chains'). Structures without a receptor are listed apart."""
    out, unusable = {}, {}
    for s in manifest.get("structures") or []:
        prep = s.get("prepared") or {}
        if not prep.get("ok"):
            unusable[s.get("pdb_id")] = prep.get("error") or (s.get("file") or {}).get("error") or "no cleaned receptor"
            continue
        out[s["pdb_id"]] = {
            "chains": [c["auth_asym_id"] for c in prep["metrics"]["chains"]],
            "pockets": [f"{_pocket_name(p)}: chain{'s' if len(p['residues_by_chain']) > 1 else ''} "
                        + ", ".join(p["residues_by_chain"]) for p in prep["pockets"]
                        if not p.get("error") and p.get("residues_by_chain")],
        }
    af = manifest.get("alphafold") or {}
    if (af.get("file") or {}).get("ok"):
        out["AlphaFold"] = {"entry_id": af.get("entry_id"), "chains": [], "pockets": []}
    return {"structures": out, "not_selectable": unusable}


def _find_pocket(structure: dict, text: str):
    """Pocket from 'comp label', 'comp' (if unique) or '<PDB> comp label' (the summary format)."""
    parts = text.split()
    if len(parts) == 3 and parts[0].upper() == structure["pdb_id"]:
        parts = parts[1:]
    if not parts or len(parts) > 2:
        return None, f"pocket '{text}' is not '<comp_id> <label_asym_id>' (e.g. '032 C')"
    comp, inst = parts[0].upper(), (parts[1] if len(parts) == 2 else None)
    pockets = structure["prepared"]["pockets"]
    hits = [p for p in pockets if p["comp_id"].upper() == comp and (inst is None or p["label_asym_id"] == inst)]
    if not hits:
        return None, f"no pocket '{text}' in {structure['pdb_id']}"
    if len(hits) > 1:
        return None, (f"'{text}' matches {len(hits)} pockets in {structure['pdb_id']} ("
                      + ", ".join(_pocket_name(p) for p in hits) + "): give the instance too")
    pk = hits[0]
    if pk.get("error") or not pk.get("residues_by_chain"):
        return None, f"pocket {_pocket_name(pk)} of {structure['pdb_id']} has no residues: {pk.get('error') or 'none within the cutoff'}"
    return pk, None


def resolve(manifest: dict, structure: str, pocket: str = "", chains=None):
    """(snapshot, None) or (None, error text). The snapshot is copied from the manifest."""
    structure = (structure or "").strip()
    pocket = (pocket or "").strip()
    wanted = []
    for c in chains or []:
        c = str(c).strip()
        if c and c not in wanted:
            wanted.append(c)
    af = manifest.get("alphafold") or {}
    if structure.lower() in ("alphafold", str(af.get("entry_id") or "").lower()) and structure:
        if not (af.get("file") or {}).get("ok"):
            return None, "the manifest has no AlphaFold file (not requested or its download failed)"
        if pocket or wanted:
            return None, ("the AlphaFold model has no ligand, so no pocket; and its chain is not chosen here: "
                          "choose AlphaFold without pocket and chains")
        return {
            "structure": "AlphaFold", "kind": "alphafold", "entry_id": af.get("entry_id"),
            "is_canonical_sequence": af.get("is_canonical_sequence"), "mean_plddt": af.get("mean_plddt"),
            "file": af["file"].get("path"), "file_sha256": af["file"].get("sha256"),
            "variant_position": af.get("variant_position"), "chains": [], "pocket": None,
            "warnings": ["AlphaFold model chosen: no ligand, no pocket and no box (a box would have to be "
                         "defined in Phase 3); the file is the raw prediction, not cleaned."],
        }, None

    s = next((x for x in manifest.get("structures") or [] if x.get("pdb_id") == structure.upper()), None)
    if s is None:
        return None, f"structure '{structure}' is not in the manifest"
    prep = s.get("prepared") or {}
    if not prep.get("ok"):
        return None, (f"{s['pdb_id']} has no cleaned receptor in the manifest: "
                      + (prep.get("error") or (s.get("file") or {}).get("error") or "not prepared"))
    target = {c["auth_asym_id"]: c for c in prep["metrics"]["chains"]}
    bad = [c for c in wanted if c not in target]
    if bad:
        return None, (f"chain(s) {', '.join(bad)} not a target chain of {s['pdb_id']} (auth ids: "
                      f"{', '.join(target)}; chain ids are case-sensitive)")
    warnings, pk = [], None
    if pocket:
        pk, err = _find_pocket(s, pocket)
        if err:
            return None, err
        pocket_chains = list(pk["residues_by_chain"])
        outside = [c for c in wanted if c not in pocket_chains]
        if outside:
            return None, (f"chain(s) {', '.join(outside)} hold no residue of pocket {_pocket_name(pk)} "
                          f"(its chains: {', '.join(pocket_chains)})")
        wanted = wanted or pocket_chains
        if pk.get("spans_target_chains") and len(wanted) < len(pocket_chains):
            warnings.append(f"Pocket {_pocket_name(pk)} spans chains {', '.join(pocket_chains)}: its residues and box "
                            f"cover all of them and are NOT reduced to chain(s) {', '.join(wanted)}.")
    elif not wanted:
        if len(target) != 1:
            return None, f"{s['pdb_id']} has {len(target)} target chains ({', '.join(target)}): name the chain(s) or a pocket"
        wanted = list(target)
    if not pk:
        warnings.append("No pocket chosen: no ligand pose and no box for Phase 3.")

    ligand = None
    if pk:
        inst = next((i for i in s.get("nonpolymer_instances") or [] if i.get("label_asym_id") == pk["label_asym_id"]), {})
        sdf = inst.get("sdf") or {}
        ligand = {"name": inst.get("name"), "synonyms": inst.get("synonyms") or [],
                  "formula_weight": inst.get("formula_weight"), "sdf": sdf.get("path") if sdf.get("ok") else None}
        if not sdf.get("ok"):
            warnings.append(f"Ligand {_pocket_name(pk)} has no SDF in the manifest"
                            + (f" ({sdf['error']})" if sdf.get("error") else "") + ": its crystal pose is not available.")
    return {
        "structure": s["pdb_id"], "kind": "experimental",
        "resolution_angstrom": s.get("resolution_angstrom"), "r_free": prep["metrics"].get("r_free"),
        "raw_file": s["file"].get("path"), "raw_sha256": s["file"].get("sha256"),
        "receptor": prep["receptor"]["path"], "receptor_format": prep["receptor"].get("format"),
        "receptor_chains_auth": prep["receptor"].get("chains_auth"),
        # atoms and kept nucleotides/metals: a change in the cleaning shows up in changed_fields
        "receptor_atoms": prep["receptor"].get("atoms"),
        "receptor_kept_nonpolymer": prep["receptor"].get("kept_nonpolymer_label_asym_ids"),
        # all three lists: an empty `mutations` does not prove a native sequence (4CQE: title + "YES")
        "construct": {k: (s.get("construct") or {}).get(k) or []
                      for k in ("mutations", "mutations_from_title", "mutations_unparsed")},
        "chains": [{
            "auth_asym_id": c, "label_asym_id": target[c]["label_asym_id"],
            "modelled_uniprot_range": target[c].get("modelled_uniprot_range"),
            "missing_uniprot_ranges": target[c].get("missing_uniprot_ranges"),
            **({"variant_status": target[c]["variant_check"].get("status"),
                "variant_observed": target[c]["variant_check"].get("observed")} if "variant_check" in target[c] else {}),
        } for c in wanted],
        "pocket": None if not pk else {
            "comp_id": pk["comp_id"], "label_asym_id": pk["label_asym_id"], "auth_asym_id": pk["auth_asym_id"],
            **ligand,
            "n_residues": len(pk["residues"]), "residues_by_chain": pk["residues_by_chain"],
            "spans_target_chains": bool(pk.get("spans_target_chains")), "box": pk["box"],
            "large_box": bool(pk.get("large_box")), "adjacent_missing_ranges": pk.get("adjacent_missing_ranges") or [],
            "alphafold_mean_plddt": (pk.get("alphafold_plddt") or {}).get("mean_plddt"),
        },
        "warnings": warnings,
    }, None


def _key(choice: dict):
    """The parameters that reproduce a choice: (structure, pocket, chains)."""
    pk = choice.get("pocket")
    return (choice["structure"], _pocket_name(pk) if pk else "",
            [c["auth_asym_id"] for c in choice.get("chains") or []])


def carry_over(old: dict, new: dict):
    """Called by prepare_structures before it writes `new` over `old`. Sets new['user_choice'],
    new['user_choice_history'] and new['pending']. Returns what happened, for the reply (None if
    there was no choice)."""
    history = list((old or {}).get("user_choice_history") or [])
    previous = (old or {}).get("user_choice")
    new["user_choice"], new["user_choice_history"] = None, history
    status = None
    if previous:
        snap, err = resolve(new, *_key(previous))
        if not err and snap.get("raw_sha256", snap.get("file_sha256")) != previous.get(
                "raw_sha256", previous.get("file_sha256")):
            err = "the structure file changed (different SHA-256)"
        if err:
            history.append({**previous, "dropped_at": new.get("created_at"), "dropped_reason": err})
            status = {"dropped": label(previous), "reason": err}
        else:
            changed = [k for k in snap if k not in _META and snap[k] != previous.get(k)]
            meta = {k: previous[k] for k in ("reason", "chosen_at", "recorded_by", "manifest_created_at") if k in previous}
            new["user_choice"] = {**snap, **meta, "revalidated_at": new.get("created_at"), "changed_fields": changed}
            status = {"kept": label(previous), "changed_fields": changed}
    new["pending"] = pending(new.get("pending"), new["user_choice"] is not None)
    return status


def select_structure(uniprot_id: str, structure: str, pocket: str = "", chains=None, reason: str = "") -> dict:
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    reason = (reason or "").strip()
    if not (structure or "").strip():
        return {"ok": False, "uniprot_id": accession, "error": "structure is required (a PDB ID of the manifest or 'AlphaFold')."}
    path = manifest_path(accession)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {"ok": False, "uniprot_id": accession,
                "error": f"No preparation manifest for {accession}: run prepare_structures first."}
    except (OSError, ValueError) as exc:
        log.exception("Could not read %s", path)
        return {"ok": False, "uniprot_id": accession, "error": f"Could not read {path}: {exc}"}

    snap, error = resolve(manifest, structure, pocket, chains)
    if error:
        return {"ok": False, "uniprot_id": accession, "error": f"Not recorded: {error}.",
                "options": options(manifest)}
    # the choice is checked first (Desktop test 2026-10-28: a missing reason hid the chain error)
    if not reason:
        return {"ok": False, "uniprot_id": accession, "valid_choice": label(snap),
                "error": f"Not recorded: reason is required. The choice {label(snap)} is valid; ask the user "
                         "for their own reason, do not write one."}
    if len(reason) > MAX_REASON:
        return {"ok": False, "uniprot_id": accession, "error": f"reason is longer than {MAX_REASON} characters."}
    now = utc_now()
    choice = {**snap, "reason": reason, "chosen_at": now, "recorded_by": "select_structure",
              "manifest_created_at": manifest.get("created_at")}
    previous = manifest.get("user_choice")
    history = list(manifest.get("user_choice_history") or [])
    if previous:
        history.append({**previous, "superseded_at": now})
    manifest.update(user_choice=choice, user_choice_history=history,
                    pending=pending(manifest.get("pending"), True))
    try:
        write_manifest(path, manifest)
    except OSError as exc:
        log.exception("Could not write %s", path)
        return {"ok": False, "uniprot_id": accession, "error": f"The manifest could not be written: {exc!r}"}
    return {
        "ok": True,
        "uniprot_id": accession,
        "recorded": label(choice),
        "manifest_path": str(path),
        "user_choice": choice,
        "replaced_choice": label(previous) if previous else None,
        "previous_choices": len(history),
        "notes": [
            "Copied from the manifest: nothing was recomputed. The receptor is not protonated or completed.",
            "If prepare_structures is run again, the choice is kept only if the structure, chains and pocket "
            "are still in the new manifest with the same file; otherwise it moves to user_choice_history.",
        ],
    }
