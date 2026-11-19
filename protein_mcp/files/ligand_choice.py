"""select_ligands: the researcher's compound selection, written to prep/manifest.json (Phase 3, step 6b).

Mirror of files.choice.select_structure, for compounds instead of structures. The decision of
2026-11-15 stands: the selection is the researcher's, and this tool supports and RECORDS it, it
does not replace it with a ranking. ESTADO 7.4, decision (b) 1: the choice and its reason are
written BEFORE anything is docked, so the record cannot be assembled after seeing the scores.

No network: a compound arrives as an id plus its SMILES (whoever retrieved it is who has the
network; this module only validates and stores). What is checked, and refused with the reason:
the SMILES must be readable by RDKit, have one component (no '.'), and have at least one heavy
atom; ids must be unique; the level must be T1 or T3 (T2, a declared inactive sample, is out of
step 6b by the user's decision (a) 2); and a reason is required, as in select_structure.

What is stored per compound: the SMILES as given, RDKit's canonical SMILES, the InChIKey, heavy
atoms, formal charge, molecular weight, rotatable bonds and the slug that dock.run_id_of uses.
Nothing is computed about activity, and no property here orders anything.
Imports common, dossier and files.choice (its atomic write). Pure function + CLI.
"""
import argparse
import json
import re
import sys

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice

LEVELS = ("T1", "T3")
MAX_REASON = 1000
MAX_COMPOUNDS = 50
STEP6B_PENDING = "step 6b: researcher's compound selection (select_ligands)"
LEVEL_MEANING = {
    "T1": ("Cross-docking: the compound has a crystallographic pose in another structure of the same target, "
           "so an RMSD against that pose is possible once the superposition exists."),
    "T3": ("Researcher's compound: no crystallographic pose in this target, so no RMSD is possible and the "
           "score is all there is."),
}
NOT_A_RANKING = ("This is the researcher's selection with their own reason, recorded before docking. It is not a "
                 "ranking and no property stored here orders the compounds.")


def slug_of(compound_id: str, inchikey: str) -> str:
    """A short file-name-safe identity for the run id: the id (ChEMBL ids are ideal) or, when it has
    nothing usable, the head of the InChIKey. Never empty, or two compounds would share a run id."""
    s = re.sub(r"[^a-z0-9]+", "", (compound_id or "").lower())[:20]
    return s or re.sub(r"[^a-z0-9]+", "", (inchikey or "").lower())[:12]


def _one(raw, seen: set) -> dict:
    """One validated compound record, or ValueError with the reason (the caller names the row)."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    if not isinstance(raw, dict):
        raise ValueError("each compound must be an object with at least 'id', 'smiles' and 'level'")
    cid = str(raw.get("id") or "").strip()
    smiles = str(raw.get("smiles") or "").strip()
    level = str(raw.get("level") or "").strip().upper()
    name = str(raw.get("name") or "").strip() or None
    if not cid:
        raise ValueError("'id' is required (a ChEMBL id, or any identifier the user uses)")
    if not smiles:
        raise ValueError("'smiles' is required: this module has no network and does not resolve ids")
    if level not in LEVELS:
        raise ValueError(f"'level' must be one of {', '.join(LEVELS)} (T2, a declared inactive sample, is "
                         f"out of step 6b by the user's decision); got {level or None!r}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not read the SMILES {smiles!r}")
    parts = smiles.split(".")
    if len(parts) > 1:
        raise ValueError(f"the SMILES has {len(parts)} components (a salt or a mixture); give the single "
                         "component that is docked")
    if mol.GetNumHeavyAtoms() == 0:
        raise ValueError("the SMILES has no heavy atom")
    try:                                  # the InChI library is optional in some RDKit builds
        inchikey = Chem.MolToInchiKey(mol) or None
    except Exception:                     # noqa: BLE001 - absence of InChI must not refuse a valid compound
        log.warning("RDKit could not make an InChIKey for %s", cid)
        inchikey = None
    slug = slug_of(cid, inchikey or "")
    if not slug:
        raise ValueError(f"no file-safe identity could be made from id {cid!r} and its InChIKey")
    if slug in seen:
        raise ValueError(f"the identity {slug!r} is already used by another compound of this selection: two "
                         "compounds would share one docking run id")
    seen.add(slug)
    return {
        "id": cid, "name": name, "level": level, "slug": slug,
        "smiles": smiles, "canonical_smiles": Chem.MolToSmiles(mol), "inchikey": inchikey,
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge(mol),
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "level_meaning": LEVEL_MEANING[level],
        "values": "calculated",
        "calculated_by": "RDKit from the SMILES given; nothing here is a measurement.",
    }


def carry_over(old: dict, new: dict):
    """Called by prepare_structures before it writes `new` over `old`: the compound selection does
    not depend on any structure, so it is kept as it is (with the date of the revalidation). It is
    the box and the pocket that can change, and those live in the docking runs, not here."""
    previous = (old or {}).get("ligand_choice")
    new["ligand_choice"] = None
    new["ligand_choice_history"] = list((old or {}).get("ligand_choice_history") or [])
    if not previous:
        return None
    new["ligand_choice"] = {**previous, "revalidated_at": new.get("created_at")}
    return {"kept": f"{len(previous.get('compounds') or [])} compounds",
            "note": "The compound selection does not depend on the structures, so it was kept unchanged."}


def select_ligands(uniprot_id: str, compounds=None, reason: str = "") -> dict:
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    reason = (reason or "").strip()
    if isinstance(compounds, str):                 # a JSON string from a CLI or a chat client
        try:
            compounds = json.loads(compounds)
        except ValueError as exc:
            return {"ok": False, "uniprot_id": accession, "error": f"compounds is not valid JSON: {exc}"}
    if isinstance(compounds, dict):
        compounds = [compounds]
    if not isinstance(compounds, list) or not compounds:
        return {"ok": False, "uniprot_id": accession,
                "error": "compounds is required: a list of objects with 'id', 'smiles', 'level' (T1 or T3) "
                         "and optionally 'name'."}
    if len(compounds) > MAX_COMPOUNDS:
        return {"ok": False, "uniprot_id": accession,
                "error": f"{len(compounds)} compounds: more than the {MAX_COMPOUNDS} this tool accepts in one call."}
    path = choice.manifest_path(accession)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {"ok": False, "uniprot_id": accession,
                "error": f"No preparation manifest for {accession}: run prepare_structures first."}
    except (OSError, ValueError) as exc:
        log.exception("Could not read %s", path)
        return {"ok": False, "uniprot_id": accession, "error": f"Could not read {path}: {exc}"}

    records, errors, seen = [], [], set()
    for i, raw in enumerate(compounds, 1):
        try:
            records.append(_one(raw, seen))
        except ValueError as exc:
            name = (raw.get("id") if isinstance(raw, dict) else None) or f"row {i}"
            errors.append(f"{name}: {exc}")
    # nothing is recorded if one row is wrong: a selection half written is not the user's selection
    if errors:
        return {"ok": False, "uniprot_id": accession,
                "error": f"Not recorded: {len(errors)} of {len(compounds)} compounds were refused.",
                "refused": errors,
                "accepted_so_far": [r["id"] for r in records]}
    # the selection is validated first, as in select_structure (a missing reason used to hide a real error)
    if not reason:
        return {"ok": False, "uniprot_id": accession,
                "valid_selection": [f"{r['id']} ({r['name'] or 'no name'}, {r['level']})" for r in records],
                "error": "Not recorded: reason is required. The compounds above are valid; ask the user for "
                         "their own reason for choosing them, do not write one."}
    if len(reason) > MAX_REASON:
        return {"ok": False, "uniprot_id": accession, "error": f"reason is longer than {MAX_REASON} characters."}

    now = utc_now()
    selection = {
        "compounds": records, "reason": reason, "chosen_at": now, "recorded_by": "select_ligands",
        "manifest_created_at": manifest.get("created_at"),
        "levels": {lv: [r["id"] for r in records if r["level"] == lv] for lv in LEVELS},
        "note": NOT_A_RANKING,
    }
    previous = manifest.get("ligand_choice")
    history = list(manifest.get("ligand_choice_history") or [])
    if previous:
        history.append({**previous, "superseded_at": now})
    pending = [p for p in (manifest.get("pending") or []) if p != STEP6B_PENDING]
    manifest.update(ligand_choice=selection, ligand_choice_history=history, pending=pending)
    try:
        choice.write_manifest(path, manifest)
    except OSError as exc:
        log.exception("Could not write %s", path)
        return {"ok": False, "uniprot_id": accession, "error": f"The manifest could not be written: {exc!r}"}
    return {
        "ok": True, "uniprot_id": accession, "manifest_path": str(path),
        "recorded": f"{len(records)} of {len(compounds)} compounds given",
        "compounds": [{k: r[k] for k in ("id", "name", "level", "slug", "canonical_smiles", "inchikey",
                                         "heavy_atoms", "formal_charge", "molecular_weight", "rotatable_bonds")}
                      for r in records],
        "levels": selection["levels"],
        "reason": reason,
        "heavy_atom_range": [min(r["heavy_atoms"] for r in records), max(r["heavy_atoms"] for r in records)],
        "replaced_selection": (f"{len(previous.get('compounds') or [])} compounds recorded "
                               f"{previous.get('chosen_at')}") if previous else None,
        "previous_selections": len(history),
        "next": (f"dock_ligand('{accession}', compound='<id>') docks one of them into the chosen pocket, with "
                 "that pocket's box unchanged. One compound per call, one run at a time on this machine."),
        "notes": [
            NOT_A_RANKING,
            "Validated with RDKit from the SMILES given: single component, readable, with heavy atoms. "
            "Whether the SMILES is the compound that was measured is not something this tool can check; a "
            "salt or a different protonation state is the expected failure mode of an exact-SMILES match.",
            "Properties stored here are calculated by RDKit, never measured, and none of them ranks anything.",
            "If prepare_structures is run again, the selection is kept: it does not depend on the structures.",
        ],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m protein_mcp.files.ligand_choice",
                                 description=__doc__.splitlines()[0])
    ap.add_argument("uniprot_id")
    ap.add_argument("--compounds", required=True, help="JSON list of objects with id, smiles, level and name, "
                                                       "or @path to a file with that JSON")
    ap.add_argument("--reason", default="", help="the user's own reason for this selection")
    a = ap.parse_args(argv)
    text = a.compounds
    if text.startswith("@"):
        with open(text[1:], encoding="utf-8") as fh:
            text = fh.read()
    result = select_ligands(a.uniprot_id, text, a.reason)
    sys.stdout.buffer.write((json.dumps(result, ensure_ascii=False, indent=1) + "\n").encode("utf-8"))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
