"""analyze_admet: ADME records, measurements and calculated properties for the ligand that was
docked and simulated (Phase 3, step 7).

Measured against the real ChEMBL API and the real files of this project on 2026-11-11
(ESTADO 8.5). What that measurement settled:

- **The compound link is an exact InChIKey, computed here with RDKit from the SDF itself.** The
  pose SDF (hydrogens added by Meeko) and the crystallographic SDF from ModelServer give the SAME
  key for both ligands of the project, so the link does not depend on which file is read:
  GPXBXXGIAQBQNI-UHFFFAOYSA-N -> CHEMBL1229517, GCVULUAJPYRJHT-UHFFFAOYSA-N -> CHEMBL5440746,
  total_count=1 each. No compound name is used to search and none is inferred.
- **Three provenances, three blocks, never one table** (rule of ESTADO 7.1): `chembl_record`
  (facts of the ChEMBL entry: phase, approval, withdrawal, routes, ATC), `measured_adme`
  (retrieved activities, each with its own units) and `calculated_properties` (RDKit here and
  ChEMBL's own engine, side by side and never merged).
- **"Predicted" is not used in this version.** No number here comes from a trained model; that is
  step 7.1. The word is reserved for it.
- **The rules are conventions, not measurements** (pre-registered in ESTADO 7.3). They go in
  `rule_reading` with status READING OFFERED, NOT MEASURED, their cut-offs written out, and no
  property carries a verdict of its own.
- **The measured block is usually almost empty, and that is the source, not a bug.** Vemurafenib,
  an approved drug, has 1 CL, 1 Papp, 1 Solubility and 1 permeability in ChEMBL; TXV has none.
  The output says how many rows were read and how many survived the allowlist, so "no data" is
  never confused with "we filtered it out".

hERG is out of this step on purpose (decision of 2026-11-11): it lives in assay_type=F against
target CHEMBL240 and is off-target safety, not ADME.

Pure function + CLI: python -m protein_mcp.files.admet P15056 --structure 3OG7 --pocket "032 C"
Imports common, files.choice, files.dock (the same manifest resolution as docking, so the ADMET
record and the docking run cannot disagree) and sources.chembl (+ rdkit).
"""
import argparse
import asyncio
import json
import sys
from pathlib import Path

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice, dock
from protein_mcp.sources import chembl

# ---- rules, pre-registered in ESTADO 7.3 on 2026-11-11 -------------------------------------
# Common-use conventions. The names are LABELS, not citations: no publication was retrieved in
# this project, so the cut-offs are the ones written here and nothing is attributed to "the
# literature". Breaking them is not a prediction of failure: vemurafenib breaks Lipinski by its
# logP and was approved in 2011, which is a fact of the chembl_record block of this same output.
#
# ONE SOURCE for each cut-off. The numbers that are APPLIED and the strings that are PUBLISHED
# are generated from this single table, so they cannot drift apart. They were two separate
# literals until the reviewer pointed out (2026-11-11) that the pre-registration test only
# checked the published strings: the applied number could then be changed silently while the
# report kept printing the pre-registered one, in a project whose whole argument is that the
# criterion was fixed in advance.
RULE_LIMITS = {
    "Lipinski": [("molecular_weight", 500), ("clogp", 5), ("hbd", 5), ("hba", 10)],
    "Veber": [("rotatable_bonds", 10), ("tpsa", 140)],
    "Egan": [("tpsa", 131.6), ("clogp", 5.88)],
}
_REPORTED_AS = {
    "Lipinski": "number of violations, not pass/fail: the rule itself allows one",
    "Veber": "number of violations",
    "Egan": "number of violations",
}
# A value is a violation only when it is STRICTLY above the limit: a molecule exactly at the
# cut-off satisfies it ("<= 500"), which is what the pre-registered strings say.
RULES = {
    name: {"cutoffs": [f"{field} <= {limit}" for field, limit in items],
           "reported_as": _REPORTED_AS[name]}
    for name, items in RULE_LIMITS.items()
}
RULE_DISCLAIMER = (
    "These rules are conventions of common use, pre-registered in this project on 2026-11-11 "
    "before any code was written. They are NOT measurements, NOT predictions and NOT derived "
    "from any source retrieved by this project; the names are labels, not citations. They are "
    "property filters on one molecule: they do not predict bioavailability, absorption or "
    "toxicity, and an approved drug can break them. The researcher's own judgement decides.")
CALCULATED_NOTE = (
    "Two engines calculate some of the same quantities: RDKit here, from the exact SDF that was "
    "docked, and ChEMBL's own engine on its entry. They are reported side by side, never averaged "
    "and never presented as confirming each other: two calculations agreeing is not a "
    "measurement, and if they disagree neither can be called correct from here.")
RECORD_DEFAULTS_NOTE = (
    "Two kinds of value in this block are NOT findings, and ChEMBL's API does not distinguish "
    "them from findings. `fields_minus_one` lists the fields that came back as -1: this project "
    "has not read ChEMBL's schema, so it does not say what -1 means and never reads it as a "
    "'no'. `boolean_fields_false` lists the flags that came back False: on an entry with nothing "
    "annotated they are all False, so a False here must not be reported as 'this compound is not "
    "oral' or 'was not withdrawn'. Both lists keep the raw values; neither is reinterpreted.")
MEASURED_NOTE = (
    "Retrieved measurements. Each row keeps the units, relation and assay description ChEMBL "
    "returned; nothing is averaged, aggregated or converted, because one standard_type arrives in "
    "incompatible units (CL as mL.min-1.kg-1 and as 'ml kg-1' in the same sample). A row with a "
    "relation ('>' or '<') is a censored value and the relation is part of the result.")
# F2, from the Desktop test of 2026-11-12. With the ChEMBL lookup down (HTTP 500) and the
# manifest pocket name None -- so no retrieved source held any name at all -- the compound was
# reported as "the chemically linked vemurafenib inhibitor (VEM-3-VEM)". Both the code and the
# description came from nowhere. The -1 and False fields each carry their own warning, so the
# absence of a name now carries one too, in the block where the absence shows up.
NAME_RULE = (
    "Refer to this compound by its PDB chemical component id and its InChIKey. ChEMBL gives no "
    "`pref_name` for it here, so THIS OUTPUT holds no name: do not invent one, do not give it a "
    "code of your own, and do not tie it to a named drug -- 'an analogue of X', 'a linked X', "
    "'X dimer' -- because nothing retrieved here says that. One other name does exist outside "
    "this output: wwPDB's systematic chemical name of the component, retrieved by get_structures "
    "and kept in the manifest (the research report prints it). It is a systematic name, not a "
    "compound or trade name, and this output deliberately does not carry it; if a structural "
    "claim is made from it, say it comes from that systematic name. Size, formula and "
    "descriptors are in `calculated_properties` and can be reported as what they are.")


class AdmetError(Exception):
    def __init__(self, message: str, **extra):
        super().__init__(message)
        self.extra = extra


# ---- local calculation ---------------------------------------------------------------------
def descriptors(mol) -> dict:
    """RDKit descriptors, each labelled with the exact function that produced it.

    The function names are part of the output on purpose: 'hbd' means different things in
    different tools, so the only auditable statement is which call returned the number.
    """
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

    values = {
        "molecular_weight": (Descriptors.MolWt(mol), "rdkit.Chem.Descriptors.MolWt"),
        "clogp": (Crippen.MolLogP(mol), "rdkit.Chem.Crippen.MolLogP"),
        "tpsa": (rdMolDescriptors.CalcTPSA(mol), "rdkit.Chem.rdMolDescriptors.CalcTPSA"),
        "hbd": (Lipinski.NumHDonors(mol), "rdkit.Chem.Lipinski.NumHDonors"),
        "hba": (Lipinski.NumHAcceptors(mol), "rdkit.Chem.Lipinski.NumHAcceptors"),
        "rotatable_bonds": (Lipinski.NumRotatableBonds(mol),
                            "rdkit.Chem.Lipinski.NumRotatableBonds"),
        "aromatic_rings": (Lipinski.NumAromaticRings(mol), "rdkit.Chem.Lipinski.NumAromaticRings"),
        "rings": (rdMolDescriptors.CalcNumRings(mol), "rdkit.Chem.rdMolDescriptors.CalcNumRings"),
        "fraction_csp3": (rdMolDescriptors.CalcFractionCSP3(mol),
                          "rdkit.Chem.rdMolDescriptors.CalcFractionCSP3"),
        "heavy_atoms": (mol.GetNumHeavyAtoms(), "rdkit Mol.GetNumHeavyAtoms"),
        "formal_charge": (Chem.GetFormalCharge(mol), "rdkit.Chem.GetFormalCharge"),
        "qed": (QED.qed(mol), "rdkit.Chem.QED.qed"),
    }
    out = {k: (round(v, 3) if isinstance(v, float) else v) for k, (v, _) in values.items()}
    out["molecular_formula"] = rdMolDescriptors.CalcMolFormula(mol)
    out["engine"] = f"RDKit {Chem.rdBase.rdkitVersion}"
    out["functions"] = {k: fn for k, (_, fn) in values.items()}
    return out


def rule_reading(desc: dict) -> dict:
    """The three rules applied to the RDKit descriptors. Conventions, offered, never measured."""
    out = {}
    for name, items in RULE_LIMITS.items():     # the same table the published strings come from
        broken, unknown = [], []
        for field, limit in items:
            value = desc.get(field)
            if value is None:
                unknown.append(field)
            elif value > limit:
                broken.append(f"{field} {value} > {limit}")
        out[name] = {
            "violations": None if unknown else len(broken),
            "violated": broken,
            "not_evaluated": unknown or None,
            "cutoffs": RULES[name]["cutoffs"],
        }
    return {
        "status": "READING OFFERED, NOT MEASURED",
        "computed_on": "the RDKit descriptors of this output, not on ChEMBL's",
        "rules_used": RULES,
        "results": out,
        "disclaimer": RULE_DISCLAIMER,
    }


# ---- which molecule ------------------------------------------------------------------------
def _inside(prep: Path, rel) -> Path:
    """The two causes are reported apart: a missing file and a path escaping prep/ are not the
    same event, and one message for both made the path guard untestable on its own."""
    path = (prep / str(rel)).resolve()
    if prep.resolve() not in path.parents:
        raise AdmetError(f"{rel} resolves outside prep/ and was refused.")
    if not path.is_file():
        raise AdmetError(f"{rel} is inside prep/ but the file is missing.")
    return path


def _first_mol(path: Path):
    """First record of an SDF, read through a file object (non-ASCII paths; ESTADO 8.3)."""
    from rdkit import Chem
    with open(path, "rb") as fh:
        for mol in Chem.ForwardSDMolSupplier(fh, removeHs=False):
            if mol is not None:
                return mol
    raise AdmetError(f"RDKit could not read a molecule from {path.name}.")


def _inchikey(mol) -> str:
    from rdkit.Chem import inchi
    key = inchi.MolToInchiKey(mol)
    if not key:
        raise AdmetError("RDKit returned no InChIKey for this molecule, so ChEMBL was not queried.")
    return key


def _pick_run(s: dict, pocket_key: str, conformer: str):
    """(run_id, run) for this pocket, or (None, None) when no docking run exists.

    Not imported from md_package on purpose: there, no run is a hard error (nothing can be
    delivered to CHARMM-GUI without a pose), while here the crystallographic ligand is a valid
    fallback, because ADMET describes a compound and not a pose.
    """
    runs = ((s.get("docking") or {}).get(pocket_key) or {}).get("runs") or {}
    if not runs:
        return None, None
    if conformer:
        run_id = dock.run_id_of(pocket_key.split()[0], conformer.strip().lower())
        if run_id not in runs:
            raise AdmetError(f"No run '{run_id}'. Available: {', '.join(sorted(runs))}.")
    elif len(runs) == 1:
        run_id = next(iter(runs))
    else:
        raise AdmetError(f"{len(runs)} docking runs for this pocket; name one with conformer "
                         f"({', '.join(sorted(r.rsplit('_', 1)[-1] for r in runs))}).")
    return run_id, runs[run_id]


def _crystal_sdf(s: dict, comp: str, label: str):
    """The crystallographic SDF of this instance, from nonpolymer_instances (the pocket block
    does not carry it)."""
    for ni in s.get("nonpolymer_instances") or []:
        if ni.get("comp_id") == comp and ni.get("label_asym_id") == label:
            sdf = ni.get("sdf") or {}
            return (sdf.get("path") if sdf.get("ok") else None), sdf.get("sha256")
    return None, None


def identity(prep: Path, s: dict, snap: dict, pocket_key: str, conformer: str, pose: int):
    """(identity record, the RDKit molecule) for the compound ADMET describes.

    The docked pose is preferred: it is the molecule that went through docking and MD. The
    crystallographic SDF is the fallback when no run exists, and is marked as such. When BOTH
    exist their InChIKeys are compared -- equal for both ligands of this project, measured on
    2026-11-11 -- and a mismatch means the run and the pocket have drifted apart, so it is
    reported and the ChEMBL lookup does not happen.
    """
    from rdkit import Chem

    comp, label = snap["pocket"]["comp_id"], snap["pocket"]["label_asym_id"]
    out = {"comp_id": comp, "pdb_id": snap["structure"], "pocket": pocket_key}
    keys, mols, files = {}, {}, {}

    run_id, run = _pick_run(s, pocket_key, conformer)
    if run:
        rel = next((v["path"] for k, v in (run.get("files") or {}).items()
                    if k.endswith(".sdf")), None)
        if not rel:
            raise AdmetError(f"Docking run {run_id} records no poses SDF.")
        path = _inside(prep, rel)
        recorded = (run["files"][Path(rel).name] or {}).get("sha256")
        actual = dock.sha256(path)
        if recorded and actual != recorded:
            raise AdmetError(f"{Path(rel).name} changed since the manifest (SHA-256): dock again "
                             "before recording ADMET.")
        with open(path, "rb") as fh:
            poses = list(Chem.ForwardSDMolSupplier(fh, removeHs=False))
        if not 1 <= pose <= len(poses):
            raise AdmetError(f"pose {pose} is outside 1..{len(poses)}.")
        if poses[pose - 1] is None:
            raise AdmetError(f"RDKit could not read pose {pose} of {Path(rel).name}.")
        mols["docked_pose"] = poses[pose - 1]
        keys["docked_pose"] = _inchikey(poses[pose - 1])
        files["docked_pose"] = (rel, actual)

    crystal_rel, _ = _crystal_sdf(s, comp, label)
    if crystal_rel:
        try:
            path = _inside(prep, crystal_rel)
            mol = _first_mol(path)
            mols["crystallographic"] = mol
            keys["crystallographic"] = _inchikey(mol)
            files["crystallographic"] = (crystal_rel, dock.sha256(path))
        except AdmetError as exc:
            out["crystallographic_not_read"] = str(exc)

    which = "docked_pose" if "docked_pose" in mols else (
        "crystallographic" if "crystallographic" in mols else None)
    if which is None:
        # the reason the crystallographic file was rejected must survive: without it the answer
        # is "no SDF" for a file that exists and was refused, which hides a refusal
        why = out.get("crystallographic_not_read")
        raise AdmetError(f"No readable SDF for {comp} {label} in {snap['structure']}: neither a "
                         "docked pose nor the crystallographic ligand file"
                         + (f". The crystallographic file was refused: {why}" if why else "."))

    out["source"] = ("docked pose" if which == "docked_pose" else
                     "crystallographic SDF (no docking run for this pocket)")
    out["source_file"], out["source_sha256"] = files[which]
    out["run_id"] = run_id if which == "docked_pose" else None
    out["pose"] = pose if which == "docked_pose" else None
    out["inchikey"] = keys[which]
    out["inchikey_by_file"] = keys

    if len(set(keys.values())) > 1:
        out["identity_check"] = "MISMATCH"
        raise AdmetError(
            "The docked pose and the crystallographic ligand have different InChIKeys "
            f"({keys}), so they are not the same compound. ChEMBL was not queried: re-run "
            "dock_ligand for this pocket before recording ADMET.", compound_identity=out)
    # "nothing to compare" must not cover the case where a second SDF EXISTED and was refused:
    # the comparison is this step's own guarantee, and saying "not applicable" would hide that
    # it was not performed (reviewer finding, 2026-11-11)
    if len(keys) > 1:
        out["identity_check"] = "the docked pose and the crystallographic SDF give the same InChIKey"
    elif out.get("crystallographic_not_read"):
        out["identity_check"] = ("NOT COMPARED: the crystallographic SDF exists in the manifest but "
                                 f"could not be read ({out['crystallographic_not_read']}), so the "
                                 "pose was not checked against a second file")
    else:
        out["identity_check"] = "only one SDF was available, so there was nothing to compare"
    out["molecule_source_note"] = (
        "ADMET describes a compound, not a pose, so the geometry is not used here: the SDF is read "
        "only for its atoms and bonds. The file and its SHA-256 are named so the compound can be "
        "traced to the run that used it.")
    return out, mols[which]


# ---- the tool ------------------------------------------------------------------------------
async def analyze_admet(uniprot_id: str, structure: str = "", pocket: str = "",
                        conformer: str = "", pose: int = 1) -> dict:
    """ADME records, measurements and calculated properties. See the module docstring."""
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    path = choice.manifest_path(accession)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"ok": False, "uniprot_id": accession,
                "error": f"No readable manifest at {path}: run prepare_structures first."}
    prep = path.parent
    try:
        s, snap, origin = dock._resolve(manifest, structure, pocket)
        pocket_key = f"{snap['pocket']['comp_id']} {snap['pocket']['label_asym_id']}"
        ident, mol = identity(prep, s, snap, pocket_key, conformer, pose)
    except (dock.DockError, AdmetError) as exc:
        out = {"ok": False, "uniprot_id": accession, "error": str(exc)}
        out.update(getattr(exc, "extra", {}) or {})
        return out

    desc = descriptors(mol)
    # identity check, not a property check: the pocket block records the heavy atoms of the
    # crystallographic instance, so a disagreement means the SDF is not that ligand
    in_manifest = next((p.get("heavy_atoms") for p in
                        ((s.get("prepared") or {}).get("pockets") or [])
                        if p.get("comp_id") == snap["pocket"]["comp_id"]
                        and p.get("label_asym_id") == snap["pocket"]["label_asym_id"]), None)
    ident["heavy_atoms_in_manifest_pocket"] = in_manifest
    ident["heavy_atoms_agree"] = (None if in_manifest is None
                                  else in_manifest == desc["heavy_atoms"])
    if ident["heavy_atoms_agree"] is False:
        log.warning("heavy atoms differ for %s: SDF %s, manifest pocket %s",
                    pocket_key, desc["heavy_atoms"], in_manifest)

    # the molecule lookup failing is UNKNOWN too, and it used to escape as a RuntimeError
    # through the MCP boundary while ESTADO claimed the block went to UNKNOWN (reviewer finding)
    # one budget for every ChEMBL request of this call: a slow source must not outlive the
    # caller, because a tool that never answers is worse than one that answers UNKNOWN
    # (measured in the Desktop test of 2026-11-12: two timeouts, nothing shown)
    deadline = chembl.Deadline(chembl.ADME_BUDGET)
    try:
        found = await chembl.fetch_molecule_by_inchikey(ident["inchikey"], deadline=deadline)
    except Exception as exc:
        log.warning("ChEMBL molecule lookup failed: %s", exc)
        found = {"found": False, "lookup_failed": True,
                 "inchi_key": ident["inchikey"], "error": f"{type(exc).__name__}: {exc}",
                 "note": "The ChEMBL molecule lookup FAILED, so whether this compound has a "
                         "ChEMBL entry is UNKNOWN. This is not 'no record'."}
    if found.get("lookup_failed"):
        measured = {"retrieved": False, "lookup_reached_chembl": False,
                    "note": "No ADME query was made: the molecule lookup failed first, so "
                            "whether measurements exist is UNKNOWN. This is not 'no data'."}
    elif not found.get("found"):
        measured = {"retrieved": False, "lookup_reached_chembl": True, "no_chembl_entry": True,
                    "note": "No ADME query was made because no ChEMBL molecule has this exact "
                            "InChIKey. ChEMBL holds no measurements for a compound it does not "
                            "have; a measurement published elsewhere would not be here."}
    else:
        try:
            measured = await chembl.fetch_adme_activities(found["chembl_id"], deadline=deadline)
        except Exception as exc:  # the lookup failing is not "no data"
            log.warning("ADME activities failed: %s", exc)
            measured = {"retrieved": False, "lookup_reached_chembl": True,
                        "error": f"{type(exc).__name__}: {exc}",
                        "note": "The ADME activity lookup FAILED, so whether measurements exist "
                                "is UNKNOWN. This is not 'no data'."}

    result = {
        "ok": True,
        "uniprot_id": accession,
        "generated_by": "analyze_admet",
        "generated_at": utc_now(),
        "structure_origin": origin,
        "compound_identity": ident,
        "chembl_lookup": {k: v for k, v in found.items()
                          if k not in ("record", "chembl_calculated_properties",
                                       "molecule_hierarchy")},
        "chembl_record": ({"status": "RECORD (facts of the ChEMBL entry)",
                           "chembl_id": found.get("chembl_id"),
                           "link": found.get("link"),
                           "pref_name": found.get("pref_name"),
                           "name_status": _name_status(found),
                           "fields": found.get("record"),
                           "fields_minus_one": found.get("record_fields_minus_one"),
                           "boolean_fields_false": found.get("record_boolean_fields_false"),
                           "defaults_note": RECORD_DEFAULTS_NOTE,
                           "atc_classifications": found.get("atc_classifications"),
                           "cross_references": found.get("cross_references"),
                           "molecule_hierarchy": found.get("molecule_hierarchy")}
                          if found.get("found") else
                          {"status": ("UNKNOWN (the ChEMBL lookup failed)"
                                      if found.get("lookup_failed") else "NO RECORD"),
                           "name_status": _name_status(found),
                           "note": found.get("note") or found.get("error")}),
        # the status must not say MEASURED when nothing was retrieved (reviewer finding): the
        # status is the field a summary reads first, and it contradicted its own note
        "measured_adme": {"status": ("MEASURED (retrieved from ChEMBL)" if measured.get("retrieved")
                                     else "UNKNOWN (nothing was retrieved)"),
                          "note": MEASURED_NOTE if measured.get("retrieved") else None,
                          **measured},
        "calculated_properties": {
            "status": "CALCULATED (not measured, not predicted)",
            "note": CALCULATED_NOTE,
            "rdkit": desc,
            "chembl": (found.get("chembl_calculated_properties") if found.get("found") else None),
            "chembl_engine": ("ChEMBL's own calculation on its entry, engine not reported by the API"
                              if found.get("found") else None),
        },
        "rule_reading": rule_reading(desc),
        "selection": _selection(ident, found, measured, desc),
        "notes": [
            "No value in this output is PREDICTED: nothing here comes from a trained model "
            "(Phase 3 step 7.1).",
            *([] if found.get("pref_name") else
              [f"This output holds NO name for this compound. {NAME_RULE}"]),
            "hERG and other off-target safety data are NOT in this output: in ChEMBL they are "
            "assay_type=F against a target, not ADME (decision of 2026-11-11).",
        ],
    }
    warn = _record(path, manifest, snap["structure"], pocket_key, result, prep)
    if warn:
        result["manifest_warning"] = warn
    else:
        result["recorded_in"] = str(path)
    return result


def _name_status(found: dict) -> str:
    """Whether this compound has a name from a retrieved source, said out loud. See NAME_RULE."""
    if found.get("lookup_failed"):
        return ("UNKNOWN, and NO NAME IS AVAILABLE HERE: the ChEMBL lookup did not complete, so "
                "whether ChEMBL names this compound is UNKNOWN and this output holds no name. "
                + NAME_RULE)
    if not found.get("found"):
        return ("NO NAME IN CHEMBL: no ChEMBL molecule has this exact InChIKey, so ChEMBL holds "
                "no name for it. " + NAME_RULE)
    if not found.get("pref_name"):
        return (f"UNNAMED IN CHEMBL: entry {found.get('chembl_id')} exists but its `pref_name` is "
                "None, which is the normal state of a research compound. " + NAME_RULE)
    return (f"NAMED BY THE SOURCE: `pref_name` of {found.get('chembl_id')} is "
            f"\"{found['pref_name']}\". Use that name. Any other name must come from a field of "
            "this output -- `cross_references` can carry trade names -- and must be attributed to "
            "it; this output returns no synonym list of its own.")


def _selection(ident, found, measured, desc) -> str:
    """Always say how the compound was selected and how much was found ('N of total')."""
    who = f"{ident['comp_id']} of {ident['pdb_id']} pocket {ident['pocket']}, read from the " \
          f"{ident['source']}"
    tail = (f" RDKit: MW {desc['molecular_weight']}, cLogP {desc['clogp']}, "
            f"TPSA {desc['tpsa']}.")
    if found.get("lookup_failed"):
        return (f"{who}. The ChEMBL molecule lookup FAILED, so whether InChIKey "
                f"{ident['inchikey']} has an entry, a record or any measurement is UNKNOWN; only "
                f"the RDKit calculation and the rule reading apply.{tail}")
    if not found.get("found"):
        return (f"{who}. InChIKey {ident['inchikey']} has no ChEMBL entry, so there are no records "
                f"and no measurements; only the RDKit calculation and the rule reading apply.{tail}")
    rows = measured.get("adme_rows")
    if rows is None:
        return (f"{who}. ChEMBL entry {found['chembl_id']}; whether ADME measurements exist is "
                f"UNKNOWN (the lookup failed).{tail}")
    read = measured.get("rows_read") or 0
    held = measured.get("rows_total")
    classes = sorted({r["adme_class"] for r in rows})
    calculated = sum(1 for r in rows if r.get("description_says_value_was_calculated"))
    return (f"{who}. InChIKey {ident['inchikey']} matches exactly 1 ChEMBL molecule "
            f"({found['chembl_id']}). {len(rows)} of {read} rows read are on this project's ADME "
            "allowlist"
            + (f" ({', '.join(classes)})" if classes else "")
            + (f"; the other {read - len(rows)} are not on the allowlist, which is narrower than "
               "ADME, so they are not 'not ADME'" if read - len(rows) else "")
            # the cap is a real one: CHEMBL53463 has 1143 rows of assay_type=A
            + (f". TRUNCATED: ChEMBL holds {held} rows in the assay types queried and only {read} "
               f"were read, so {held - read} were never examined and the count above is a LOWER "
               "BOUND" if measured.get("truncated") else
               f". All {held} rows ChEMBL holds in the assay types queried were read"
               if held is not None else "")
            + (f". {calculated} of the allowlisted rows say in their own assay description that "
               "the value was CALCULATED by the depositor, not measured" if calculated else "")
            + "." + tail)


def _record(path: Path, manifest: dict, pdb_id: str, pocket_key: str, result: dict, prep: Path):
    """Write an `admet` block next to the docking block. Re-reads the manifest first, as
    build_md_package does: if prepare_structures rebuilt it meanwhile, say so."""
    try:
        current = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return f"The result could not be recorded: the manifest could not be re-read ({exc!r})."
    if current.get("created_at") != manifest.get("created_at"):
        return ("prepare_structures rewrote the manifest while this ran, so the result was not "
                "recorded there. Run analyze_admet again.")
    try:
        s = next(x for x in current["structures"] if x.get("pdb_id") == pdb_id)
        block = s.setdefault("admet", {})
        block[pocket_key] = {k: result[k] for k in
                             ("generated_by", "generated_at", "compound_identity", "chembl_lookup",
                              "chembl_record", "measured_adme", "calculated_properties",
                              "rule_reading", "selection", "notes")}
        choice.write_manifest(path, current)
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        log.exception("Could not write %s", path)
        return f"The result could not be recorded in the manifest ({exc!r})."
    return None


def compact(result: dict) -> dict:
    """Chat view: drops the per-row assay descriptions' bulk and the RDKit function map."""
    if not result.get("ok"):
        return result
    out = json.loads(json.dumps(result))
    out["calculated_properties"]["rdkit"].pop("functions", None)
    rows = (out.get("measured_adme") or {}).get("adme_rows")
    if isinstance(rows, list) and len(rows) > 20:
        out["measured_adme"]["adme_rows"] = rows[:20]
        out["measured_adme"]["adme_rows_shown"] = f"20 of {len(rows)}"
    # the excluded rows carry 300-character descriptions and were NOT trimmed, so in the chat
    # view they could outnumber the measured ones two to one -- which pushes towards exactly the
    # reading `rows_excluded_note` forbids (reviewer finding 2026-11-12)
    excl = (out.get("measured_adme") or {}).get("rows_excluded_by_allowlist")
    if isinstance(excl, list) and len(excl) > 10:
        kept = excl[:10]
        out["measured_adme"]["rows_excluded_by_allowlist"] = kept
        out["measured_adme"]["rows_excluded_shown"] = (
            f"{len(kept)} of {sum((result.get('measured_adme') or {}).get(
                'rows_excluded_by_allowlist_by_standard_type', {}).values())} (chat view)")
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m protein_mcp.files.admet",
                               description="ADME records, measurements and calculated properties "
                                           "for the ligand of one pocket.")
    p.add_argument("uniprot_id")
    p.add_argument("--structure", default="")
    p.add_argument("--pocket", default="")
    p.add_argument("--conformer", default="")
    p.add_argument("--pose", type=int, default=1)
    a = p.parse_args(argv)
    res = asyncio.run(analyze_admet(a.uniprot_id, a.structure, a.pocket, a.conformer, a.pose))
    json.dump(res, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if res.get("ok") else 1


if __name__ == "__main__":
    raise SystemExit(main())
