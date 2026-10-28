"""Phase 2 step 2: cleaned receptor, structure metrics, pocket and box, computed with gemmi from
the ORIGINAL mmCIF downloaded by prepare_structures. Pure computation, no network.
Imports only common (and gemmi).

Validated on the real BRAF V600E files (2026-10-22):
- 3OG7: the 22 residues within 5 Å of vemurafenib (032, chain A) include K483, T529,
  Q530-C532 and D594-F595-G596; chain A has residues 597-614 (activation segment, V600)
  NOT modelled, chain B has GLU 600.
- 7ZR0: "residue 600" of chains A/B is HSP90 (ALA); the BRAF chain (label D, auth K) has 186
  of 828 residues modelled. So numbering is only read inside the target chains, and UniProt
  positions come from label_seq_id + RCSB's alignment (label 583 -> UniProt 521 = auth 521).
- gemmi's remove_ligands_and_waters also removes nucleotides/cofactors (e.g. KRAS GDP): the
  receptor is selected by hand. 5JSM has hydrogens and alternate conformations. 7ZR0 has no
  R-free ("?"). AlphaFold's B-factor column equals _ma_qa_metric_local (pLDDT), max diff 0.0.
"""
import re
from pathlib import Path

import gemmi

from protein_mcp.common import log

# >>> EDIT THESE LINES: pocket = target residues with any atom within POCKET_CUTOFF of the ligand;
# >>> box = ligand bounding box + BOX_PADDING on every side. <<<
POCKET_CUTOFF = 5.0
BOX_PADDING = 5.0
MIN_POCKET_HEAVY_ATOMS = 10  # smaller instances (EOH 3, BEN 9 in 5JSM: additives) get no pocket
LARGE_BOX_SIDE = 30.0        # a box side above this is flagged (docking boxes are usually smaller)
AF_LOW_PLDDT = 70.0          # AlphaFold "confident" threshold
KEEP_IN_RECEPTOR = {"nucleotide", "cofactor"}   # kept when they touch a target chain
# Metal ions usually functional (KRAS: Mg2+ with GDP, 8TVK); kept when they touch a target chain.
# NA, K, CL, CA, CD... are left out as likely crystallisation components.
KEEP_METALS = {"MG", "MN", "ZN", "FE", "CO", "NI", "CU"}
_SUBSTITUTION_RE = re.compile(r"([A-Z])(\d{1,5})([A-Z])")


def _one_letter(resname: str):
    info = gemmi.find_tabulated_residue(resname)
    code = info.one_letter_code.upper() if info else ""
    return code if code.isalpha() else None


def _uniprot_pos(label_seq, alignment):
    """label_seq_id -> UniProt position through RCSB's aligned regions; None outside them (tags, linkers)."""
    if label_seq is None:
        return None
    for reg in alignment or []:
        beg, ubeg, length = reg.get("entity_begin"), reg.get("uniprot_begin"), reg.get("length")
        if isinstance(beg, int) and isinstance(ubeg, int) and isinstance(length, int) and beg <= label_seq < beg + length:
            return ubeg + (label_seq - beg)
    return None


def _ranges(nums):
    out, start, prev = [], None, None
    for n in sorted(nums):
        if start is None:
            start = prev = n
        elif n == prev + 1:
            prev = n
        else:
            out.append([start, prev])
            start = prev = n
    if start is not None:
        out.append([start, prev])
    return out


def _r_free(path: Path):
    try:
        value = gemmi.cif.read(str(path)).sole_block().find_value("_refine.ls_R_factor_R_free")
        return round(float(value), 4) if value not in (None, "?", ".") else None
    except (RuntimeError, ValueError):
        return None


def parse_variant(mutation_filter, canonical: str):
    """'V600E' -> (wild-type letter FROM UNIPROT, position, mutant letter, warning or None).
    The wild-type letter is read from the canonical sequence, not from the user's text: with
    'E600V' a GLU 600 would otherwise be reported as wild type. None if not a substitution."""
    mo = _SUBSTITUTION_RE.fullmatch(str(mutation_filter or "").upper())
    if not mo:
        return None
    typed_wt, pos, mut = mo.group(1), int(mo.group(2)), mo.group(3)
    wt = canonical[pos - 1] if canonical and 0 < pos <= len(canonical) else None
    warning = None
    if wt is None:
        warning = "canonical UniProt sequence not available: the wild-type letter is unchecked"
        wt = typed_wt
    elif wt != typed_wt:
        warning = f"the filter says {typed_wt}{pos} but UniProt has {wt} at {pos}"
    return wt, pos, mut, warning


def _chain_metrics(model, chain, alignment, sample_length, variant):
    """Modelled residues, missing UniProt ranges inside the modelled span, and the residue at
    the requested variant position, for ONE target chain (label_asym_id)."""
    residues = [r for ch in model for r in ch if r.subchain == chain["label_asym_id"]]
    upos, shared = {}, {}
    for r in residues:
        u = _uniprot_pos(r.label_seq, alignment)
        if u is not None:
            if u in upos:                  # microheterogeneity: two residues at one position
                shared.setdefault(u, [upos[u].name]).append(r.name)
            else:
                upos[u] = r
    expected = set()
    for reg in alignment or []:
        if isinstance(reg.get("uniprot_begin"), int) and isinstance(reg.get("length"), int):
            expected.update(range(reg["uniprot_begin"], reg["uniprot_begin"] + reg["length"]))
    modelled = sorted(upos)
    span = set(range(modelled[0], modelled[-1] + 1)) if modelled else set()
    out = {
        "label_asym_id": chain["label_asym_id"],
        "auth_asym_id": chain["auth_asym_id"],
        "residues_modelled": len(residues),
        "sample_sequence_length": sample_length,
        "modelled_uniprot_range": [modelled[0], modelled[-1]] if modelled else None,
        "missing_uniprot_ranges": _ranges((expected & span) - set(modelled)),
    }
    if shared:
        out["microheterogeneity_uniprot"] = {str(u): names for u, names in sorted(shared.items())}
    if variant:
        wt, pos, mut, warning = variant
        r = upos.get(pos)
        observed = _one_letter(r.name) if r else None
        check = {
            "uniprot_position": pos,
            "uniprot_wild_type": wt,
            "observed": f"{r.name} {r.seqid.num}{r.seqid.icode.strip()}" if r else None,
            # wild type is tested first: with a wrong filter ("D13G", UniProt G13) a GLY 13 is wild type
            "status": ("not modelled" if r is None else "ambiguous (microheterogeneity)" if pos in shared
                       else "wild type" if observed == wt else "variant" if observed == mut else "other residue"),
        }
        if warning:
            check["warning"] = warning
        out["variant_check"] = check
    return out, upos


def _first_conformer(residue):
    """Ligand atoms without altloc plus those of its first altloc (the receptor keeps the first too)."""
    alts = sorted({a.altloc for a in residue if a.altloc != "\0"})
    return [a for a in residue if a.altloc == "\0" or (alts and a.altloc == alts[0])]


def _pocket(model, ligand_atoms, target_subchains, uniprot_of, ns):
    """Residues within POCKET_CUTOFF of any ligand atom: target chains vs other polymer chains.
    Asymmetric unit only: NeighborSearch also returns crystal-symmetry images (image_idx != 0),
    which put residues ~47 Å away in the pocket (8C7Y TXV, found by the reviewer)."""
    near = {}
    for atom in ligand_atoms:
        for mark in ns.find_atoms(atom.pos, "\0", radius=POCKET_CUTOFF):
            if mark.image_idx != 0:
                continue
            cra = mark.to_cra(model)
            res = cra.residue
            if res.entity_type != gemmi.EntityType.Polymer:
                continue
            key = (res.subchain, res.seqid.num, res.seqid.icode, res.name)
            d = cra.atom.pos.dist(atom.pos)
            if key not in near or d < near[key][0]:
                near[key] = (d, cra.chain.name, res)
    target, other = [], []
    for (sub, num, icode, name), (d, auth, res) in sorted(near.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        row = {"auth_asym_id": auth, "resname": name, "auth_seq_id": f"{num}{icode.strip()}",
               "min_distance": round(d, 2)}
        if sub in target_subchains:
            row["uniprot"] = uniprot_of.get(sub, {}).get(res.label_seq)
            target.append(row)
        else:
            other.append(row)
    return target, other


def _box(atoms):
    xs, ys, zs = zip(*((a.pos.x, a.pos.y, a.pos.z) for a in atoms))
    lo, hi = (min(xs), min(ys), min(zs)), (max(xs), max(ys), max(zs))
    return {
        "center": [round((a + b) / 2, 3) for a, b in zip(lo, hi)],
        "size": [round(b - a + 2 * BOX_PADDING, 2) for a, b in zip(lo, hi)],
        "ligand_extent": [round(b - a, 2) for a, b in zip(lo, hi)],
        "padding": BOX_PADDING,
        "frame": "coordinates of the original mmCIF (unchanged in the cleaned receptor)",
    }


def _write_receptor(st, out_dir: Path, pdb_id: str):
    """PDB when the format can hold it (1-character chain names, <= 99999 atoms), else mmCIF."""
    atoms = sum(len(r) for ch in st[0] for r in ch)
    if all(len(ch.name) == 1 for ch in st[0]) and atoms <= 99999:
        path = out_dir / f"{pdb_id}_receptor.pdb"
        st.write_pdb(str(path))
        return path, "pdb", atoms
    path = out_dir / f"{pdb_id}_receptor.cif"
    st.make_mmcif_document().write_file(str(path))
    return path, "mmcif", atoms


def clean_structure(cif: Path, out_dir: Path, record: dict, mutation_filter, canonical: str = "") -> dict:
    """Cleaned receptor + metrics + pocket/box per ligand instance for ONE experimental entry.
    `record` is the manifest record of prepare_structures (target chains, instances);
    `canonical` is the UniProt canonical sequence (for the wild-type letter of the variant).
    Return {ok, receptor, metrics, pockets} or {ok: False, error}; never raises."""
    try:
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        n_models = len(st)
        if n_models == 0:
            return {"ok": False, "error": "cleaning failed: the file has no model (no atoms)"}
        while len(st) > 1:            # NMR / multi-model: keep the first model
            del st[len(st) - 1]
        model = st[0]
        variant = parse_variant(mutation_filter, canonical)

        chains, uniprot_of, target_subchains = [], {}, set()
        for ent in record.get("target_entities") or []:
            for chain in ent.get("chains") or []:
                metrics, upos = _chain_metrics(model, chain, ent.get("uniprot_alignment"),
                                               ent.get("sample_sequence_length"), variant)
                chains.append(metrics)
                target_subchains.add(chain["label_asym_id"])
                uniprot_of[chain["label_asym_id"]] = {r.label_seq: u for u, r in upos.items()}

        instances = {i["label_asym_id"]: i for i in record.get("nonpolymer_instances") or []}
        ns = gemmi.NeighborSearch(model, st.cell, POCKET_CUTOFF + 1).populate()
        pockets, skipped = [], []
        for label, inst in sorted(instances.items()):
            if inst["category"] != "ligand":
                continue                    # nucleotides/cofactors stay in the receptor (DECISIONS)
            if not inst.get("neighbor_target_chains"):
                skipped.append({"comp_id": inst["comp_id"], "label_asym_id": label,
                                "reason": "RCSB lists no target-chain neighbour"})
                continue
            residue = next((r for ch in model for r in ch if r.subchain == label and r.name == inst["comp_id"]), None)
            if residue is None:
                pockets.append({"comp_id": inst["comp_id"], "label_asym_id": label, "error": "not found in the file"})
                continue
            atoms = _first_conformer(residue)
            heavy = sum(1 for a in atoms if not a.is_hydrogen())
            if heavy < MIN_POCKET_HEAVY_ATOMS:
                skipped.append({"comp_id": inst["comp_id"], "label_asym_id": label, "heavy_atoms": heavy,
                                "reason": f"fewer than {MIN_POCKET_HEAVY_ATOMS} heavy atoms"})
                continue
            target, other = _pocket(model, atoms, target_subchains, uniprot_of, ns)
            box = _box(atoms)
            by_chain = {}
            for row in target:
                by_chain[row["auth_asym_id"]] = by_chain.get(row["auth_asym_id"], 0) + 1
            pockets.append({
                "comp_id": inst["comp_id"], "label_asym_id": label, "auth_asym_id": inst.get("auth_asym_id"),
                "heavy_atoms": heavy, "cutoff": POCKET_CUTOFF,
                # 5JSM 6NB (two vemurafenib units joined by a linker, 40 Å between farthest heavy atoms): one sits in
                # the ATP sites of chains A (24 residues) and B (23), so its box covers both sites
                "residues_by_chain": dict(sorted(by_chain.items())),
                "spans_target_chains": len(by_chain) > 1,
                "residues": target, "partner_residues": other,
                "box": box, "large_box": max(box["size"]) > LARGE_BOX_SIDE,
            })
        _flag_gaps(pockets, chains)

        # Receptor: target polymer chains + nucleotides/cofactors touching them; no waters, additives,
        # ligands, hydrogens; first alternate conformation only.
        keep_np = {lab for lab, i in instances.items()
                   if (i["category"] in KEEP_IN_RECEPTOR or i["comp_id"] in KEEP_METALS)
                   and i.get("neighbor_target_chains")}
        for ch in model:
            for idx in range(len(ch) - 1, -1, -1):
                r = ch[idx]
                if r.subchain in target_subchains and r.entity_type == gemmi.EntityType.Polymer:
                    continue
                if r.subchain in keep_np:
                    continue
                del ch[idx]
        st.remove_empty_chains()
        st.remove_hydrogens()
        st.remove_alternative_conformations()
        if not any(len(r) for ch in st[0] for r in ch):
            return {"ok": False, "error": "cleaning failed: no atom of the target chains "
                                          f"({', '.join(sorted(target_subchains)) or 'none'}) is in the file"}
        out_dir.mkdir(parents=True, exist_ok=True)
        path, fmt, atoms = _write_receptor(st, out_dir, record["pdb_id"])
    except (RuntimeError, ValueError, OSError, IndexError) as exc:
        log.exception("Cleaning failed for %s", cif)
        return {"ok": False, "error": f"cleaning failed: {exc!r}"}
    return {
        "ok": True,
        "receptor": {"path": str(path), "format": fmt, "atoms": atoms,
                     "chains_auth": [ch.name for ch in st[0]],
                     "kept_nonpolymer_label_asym_ids": sorted(keep_np)},
        "metrics": {"r_free": _r_free(cif), "models_in_file": n_models, "chains": chains},
        "pockets": pockets,
        "pocket_not_computed": skipped,
    }


def _flag_gaps(pockets, chains):
    """A missing range next to a pocket residue (flank at +/-1 in UniProt numbering) is flagged:
    the pocket may be incomplete there (3OG7 chain A: 597-614 missing, G596 in the pocket)."""
    missing = {c["auth_asym_id"]: c["missing_uniprot_ranges"] for c in chains}
    for p in pockets:
        flagged = set()
        for row in p.get("residues") or []:
            u = row.get("uniprot")
            for lo, hi in missing.get(row["auth_asym_id"]) or []:
                if u is not None and (u == lo - 1 or u == hi + 1):
                    flagged.add((row["auth_asym_id"], lo, hi))
        p["adjacent_missing_ranges"] = [{"auth_asym_id": c, "missing_uniprot_range": [lo, hi]}
                                        for c, lo, hi in sorted(flagged)]


def residue_map(cif: Path, record=None, auth_chain: str = "") -> dict:
    """UniProt number -> modelled residue ("GLU", auth_seq_id with insertion code) in ONE chain, for
    the chat control of PyMOL (step 4b). Experimental: the target chain with this auth id, numbers
    from label_seq_id + RCSB's alignment (as in _chain_metrics). record None = AlphaFold model:
    its residue numbers are UniProt numbers (as in alphafold_metrics). Positions with two residues
    (microheterogeneity) are listed apart. Never raises."""
    try:
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        model = st[0]
        if record is None:
            residues = {r.seqid.num: (r.name, str(r.seqid.num)) for r in model[0] if len(r)}
            return {"ok": True, "residues": residues, "ambiguous": {}}
        for ent in record.get("target_entities") or []:
            for chain in ent.get("chains") or []:
                if chain.get("auth_asym_id") != auth_chain:
                    continue
                residues, ambiguous = {}, {}
                for r in (r for ch in model for r in ch if r.subchain == chain["label_asym_id"]):
                    u = _uniprot_pos(r.label_seq, ent.get("uniprot_alignment"))
                    if u is None:
                        continue
                    if u in residues:
                        ambiguous.setdefault(u, [residues[u][0]]).append(r.name)
                    else:
                        residues[u] = (r.name, f"{r.seqid.num}{r.seqid.icode.strip()}")
                return {"ok": True, "residues": residues, "ambiguous": ambiguous}
        return {"ok": False, "error": f"chain {auth_chain} is not a target chain of {record.get('pdb_id')}"}
    except (RuntimeError, ValueError, OSError, IndexError) as exc:
        log.exception("Residue map failed for %s", cif)
        return {"ok": False, "error": f"could not read {cif.name}: {exc!r}"}


def alphafold_metrics(cif: Path, pockets_by_pdb: dict, variant_filter, canonical: str = "") -> dict:
    """pLDDT (B-factor column, equal to _ma_qa_metric_local) of the AlphaFold model at each
    experimental pocket (same UniProt numbers) and at the variant position. Only valid for the
    canonical-sequence model: the caller checks is_canonical_sequence."""
    try:
        st = gemmi.read_structure(str(cif))
        plddt = {r.seqid.num: (r.name, r[0].b_iso) for r in st[0][0] if len(r)}
    except (RuntimeError, ValueError, OSError, IndexError) as exc:
        log.exception("AlphaFold metrics failed for %s", cif)
        return {"ok": False, "error": f"AlphaFold metrics failed: {exc!r}"}
    out = {"ok": True, "pockets": []}
    for pdb_id, pockets in pockets_by_pdb.items():
        for p in pockets:
            nums = sorted({r["uniprot"] for r in p.get("residues") or [] if r.get("uniprot") is not None})
            values = [plddt[n][1] for n in nums if n in plddt]
            out["pockets"].append({
                "pdb_id": pdb_id, "comp_id": p["comp_id"], "label_asym_id": p["label_asym_id"],
                "residues_compared": len(values),
                "mean_plddt": round(sum(values) / len(values), 2) if values else None,
                "residues_below_70": [n for n in nums if n in plddt and plddt[n][1] < AF_LOW_PLDDT],
            })
    variant = parse_variant(variant_filter, canonical)
    if variant and variant[1] in plddt:
        wt, pos = variant[0], variant[1]
        name, value = plddt[pos]
        row = {"uniprot_position": pos, "residue": name, "plddt": round(value, 2)}
        if not canonical:
            row["note"] = "Not compared with UniProt: the canonical sequence is not available."
        elif _one_letter(name) == wt:
            row["note"] = "Wild-type residue (matches the UniProt canonical sequence)."
        else:
            row["note"] = f"Differs from the UniProt canonical residue ({wt} at {pos})."
        out["variant_position"] = row
    return out
