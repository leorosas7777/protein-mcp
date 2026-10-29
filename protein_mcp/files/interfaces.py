"""analyze_interfaces: partner chains, sequences, assemblies and protein-protein interfaces of the
prepared structures (Phase 2, step 5), written per structure as `interfaces` in prep/manifest.json.

Validated 2026-10-28 (ESTADO_PROYECTO.md 8.2) on 3WD5 (TNF + adalimumab Fab), 2W2N, 3H42, 6MV5, 6XIC:
- RCSB's biological assembly can be just the asymmetric unit (3WD5: 1 TNF + Fab, operator "1"),
  while the TNF trimer is formed by crystal symmetry; PDBe PISA calls the trimer stable and gives
  css 0 to every TNF-Fab interface. The sources disagree, so both are reported with their origin
  and NO interface is called biological here; no oligomer is built.
- Contacts are computed with gemmi, as the pockets: heavy atoms of polymer residues within
  CONTACT_CUTOFF, including crystal-symmetry copies, each contact labelled with its symmetry
  operator (gemmi code "5_555" = "z,x,y", the same operator PISA writes "Z,X,Y"). The target
  residues found this way are all in PISA's buried-residue list for 3WD5 (PISA adds 1-9 per
  interface buried farther than 5 Å). A Fab of 3WD5 touches two TNF protomers (own copy and the
  3-fold copy): visible only with symmetry.
- Partner chains are described, not classified: "no UniProt = antibody" fails (6MV5 light chain
  has one). Sequences are those of the crystallised construct (tags, single domains, "X").
Nothing is computed beyond contacts; coordinates are never changed. Imports common, dossier,
files.{choice, clean} and sources.{rcsb, pisa}.
"""
import asyncio
import hashlib
import json
import re

import gemmi
import httpx

from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files import choice, clean
from protein_mcp.sources import pisa, rcsb

CONTACT_CUTOFF = 5.0      # Å, heavy atoms, as POCKET_CUTOFF
PARALLEL_PISA = 3
_HIS_TAG = re.compile(r"H{6,}")


def sequence_flags(seq: str) -> list:
    """Facts read from a construct sequence that matter for a complex prediction (Phase 3). A run
    of histidines is often an affinity tag, but RCSB does not say so and native poly-His exists:
    it is reported as a run, not as a tag. Other tags (Strep, FLAG, fusions) are not detected."""
    flags = ["construct sequence as deposited (may differ from the native protein)"]
    if not seq:
        return flags + ["no sequence from RCSB"]
    n_x = seq.count("X")
    if n_x:
        flags.append(f"{n_x} non-standard or unknown residue{'s' if n_x != 1 else ''} (X)")
    for m in _HIS_TAG.finditer(seq):
        flags.append(f"{len(m.group(0))} consecutive histidines at residues {m.start() + 1}-{m.end()} "
                     "(often an affinity tag; not stated by RCSB)")
    return flags


def _ops(st):
    sg = st.find_spacegroup()
    return list(sg.operations()) if sg else []


def compute_contacts(cif, record: dict, cutoff: float = CONTACT_CUTOFF) -> dict:
    """Target-chain residues within `cutoff` of any other polymer chain or symmetry copy, grouped
    by (target chain, partner chain, symmetry operator). First model only. Never raises."""
    try:
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        while len(st) > 1:
            del st[len(st) - 1]
        model, cell = st[0], st.cell
        ops = _ops(st)
        targets, uniprot_aln, role, auth_of = set(), {}, {}, {}
        for ent in record.get("target_entities") or []:
            for c in ent.get("chains") or []:
                targets.add(c["label_asym_id"])
                uniprot_aln[c["label_asym_id"]] = ent.get("uniprot_alignment")
        for ent in (record.get("target_entities") or []) + (record.get("other_polymer_entities") or []):
            for c in ent.get("chains") or []:
                auth_of[c["label_asym_id"]] = c["auth_asym_id"]
                role[c["label_asym_id"]] = "target" if c["label_asym_id"] in targets else "partner"
        ns = gemmi.NeighborSearch(model, cell, cutoff + 1).populate(include_h=False)
        groups = {}
        for ch in model:
            for res in ch:
                if res.subchain not in targets or res.entity_type != gemmi.EntityType.Polymer:
                    continue
                for atom in res:
                    if atom.is_hydrogen():
                        continue
                    for mark in ns.find_atoms(atom.pos, "\0", radius=cutoff):
                        cra = mark.to_cra(model)
                        if cra.residue.entity_type != gemmi.EntityType.Polymer or cra.atom.is_hydrogen():
                            continue
                        # the image of THIS mark (find_nearest_image would give one image per pair and
                        # lose residues that touch two copies: 27 instead of 31 for 3WD5 A-A z,x,y)
                        im = cell.find_nearest_pbc_image(atom.pos, cra.atom.pos, mark.image_idx)
                        if im.dist() > cutoff:
                            continue
                        code = im.symmetry_code(True)
                        if code == "1_555" and cra.residue.subchain == res.subchain:
                            continue                       # same chain, same copy
                        key = (res.subchain, cra.residue.subchain, code)
                        g = groups.setdefault(key, {"im": im, "target": {}, "partner": set()})
                        rid = f"{res.seqid.num}{res.seqid.icode.strip()}"
                        g["target"][rid] = (res.name, clean._uniprot_pos(res.label_seq, uniprot_aln.get(res.subchain)))
                        g["partner"].add(f"{cra.residue.seqid.num}{cra.residue.seqid.icode.strip()}")
    except (RuntimeError, ValueError, OSError, IndexError) as exc:
        log.exception("Contacts failed for %s", cif)
        return {"ok": False, "error": f"contacts failed: {exc!r}"}

    contacts = []
    for (tl, pl, code), g in sorted(groups.items(), key=lambda kv: (auth_of.get(kv[0][0]) or "", auth_of.get(kv[0][1]) or "", kv[0][2])):
        im = g["im"]
        op = None
        if im.sym_idx < len(ops):
            # lattice shift folded in, as PISA writes it (2W2N "2_544": -y+1/2,x-1/2,z-1/4)
            shift = gemmi.Op("x{:+d},y{:+d},z{:+d}".format(*im.pbc_shift)) if any(im.pbc_shift) else None
            op = (shift.combine(ops[im.sym_idx]) if shift else ops[im.sym_idx]).triplet()
        elif code == "1_555":
            op = "x,y,z"

        def order(rid):
            m = re.match(r"-?\d+", rid)
            return (int(m.group(0)) if m else 0, rid)
        residues = [{"auth_seq_id": rid, "resname": g["target"][rid][0], "uniprot": g["target"][rid][1]}
                    for rid in sorted(g["target"], key=order)]
        contacts.append({
            "target_chain": auth_of.get(tl, tl), "target_label": tl,
            "partner_chain": auth_of.get(pl, pl), "partner_label": pl, "partner_role": role.get(pl, "partner"),
            "symmetry": {"code": code, "operator": op, "lattice_shift": list(im.pbc_shift)},
            "in_asymmetric_unit": code == "1_555",
            "n_target_residues": len(residues), "target_residues": residues,
            "n_partner_residues": len(g["partner"]),
            "_partner_ids": sorted(g["partner"], key=order),     # for the PISA match only, not stored
        })
    return {"ok": True, "contacts": contacts, "spacegroup": st.spacegroup_hm or None,
            "crystal": bool(cell.is_crystal())}


def _op(text):
    try:
        return gemmi.Op(str(text).lower().replace(" ", ""))
    except (RuntimeError, ValueError):
        return None


def match_pisa(contact: dict, pisa_interfaces: list):
    """The PISA interface with the same chain pair AND the same relative symmetry operator, or None.
    PISA writes an operator for each molecule (3WD5 #4: H x,y,z, A Z,X,Y), so the operator that
    takes the target to the partner is inv(op_target) * op_partner; for a chain with itself the
    inverse is the same interface (3WD5 #2 is both z,x,y and y,z,x). Reviewer 2026-10-28: matching by
    shared residues alone gave a 1-residue lattice contact of 4Y6O (C-A x,y-1,z) the css 1 of the
    large interface; PISA's own interface for it is #13, css 0, 16 Å²."""
    ours = _op(contact["symmetry"]["operator"]) if contact["symmetry"]["operator"] else None
    if ours is None:
        return None
    same_chain = contact["target_chain"] == contact["partner_chain"]
    best = None
    for itf in pisa_interfaces or []:
        mols = itf["molecules"]
        for i, t in enumerate(mols):
            p = mols[1 - i]
            if t["chain"] != contact["target_chain"] or p["chain"] != contact["partner_chain"]:
                continue
            op_t, op_p = _op(t["symop"]), _op(p["symop"])
            if op_t is None or op_p is None:
                continue
            rel = op_t.inverse().combine(op_p)
            direct = rel.triplet() == ours.triplet()
            if direct or (same_chain and rel.inverse().triplet() == ours.triplet()):
                # a direct match names the target-side molecule; an inverse one (same chain) the other
                # side, whose residues are not ours (reviewer 2026-10-28: 3WD5 A-A y,z,x said 5 of 29,
                # 28 of 29 on the right side): prefer direct, then the larger area
                rank = (direct, itf["area"] or 0)
                if best is None or rank > best[0]:
                    best = (rank, itf, t if direct else p)
    if not best:
        return None
    _, itf, t = best
    ours_res = {r["auth_seq_id"] for r in contact["target_residues"]}
    return {"id": itf["id"], "css": itf["css"], "area": itf["area"],
            "operators": [f"{m['chain']} {m['symop']}" for m in itf["molecules"]],
            # agreement check, not the matching rule: PISA lists residues with buried surface
            "target_residues_in_pisa": f"{len(ours_res & set(t['residues']))} of {contact['n_target_residues']}"}


def _sequences(record: dict, seqs: dict, accession: str) -> list:
    out = []
    for role, ents in (("target", record.get("target_entities") or []), ("partner", record.get("other_polymer_entities") or [])):
        for ent in ents:
            s = seqs.get(str(ent.get("entity_id"))) or {}
            seq = s.get("sequence")
            out.append({
                "entity_id": ent.get("entity_id"), "role": role,
                "description": ent.get("description") if role == "partner" else f"target ({accession})",
                "uniprot_ids": ent.get("uniprot_ids") if role == "partner" else [accession],
                "chains": [{"label": c["label_asym_id"], "auth": c["auth_asym_id"]} for c in ent.get("chains") or []],
                "polymer_type": s.get("polymer_type"),
                "sequence": seq, "length": len(seq) if seq else None,
                "flags": sequence_flags(seq or ""),
            })
    return out


def _sha(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


async def analyze_interfaces(uniprot_id: str, pdb_ids=None) -> dict:
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    path = choice.manifest_path(accession)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {"ok": False, "uniprot_id": accession, "error": f"No preparation manifest for {accession}: run prepare_structures first."}
    except (OSError, ValueError) as exc:
        return {"ok": False, "uniprot_id": accession, "error": f"Could not read {path}: {exc}"}
    prep = path.parent
    listed = [s.get("pdb_id") for s in manifest.get("structures") or []]
    wanted = [str(p).strip().upper() for p in pdb_ids if str(p).strip()] if pdb_ids else None
    chosen, skipped = [], []
    for s in manifest.get("structures") or []:
        if wanted is not None and s.get("pdb_id") not in wanted:
            continue
        f = s.get("file") or {}
        raw = (prep / str(f.get("path") or "")).resolve()
        if not f.get("ok") or prep.resolve() not in raw.parents or not raw.is_file():
            skipped.append({"pdb_id": s.get("pdb_id"), "reason": f.get("error") or "raw file not in prep/"})
        elif not s.get("target_entities"):
            skipped.append({"pdb_id": s.get("pdb_id"), "reason": s.get("chains_error") or "no target chain"})
        elif _sha(raw) != f.get("sha256"):
            skipped.append({"pdb_id": s.get("pdb_id"), "reason": "raw file changed since the manifest (SHA-256): run prepare_structures again"})
        else:
            chosen.append((s, raw))
    not_in_manifest = [p for p in wanted if p not in listed] if wanted else []
    if not chosen:
        return {"ok": False, "uniprot_id": accession, "error": "No structure to analyse.", "skipped": skipped,
                "pdb_ids_not_in_manifest": not_in_manifest}

    ids = [s["pdb_id"] for s, _ in chosen]
    gate = asyncio.Semaphore(PARALLEL_PISA)

    async def one_pisa(client, pdb_id):
        async with gate:
            return await pisa.fetch_pisa(client, pdb_id)
    async with httpx.AsyncClient(timeout=120.0, follow_redirects=True) as client:
        rcsb_result, *pisa_results = await asyncio.gather(
            rcsb.fetch_assemblies(client, ids), *(one_pisa(client, i) for i in ids))
    computed = await asyncio.to_thread(lambda: [compute_contacts(raw, s) for s, raw in chosen])

    warnings, now = [], utc_now()
    if not rcsb_result.get("ok"):
        warnings.append(f"RCSB assemblies and sequences unknown: {rcsb_result.get('error')}")
    blocks = {}
    for (s, raw), p, c in zip(chosen, pisa_results, computed):
        pdb = s["pdb_id"]
        entry = (rcsb_result.get("entries") or {}).get(pdb) if rcsb_result.get("ok") else None
        if rcsb_result.get("ok") and entry is None:
            warnings.append(f"{pdb}: RCSB returned no assembly data.")
        for key in ("interfaces", "multimers"):
            if p.get(f"{key}_error"):
                warnings.append(f"{pdb}: PISA {key} not available ({p[f'{key}_error']}); contacts are still computed here.")
        if not c.get("ok"):
            warnings.append(f"{pdb}: {c['error']}")
        pisa_checked = p.get("interfaces") is not None
        contacts = []
        for ct in c.get("contacts") or []:
            # None + pisa_checked False = PISA not available; None + True = PISA has no interface there
            ct["pisa_interface"] = match_pisa(ct, p["interfaces"]) if pisa_checked else None
            ct["pisa_checked"] = pisa_checked
            ct.pop("_partner_ids", None)
            contacts.append(ct)
        blocks[pdb] = {
            "generated_by": "analyze_interfaces", "computed_at": now, "cutoff": CONTACT_CUTOFF,
            "raw_sha256": s["file"]["sha256"],
            "sequences": _sequences(s, (entry or {}).get("sequences") or {}, accession) if entry else None,
            "rcsb_assemblies": (entry or {}).get("assemblies") if entry else None,
            "pisa": {
                "sources": p.get("sources"),
                "interfaces": [{k: i[k] for k in ("id", "area", "css", "n_occ")} | {
                    "molecules": [f"{m['chain']} {m['symop']}" for m in i["molecules"]]} for i in p["interfaces"]]
                if pisa_checked else None,
                "interfaces_error": p.get("interfaces_error"), "interfaces_note": p.get("interfaces_note"),
                # each set divides all chains into assemblies; the sets are PISA's alternatives
                "multimer_sets": p.get("sets"), "total_assemblies": p.get("total_assemblies"),
                "multimers_error": p.get("multimers_error"),
            },
            "contacts": contacts if c.get("ok") else None,
            "contacts_error": None if c.get("ok") else c["error"],
            "spacegroup": c.get("spacegroup"),
            "symmetry_copies_checked": c.get("crystal"),
        }

    # write: re-read so that a choice recorded meanwhile is not lost; refuse if the manifest was rebuilt
    try:
        current = json.loads(path.read_text(encoding="utf-8"))
        if current.get("created_at") != manifest.get("created_at"):
            return {"ok": False, "uniprot_id": accession,
                    "error": "prepare_structures rewrote the manifest during the analysis; run analyze_interfaces again."}
        for s in current.get("structures") or []:
            if s.get("pdb_id") in blocks:
                s["interfaces"] = blocks[s["pdb_id"]]
        choice.write_manifest(path, current)
    except (OSError, ValueError) as exc:
        log.exception("Could not write %s", path)
        return {"ok": False, "uniprot_id": accession, "error": f"The manifest could not be written: {exc!r}"}

    return {
        "ok": True, "uniprot_id": accession, "manifest_path": str(path),
        "selection": f"{len(blocks)} of {len(listed)} structures in the manifest"
                     + (f" (requested: {', '.join(wanted)})" if wanted else ""),
        "structures": [compact(pdb, b) for pdb, b in blocks.items()],
        "skipped": skipped, "pdb_ids_not_in_manifest": not_in_manifest, "warnings": warnings,
        "notes": [
            f"Contacts: target residues with a heavy atom within {CONTACT_CUTOFF} Å of another polymer chain, "
            "including crystal-symmetry copies (symmetry code 1_555 = the deposited coordinates). A contact is "
            "geometry only: it is not called biological or crystal packing here.",
            "rcsb_assemblies: the assembly as RCSB gives it, with who defined it; equals_asymmetric_unit true means "
            "no symmetry copy is in it. PISA (PDBe) gives its own opinion (css, multimers); the sources can disagree. "
            "A contact is matched to the PISA interface with the same chain pair and the same symmetry operator.",
            "PISA multimer sets: each set divides all chains of the crystal into assemblies that coexist; the sets "
            "are PISA's alternative solutions.",
            "Sequences are those of the crystallised constructs (tags, single domains, non-standard residues).",
        ],
    }


def compact(pdb: str, b: dict) -> dict:
    """What the chat needs per structure; the manifest keeps residues and sequences in full."""
    def contact_line(c):
        pi = c.get("pisa_interface")
        nums = [str(r["uniprot"]) if r["uniprot"] is not None else f"{r['auth_seq_id']}(auth)" for r in c["target_residues"]]
        return (f"{c['target_chain']}-{c['partner_chain']} ({c['partner_role']}) {c['symmetry']['code']} "
                f"[{c['symmetry']['operator']}]" + (" asymmetric unit" if c["in_asymmetric_unit"] else " symmetry copy")
                + f": {c['n_target_residues']} target residues (UniProt {', '.join(nums)})"
                + (f"; PISA interface {pi['id']} css {pi['css']} area {pi['area']} ({pi['target_residues_in_pisa']} in PISA)"
                   if pi else "; PISA lists no interface for this chain pair and operator" if c.get("pisa_checked")
                   else "; PISA not available"))
    pisa_b = b.get("pisa") or {}
    return {
        "pdb_id": pdb,
        "partners": [f"{s['role']} entity {s['entity_id']}: {s['description']}; chains "
                     + ", ".join(c["auth"] for c in s["chains"]) + f"; UniProt {', '.join(s['uniprot_ids'] or []) or 'none'}"
                     + f"; length {s['length']}" + "".join(f"; {f}" for f in s["flags"][1:])
                     for s in b.get("sequences") or []],
        "rcsb_assemblies": [f"{a['id']}: {a['defined_by']}" + (f" ({a['method']})" if a.get("method") else "")
                            + f", polymer chains {', '.join(a['polymer_chains'])}, operators {', '.join(a['operators'])}"
                            + ("; same as the asymmetric unit" if a["equals_asymmetric_unit"] else "")
                            for a in b.get("rcsb_assemblies") or []] if b.get("rcsb_assemblies") is not None else None,
        # one line per set: its assemblies coexist; the sets are PISA's alternatives
        "pisa_multimer_sets": [f"set {s['set']}: " + " + ".join(f"{a['composition']} ({a['verdict']})" for a in s["assemblies"])
                               for s in pisa_b.get("multimer_sets") or []] if pisa_b.get("multimer_sets") is not None else None,
        "pisa_errors": [e for e in (pisa_b.get("interfaces_error"), pisa_b.get("multimers_error")) if e] or None,
        "contacts": [contact_line(c) for c in b.get("contacts") or []],
        "contacts_error": b.get("contacts_error"),
        "symmetry_copies_checked": b.get("symmetry_copies_checked"),
    }


def carry_over(old: dict, new: dict) -> dict:
    """Called by prepare_structures: keep a structure's `interfaces` if its raw file is the same."""
    before = {s.get("pdb_id"): s for s in (old or {}).get("structures") or [] if s.get("interfaces")}
    kept, dropped = [], []
    for s in new.get("structures") or []:
        o = before.pop(s.get("pdb_id"), None)
        if not o:
            continue
        if (s.get("file") or {}).get("ok") and o["interfaces"].get("raw_sha256") == s["file"].get("sha256"):
            s["interfaces"] = o["interfaces"]
            kept.append(s["pdb_id"])
        else:
            dropped.append(f"{s['pdb_id']}: structure file changed or not downloaded")
    dropped += [f"{p}: not in the new manifest" for p in before]
    return {"kept": kept, "dropped": dropped}
