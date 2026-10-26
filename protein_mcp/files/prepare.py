"""prepare_structures: download the structures chosen by get_structures and write the
manifest, the contract with Phase 3 (dossiers/<acc>/prep/manifest.json).

Step 1 of Phase 2: original files (mmCIF, AlphaFold mmCIF, UniProt FASTA, crystal ligand
poses as SDF). Step 2 (files.clean, gemmi): cleaned receptor, chain metrics, pocket and box per
ligand instance, AlphaFold pLDDT at each pocket. Nothing is protonated, completed or mutated;
the manifest says so in `pending`. Imports common, dossier, files.{download, clean} and
sources.{rcsb, uniprot}.
"""
import asyncio
import json
import re
from pathlib import Path

import httpx

from protein_mcp import dossier
from protein_mcp.common import log, parse_accession, short, utc_now
from protein_mcp.files import clean
from protein_mcp.files.download import check_fasta, check_mmcif, check_sdf, download
from protein_mcp.sources import rcsb
from protein_mcp.sources.uniprot import UNIPROT_FASTA

MANIFEST_VERSION = 1
PARALLEL_DOWNLOADS = 4
_PDB_ID_RE = re.compile(r"[0-9][A-Z0-9]{3}")
_AF_FILE_RE = re.compile(r"https://alphafold\.ebi\.ac\.uk/files/(AF-[A-Za-z0-9-]+-model_v\d+\.cif)")

DECISIONS = [
    "Original files are kept unmodified in raw/. Loop completion, protonation and charges are Phase 3.",
    "Receptor (clean/<PDB>_receptor.pdb): target chains only (RCSB mapping), first model, first "
    "alternate conformation, no hydrogens, no waters, no additives, no ligands; nucleotides, "
    "cofactors and metal ions MG/MN/ZN/FE/CO/NI/CU that touch a target chain are KEPT (e.g. KRAS "
    "GDP + Mg, 8TVK); NA, K, CL, CA, CD are removed. Coordinates are unchanged. PDB "
    "format unless chain names or atom count do not fit, then mmCIF.",
    "Pocket: target-chain residues with any atom within 5 Å of the ligand instance (gemmi, first "
    "alternate conformation of the ligand), in the asymmetric unit only: crystal-symmetry copies are "
    "excluded (they put residues ~47 Å away in 8C7Y). UniProt numbers come from label_seq_id + RCSB "
    "alignment; residues of other polymer chains are listed apart. No pocket for nucleotides and "
    "cofactors (kept in the receptor), for ligands with fewer than 10 heavy atoms, or without a "
    "target-chain neighbour in RCSB: each is listed in pocket_not_computed. Box: ligand bounding box "
    "+ 5 Å per side, same frame as the mmCIF; sides above 30 Å are flagged. A ligand instance with "
    "pocket residues in more than one target chain is flagged spans_target_chains (e.g. a bivalent "
    "compound bridging two sites): its pocket and box cover every chain it touches and are NOT split "
    "into one box per site; residues_by_chain gives the count per chain.",
    "Chain metrics: missing_uniprot_ranges are residues of the sample sequence not modelled inside "
    "the modelled span; a missing range next to a pocket residue is flagged. variant_check reads "
    "the residue actually modelled at the variant position in each target chain and compares it "
    "with the UniProt canonical letter (not the letter typed in the filter): the variant can be "
    "named in the title and still be unmodelled in the chain that holds the ligand (3OG7 A). Two "
    "residues at one position (microheterogeneity) give status 'ambiguous'.",
    "AlphaFold pLDDT (B-factor column, equal to _ma_qa_metric_local) is read at the same UniProt "
    "numbers as each experimental pocket, only for the canonical-sequence model; residues without "
    "a UniProt number are not compared (alphafold_residues_compared vs n_residues).",
    "Target chains come from RCSB's chain-to-UniProt mapping, not from the mmCIF _struct_ref, "
    "which is unreliable (3OG7 cites UniProt Q5IBP5; antibody chains often cite no UniProt).",
    "Chains are given as label_asym_id (mmCIF) and auth_asym_id (author, used by PyMOL and the "
    "PDB format); they can differ (7ZR0: BRAF is label D, auth K).",
    "Ligand poses are the crystal coordinates of each instance (RCSB ModelServer, by "
    "label_asym_id). A ligand can sit on one chain only (3OG7: chain A) or touch several.",
    "neighbor_*_chains come from RCSB rcsb_target_neighbors with the closest distance per "
    "chain; RCSB does not document the distance cutoff.",
    "construct.mutations are the depositor's declared mutations of the target entity; they can "
    "include solubility/engineering mutations besides the variant (3OG7 declares 13). Any declared "
    "mutation means the structure is not the native sequence; an empty list does NOT prove it is "
    "native (the variant may be named only in the title, or the text may be unparsed, e.g. 'YES').",
    "Additives (ions, buffers, cryoprotectants) get no SDF. The classification is the one of "
    "get_structures; some additives (e.g. BEN, EOH in 5JSM) are not in its list and count as ligands.",
    "Glycans (RCSB branched entities, e.g. NAG chains) are not queried and are absent from "
    "nonpolymer_instances.",
]
PENDING = [
    "step 3: PyMOL session (.pml)",
    "step 4: user's structure choice (select_structure)",
]


def prep_dir(accession: str) -> Path:
    return dossier.DOSSIER_DIR / accession / "prep"


def _rel(path, base: Path):
    return Path(path).relative_to(base).as_posix() if path else None


def _file_record(dl: dict, base: Path) -> dict:
    """Provenance of one download, with the path relative to prep/."""
    if not dl.get("ok"):
        return {"ok": False, "url": dl.get("url"), "error": dl.get("error")}
    out = {"ok": True, "path": _rel(dl["path"], base), "url": dl["url"], "bytes": dl["bytes"],
           "sha256": dl["sha256"], "retrieved_at": dl["retrieved_at"]}
    if dl.get("attempts"):
        out["attempts"] = dl["attempts"]
    return out


async def _run(jobs: dict, client) -> dict:
    """{key: (url, path, check)} -> {key: download result}, at most PARALLEL_DOWNLOADS at once."""
    gate = asyncio.Semaphore(PARALLEL_DOWNLOADS)

    async def one(url, path, check):
        async with gate:
            return await download(client, url, path, check)

    keys = list(jobs)
    results = await asyncio.gather(*(one(*jobs[k]) for k in keys))
    return dict(zip(keys, results))


def _available_structure_files(accession: str) -> list:
    folder = dossier.DOSSIER_DIR / accession
    return sorted(p.name for p in folder.glob("structures*.json")) if folder.is_dir() else []


async def prepare_structures(uniprot_id: str, mutation: str = "", pdb_ids=None,
                             include_alphafold: bool = True) -> dict:
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    mutation = (mutation or "").strip()
    name = dossier.dossier_name("structures", mutation)
    saved = dossier.read_dossier(accession, name)
    if saved is None:
        call = f'get_structures("{accession}"' + (f', mutation="{mutation}")' if mutation else ")")
        return {"ok": False, "error": f"No dossier file {name}.json for {accession}. Run {call} first.",
                "available_structure_files": _available_structure_files(accession)}
    result = saved.get("result") or {}
    exp = result.get("experimental") or {}
    in_dossier = [s for s in exp.get("structures") or [] if _PDB_ID_RE.fullmatch(str(s.get("pdb_id") or ""))]
    dossier_ids = [s["pdb_id"] for s in in_dossier]

    not_in_dossier = []
    if pdb_ids:
        wanted = [str(p).strip().upper() for p in pdb_ids if str(p).strip()]
        not_in_dossier = [p for p in wanted if p not in dossier_ids]
        chosen = [s for s in in_dossier if s["pdb_id"] in wanted]
    else:
        chosen = in_dossier
    warnings = []
    if not exp.get("ok"):
        warnings.append(f"RCSB failed when {name}.json was made ({exp.get('error') or 'no experimental result'}): "
                        "its experimental structures are unknown, not zero. Rerun get_structures.")
    af = result.get("alphafold") or {}
    af_url = str((af.get("links") or {}).get("cif") or "")
    af_match = _AF_FILE_RE.fullmatch(af_url) if include_alphafold else None
    if include_alphafold and not af_match:
        if not af.get("ok"):
            warnings.append(f"AlphaFold not used: it failed when {name}.json was made "
                            f"({af.get('error') or 'no AlphaFold result'}). Rerun get_structures.")
        elif af.get("available") is False:
            warnings.append("AlphaFold not used: AlphaFold DB has no model for this accession (per the dossier).")
        else:
            warnings.append("AlphaFold not used: the dossier's model URL is not a recognised AlphaFold DB "
                            f"mmCIF URL ({af_url or 'none'}).")
    if not chosen and not af_match:
        return {"ok": False, "uniprot_id": accession, "source_dossier": f"{name}.json",
                "error": "Nothing to prepare: none of the requested PDB IDs is in this dossier file"
                         + ("" if pdb_ids else " (it has no structures)") + " and no AlphaFold model is used.",
                "warnings": warnings, "pdb_ids_in_dossier": dossier_ids, "pdb_ids_not_in_dossier": not_in_dossier}

    base = prep_dir(accession)
    ids = [s["pdb_id"] for s in chosen]
    async with httpx.AsyncClient(timeout=120.0, follow_redirects=True) as client:
        chains = await rcsb.fetch_chains(client, ids, accession) if ids else {"ok": True, "entries": {}}
        if not chains.get("ok"):
            return {"ok": False, "uniprot_id": accession, "source_dossier": f"{name}.json",
                    "error": f"{chains.get('error')}. Nothing was downloaded and the previous manifest "
                             "(if any) was kept.", "warnings": warnings}
        entries = chains.get("entries") or {}
        jobs = {("cif", i): (rcsb.cif_url(i), base / "raw" / f"{i}.cif", check_mmcif(i)) for i in ids}
        jobs["fasta"] = (UNIPROT_FASTA.format(accession=accession), base / "raw" / f"{accession}.fasta",
                         check_fasta(accession))
        if af_match:
            jobs["alphafold"] = (af_match.group(0), base / "raw" / af_match.group(1), check_mmcif())
        for i in ids:
            for inst in (entries.get(i) or {}).get("nonpolymer_instances") or []:
                if inst["category"] != "additive":
                    label = inst["label_asym_id"]
                    jobs[("sdf", i, label)] = (rcsb.ligand_sdf_url(i, label),
                                               base / "ligands" / f"{i}_{inst['comp_id']}_{label}.sdf", check_sdf)
        downloads = await _run(jobs, client)

    structures = []
    for s in chosen:
        i = s["pdb_id"]
        cif = _file_record(downloads[("cif", i)], base)
        chain_info = entries.get(i)
        record = {
            "pdb_id": i,
            "title": s.get("title"),
            "method": s.get("method"),
            "resolution_angstrom": s.get("resolution_angstrom"),
            "release_date": s.get("release_date"),
            "links": s.get("links") or {},
            "file": cif,
            "construct": {
                "mutations": s.get("mutations") or [],
                "mutations_from_title": s.get("mutations_from_title") or [],
                "mutations_unparsed": s.get("mutations_unparsed") or [],
            },
        }
        if chain_info is None:
            record["chains_error"] = "RCSB did not return chain data for this entry."
            record.update(target_entities=None, other_polymer_entities=None, nonpolymer_instances=None)
        else:
            instances = []
            for inst in chain_info["nonpolymer_instances"]:
                row = dict(inst)
                dl = downloads.get(("sdf", i, inst["label_asym_id"]))
                row["sdf"] = _file_record(dl, base) if dl else None
                instances.append(row)
            record.update(target_entities=chain_info["target_entities"],
                          other_polymer_entities=chain_info["other_polymer_entities"],
                          nonpolymer_instances=instances)
            if not chain_info["target_entities"]:
                record["chains_error"] = f"RCSB maps no chain of this entry to {accession}."
        if record.get("chains_error"):
            warnings.append(f"{i}: {record['chains_error']}")
        structures.append(record)

    failed = [f"{k if isinstance(k, str) else '/'.join(k)}: {d.get('error')}"
              for k, d in downloads.items() if not d.get("ok")]
    if not any(r["file"]["ok"] for r in structures) and not (af_match and downloads["alphafold"].get("ok")):
        return {"ok": False, "uniprot_id": accession, "source_dossier": f"{name}.json",
                "error": "Every structure download failed; the previous manifest (if any) was kept.",
                "failed_downloads": failed, "warnings": warnings}

    fasta = _file_record(downloads["fasta"], base)
    canonical = ""
    if fasta["ok"]:
        text = (base / fasta["path"]).read_text(encoding="ascii", errors="replace").splitlines()
        fasta["header"] = text[0] if text else None
        canonical = "".join(line.strip() for line in text[1:])
        fasta["length"] = len(canonical)

    to_clean = [r for r in structures if r["file"]["ok"] and r.get("target_entities")]
    cleaned = await asyncio.to_thread(
        lambda: [clean.clean_structure(base / r["file"]["path"], base / "clean", r, exp.get("mutation_filter"),
                                       canonical)
                 for r in to_clean])
    for r in structures:
        r["prepared"] = None
    for r, c in zip(to_clean, cleaned):
        if c["ok"]:
            c["receptor"]["path"] = _rel(c["receptor"]["path"], base)
        else:
            warnings.append(f"{r['pdb_id']}: {c['error']}")
        r["prepared"] = c

    alphafold = None
    if af_match:
        alphafold = {key: af.get(key) for key in (
            "entry_id", "isoform", "is_canonical_sequence", "selection_reason", "mean_plddt",
            "confidence_band", "residues", "model_version", "created")}
        alphafold["file"] = _file_record(downloads["alphafold"], base)
        pockets = {r["pdb_id"]: r["prepared"]["pockets"] for r in structures if (r["prepared"] or {}).get("ok")}
        alphafold["pocket_plddt"] = None
        if alphafold["file"]["ok"] and af.get("is_canonical_sequence") is not True:
            warnings.append("AlphaFold pLDDT not compared with the pockets: the model is not the canonical sequence.")
        elif alphafold["file"]["ok"]:
            m = await asyncio.to_thread(clean.alphafold_metrics, base / alphafold["file"]["path"], pockets,
                                        exp.get("mutation_filter"), canonical)
            if not m["ok"]:
                warnings.append(m["error"])
            alphafold["pocket_plddt"] = m.get("pockets")
            alphafold["variant_position"] = m.get("variant_position")
            by_key = {(x["pdb_id"], x["label_asym_id"]): x for x in m.get("pockets") or []}
            for pdb_id, plist in pockets.items():
                for pk in plist:
                    hit = by_key.get((pdb_id, pk["label_asym_id"]))
                    if hit:
                        pk["alphafold_plddt"] = {k: hit[k] for k in ("mean_plddt", "residues_below_70", "residues_compared")}
    summary = summarize(structures)
    variant = clean.parse_variant(exp.get("mutation_filter"), canonical)
    if variant and variant[3]:
        warnings.append(f"Variant check: {variant[3]}.")

    manifest_path = base / "manifest.json"
    previous = None
    try:
        previous = json.loads(manifest_path.read_text(encoding="utf-8")).get("source_dossier")
    except FileNotFoundError:
        pass
    except (OSError, ValueError, AttributeError):
        previous = {"file": "unreadable previous manifest"}

    selection = (f"{len(chosen)} of {len(dossier_ids)} experimental structures in {name}.json"
                 + (f" (requested: {', '.join(str(p) for p in pdb_ids)})" if pdb_ids else "")
                 + (" + the AlphaFold model" if af_match else "")
                 + (f". The dossier file itself holds {exp.get('returned', len(dossier_ids))} of "
                    f"{exp.get('total_entries')} RCSB entries (see its `selection`)."
                    if exp.get("ok") else ". RCSB failed when the dossier file was made (see warnings)."))
    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "uniprot_id": accession,
        "created_at": utc_now(),
        "generated_by": "prepare_structures",
        "source_dossier": {"file": f"{name}.json", "tool": saved.get("tool"), "params": saved.get("params"),
                           "retrieved_at": saved.get("retrieved_at")},
        "mutation_filter": exp.get("mutation_filter"),
        "selection": selection,
        "summary": summary,
        "sequence": fasta,
        "structures": structures,
        "alphafold": alphafold,
        # the same list the chat gets: a failure is also in its own file record, but a reader of
        # the manifest must not have to scan every ligand instance to find one (6V34 QOP C)
        "failed_downloads": failed,
        "warnings": warnings,
        "decisions": DECISIONS,
        "pending": PENDING,
        "user_choice": None,
    }
    try:
        base.mkdir(parents=True, exist_ok=True)
        part = manifest_path.with_name("manifest.json.part")
        part.write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
        part.replace(manifest_path)
    except OSError as exc:
        log.exception("Could not write manifest %s", manifest_path)
        return {"ok": False, "error": f"Files downloaded to {base} but the manifest could not be written: {exc!r}"}

    return {
        "ok": True,
        "uniprot_id": accession,
        "manifest_path": str(manifest_path),
        "prep_dir": str(base),
        "source_dossier": f"{name}.json",
        "selection": selection,
        "summary": summary,
        "pdb_ids_not_in_dossier": not_in_dossier,
        "replaced_manifest_from": previous,
        "sequence": {k: fasta.get(k) for k in ("ok", "path", "header", "length", "error") if k in fasta},
        "structures": [compact_structure(r) for r in structures],
        "alphafold": ({"entry_id": alphafold["entry_id"], "file": alphafold["file"].get("path"),
                       "ok": alphafold["file"]["ok"], "mean_plddt": alphafold["mean_plddt"],
                       "variant_position": alphafold.get("variant_position")}
                      if alphafold else None),
        "failed_downloads": failed,
        "warnings": warnings,
        "notes": [
            "raw/ holds the ORIGINAL downloads; clean/ holds receptors with other chains, waters, additives "
            "and ligands removed. Nothing was protonated, completed or mutated.",
            f"Pocket = target residues within {clean.POCKET_CUTOFF} Å of the ligand; box = ligand + "
            f"{clean.BOX_PADDING} Å per side, in the mmCIF frame.",
            "Chains: label_asym_id/auth_asym_id from RCSB; ligand neighbours from RCSB rcsb_target_neighbors.",
            "Full provenance (URL, size, SHA-256, time) and decisions are in the manifest.",
        ],
    }


def summarize(structures: list) -> dict:
    """Facts ACROSS entries, computed here so the chat does not have to cross them by hand.
    Desktop tests 2026-10-26: Claude twice said "the 597-614 gap is in all 7 entries" (8C7Y A
    has no gap) and that no pocket has V600E modelled in its own chain (8C7Y TXV A and 5JSM 6NB
    have). Sorting and grouping only: no score and no docking ranking. File-derived lists cover
    structures_compared only; pocket errors and ligands without a pocket are listed, not dropped."""
    done = [r for r in structures if (r.get("prepared") or {}).get("ok")]
    not_done = [r["pdb_id"] for r in structures if r not in done]

    def label(r, value):
        text = f"{r['pdb_id']} {value}" if value is not None else f"{r['pdb_id']} not reported"
        return text if r in done else text + " (not compared)"

    out = {
        "scope": (f"{len(done)} of {len(structures)} structures compared (downloaded and cleaned). "
                  "Every list below except by_resolution covers only structures_compared: an empty "
                  "list says nothing about structures_not_compared."),
        "structures_compared": [r["pdb_id"] for r in done],
        "structures_not_compared": not_done,
        # resolution comes from the dossier, so every entry is ranked, compared or not; best first
        "by_resolution": [label(r, r.get("resolution_angstrom")) for r in sorted(
            structures, key=lambda r: (r.get("resolution_angstrom") is None, r.get("resolution_angstrom") or 0))],
        "by_r_free": [label(r, v) for r, v in sorted(((r, r["prepared"]["metrics"]["r_free"]) for r in done),
                                                     key=lambda rv: (rv[1] is None, rv[1] or 0))],
        # no gap INSIDE the modelled span; unmodelled termini are not gaps (8C7Y A: sample 444-721,
        # modelled 448-721), so the modelled range is given with each chain
        "chains_without_internal_gaps": [
            f"{r['pdb_id']}:{c['auth_asym_id']} (modelled {c['modelled_uniprot_range'][0]}-{c['modelled_uniprot_range'][1]})"
            for r in done for c in r["prepared"]["metrics"]["chains"]
            if c["modelled_uniprot_range"] and not c["missing_uniprot_ranges"]],
        # ligand instances without a pocket (too small, no target neighbour in RCSB), so that
        # "every ligand" statements can see them
        "ligands_without_pocket": [f"{r['pdb_id']} {x['comp_id']} {x['label_asym_id']}: {x['reason']}"
                                   for r in done for x in r["prepared"]["pocket_not_computed"]],
    }
    chains = [(r, c) for r in done for c in r["prepared"]["metrics"]["chains"]]
    if not any("variant_check" in c for _, c in chains):
        return out                                 # no substitution filter: nothing to check
    by_status = {}
    for r, c in chains:
        status = (c.get("variant_check") or {}).get("status", "not checked")
        by_status.setdefault(status, []).append(f"{r['pdb_id']}:{c['auth_asym_id']}")
    out["variant_position_by_chain_status"] = by_status
    with_variant, ambiguous, without, unclassified = [], [], [], []
    for r in done:
        status = {c["auth_asym_id"]: (c.get("variant_check") or {}).get("status") for c in r["prepared"]["metrics"]["chains"]}
        for pk in r["prepared"]["pockets"]:
            name = f"{r['pdb_id']} {pk['comp_id']} {pk['label_asym_id']}"
            pocket_chains = list(pk.get("residues_by_chain") or {})
            if pk.get("error") or not pocket_chains:
                unclassified.append(f"{name}: {pk.get('error') or 'no target residue within the cutoff'}")
                continue
            row = {"pdb_id": r["pdb_id"], "comp_id": pk["comp_id"], "label_asym_id": pk["label_asym_id"],
                   "pocket_chains": pocket_chains,
                   "residues_by_chain": pk["residues_by_chain"],
                   # at least one pocket chain; the variant residue itself need not be a pocket residue
                   "variant_modelled_in": [c for c in pocket_chains if status.get(c) == "variant"]}
            if pk.get("spans_target_chains"):
                row["spans_target_chains"] = True
            if row["variant_modelled_in"]:
                with_variant.append(row)
            elif any(str(status.get(c)).startswith("ambiguous") for c in pocket_chains):
                ambiguous.append(f"{name} ({', '.join(f'{c}: {status.get(c)}' for c in pocket_chains)})")
            else:
                without.append(f"{name} ({', '.join(f'{c}: {status.get(c)}' for c in pocket_chains)})")
    out["pockets_with_variant_in_pocket_chain"] = with_variant
    out["pockets_with_ambiguous_variant_position"] = ambiguous
    out["pockets_without_variant_in_pocket_chain"] = without
    out["pockets_not_classified"] = unclassified
    return out


def compact_structure(r: dict) -> dict:
    """What the chat needs per structure; the manifest keeps the rest."""
    targets = r.get("target_entities")
    ligands = None
    if r.get("nonpolymer_instances") is not None:
        ligands = [{
            "comp_id": x["comp_id"], "category": x["category"],
            "name": short(x.get("name"), 100), "synonyms": x.get("synonyms") or [],
            "formula_weight": x.get("formula_weight"),
            "label_asym_id": x["label_asym_id"], "auth_asym_id": x["auth_asym_id"],
            "neighbor_target_auth_chains": sorted({v["auth_asym_id"] or k for k, v in x["neighbor_target_chains"].items()}),
            "touches_other_chains": bool(x["neighbor_other_chains"]),
            "sdf": (x["sdf"] or {}).get("path"),
        } for x in r["nonpolymer_instances"] if x["category"] != "additive"]
    out = {
        "pdb_id": r["pdb_id"],
        "title": r.get("title"),        # 5JSM: "... Chemically Linked Vemurafenib Inhibitor VEM-3-VEM"
        "resolution_angstrom": r.get("resolution_angstrom"),
        "file": r["file"].get("path"),
        "download_ok": r["file"]["ok"],
        "target_chains": None if targets is None else [
            {"label_asym_id": c["label_asym_id"], "auth_asym_id": c["auth_asym_id"]}
            for t in targets for c in t["chains"]],
        # UniProt range of the SAMPLE sequence (construct), not of the residues resolved in the model
        "sample_uniprot_ranges": None if targets is None else sorted({
            f"{a['uniprot_begin']}-{a['uniprot_begin'] + a['length'] - 1}"
            for t in targets for a in t.get("uniprot_alignment") or []
            if isinstance(a.get("uniprot_begin"), int) and isinstance(a.get("length"), int)}),
        "declared_mutations": r["construct"]["mutations"],
        "mutations_unparsed": r["construct"]["mutations_unparsed"],
        "other_polymer_entities": None if r.get("other_polymer_entities") is None else [
            o.get("description") for o in r["other_polymer_entities"]],
        "ligand_instances": ligands,
    }
    prep = r.get("prepared")
    if prep and prep.get("ok"):
        out["r_free"] = prep["metrics"]["r_free"]
        out["receptor"] = prep["receptor"]["path"]
        out["chain_checks"] = [{
            "auth_asym_id": c["auth_asym_id"], "modelled_uniprot_range": c["modelled_uniprot_range"],
            "missing_uniprot_ranges": c["missing_uniprot_ranges"],
            **({"variant": c["variant_check"]["status"]} if "variant_check" in c else {}),
            **({"variant_warning": c["variant_check"]["warning"]} if "warning" in c.get("variant_check", {}) else {}),
        } for c in prep["metrics"]["chains"]]
        out["pockets"] = [_compact_pocket(p) for p in prep["pockets"]]
        if prep["pocket_not_computed"]:
            out["pocket_not_computed"] = [f"{x['comp_id']} {x['label_asym_id']}: {x['reason']}"
                                          for x in prep["pocket_not_computed"]]
    elif prep:
        out["clean_error"] = prep.get("error")
    if r.get("chains_error"):
        out["chains_error"] = r["chains_error"]
    if not r["file"]["ok"]:
        out["download_error"] = r["file"].get("error")
    return out


def _compact_pocket(p: dict) -> dict:
    """Pocket for the chat: residues as 'chain:<one letter><UniProt number>' (auth number if unmapped)."""
    if p.get("error"):
        return {"comp_id": p["comp_id"], "label_asym_id": p["label_asym_id"], "error": p["error"]}
    names = []
    for x in p["residues"]:
        code = clean._one_letter(x["resname"]) or x["resname"]
        number = x["uniprot"] if x.get("uniprot") is not None else f"{x['auth_seq_id']}(auth)"
        names.append(f"{x['auth_asym_id']}:{code}{number}")
    out = {
        "comp_id": p["comp_id"], "label_asym_id": p["label_asym_id"], "auth_asym_id": p["auth_asym_id"],
        "heavy_atoms": p["heavy_atoms"],
        "n_residues": len(p["residues"]),
        "residues_by_chain": p.get("residues_by_chain"),
        "residues": " ".join(names),
        "partner_residues": len(p["partner_residues"]),
        "box_center": p["box"]["center"], "box_size": p["box"]["size"],
        "adjacent_missing_ranges": p["adjacent_missing_ranges"],
    }
    if p.get("large_box"):
        out["large_box"] = True
    if p.get("spans_target_chains"):
        out["spans_target_chains"] = True
    if p.get("alphafold_plddt"):
        out["alphafold_mean_plddt"] = p["alphafold_plddt"]["mean_plddt"]
        out["alphafold_residues_compared"] = p["alphafold_plddt"]["residues_compared"]   # of n_residues
    return out
