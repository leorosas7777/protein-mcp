"""prepare_structures: download the structures chosen by get_structures and write the
manifest, the contract with Phase 3 (dossiers/<acc>/prep/manifest.json).

Step 1 of Phase 2: original files only (mmCIF, AlphaFold mmCIF, UniProt FASTA, crystal
ligand poses as SDF). Nothing is cleaned, protonated or completed here; the manifest
says so in `pending`. Imports common, dossier, files.download and sources.{rcsb, uniprot}.
"""
import asyncio
import json
import re
from pathlib import Path

import httpx

from protein_mcp import dossier
from protein_mcp.common import log, parse_accession, utc_now
from protein_mcp.files.download import check_fasta, check_mmcif, check_sdf, download
from protein_mcp.sources import rcsb
from protein_mcp.sources.uniprot import UNIPROT_FASTA

MANIFEST_VERSION = 1
PARALLEL_DOWNLOADS = 4
_PDB_ID_RE = re.compile(r"[0-9][A-Z0-9]{3}")
_AF_FILE_RE = re.compile(r"https://alphafold\.ebi\.ac\.uk/files/(AF-[A-Za-z0-9-]+-model_v\d+\.cif)")

DECISIONS = [
    "Original files are kept unmodified in raw/. Cleaning (target chains only, no waters or "
    "additives), pocket and box are step 2; loop completion and protonation are Phase 3.",
    "Target chains come from RCSB's chain-to-UniProt mapping, not from the mmCIF _struct_ref, "
    "which is unreliable (3OG7 cites UniProt Q5IBP5; antibody chains often cite no UniProt).",
    "Chains are given as label_asym_id (mmCIF) and auth_asym_id (author, used by PyMOL and the "
    "PDB format); they can differ (7ZR0: BRAF is label D, auth K).",
    "Ligand poses are the crystal coordinates of each instance (RCSB ModelServer, by "
    "label_asym_id). A ligand can sit on one chain only (3OG7: chain A) or touch several.",
    "neighbor_*_chains come from RCSB rcsb_target_neighbors with the closest distance per "
    "chain; RCSB does not document the distance cutoff. The pocket is computed in step 2.",
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
    "step 2: cleaned PDB of the target chains, structure metrics, pocket and box",
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
    return {"ok": True, "path": _rel(dl["path"], base), "url": dl["url"], "bytes": dl["bytes"],
            "sha256": dl["sha256"], "retrieved_at": dl["retrieved_at"]}


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

    alphafold = None
    if af_match:
        alphafold = {key: af.get(key) for key in (
            "entry_id", "isoform", "is_canonical_sequence", "selection_reason", "mean_plddt",
            "confidence_band", "residues", "model_version", "created")}
        alphafold["file"] = _file_record(downloads["alphafold"], base)
    fasta = _file_record(downloads["fasta"], base)
    if fasta["ok"]:
        text = (base / fasta["path"]).read_text(encoding="ascii", errors="replace").splitlines()
        fasta["header"] = text[0] if text else None
        fasta["length"] = sum(len(line.strip()) for line in text[1:])

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
        "sequence": fasta,
        "structures": structures,
        "alphafold": alphafold,
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
        "pdb_ids_not_in_dossier": not_in_dossier,
        "replaced_manifest_from": previous,
        "sequence": {k: fasta.get(k) for k in ("ok", "path", "header", "length", "error") if k in fasta},
        "structures": [compact_structure(r) for r in structures],
        "alphafold": ({"entry_id": alphafold["entry_id"], "file": alphafold["file"].get("path"),
                       "ok": alphafold["file"]["ok"], "mean_plddt": alphafold["mean_plddt"]}
                      if alphafold else None),
        "failed_downloads": failed,
        "warnings": warnings,
        "notes": [
            "Files are the ORIGINAL downloads: nothing was cleaned, protonated or completed.",
            "Chains: label_asym_id/auth_asym_id from RCSB; ligand neighbours from RCSB rcsb_target_neighbors.",
            "Full provenance (URL, size, SHA-256, time) and decisions are in the manifest.",
        ],
    }


def compact_structure(r: dict) -> dict:
    """What the chat needs per structure; the manifest keeps the rest."""
    targets = r.get("target_entities")
    ligands = None
    if r.get("nonpolymer_instances") is not None:
        ligands = [{
            "comp_id": x["comp_id"], "category": x["category"], "label_asym_id": x["label_asym_id"],
            "auth_asym_id": x["auth_asym_id"],
            "neighbor_target_auth_chains": sorted({v["auth_asym_id"] or k for k, v in x["neighbor_target_chains"].items()}),
            "touches_other_chains": bool(x["neighbor_other_chains"]),
            "sdf": (x["sdf"] or {}).get("path"),
        } for x in r["nonpolymer_instances"] if x["category"] != "additive"]
    out = {
        "pdb_id": r["pdb_id"],
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
    if r.get("chains_error"):
        out["chains_error"] = r["chains_error"]
    if not r["file"]["ok"]:
        out["download_error"] = r["file"].get("error")
    return out
