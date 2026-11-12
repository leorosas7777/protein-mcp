"""Protein research MCP server. This file only registers the tools.

Tools:
  - resolve_protein: protein/gene name -> UniProt accession candidates.
  - search_literature: UniProt accession -> papers from Europe PMC
    (PubMed/MEDLINE records, PMC full-text records and, optionally, preprints).
  - get_structures: UniProt accession -> a few experimental PDB structures
    (RCSB) plus the AlphaFold DB model with its confidence.
  - get_ligands: UniProt accession -> ChEMBL ligands ranked by potency.
  - get_protein_profile: UniProt accession -> UniProt annotation + Open Targets
    associations and tractability (full result saved to the dossier on disk).
  - build_report: dossier on disk -> PDF + Markdown data report.
  - prepare_structures (Phase 2): structures dossier file -> original files on disk
    (mmCIF, AlphaFold, FASTA, ligand poses), cleaned receptors, chain checks,
    pocket and box per ligand + prep/manifest.json.
  - open_in_pymol (Phase 2): prep/manifest.json -> PyMOL session (.pml); starts
    PyMOL only if it is found, otherwise returns the path of the .pml.
  - select_structure (Phase 2): records the user's structure, chain(s) and pocket
    in prep/manifest.json, checked against the manifest.
  - pymol_view (Phase 2): changes the view of the PyMOL window opened with
    open_in_pymol (fixed actions only, never PyMOL text from the chat).
  - analyze_interfaces (Phase 2): partner chains, construct sequences, RCSB
    assemblies, PDBe PISA interfaces and protein-protein contacts.
  - dock_ligand (Phase 3): starts, in a detached process, a re-docking of the
    crystal ligand of a pocket with AutoDock Vina.
  - dock_status (Phase 3): the state of those runs and their result, the RMSD to
    the crystal pose calculated in place (calculated values); can stop a run.
  - build_experiment_report (Phase 3): prep/manifest.json -> separate PDF +
    Markdown record of the experiments (docking now; molecular dynamics and
    ADMET declared as not run), with the provenance of every run.

Logic lives in one module per source in protein_mcp/sources/: uniprot,
europepmc, rcsb, alphafold, chembl, opentargets. target_profile combines
uniprot + opentargets; dossier saves results to disk; report builds the PDF.
Shared helpers are in common. Phase 2 (files for in-silico experiments) is in protein_mcp/files/.

IMPORTANT: this server talks to Claude Desktop over stdio.
Never use print() anywhere: anything written to stdout corrupts the protocol.
Use logging (goes to stderr) instead.
"""
import asyncio
import logging
import re

import httpx
from mcp.server.fastmcp import FastMCP

logging.basicConfig(level=logging.INFO)

from protein_mcp import experiment_report, report, target_profile  # noqa: E402
from protein_mcp.common import parse_accession, safe  # noqa: E402
from protein_mcp.dossier import attach_dossier, dossier_name  # noqa: E402
from protein_mcp.files import (  # noqa: E402
    admet, choice, control, dock, interfaces, md_analyze, md_package, md_run, prepare, viewer,
)
from protein_mcp.sources import chembl, europepmc, uniprot  # noqa: E402
from protein_mcp.sources.alphafold import fetch_alphafold  # noqa: E402
from protein_mcp.sources.rcsb import MAX_PDB_MODELS, fetch_pdb  # noqa: E402

mcp = FastMCP("protein-research")


def _add_accession_warning(result: dict, check: dict) -> dict:
    """UniProt could not answer: keep the result and say so, in the chat and in the report's notes."""
    if check.get("warning"):
        result["accession_warning"] = check["warning"]
        result["notes"] = list(result.get("notes") or []) + [check["warning"]]
    return result


@mcp.tool()
async def resolve_protein(
    name: str, organism_id: int = 9606, max_results: int = 5
) -> dict:
    """Resolve a protein or gene name to UniProt accession candidates.

    Args:
        name: Protein or gene name, e.g. "KRAS", "BACE1", "insulin receptor".
        organism_id: NCBI taxonomy ID. Default 9606 (human). Use 10090 for mouse.
        max_results: How many candidates to return (1-10).

    Returns a list of candidates. If more than one plausible match exists,
    ask the user which one they mean before continuing. Never guess.
    """
    return await uniprot.resolve(name, organism_id, max_results)


@mcp.tool()
async def search_literature(
    uniprot_id: str,
    topic: str = "",
    max_results: int = 10,
    sort: str = "relevance",
    include_preprints: bool = False,
) -> dict:
    """Search the scientific literature (Europe PMC) for a protein.

    Args:
        uniprot_id: UniProt accession from resolve_protein, e.g. "P01116".
        topic: Optional extra words to narrow the search, e.g. "inhibitor",
            "crystal structure", "allosteric". Plain words must ALL appear in the
            title or abstract. Advanced Europe PMC syntax (with ":", "(" or
            AND/OR/NOT) is passed through unchanged.
        max_results: Papers to return (1-25). Keep it small: abstracts are long.
        sort: "relevance" (default), "recent" or "cited".
        include_preprints: Include bioRxiv/medRxiv etc. Default False. Preprints
            are NOT peer reviewed; they are flagged with is_preprint=true.

    Rules for using the result:
    - Base every statement ONLY on the returned titles and abstracts, and cite
      each one with the links provided. Never add facts that are not in them.
    - Only abstracts are returned, not full text. open_access/has_pdf only say
      whether a legal free full text exists.
    - If total_hits is much larger than returned, say the review is partial.
    The same result is also saved to `dossier_path` for build_report.
    """
    result = await europepmc.search(uniprot_id, topic, max_results, sort, include_preprints)
    accession, _ = parse_accession(uniprot_id)
    name = dossier_name("literature", topic, sort if sort != "relevance" else "",
                        "preprints" if include_preprints else "")
    params = {"topic": topic, "max_results": max_results, "sort": sort, "include_preprints": include_preprints}
    return attach_dossier(result, accession, name, "search_literature", params)


@mcp.tool()
async def get_structures(uniprot_id: str, max_models: int = 0, mutation: str = "") -> dict:
    """Get 3D structures for a protein: experimental PDB entries + AlphaFold model.

    Args:
        uniprot_id: UniProt accession from resolve_protein, e.g. "P01116".
        max_models: How many experimental PDB structures to return (1-25).
            0 (default) uses the MAX_PDB_MODELS constant in sources/rcsb.py.
        mutation: Optional text to restrict the PDB entries to one variant,
            e.g. "G12D". Use it when the user cares about a specific mutation;
            without it the picks come from the best-resolution entries overall,
            which may not include rare variants. RCSB finds candidates by
            full-text search; each one is then CONFIRMED to name the variant in
            its mutation field or title, and the rest are listed in
            `discarded_not_confirmed`. A bare position ("G12") matches any
            substitution there.

    The PDB structures are a SELECTION, not the full list: state "N of
    total_entries" and read `selection` to explain how they were chosen. With a
    mutation filter, total_entries counts RCSB text matches, not confirmed ones:
    report `confirmed_in_examined` and `discarded_count` too.
    The AlphaFold model is the canonical UniProt sequence unless an isoform
    model beats it by AF_ISOFORM_MARGIN pLDDT; `selection_reason` says why.

    Rules for using the result:
    - Report only what is in the result; never invent ligands, mutations or
      resolutions. Link every structure with its `links`.
    - Mutation data is depositor-written and unreliable. `mutations` is the
      official field, `mutations_from_title` are mutations named only in the
      title, `mutations_unparsed` is non-mutation text. An empty `mutations`
      does NOT prove wild type: check the title and say it is unconfirmed.
    - Resolution is not binding-site quality: check the ligand density before
      trusting a pose.
    - The AlphaFold model may be a non-canonical isoform (see `isoform`,
      `is_canonical_sequence`, `other_models`): always say which one was chosen.
    - AlphaFold gives one static prediction with no ligands, nucleotides or ions;
      do not use it to describe inhibitor-bound conformations.
    - The accession is checked in UniProt first: nonexistent, deleted or merged
      accessions return ok=false and save nothing. If merged, ask the user before
      using the accession named in the error. If `accession_warning` is present,
      UniProt could not answer: say the accession is unconfirmed and do not read
      an empty result as "no structures".
    The same result is also saved to `dossier_path` for build_report.
    """
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    check = await uniprot.check_accession(accession)
    if not check["ok"]:
        return {**check, "uniprot_id": accession}
    limit = max_models if max_models > 0 else MAX_PDB_MODELS
    limit = max(1, min(limit, 25))
    mutation = (mutation or "").strip()
    if mutation and not re.fullmatch(r"[A-Za-z0-9 _\-/.]{1,40}", mutation):
        return {
            "ok": False,
            "error": f"Invalid mutation filter '{mutation}': use letters, digits and - _ / . only (max 40 chars).",
        }

    async with httpx.AsyncClient(timeout=30.0, follow_redirects=True) as client:
        pdb, af = await asyncio.gather(
            safe(fetch_pdb(client, accession, limit, mutation), "RCSB PDB"),
            safe(fetch_alphafold(client, accession), "AlphaFold DB"),
        )

    result = {
        "ok": bool(pdb.get("ok") or af.get("ok")),
        "uniprot_id": accession,
        "experimental": pdb,
        "alphafold": af,
        "notes": [
            "PDB structures are a diverse selection, not the complete list.",
            "AlphaFold is a single static prediction without ligands, nucleotides or ions.",
            "All links are built from IDs/URLs returned by RCSB and AlphaFold DB.",
        ],
    }
    _add_accession_warning(result, check)
    return attach_dossier(result, accession, dossier_name("structures", mutation), "get_structures",
                          {"max_models": limit, "mutation": mutation})


@mcp.tool()
async def get_ligands(uniprot_id: str, mutation: str = "", max_ligands: int = 15):
    """Get ChEMBL bioactivity for a protein: ligands ranked by potency, with development status.

    Args:
        uniprot_id: UniProt accession from resolve_protein, e.g. "P01116".
        mutation: Optional variant, e.g. "G12D" ("G12" matches any substitution).
            Potency is reported per mutation; 'not stated' does NOT prove wild type.
        max_ligands: How many top ligands to return (1-50).

    This reply is a COMPACT view: per ligand it keeps the 3 best potency groups
    (`potency_groups_total` says how many exist) and per row one `source` link and
    `n_documents`. SMILES, properties, assay descriptions, every source document
    and the label-conflict details are saved to `dossier_path`.

    Report only what is in the result. Say "N of ligands_with_qualifying_data";
    if ligands_with_qualifying_data_is_lower_bound, say "at least N".
    Ranking is by best pChEMBL: potency is not selectivity or developability.
    ALWAYS mention the three mechanism lists, these are drugs too:
    - mechanism_compounds_without_qualifying_activity: no IC50/Ki/Kd (often
      antibodies); never say "no approved drugs" without checking their max_phase.
    - mechanism_compounds_other_variants_only: data exist only for other
      variants; say "no qualifying measurement on this variant", not "untested".
    - mechanism_compounds_unchecked: their activity download failed; unknown.
    If mechanisms.retrieved is false the lists are None: say drugs were NOT
    checked, never "none found". If activities.truncated, say that only the
    most potent activities (ranked_download_cutoff_pchembl) were examined, and
    unless activities.mechanism_compounds_download is "complete", that drug
    statistics may be partial.
    name is None for most research compounds: never invent names.
    Cite potency with its `source` link. If target is None, say "no record in ChEMBL".
    The accession is checked in UniProt first: nonexistent, deleted or merged
    accessions return ok=false and save nothing. If merged, ask the user before
    using the accession named in the error. If `accession_warning` is present,
    UniProt could not answer: say the accession is unconfirmed and do not read an
    empty result as "no ligands".
    """
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    check = await uniprot.check_accession(accession)
    if not check["ok"]:
        return {**check, "uniprot_id": accession}
    result = await chembl.get_ligands(accession, mutation, max_ligands)
    _add_accession_warning(result, check)
    return attach_dossier(result, accession, dossier_name("ligands", mutation), "get_ligands",
                          {"mutation": mutation, "max_ligands": max_ligands}, compact=chembl.compact)


@mcp.tool()
async def get_protein_profile(uniprot_id: str, max_diseases: int = 10) -> dict:
    """Get the target profile of a protein: UniProt annotation + Open Targets.

    Args:
        uniprot_id: UniProt accession from resolve_protein, e.g. "P01116".
        max_diseases: How many Open Targets disease associations to list in each
            of top_overall and top_oncology (1-50).

    This reply is a COMPACT summary. The full result (every variant with its
    description and PubMed IDs, disease descriptions, all examined associations
    with datatype scores, all tractability flags) is saved to `dossier_path`.

    Rules for using the result:
    - Keep the two sources separate: UniProt = curated annotation; Open Targets =
      association scores (0-1) computed from evidence, not proof of causality.
    - Say "N of total" for associations: they are the top `examined` of `total`
      by overall score. top_oncology keeps those in Open Targets' therapeutic
      area "cancer or benign tumor" (its classification; it includes some
      syndromes). Overall score can rank germline syndromes first (e.g. KRAS):
      always say which list you describe.
    - A UniProt disease with name None only has a note: quote the note, never
      give it a name.
    - Function evidence "by similarity" means inferred from another protein.
    - Tractability gives Open Targets' yes/no flags per modality; it is not a
      druggability verdict. Cite data_version.
    - When UniProt lists several Ensembl genes (e.g. alternate haplotypes), they
      are queried in order until one has an Open Targets target; the others are
      in ensembl_genes_without_target / other_ensembl_genes_not_queried.
    - If open_targets.found is false, say "no record in Open Targets". If
      uniprot_crosscheck is false, warn that Open Targets maps the gene to
      another protein.
    """
    full = await target_profile.get_protein_profile(uniprot_id)
    accession, _ = parse_accession(uniprot_id)
    return attach_dossier(full, accession, "profile", "get_protein_profile", {"max_diseases": max_diseases},
                          compact=lambda r: target_profile.compact(r, max_diseases))


@mcp.tool()
async def build_report(uniprot_id: str) -> dict:
    """Build the PDF + Markdown data report for a protein from its dossier on disk.

    Args:
        uniprot_id: UniProt accession, e.g. "P01116".

    The report uses ONLY the files that the other tools saved to the dossier
    (get_protein_profile, get_structures, get_ligands, search_literature; every
    call with different parameters adds its own file). Run the ones the user
    wants first. Sections without a file are listed in `missing` and left empty
    in the report, never filled from the conversation. If prepare_structures
    was run, section 6 "Structure preparation" shows prep/manifest.json:
    structures, chains, pockets and boxes, and the structure chosen with
    select_structure (or that none was chosen).

    Rules for using the result:
    - Give the user `pdf_path` (and `markdown_path`). If `missing` is not empty,
      say which sections are absent and offer to run those tools and rebuild.
    - `preparation` names the manifest used (None = no preparation section).
    - The report holds data and links only, no interpretation. If you add an
      interpretation in the chat, base it only on tool results and cite them.
    - If it fails because the PDF is open in a viewer, ask the user to close it.
    """
    return await asyncio.to_thread(report.build_report, uniprot_id)


@mcp.tool()
async def build_experiment_report(uniprot_id: str) -> dict:
    """Build the separate PDF + Markdown record of the in-silico experiments from prep/manifest.json (Phase 3).

    Args:
        uniprot_id: UniProt accession, e.g. "P15056".

    A document of its own, NOT a section of build_report: it holds the docking
    runs recorded in dossiers/<accession>/prep/manifest.json, with the
    provenance of each one (PDB entry, resolution, pocket, box, SHA-256 of the
    mmCIF, receptor chains and, where it applies, the literal reason recorded by
    select_structure), and the ADMET results recorded by analyze_admet. The
    molecular dynamics section is present and says "not run": no tool writes it
    yet.

    Rules for using the result:
    - Give the user `pdf_path` (and `markdown_path`), and name the runs in
      `docking_runs` and the compounds in `admet_compounds` as they come.
    - `sections` is the document's own inventory, read back from what was
      written: each entry has its `number`, `title` and `status`, which is
      either "has content" or "DECLARED HOLE (not run)" with the sentence the
      document prints. Answer "what does section N say" from it, and do not
      claim a section has content when its status says otherwise. It is an
      inventory, not the text: the content of a section is in the file, so to
      quote a section say that the file has to be opened.
    - The docking results (Vina scores, RMSD) are CALCULATED. The provenance
      quotes other kinds of value: the resolution is the experimental value of
      the PDB entry, and run times are measurements of this machine. Do not call
      a score an affinity and do not compare it with the measured IC50/Ki/Kd
      that build_report shows from ChEMBL: those stay in report_<accession>.pdf.
    - If `docking_runs` is empty, say so and offer dock_ligand + dock_status;
      never describe a run that is not listed.
    - If it fails because the PDF is open in a viewer, ask the user to close it.
    """
    return await asyncio.to_thread(experiment_report.build_experiment_report, uniprot_id)


@mcp.tool()
async def prepare_structures(
    uniprot_id: str, mutation: str = "", pdb_ids: list[str] | None = None, include_alphafold: bool = True
) -> dict:
    """Download the structures chosen by get_structures, clean them, find each ligand's pocket and box, and write the preparation manifest (Phase 2, steps 1-2).

    Args:
        uniprot_id: UniProt accession, e.g. "P15056".
        mutation: The same mutation used in get_structures, e.g. "V600E"; it picks
            the dossier file (structures_v600e.json). Empty = structures.json.
            Run get_structures with that mutation first.
        pdb_ids: Optional subset of the PDB IDs in that file, e.g. ["3OG7"]. IDs
            that are not in the file are NOT downloaded (`pdb_ids_not_in_dossier`).
        include_alphafold: Also download the AlphaFold model of the dossier file.

    Downloads into dossiers/<acc>/prep/: raw/<PDB>.cif (original mmCIF), the
    AlphaFold mmCIF, raw/<acc>.fasta (UniProt canonical sequence) and
    ligands/<PDB>_<comp>_<label chain>.sdf (crystal pose of every non-additive
    ligand instance). Then, with gemmi: clean/<PDB>_receptor.pdb (target chains
    only, first model and altloc, no hydrogens, waters, additives or ligands;
    nucleotides, cofactors and Mg/Mn/Zn/Fe/Co/Ni/Cu touching the target are
    kept), per-chain metrics (modelled and missing UniProt ranges, the residue
    actually modelled at the variant position), and for every ligand instance
    with >= 10 heavy atoms on a target chain: pocket residues within 5 Å and a
    box (ligand + 5 Å per side, mmCIF frame), plus AlphaFold pLDDT at the same
    residues. Writes prep/manifest.json (provenance with SHA-256, chains, ligand
    instances, construct mutations, all of the above, decisions, pending steps).
    Requested IDs that are in the file are prepared even when others are not.
    Running it again replaces the manifest (`replaced_manifest_from`), also with a
    subset of pdb_ids: then the manifest holds only that subset. If the RCSB
    chain query or every structure download fails, it returns ok false and keeps
    the previous manifest. A structure choice made with select_structure is kept
    only if its structure, chains and pocket are still in the new manifest with
    the same file: `user_choice` says "kept" (with `changed_fields`) or
    "dropped" (with the reason; tell the user it must be chosen again).
    Results of analyze_interfaces are kept per structure on the same rule
    (`interfaces`: kept / dropped).

    Rules for using the result:
    - Say "N of M" from `selection`, and list `failed_downloads` and `warnings`
      if any (e.g. RCSB or AlphaFold failed when the dossier was made: then the
      structures are unknown, not absent).
    - raw/ files are the ORIGINAL downloads; the receptor is only stripped as
      above. Nothing is protonated, completed (missing loops stay missing) or
      mutated: never say the structure is ready to simulate.
    - Per chain, report `chain_checks.variant`: "variant", "wild type" (the
      UniProt canonical residue), "other residue", "not modelled" or "ambiguous
      (microheterogeneity)" (two residues modelled at that position). A structure titled with the variant can have it
      unmodelled in the chain that holds the ligand (3OG7 chain A): say so.
    - Start from `summary`: facts across entries computed by the tool. Any
      statement about several entries ("all", "none", "the best", "only")
      must come from it, within its `scope`: every list except `by_resolution`
      covers only `structures_compared`. If `structures_not_compared` is not
      empty, say those entries were not checked; never turn an empty list
      into "none exists in the selection".
      - `by_resolution` (all entries, "(not compared)" marked) and
        `by_r_free` (compared only), sorted best first.
      - `chains_without_internal_gaps`: no gap inside the modelled span;
        the termini may still be unmodelled (the modelled range is given),
        so do not call such a chain complete.
      - `ligands_without_pocket`: ligand instances that got no pocket, why.
      - Only with a substitution filter such as "V600E" (not "G12"):
        `variant_position_by_chain_status` (each target chain of the
        compared entries by "variant", "not modelled", ...);
        `pocket_counts`, then `pockets_with_variant_in_pocket_chain` = the
        variant is modelled in at least one chain of the pocket (each line
        gives every pocket chain's status and residue count, and whether the
        pocket is single-chain; the variant residue need not be a pocket
        residue). Name every pocket of this list when asked where the variant
        and a ligand coexist;
        `pockets_with_ambiguous_variant_position`;
        `pockets_without_variant_in_pocket_chain`; and
        `pockets_not_classified` (pocket errors). If a list is not empty,
        never say that no such structure or chain exists.
    - Quote every number per structure and per chain, as given. Missing ranges,
      variant status, R-free and resolution differ between entries and between
      chains of one entry (8C7Y chain A has no missing range, 3OG7 chain A misses
      597-614): never generalise one entry's value to "all structures". To rank
      entries by a metric (resolution, r_free), sort the returned values; r_free
      None means not reported, not good.
    - If asked which structure to use, a recommendation is your judgement, not
      a result: compare all entries on the returned values, say which ones it
      rests on, and that the protein conformation (DFG, αC helix) and the fit
      to a different ligand are NOT measured by this tool.
    - For each pocket, give its chain, n_residues, box and any
      `adjacent_missing_ranges` (the pocket may be incomplete there). Flag
      `large_box`. `spans_target_chains` means ONE ligand instance has pocket
      residues in several target chains (`residues_by_chain`), e.g. a bivalent
      compound bridging two sites: its box covers all of them and is not a
      single-site box. Do not invent a reduced box; say that
      one would have to be defined. `alphafold_mean_plddt` is AlphaFold's confidence at those
      residues (`alphafold_residues_compared` of n_residues); AlphaFold has no
      ligand, and its residue at the variant position is in
      alphafold.variant_position. `pocket_not_computed` lists ligand instances
      without a pocket and why (too small: often additives). Pockets use the
      asymmetric unit only (no crystal-symmetry contacts).
    - If `warnings` has a "Variant check" line, the filter's wild-type letter
      does not match UniProt: say so before using chain_checks.
    - Name chains by auth_asym_id (what PyMOL shows) and say when it differs from
      label_asym_id. Target chains come from RCSB's UniProt mapping.
    - Identify ligands only by `name`, `synonyms` (wwPDB-curated drug or code
      names), `formula_weight` and the entry `title` of THIS result. If none
      names the compound, say it is unnamed; do not identify it from memory or
      from the examples in this text.
      A comp_id is not a ChEMBL ID: no link to get_ligands is made.
    - Every failed download is in `failed_downloads` (the manifest has the
      same list); transient errors were retried first, "(after N attempts)".
      Name each missing file.
    - Say which chain(s) each ligand touches (`neighbor_target_auth_chains`): a
      ligand can be bound to one chain only. An empty list with
      touches_other_chains true means it sits on a partner protein only: do not
      call it a target ligand. Empty with touches_other_chains false means RCSB
      gave no neighbours: say it is unknown.
    - `declared_mutations` lists the depositor's mutations of the target chain,
      including engineered ones: if it is not empty, the construct is not the
      native sequence. Empty does NOT prove a native sequence (the variant may be
      named only in the title, or `mutations_unparsed` may hold text like "YES"):
      say it is unconfirmed.
    - `sample_uniprot_ranges` is the UniProt range of the crystallised/sampled
      construct, not the residues resolved in the model (7ZR0: 1-766 sampled).
    """
    return await prepare.prepare_structures(uniprot_id, mutation, pdb_ids, include_alphafold)


@mcp.tool()
async def open_in_pymol(
    uniprot_id: str, pdb_ids: list[str] | None = None, include_alphafold: bool = True, start: bool = True
) -> dict:
    """Show the prepared structures in PyMOL: write a session file (.pml) from prep/manifest.json and start PyMOL if it is installed (Phase 2, step 3). Use it ONLY when the user asks to see or open the structures in PyMOL.

    Args:
        uniprot_id: UniProt accession, e.g. "P15056". Run prepare_structures first.
        pdb_ids: Optional subset of the structures in the manifest, e.g. ["3OG7"].
            Empty = all of them. IDs not in the manifest are listed in
            `pdb_ids_not_in_manifest`.
        include_alphafold: Also load the AlphaFold model, coloured by pLDDT.
        start: False = only write the .pml (when the user wants the file).

    The session uses only files already in dossiers/<acc>/prep/ (no network):
    each cleaned receptor, the crystal ligand poses (SDF), one selection per
    pocket of the manifest, the box of each pocket (lines), the residue modelled
    at the variant position, and the AlphaFold model superposed on the first
    structure (in the session only). Only the first structure is shown; the
    others are hidden in the object panel. The file is
    prep/pymol/session.pml and is replaced on every call.

    Rules for using the result:
    - If `pymol.launched` is true, say PyMOL was started (the server cannot see
      the window). List every `pymol.session_errors` line (a file that failed
      to load); if `session_loaded` is false, say loading was still running at
      the check and errors after it are only in `pymol.log_path`. If launched
      is false, give `pml_path` and `how_to_open`, plus `pymol.error` and
      `pymol.hint` if present (e.g. PyMOL without Qt).
    - Explain what the user will see with `legend` and the names in
      `structures` (group, receptor, ligands, pockets, variant_position) and
      `alphafold`; say which structure is shown first (`enabled`). Mention
      `skipped`, `pdb_ids_not_in_manifest` and `notes` (frames: structures are
      not superposed with each other).
    - Name what is NOT shown: per structure `ligands_not_loaded` (failed or
      missing SDF), `pockets_not_shown` and, per pocket,
      `residues_not_selectable` (then `n_residues_selected` < `n_residues`).
    - A ligand with `on_target_chain` false is drawn but touches no target
      chain (e.g. ATP on a partner protein that is not in the receptor): never
      call it bound to the target. Name the protein it sits on only from its
      `partner_chains` ("auth chain: RCSB entity description"), never from the
      entry title. `has_pocket` false means the manifest has
      no pocket for it (often a small additive-like compound).
    - The magenta residue (`variant_position`) is whatever is modelled at the
      variant position: give each chain's `observed` and `status` from
      `by_chain` ("variant", "wild type", ...); do not call it the mutation
      unless the status is "variant".
    - The view shows the manifest's data; it measures nothing new. Numbers
      (pocket residues, pLDDT, variant status) come from prepare_structures,
      not from this tool.
    - The session file also starts the chat control: while this window is open,
      use pymol_view to change the view (`control.session_id`). If
      `pymol.control_listening` is false after a launch, say the chat control
      may not be available.
    """
    return await viewer.open_in_pymol(uniprot_id, pdb_ids, include_alphafold, start)


@mcp.tool()
async def select_structure(
    uniprot_id: str, structure: str, reason: str, pocket: str = "", chains: list[str] | None = None
) -> dict:
    """Record the USER's choice of structure (and pocket) for the next phase in prep/manifest.json (Phase 2, step 4). Call it ONLY when the user has explicitly chosen; never to record your own recommendation without the user confirming it.

    Args:
        uniprot_id: UniProt accession, e.g. "P15056". Run prepare_structures first.
        structure: A PDB ID of the manifest, e.g. "3OG7", or "AlphaFold".
        reason: The user's own reason, in their words (quote them). If the user
            gave none, ask for it; do not write one yourself.
        pocket: Optional ligand pocket as "<comp_id> <label_asym_id>", e.g.
            "032 C", exactly as in prepare_structures (`pockets`, or the summary
            lines "8C7Y TXV U: ..."). "032" alone works if it is unique in that
            structure. AlphaFold has no pocket.
        chains: Optional target chains by auth_asym_id (what PyMOL shows), e.g.
            ["A"]. With a pocket, the default is the pocket's chains; without a
            pocket, required if the structure has several target chains.

    Everything is checked against the manifest: a structure, chain or pocket
    that is not there, or a structure without a cleaned receptor, is rejected
    with `options` (what can be chosen) and nothing is written. The choice
    stores a copy of the manifest data (receptor file, chains with modelled
    ranges and variant status, pocket residues per chain, box, SDF, construct
    mutations), the reason and the time. Choosing again replaces it; earlier
    choices stay in the manifest's user_choice_history.

    Rules for using the result:
    - The choice is checked before the reason: if only the reason is missing,
      `valid_choice` names the choice; ask the user for their reason.
    - If ok is false, say why and offer the `options`; do not retry with a
      structure, pocket or chain the user did not name (e.g. "has 2 target
      chains": ask which one, do not pick it).
    - Confirm with `recorded` and list every `user_choice.warnings` line (e.g.
      no pocket = no box; a pocket spanning several chains keeps its full box;
      a missing SDF).
    - Say that nothing was recomputed and the receptor is not protonated or
      completed: it is a recorded choice, not a ready-to-simulate system.
    - If `replaced_choice` is set, name the choice it replaced.
    """
    return await asyncio.to_thread(choice.select_structure, uniprot_id, structure, pocket, chains, reason)


@mcp.tool()
async def pymol_view(
    uniprot_id: str, action: str, structure: str = "", pocket: str = "", chain: str = "",
    residues: str = "", color: str = "", target: str = ""
) -> dict:
    """Change the view of the PyMOL window opened with open_in_pymol, from the chat (Phase 2, step 4b). Use it when the user asks to show, hide, colour or focus something in the open PyMOL window.

    Args:
        uniprot_id: UniProt accession, e.g. "P15056".
        action: One of:
            "status" (what is loaded and shown),
            "show_structure" (show `structure`, hide the others),
            "focus_pocket" (show `structure` and zoom on `pocket`, with its box),
            "color_chain" (colour `chain` of `structure` with `color`),
            "show_residues" (sticks + `color`, default green, and zoom on
                `residues` of `structure`/`chain`),
            "hide" / "show" (`target`: "ligands", "boxes" or "alphafold"; for
                ligands and boxes, only of `structure` if given),
            "reset" (draw the open_in_pymol scene again).
        structure: A PDB ID loaded in the session, or "AlphaFold".
        pocket: "<comp_id> <label_asym_id>", e.g. "032 C" or "TXV U" ("032"
            alone if unique), as in prepare_structures.
        chain: Target chain by auth_asym_id (what PyMOL shows), e.g. "B".
            Needed when the structure has several target chains.
        residues: UniProt numbers, e.g. "594-596, 600". They are translated to
            the structure's own numbering; residues not modelled in that chain
            are listed in `not_modelled`, never guessed. On "AlphaFold" only
            for a single canonical model covering the whole sequence (F1).
        color: red, green, blue, yellow, orange, magenta, cyan, white, grey,
            pink, purple, salmon, wheat or slate.
        target: For hide/show: "ligands", "boxes" or "alphafold".

    Only these actions exist: the server builds every PyMOL command from the
    manifest and the session, and never runs PyMOL text from the chat. It does
    not compute anything about the structure, and changes are not saved
    (reopening the session restores the open_in_pymol scene).

    Rules for using the result:
    - The server cannot see the window: say the command was sent and PyMOL
      accepted it (`applied` of `sent`), not what the picture looks like. List
      every `failed` line.
    - For status, report per structure `shown`, `ligands_hidden` and
      `boxes_hidden`, and `alphafold.shown`.
    - If ok is false with `pymol_open` false, PyMOL is closed or still loading
      (the error says which is likely): offer to try again in a few seconds or
      to open it again with open_in_pymol (that restores the initial view). If
      the error says the commands were sent without an answer, do not say
      PyMOL is closed: they may have been applied; ask the user. Otherwise say
      the error and offer the `options`; do not retry with a structure,
      chain, pocket, residue or colour the user did not name.
    - For show_residues, give `residues` (UniProt number: residue, author
      number) and name every `not_modelled` number: those are absent from the
      model, not hidden. Mention `ambiguous` (two residues at one position).
    - If the user asks for something not in the actions (a label, a surface,
      measuring a distance, a new file), say it is not available from the chat
      and how to do it in the PyMOL menus.
    """
    return await asyncio.to_thread(control.pymol_view, uniprot_id, action, structure, pocket, chain,
                                   residues, color, target)


@mcp.tool()
async def analyze_interfaces(uniprot_id: str, pdb_ids: list[str] | None = None) -> dict:
    """Partner chains, their sequences, the assemblies and the protein-protein contacts of the prepared structures (Phase 2, step 5). Use it for targets of biologics (antibodies, receptors, peptides) or when the user asks what touches the target protein. Run prepare_structures first.

    Args:
        uniprot_id: UniProt accession, e.g. "P01375".
        pdb_ids: Optional subset of the structures in the manifest, e.g. ["3WD5"].

    For each structure (from its raw mmCIF and two sources) it writes an
    `interfaces` block in prep/manifest.json:
    - sequences: every polymer entity (target and partners) with description,
      UniProt (if any), chains and the CONSTRUCT sequence with flags (runs of
      6+ histidines, often a tag but not stated by RCSB; non-standard
      residues X; other tags are not detected);
    - rcsb_assemblies: RCSB's biological assembly, who defined it (author or
      software) and whether it is just the asymmetric unit;
    - pisa: PDBe PISA's interfaces (area, css) and multimer sets with PISA's
      verdict; each PISA service can fail alone: that is said and the rest
      goes on;
    - contacts: target residues within 5 Å of another polymer chain, also of
      crystal-symmetry copies, per partner chain and symmetry operator, in
      UniProt numbering, with the matching PISA interface if any.
    Nothing is built or predicted: no oligomer, no epitope, no antibody label.

    Rules for using the result:
    - Say "N of M" from `selection`; list `skipped`, `warnings` and `pisa_errors`.
    - A contact line ends with its PISA interface (same chain pair and
      operator), "PISA lists no interface for this chain pair and operator",
      or "PISA not available" (then nothing is known from PISA about it).
    - `pisa_multimer_sets`: each set divides all chains of the crystal into
      assemblies that coexist (joined by "+"); the sets are PISA's
      alternative solutions. Do not mix assemblies of different sets.
    - Never call a contact or assembly "biological" or "crystal packing" on
      your own. Quote the sources: RCSB's assembly with its `defined_by`, PISA's
      verdict for each multimer and the css of each interface. If they
      disagree (e.g. RCSB's assembly is the asymmetric unit while PISA calls a
      trimer of the target stable), say both. A contact at a symmetry copy
      ("symmetry copy", code other than 1_555) exists in the crystal; whether
      it exists in solution is not measured here.
    - Name partners only from `partners` (RCSB entity description); do not
      call a chain an antibody, Fab or receptor unless its description says
      so. Partner sequences are crystallised constructs: mention tags and X
      residues before using them for a complex prediction.
    - Give contact residues with their UniProt numbers as returned; "(auth)"
      means no UniProt number. A partner chain can touch the target in its own
      copy and in a neighbouring copy: report each line separately.
    """
    return await interfaces.analyze_interfaces(uniprot_id, pdb_ids)


@mcp.tool()
async def dock_ligand(
    uniprot_id: str,
    structure: str = "",
    pocket: str = "",
    conformer: str = "regenerated",
    exhaustiveness: int = 32,
    seed: int = 42,
) -> dict:
    """Phase 3: START a re-docking of the crystal ligand of one pocket into its own receptor with
    AutoDock Vina, in the background. Answers in seconds; the result comes from dock_status.

    Docking takes minutes (Vina alone: 61-124 s measured for vemurafenib in 3OG7 with
    exhaustiveness 32) and a chat request is cut long before that (60 s or 240 s, measured in
    Claude Desktop), so this tool validates everything, launches a separate process and returns.
    The run is not cut when this request or this chat ends; whether it also survives Claude
    Desktop closing is not known, so do not promise that. NEVER promise a result in this reply
    and never say the docking finished: `launched: true` only means it started. Tell the user it is running, that one run at a time is allowed on this machine, and
    that you will read the result with dock_status (wait at least a minute before the first call).

    Reads dossiers/<acc>/prep/manifest.json: run prepare_structures first. Receptor, box and
    ligand all come from the manifest; nothing is taken from the chat. With no `structure`, the
    user's choice (select_structure) is used and must have a pocket. `pocket` as in
    select_structure: "032 C", "032" if unique, or "3OG7 032 C". One ligand per call.
    - conformer "regenerated" (default): a new 3D conformer built from the crystal ligand's own
      atoms and bonds, moved to the box centre. With exhaustiveness 32 and seed 42 this is the
      pre-registered control (`is_preregistered_control: true`); other values are reported
      against the same threshold but are not the control. "crystal" starts from the crystal
      conformation: biased in favour; say so whenever you report it.
    - The receptor is the whole cleaned receptor (all target chains), even if the user's choice
      names fewer chains: report `receptor_preparation.chains_note` when present.
    - exhaustiveness (1-64) and seed are recorded; num_modes 9, energy_range 3, min_rmsd 1 and
      scoring "vina" are fixed.

    Refused before anything is launched (`ok: false`, `launched: false`, with `options` when the
    pocket is not found): pockets marked large_box or spans_target_chains, receptors that kept a
    nucleotide/cofactor/metal (named), changed files, Vina not found (PROTEIN_MCP_VINA), and a
    docking run already in progress on this machine (the reply names it; offer dock_status, and
    do not launch anything else until it ends). The checks that need Meeko and Vina (the box
    guard, residues without a template) happen inside the run and appear in dock_status.

    Say only what this reply contains: the pocket, the box, the receptor chains (quote
    `receptor.chains_note` when present), the ligand by `ligand.name`/`ligand.synonyms` only, and
    that it is running. Do not call dock_ligand again for the same pocket while it is live.
    """
    return await dock.start_dock(uniprot_id, structure, pocket, conformer, exhaustiveness, seed)


@mcp.tool()
async def dock_status(
    uniprot_id: str,
    structure: str = "",
    pocket: str = "",
    conformer: str = "",
    stop: bool = False,
) -> dict:
    """Phase 3: read what the background docking runs of a protein are doing, and the result of
    the ones that finished. With stop=True, stop the run that is live.

    Each run reports a `state`: "running" (or "queued") with its `phase` and `elapsed_seconds`;
    "finished" with the full result; "failed" or "stopped" with the reason in `error`;
    "vanished", which means its process stopped reporting (the computer was shut down or it was
    killed); or "unknown", which means its state file could not be read and nothing can be said
    about it. A "running" run is not a result: say it is still running, give the elapsed time, and
    offer to look again; never guess an RMSD. Use stop=True only when the user asks to stop it.

    Two fields must be quoted whenever they are present, because they change what the numbers
    mean:
    - `warning` on a run: it has not reported for a while, so `phase` and `elapsed_seconds` are
      from then and not from now, and the run may already be dead. Never present them as current.
    - `note_on_state`: the result was recovered from the manifest after the run stopped reporting,
      or the run was not started by dock_ligand (the command line, or an older version). Say so;
      do not present it as a run this chat started and watched.
    Report every run by its own structure, pocket and conformer, and `selection` as given.

    Rules for using a finished result:
    - A Vina score (kcal/mol) is CALCULATED, not a measured affinity. Never call it an affinity
      or a binding energy, never put it next to ChEMBL IC50/Ki/Kd as if comparable, and never
      rank compounds by it here.
    - Report `redocking.top_pose_rmsd_angstrom` with the criterion and `within_criterion`, and
      quote `scope`: it is one re-docking of one ligand in one structure, not the accuracy of
      docking for this target. If `is_preregistered_control` is false, quote `redocking.note`:
      "biased" only for conformer "crystal"; other exhaustiveness/seed is "not the
      pre-registered control", not biased.
    - Name every residue in `receptor_preparation.deleted_residues` (Meeko removed them because
      they match no template and lie outside the box; auth numbering). If
      `receptor_preparation.warning` is present, quote it: the list may be incomplete.
    - Give `poses_returned` as is: fewer than 9 poses means the others were outside energy_range.
    - Identify the ligand only by `ligand.name` and `ligand.synonyms` (wwPDB); without a
      synonym, do not give it a drug name.
    - The receptor is rigid and was not completed: do not call the result a binding mode
      prediction for a new compound.
    - `vina.cpu_seconds` and `vina.cpu_per_wall_second` are measurements of this machine, not of
      the docking: report them only if asked about speed, and if `vina.cpu_seconds_note` is
      present say nothing about speed at all (it could not be measured; it is not a slow run).
    """
    return await dock.dock_status(uniprot_id, structure, pocket, conformer, stop)


@mcp.tool()
async def build_md_package(
    uniprot_id: str,
    structure: str = "",
    pocket: str = "",
    conformer: str = "",
    pose: int = 1,
) -> dict:
    """Phase 3: write the folder the USER uploads to CHARMM-GUI to build a molecular dynamics
    system for a docked complex (receptor + one pose, in solution).

    CHARMM-GUI has no public API for this, so the delivery is manual: this tool only writes
    files and checks them. It does not contact CHARMM-GUI and starts no simulation. Run
    dock_ligand first; the pose comes from the manifest, never from the chat.

    `structure` and `pocket` as in dock_ligand (with no `structure`, the user's choice is used).
    `conformer` names the run when the pocket has more than one. `pose` is 1-based over the
    poses of that run, 1 (best scored) by default.

    Tell the user, from this reply and nothing else:
    - the two files of `upload` and the folder `package_dir`, and that the PDB must keep its
      name (CHARMM-GUI reads the entry id from the file name);
    - which runs to tick under "Model missing residues" (`tick_model_missing_residues`) and
      which to leave unticked (`leave_unticked`: disordered terminal tails, not gaps);
    - that README.md in that folder has the full click-by-click sequence;
    - that when CHARMM-GUI returns the .tgz they should come back, because
      check_charmm_gui_output has to verify it before anything is simulated. Of the four jobs of
      the validation two finished with a broken system, so finishing the wizard proves nothing.

    Always repeat `modelled_not_experimental`: the residues CHARMM-GUI builds into the gaps are
    modelled, not experimental, and any later report must say so where they are discussed. Give
    `missing_residues_source` as it comes: when the mmCIF has no record of unobserved residues,
    say that the tool cannot tell, never that the structure has no gaps. If
    `internal_runs_partly_in_the_file` is present, pass it on: those runs can be neither modelled
    nor ignored, and the user has to look at the entry.

    Refused (`ok: false`): no manifest, no docking run for that pocket, several runs without
    `conformer`, a pose index outside the range, files changed since the manifest, or a written
    PDB and SDF that do not agree atom by atom.
    """
    return md_package.build_md_package(uniprot_id, structure, pocket, conformer, pose)


@mcp.tool()
async def check_charmm_gui_output(
    uniprot_id: str,
    archive: str,
    structure: str = "",
    pocket: str = "",
) -> dict:
    """Phase 3: verify the .tgz CHARMM-GUI returns against the package build_md_package handed
    over, BEFORE any simulation is submitted.

    `archive` is the path to the .tgz on this machine. Nothing is sent anywhere; the archive is
    read into a temporary folder and deleted.

    Read the answer like this:
    - `passes: true` means every check that COULD run passed. It is not "everything is fine":
      always read `not_checked` and repeat it to the user. A check that could not run is a real
      gap, not a detail (without the Modeller intermediate the ANISOU trap that broke one of the
      validation jobs is not looked for at all).
    - `passes: false` means it must NOT be submitted: list `failures` verbatim. They are
      diagnoses, not guesses: a system without the ligand, a protein that came back distorted,
      gaps bridged by a long artificial bond instead of filled, a pose that changed shape, or an
      archive from another job.
    - `checks` holds the measured values (geometry, ligand, box, ions, temperature, steps,
      CGenFF penalties). Quote numbers from there, never from memory.
    - Repeat `modelled_not_experimental` whenever the filled gaps are discussed.

    Say plainly that reaching the end of the CHARMM-GUI wizard proves nothing: it reports no
    error for any of these failures. What is NOT measured: the absolute displacement of the pose
    (shape, pocket contacts and clashes are, so a rigid shift of about 1.5 A or less would not be
    caught) and the values under `system`, which are reported without any threshold.
    """
    return md_package.check_charmm_gui_output(uniprot_id, archive, structure, pocket)


@mcp.tool()
async def build_md_job(
    uniprot_id: str,
    archive: str,
    structure: str = "",
    pocket: str = "",
    gpu_resident: bool = False,
    site_report: str = "",
    partition: str = "",
    account: str = "",
    time_limit: str = "",
    cpus_per_task: int = 0,
    gpus: int = 1,
    gpu_type: str = "",
    namd_module: str = "",
    namd_binary: str = "",
) -> dict:
    """Phase 3: write the SLURM job the USER submits on THEIR cluster to run the molecular
    dynamics of a CHARMM-GUI system.

    `archive` is the path to the .tgz on this machine. The archive is re-checked here with the
    same measurements as check_charmm_gui_output and the job is REFUSED if they fail, so a
    `ok: false` with `failures` means the system must not be simulated, not that the tool broke.
    Nothing is sent anywhere and no simulation is started: the transfer and `sbatch` are manual,
    exactly like the CHARMM-GUI upload. `structure` and `pocket` as in build_md_package.

    THIS TOOL KNOWS NOTHING ABOUT ANY CLUSTER, and will not guess. Called without `site_report`
    it writes no script and comes back with `site_unknown` plus the commands in
    `paste_the_output_of`. Give those commands to the user to run on their own machine and ask
    them to paste everything back; then call again with all of it in `site_report`. Never fill a
    partition, account, module, binary or time limit from memory, from another session, or from
    what some other cluster uses, and never invent the output of those commands yourself.

    The reply then separates two things that must not be mixed, and you should keep them apart
    when you talk about them:
    - `simulation`: read from the archive. The chunk count is `cntmax` from its own namd/README
      and the chunk length comes from the `run` line of the production input. The equilibration
      runs as CHARMM-GUI wrote it.
    - `site`: one entry per `#SBATCH` field, each with a `mark` and a `source`. REPORTED means it
      was read from the pasted output and the source names the command; CHOSEN means the user
      said so; CONFIRM means nothing could fill it and the script will abort until they do;
      EDITABLE is a knob whose default comes from a stated rule. Quote the `source` when the
      user asks where a value came from, and say plainly when something is EDITABLE rather than
      measured. `site_report_unparsed` holds the pasted lines that were not understood: pass them
      on instead of pretending they were used.

    Tell the user, from this reply and nothing else: `run_dir`; that README.md there has the
    transfer and submit lines; whatever is in `to_confirm`; `run_several_at_once` when present
    (one SLURM array task per system, rebuilt from the manifest with every system whose archive
    passed the checks); `verified.not_checked`, repeated as it comes, because a check that could
    not run is a real gap; and `modelled_not_experimental` wherever the filled gaps are discussed.

    For a longer run, give both places in `edit_to_change_the_length`: CHUNKS in run_md.slurm, or
    the chunk length in the production input. For that second one give BOTH lines it names, with
    their numbers and current text: the real input carries `numsteps` AND `run`, NAMD executes the
    steps of `run`, and raising only `numsteps` would leave every chunk at its current length.
    Those are STEPS, not seconds, and editing the input means the system is no longer the one
    CHARMM-GUI returned and verified.

    Do NOT estimate how long the run will take: no wall time per nanosecond has been measured for
    this system on any machine. `--time` is the limit the cluster REPORTED, used as a ceiling; the
    honest way to replace it is to measure the first chunk with the `sacct` line the reply and the
    README give. The script times each production chunk and keeps the last measurement, so a
    re-submitted job is guarded from its first chunk and continues where it stopped; the very
    first chunk of the first submission runs unguarded.

    What the cluster reports is not a promise that the job runs: whether NAMD works there, and
    whether the site wants `--gres=gpu` or `--gpus`, only the first `sbatch` will tell. Say so.

    A partition can be several kinds of node: on a real cluster `gpu` had six rows with five GPU
    models (rtx_2080, rtx_3090, a100, V100). When that happens nothing is sized from them and
    `site.partition.gpu_types_reported` lists the models: tell the user, because without
    `gpu_type` the job lands on whichever node the scheduler gives it, and for molecular dynamics
    an rtx_2080 and an a100 are not the same machine. With `gpu_type` the request is sized against
    that node and the model is pinned in the GRES.

    The overrides (`partition`, `account`, `time_limit`, `cpus_per_task`, `gpus`, `gpu_type`,
    `namd_module`, `namd_binary`) are for values the USER states. A partition that the pasted output does not
    list, or a time over the limit it reports, is refused rather than written.

    `gpu_resident: true` inserts `CUDASOAintegrate on` before the line that executes the run. It
    is an edit to CHARMM-GUI's output, reported in `inputs_edited`: only use it if the user asks,
    and then say the input was changed. It has not been tested on any cluster, and the
    equilibration is deliberately left untouched.
    """
    return md_run.build_md_job(uniprot_id, archive, structure, pocket, gpu_resident, site_report,
                               partition, account, time_limit, cpus_per_task, gpus, gpu_type,
                               namd_module, namd_binary)


@mcp.tool()
async def analyze_md(
    folder: str,
    topology: str = "",
    ligand_resname: str = "LIG",
    ligand_segid: str = "",
    pocket_cutoff: float = 5.0,
    skip_frames: int = -1,
    reference_frame: int = -1,
    include_series: bool = False,
    stage: str = "",
) -> dict:
    """Phase 3: measure a FINISHED molecular dynamics trajectory that the user has brought back
    from their cluster, and report five metrics with the selection and the criterion behind each.

    `folder` is a path ON THIS MACHINE holding the PSF and the DCD files (the folder itself, or
    its `system/namd/`, which is where a build_md_job package keeps them). Nothing is downloaded:
    the user copies the files back themselves, with scp, rsync or a client such as WinSCP. If the
    folder has no .dcd yet, the run has not written one and the reply says so; do not retry.

    What each metric answers, and in this order, because the first one gates the rest:
    - `protein_rmsd`: the backbone against the reference frame, ONE SERIES PER PROTEIN SEGMENT
      (`protein_rmsd.per_segment`); the headline numbers are the segment with the LARGEST DRIFT,
      named in `protein_rmsd.segid` and with the rule in `protein_rmsd.headline_rule`. That is not
      the segment that moved most: a chain swinging widely about one position has a drift near
      zero, so check `per_segment` ranges before saying anything about how much the protein moved.
      This is the EQUILIBRATION CONTROL. If it has not reached a plateau the
      protein is still settling and every other number is premature; say that before discussing
      the ligand. Chains are never fitted together, so motion of one chain relative to another is
      NOT measured here: do not describe the whole protein from one segment's number.
      BEFORE quoting it, read `trajectory.restraints`. The tool READS the NAMD restraint reference
      file and reports, per group, how many of the atoms it measured were held: `restraints.held.
      protein_backbone` and `restraints.held.ligand_heavy`. A CHARMM-GUI equilibration holds the
      ligand's heavy atoms too, so BOTH the plateau and the ligand figures can be meaningless on
      the same trajectory. Where `held` is above zero, say the figure passes by construction and
      is not evidence; the offered reading marks those lines itself.
      Each group carries `known`. When it is FALSE there is no number, and you report it as
      UNKNOWN rather than as "not held": it means the reference file was absent, unreadable, named
      a force-constant column that is not X/Y/Z/O/B, described a different number of atoms, or had
      the right count in another order (its atom and residue names are checked against the
      topology at every index this tool measures). Read `restraints.force_constant_column_source`
      too: when no column was declared the tool ASSUMES B, which is what CHARMM-GUI writes, and
      the counts are then unverified - say so. `restraints.declarations` lists every hold the
      input declares; with both a spring and fixed atoms, the counts come from the frozen one.
      And read `protein_rmsd.drift_is_meaningful`: with few analysed frames the drift is a
      difference between single frames and cannot show a plateau at all. When it is false, pass on
      `drift_warning` and do not call anything equilibrated.
    - `ligand_rmsd`: the ligand's heavy atoms, superposed on the POCKET backbone. Report that it
      is fitted on the pocket and not on the ligand, because the second would show a ligand
      sitting still while it leaves the site.
    - `com_distance`: ligand centre of mass against the pocket's. This is what separates a ligand
      that reorients inside the site from one that dissociates. `pbc_corrected_frames` says
      whether any frame needed a periodic-image correction.
    - `hbonds`: per donor-H-acceptor pair, `frames_present` out of `frames_analysed`.
      `occupancy_sampled` is a SHARE OF FRAMES, never a share of elapsed time: hydrogen bonds
      break and re-form faster than a DCD stores frames. Repeat that whenever you quote it.
    - `rmsf`: per pocket residue. With no apo control, a flexible residue cannot be attributed to
      the ligand; say so rather than implying the ligand caused it.

    Always report, from the reply and nothing else: `trajectory.frames_note` ("N of M frames
    analysed"), `trajectory.minimization_note`, `trajectory.restraints.note`,
    `trajectory.span_ns`, `ligand.selection` with
    `ligand.dropped_note`, and `pocket.selection` with `pocket.segments_note`. Name the pocket
    residues from `pocket.residues` when they are discussed; they are measured by distance in the
    reference frame, so they are evidence, and a dossier label such as `3OG7_032_C` names a
    POCKET, not a segment of the simulated system.

    ONE STAGE is analysed, never two. `trajectory.stage` is the CHARMM-GUI stage (`step5`,
    `step4`, ...), taken from the trajectory file names, and by default it is the LATEST stage
    present, because CHARMM-GUI numbers stages in protocol order. Every other trajectory is named
    in `trajectory.left_out` with its own stage, and `trajectory.stages_present` lists them all:
    they were simulated at another point of the protocol, usually at another `dcdfreq`, so
    chaining them would run the time axis backwards. The `stage` argument analyses a different one
    instead, and naming a stage that is not there is a refusal that lists the ones that are.

    THE TOOL NEVER DECIDES THAT A TRAJECTORY IS "PRODUCTION" OR "EQUILIBRATION", and neither may
    you from its file name. Those words come only from the NAMD input of that stage, named in
    `trajectory.stage_inputs`, which whoever built the job wrote; say where the word comes from.
    What the input actually declares is in `trajectory.restraints` (were atoms held in place?) and
    `trajectory.minimization` (are the first frames minimization rather than dynamics?), and those
    two are the facts that decide whether a number means anything. When `stage_inputs` is empty, no
    input was found and the part of the protocol is UNKNOWN: report it as unknown.
    Always report `trajectory.stage`, `trajectory.stage_selection` and `trajectory.stage_note`.

    Refusals are results, not failures, and come back as `ok: false` with the reason in `error`.
    Trajectories whose frame spacing differs are refused by name instead of averaged, which
    happens when `dcdfreq` was changed between chunks: pass on which files differ. A reference
    frame with no residue near the ligand is refused instead of producing a residue list nothing
    supports. Never retry a refusal unchanged, and never present one as a measurement.

    The optional arguments are for the cases the defaults do not cover, and all of them are the
    USER's choice, never yours to guess: `topology` names the PSF when it is not in the folder or
    not the one the search would pick; `ligand_resname` defaults to `LIG`, which is what
    CHARMM-GUI names the ligand it parameterised, and a refusal lists the non-standard residue
    names actually present, so pass one of THOSE rather than inventing one; `ligand_segid` picks
    ONE copy when the system holds several, which is refused otherwise - two copies measured
    together put the centre of mass between them, where there is no ligand - and the refusal lists
    the copies, so pass one of those segments and run the others separately, reporting each on its
    own; `pocket_cutoff`
    (angstrom, default 5.0) changes how far from the ligand a residue still counts as pocket, and
    changing it changes the residue list, the superposition and the centre-of-mass reference
    together, so say that it was changed whenever it was.

    `com_distance.pbc_corrected_frames` is not a warning about the ligand. CHARMM-GUI writes
    `wrapAll on`, so NAMD re-images whole molecules back into the box; each frame's ligand is put
    in the image nearest where it was in the PREVIOUS frame, and the same corrected coordinates
    feed the RMSD. Report the count when it is above zero, as re-imaging and not as motion -
    but ONLY as a count of corrections, never as a reason to discount a large RMSD: a ligand that
    really leaves the site produces a large RMSD and no correction at all, and the two are
    different facts. Three fields say when the numbers above are not trustworthy and each must be
    passed on plainly when it is above zero: `non_orthorhombic_frames` (the correction rule
    assumes 90 deg angles and was NOT applied), `frames_without_box` (the frame declares no unit
    cell, so nothing could be corrected) and `pocket_split_frames` (the pocket's own atoms span
    more than half a box, so its centre of mass and the superposition are both unreliable and
    nothing was corrected for it).

    `skip_frames` defaults to the minimization frames that the NAMD input itself declares, read
    from its `minimize` and `dcdfreq` lines. Those frames are NOT dynamics: minimization carries
    a step number but no time, so counting them as simulated time is wrong. Pass 0 only if the
    user asks, and then say the time axis is fiction for them. `trajectory.minimization_note` can
    also say the answer is UNKNOWN, which is not the same as "none": it means no NAMD input of
    this stage was found, so nothing was dropped and nobody knows whether something should have
    been. Pass that on as unknown.

    `reference_frame` is the frame everything is measured AGAINST: the pocket is selected around
    the ligand there, both RMSDs take their zero there, and the series need not start at 0 if it
    lies outside the analysed range. It is an index into the whole chained trajectory, counted
    from 0, and it defaults to the first analysed frame. It is the USER's choice; `ligand_rmsd.
    reference` says which frame was used, and that is what to report. A frame the trajectory does
    not have comes back as `ok: false` with the frame count, and so does any negative value other
    than the default: negatives are NOT counted from the end here. Pass the refusal on and ask the
    user for a frame in range rather than guessing one.

    `trajectory.span_ps` and `span_ns` are the elapsed time between the FIRST and the LAST
    analysed frame, which is one frame spacing less than frames x dt. Quote the field, do not
    recompute it.

    `stage` names the CHARMM-GUI stage to analyse (`step5`, `step4`, ...) when the folder holds
    more than one and the user does not want the latest. It is the USER's choice: offer it when
    `trajectory.stages_present` shows another stage they might have meant, and never pass it on a
    guess.

    `reading` is an OFFERED INTERPRETATION and is the only part of the reply that is not a
    measurement. Present it as such, always with `reading.thresholds_used`, and pass on
    `reading.disclaimer`: those thresholds are conventions, not measurements of this system, and
    the researcher's own judgement decides. Never upgrade it to a conclusion of your own and
    never call a pose stable on your own authority.

    What is NOT measured, and is listed in `not_measured`: binding affinity or any free energy
    (none of these metrics is one, and a stable pose is not a potent ligand), bridging waters,
    and whether the pose is correct. Repeat that list when the user asks what the simulation
    proves.

    `include_series` defaults to false and returns only the summary of each series. Pass true
    when the user wants the per-frame values, for a plot or to see when something changed.
    """
    return md_analyze.analyze_md(
        folder, topology=topology or None, ligand_resname=ligand_resname,
        ligand_segid=ligand_segid, pocket_cutoff=pocket_cutoff, skip_frames=skip_frames,
        reference_frame=reference_frame, include_series=include_series,
        stage=stage or None)


@mcp.tool()
async def analyze_admet(
    uniprot_id: str,
    structure: str = "",
    pocket: str = "",
    conformer: str = "",
    pose: int = 1,
) -> dict:
    """Phase 3: retrieved ADME records and measurements plus locally calculated properties for the
    ligand of one pocket. NOTHING HERE IS PREDICTED and NO RISK IS ASSESSED: there is no hERG
    result at all, and no classification of anything. If a prediction, a risk or a class was
    asked for, say this tool does not produce one rather than estimating it. (A CYP or
    cytotoxicity row can appear, but only inside `rows_excluded_by_allowlist`, unclassified and
    quoted from its own assay description -- never as an endpoint of this tool.)

    The compound is the one that was docked and simulated, identified by its InChIKey.

    `structure` and `pocket` as in build_md_package; with neither, the user's choice recorded by
    select_structure is used. The compound is read from the docked pose when a run exists, and
    from the crystallographic SDF otherwise; `compound_identity.source` says which, and
    `identity_check` says whether the two files agree. Nothing is downloaded except the ChEMBL
    query.

    THE THREE BLOCKS ARE THREE DIFFERENT KINDS OF STATEMENT AND MUST NOT BE MIXED IN ONE
    SENTENCE, ONE TABLE OR ONE LIST:
    - `chembl_record`: facts of the ChEMBL entry (`max_phase`, `first_approval`,
      `withdrawn_flag`, `black_box_warning`, the route flags, ATC). These are RECORDS, not
      measurements and not predictions. `pref_name` is None for most research compounds: report
      it as unnamed, never substitute a synonym, a PDB chemical-component name or a guess.
      `name_status` says in words whether any retrieved source names this compound, and when
      none does you must not name it, give it a code, or say what it IS (an analogue, dimer or
      linked version of some named drug): there is no source for that. Pass `name_status` on.
      `fields_minus_one` and `boolean_fields_false` name values that are NOT findings: a field
      returned as -1 means nothing this project can state, and a flag returned as False is the
      value an unannotated entry carries, so never report "this compound is not oral" or "was
      not withdrawn" from one. Pass on `defaults_note`.
    - `measured_adme`: retrieved measurements, each with ITS OWN units, relation and assay
      description. Quote every value with its units exactly as given and NEVER convert, average
      or pool them: one standard_type arrives in incompatible units across rows. A row with a
      `relation` of '>' or '<' is a CENSORED value and the relation is part of the number -- "> 30
      uM" is not "30 uM". Cite `assay_description` when you quote a row: an in-vitro clearance in
      hepatocytes and an in-vivo clearance are both "CL" and are not comparable. A row with
      `description_says_value_was_calculated` is one whose own assay description says the
      depositor calculated or estimated the value: report it as such, not as a measurement.
      When `truncated` is true, every count is a LOWER BOUND and you must say so: ChEMBL held
      `rows_total` rows and only `rows_read` were read.
    - `calculated_properties`: `rdkit` (computed here from the SDF) and `chembl` (ChEMBL's own
      engine) side by side. Report them as calculated, name the engine, and do NOT present their
      agreement as confirmation or a measurement: if they disagree, neither can be called correct
      from here.

    `rule_reading` (Lipinski, Veber, Egan) carries `status: "READING OFFERED, NOT MEASURED"`.
    Those cut-offs are CONVENTIONS pre-registered by this project, not measurements and not
    citations: pass on `rule_reading.disclaimer`, give the violation COUNT and which cut-off was
    exceeded, and never turn it into a verdict of your own. Breaking a rule is not a prediction of
    failure -- vemurafenib breaks Lipinski on logP and `chembl_record` shows it was approved in
    2011 -- so say both without implying a contradiction.

    NOTHING IN THIS OUTPUT IS PREDICTED. No value comes from a trained model, so do not call any
    of it a prediction, and do not add endpoints that are absent (hERG, Caco-2 class, CYP
    inhibition, clearance class, toxicity). hERG and other off-target safety data are NOT here by
    design.

    Always report the selection with `selection`, which gives "N of total" for the measurements.
    When `measured_adme.retrieved` is false, distinguish the two cases it reports: with
    `no_chembl_entry` no query was made because the compound is not in ChEMBL; otherwise the
    lookup FAILED and the answer is UNKNOWN, never "no ADME data". The same applies to
    `chembl_record.status`, which is "UNKNOWN (the ChEMBL lookup failed)" when even the molecule
    lookup did not complete.

    `rows_excluded_by_allowlist_by_standard_type` counts the rows left out by THIS PROJECT'S
    allowlist, and `rows_excluded_by_allowlist` gives them one by one with their own assay
    description so they can be read. That list is CAPPED: report `rows_excluded_shown` ("N of
    total") with it and never present the rows shown as all there are. Do NOT call them "not ADME": the allowlist is narrower than
    ADME on purpose, and the generic types `Activity` and `Ratio` carry real ADME in some rows
    (efflux ratios, CYP induction) mixed with cytotoxicity in others. Quote one only as what its
    own description says, never as an allowlisted measurement, and never add them to the ADME
    count. `assay_types_not_queried` names what was never looked at, while
    `assay_types_skipped_at_the_row_cap` and `assay_types_skipped_at_the_time_budget` name types
    whose row counts are UNKNOWN for those two different reasons -- a slow source, not an empty
    one. Any count that comes with `truncated` true is a LOWER BOUND.
    """
    return admet.compact(await admet.analyze_admet(uniprot_id, structure, pocket, conformer, pose))


def main():
    """Entry point of the `protein-mcp` script (pyproject.toml)."""
    mcp.run()  # stdio transport


if __name__ == "__main__":
    main()