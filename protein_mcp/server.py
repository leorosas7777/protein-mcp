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

from protein_mcp import report, target_profile  # noqa: E402
from protein_mcp.common import parse_accession, safe  # noqa: E402
from protein_mcp.dossier import attach_dossier, dossier_name  # noqa: E402
from protein_mcp.files import choice, control, prepare, viewer  # noqa: E402
from protein_mcp.sources import chembl, europepmc, uniprot  # noqa: E402
from protein_mcp.sources.alphafold import fetch_alphafold  # noqa: E402
from protein_mcp.sources.rcsb import MAX_PDB_MODELS, fetch_pdb  # noqa: E402

mcp = FastMCP("protein-research")


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
    The same result is also saved to `dossier_path` for build_report.
    """
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
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
    Cite potency with its `source` link. If target is None, say "sin registro en ChEMBL".
    """
    result = await chembl.get_ligands(uniprot_id, mutation, max_ligands)
    accession, _ = parse_accession(uniprot_id)
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
    - If open_targets.found is false, say "sin registro en Open Targets". If
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



def main():
    """Entry point of the `protein-mcp` script (pyproject.toml)."""
    mcp.run()  # stdio transport


if __name__ == "__main__":
    main()