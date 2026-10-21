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

Logic lives in one module per source: uniprot, europepmc, rcsb, alphafold,
chembl, opentargets; target_profile combines uniprot + opentargets.
Shared helpers are in common.

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

import europepmc  # noqa: E402
import report  # noqa: E402
import target_profile  # noqa: E402
import uniprot  # noqa: E402
import chembl
from alphafold import fetch_alphafold  # noqa: E402
from common import attach_dossier, dossier_name, parse_accession, safe  # noqa: E402
from rcsb import MAX_PDB_MODELS, fetch_pdb  # noqa: E402

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
            0 (default) uses the MAX_PDB_MODELS constant in rcsb.py.
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
    in the report, never filled from the conversation.

    Rules for using the result:
    - Give the user `pdf_path` (and `markdown_path`). If `missing` is not empty,
      say which sections are absent and offer to run those tools and rebuild.
    - The report holds data and links only, no interpretation. If you add an
      interpretation in the chat, base it only on tool results and cite them.
    - If it fails because the PDF is open in a viewer, ask the user to close it.
    """
    return await asyncio.to_thread(report.build_report, uniprot_id)



if __name__ == "__main__":
    mcp.run()  # stdio transport