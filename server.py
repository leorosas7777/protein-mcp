"""Protein research MCP server. This file only registers the tools.

Tools:
  - resolve_protein: protein/gene name -> UniProt accession candidates.
  - search_literature: UniProt accession -> papers from Europe PMC
    (PubMed/MEDLINE records, PMC full-text records and, optionally, preprints).
  - get_structures: UniProt accession -> a few experimental PDB structures
    (RCSB) plus the AlphaFold DB model with its confidence.

Logic lives in one module per source: uniprot, europepmc, rcsb, alphafold.
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
import uniprot  # noqa: E402
from alphafold import fetch_alphafold  # noqa: E402
from common import parse_accession, safe  # noqa: E402
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
    """
    return await europepmc.search(uniprot_id, topic, max_results, sort, include_preprints)


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

    return {
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


if __name__ == "__main__":
    mcp.run()  # stdio transport