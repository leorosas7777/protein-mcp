# protein-research MCP

Local MCP server (Python, `uv`, `mcp[cli]<2`/FastMCP) for Claude Desktop. The user says "review protein X" and the server resolves the protein, retrieves literature, structures and ligands, and produces a report with hyperlinks. Portfolio project for a European company. Intended reader of the report: a drug-discovery analyst.

**At the start of a session, read all of `ESTADO_PROYECTO.md`.** It is the real state of the project, the decisions taken and the known limitations. Do not reopen decisions in section 5 without a reason.

`HISTORIAL_PROYECTO.md` (in `.gitignore`) is the same document before it was split on 2026-11-12, with the whole history unsummarised and the same headings: tests, reviewer rounds, closed steps and validations 8.1 to 8.5. **It is not read at the start of a session**: search it when the detail of something already closed is needed, instead of assuming it.

## Rules that are not negotiable
- **Never `print()` to stdout** (it breaks the MCP stdio channel). Use `logging`.
- **Everything in this repository is written in English** (the user decided on 2026-11-04 that everything the project produces goes in English; extended on 2026-11-23 to the whole repository, because it is public and the intended reader does not read Spanish): tool output, reports, docstrings, `README.md`, `docs/`, this file, `ESTADO_PROYECTO.md` and `.claude/agents/`, in any phase. Spanish is only for the conversation.
  - **Three exceptions, which are DATA and not prose: do not translate them.** (1) The `LEGACY_NOTES` keys and `_LEGACY_GENES_RE` in `protein_mcp/report.py`: they must stay the exact Spanish text that was written into dossiers before 2026-11-04, because they are the lookup keys that translate those old records on read. Translating them makes every pre-2026-11-04 dossier stop being translated, silently, with every test still passing. (2) The fixtures that feed them in `tests/test_report.py`. (3) The non-ASCII path fixtures in `tests/test_dock.py` (`"carpeta ñá"`, `C:\Users\Jesús\...`), which exist to test Unicode paths. Still pending and not an exception: the Spanish `select_structure` reason quoted inside the public English report (`ESTADO_PROYECTO.md` section 8).
- **Anti-hallucination:** every link is built only from IDs the source returned; every claim in the report must have a retrieved source; fields with no data stay `None` or "no record", and are never filled in or inferred.
- **Always report the selection** ("N of total").
- **Validate against the real API before writing code** for a new field. A WebFetch summary is not evidence. Simulated tests must copy the real shapes and real defects of the APIs.
- **Before recording anything as implemented in `ESTADO_PROYECTO.md`, verify it in the code and in a real test** (this document has already described functions as done that did not exist).
- Question the user's decisions if they get in the way of the project succeeding; do not confirm biases.

## Architecture
- Package `protein_mcp/` (structure in `ESTADO_PROYECTO.md` section 3). `protein_mcp/server.py` only registers tools (`@mcp.tool()`) and delegates. The docstrings stay there (Claude Desktop reads them in order to use each tool).
- One module per source in `protein_mcp/sources/`: `uniprot.py`, `europepmc.py`, `rcsb.py`, `alphafold.py`, `chembl.py`, `opentargets.py`. `common.py` holds helpers and imports nothing from the project (this avoids import cycles); `dossier.py` writes and reads dossiers. Absolute imports (`from protein_mcp.sources import uniprot`). A single server for every phase; Phase 2 lives in `protein_mcp/files/`.
- Modules do not import each other except where `ESTADO_PROYECTO.md` section 3 says so.

## Commands
- Run the server: `uv --directory C:\bioinformatics\protein-mcp run protein-mcp` (after changing `pyproject.toml`: `uv sync`)
- Inspector: `uv run mcp dev protein_mcp/server.py`
- Offline tests: `uv run python tests/test_chembl.py` (use the same pattern for the rest).
- Real tests: `uv run python tests/smoke_live.py [--skip-slow | kras pcsk9 egfr notfound tnf]` (real APIs, dossiers in a temporary folder; external, regression and anti-hallucination invariants).
- Review: after changing a tool, run the `revisor` agent (`.claude/agents/revisor.md`) in a clean context before the commit.
- To make Claude Desktop reload the code: full Quit from the tray, then open a new chat. `query_used` (resolve_protein) and `selection_reason` (AlphaFold) tell which version is running.
- Environment: Windows, PowerShell. Working language of the conversation: Spanish.

## Workflow
- One task per session. When it closes, update `ESTADO_PROYECTO.md` (only with what has been verified).
- The user makes the commits. Do not `git push`, do not delete files, do not version backups (`server_old.py`), `dossiers/`, or test output.
- If a change goes wrong and going back to an earlier commit to rebuild from there is cleaner than carrying on patching, say so explicitly: which commit, why, and which uncommitted changes would be lost. Do not run `git reset`/`git checkout`/`git restore` without the user's approval.
- Keep replies short; show changes as diffs except for large refactors.
- If you change the behaviour of a tool, adjust its docstring in `protein_mcp/server.py` and its test.
