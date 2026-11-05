# Demo video script (about 2-3 min, no voice)

Record with `Win + Alt + R` (Xbox Game Bar) or with Snipping Tool → Video. Keep the Claude Desktop window maximized, in a new chat. Add titles and subtitles afterwards.

Before recording, quit Claude Desktop from the system tray, open it, and check that the 11 `protein-research` tools appear. Use a protein that is **not** in `dossiers\` (e.g. **EGFR**) so that everything is queried live.

| # | On-screen title | What you do | Subtitle |
|---|---|---|---|
| 0 | **protein-research MCP** | Empty chat screen (3 s) | Target assessment from public databases, inside Claude Desktop |
| 1 | **1 · Resolve the protein** | Type: `Review the protein EGFR in human` | The server resolves the name to a UniProt accession and asks if it is ambiguous |
| 2 | **2 · Target profile** | `Give me the target profile: function, diseases and tractability` | UniProt annotation + Open Targets evidence, each value with its source |
| 3 | **3 · Structures** | `Show experimental structures with the mutation T790M and the AlphaFold model` | RCSB PDB entries are confirmed for the variant; AlphaFold confidence is reported |
| 4 | **4 · Ligands** | `Which ligands and drugs act on EGFR?` | ChEMBL potency per mutation and assay type; approved drugs first |
| 5 | **5 · Literature** | `Find recent papers about EGFR inhibitor resistance` | Europe PMC search; "N of total" is always stated |
| 6 | **6 · Report** | `Build the report`, then open the PDF at the returned path | Deterministic PDF: data and links only, no generated text |
| 7 | **Every claim has a source** | Scroll through the PDF and click a link (PubMed or ChEMBL) | Every link comes from identifiers returned by the sources |
| 8 | **github.com/leorosas7777/protein-mcp** | Final screen (3 s) | Python · MCP · UniProt · Open Targets · RCSB · AlphaFold · ChEMBL · Europe PMC |

Tips:
- Let each reply finish before typing the next prompt, and wait 2 more seconds.
- If a reply is very long, scroll quickly; it can be sped up in editing.
- If something fails, do not stop: note the minute. It also serves as evidence.
