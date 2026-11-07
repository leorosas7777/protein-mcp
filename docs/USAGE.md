# How to use protein-research MCP

This is a local server that Claude Desktop uses as a set of tools. You write in the chat in plain language ("review the protein PCSK9") and Claude decides which tools to call. Everything the server retrieves is saved to disk (`dossiers\<accession>\`) and can then become a PDF report with links to every source.

## 1. Installation (once)

Requirements: Windows, [uv](https://docs.astral.sh/uv/), Claude Desktop and internet access. The server queries UniProt, Europe PMC, RCSB PDB, AlphaFold DB, ChEMBL, Open Targets and PDBe PISA, and needs no API keys.

1. Clone the repository to `C:\bioinformatics\protein-mcp` and run `uv sync` inside the folder.
2. In Claude Desktop, go to *Settings → Developer → Edit Config* and add this to `claude_desktop_config.json`:

   ```json
   {
     "mcpServers": {
       "protein-research": {
         "command": "uv",
         "args": ["--directory", "C:\\bioinformatics\\protein-mcp", "run", "protein-mcp"]
       }
     }
   }
   ```
3. Close Claude Desktop **from the system tray** (Quit, not just the window) and open it again. In a new chat, 11 `protein-research` tools should appear.

Repeat step 3 every time the code changes.

### PyMOL (optional)

`open_in_pymol` always writes the session `dossiers\<accession>\prep\pymol\session.pml`. It only opens PyMOL if it finds it: first through the `PROTEIN_MCP_PYMOL` variable, then as `pymol` on the PATH. Without PyMOL, it gives you the path of the `.pml` to open with *File → Run Script...*.

While that PyMOL window is open, you can ask for view changes in the chat (`pymol_view`): "show only 3OG7", "focus on the 032 pocket", "colour chain B red", "show residues 594-596" (UniProt numbering), "hide the boxes", "go back to the initial view", "what is loaded?". Only these actions exist, and the chat never sends PyMOL commands. Labels, surfaces and distance measurements are done from the PyMOL menus. Changes are not saved: reopening the session restores the initial view. The session opens a port on your computer only (127.0.0.1) with a random key, and it closes when PyMOL closes.

Open-source PyMOL in a separate environment (tested on Windows 11 with Python 3.13, 2026-10-26; about 207 MB):

```
uv venv C:\pymol-env --python 3.13
uv pip install --python C:\pymol-env\Scripts\python.exe pymol-open-source-whl PyQt5
```

The wheel has no graphical interface of its own: without `PyQt5` (or `PySide6`), PyMOL exits on startup and the tool tells you so. Then add this to the server entry in `claude_desktop_config.json`:

```json
"env": { "PROTEIN_MCP_PYMOL": "C:\\pymol-env\\Scripts\\python.exe" }
```

If you already have PyMOL (Schrödinger or conda), point the variable to its `pymol.exe`. Those two distributions have not been tested yet.

### AutoDock Vina (optional, for `dock_ligand`)

`dock_ligand` runs the official AutoDock Vina 1.2.7 executable. The `vina` package on PyPI has no Windows wheel. Download `vina_1.2.7_win.exe` from the v1.2.7 release of `forlilab/AutoDock-Vina` on GitHub and put it in a folder whose path has only ASCII characters, such as `C:\vina\`. Check its SHA-256: `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5`. Then add it to the same `env`:

```json
"env": { "PROTEIN_MCP_PYMOL": "C:\\pymol-env\\Scripts\\python.exe",
         "PROTEIN_MCP_VINA": "C:\\vina\\vina_1.2.7_win.exe" }
```

RDKit and Meeko are installed with the project. One call docks one ligand. On a 12-thread laptop, Vina itself took 61 to 124 seconds across the runs measured on 2026-11-05, plus a few seconds for Meeko. With little free memory it took much longer (one run did not finish in 15 minutes and was stopped). Close memory-hungry programs before docking. Because a chat request is cut long before that, `dock_ligand` only starts the run and `dock_status` reads its result.

## 2. Recommended flow

| Step | What you ask | Tool |
|---|---|---|
| 1. Identify | "Find the protein PCSK9" | `resolve_protein` |
| 2. Collect data | "Give me the target profile", "structures", "ligands", "literature" | `get_protein_profile`, `get_structures`, `get_ligands`, `search_literature` |
| 3. Report | "Build the report" | `build_report` |
| 4. Prepare structures (optional) | "Prepare the BRAF V600E structures" | `prepare_structures` |
| 5. View them (optional) | "Open them in PyMOL" | `open_in_pymol` |
| 6. Choose (optional) | "I choose 8C7Y, pocket TXV U, because ..." | `select_structure` (recorded in the manifest and in section 6 of the report) |
| 6b. Interfaces (optional, targets of biologics) | "Which chains touch the protein in 3WD5?", "analyze the interfaces" | `analyze_interfaces` (partners, construct sequences, RCSB assembly, PISA's view, contacts; section 6 of the report) |
| 7. Change the view (optional, with PyMOL open) | "Colour chain B red", "show residues 594-596" | `pymol_view` |
| 8. Re-docking control (optional, needs Vina) | "Re-dock the ligand of the pocket I chose" | `dock_ligand` starts it in the background and answers in seconds (crystal ligand only; the score is calculated, not an affinity) |
| 8b. Its result | "Is the docking done?" | `dock_status` (also stops a run); docking takes about a minute or two, and only one runs at a time |
| 8c. Experiment record | "Write up the experiments" | `build_experiment_report`: a separate PDF/Markdown record of the docking runs with their provenance; molecular dynamics and ADMET say "not run". Its calculated values are not mixed with the measured potency of the research report |

- If several proteins are possible, Claude should ask you which one you mean. Always confirm the **accession** (e.g. `Q8NBP7`): everything else is saved under that code.
- The report only includes what was collected in step 2. If something is missing, the PDF says so ("Missing sections") instead of making it up.
- Each query with a different mutation, topic, sort order or preprint setting is saved as a separate file (`ligands.json`, `ligands_g12c.json`...), and the report includes them all. Size parameters (`max_ligands`, `max_results`, `max_models`, `max_diseases`) do **not** change the file name: asking for 50 ligands overwrites the `ligands.json` of the default call.

## 3. Prompts that work well

**Identify the protein**
- "Review the protein KRAS" / "Find BACE1 in human"
- "Find Kras in mouse" (uses the organism; Open Targets only covers human)

**Target profile** (UniProt + Open Targets)
- "What does PCSK9 do and which diseases is it associated with?"
- "Is this target tractable with small molecules or antibodies?"
- "Show the 20 most associated diseases" (adjusts `max_diseases`)

**Structures** (RCSB PDB + AlphaFold)
- "Give me KRAS structures with the G12D mutation"
- "Which KRAS structures have mutations at position 12?" (`G12` accepts any substitution)
- "Show 10 structures" (maximum 25)

**Ligands** (ChEMBL)
- "Which ligands are there for KRAS G12C?"
- "List the 30 most potent compounds against PCSK9" (maximum 50)

**Literature** (Europe PMC)
- "Find recent papers on KRAS inhibitors" (topic + sort by date)
- "The most cited papers on the crystal structure of BACE1"
- "Include preprints" (not peer reviewed; they are flagged)

**Report**
- "Build the KRAS report" returns the paths of the PDF and the Markdown file.

## 4. What it does not do (yet) or does poorly

| Request | What happens | Alternative |
|---|---|---|
| "Find targets for Alzheimer's" (disease to targets) | There is no tool for that | Start from a specific protein |
| "What does sotorasib do?" (compound to data) | There is no search by compound (`lookup_compound` is pending) | Ask about the target protein |
| Drugs that are **antibodies** or proteins (e.g. evolocumab on PCSK9) | They have no IC50/Ki/Kd potency, so they are left out of the ranking; they appear in a separate table with name, type and phase | Read that table before concluding there are no drugs |
| Selectivity, toxicity, ADME, negative data | Not retrieved | Out of current scope |
| Clinical trials (ClinicalTrials.gov) | Not connected | Pending |
| Labels, surfaces or distances in PyMOL from the chat | Not among the actions of `pymol_view` | PyMOL menus |
| Docking a ChEMBL compound or your own molecule, molecule design, predictions | `dock_ligand` only re-docks the crystal ligand of a pocket | Planned for Phase 3 |
| Docking into a pocket that spans two chains (e.g. 5JSM 6NB), or a receptor with GDP/ATP kept | Refused before Vina starts (too slow; no template in Meeko) | Choose another pocket |
| Ambiguous names ("insulin receptor") | Name search returns noise | Use the gene symbol (INSR) |
| Free-text mutations ("the most common mutation") | The filter needs notation such as `G12C` or `G12` | Specify the variant |
| "Write the conclusion in the PDF" | The PDF only contains sourced data; interpretation stays in the chat | Ask Claude to interpret it, citing the report |

## 5. Where files are saved

```
dossiers\Q8NBP7\
  profile.json              full result of get_protein_profile
  structures.json           get_structures (default call)
  ligands.json              get_ligands (default call)
  literature.json           search_literature (default call)
  report_Q8NBP7.pdf         report
  report_Q8NBP7.md          the same report in Markdown
  prep\                     Phase 2: manifest.json, raw\, clean\, ligands\, pymol\
```

- `dossiers/` is not committed to git. To start over with a protein, delete its folder.
- Repeating a query with the same parameters overwrites its file with new data.

## 6. Common problems

- **The tools do not appear**: quit Claude Desktop from the system tray and open a new chat. Check that the path in the JSON is correct.
- **"Could not write the report... close it"**: the PDF is open in another program. Close it and try again.
- **Very long or cut-off reply**: ask for fewer results (`max_results`, `max_ligands`). The full detail is in the dossier.
- **"inactive" or "not found"**: the accession was deleted or does not exist in UniProt. Go back to `resolve_protein`.
- **"is inactive (MERGED...)" or "secondary accession"**: UniProt merged that accession into another. The error names the new one; confirm that it is the protein you mean and use it.
- **`accession_warning` in a result**: UniProt did not answer, so the accession is unconfirmed. An empty result does not mean there are no data; repeat the query.
