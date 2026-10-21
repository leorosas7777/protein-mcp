# Guion del video de demostración (≈ 2-3 min, sin voz)

Grábalo con `Win + Alt + R` (Xbox Game Bar) o con Recortes → Video. Ventana de Claude Desktop maximizada, chat nuevo. Los prompts van en inglés porque el video es para el portafolio. Los títulos y subtítulos se añaden después.

Antes de grabar: cierra Claude Desktop desde la bandeja, ábrelo y comprueba que aparecen las 6 herramientas de `protein-research`. Usa una proteína que **no** esté en `dossiers\` (p. ej. **EGFR**) para que todo se consulte en vivo.

| # | Título en pantalla | Qué haces | Subtítulo |
|---|---|---|---|
| 0 | **protein-research MCP** | Pantalla del chat vacío (3 s) | Target assessment from public databases, inside Claude Desktop |
| 1 | **1 · Resolve the protein** | Escribe: `Review the protein EGFR in human` | The server resolves the name to a UniProt accession and asks if it is ambiguous |
| 2 | **2 · Target profile** | `Give me the target profile: function, diseases and tractability` | UniProt annotation + Open Targets evidence, each value with its source |
| 3 | **3 · Structures** | `Show experimental structures with the mutation T790M and the AlphaFold model` | RCSB PDB entries are confirmed for the variant; AlphaFold confidence is reported |
| 4 | **4 · Ligands** | `Which ligands and drugs act on EGFR?` | ChEMBL potency per mutation and assay type; approved drugs first |
| 5 | **5 · Literature** | `Find recent papers about EGFR inhibitor resistance` | Europe PMC search; "N of total" is always stated |
| 6 | **6 · Report** | `Build the report` y abre el PDF de la ruta que devuelve | Deterministic PDF: data and links only, no generated text |
| 7 | **Every claim has a source** | Desplázate por el PDF; haz clic en un enlace (PubMed o ChEMBL) | 100 % of the links come from identifiers returned by the sources |
| 8 | **github.com/leorosas7777/protein-mcp** | Pantalla final (3 s) | Python · MCP · UniProt · Open Targets · RCSB · AlphaFold · ChEMBL · Europe PMC |

Consejos:
- Deja que cada respuesta termine antes de escribir la siguiente y espera 2 s más.
- Si una respuesta es muy larga, desplázate rápido; se acelera en la edición.
- Si algo falla, no pares: anota el minuto. También sirve como prueba.
