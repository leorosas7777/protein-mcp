# Cómo usar protein-research MCP

Servidor local que Claude Desktop usa como "herramientas". Tú escribes en el chat en lenguaje normal ("revisa la proteína PCSK9") y Claude decide qué herramientas llamar. Todo lo que el servidor recupera se guarda en disco (`dossiers\<accession>\`) y al final puede convertirse en un reporte PDF con enlaces a cada fuente.

## 1. Instalación (una vez)

Requisitos: Windows, [uv](https://docs.astral.sh/uv/), Claude Desktop y conexión a internet (el servidor consulta UniProt, Europe PMC, RCSB PDB, AlphaFold DB, ChEMBL y Open Targets; no necesita claves).

1. Clona el repositorio en `C:\bioinformatics\protein-mcp` y ejecuta `uv sync` dentro de la carpeta.
2. En Claude Desktop: *Settings → Developer → Edit Config* y añade en `claude_desktop_config.json`:

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
3. Cierra Claude Desktop **desde la bandeja del sistema** (Quit, no solo la ventana) y ábrelo otra vez. En un chat nuevo deben aparecer 11 herramientas de `protein-research`.

Cada vez que cambie el código hay que repetir el paso 3.

### PyMOL (opcional)

`open_in_pymol` siempre escribe la sesión `dossiers\<accession>\prep\pymol\session.pml`. Solo abre PyMOL si lo encuentra: primero en la variable `PROTEIN_MCP_PYMOL`, luego como `pymol` en el PATH. Sin PyMOL te da la ruta del `.pml` para abrirlo con *File → Run Script...*.

Mientras esa ventana de PyMOL esté abierta, puedes pedir cambios de vista en el chat (`pymol_view`): "muestra solo 3OG7", "enfoca el bolsillo 032", "colorea la cadena B de rojo", "muestra los residuos 594-596" (numeración UniProt), "oculta las cajas", "vuelve a la vista inicial", "¿qué está cargado?". Solo existen esas acciones: el chat nunca envía comandos de PyMOL. Etiquetas, superficies o medir distancias se hacen en los menús de PyMOL. Los cambios no se guardan: al reabrir la sesión vuelve la vista inicial. La sesión abre un puerto solo en tu equipo (127.0.0.1) con una clave aleatoria; se cierra al cerrar PyMOL.

PyMOL open-source en un entorno aparte (probado en Windows 11 con Python 3.13, 2026-10-26; ~207 MB):

```
uv venv C:\pymol-env --python 3.13
uv pip install --python C:\pymol-env\Scripts\python.exe pymol-open-source-whl PyQt5
```

La rueda no trae interfaz gráfica: sin `PyQt5` (o `PySide6`) PyMOL se cierra al arrancar y la herramienta lo avisa. Después añade al servidor en `claude_desktop_config.json`:

```json
"env": { "PROTEIN_MCP_PYMOL": "C:\\pymol-env\\Scripts\\python.exe" }
```

Si ya tienes PyMOL (Schrödinger o conda), apunta la variable a su `pymol.exe`. Ni el arranque desde Claude Desktop con esta variable ni un PyMOL de Schrödinger están probados todavía.

## 2. Flujo recomendado

| Paso | Qué pides | Herramienta |
|---|---|---|
| 1. Identificar | "Busca la proteína PCSK9" | `resolve_protein` |
| 2. Recoger datos | "Dame el perfil de la diana", "estructuras", "ligandos", "literatura" | `get_protein_profile`, `get_structures`, `get_ligands`, `search_literature` |
| 3. Reporte | "Genera el reporte" | `build_report` |
| 4. Preparar estructuras (opcional) | "Prepara las estructuras de BRAF V600E" | `prepare_structures` |
| 5. Verlas (opcional) | "Ábrelas en PyMOL" | `open_in_pymol` |
| 6. Elegir (opcional) | "Elijo 8C7Y, bolsillo TXV U, porque ..." | `select_structure` (queda en el manifiesto y en la sección 6 del reporte) |
| 6b. Interfaces (opcional, dianas de biológicos) | "¿Qué cadenas tocan la proteína en 3WD5?", "analiza las interfaces" | `analyze_interfaces` (parejas, secuencias del constructo, ensamblaje de RCSB, opinión de PISA, contactos; sección 6 del reporte) |
| 7. Cambiar la vista (opcional, con PyMOL abierto) | "Colorea la cadena B de rojo", "muestra los residuos 594-596" | `pymol_view` |

- Si hay varias proteínas posibles, Claude debe preguntarte cuál es. Confirma siempre el **accession** (p. ej. `Q8NBP7`): todo lo demás se guarda bajo ese código.
- El reporte solo incluye lo que ya se recogió en el paso 2. Si falta algo, el PDF lo dice ("Missing sections") en vez de inventarlo.
- Cada consulta con distinta mutación, tema, orden o preprints se guarda como archivo aparte (`ligands.json`, `ligands_g12c.json`...). El reporte las incluye todas. Los parámetros de tamaño (`max_ligands`, `max_results`, `max_models`, `max_diseases`) **no** cambian el nombre: pedir 50 ligandos sobrescribe `ligands.json` de la llamada por defecto.

## 3. Prompts que funcionan bien

**Identificar la proteína**
- "Revisa la proteína KRAS" / "Busca BACE1 en humano"
- "Busca Kras en ratón" (usa el organismo; Open Targets solo cubre humano)

**Perfil de la diana** (UniProt + Open Targets)
- "¿Qué función tiene PCSK9 y con qué enfermedades se asocia?"
- "¿Es una diana tratable con molécula pequeña o con anticuerpos?"
- "Muestra las 20 enfermedades más asociadas" (ajusta `max_diseases`)

**Estructuras** (RCSB PDB + AlphaFold)
- "Dame estructuras de KRAS con la mutación G12D"
- "¿Qué estructuras hay de KRAS con mutaciones en la posición 12?" (`G12` acepta cualquier sustitución)
- "Muestra 10 estructuras" (máximo 25)

**Ligandos** (ChEMBL)
- "¿Qué ligandos hay para KRAS G12C?"
- "Lista los 30 compuestos más potentes contra PCSK9" (máximo 50)

**Literatura** (Europe PMC)
- "Busca artículos recientes sobre inhibidores de KRAS" (tema + orden por fecha)
- "Los artículos más citados sobre la estructura cristalina de BACE1"
- "Incluye preprints" (no revisados por pares; se marcan)

**Reporte**
- "Genera el reporte de KRAS" → devuelve la ruta del PDF y del Markdown.

## 4. Lo que NO hace (todavía) o hace mal

| Petición | Qué pasa | Alternativa |
|---|---|---|
| "Busca dianas para el Alzheimer" (de enfermedad a dianas) | No hay herramienta para eso | Parte de una proteína concreta |
| "¿Qué hace el sotorasib?" (de compuesto a datos) | No hay búsqueda por compuesto (`lookup_compound` está pendiente) | Pregunta por la proteína diana |
| Fármacos que son **anticuerpos** o proteínas (p. ej. evolocumab en PCSK9) | No tienen potencia IC50/Ki/Kd: no entran al ranking; salen en una tabla aparte con nombre, tipo y fase | Lee esa tabla antes de concluir que no hay fármacos |
| Selectividad, toxicidad, ADME, datos negativos | No se recuperan | Fuera del alcance actual |
| Ensayos clínicos (ClinicalTrials.gov) | No conectado | Pendiente |
| Cambiar la vista de PyMOL desde el chat ("colorea la cadena B") | Todavía no: la sesión se abre ya preparada | Menús y panel de objetos de PyMOL |
| Docking, diseño de moléculas, predicciones | No es una herramienta de cálculo | — |
| Nombres ambiguos ("insulin receptor") | La búsqueda por nombre devuelve ruido | Usa el símbolo del gen (INSR) |
| Mutaciones en texto libre ("la mutación más común") | El filtro necesita notación tipo `G12C` o `G12` | Especifica la variante |
| "Escribe la conclusión en el PDF" | El PDF solo contiene datos con fuente; la interpretación queda en el chat | Pide a Claude que interprete citando el reporte |

## 5. Dónde quedan los archivos

```
dossiers\Q8NBP7\
  profile.json              resultado completo de get_protein_profile
  structures.json           get_structures (por defecto)
  ligands.json              get_ligands (por defecto)
  literature.json           search_literature (por defecto)
  report_Q8NBP7.pdf         reporte
  report_Q8NBP7.md          el mismo reporte en Markdown
```

- `dossiers/` no se sube a git. Para empezar de cero con una proteína, borra su carpeta.
- Repetir una consulta con los mismos parámetros sobrescribe su archivo con datos nuevos.

## 6. Problemas frecuentes

- **No aparecen las herramientas**: cierra Claude Desktop desde la bandeja y abre un chat nuevo. Revisa que la ruta del JSON sea correcta.
- **"Could not write the report... close it"**: el PDF está abierto en otro programa; ciérralo y repite.
- **Respuesta enorme o cortada**: pide menos resultados (`max_results`, `max_ligands`); el detalle completo está en el expediente.
- **"inactive" o "not found"**: el accession está eliminado o no existe en UniProt; vuelve a `resolve_protein`.
