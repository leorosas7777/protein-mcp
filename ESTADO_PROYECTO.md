# ESTADO DEL PROYECTO: protein-research MCP
Última actualización: 2026-10-20. Leer completo al iniciar cualquier sesión o chat nuevo.

## 1. Objetivo
Servidor MCP local para Claude Desktop. El investigador dice "revisa la proteína X" y el servidor: resuelve la proteína, recupera literatura, estructuras y ligandos, y produce un reporte con hipervínculos (chat o PDF). PyMOL se abre **solo a petición del usuario**.
Meta final: proyecto de portafolio para ser admitido en una empresa europea.
**Lector objetivo del reporte: analista de descubrimiento de fármacos** (decidido por el usuario el 2026-10-17).

## 2. Entorno
- Windows, PowerShell, PyMOL gratuito. Idioma: español. Plan de Claude: gratis hasta 2026-10-20; decisión de pasar a Pro (incluye Claude Code) tomada ese día.
- Carpeta del proyecto: `C:\bioinformatics\protein-mcp\` (todos los .py en la raíz). Python con `uv`, `mcp[cli]<2` (1.30.0, FastMCP). Migrar a mcp v2 es una tarea aparte y diferida.
- Claude Desktop: `claude_desktop_config.json` → `uv --directory C:\bioinformatics\protein-mcp run server.py`.
- Para recargar código: Quit completo desde la bandeja + chat nuevo. El campo `query_used` de `resolve_protein` y `selection_reason` de AlphaFold indican qué versión corre.
- Node.js solo para el MCP Inspector (`uv run mcp dev server.py`).
- Nunca `print()` a stdout (rompe stdio); usar `logging`.
- Git/GitHub: repo privado, un commit por paso, sin secretos. El usuario hace los commits. (Verificar estado de commits: división, correcciones, `chembl.py`.)
- Chat de claude.ai: su shell NO alcanza ChEMBL, UniProt ni Open Targets (proxy rechaza); WebFetch devuelve un resumen hecho por otro modelo, no sirve como evidencia de formas JSON. Desde Claude Code en la máquina local se espera tener red completa (verificar con un `curl` en la primera sesión).

## 3. Arquitectura (decidida 2026-10-17)
- `server.py` es el único punto de entrada: solo registra herramientas (`@mcp.tool()`) y delega. Los docstrings se quedan ahí (es lo que Claude lee para usar cada herramienta).
- Un módulo por fuente: `uniprot.py`, `europepmc.py`, `rcsb.py` (no `pdb.py`: taparía el módulo `pdb` de Python), `alphafold.py`, `chembl.py`. `common.py` guarda helpers compartidos (`log`, `parse_accession`, `clean_text`, `short`, `safe`) y no importa nada del proyecto (evita ciclos). `chembl.py` tampoco importa módulos del proyecto (solo httpx).
- Grafo: server → uniprot, europepmc, rcsb, alphafold, chembl, common; europepmc → uniprot, common; rcsb/alphafold → common.
- **Un solo servidor para todas las fases** (config sin cambios). Al empezar la Fase 2: migrar a paquetes por fase con `register(mcp)` y mover `common.py` a un lugar compartido. Mover la Fase 3 a su propio servidor solo si necesita dependencias distintas, cómputo largo o subprocesos que bloqueen. No crear carpetas de fase antes de tiempo.

## 4. Herramientas implementadas
Probadas contra APIs reales, KRAS/P01116.
1. `resolve_protein(name, organism_id=9606, max_results=5)`: UniProt. Sin espacios: `gene_exact` + reviewed primero; si no, texto libre. Devuelve `candidates`, `query_used`, `needs_user_confirmation`.
2. `search_literature(uniprot_id, topic="", max_results=10, sort="relevance", include_preprints=False)`: Europe PMC (`resultType=core`, `TITLE_ABS`, `NOT SRC:PPR`). Limpia HTML, campos anidados de revista y tipo de publicación; enlaces pubmed/doi/pmc/europepmc.
3. `get_structures(uniprot_id, max_models=0, mutation="")`: RCSB (Search v2 + GraphQL) y AlphaFold DB en paralelo con `asyncio.gather`.
   - Constantes editables: `MAX_PDB_MODELS = 7` y `CANDIDATE_POOL = 100` (rcsb.py); `AF_ISOFORM_MARGIN = 5.0` (alphafold.py).
   - Selección diversa de PDB (mutaciones, nucleótido, con/sin ligando); informa "N de total".
   - Filtro de mutación: RCSB busca por texto completo y cada resultado se CONFIRMA (mutación nombrada en el campo de mutación o en el título; una posición sola como "G12" acepta cualquier sustitución). Los no confirmados se listan en `discarded_not_confirmed` (hasta 10) con `discarded_count`; `confirmed_in_examined` cuenta los confirmados. En el título se acepta el gen pegado a la variante ("KRASG12D", "KRasG12C", "NRASQ61K": 3+ letras delante) y solo letras de aminoácido reales, para que "KRAS4B" no dé "S4B". Prueba offline: `tests/test_rcsb.py`.
   - AlphaFold: se prefiere el modelo canónico salvo que una isoforma lo supere por ≥ `AF_ISOFORM_MARGIN` pLDDT; `selection_reason` explica la elección.
   - Ligandos filtrados: nucleótidos, aditivos y cofactores se separan.
   - Fallos aislados por fuente con `safe`.
   - Pruebas reales 2026-10-17 y 2026-10-20 (`get_structures("P01116")`, `mutation="G12D", max_models=4`, `mutation="Q61H"`): OK.
4. `get_ligands(uniprot_id, mutation="", max_ligands=15)`: ChEMBL, módulo `chembl.py`. Prueba real P01116 + G12C (v2): 5.070 actividades descargadas, 465 usadas, 348 ligandos con datos válidos, 17 devueltos; sotorasib y adagrasib (fase 4) y glecirasib (fase 2) salen primero. Prueba offline: `tests/test_chembl.py` (usa formas y defectos reales de KRAS).
   - Diana: solo el target `SINGLE PROTEIN` con un único componente igual al accession (KRAS: CHEMBL2189121). Los demás targets (PPI, familia, complejos) se listan en `other_targets_ignored`.
   - Actividades: `standard_type` IC50/Ki/Kd con `pchembl_value`; unidades nM; solo relación "=" entra a las estadísticas; excluye `potential_duplicate` y `data_validity_comment`.
   - Confianza: `confidence_score` se lee del ensayo (endpoint `assay`) y se aplica en el cliente (`MIN_CONFIDENCE = 8`); en KRAS todo vale 9, así que no filtra.
   - Mutación: `assay_variant_mutation` del registro; se agrupa por (molécula, mutación, tipo, `assay_class`) y nunca se mezclan mutaciones ni ensayos celulares con bioquímicos. "not stated" NO prueba tipo silvestre.
   - **Conflicto etiqueta/texto:** si la mutación nombrada en la descripción del ensayo no coincide con la etiqueta curada, la actividad se excluye y se lista en `label_conflicts_flagged` (hasta 10). Caso real: 23 actividades etiquetadas G12C cuyo texto dice "KRAS G12D" (ensayo de intercambio con bodipy).
   - Orden de salida: primero por `max_phase`, luego por mejor pChEMBL. Se añaden siempre los compuestos con mecanismo registrado (`mechanism` endpoint).
   - Identidad, propiedades (None si ChEMBL no las tiene), estado de desarrollo y documento fuente con enlaces PubMed/DOI. El nombre es el `pref_name` de ChEMBL y es None en la mayoría de los compuestos de investigación.
   - Constantes: `MAX_ACTIVITIES = 6000`, `MAX_LIGANDS = 15`, `PAGE = 1000`.
   - Alcance v1: bloques 1, 2, 3, 6 y 8 de la lista original (identidad, potencia, ensayo, desarrollo, fuente). Fuera: selectividad (4), datos negativos (7), cruce PDB↔ChEMBL.

## 5. Decisiones tomadas (no reabrir sin motivo)
- **Fuente de ligandos: solo ChEMBL** (licencia CC BY-SA 3.0, API gratuita; ChEMBL_37, ~24,5 M de actividades; verificado en `ebi.ac.uk/chembl/api/data/status.json`).
- **DrugBank se descarta por ahora.** Plan académico no comercial y sin API; uso comercial requiere licencia de pago. La cifra "500.000 fármacos" NO se sostiene (19.872 entradas en su página). Posible módulo opcional futuro en una rama, solo si la empresa destino tiene licencia. No se leyeron los términos legales completos.
- ClinicalTrials.gov: después y limitado (solo ligandos con max_phase ≥ 1, buscando por el nombre exacto de ChEMBL). Open Targets: irá en el perfil de proteína (validación de diana), no en ligandos.
- **Ligandos desde la literatura:** NO extraer nombres de compuestos de los abstracts (frágil y contra la regla anti-alucinación). `get_ligands` parte de la diana (UniProt → target ChEMBL → actividades). Cruce con la revisión de literatura por PMID/DOI. Más adelante, `lookup_compound(nombre)` con coincidencia exacta de nombre/sinónimo de ChEMBL; si no hay, "sin registro".
- Reglas de diseño: todo enlace se construye solo con IDs devueltos; toda afirmación del reporte debe tener una fuente recuperada; se informa cómo se seleccionó ("N de total"); el reporte y la recuperación de archivos son fases separadas; PyMOL a petición.
- **Claude Code (decidido 2026-10-20):** el trabajo de código pasa a Claude Code (incluido en Pro) porque ve el repo y puede llamar a las APIs desde la máquina local. La API de Anthropic NO se necesita: el servidor MCP no llama a ningún modelo. Los chats de claude.ai quedan para diseño y revisión. Límites compartidos entre Claude y Claude Code.
- Instrucción del proyecto: no confirmar sesgos del usuario y cuestionar decisiones que obstaculicen el éxito.
- **Lección:** este documento llegó a describir como hechas funciones que no existían en el código. Antes de anotar algo como implementado, verificarlo en el código (por ejemplo `Select-String -Pattern "NOMBRE"`) y en una prueba real. Un resumen de WebFetch no cuenta como evidencia de la forma de una API.

## 6. Limitaciones conocidas y bugs abiertos
- **RCSB, falso negativo del lado del servidor:** la búsqueda `full_text` de RCSB no separa el gen de la variante, así que entradas con el campo de mutación vacío y la variante solo pegada en el título no llegan al conjunto (verificado 2026-10-20: 9GLU "KRasG12D-GDP…" y 5UFQ no aparecen al buscar "G12D"). Nuestra confirmación sí las reconocería. Posible arreglo: añadir al OR una búsqueda por `pdbx_mutation` o por título con prefijos de gen; no validado.
- **Salida demasiado grande:** `get_ligands` con P01116/G12C devuelve ~55 KB y la herramienta la guarda en un archivo en lugar de mostrarla. Encadenar herramientas por el contexto del modelo no escala: hace falta el expediente en disco y salidas compactas.
- `get_ligands`: el ranking por mejor pChEMBL favorece mediciones únicas; solo se comprueba `max_phase` en los compuestos devueltos y en los de mecanismo registrado, así que otros compuestos clínicos sin mecanismo pueden no aparecer (el filtro `max_phase__gte` sobre molecule NO está validado). Los ligandos del PDB no se cruzan con ChEMBL. Un ensayo de ChEMBL puede estar mal etiquetado en el sentido contrario (texto sin mutación nombrada y etiqueta equivocada): no se detecta.
- Filtro de mutación de RCSB: solo se examinan los 100 mejores por resolución de los que RCSB encuentra por texto; la confirmación puede descartar entradas que nombran la variante solo en otros campos (se listan en `discarded_not_confirmed`).
- Con filtro de mutación, varias estructuras salen sin ligando de molécula pequeña; los análogos de nucleótido (GDP/GCP/GNP) ocupan slots (pendiente: agruparlos GDP vs análogo de GTP).
- 7ROV y 6WGN están unidas a péptidos (aparecen en `partner_chains`, `ligands: []`); pendiente contar binders poliméricos como ligandos.
- El DOI de 7T47 es un preprint (Research Square).
- Los caminos "no encontrado" no se han probado contra APIs reales (usar un accession bien formado que no exista, verificado en uniprot.org; P99999 NO sirve: es citocromo c humano).
- `resolve_protein("insulin receptor")` devuelve ruido; mejora diferida (búsqueda por nombre entre comillas).
- La relevancia de `search_literature` con `topic="inhibitor"` incluye artículos poco ligados a KRAS (p. ej. PP2A/EGFR); es el comportamiento de `TITLE_ABS` y se aceptó.

## 7. Historial reciente
- 2026-10-15/17: división en módulos (sin cambiar lógica). Pruebas reales OK, con dos fallos preexistentes detectados: AlphaFold elegía la isoforma P01116-2 (pLDDT 92,56 vs 91,5 del canónico) y el filtro `mutation="G12D"` dejaba pasar 9E3S (G12N). Corregidos en un commit aparte.
- 2026-10-17: `get_ligands` v1 implementado y probado contra ChEMBL real. La prueba con G12C reveló tres defectos que se corrigieron en v2: un ligando #1 con etiqueta G12C pero ensayo G12D, mezcla de ensayos celulares y bioquímicos, y fármacos enterrados bajo preclínicos de una sola medición. `get_structures` repasado: nuevo falso negativo (9TBR), sin corregir.
- 2026-10-20: se decide pasar a Claude Code. Se prepara `CLAUDE.md`. Se planeó `probe_profile.ps1` para validar campos reales de UniProt y Open Targets, pero el archivo NO existe en el repo (verificado 2026-10-20).
- 2026-10-20 (Claude Code): `.gitignore` creado; `test_chembl.py` movido a `tests/` (pasa: ALL OK); commit y push del paso de ChEMBL.
- 2026-10-20 (Claude Code): validados los campos de UniProt y Open Targets para el perfil (sección 8.1). `pyproject.toml`: descripción real, se quitan `[project.scripts]` y `[build-system]` (apuntaban a la plantilla `src/protein_mcp/`, sin uso; el proyecto pasa a `virtual` en `uv.lock`). `uv sync` quitó del `.venv` paquetes que no estaban en el lock (`httpx2`, `mcp-types`, `opentelemetry-api`, `truststore`); el servidor carga sus 4 herramientas y las pruebas pasan.
- 2026-10-20 (Claude Code): **corregido 9TBR** en `rcsb.py` (`_LEFT_EDGE` compartido por `_TITLE_MUTATION_RE` y `_mutation_matches`). El defecto era más amplio: en los 552 títulos de KRAS, 34 entradas ganan su mutación de título y ninguna pierde la que tenía. Prueba real: `G12D` 96 confirmados / 4 descartados de 100 examinados (127 encontrados), 9TBR (puesto 69) confirmado; 9E3S (G12N) sigue descartado; `G12C` 95/5; `Q61H` 12/0.

## 8. Siguientes tareas (en este orden)
1. ~~Corregir 9TBR~~ (hecho 2026-10-20).
2. **Expediente en disco** (`dossiers\P01116\`, en .gitignore): cada herramienta guarda su JSON completo y devuelve al chat solo un resumen compacto con la ruta. `build_report` leerá del disco, no de la conversación.
3. **`get_protein_profile`** (`profile.py`): UniProt (función con evidencia, localización, comentarios de enfermedad, variantes) y Open Targets (asociaciones diana-enfermedad con score, tractabilidad por modalidad). Mantener separadas las dos fuentes. Campos ya validados contra las APIs reales (ver 8.1); `probe_profile.ps1` ya no hace falta.
4. **`build_report`** (Markdown → PDF, ~11 secciones, regla anti-alucinación). Sección de métodos con las limitaciones de la sección 6.

### 8.1 Campos validados para `get_protein_profile` (2026-10-20, KRAS/P01116, llamadas reales desde Claude Code)
Fixtures reales (recortados, valores sin tocar) en `tests/fixtures/`: `uniprot_P01116_profile.json`, `opentargets_ENSG00000133703.json`, `opentargets_versioned_id_null.json`.
- **UniProt** `GET https://rest.uniprot.org/uniprotkb/{acc}.json` (el mismo que usa `uniprot.fetch_entry`).
  - `comments[]` por `commentType`. `FUNCTION`: `texts[].value` + `texts[].evidences[]` (`evidenceCode` ECO, `source` PubMed/UniProtKB, `id`). ECO:0000269 = experimental; ECO:0000250 = "By similarity" (inferida de otra proteína, `source: UniProtKB`); ECO:0000305 = inferencia del curador.
  - `SUBCELLULAR LOCATION` (2 en KRAS, una con `molecule: "Isoform 2B"`): `subcellularLocations[].location/topology/orientation` con `value`, `id` (SL-xxxx) y `evidences`.
  - `DISEASE` (10 en KRAS): `disease.diseaseId` (nombre), `diseaseAccession` (DI-xxxxx), `acronym`, `description`, `diseaseCrossReference` (MIM), `evidences`; `note.texts[]`. **Defecto real: 2 de 10 NO tienen `disease`, solo `note`** (astrocitoma pilocítico y otro): no inventar nombre, usar la nota tal cual.
  - `features[]` `type == "Natural variant"` (35): `location.start.value`, `alternativeSequence.originalSequence` / `alternativeSequences[]` (puede ser inserción: G10 -> "GG"), `description` en texto libre ("in NS3; dbSNP:rs…"; usa los `acronym` de DISEASE), `featureId` (VAR_…), `featureCrossReferences` (dbSNP), `evidences`. `Mutagenesis` (8): misma forma, sin featureId.
  - Ensembl en `uniProtKBCrossReferences` (`database == "Ensembl"`, 5 entradas, una por transcrito): `properties[key=GeneId].value` = `ENSG00000133703.15` **con versión**; deduplicar. También `ChEMBL` (CHEMBL2189121, coincide con `chembl.py`), `MIM`, `Orphanet` (`properties[key=Disease]`).
  - `annotationScore` (5.0), `proteinExistence`, `entryAudit.lastAnnotationUpdateDate` (útil para "fecha del dato").
- **Open Targets** GraphQL `POST https://api.platform.opentargets.org/api/v4/graphql`, `target(ensemblId:)`. Esquema verificado por introspección.
  - **El ID con versión devuelve `target: null` con HTTP 200 y sin `errors`** (igual que un ID inexistente): quitar `.NN` antes de consultar y distinguir "no encontrado" de "ID mal formado".
  - `meta.dataVersion {year month iteration}` = 26/09/null; `apiVersion` 26.9.0 (citar versión en el reporte).
  - `approvedSymbol`, `approvedName`, `biotype`, `proteinIds[] {id source}` (uniprot_swissprot P01116: sirve para comprobar el cruce).
  - `tractability[] {modality label value}`: 28 filas, modalidades SM/AB/PR/OC, incluye las `false`. KRAS SM: Approved Drug true, Structure with Ligand true, High-Quality Ligand true, Druggable Family false, pockets false.
  - `associatedDiseases(page:{index,size})` -> `count` (**2262 en KRAS**), `rows[] {score, disease{id name therapeuticAreas{id name}}, datatypeScores[]{id score}}`, ordenadas por score descendente. `datatypeScores` solo trae los tipos con dato (genetic_association, somatic_mutation, literature, clinical, animal_model, affected_pathway, genetic_literature…).
  - Por score global, el top 4 de KRAS son síndromes germinales (Noonan, CFC); NSCLC aparece en el 5.º puesto. Para el analista de fármacos hará falta decidir cómo presentar oncología vs. germinal (por `therapeuticAreas` o `datatypeScores`).
  - Open Targets solo cubre diana humana: con `organism_id != 9606` debe quedar "sin registro".

## 9. Hoja de ruta posterior
- Fase 2: recuperación de archivos (PDB/FASTA, CIF→PDB) y PyMOL a petición (`.pml` + subprocess). Aquí se migra a paquetes por fase.
- Opcional: validar `max_phase__gte` para encontrar compuestos clínicos sin mecanismo; cruce PDB↔ChEMBL por InChIKey; clave NCBI/PubMed; migración a mcp v2; ClinicalTrials.gov limitado; rama de DrugBank; `lookup_compound`; selectividad y datos negativos de ChEMBL.
- Fila no cubierta por ninguna fuente: disponibilidad comercial / calidad de sonda química.

## 10. Protocolo de trabajo
- Una tarea por sesión (en Claude Code) o un chat por paso (en claude.ai); el rol se indica al inicio.
- Probar siempre contra APIs reales; los tests simulados deben reflejar defectos reales.
- Respuestas breves; cambios como diffs salvo en refactors grandes; no subir `server.py` completo salvo que haga falta.
- Commit por paso, hecho por el usuario; no versionar archivos de respaldo (`server_old.py`), `dossiers/` ni salidas de prueba (`probe.txt`).
- Pruebas: carpeta `tests/`, ejecutadas con `uv run python tests/<archivo>.py`.
