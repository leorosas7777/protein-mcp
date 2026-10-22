# protein-research MCP

Servidor MCP local (Python, `uv`, `mcp[cli]<2`/FastMCP) para Claude Desktop. El usuario dice "revisa la proteína X" y el servidor resuelve la proteína, recupera literatura, estructuras y ligandos, y produce un reporte con hipervínculos. Proyecto de portafolio para una empresa europea. Lector objetivo del reporte: analista de descubrimiento de fármacos.

**Al iniciar una sesión, lee `ESTADO_PROYECTO.md` completo.** Es el estado real del proyecto, las decisiones tomadas y las limitaciones conocidas. No reabras decisiones de la sección 5 sin motivo.

## Reglas que no se negocian
- **Nunca `print()` a stdout** (rompe stdio del MCP). Usa `logging`.
- **Anti-alucinación:** todo enlace se construye solo con IDs devueltos por la fuente; toda afirmación del reporte debe tener fuente recuperada; los campos sin dato quedan en `None` o "sin registro", nunca se rellenan ni se infieren nombres.
- **Informa siempre la selección** ("N de total").
- **Valida contra la API real antes de codificar** un campo nuevo. Un resumen de WebFetch no es evidencia. Los tests simulados deben copiar formas y defectos reales de las APIs.
- **Antes de anotar algo como implementado en `ESTADO_PROYECTO.md`, verifícalo en el código y en una prueba real** (este documento ya describió como hechas funciones que no existían).
- Cuestiona las decisiones del usuario si obstaculizan el éxito del proyecto; no confirmes sesgos.

## Arquitectura
- Paquete `protein_mcp/` (estructura en `ESTADO_PROYECTO.md` sección 3). `protein_mcp/server.py` solo registra herramientas (`@mcp.tool()`) y delega. Los docstrings se quedan ahí (Claude Desktop los lee para usar cada herramienta).
- Un módulo por fuente en `protein_mcp/sources/`: `uniprot.py`, `europepmc.py`, `rcsb.py`, `alphafold.py`, `chembl.py`, `opentargets.py`. `common.py` guarda helpers y no importa nada del proyecto (evita ciclos); `dossier.py` guarda y lee expedientes. Imports absolutos (`from protein_mcp.sources import uniprot`). Un solo servidor para todas las fases; la Fase 2 va en `protein_mcp/files/`.
- Los módulos no se importan entre sí salvo lo indicado en `ESTADO_PROYECTO.md` sección 3.

## Comandos
- Ejecutar el servidor: `uv --directory C:\bioinformatics\protein-mcp run protein-mcp` (tras cambiar `pyproject.toml`: `uv sync`)
- Inspector: `uv run mcp dev protein_mcp/server.py`
- Pruebas offline: `uv run python tests/test_chembl.py` (usar el mismo patrón para las demás).
- Pruebas reales: `uv run python tests/smoke_live.py [--skip-slow | kras pcsk9 egfr notfound tnf]` (APIs reales, expedientes en carpeta temporal; invariantes externas, de regresión y anti-alucinación).
- Revisión: tras cambiar una herramienta, lanzar el agente `revisor` (`.claude/agents/revisor.md`) en contexto limpio antes del commit.
- Para que Claude Desktop recargue el código: Quit completo desde la bandeja y abrir un chat nuevo. `query_used` (resolve_protein) y `selection_reason` (AlphaFold) indican qué versión corre.
- Entorno: Windows, PowerShell. Idioma de trabajo: español.

## Flujo de trabajo
- Una tarea por sesión. Al cerrarla, actualiza `ESTADO_PROYECTO.md` (solo con lo verificado).
- Commits los hace el usuario. No hagas `git push`, no borres archivos, no versiones respaldos (`server_old.py`), `dossiers/`, ni salidas de prueba.
- Si un cambio sale mal y volver a un commit anterior para reconstruir desde ahí es más limpio que seguir parchando, dilo explícitamente: qué commit, por qué, y qué cambios sin commitear se perderían. No ejecutes `git reset`/`git checkout`/`git restore` sin aprobación del usuario.
- Respuestas breves; cambios como diffs salvo refactors grandes.
- Si tocas el comportamiento de una herramienta, ajusta su docstring en `protein_mcp/server.py` y su prueba.
