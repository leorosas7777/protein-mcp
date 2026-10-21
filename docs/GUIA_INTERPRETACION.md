# Guía para interpretar el reporte

Para quien no es farmacólogo. Explica qué responde cada sección del reporte, cómo leer sus números y qué **no** se puede concluir de ellos. Los ejemplos son de reportes reales de KRAS (`P01116`) y PCSK9 (`Q8NBP7`) generados el 2026-10-20.

## La idea de fondo

En descubrimiento de fármacos, una **diana** es la proteína sobre la que se quiere actuar. Antes de invertir en ella, un analista se hace cinco preguntas, y cada sección del reporte responde a una:

| Pregunta | Sección |
|---|---|
| ¿Qué hace esta proteína y dónde está? | 1. Target profile (UniProt) |
| ¿Hay pruebas de que tenga que ver con una enfermedad? | 1 (enfermedades, variantes) y 2 (Open Targets) |
| ¿Se puede "atacar" con un fármaco? | 2 (tractability) y 4 (ligandos) |
| ¿Conocemos su forma 3D? | 3. Structures |
| ¿Qué se está publicando? | 5. Literature |

El reporte **no da un veredicto**: reúne la evidencia con su fuente. La decisión la toma una persona.

## Reglas de lectura que valen para todo

- **Cada dato tiene un enlace.** Si algo te sorprende, haz clic y mira la fuente original.
- **"no record"** significa que la fuente no tiene ese dato. No significa "no existe" ni "es cero".
- **"N of total"** (p. ej. "7 of 552 entries shown") indica que es una selección. El reporte nunca es la lista completa.
- **La fecha de consulta** aparece en cada sección. Las bases de datos cambian cada mes.

---

## Contents of the dossier

Tabla de los archivos que alimentan el reporte, con sus parámetros (p. ej. `mutation='G12C'`) y la fecha. Si una sección aparece en **Missing sections**, nadie la consultó; el reporte la deja vacía a propósito.

## 1. Target profile (UniProt)

UniProt es la enciclopedia curada de proteínas: expertos leen la literatura y resumen lo demostrado.

- **reviewed (Swiss-Prot)**: la entrada la revisó un curador humano. *unreviewed (TrEMBL)* es automática y menos fiable.
- **annotation score x/5**: qué tan completa es la ficha. 5/5 (KRAS, muy estudiada) no significa que la proteína sea "mejor" diana.
- **Function**: resumen en texto. Fíjate en la línea **Evidence**:
  - *experimental*: demostrado en un experimento con esta proteína (con su PMID).
  - *by similarity*: deducido de una proteína parecida, no demostrado en esta.
  - *curator inference*: deducción del curador.
- **Subcellular location**: dónde vive la proteína. Importa para el tipo de fármaco: un anticuerpo no entra en la célula, así que solo sirve contra proteínas de la membrana o de fuera. PCSK9 se secreta, por eso funcionan los anticuerpos.
- **Disease involvement**: enfermedades en las que UniProt documenta un papel de la proteína, con enlace a OMIM (catálogo de enfermedades genéticas). Si aparece *"Unnamed in UniProt. Note: ..."*, UniProt solo tiene una nota sin nombre de enfermedad y el reporte la cita tal cual.
- **Natural variants**: cambios en la proteína encontrados en personas. `G12C` significa que en la posición 12 una glicina (G) pasa a cisteína (C). En la descripción:
  - *somatic mutation*: aparece en el tumor, no se hereda (típico en cáncer: KRAS G12C en pulmón).
  - las siglas como NS3, CFC2 o RALD2 son enfermedades hereditarias definidas en la tabla de enfermedades.
- **Mutagenesis**: cambios hechos a propósito en el laboratorio para estudiar la función. No son variantes de pacientes.

**No concluyas:** que una variante "causa" una enfermedad solo porque aparece en la tabla. Lee la descripción y la evidencia.

## 2. Target-disease evidence and tractability (Open Targets)

Open Targets cruza más de 20 fuentes (genética humana, mutaciones en tumores, ensayos clínicos, modelos animales, literatura) y calcula cuánto respaldo hay para cada par diana-enfermedad.

**Tractability flags**: ¿qué tipo de fármaco podría funcionar?

| Modalidad | Qué es | Señales a favor (`yes`) |
|---|---|---|
| small molecule | Pastilla química clásica | *Approved Drug* (ya hay uno), *Structure with Ligand* (hay estructura 3D con un compuesto unido), *High-Quality Pocket* (tiene un hueco donde encajar), *Druggable Family* (pertenece a una familia que suele funcionar) |
| antibody | Proteína que se une desde fuera de la célula | *UniProt loc high conf* (localización accesible), *Approved Drug* |
| PROTAC | Molécula que hace que la célula destruya la proteína | *Small Molecule Binder*, datos de ubiquitinación |
| other clinical | Otras modalidades en clínica | *Approved Drug*, etc. |

Ejemplos reales:
- **KRAS**: *Approved Drug = yes* pero *Druggable Family = no* y pocket = no. Es la historia de KRAS, considerada "imposible de atacar" durante décadas hasta que llegaron los inhibidores de G12C.
- **PCSK9**: *antibody → Approved Drug = yes* (evolocumab, alirocumab).

Son indicadores, **no un veredicto** de que se pueda hacer un fármaco.

**Top associated diseases**:
- **Score (0 a 1)**: cuánta evidencia combinada hay. 0,8 es mucha; 0,1 es poca. **No es una probabilidad ni prueba causalidad.**
- **Evidence types**: de dónde viene el score.
  - *genetic association*: genética humana (variantes asociadas a la enfermedad). Suele ser la evidencia más valiosa.
  - *somatic mutation*: mutaciones en tumores.
  - *clinical*: hay fármacos en ensayos clínicos para esa enfermedad.
  - *animal model*: modelos en animales.
  - *literature*: menciones en artículos (la más débil por sí sola).
  - *affected pathway*: la proteína participa en una vía alterada en la enfermedad.
- **Top overall vs. Top oncology**: el ranking general puede estar dominado por enfermedades hereditarias raras. En KRAS, Noonan (0,834) sale antes que el cáncer de pulmón (0,804). El bloque de oncología filtra por la clasificación de Open Targets, que incluye algunos síndromes.

**No concluyas:** que la enfermedad con mayor score es "la indicación" del fármaco. Mira el tipo de evidencia.

## 3. Structures

- **Experimental (RCSB PDB)**: estructuras 3D medidas en el laboratorio.
  - **Resolución (Å)**: menor es mejor. Menos de 2 Å es muy buena; más de 3 Å, borrosa. No garantiza que el compuesto unido esté bien resuelto.
  - **Method**: *X-RAY DIFFRACTION* (cristal), *ELECTRON MICROSCOPY* (crio-EM), *SOLUTION NMR*.
  - **Mutations**: `field` es el campo oficial; `title` indica que solo se menciona en el título. *"none listed"* **no** prueba que sea la proteína normal.
  - **Ligands**: compuestos unidos en esa estructura (el código de 3 letras es del PDB). Son la prueba física de que algo se une ahí.
  - **Nucl.**: GDP/GTP o análogos (GNP, GCP). En KRAS indican si está "apagada" (GDP) o "encendida" (análogo de GTP).
  - Con filtro de mutación, se muestran las confirmadas y una muestra de las descartadas.
- **Predicted (AlphaFold DB)**: estructura predicha por IA.
  - **pLDDT**: confianza de la predicción. Más de 90 es muy alta, entre 70 y 90 es fiable, por debajo de 50 no fiable (suelen ser zonas flexibles).
  - Es **una sola forma estática, sin fármacos ni iones**. No sirve para decir cómo se une un compuesto.

## 4. Ligands and bioactivity (ChEMBL)

ChEMBL recopila mediciones de potencia publicadas: cuánto compuesto hace falta para actuar sobre la proteína.

- **Max phase**: lo más lejos que ha llegado el compuesto. **4** = aprobado (con año); **3/2/1** = fase de ensayo clínico; *none recorded* = investigación.
- **Potency**: `G12C · IC50 · cellular: 5 nM (median 10, n=3)` significa:
  - **IC50**: concentración que inhibe el 50 % de la actividad. **Ki/Kd**: afinidad de unión. En todos, **menor = más potente**.
  - **nM** (nanomolar): 1 nM es muy potente, 1.000 nM (1 µM) es débil.
  - **binding** (proteína aislada) vs **cellular** (en células). El valor celular suele ser peor pero es más realista. El reporte nunca mezcla ambos.
  - **n**: cuántas mediciones. Un valor con n=1 es una sola medición; desconfía.
- **pChEMBL (best)**: la potencia en escala logarítmica. **9 = 1 nM, 8 = 10 nM, 7 = 100 nM, 6 = 1 µM**. Más alto es más potente.
- **no name in ChEMBL**: la mayoría de compuestos de investigación no tiene nombre. El reporte no lo inventa.
- **Orden**: primero los que llegaron a clínica, luego por potencia. **Potencia no es selectividad** (que no afecte a otras proteínas) **ni que sea un buen fármaco** (absorción, toxicidad).

⚠️ **Punto ciego conocido**: los fármacos que son **anticuerpos** casi nunca tienen IC50 en ChEMBL. En PCSK9, evolocumab y alirocumab (aprobados en 2015) solo aparecen como IDs en la nota *"Compounds with a recorded mechanism but no qualifying activity"*. Lee esa nota y la tractabilidad de Open Targets antes de concluir que "no hay fármacos".

## 5. Literature (Europe PMC)

Artículos que mencionan la proteína en título o resumen. "10 of 25779" significa que es una muestra. *Cited by* indica cuántas veces se citó; los artículos recientes tienen 0 aunque sean buenos. *PREPRINT* = no revisado por pares.

## 6. Methods and limitations

Cómo se seleccionó cada dato y sus advertencias, tal como las declaran las herramientas. Léelo antes de citar el reporte.

---

## Ejemplo: leer KRAS en dos minutos

1. **Función**: interruptor GDP/GTP de la vía Ras-MAPK (proliferación celular). La evidencia es experimental.
2. **Enfermedad**: variantes somáticas G12/G13/Q61 en cánceres de pulmón, colon y páncreas; variantes hereditarias en síndromes (Noonan, CFC).
3. **Tratabilidad**: molécula pequeña aprobada, pero no pertenece a una familia "fácil" y no tiene un bolsillo clásico. Es una diana difícil que se logró atacar.
4. **Estructura**: más de 550 estructuras en el PDB, muchas con inhibidores unidos a G12C/G12D. AlphaFold muy confiable (pLDDT 91,5).
5. **Ligandos (G12C)**: adagrasib y sotorasib aprobados (2022 y 2021), glecirasib en fase 2, cientos de compuestos de investigación con potencia nanomolar.

Conclusión que **sí** permite el reporte: hay evidencia fuerte de que KRAS mutada impulsa ciertos cánceres, y fármacos aprobados para G12C. Conclusión que **no** permite: que un compuesto concreto sea mejor que otro, o que G12D tenga la misma tratabilidad que G12C.

## Glosario rápido

| Término | Significado |
|---|---|
| Accession (P01116) | Código único de la proteína en UniProt |
| Diana (target) | Proteína sobre la que actúa un fármaco |
| Ligando | Molécula que se une a la proteína |
| IC50 / Ki / Kd | Medidas de potencia o afinidad; menor es más potente |
| pChEMBL | −log10 de la potencia molar; 9 = 1 nM |
| nM, µM | Nanomolar, micromolar (1 µM = 1.000 nM) |
| Germinal vs somática | Heredada vs adquirida (p. ej. en el tumor) |
| Å (angstrom) | Unidad de resolución; menor = más nítido |
| pLDDT | Confianza de AlphaFold (0-100) |
| Fase clínica | 1 seguridad, 2 eficacia inicial, 3 eficacia a gran escala, 4 aprobado |
| PMID | Identificador de un artículo en PubMed |
