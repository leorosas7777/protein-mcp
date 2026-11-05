# Guide to interpreting the report

For readers who are not pharmacologists. It explains what each section of the report answers, how to read its numbers, and what **cannot** be concluded from them. The examples come from real KRAS (`P01116`) and PCSK9 (`Q8NBP7`) reports generated on 2026-10-20, except section 6, which uses BRAF (`P15056`).

## The underlying idea

In drug discovery, a **target** is the protein a drug is meant to act on. Before investing in one, an analyst asks five questions, and each section of the report answers one:

| Question | Section |
|---|---|
| What does this protein do and where is it? | 1. Target profile (UniProt) |
| Is there evidence that it is involved in a disease? | 1 (diseases, variants) and 2 (Open Targets) |
| Can a drug act on it? | 2 (tractability) and 4 (ligands) |
| Do we know its 3D shape? | 3. Structures |
| What is being published? | 5. Literature |

The report **gives no verdict**: it gathers the evidence with its source, and a person makes the decision.

## Reading rules that apply throughout

- **Every value has a link.** If something surprises you, click it and check the original source.
- **"no record"** means the source does not have that value. It does not mean "does not exist" or "is zero".
- **"N of total"** (e.g. "7 of 552 entries shown") indicates a selection. The report is never the complete list.
- **The retrieval date** appears in each section. Databases change every month.

---

## Contents of the dossier

A table of the files that feed the report, with their parameters (e.g. `mutation='G12C'`) and date. If a section appears under **Missing sections**, nobody queried it, and the report leaves it empty on purpose.

## 1. Target profile (UniProt)

UniProt is the curated encyclopedia of proteins: experts read the literature and summarize what has been shown.

- **reviewed (Swiss-Prot)**: a human curator reviewed the entry. *unreviewed (TrEMBL)* is automatic and less reliable.
- **annotation score x/5**: how complete the entry is. A score of 5/5 (KRAS, heavily studied) does not make the protein a "better" target.
- **Function**: a text summary. Pay attention to the **Evidence** line:
  - *experimental*: shown in an experiment with this protein (with its PMID).
  - *by similarity*: inferred from a similar protein, not shown in this one.
  - *curator inference*: deduced by the curator.
- **Subcellular location**: where the protein lives. This matters for the type of drug: an antibody does not enter the cell, so it only works against membrane or extracellular proteins. PCSK9 is secreted, which is why antibodies work.
- **Disease involvement**: diseases in which UniProt documents a role for the protein, with a link to OMIM (a catalogue of genetic diseases). If you see *"Unnamed in UniProt. Note: ..."*, UniProt only has a note without a disease name, and the report quotes it as it is.
- **Natural variants**: changes in the protein found in people. `G12C` means that position 12 changes from glycine (G) to cysteine (C). In the description:
  - *somatic mutation*: appears in the tumour and is not inherited (typical in cancer: KRAS G12C in lung cancer).
  - abbreviations such as NS3, CFC2 or RALD2 are hereditary diseases defined in the disease table.
- **Mutagenesis**: changes made on purpose in the laboratory to study function. They are not patient variants.

**Do not conclude** that a variant "causes" a disease just because it appears in the table. Read the description and the evidence.

## 2. Target-disease evidence and tractability (Open Targets)

Open Targets combines more than 20 sources (human genetics, tumour mutations, clinical trials, animal models, literature) and calculates how much support there is for each target-disease pair.

**Tractability flags**: what type of drug could work?

| Modality | What it is | Signals in favour (`yes`) |
|---|---|---|
| small molecule | Classic chemical drug | *Approved Drug* (one already exists), *Structure with Ligand* (a 3D structure with a bound compound exists), *High-Quality Pocket* (it has a cavity to fit into), *Druggable Family* (it belongs to a family that often works) |
| antibody | Protein that binds from outside the cell | *UniProt loc high conf* (accessible location), *Approved Drug* |
| PROTAC | Molecule that makes the cell destroy the protein | *Small Molecule Binder*, ubiquitination data |
| other clinical | Other modalities in the clinic | *Approved Drug*, etc. |

Real examples:
- **KRAS**: *Approved Drug = yes* but *Druggable Family = no* and pocket = no. This is the story of KRAS, considered "undruggable" for decades until G12C inhibitors arrived.
- **PCSK9**: *antibody → Approved Drug = yes* (evolocumab, alirocumab).

These are indicators, **not a verdict** that a drug can be made.

**Top associated diseases**:
- **Score (0 to 1)**: how much combined evidence exists. 0.8 is a lot; 0.1 is little. **It is not a probability and does not prove causality.**
- **Evidence types**: where the score comes from.
  - *genetic association*: human genetics (variants associated with the disease). Usually the most valuable evidence.
  - *somatic mutation*: mutations in tumours.
  - *clinical*: there are drugs in clinical trials for that disease.
  - *animal model*: models in animals.
  - *literature*: mentions in papers (the weakest on its own).
  - *affected pathway*: the protein takes part in a pathway altered in the disease.
- **Top overall vs. Top oncology**: the overall ranking can be dominated by rare hereditary diseases. For KRAS, Noonan syndrome (0.834) comes before lung cancer (0.804). The oncology block filters by Open Targets' classification, which includes some syndromes.

**Do not conclude** that the disease with the highest score is "the indication" for a drug. Look at the type of evidence.

## 3. Structures

- **Experimental (RCSB PDB)**: 3D structures measured in the laboratory.
  - **Resolution (Å)**: lower is better. Below 2 Å is very good; above 3 Å is blurry. It does not guarantee that the bound compound is well resolved.
  - **Method**: *X-RAY DIFFRACTION* (crystal), *ELECTRON MICROSCOPY* (cryo-EM), *SOLUTION NMR*.
  - **Mutations**: `field` is the official field; `title` means it is only mentioned in the title. *"none listed"* does **not** prove it is the normal protein.
  - **Ligands**: compounds bound in that structure (the 3-letter code comes from the PDB). They are physical proof that something binds there.
  - **Nucl.**: GDP/GTP or analogues (GNP, GCP). In KRAS, they indicate whether it is "off" (GDP) or "on" (GTP analogue).
  - With a mutation filter, the confirmed entries are shown along with a sample of the discarded ones.
- **Predicted (AlphaFold DB)**: a structure predicted by AI.
  - **pLDDT**: prediction confidence. Above 90 is very high, 70 to 90 is reliable, and below 50 is unreliable (usually flexible regions).
  - It is **a single static shape, with no drugs or ions**. It cannot tell you how a compound binds.

## 4. Ligands and bioactivity (ChEMBL)

ChEMBL collects published potency measurements: how much compound it takes to act on the protein.

- **Max phase**: the furthest the compound has reached. **4** = approved (with year); **3/2/1** = clinical trial phase; **early phase 1 (0.5)**; **clinical, phase unknown (ChEMBL -1)** = ChEMBL knows it reached the clinic but not which phase; *none recorded* = research.
- **molecule record not retrieved**: the query for the molecule record failed, so name, type and phase are unknown. It does not mean the compound has no name.
- **Potency**: `G12C · IC50 · cellular: 5 nM (median 10, n=3)` means:
  - **IC50**: the concentration that inhibits 50% of the activity. **Ki/Kd**: binding affinity. In all of them, **lower = more potent**.
  - **nM** (nanomolar): 1 nM is very potent; 1,000 nM (1 µM) is weak.
  - **binding** (isolated protein) vs **cellular** (in cells). The cellular value is usually worse but more realistic. The report never mixes the two.
  - **n**: how many measurements. A value with n=1 is a single measurement; be cautious.
- **pChEMBL (best)**: potency on a logarithmic scale. **9 = 1 nM, 8 = 10 nM, 7 = 100 nM, 6 = 1 µM**. Higher is more potent.
- **no name in ChEMBL**: most research compounds have no name. The report does not invent one.
- **Order**: first those that reached the clinic, then by potency. **Potency is not selectivity** (not affecting other proteins) **and does not make a good drug** (absorption, toxicity).

**Antibodies and proteins** almost never have an IC50 in ChEMBL, so they do not enter the potency table. They appear in the table *"Drugs and clinical candidates with a recorded mechanism but no qualifying IC50/Ki/Kd"* with name, type and phase (for PCSK9: evolocumab and alirocumab, approved in 2015). Read it before concluding that "there are no drugs".

Two more tables may appear:
- *"...qualifying IC50/Ki/Kd only for variants other than 'X'"*: the drug has data, but none valid for the requested variant (e.g. lapatinib with EGFR C797S). This does not mean "untested": a measurement may have been excluded by the filters (censored value, low confidence, label conflict).
- *"...could not be checked because their activity download failed"*: the download failed, so it is unknown whether they have potency data.

If **"Warning: the ChEMBL mechanism lookup failed"** appears, the mechanism query failed and drugs without IC50 (antibodies) were not checked. Their absence does not mean they do not exist; repeat the query.

For heavily studied targets (EGFR), the *Activities* line says *truncated*: only the most potent activities were examined (down to the lowest pChEMBL shown, and ties at that value may be partly left out), plus those of drugs with a mechanism (if that extra download failed or was cut off, the line says so and the drug statistics are partial too). That is why the selection says *"At least N ligands"*: it is a minimum, not ChEMBL's total. For other compounds, *n* and the median only count those measurements.

## 5. Literature (Europe PMC)

Papers that mention the protein in the title or abstract. "10 of 25779" means it is a sample. *Cited by* shows how many times a paper was cited; recent papers have 0 even if they are good. *PREPRINT* = not peer reviewed.

## 6. Structure preparation (Phase 2)

This section appears only if `prepare_structures` was run. It describes the files prepared for in-silico experiments.

- **Structures in the manifest**: resolution, R-free (*not reported* when the entry has no value; lower is better), and the cleaned receptor file. Per chain, it gives the modelled UniProt range, the missing residues and what is actually modelled at the variant position.
  - **Variant position** is read per chain. In BRAF 3OG7, position 600 is not modelled in chain A and is GLU in chain B: the title's "V600E" does not hold for every chain.
  - **Mutations** can include solubility mutations added to make the protein crystallize. A structure with them is not the native sequence.
- **Ligand pockets**: target residues within 5 Å of each ligand instance, with a box (ligand + 5 Å per side) for docking programs.
  - *missing ... next to the pocket*: a stretch without coordinates touches the pocket, so the pocket is incomplete.
  - *spans several target chains*: the ligand touches more than one target chain, so the pocket and the box cover all of them and were not split into one box per site.
  - *large box*: one side of the box is longer than 30 Å. A big box slows docking down and lets a compound drift away from the real site.
  - **AlphaFold mean pLDDT**: AlphaFold's confidence at the same UniProt positions. It says nothing about whether the model matches the ligand-bound conformation.
- **Protein-protein interfaces** (only if `analyze_interfaces` was run): partner chains, assemblies according to RCSB and to PDBe PISA, and target residues within 5 Å of another chain.
  - **css** is PISA's score for how significant an interface is (0 to 1). The sources can disagree, and the report shows both without choosing.
  - A contact in a crystal does not prove a contact in solution, and vice versa.
- **Structure chosen by the user**: the structure, chain and pocket recorded with `select_structure`, with the user's reason in quotation marks. "No structure chosen" means nobody has chosen one yet.

**Do not conclude** that the first structure is the best one for docking. The report does not measure conformation or fit to another ligand.

## 7. Methods and limitations

How each value was selected and its caveats, as stated by the tools. Read this before citing the report.

---

## Example: reading KRAS in two minutes

1. **Function**: GDP/GTP switch of the Ras-MAPK pathway (cell proliferation). The evidence is experimental.
2. **Disease**: somatic G12/G13/Q61 variants in lung, colon and pancreatic cancers; hereditary variants in syndromes (Noonan, CFC).
3. **Tractability**: approved small molecule, but it does not belong to an "easy" family and has no classic pocket. It is a hard target that was eventually drugged.
4. **Structure**: more than 550 structures in the PDB, many with inhibitors bound to G12C/G12D. AlphaFold is highly confident (pLDDT 91.5).
5. **Ligands (G12C)**: adagrasib and sotorasib approved (2022 and 2021), glecirasib in phase 2, hundreds of research compounds with nanomolar potency.

A conclusion the report **does** support: there is strong evidence that mutated KRAS drives certain cancers, and there are approved drugs for G12C. A conclusion it does **not** support: that one specific compound is better than another, or that G12D has the same tractability as G12C.

## Quick glossary

| Term | Meaning |
|---|---|
| Accession (P01116) | Unique code of the protein in UniProt |
| Target | Protein a drug acts on |
| Ligand | Molecule that binds to the protein |
| IC50 / Ki / Kd | Measures of potency or affinity; lower is more potent |
| pChEMBL | −log10 of molar potency; 9 = 1 nM |
| nM, µM | Nanomolar, micromolar (1 µM = 1,000 nM) |
| Germline vs somatic | Inherited vs acquired (e.g. in the tumour) |
| Å (angstrom) | Unit of resolution; lower = sharper |
| pLDDT | AlphaFold confidence (0-100) |
| Clinical phase | 1 safety, 2 initial efficacy, 3 large-scale efficacy, 4 approved |
| PMID | Identifier of a paper in PubMed |
| Pocket | Target residues around a bound ligand |
| Box | Search region for docking, centred on the ligand |
