# In-silico experiment record: P15056

<small>Generated 2026-11-15T05:01:17Z from prep/manifest.json (created 2026-10-28T19:25:20Z by prepare_structures). Variant filter: V600E. UniProt sequence used for the numbering: >sp\|P15056\|BRAF_HUMAN Serine/threonine-protein kinase B-raf OS=Homo sapiens OX=9606 GN=BRAF PE=1 SV=4</small>

<small>Results of in-silico experiments (calculated docking scores and RMSD values; ADME records, retrieved measurements and calculated properties of the ligand), each with the provenance of the structure it came from, quoted from the manifest. No value in this document is a prediction: nothing here comes from a trained model. Measured potency against the target (ChEMBL activities) is not part of this document; it is in the research report, report_P15056.pdf. Interpretation is not part of either.</small>

## 1. Molecular docking (AutoDock Vina)

<small>2 docking runs recorded, on 2 of the 10 pockets computed in the manifest (2 of 7 structures). Re-docking of the crystal ligand of a pocket into its own receptor: one ligand per run. The docking results in section 1 (Vina scores and RMSD values) are CALCULATED by AutoDock Vina and RDKit, as recorded in the manifest. The provenance lines quote other kinds of value: the resolution is the experimental value deposited for the PDB entry, the box was computed when the structure was prepared, and run times are measurements of this machine. A docking score is not a measured affinity and is not comparable with the IC50/Ki/Kd values that ChEMBL reports for the same target, which stay in the research report. The poses of one run are listed in the order Vina returned them; no compound is ranked against another here.</small>

| PDB | Pocket | Run | Recorded (UTC) | Written by | Started from |
|---|---|---|---|---|---|
| 8C7Y | TXV U | TXV_regenerated | 2026-11-07T01:35:58Z | dock_ligand | arguments |
| 3OG7 | 032 C | 032_regenerated | 2026-11-06T01:20:59Z | dock_ligand | arguments |

### 8C7Y · pocket TXV U · run TXV_regenerated

#### Provenance

- [8C7Y](https://www.rcsb.org/structure/8C7Y) · X-RAY DIFFRACTION · resolution 1.65 Å · Crystal structure of BRAF V600E in complex with a hybrid compound 6
- mmCIF raw/8C7Y.cif, SHA-256 b9ee9b36bf67b4a2d68beaa74aeeec9e03b4978e4487b5f8afeba82492673195 (retrieved 2026-10-28T19:25:14Z)
- Receptor clean/8C7Y_receptor.pdb, chains A, B
- Pocket TXV U: box 21.8 × 16.27 × 22 Å, centre (0.02, 110.789, 15.175); frame: coordinates of the original mmCIF (unchanged in the cleaned receptor)
- Crystal ligand pose ligands/8C7Y_TXV_U.sdf, SHA-256 60df7b58fd9e0fc8242f32dbc30ea6dd1817922418a14c577a5249c71b6d9112
- Chosen by the user: 8C7Y chain A, pocket TXV U; reason recorded by select_structure: "I want to work on the structure with the most recent inhibitor" (recorded 2026-11-15T05:01:17Z); this run was started with explicit arguments
- Construct mutations recorded by get_structures: title: V600E

#### Ligand

- TXV (instance U, auth chain A); wwPDB synonyms: no wwPDB synonym
- Chemical name: ~{N}-[3-[(5-chloranyl-1~{H}-pyrrolo[2,3-b]pyridin-3-yl)carbonyl]-2,4-bis(fluoranyl)phenyl]-3-(2-cyanopropan-2-yl)benzamide
- SMILES: CC(C)(C#N)c1cccc(C(=O)Nc2ccc(F)c(C(=O)c3c[nH]c4ncc(Cl)cc34)c2F)c1; heavy atoms: 34; rotatable torsions: 5
- Source: crystal pose of this pocket; conformer 'regenerated'; RDKit ETKDGv3 (randomSeed 42) from the crystal ligand's atoms and bonds, then MMFF94 minimisation (not converged in the default iterations)
- Translation applied to the conformer: (0.02, 110.789, 15.073) Å
- Crystal coordinates discarded; the conformer was moved to the box centre and the translation recorded.

#### Receptor preparation

- Meeko mk_prepare_receptor --read_pdb --delete_bad_res_from_box_radius 8; exit code 0; 6.7 s
- Residues matching no Meeko template: A:607, A:608, A:609
- Deleted by Meeko (outside the box, so the pocket is unchanged): A:607 (17.487 Å outside the box), A:608 (14.981 Å outside the box), A:609 (11.859 Å outside the box)
- Refused by Meeko (too close to the box, not deleted): none
- Residue identifiers: chain:number as in the cleaned receptor PDB (auth numbering)

#### Checks run before Vina was started

- Ligand atoms inside the box: 51 of 51
- Ligand centroid inside the box: True
- Receptor atoms within 8 Å of the box centre: 93; closest at 2.78 Å
- The conformer was moved to the box centre, so the ligand checks confirm the move; the receptor count is the check that can fail.

#### AutoDock Vina run

- Vina 1.2.7 (executable SHA-256 e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5)
- Parameters: exhaustiveness 32, seed 42, num_modes 9, energy_range 3, min_rmsd 1, scoring vina; CPUs: all (Vina default)
- Search box: 21.8 × 16.27 × 22 Å, centre (0.02, 110.789, 15.175)
- Wall-clock time: 169.8 s
- CPU time: 1108.4 s on 12 logical CPUs (6.53 s of CPU per second of wall clock); a measurement of this machine, not of the docking

#### Poses (calculated)

<small>9 poses returned (num_modes 9). Both columns are calculated values of this run; the RMSD is against the crystal pose of this same pocket (method: rdMolAlign.CalcRMS, heavy atoms, symmetry-corrected, in place (no superposition)).</small>

| Mode | Vina score (kcal/mol, calculated) | RMSD to the crystal pose (Å, calculated) |
|---|---|---|
| 1 | -13.15 | 1.0530 |
| 2 | -12.28 | 2.3663 |
| 3 | -11.51 | 2.3975 |
| 4 | -11.08 | 9.7305 |
| 5 | -10.99 | 11.0949 |
| 6 | -10.68 | 1.9501 |
| 7 | -10.52 | 11.8667 |
| 8 | -10.31 | 9.1127 |
| 9 | -10.20 | 10.0447 |

- Top-pose RMSD: 1.0530 Å; criterion <= 2 Å; within the criterion: True
- Pre-registered control run: True
- RMSD method: rdMolAlign.CalcRMS, heavy atoms, symmetry-corrected, in place (no superposition)
- Negative control of this run, measured with the same function as the poses: the crystal pose against itself gave 0 Å and the crystal pose moved 5 Å gave 5 Å; passed: True
- Scope: One re-docking of one ligand into one structure (8C7Y TXV U); not the accuracy of the method.

#### Files written (paths relative to the prep folder)

| Name | Path | SHA-256 |
|---|---|---|
| ligand_TXV_regenerated.pdbqt | dock/8C7Y_TXV_U/ligand_TXV_regenerated.pdbqt | b34d22de82f05ab18439842dffe195a235b9661fd239577763fd465de6fa639e |
| poses_TXV_regenerated.pdbqt | dock/8C7Y_TXV_U/poses_TXV_regenerated.pdbqt | 2f8fb134b64f3f585329883806795fe09e77decfc0ed698778d32034284a84ed |
| poses_TXV_regenerated.sdf | dock/8C7Y_TXV_U/poses_TXV_regenerated.sdf | a4f0505f979d63fcb9e98c8c78d3c941ac4fa490443690b7ddcb29d8c3c69217 |
| receptor.pdbqt | dock/8C7Y_TXV_U/receptor.pdbqt | ba4ac612c9388e58c169c59f21696bf757758cc6f8b3c6d6998a459e394d12b3 |
| receptor_meeko.log | dock/8C7Y_TXV_U/receptor_meeko.log | 6ecd45295a62e45649afe4cc423cb418bc7a6c600c17ac7e94c71901ec0ba986 |
| vina_TXV_regenerated.log | dock/8C7Y_TXV_U/vina_TXV_regenerated.log | 4ab041e22db15ddedb1570b25f84db0674c0f1b57bc47f4c98d20668ea4f0c37 |

<small>calculated: Calculated by AutoDock Vina (scoring function 'vina'): a docking score is not a measured affinity and is not comparable with ChEMBL IC50/Ki/Kd values of the same dossier.</small>

### 3OG7 · pocket 032 C · run 032_regenerated

#### Provenance

- [3OG7](https://www.rcsb.org/structure/3OG7) · X-RAY DIFFRACTION · resolution 2.45 Å · B-Raf Kinase V600E oncogenic mutant in complex with PLX4032
- mmCIF raw/3OG7.cif, SHA-256 f59e2f517398236837cdc228f46fe2e2b29a387977004e4a3b4e453c6193020c (retrieved 2026-10-28T19:25:14Z)
- Receptor clean/3OG7_receptor.pdb, chains A, B
- Pocket 032 C: box 26.3 × 15.13 × 16.4 Å, centre (2.643, -2.28, -19.403); frame: coordinates of the original mmCIF (unchanged in the cleaned receptor)
- Crystal ligand pose ligands/3OG7_032_C.sdf, SHA-256 21501fac3aac6535e3e97a564bd1fede32a8466121d642a2a8cc3dede502ceff
- Not the structure and pocket recorded by select_structure (8C7Y chain A, pocket TXV U): the reason the user gave there is not quoted here.
- Construct mutations recorded by get_structures: field: A688R, F667E, I544A, I551K, K630S, L588N, L706S, L716E, Q562R, Q709R, S713E, S720E, Y673S; title: V600E

#### Ligand

- 032 (instance C, auth chain A); wwPDB synonyms: Vemurafenib, PLX4032
- Chemical name: N-(3-{[5-(4-chlorophenyl)-1H-pyrrolo[2,3-b]pyridin-3-yl]carbonyl}-2,4-difluorophenyl)propane-1-sulfonamide
- SMILES: CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(-c4ccc(Cl)cc4)cc23)c1F; heavy atoms: 33; rotatable torsions: 7
- Source: crystal pose of this pocket; conformer 'regenerated'; RDKit ETKDGv3 (randomSeed 42) from the crystal ligand's atoms and bonds, then MMFF94 minimisation (not converged in the default iterations)
- Translation applied to the conformer: (2.643, -2.28, -19.327) Å
- Crystal coordinates discarded; the conformer was moved to the box centre and the translation recorded.

#### Receptor preparation

- Meeko mk_prepare_receptor --read_pdb --delete_bad_res_from_box_radius 8; exit code 0; 5.4 s
- Residues matching no Meeko template: none reported
- Deleted by Meeko (outside the box, so the pocket is unchanged): none
- Refused by Meeko (too close to the box, not deleted): none
- Residue identifiers: chain:number as in the cleaned receptor PDB (auth numbering)

#### Checks run before Vina was started

- Ligand atoms inside the box: 51 of 51
- Ligand centroid inside the box: True
- Receptor atoms within 8 Å of the box centre: 97; closest at 4.06 Å
- The conformer was moved to the box centre, so the ligand checks confirm the move; the receptor count is the check that can fail.

#### AutoDock Vina run

- Vina 1.2.7 (executable SHA-256 e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5)
- Parameters: exhaustiveness 32, seed 42, num_modes 9, energy_range 3, min_rmsd 1, scoring vina; CPUs: all (Vina default)
- Search box: 26.3 × 15.13 × 16.4 Å, centre (2.643, -2.28, -19.403)
- Wall-clock time: 214.8 s
- CPU time: 1532 s on 12 logical CPUs (7.13 s of CPU per second of wall clock); a measurement of this machine, not of the docking

#### Poses (calculated)

<small>9 poses returned (num_modes 9). Both columns are calculated values of this run; the RMSD is against the crystal pose of this same pocket (method: rdMolAlign.CalcRMS, heavy atoms, symmetry-corrected, in place (no superposition)).</small>

| Mode | Vina score (kcal/mol, calculated) | RMSD to the crystal pose (Å, calculated) |
|---|---|---|
| 1 | -11.35 | 0.5396 |
| 2 | -10.96 | 2.8565 |
| 3 | -10.94 | 6.0683 |
| 4 | -10.61 | 8.3492 |
| 5 | -10.61 | 8.7793 |
| 6 | -10.52 | 2.1040 |
| 7 | -10.51 | 6.6953 |
| 8 | -10.44 | 9.3754 |
| 9 | -10.35 | 8.2993 |

- Top-pose RMSD: 0.5396 Å; criterion <= 2 Å; within the criterion: True
- Pre-registered control run: True
- RMSD method: rdMolAlign.CalcRMS, heavy atoms, symmetry-corrected, in place (no superposition)
- Negative control of this run, measured with the same function as the poses: the crystal pose against itself gave 0 Å and the crystal pose moved 5 Å gave 5 Å; passed: True
- Scope: One re-docking of one ligand into one structure (3OG7 032 C); not the accuracy of the method.

#### Files written (paths relative to the prep folder)

| Name | Path | SHA-256 |
|---|---|---|
| ligand_032_regenerated.pdbqt | dock/3OG7_032_C/ligand_032_regenerated.pdbqt | a0496e9eb094553c7bd0c9026637014799ec0340d3043bcea261f646a1f702e3 |
| poses_032_regenerated.pdbqt | dock/3OG7_032_C/poses_032_regenerated.pdbqt | ff5cc704469b2f594a672eb45d42b39b8c7cafc2f6b1fce46d93c57c095f0b74 |
| poses_032_regenerated.sdf | dock/3OG7_032_C/poses_032_regenerated.sdf | b708374495ef4eaa5257cfa931007749f7964f3e8fbf590f2cab86d00a5eddf4 |
| receptor.pdbqt | dock/3OG7_032_C/receptor.pdbqt | 988a8c0620c234c5825656f8a8de486e1dd2e36474d953c7005e409b073103cf |
| receptor_meeko.log | dock/3OG7_032_C/receptor_meeko.log | 12a7b0d157af7f90889e3819343ced96e3e3b47f2bc413e96a648d3b241f91ec |
| vina_032_regenerated.log | dock/3OG7_032_C/vina_032_regenerated.log | d1cb21cf090b6f419e3d19f4045bbcf1f730d7238669fb6842999fdd436a8163 |

<small>calculated: Calculated by AutoDock Vina (scoring function 'vina'): a docking score is not a measured affinity and is not comparable with ChEMBL IC50/Ki/Kd values of the same dossier.</small>

## 2. Molecular dynamics

<small>2 pockets with a molecular dynamics record; 2 analysed runs (3OG7 032 C production, 3OG7 032 C equilibration). Section 2 holds CALCULATED values: a molecular dynamics simulation run with NAMD and measured with MDAnalysis by analyze_md, quoted from the manifest. Not one of them is a binding affinity or a free energy, and none is comparable with the docking scores of section 1 or the measured potencies of the research report. The thresholds in the control line and in the offered reading are conventions, not measurements of this system: they were written down before the first production trajectory existed, but after the restrained equilibration of the same system had been measured, so they cannot have been tuned to a production result and could have been to that equilibration.</small>

### 8C7Y · pocket TXV U

#### Configuration (as recorded when the package and the job were written)

| Field | Value | Recorded by |
|---|---|---|
| Starting pose | pose 1 (Vina mode 1) of docking run TXV_regenerated; Vina score -13.154 kcal/mol and RMSD to the crystal pose 1.053 Å, both calculated | md_package |
| Gaps left for CHARMM-GUI to build | B 600-614 (internal) | md_package (gaps read from mmCIF) |
| Simulation job | no md_job recorded for this pocket | — |

<small>Quoted from the md_package and md_job blocks of the manifest. They record what was WRITTEN; whether the cluster ran these inputs unchanged is not recorded there. Each analysed run below names its own NAMD input, and its time axis comes from the trajectory, not from this table.</small>

<small>No trajectory analysis is recorded for this pocket. Run analyze_md on the folder brought back from the cluster with uniprot_id="P15056" to record one.</small>

### 3OG7 · pocket 032 C

#### Configuration (as recorded when the package and the job were written)

| Field | Value | Recorded by |
|---|---|---|
| Starting pose | pose 1 (Vina mode 1) of docking run 032_regenerated; Vina score -11.354 kcal/mol and RMSD to the crystal pose 0.5396 Å, both calculated | md_package |
| Gaps left for CHARMM-GUI to build | A 545-547 (internal), A 597-614 (internal), A 627-630 (internal), B 601-614 (internal), B 627-630 (internal) | md_package (gaps read from mmCIF) |
| Residues built in the system CHARMM-GUI returned | A 545-547, A 597-614, A 627-630, B 601-614, B 627-630; closest to the ligand: A 597-614 at 4.76 Å | md_job (archive checks) |
| CHARMM-GUI job | 9418739369 (folder charmm-gui-9418739369); archive SHA-256 f3160bd96a45fc78019378255ff431282c0176ddc6d968be5aa26ca98549e9b5 | md_job |
| Topology (step3_input.psf) SHA-256 | not recorded by the version of build_md_job that wrote this job | md_job |
| System | box side 114.0 Å; 43655 waters; ions POT 123, CLA 133; 303.15 K; time step 2.0 fs; cutoff 12.0 Å | md_job (archive checks) |
| Equilibration as written | step4_equilibration.inp: minimize 10000 steps, then run 125000 steps | md_job |
| Production as written | step5_production.inp: 10 chunk(s) × 1.0 ns = 10.0 ns; dcdfreq 50000 steps | md_job |
| Inputs edited by build_md_job | none | md_job |
| Cluster job | partition gpu; GPUs gpu:a100:1; time limit 1-00:00:00; NAMD .../NAMD_3.0.3_Linux-x86_64-multicore-CUDA/namd3 (full path in the manifest) | md_job (each value's mark and source are in the manifest) |

<small>Quoted from the md_package and md_job blocks of the manifest. They record what was WRITTEN; whether the cluster ran these inputs unchanged is not recorded there. Each analysed run below names its own NAMD input, and its time axis comes from the trajectory, not from this table.</small>

#### Production run (stage step5; input step5_production.inp)

EQUILIBRATION CONTROL: NOT MET. The full backbone of segment PROA has an absolute drift of 1.3534 Å (\|last-quarter mean minus first-quarter mean\|; signed 1.3534 Å) against a threshold of 0.5 Å. The threshold is a convention, not a measurement, written down before any production trajectory existed but after the restrained equilibration of the same system had been measured. Every figure below is premature.

<small>Recorded 2026-11-15T04:54:58Z by analyze_md. Every figure under this heading comes from the production run only. The word 'production' is the label the recorded job gives this stage: md_job.scheme names 'step5' as production and 'step4' as equilibration (namd/README of the archive).</small>

- Trajectory md/3OG7_032_C/cluster/system/namd/step5_production.dcd: 100 frames, 167995076 bytes, SHA-256 e38ed45e6910471a5f3871e56d195012ea8a9006d04db5179790ee8b877df3c4
- Topology md/3OG7_032_C/cluster/system/namd/step3_input.psf, SHA-256 b619657fe4407af7d90414d6dc5d956373a5d37ea3634032589f7d8fbb7cb9e8; equal to the PSF of the registered package (charmm-gui-9418739369/namd/step3_input.psf inside charmm-gui_D.tgz, whose SHA-256 matches md_job.archive (the md_job block predates the PSF hash))
- 100 of 100 frames analysed; 0 dropped from the start; frame spacing 100.0 ps; span 9.9 ns between the first and the last analysed frame
- Minimization: 1 NAMD input(s) of this stage were read (step5_production.inp) and none declares a minimization, so no leading frame was treated as one
- Restraints: no NAMD input of this stage declares `constraints on` or `fixedAtoms on` (both were looked for), so nothing read here says the atoms were held down
- Held atoms among those measured: ligand heavy atoms: whether the ligand heavy was held is UNKNOWN: no restraint reference file was read; protein backbone: whether the protein backbone was held is UNKNOWN: no restraint reference file was read
- Ligand: 33 heavy of 52 atoms in resname LIG; dropped by mass: H1, H10, H11, H12, H13, H14, H15, H16, H17, H18, H2, H3, H4, H5, H6, H7, H8, H9, LP1
- Pocket: 23 of 544 protein residues (protein residues with any atom within 5.0 A of a ligand heavy atom, in the reference frame): PROA:ILE463, PROA:VAL471, PROA:ALA481, PROA:VAL482, PROA:LYS483, PROA:LEU505, PROA:THR508, PROA:ILE513, PROA:LEU514, PROA:LEU515, PROA:PHE516, PROA:ILE527, PROA:THR529, PROA:GLN530, PROA:TRP531, PROA:CYS532, PROA:SER535, PROA:PHE583, PROA:ILE592, PROA:GLY593, PROA:ASP594, PROA:PHE595, PROA:GLY596
- Built residues in the pocket: 0 of 23 pocket residues were built
- Periodic images: 0 of 100 frames re-imaged; non-rectangular boxes 0; frames without a box 0; frames with the pocket split across the box 0
- Arguments: ligand_resname LIG, pocket_cutoff 5.0, skip_frames -1, reference_frame -1

#### Measured over the production run (calculated)

| Metric | Selection | Min | Max | Mean | First-quarter mean | Last-quarter mean | Drift | Last-quarter range |
|---|---|---|---|---|---|---|---|---|
| Protein backbone RMSD, PROA, full segment (pre-registered control) | 1087 backbone atoms | 0.0000 | 3.9339 | 3.1091 | 2.2517 | 3.6051 | 1.3534 | 0.8843 |
| Protein backbone RMSD, PROA, observed residues only (post hoc, NOT pre-registered) | 987 backbone atoms of the 247 OBSERVED residues of PROA, superposed on themselves; the 25 built residues are left out | 0.0000 | 2.2518 | 1.7565 | 1.5357 | 1.9611 | 0.4254 | 0.4764 |
| Protein backbone RMSD, PROA, residues built by CHARMM-GUI only (post hoc, NOT pre-registered) | 100 backbone atoms of the 25 BUILT residues of PROA, carried by the superposition on the observed ones (mass-weighted) | 0.0000 | 12.5750 | 9.0531 | 5.8733 | 10.7228 | 4.8495 | 3.7877 |
| Protein backbone RMSD, PROB, full segment (pre-registered control) | 1087 backbone atoms | 0.0000 | 2.5656 | 2.2644 | 1.9185 | 2.4201 | 0.5016 | 0.3089 |
| Protein backbone RMSD, PROB, observed residues only (post hoc, NOT pre-registered) | 1015 backbone atoms of the 254 OBSERVED residues of PROB, superposed on themselves; the 18 built residues are left out | 0.0000 | 2.3391 | 1.9981 | 1.6330 | 2.1027 | 0.4697 | 0.3255 |
| Protein backbone RMSD, PROB, residues built by CHARMM-GUI only (post hoc, NOT pre-registered) | 72 backbone atoms of the 18 BUILT residues of PROB, carried by the superposition on the observed ones (mass-weighted) | 0.0000 | 6.4930 | 4.6938 | 4.3002 | 5.2672 | 0.9669 | 2.6183 |
| Ligand heavy-atom RMSD, fitted on the pocket backbone (never on the ligand) | 33 heavy of 52 atoms in resname LIG; superposed on 92 pocket backbone atoms | 0.0000 | 2.0963 | 1.2092 | 1.0210 | 1.3175 | 0.2965 | 1.0240 |
| Ligand–pocket centre-of-mass distance | centre of mass of 33 ligand heavy atoms against 374 pocket atoms | 3.0965 | 4.0341 | 3.6427 | 3.5594 | 3.7553 | 0.1958 | 0.6814 |
| Ligand–protein hydrogen bonds per frame (a count, not Å) | pairs meeting the criterion below | 0.0000 | 3.0000 | 1.9100 | 1.9600 | 1.9600 | 0.0000 | 2.0000 |

<small>Values in Å except the hydrogen-bond count. Drift = last-quarter mean minus first-quarter mean. Each protein segment is superposed on itself; segments are never fitted together. Reference: frame 0 of this trajectory, the first analysed frame.</small>

<small>The 'observed only' and 'built only' rows are an analysis added on 2026-11-12 AFTER the pre-registered full-backbone control failed on this project's first production run. Its reason (residues CHARMM-GUI built are a guess, not a measurement) does not depend on the result, but the decision to look came after seeing it: it was NOT pre-registered, and the full-segment row stays the control. 'Observed' means present in the structure uploaded to CHARMM-GUI.</small>

#### Hydrogen bonds, production run

<small>3 ligand–protein pair(s) met the criterion in at least one analysed frame. Criterion: acceptor = N or O by mass; donor = the N or O that the hydrogen is BONDED to in the topology; geometry = donor-acceptor distance <= 3.5 A and donor-H-acceptor angle >= 150.0 deg, counted per frame. No partial charge takes part in deciding which atoms are eligible. The share is a PROPORTION OF ANALYSED FRAMES, never a fraction of time: hydrogen bonds break and re-form faster than a frame is stored.</small>

| Donor | Hydrogen | Acceptor | Frames present / analysed | Share of analysed frames |
|---|---|---|---|---|
| HETA:LIG1:N1 | H8 | PROA:ASP594:OD1 | 97 / 100 | 0.97 |
| HETA:LIG1:N2 | H12 | PROA:GLN530:O | 47 / 100 | 0.47 |
| PROA:CYS532:N | HN | HETA:LIG1:N3 | 47 / 100 | 0.47 |

#### Pocket flexibility (RMSF), production run

<small>the CA of 23 of 23 pocket residues. deviation about each atom's own mean position after the same pocket backbone superposition; each value is labelled from the CA atom it was computed from, not from a parallel list. With no apo simulation, a flexible residue cannot be attributed to the ligand.</small>

| Residue | RMSF (Å) |
|---|---|
| PROA:ILE463 | 0.8792 |
| PROA:VAL471 | 0.6020 |
| PROA:ALA481 | 0.3219 |
| PROA:VAL482 | 0.3441 |
| PROA:LYS483 | 0.3742 |
| PROA:LEU505 | 0.5566 |
| PROA:THR508 | 0.5150 |
| PROA:ILE513 | 0.5008 |
| PROA:LEU514 | 0.3699 |
| PROA:LEU515 | 0.3350 |
| PROA:PHE516 | 0.3567 |
| PROA:ILE527 | 0.3740 |
| PROA:THR529 | 0.3018 |
| PROA:GLN530 | 0.3179 |
| PROA:TRP531 | 0.3720 |
| PROA:CYS532 | 0.6427 |
| PROA:SER535 | 0.8036 |
| PROA:PHE583 | 0.5991 |
| PROA:ILE592 | 0.5189 |
| PROA:GLY593 | 0.6242 |
| PROA:ASP594 | 0.4848 |
| PROA:PHE595 | 0.5833 |
| PROA:GLY596 | 1.2151 |

<small>Figure (drawn in the PDF only): Production run, stage step5: drawn from md/3OG7_032_C/analysis_production.json, whose SHA-256 matches the manifest; nothing is recomputed. All 100 analysed frames plotted. Solid lines are full segments, dashed lines the observed residues only (post hoc, not pre-registered). The hydrogen-bond panel counts bonds per stored frame, a sample, not a fraction of time; its count is recorded without times, so it is placed on the ligand series' time axis (same analysed frames, same length), and left out if the lengths differ.</small>

#### Offered reading of the production run (not measured)

<small>READING OFFERED, NOT MEASURED. These lines are an interpretation built from the numbers in this same result and from the thresholds listed here, which are conventions and not measurements of this system. The measurements above stand on their own; this reading does not. Apply your own judgement. Thresholds: protein_plateau_angstrom 0.5, ligand_stable_angstrom 2.0, com_stable_range_angstrom 2.0, hbond_persistent_occupancy 0.5. Source of the thresholds: common conventions, written down with a reason for each on 2026-11-11, before any production trajectory existed but after a restrained equilibration of this project's first simulated system (BRAF with vemurafenib, PDB 3OG7) had been measured, so they were chosen with those numbers in view. No literature is cited because none was consulted. Backbone plateau: drift (last-quarter mean minus first-quarter mean) of at most 0.5 A, of the order of the thermal noise of a protein backbone at 303 K. Ligand RMSD: at most 2.0 A over the last quarter; the same number as this project's docking validation threshold BY CONVENTION, not by derivation: that one compares a calculated pose with a crystal pose, this one asks whether a pose persists over time. Range of the ligand-pocket centre-of-mass distance: at most 2.0 A, to separate a reorientation within the site from a dissociation; the scale comes from the author's earlier simulations of a different system (neuropeptide-GPCR complexes, 10 ns) and is not a value transferable to this one. Hydrogen-bond occupancy: at least 0.5 of the analysed frames, the conventional cut between a contact that holds and an occasional one; it is a sample of frames, not a fraction of time.</small>

- Protein backbone RMSD of segment PROA (the one of 2 with the largest drift) averages 3.6051 A over the last quarter against 2.2517 A over the first (a drift of 1.3534 A, compared with a chosen 0.5 A). That is a drift rather than a plateau, so the protein may not be equilibrated and everything below is premature.
- Split into observed and BUILT residues (built = absent from, or without a full backbone in, CHARMM-GUI's .orig.pdb; observed = present in the structure uploaded to CHARMM-GUI, which is experimental only if that was a deposited structure - not checked here): PROA: the observed residues drift 0.4254 A over 987 atoms (within the chosen 0.5 A; full segment 1.3534 A), with 25 of 272 residues built and left out (545-547, 597-614, 627-630); the built residues alone drift 4.8495 A and average 10.7228 A over the last quarter; PROB: the observed residues drift 0.4697 A over 1015 atoms (within the chosen 0.5 A; full segment 0.5016 A), with 18 of 272 residues built and left out (601-614, 627-630); the built residues alone drift 0.9669 A and average 5.2672 A over the last quarter. For PROA, the segment of the line above, the observed residues are within the chosen figure, so what exceeds it comes from the residues CHARMM-GUI built. THIS SPLIT WAS NOT PRE-REGISTERED: it was added on 2026-11-12, after the pre-registered full-backbone control failed on this project's first production run (the 10 ns production of BRAF with vemurafenib, PDB 3OG7). Its reason - built residues are a guess, not a measurement - does not depend on the result, but the decision to look came after seeing it. The line above remains the pre-registered control.
- Ligand RMSD reaches 2.0963 A at most and averages 1.3175 A over the last quarter (chosen comparison: 2.0 A). The pose held within that figure.
- Centre-of-mass distance stays in a 0.9376 A range (chosen comparison: 2.0 A), mean 3.6427 A. No dissociation is visible in this window.
- 1 of 3 hydrogen-bond pairs held in at least 50% of the analysed frames: HETA:LIG1:N1 -> PROA:ASP594:OD1 (97%).

<small>Not measured: Binding affinity. Nothing here computes a free energy, and none of these metrics is one: a stable pose is not a potent ligand. Bridging waters. Left out of this version on purpose. Whether the pose is correct. With no apo control, a change in the protein cannot be attributed to the ligand.</small>

#### Equilibration run (stage step4; input step4_equilibration.inp)

EQUILIBRATION CONTROL: NO EVIDENCE. The full backbone of segment PROB has an absolute drift of 0.0694 Å (\|last-quarter mean minus first-quarter mean\|; signed 0.0694 Å) against a threshold of 0.5 Å. The threshold is a convention, not a measurement, written down before any production trajectory existed but after the restrained equilibration of the same system had been measured, but the protein backbone was held by the restraints of this stage (2174 of 2174 protein backbone atoms are held in place by the restraint file, so their motion is damped and NOT free. Any figure built on those atoms passes by construction), so a drift within the threshold is what the restraints produce by construction and is not evidence of equilibration. A figure below that concerns held atoms carries no evidence either.

<small>Recorded 2026-11-15T04:55:33Z by analyze_md. Every figure under this heading comes from the equilibration run only. The word 'equilibration' is the label the recorded job gives this stage: md_job.scheme names 'step5' as production and 'step4' as equilibration (namd/README of the archive).</small>

- Trajectory md/3OG7_032_C/cluster/system/namd/step4_equilibration.dcd: 27 frames, 45358872 bytes, SHA-256 e36b3f30ddb4239dcebee827413c913ef34f553c8b0a5530218a51e6aedfac55
- Topology md/3OG7_032_C/cluster/system/namd/step3_input.psf, SHA-256 b619657fe4407af7d90414d6dc5d956373a5d37ea3634032589f7d8fbb7cb9e8; equal to the PSF of the registered package (charmm-gui-9418739369/namd/step3_input.psf inside charmm-gui_D.tgz, whose SHA-256 matches md_job.archive (the md_job block predates the PSF hash))
- 25 of 27 frames analysed; 2 dropped from the start; frame spacing 10.0 ps; span 0.24 ns between the first and the last analysed frame
- Minimization: 2 leading frame(s) come from `minimize` in step4_equilibration.inp and are NOT dynamics: minimization has no time axis. They were dropped.
- Restraints: step4_equilibration.inp declares `constraints` on: the atoms its reference file marks are held in place. On such a trajectory an RMSD of those atoms CANNOT move much, so any control built on it passes by construction and is not evidence of anything. Which atoms were held is in `held`, read from the file itself
- Held atoms among those measured: ligand heavy atoms: 33 of 33 ligand heavy atoms are held in place by the restraint file, so their motion is damped and NOT free. Any figure built on those atoms passes by construction; protein backbone: 2174 of 2174 protein backbone atoms are held in place by the restraint file, so their motion is damped and NOT free. Any figure built on those atoms passes by construction
- Ligand: 33 heavy of 52 atoms in resname LIG; dropped by mass: H1, H10, H11, H12, H13, H14, H15, H16, H17, H18, H2, H3, H4, H5, H6, H7, H8, H9, LP1
- Pocket: 27 of 544 protein residues (protein residues with any atom within 5.0 A of a ligand heavy atom, in the reference frame): PROA:ILE463, PROA:VAL471, PROA:TYR472, PROA:ALA481, PROA:VAL482, PROA:LYS483, PROA:LEU505, PROA:THR508, PROA:ILE513, PROA:LEU514, PROA:LEU515, PROA:PHE516, PROA:ILE527, PROA:THR529, PROA:GLN530, PROA:TRP531, PROA:CYS532, PROA:GLY534, PROA:SER535, PROA:SER536, PROA:PHE583, PROA:ILE592, PROA:GLY593, PROA:ASP594, PROA:PHE595, PROA:GLY596, PROA:LEU597
- Built residues in the pocket: 1 of 27 pocket residues were built, not observed: PROA:LEU597
- Periodic images: 0 of 25 frames re-imaged; non-rectangular boxes 0; frames without a box 0; frames with the pocket split across the box 0
- Arguments: ligand_resname LIG, pocket_cutoff 5.0, skip_frames -1, reference_frame -1, stage step4

#### Measured over the equilibration run (calculated)

| Metric | Selection | Min | Max | Mean | First-quarter mean | Last-quarter mean | Drift | Last-quarter range |
|---|---|---|---|---|---|---|---|---|
| Protein backbone RMSD, PROA, full segment (pre-registered control) | 1087 backbone atoms | 0.0000 | 0.3794 | 0.3498 | 0.2998 | 0.3681 | 0.0683 | 0.0144 |
| Protein backbone RMSD, PROA, observed residues only (post hoc, NOT pre-registered) | 987 backbone atoms of the 247 OBSERVED residues of PROA, superposed on themselves; the 25 built residues are left out | 0.0000 | 0.3739 | 0.3437 | 0.2942 | 0.3602 | 0.0660 | 0.0162 |
| Protein backbone RMSD, PROA, residues built by CHARMM-GUI only (post hoc, NOT pre-registered) | 100 backbone atoms of the 25 BUILT residues of PROA, carried by the superposition on the observed ones (mass-weighted) | 0.0000 | 0.5053 | 0.4089 | 0.3531 | 0.4432 | 0.0901 | 0.0932 |
| Protein backbone RMSD, PROB, full segment (pre-registered control) | 1087 backbone atoms | 0.0000 | 0.3829 | 0.3536 | 0.3028 | 0.3722 | 0.0694 | 0.0082 |
| Protein backbone RMSD, PROB, observed residues only (post hoc, NOT pre-registered) | 1015 backbone atoms of the 254 OBSERVED residues of PROB, superposed on themselves; the 18 built residues are left out | 0.0000 | 0.3758 | 0.3444 | 0.2951 | 0.3626 | 0.0675 | 0.0199 |
| Protein backbone RMSD, PROB, residues built by CHARMM-GUI only (post hoc, NOT pre-registered) | 72 backbone atoms of the 18 BUILT residues of PROB, carried by the superposition on the observed ones (mass-weighted) | 0.0000 | 0.5647 | 0.4650 | 0.3962 | 0.4890 | 0.0928 | 0.1524 |
| Ligand heavy-atom RMSD, fitted on the pocket backbone (never on the ligand) | 33 heavy of 52 atoms in resname LIG; superposed on 108 pocket backbone atoms | 0.0000 | 0.4184 | 0.3385 | 0.3125 | 0.3505 | 0.0381 | 0.1001 |
| Ligand–pocket centre-of-mass distance | centre of mass of 33 ligand heavy atoms against 432 pocket atoms | 2.6135 | 2.9255 | 2.7948 | 2.7834 | 2.8044 | 0.0210 | 0.2151 |
| Ligand–protein hydrogen bonds per frame (a count, not Å) | pairs meeting the criterion below | 1.0000 | 3.0000 | 2.3200 | 2.5000 | 2.8333 | 0.3333 | 1.0000 |

<small>Values in Å except the hydrogen-bond count. Drift = last-quarter mean minus first-quarter mean. Each protein segment is superposed on itself; segments are never fitted together. Reference: frame 2 of this trajectory, the first analysed frame.</small>

<small>The 'observed only' and 'built only' rows are an analysis added on 2026-11-12 AFTER the pre-registered full-backbone control failed on this project's first production run. Its reason (residues CHARMM-GUI built are a guess, not a measurement) does not depend on the result, but the decision to look came after seeing it: it was NOT pre-registered, and the full-segment row stays the control. 'Observed' means present in the structure uploaded to CHARMM-GUI.</small>

#### Hydrogen bonds, equilibration run

<small>5 ligand–protein pair(s) met the criterion in at least one analysed frame. Criterion: acceptor = N or O by mass; donor = the N or O that the hydrogen is BONDED to in the topology; geometry = donor-acceptor distance <= 3.5 A and donor-H-acceptor angle >= 150.0 deg, counted per frame. No partial charge takes part in deciding which atoms are eligible. The share is a PROPORTION OF ANALYSED FRAMES, never a fraction of time: hydrogen bonds break and re-form faster than a frame is stored.</small>

| Donor | Hydrogen | Acceptor | Frames present / analysed | Share of analysed frames |
|---|---|---|---|---|
| HETA:LIG1:N2 | H12 | PROA:GLN530:O | 22 / 25 | 0.88 |
| HETA:LIG1:N1 | H8 | PROA:ASP594:OD1 | 12 / 25 | 0.48 |
| PROA:GLY596:N | HN | HETA:LIG1:O1 | 11 / 25 | 0.44 |
| PROA:CYS532:N | HN | HETA:LIG1:N3 | 10 / 25 | 0.40 |
| PROA:ASP594:N | HN | HETA:LIG1:N1 | 3 / 25 | 0.12 |

#### Pocket flexibility (RMSF), equilibration run

<small>the CA of 27 of 27 pocket residues. deviation about each atom's own mean position after the same pocket backbone superposition; each value is labelled from the CA atom it was computed from, not from a parallel list. With no apo simulation, a flexible residue cannot be attributed to the ligand.</small>

| Residue | RMSF (Å) |
|---|---|
| PROA:ILE463 | 0.2328 |
| PROA:VAL471 | 0.1708 |
| PROA:TYR472 | 0.1838 |
| PROA:ALA481 | 0.1754 |
| PROA:VAL482 | 0.2085 |
| PROA:LYS483 | 0.1715 |
| PROA:LEU505 | 0.2082 |
| PROA:THR508 | 0.2206 |
| PROA:ILE513 | 0.1804 |
| PROA:LEU514 | 0.1857 |
| PROA:LEU515 | 0.1722 |
| PROA:PHE516 | 0.1800 |
| PROA:ILE527 | 0.1820 |
| PROA:THR529 | 0.1791 |
| PROA:GLN530 | 0.1896 |
| PROA:TRP531 | 0.1960 |
| PROA:CYS532 | 0.2175 |
| PROA:GLY534 | 0.2405 |
| PROA:SER535 | 0.2141 |
| PROA:SER536 | 0.1831 |
| PROA:PHE583 | 0.2050 |
| PROA:ILE592 | 0.1834 |
| PROA:GLY593 | 0.2708 |
| PROA:ASP594 | 0.1985 |
| PROA:PHE595 | 0.2086 |
| PROA:GLY596 | 0.2314 |
| PROA:LEU597 | 0.1928 |

<small>Figure (drawn in the PDF only): Equilibration run, stage step4: drawn from md/3OG7_032_C/analysis_equilibration.json, whose SHA-256 matches the manifest; nothing is recomputed. All 25 analysed frames plotted. Solid lines are full segments, dashed lines the observed residues only (post hoc, not pre-registered). The hydrogen-bond panel counts bonds per stored frame, a sample, not a fraction of time; its count is recorded without times, so it is placed on the ligand series' time axis (same analysed frames, same length), and left out if the lengths differ.</small>

#### Offered reading of the equilibration run (not measured)

<small>READING OFFERED, NOT MEASURED. These lines are an interpretation built from the numbers in this same result and from the thresholds listed here, which are conventions and not measurements of this system. The measurements above stand on their own; this reading does not. Apply your own judgement. Thresholds: protein_plateau_angstrom 0.5, ligand_stable_angstrom 2.0, com_stable_range_angstrom 2.0, hbond_persistent_occupancy 0.5. Source of the thresholds: common conventions, written down with a reason for each on 2026-11-11, before any production trajectory existed but after a restrained equilibration of this project's first simulated system (BRAF with vemurafenib, PDB 3OG7) had been measured, so they were chosen with those numbers in view. No literature is cited because none was consulted. Backbone plateau: drift (last-quarter mean minus first-quarter mean) of at most 0.5 A, of the order of the thermal noise of a protein backbone at 303 K. Ligand RMSD: at most 2.0 A over the last quarter; the same number as this project's docking validation threshold BY CONVENTION, not by derivation: that one compares a calculated pose with a crystal pose, this one asks whether a pose persists over time. Range of the ligand-pocket centre-of-mass distance: at most 2.0 A, to separate a reorientation within the site from a dissociation; the scale comes from the author's earlier simulations of a different system (neuropeptide-GPCR complexes, 10 ns) and is not a value transferable to this one. Hydrogen-bond occupancy: at least 0.5 of the analysed frames, the conventional cut between a contact that holds and an occasional one; it is a sample of frames, not a fraction of time.</small>

- Protein backbone RMSD of segment PROB (the one of 2 with the largest drift) averages 0.3722 A over the last quarter against 0.3028 A over the first (a drift of 0.0694 A, compared with a chosen 0.5 A). That is consistent with a plateau. THIS LINE CARRIES NO EVIDENCE: the NAMD input of this stage holds these atoms in place (2174 of 2174 protein backbone atoms are held in place by the restraint file, so their motion is damped and NOT free. Any figure built on those atoms passes by construction), so the figure cannot move and the comparison passes by construction. Read it again on an unrestrained production trajectory.
- Split into observed and BUILT residues (built = absent from, or without a full backbone in, CHARMM-GUI's .orig.pdb; observed = present in the structure uploaded to CHARMM-GUI, which is experimental only if that was a deposited structure - not checked here): PROA: the observed residues drift 0.066 A over 987 atoms (within the chosen 0.5 A; full segment 0.0683 A), with 25 of 272 residues built and left out (545-547, 597-614, 627-630); the built residues alone drift 0.0901 A and average 0.4432 A over the last quarter; PROB: the observed residues drift 0.0675 A over 1015 atoms (within the chosen 0.5 A; full segment 0.0694 A), with 18 of 272 residues built and left out (601-614, 627-630); the built residues alone drift 0.0928 A and average 0.489 A over the last quarter. For PROB, the segment of the line above, the full backbone is within the chosen figure already, so the split changes nothing there. THIS SPLIT WAS NOT PRE-REGISTERED: it was added on 2026-11-12, after the pre-registered full-backbone control failed on this project's first production run (the 10 ns production of BRAF with vemurafenib, PDB 3OG7). Its reason - built residues are a guess, not a measurement - does not depend on the result, but the decision to look came after seeing it. The line above remains the pre-registered control. THIS LINE CARRIES NO EVIDENCE: the NAMD input of this stage holds these atoms in place (2174 of 2174 protein backbone atoms are held in place by the restraint file, so their motion is damped and NOT free. Any figure built on those atoms passes by construction), so the figure cannot move and the comparison passes by construction. Read it again on an unrestrained production trajectory.
- 1 pocket residue(s) were built, not observed (PROA:LEU597). The ligand RMSD is superposed on the pocket backbone, theirs included, and the pocket itself was selected around a ligand placed next to a guessed geometry.
- Ligand RMSD reaches 0.4184 A at most and averages 0.3505 A over the last quarter (chosen comparison: 2.0 A). The pose held within that figure. THIS LINE CARRIES NO EVIDENCE: the NAMD input of this stage holds these atoms in place (33 of 33 ligand heavy atoms are held in place by the restraint file, so their motion is damped and NOT free. Any figure built on those atoms passes by construction), so the figure cannot move and the comparison passes by construction. Read it again on an unrestrained production trajectory.
- Centre-of-mass distance stays in a 0.312 A range (chosen comparison: 2.0 A), mean 2.7948 A. No dissociation is visible in this window. THIS LINE CARRIES NO EVIDENCE: the NAMD input of this stage holds these atoms in place (33 of 33 ligand heavy atoms are held in place by the restraint file, so their motion is damped and NOT free. Any figure built on those atoms passes by construction), so the figure cannot move and the comparison passes by construction. Read it again on an unrestrained production trajectory.
- 1 of 5 hydrogen-bond pairs held in at least 50% of the analysed frames: HETA:LIG1:N2 -> PROA:GLN530:O (88%). This figure is not marked as carrying no evidence, because an occupancy is not a displacement - but the atoms on both sides of these bonds are held in place (33 of 33 ligand heavy atoms are held in place by the restraint file, so their motion is damped and NOT free. Any figure built on those atoms passes by construction), so the geometry was kept near its starting value and the occupancy is an upper bound rather than a result.

<small>Not measured: Binding affinity. Nothing here computes a free energy, and none of these metrics is one: a stable pose is not a potent ligand. Bridging waters. Left out of this version on purpose. Whether the pose is correct. With no apo control, a change in the protein cannot be attributed to the ligand.</small>

## 3. ADMET

<small>3 compounds recorded. Section 3 holds three different kinds of statement, in three separate tables, and they are not comparable with each other or with the calculated docking scores of section 1: the ChEMBL record (facts of a database entry), the retrieved measurements (each with the units and the assay it was measured in) and the calculated properties (computed from the ligand's atoms and bonds). NOTHING in this section is a prediction: no value here comes from a trained model.</small>

### 8C7Y · pocket TXV U · TXV

<small>Recorded 2026-11-12T18:59:15Z by analyze_admet. TXV of 8C7Y pocket TXV U, read from the docked pose. InChIKey GCVULUAJPYRJHT-UHFFFAOYSA-N matches exactly 1 ChEMBL molecule (CHEMBL5440746). 0 of 0 rows read are on this project's ADME allowlist. All 0 rows ChEMBL holds in the assay types queried were read. RDKit: MW 478.886, cLogP 5.779, TPSA 98.64.</small>

#### Compound identity

| Field | Value |
|---|---|
| InChIKey (computed here with RDKit) | GCVULUAJPYRJHT-UHFFFAOYSA-N |
| Read from | docked pose |
| File | dock/8C7Y_TXV_U/poses_TXV_regenerated.sdf |
| SHA-256 of that file | a4f0505f979d63fcb9e98c8c78d3c941ac4fa490443690b7ddcb29d8c3c69217 |
| Docking run | TXV_regenerated |
| Agreement between the files | the docked pose and the crystallographic SDF give the same InChIKey |
| Heavy atoms (SDF / manifest pocket) | 34 / 34 |
| ChEMBL molecules with this exact InChIKey | 1 |

#### ChEMBL record (facts of the entry, not measurements)

<small>[ChEMBL CHEMBL5440746](https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL5440746/)<br> · name in ChEMBL: unnamed (no pref_name)</small>

| Field | Value |
|---|---|
| max_phase | no record |
| first_approval | no record |
| withdrawn_flag | False (default, see below) |
| black_box_warning | 0 |
| oral | False (default, see below) |
| parenteral | False (default, see below) |
| topical | False (default, see below) |
| molecule_type | no record |
| therapeutic_flag | False (default, see below) |
| usan_stem | no record |
| first_in_class | -1 (-1, see below) |
| ATC | no record |

<small>Two kinds of value in this block are NOT findings, and ChEMBL's API does not distinguish them from findings. `fields_minus_one` lists the fields that came back as -1: this project has not read ChEMBL's schema, so it does not say what -1 means and never reads it as a 'no'. `boolean_fields_false` lists the flags that came back False: on an entry with nothing annotated they are all False, so a False here must not be reported as 'this compound is not oral' or 'was not withdrawn'. Both lists keep the raw values; neither is reinterpreted.</small>

<small>Returned as -1: availability_type, chirality, first_in_class, inorganic_flag, orphan, prodrug, veterinary.</small>

<small>Returned as False: dosed_ingredient, oral, parenteral, therapeutic_flag, topical, withdrawn_flag. On an entry with nothing annotated every one of these is False, so none of them is evidence that the compound lacks that property.</small>

#### Measured ADME (retrieved from ChEMBL)

<small>0 of 0 rows read from ChEMBL assay types A, P are on this project's ADME allowlist. Assay types B, F, T, U were not queried, so a measurement filed only there is not here. Values are quoted with their own units and are never converted or averaged; a relation ('>' or '<') is part of the value.</small>

<small>No row in the assay types queried is on this project's ADME allowlist.</small>

#### Calculated properties (not measured, not predicted)

<small>Two engines calculate some of the same quantities: RDKit here, from the exact SDF that was docked, and ChEMBL's own engine on its entry. They are reported side by side, never averaged and never presented as confirming each other: two calculations agreeing is not a measurement, and if they disagree neither can be called correct from here.</small>

| Property | This project (RDKit 2026.09.1) | ChEMBL's own engine |
|---|---|---|
| Molecular weight | 478.886 | 478.89 |
| cLogP | 5.779 | 5.78 |
| TPSA | 98.64 | 98.64 |
| H-bond donors | 2 | 2 |
| H-bond acceptors | 4 | 4 |
| Rotatable bonds | 5 | 5 |
| Aromatic rings | 4 | 4 |
| Heavy atoms | 34 | 34 |
| Fraction Csp3 | 0.12 | no record |
| QED | 0.355 | 0.36 |
| Formula | C25H17ClF2N4O2 | C25H17ClF2N4O2 |

#### Property rules (offered, not measured)

<small>READING OFFERED, NOT MEASURED. These rules are conventions of common use, pre-registered in this project on 2026-11-11 before any code was written. They are NOT measurements, NOT predictions and NOT derived from any source retrieved by this project; the names are labels, not citations. They are property filters on one molecule: they do not predict bioavailability, absorption or toxicity, and an approved drug can break them. The researcher's own judgement decides.</small>

| Rule | Cut-offs | Violations | Exceeded |
|---|---|---|---|
| Lipinski | molecular_weight <= 500; clogp <= 5; hbd <= 5; hba <= 10 | 1 | clogp 5.779 > 5 |
| Veber | rotatable_bonds <= 10; tpsa <= 140 | 0 | none |
| Egan | tpsa <= 131.6; clogp <= 5.88 | 0 | none |

<small>No value in this output is PREDICTED: nothing here comes from a trained model (Phase 3 step 7.1).</small>

<small>This output holds NO name for this compound. Refer to this compound by its PDB chemical component id and its InChIKey. ChEMBL gives no `pref_name` for it here, so THIS OUTPUT holds no name: do not invent one, do not give it a code of your own, and do not tie it to a named drug -- 'an analogue of X', 'a linked X', 'X dimer' -- because nothing retrieved here says that. One other name does exist outside this output: wwPDB's systematic chemical name of the component, retrieved by get_structures and kept in the manifest (the research report prints it). It is a systematic name, not a compound or trade name, and this output deliberately does not carry it; if a structural claim is made from it, say it comes from that systematic name. Size, formula and descriptors are in `calculated_properties` and can be reported as what they are.</small>

<small>hERG and other off-target safety data are NOT in this output: in ChEMBL they are assay_type=F against a target, not ADME (decision of 2026-11-11).</small>

### 5JSM · pocket 6NB J · 6NB

<small>Recorded 2026-11-12T19:00:55Z by analyze_admet. 6NB of 5JSM pocket 6NB J, read from the crystallographic SDF (no docking run for this pocket). InChIKey VESCTIQVMMXAQY-UHFFFAOYSA-N has no ChEMBL entry, so there are no records and no measurements; only the RDKit calculation and the rule reading apply. RDKit: MW 1057.114, cLogP 9.586, TPSA 220.76.</small>

#### Compound identity

| Field | Value |
|---|---|
| InChIKey (computed here with RDKit) | VESCTIQVMMXAQY-UHFFFAOYSA-N |
| Read from | crystallographic SDF (no docking run for this pocket) |
| File | ligands/5JSM_6NB_J.sdf |
| SHA-256 of that file | 576fe79b16d56f1e60cde44d15bc5bac778bc8cc84fca9e6acd03a38893889c0 |
| Docking run | none (crystallographic ligand) |
| Agreement between the files | only one SDF was available, so there was nothing to compare |
| Heavy atoms (SDF / manifest pocket) | 74 / 74 |
| ChEMBL molecules with this exact InChIKey | 0 |

#### ChEMBL record (facts of the entry, not measurements)

<small>No record. No ChEMBL molecule has this exact InChIKey (no record). NO NAME IN CHEMBL: no ChEMBL molecule has this exact InChIKey, so ChEMBL holds no name for it. Refer to this compound by its PDB chemical component id and its InChIKey. ChEMBL gives no `pref_name` for it here, so THIS OUTPUT holds no name: do not invent one, do not give it a code of your own, and do not tie it to a named drug -- 'an analogue of X', 'a linked X', 'X dimer' -- because nothing retrieved here says that. One other name does exist outside this output: wwPDB's systematic chemical name of the component, retrieved by get_structures and kept in the manifest (the research report prints it). It is a systematic name, not a compound or trade name, and this output deliberately does not carry it; if a structural claim is made from it, say it comes from that systematic name. Size, formula and descriptors are in `calculated_properties` and can be reported as what they are. Only the calculated properties and the rule reading below apply.</small>

#### Measured ADME (retrieved from ChEMBL)

<small>No ADME measurement is reported because no ChEMBL molecule has this exact InChIKey, so no ADME query was made. ChEMBL holds nothing for a compound it does not have; a measurement published elsewhere would not be here.</small>

#### Calculated properties (not measured, not predicted)

<small>Two engines calculate some of the same quantities: RDKit here, from the exact SDF that was docked, and ChEMBL's own engine on its entry. They are reported side by side, never averaged and never presented as confirming each other: two calculations agreeing is not a measurement, and if they disagree neither can be called correct from here.</small>

| Property | This project (RDKit 2026.09.1) | ChEMBL's own engine |
|---|---|---|
| Molecular weight | 1057.114 | no record |
| cLogP | 9.586 | no record |
| TPSA | 220.76 | no record |
| H-bond donors | 4 | no record |
| H-bond acceptors | 12 | no record |
| Rotatable bonds | 25 | no record |
| Aromatic rings | 8 | no record |
| Heavy atoms | 74 | no record |
| Fraction Csp3 | 0.231 | no record |
| QED | 0.024 | no record |
| Formula | C52H48F4N6O10S2 | no record |

#### Property rules (offered, not measured)

<small>READING OFFERED, NOT MEASURED. These rules are conventions of common use, pre-registered in this project on 2026-11-11 before any code was written. They are NOT measurements, NOT predictions and NOT derived from any source retrieved by this project; the names are labels, not citations. They are property filters on one molecule: they do not predict bioavailability, absorption or toxicity, and an approved drug can break them. The researcher's own judgement decides.</small>

| Rule | Cut-offs | Violations | Exceeded |
|---|---|---|---|
| Lipinski | molecular_weight <= 500; clogp <= 5; hbd <= 5; hba <= 10 | 3 | molecular_weight 1057.114 > 500; clogp 9.586 > 5; hba 12 > 10 |
| Veber | rotatable_bonds <= 10; tpsa <= 140 | 2 | rotatable_bonds 25 > 10; tpsa 220.76 > 140 |
| Egan | tpsa <= 131.6; clogp <= 5.88 | 2 | tpsa 220.76 > 131.6; clogp 9.586 > 5.88 |

<small>No value in this output is PREDICTED: nothing here comes from a trained model (Phase 3 step 7.1).</small>

<small>This output holds NO name for this compound. Refer to this compound by its PDB chemical component id and its InChIKey. ChEMBL gives no `pref_name` for it here, so THIS OUTPUT holds no name: do not invent one, do not give it a code of your own, and do not tie it to a named drug -- 'an analogue of X', 'a linked X', 'X dimer' -- because nothing retrieved here says that. One other name does exist outside this output: wwPDB's systematic chemical name of the component, retrieved by get_structures and kept in the manifest (the research report prints it). It is a systematic name, not a compound or trade name, and this output deliberately does not carry it; if a structural claim is made from it, say it comes from that systematic name. Size, formula and descriptors are in `calculated_properties` and can be reported as what they are.</small>

<small>hERG and other off-target safety data are NOT in this output: in ChEMBL they are assay_type=F against a target, not ADME (decision of 2026-11-11).</small>

### 3OG7 · pocket 032 C · 032

<small>Recorded 2026-11-12T17:41:43Z by analyze_admet. 032 of 3OG7 pocket 032 C, read from the docked pose. InChIKey GPXBXXGIAQBQNI-UHFFFAOYSA-N matches exactly 1 ChEMBL molecule (CHEMBL1229517). 4 of 18 rows read are on this project's ADME allowlist (clearance, permeability, solubility); the other 14 are not on the allowlist, which is narrower than ADME, so they are not 'not ADME'. All 18 rows ChEMBL holds in the assay types queried were read. RDKit: MW 489.931, cLogP 5.544, TPSA 91.92.</small>

#### Compound identity

| Field | Value |
|---|---|
| InChIKey (computed here with RDKit) | GPXBXXGIAQBQNI-UHFFFAOYSA-N |
| Read from | docked pose |
| File | dock/3OG7_032_C/poses_032_regenerated.sdf |
| SHA-256 of that file | b708374495ef4eaa5257cfa931007749f7964f3e8fbf590f2cab86d00a5eddf4 |
| Docking run | 032_regenerated |
| Agreement between the files | the docked pose and the crystallographic SDF give the same InChIKey |
| Heavy atoms (SDF / manifest pocket) | 33 / 33 |
| ChEMBL molecules with this exact InChIKey | 1 |

#### ChEMBL record (facts of the entry, not measurements)

<small>[ChEMBL CHEMBL1229517](https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL1229517/)<br> · name in ChEMBL: VEMURAFENIB</small>

| Field | Value |
|---|---|
| max_phase | 4.0 |
| first_approval | 2011 |
| withdrawn_flag | False (default, see below) |
| black_box_warning | 0 |
| oral | True |
| parenteral | False (default, see below) |
| topical | False (default, see below) |
| molecule_type | Small molecule |
| therapeutic_flag | True |
| usan_stem | -rafenib |
| first_in_class | 1 |
| ATC | L01EC01 |

<small>Two kinds of value in this block are NOT findings, and ChEMBL's API does not distinguish them from findings. `fields_minus_one` lists the fields that came back as -1: this project has not read ChEMBL's schema, so it does not say what -1 means and never reads it as a 'no'. `boolean_fields_false` lists the flags that came back False: on an entry with nothing annotated they are all False, so a False here must not be reported as 'this compound is not oral' or 'was not withdrawn'. Both lists keep the raw values; neither is reinterpreted.</small>

<small>Returned as False: parenteral, topical, withdrawn_flag. On an entry with nothing annotated every one of these is False, so none of them is evidence that the compound lacks that property.</small>

#### Measured ADME (retrieved from ChEMBL)

<small>4 of 18 rows read from ChEMBL assay types A, P are on this project's ADME allowlist. Assay types B, F, T, U were not queried, so a measurement filed only there is not here. Values are quoted with their own units and are never converted or averaged; a relation ('>' or '<') is part of the value.</small>

| Property | Type | Value | Units | Assay |
|---|---|---|---|---|
| clearance | CL | 6.61 | uL.min-1.(10^6cells)-1 | ASTRAZENECA: Intrinsic clearance measured in human hepatocytes following incubation at 37C. Experimental range <3 to >150 microL/min/1E6 cells. Rapid Commun. Mass Spectrom. 2010, 24, 1730-1736. |
| permeability | permeability | = 2.9 | 10^-6 cm/s | Permeability in human Caco-2 cells |
| permeability | Papp | = 3.0 | 10^-6 cm/s | Apparent permeability across apical to basolateral side in LLC-PK1 cells |
| solubility | Solubility | = 2.0 | ug.mL-1 | Aqueous solubility of the compound in FaSSIF at pH 6.5 measured by dissolution tester |

<small>Rows in those assay types that are NOT on this project's allowlist and were therefore excluded: Activity (9), IC50 (3), Ratio (2). This does not mean they are not ADME. The allowlist is narrower than ADME on purpose: the generic types 'Activity' and 'Ratio' carry real ADME in some rows (efflux ratios, CYP induction) and cytotoxicity or target-engagement in others, and the two cannot be separated without reading each assay description. An excluded row may hold a measurement this subsection does not show.</small>

<small>The excluded rows themselves (14 of 14), unclassified and NOT part of the count above, so each can be read on its own terms:</small>

| Type | Value | Units | Assay (ChEMBL's own description) |
|---|---|---|---|
| Activity | > 75.0 | % | Cytotoxicity against human HepG2 cells assessed as cell viability at 1 uM after 4 hrs by Cell-Titer Glo assay |
| Activity | > 75.0 | % | Cytotoxicity against human HepG2 cells assessed as cell viability at 10 uM after 4 hrs by Cell-Titer Glo assay |
| Activity | > 75.0 | % | Cytotoxicity against human HepG2 cells assessed as cell viability at 100 uM after 4 hrs by Cell-Titer Glo assay |
| Activity | > 75.0 | % | Cytotoxicity against human HepG2 cells assessed as cell viability at 1 uM after 24 hrs by Cell-Titer Glo assay |
| Activity | > 75.0 | % | Cytotoxicity against human HepG2 cells assessed as cell viability at 10 uM after 24 hrs by Cell-Titer Glo assay |
| Activity | = 72.0 | % | Cytotoxicity against human HepG2 cells assessed as cell viability at 100 uM after 24 hrs by Cell-Titer Glo assay |
| Activity | no record | — | Effect on B-Raf level in human A549 cells at five times IC50 value after 12 to 24 hrs by Western blot analysis |
| Activity | no record | — | Activation of CRAF in human Calu6 cells assessed as increase in MEK phosphorylation at 3 uM after 2 hrs by FRET assay |
| IC50 | no record | — | Cytotoxicity against human HELF cells assessed as reduction in cell viability incubated for 72 hrs by CellTiter-Glo assay |
| IC50 | no record | — | Cytotoxicity against mouse C2C12 cells assessed as reduction in cell viability incubated for 72 hrs by CellTiter-Glo assay |
| IC50 | > 100000.0 | nM | Cytotoxicity against human T75 cells |
| Activity | = 39.0 | % | Induction of CYP450 in human assessed as reduction in midazolam AUC |
| Ratio | = 21.7 | — | Substrate activity at P-gp in LLC-PK1 cells assessed as efflux ratio of apparent permeability from apical to basolateral side over basolateral to apical side |
| Ratio | = 15.2 | — | Substrate activity at BCRP in LLC-PK1 cells assessed as efflux ratio of apparent permeability from apical to basolateral side over basolateral to apical side |

#### Calculated properties (not measured, not predicted)

<small>Two engines calculate some of the same quantities: RDKit here, from the exact SDF that was docked, and ChEMBL's own engine on its entry. They are reported side by side, never averaged and never presented as confirming each other: two calculations agreeing is not a measurement, and if they disagree neither can be called correct from here.</small>

| Property | This project (RDKit 2026.09.1) | ChEMBL's own engine |
|---|---|---|
| Molecular weight | 489.931 | 489.93 |
| cLogP | 5.544 | 5.54 |
| TPSA | 91.92 | 91.92 |
| H-bond donors | 2 | 2 |
| H-bond acceptors | 4 | 4 |
| Rotatable bonds | 7 | 7 |
| Aromatic rings | 4 | 4 |
| Heavy atoms | 33 | 33 |
| Fraction Csp3 | 0.13 | no record |
| QED | 0.333 | 0.33 |
| Formula | C23H18ClF2N3O3S | C23H18ClF2N3O3S |

#### Property rules (offered, not measured)

<small>READING OFFERED, NOT MEASURED. These rules are conventions of common use, pre-registered in this project on 2026-11-11 before any code was written. They are NOT measurements, NOT predictions and NOT derived from any source retrieved by this project; the names are labels, not citations. They are property filters on one molecule: they do not predict bioavailability, absorption or toxicity, and an approved drug can break them. The researcher's own judgement decides.</small>

| Rule | Cut-offs | Violations | Exceeded |
|---|---|---|---|
| Lipinski | molecular_weight <= 500; clogp <= 5; hbd <= 5; hba <= 10 | 1 | clogp 5.544 > 5 |
| Veber | rotatable_bonds <= 10; tpsa <= 140 | 0 | none |
| Egan | tpsa <= 131.6; clogp <= 5.88 | 0 | none |

<small>No value in this output is PREDICTED: nothing here comes from a trained model (Phase 3 step 7.1).</small>

<small>hERG and other off-target safety data are NOT in this output: in ChEMBL they are assay_type=F against a target, not ADME (decision of 2026-11-11).</small>

## 4. Methods and limitations

This record is generated deterministically from prep/manifest.json: every value and every identifier is quoted from the manifest, fields without data are shown as 'no record', and nothing is recomputed or inferred here. No text of a Vina log is copied into this document.

- The docking results in section 1 (Vina scores and RMSD values) are CALCULATED by AutoDock Vina and RDKit, as recorded in the manifest. The provenance lines quote other kinds of value: the resolution is the experimental value deposited for the PDB entry, the box was computed when the structure was prepared, and run times are measurements of this machine. A docking score is not a measured affinity and is not comparable with the IC50/Ki/Kd values that ChEMBL reports for the same target, which stay in the research report. The poses of one run are listed in the order Vina returned them; no compound is ranked against another here.
- Re-docking the crystal ligand of a pocket into its own receptor tests the set-up, not the method: the scope line of each run states what it covers.
- The receptor is rigid and was not completed: missing loops stay missing, side chains are not rebuilt, and residues that match no Meeko template are listed per run (deleted only when outside the box).
- Wall-clock and CPU times are measurements of the machine that ran the docking, not properties of the result.

Rules that produced the receptor, the pocket and the box used above (prep/manifest.json, decisions):

- Original files are kept unmodified in raw/. Loop completion, protonation and charges are Phase 3.
- Receptor (clean/<PDB>_receptor.pdb): target chains only (RCSB mapping), first model, first alternate conformation, no hydrogens, no waters, no additives, no ligands; nucleotides, cofactors and metal ions MG/MN/ZN/FE/CO/NI/CU that touch a target chain are KEPT (e.g. KRAS GDP + Mg, 8TVK); NA, K, CL, CA, CD are removed. Coordinates are unchanged. PDB format unless chain names or atom count do not fit, then mmCIF.
- Pocket: target-chain residues with any atom within 5 Å of the ligand instance (gemmi, first alternate conformation of the ligand), in the asymmetric unit only: crystal-symmetry copies are excluded (they put residues ~47 Å away in 8C7Y). UniProt numbers come from label_seq_id + RCSB alignment; residues of other polymer chains are listed apart. No pocket for nucleotides and cofactors (kept in the receptor), for ligands with fewer than 10 heavy atoms, or without a target-chain neighbour in RCSB: each is listed in pocket_not_computed. Box: ligand bounding box + 5 Å per side, same frame as the mmCIF; sides above 30 Å are flagged. A ligand instance with pocket residues in more than one target chain is flagged spans_target_chains (e.g. a bivalent compound bridging two sites): its pocket and box cover every chain it touches and are NOT split into one box per site; residues_by_chain gives the count per chain.
- Chain metrics: missing_uniprot_ranges are residues of the sample sequence not modelled inside the modelled span; a missing range next to a pocket residue is flagged. variant_check reads the residue actually modelled at the variant position in each target chain and compares it with the UniProt canonical letter (not the letter typed in the filter): the variant can be named in the title and still be unmodelled in the chain that holds the ligand (3OG7 A). Two residues at one position (microheterogeneity) give status 'ambiguous'.
- AlphaFold pLDDT (B-factor column, equal to _ma_qa_metric_local) is read at the same UniProt numbers as each experimental pocket, only for the canonical-sequence model; residues without a UniProt number are not compared (alphafold_residues_compared vs n_residues).
- Target chains come from RCSB's chain-to-UniProt mapping, not from the mmCIF _struct_ref, which is unreliable (3OG7 cites UniProt Q5IBP5; antibody chains often cite no UniProt).
- Chains are given as label_asym_id (mmCIF) and auth_asym_id (author, used by PyMOL and the PDB format); they can differ (7ZR0: BRAF is label D, auth K).
- Ligand poses are the crystal coordinates of each instance (RCSB ModelServer, by label_asym_id). A ligand can sit on one chain only (3OG7: chain A) or touch several.
- neighbor_*_chains come from RCSB rcsb_target_neighbors with the closest distance per chain; RCSB does not document the distance cutoff.
- construct.mutations are the depositor's declared mutations of the target entity; they can include solubility/engineering mutations besides the variant (3OG7 declares 13). Any declared mutation means the structure is not the native sequence; an empty list does NOT prove it is native (the variant may be named only in the title, or the text may be unparsed, e.g. 'YES').
- Additives (ions, buffers, cryoprotectants) get no SDF. The classification is the one of get_structures; some additives (e.g. BEN, EOH in 5JSM) are not in its list and count as ligands.
- Glycans (RCSB branched entities, e.g. NAG chains) are not queried and are absent from nonpolymer_instances.

Preparation steps still pending in the manifest: none listed. Molecular dynamics (section 2): the simulations ran outside this machine; their trajectories were measured locally by analyze_md and every figure, file hash and plotted series is quoted from what it recorded.
