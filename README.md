The idea: start from a large compound library (thousands–millions of
molecules), cheaply filter it down to a small, drug-like,
structurally-clean shortlist (Stages 1–3), then run the
compute-intensive docking + MD analysis (Stages 4–8) only on that
shortlist.

---

## Workflow Stages

The full process follows this sequence:

```
Stage 1: Standardize & deduplicate SMILES
   |
Stage 2: Lipinski/Veber + PAINS/BRENK structural filters
   |
Stage 3: ADMET rule-based filtering
   |
Stage 4: Ligand preparation & format conversion (PDB -> MOL2 -> PDBQT)
   |
Stage 5: Protein-ligand docking (AutoDock Vina)
   |
Stage 6: Complex generation
   |
Stage 7: OPLS-AA force-field parameter generation (LigParGen, for GROMACS)
   |
Stage 8: GROMACS MD simulation setup
   |
Stage 9: PyMOL visualization & animation generation
```

---

### Stage 1: Standardize & Deduplicate

**Script:** `stage1_standardize.py`

- Reads the raw compound library (`data/input/library.csv`).
- Canonicalizes every SMILES string with RDKit and discards molecules
  that fail to parse/sanitize.
- Removes duplicate molecules (by canonical SMILES), both within and
  across processing chunks.

**Output:** `data/intermediate/stage1/all_standardized.parquet`

```bash
python stage1_standardize.py --config config.yaml
```

---

### Stage 2: Lipinski/Veber + PAINS/BRENK Filters

**Script:** `stage2_filters.py`

- Computes 2D descriptors (MW, LogP, HBD, HBA, TPSA, rotatable bonds,
  heavy atoms, aromatic rings, QED, etc.) for every Stage 1 molecule.
- Rejects molecules violating **Lipinski's Rule of Five** /
  **Veber's rules** (drug-likeness / oral bioavailability proxies).
- Rejects molecules matching **PAINS** (assay interference) or
  **BRENK** (reactive/unstable/toxic) structural alert patterns.

**Output:** `data/intermediate/stage2/filtered.parquet`

```bash
python stage2_filters.py --config config.yaml
```

---

### Stage 3: ADMET Rule-Based Filtering

**Script:** `stage3_admet.py`

- Applies a second, tighter set of ADMET (Absorption, Distribution,
  Metabolism, Excretion, Toxicity) property thresholds to the Stage 2
  survivors, re-using their already-computed descriptors.
- Requires **all** conditions (LogP/TPSA range, MW, HBD, HBA,
  rotatable bonds, QED minimum, aromatic ring count) to pass
  simultaneously.

**Output:** `data/intermediate/stage3/admet_passed.parquet`

```bash
python stage3_admet.py --config config.yaml
```

> Optional: use `bootstrap_sample_dock.py` + `label_qsar_set.py` +
> `train_qsar.py` + `stage4_qsar.py` to further rank/shortlist the
> Stage 3 output with a QSAR model before moving on to Stage 4 below
> (see `PROTOCOL.md`).

---

### Stage 4: Ligand Preparation & Format Conversion

**Script:** `ligand_format_pdbqt.sh`
*(from Protein-Ligand-Binding-Affinity-Workflow)*

- Takes the shortlisted ligands from Stage 3 (or the QSAR/docking
  shortlist) and splits/converts them through the format chain
  **PDB → MOL2 → PDBQT**, preparing each ligand for AutoDock Vina.

---

### Stage 5: Protein–Ligand Docking (AutoDock Vina)

**Script:** `automatic_docking_ligand.sh`
*(from Protein-Ligand-Binding-Affinity-Workflow)*

- Prepares the receptor (protein target) and docking grid box.
- Batch-docks every prepared ligand (PDBQT) from Stage 4 against the
  receptor using AutoDock Vina.
- Produces a binding-affinity score (kcal/mol) per ligand pose.

---

### Stage 6: Complex Generation

- Combines the best-scoring docked ligand pose(s) with the receptor
  structure to form protein–ligand complex files, used as the starting
  structures for force-field parameterization and MD simulation.

---

### Stage 7: OPLS-AA Force-Field Parameter Generation

**Script:** `opls_ff_ligand.sh` (LigParGen)
*(from Protein-Ligand-Binding-Affinity-Workflow)*

- Submits each docked ligand to the **LigParGen** server (or local
  LigParGen tooling) to generate OPLS-AA force-field parameters
  required by GROMACS.
- Supports single-ligand parameterization, configurable output
  formats, charge models for charged molecules, and preview-only runs.

---

### Stage 8: GROMACS MD Simulation Setup

- Combines the protein–ligand complex (Stage 6) with the OPLS-AA
  ligand parameters (Stage 7) to assemble a solvated, ion-balanced
  simulation system.
- Runs energy minimization, equilibration, and production molecular
  dynamics in GROMACS.

---

### Stage 9: PyMOL Visualization & Animation

- PyMOL renders trajectory frames from the MD run (recommended
  800×600 resolution).
- ImageMagick's `convert` stitches the frames into a GIF animation
  (configurable frame delay).

---

## Installation Requirements

### Cheminformatics environment (Stages 1–3)
```bash
conda create -n molscreen python=3.11 -y
conda activate molscreen
conda install -c conda-forge rdkit pandas numpy pyarrow pyyaml -y
```

### Docking / MD environment (Stages 4–9)
```bash
conda create -n vina python=3.10 -y
conda activate vina
conda config --add channels conda-forge
conda install -c conda-forge vina numpy boost-cpp sphinx -y
conda install -c conda-forge openbabel meeko -y
```
Also required: **PyMOL**, **GROMACS**, and **ImageMagick**.

---

## Execution Steps Summary

```bash
# Stages 1-3: cheminformatics triage (this repo)
python stage1_standardize.py --config config.yaml
python stage2_filters.py --config config.yaml
python stage3_admet.py --config config.yaml

# (optional) Stages 4 (QSAR bootstrap) -- see PROTOCOL.md

# Stage 4: ligand preparation & format conversion
bash ligand_format_pdbqt.sh

# Stage 5: docking
bash automatic_docking_ligand.sh

# Stage 6: complex generation (output of Stage 5 + receptor)

# Stage 7: OPLS-AA parameters
bash opls_ff_ligand.sh

# Stage 8: GROMACS MD setup & run (gmx commands using Stage 6 + Stage 7 outputs)

# Stage 9: visualization
pymol -c render_trajectory.pml
convert -delay 10 -loop 0 frame_*.png trajectory.gif
```

---

## Prerequisites

- Python 3.10/3.11, RDKit, AutoDock Vina, Open Babel, Meeko, PyMOL,
  GROMACS, ImageMagick.
- A prepared receptor structure for your target protein.

## Important Notes

- **Verify docking box coordinates** (`docking.center_x/y/z`,
  `docking.size_x/y/z` in `config.yaml`) against your target's actual
  binding pocket before any production docking runs (Stage 5).
- **Validate OPLS-AA parameters** (Stage 7) before launching MD
  simulations (Stage 8).
- Stages 1–3 are designed for **large libraries** (thousands–millions
  of compounds); Stages 4–9 are designed for a **small shortlist**
  (tens of compounds) due to their computational cost.

# Directory Structure

```bash
.
├── ligand_format_pdbqt.sh
├── automatic_docking_ligand.sh
├── opls_ff_ligand.sh
├── ligpargen.py
│
├── ligand_prep/
│   ├── ligand_1.pdbqt
│   ├── ligand_2.pdbqt
│   └── ...
│
├── autodoc1/
├── autodoc2/
├──.........
│   ├── protein_ligand1_ad4_out.pdbqt
│   ├── protein_ligand2_ad4_out.pdbqt
│   └── ...
├── animation.gif
└── README.md
```

---

# Software Requirements

- Python 3.10
- AutoDock Vina
- Open Babel
- Meeko
- PyMOL
- GROMACS
- VMD (optional)
- ImageMagick

---

# Notes

- Verify docking box coordinates before production docking.
- Generated OPLS parameters should be validated before MD simulations.

---
