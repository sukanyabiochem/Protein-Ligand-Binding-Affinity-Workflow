# Automated Protein–Ligand Docking and MD Simulation Pipeline

This repository contains an automated workflow for:

- Protein–ligand docking using AutoDock Vina
- Ligand preparation and format conversion
- OPLS force-field generation for GROMACS
- Molecular dynamics (MD) simulation setup
- PyMOL movie generation
- GIF animation creation

The workflow supports batch processing of multiple ligands for high-throughput computational studies.

---

# Project Workflow

```text
Ligand Preparation
        ↓
PDB → MOL2 → PDBQT Conversion
        ↓
AutoDock Vina Docking
        ↓
Docked Complex Generation
        ↓
LigParGen OPLS Parameter Generation
        ↓
GROMACS MD Simulation Setup
        ↓
PyMOL Visualization & Movie Generation
```

---

# Installation

## Create a Conda Environment for AutoDock Vina

```bash
conda create -n vina python=3.10
conda activate vina
```

## Configure Conda Channels

```bash
conda config --env --add channels conda-forge
```

## Install AutoDock Vina

```bash
conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
conda install -c conda-forge vina
```

## Check Installation

```bash
vina --help
```

---

# Install Additional Packages

## Open Babel

Used for ligand format conversion.

```bash
conda install -c conda-forge openbabel
```

## Meeko

Used for generating PDBQT files.

```bash
conda install -c conda-forge meeko
```

---

## Prepare a library of Ligands ##
The PDB file Processed in Openbabel:
1. Download the PDB files from the [ZINC Database](https://zinc.docking.org/).
2. Other database source used in this work [PubChem](https://pubchem.ncbi.nlm.nih.gov/) & [RCSB PDB](https://rcsb.org/)

---
# AutoDock Workflow

---

## Step 1 — Ligand Preparation

Split multiple ligands into separate files and convert formats.

Run:

```bash
bash ligand_format_pdbqt.sh
```

This script performs:

- Ligand splitting
- PDB → MOL2 conversion
- MOL2 → PDBQT conversion

---

## Step 2 — Run AutoDock Vina Docking

Run automated docking for all ligands:

```bash
bash automatic_docking_ligand.sh
```

This script performs:

- Receptor preparation
- Grid generation
- Batch docking using AutoDock Vina
- Docked complex generation

---

# GROMACS MD Simulation Workflow

---

## Generate OPLS Force-Field Parameters

Prepare ligand topology and parameter files for GROMACS.

The `ligpargen.py` script requires several Python libraries for:

- Web requests
- HTML parsing
- Automated browser interaction
- Download management

# Python Package Dependencies

| Package | Purpose |
|---|---|
| `requests` | HTTP requests and file downloads |
| `beautifulsoup4` | HTML parsing and response extraction |
| `certifi` | SSL certificate validation |

Install dependencies:

```bash
conda install -c conda-forge mechanize requests beautifulsoup4 html5lib
```


Run:

```bash
bash opls_ff_ligand.sh
```

The script obtains OPLS-AA force-field parameters from the LigParGen server.

---

# LigParGen Python Workflow

Examples for generating ligand parameters.

---

## Generate OPLS Parameters for One Ligand

```bash
python ligpargen.py -i ligand_1.pdb -o ligpargen_output
```

---

## Generate Only Specific File Types

```bash
python ligpargen.py \
    -i ligand_1.pdb \
    -o ligpargen_output \
    --only prm,rtf,itp,gro,zip
```

---

## Generate Parameters for Charged Molecules

```bash
python ligpargen.py \
    -i ligand_1.pdb \
    --charge-model cm1a \
    --charge -1
```

---

## List Download Targets Without Saving

```bash
python ligpargen.py -i ligand_1.pdb --no-download
```

---

# PyMOL Movie Generation

Generate movies for protein–ligand complexes.

---

## Create Frame Directory in PyMOL

```python
import os
os.makedirs("frames", exist_ok=True)
```

---

## Generate Frames from Loaded Trajectory

```python
mplay
png frame, width=800, height=600, dpi=300, ray=1
mpng frames/frame
```

---

## Higher Quality Rendering

```python
set ray_trace_frames, 1
mpng frames/frame
```

---

# Convert PNG Frames to GIF Animation

## Install ImageMagick

### Ubuntu/Debian

```bash
sudo apt update
sudo apt install imagemagick
```

---

## Create GIF Animation

```bash
convert -delay 5 -loop 0 frames/frame*.png animation.gif
```

Alternative slower animation:

```bash
convert -delay 10 frames/frame*.png animation.gif
```

---

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
├── frames/
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
