"""Stage 5: Prepare ligands and run AutoDock Vina on top QSAR hits.

PURPOSE
-------
After Stage 4 (QSAR scoring) has ranked all 5,825 ADMET-passing
molecules by their predicted binding probability, Stage 5 takes the
top `docking.top_n` molecules (default 100) and performs actual
physics-based molecular docking against the target receptor (FVIIa)
using AutoDock Vina.

Docking is the most computationally expensive step in this pipeline
(each molecule requires seconds to minutes), which is why Stages 1-4
exist: they cheaply filter ~10,000 molecules down to ~100 high-priority
candidates so that docking effort is focused where it is most likely
to succeed.

The best `docking.final_shortlist` molecules (default 34) sorted by
Vina binding-affinity score are written to data/output/final_shortlist.csv
-- the pipeline's final deliverable.

WHAT IS AUTODOCK VINA?
-----------------------
AutoDock Vina is a widely-used open-source molecular docking program.
Given a receptor (protein) and a ligand (small molecule), it:
  1. Explores many possible orientations and conformations (poses) of
     the ligand inside a defined 3D search box (the binding pocket).
  2. Scores each pose using a hybrid empirical/knowledge-based scoring
     function that approximates binding free energy in kcal/mol.
  3. Returns the best poses ranked by score -- MORE NEGATIVE = stronger
     predicted binding.

WHAT IS PDBQT FORMAT?
----------------------
PDBQT (Protein Data Bank, Partial Charge, Atom Type) is the file
format AutoDock Vina requires for both receptor and ligand:
  - Extends the standard PDB format with two extra columns: partial
    atomic charges (from the Gasteiger method) and AutoDock atom types
    (e.g. A=aromatic C, C=aliphatic C, OA=H-bond acceptor O, HD=H-bond
    donor H, N=nitrogen, NA=H-bond acceptor N, SA=S acceptor ...).
  - For ligands, it also encodes the "torsion tree" -- which bonds are
    rotatable and how they connect rigid fragments -- so Vina can
    explore different conformations during docking.

HOW THIS STAGE WORKS (overview)
---------------------------------
For each top-QSAR-ranked molecule:
  1. SMILES -> 3D structure (RDKit: add H, ETKDGv3 embedding, MMFF94
     force-field optimization) -> SDF file.
  2. SDF -> PDBQT (Open Babel: adds Gasteiger charges + torsion tree).
  3. PDBQT ligand + receptor PDBQT -> AutoDock Vina -> best pose score.

Final results are sorted by Vina score (ascending) and the top
`docking.final_shortlist` are written to data/output/final_shortlist.csv.

KEY CONFIG VALUES (config.yaml docking: section)
--------------------------------------------------
  receptor        path to prepared receptor PDBQT (FVIIa)
  center_x/y/z   centroid of the FVIIa binding pocket (Angstroms)
  size_x/y/z     dimensions of the 3D search box (Angstroms)
  exhaustiveness  Vina search effort (higher = slower, more poses tried)
  top_n           how many QSAR-ranked molecules to dock (default 100)
  final_shortlist how many to keep in final output (default 34)

SHARED FUNCTIONS
-----------------
smiles_to_pdbqt() and run_vina() are also imported and reused by
bootstrap_sample_dock.py for the earlier bootstrap docking step (Stage
3 -> Stage 4 training data generation), so any fix here applies there
too.

REQUIREMENTS
------------
  - `vina` (AutoDock Vina) must be on PATH.
  - `obabel` (Open Babel) must be on PATH.
  - `docking.receptor` must point to a valid prepared receptor PDBQT
    (protein only, with H added and Gasteiger charges assigned -- e.g.
    via Meeko: mk_prepare_receptor.py -i receptor.pdb -o receptor -p).
  - `docking.center_x/y/z` and `docking.size_x/y/z` must describe your
    target's ACTUAL binding pocket (not the coordinate origin).
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from utils import load_config, log, mol_from_smiles


# ---------------------------------------------------------------------------
# Function 1: SMILES -> 3D PDBQT
# ---------------------------------------------------------------------------

def smiles_to_pdbqt(smiles: str, out_path: Path) -> bool:
    """Convert a SMILES string into a 3D-conformer PDBQT file for Vina.

    A SMILES string is a flat, 2D text representation of a molecule's
    connectivity (e.g. "CC(=O)Oc1ccccc1C(=O)O" for aspirin). AutoDock
    Vina needs a full 3D structure with explicit hydrogens, partial
    charges, and rotatable-bond information -- this function builds all
    of that from scratch.

    Pipeline inside this function
    ------------------------------
    SMILES
      |
      v  mol_from_smiles()
    RDKit Mol (2D, no H)          -- parses + sanitizes connectivity
      |
      v  Chem.AddHs()
    RDKit Mol (2D, with H)        -- adds all implicit H as explicit atoms
      |
      v  AllChem.EmbedMolecule(ETKDGv3)
    RDKit Mol (3D, with H)        -- places every atom in 3D space
      |
      v  AllChem.MMFFOptimizeMolecule()
    RDKit Mol (3D, optimized)     -- relaxes geometry via MMFF94 force field
      |
      v  Chem.MolToMolFile()
    Temporary .sdf file           -- standard chemistry 3D file format
      |
      v  obabel ... --gen3d
    Final .pdbqt file             -- AutoDock atom types + charges + torsion tree
      |
      (delete .sdf)

    Parameters
    ----------
    smiles : str
        A canonical SMILES string (already sanitized by Stages 1-3).
    out_path : Path
        Where to write the output .pdbqt file.

    Returns
    -------
    bool
        True on success.
        False if:
          - SMILES could not be parsed (mol_from_smiles returns None).
          - 3D embedding failed (EmbedMolecule returns -1) -- can happen
            for molecules with unusual ring systems or very large flexible
            chains.
    """
    # --- Step 1: Parse the SMILES into an RDKit Mol object ---
    # mol_from_smiles() calls Chem.MolFromSmiles + Chem.SanitizeMol,
    # returning None if the SMILES is malformed. At this point in the
    # pipeline the SMILES have already been canonicalized by Stage 1,
    # but we guard anyway for robustness.
    mol = mol_from_smiles(smiles)
    if mol is None:
        return False

    # --- Step 2: Add explicit hydrogen atoms ---
    # RDKit normally stores H implicitly (just counts them per heavy atom).
    # Chem.AddHs() makes every H an explicit atom with its own position,
    # which is required for:
    #   (a) realistic 3D geometry (H-H and X-H bond lengths/angles)
    #   (b) correct H-bond donor/acceptor assignment in Vina
    #   (c) Gasteiger charge computation by Open Babel
    mol = Chem.AddHs(mol)

    # --- Step 3: Generate a 3D conformation (ETKDGv3) ---
    # EmbedMolecule() places every atom in 3D space using the
    # ETKDGv3 (Experimental-Torsion basic Knowledge Distance Geometry
    # version 3) algorithm, which uses experimental torsion-angle
    # preferences from the CSD and PDB to produce realistic
    # starting geometries.
    #
    # Return value: the conformer's integer ID (>= 0) on success,
    # or -1 if embedding failed. We MUST check for -1 before calling
    # MMFFOptimizeMolecule -- passing a mol with no conformer raises
    # "ValueError: Bad Conformer Id".
    conf_id = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if conf_id < 0:
        # Embedding can fail for macrocycles, highly strained systems,
        # or very large flexible molecules where distance-geometry
        # constraints cannot be satisfied.
        return False

    # --- Step 4: Optimize the 3D geometry (MMFF94 force field) ---
    # The ETKDGv3 geometry is a reasonable starting point but not a
    # local energy minimum. MMFFOptimizeMolecule() runs a short
    # molecular mechanics minimization using the MMFF94 force field,
    # which brings bond lengths, angles, and torsions to more realistic
    # values and removes any steric clashes introduced by embedding.
    #
    # Wrapped in try/except because MMFF94 has no parameters for some
    # exotic atom types (certain metals, unusual halogens) -- in those
    # cases we keep the unoptimized ETKDGv3 geometry rather than failing.
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except ValueError:
        pass

    # --- Step 5: Write the 3D structure to a temporary SDF file ---
    # SDF (Structure Data File) is a standard chemistry format that
    # stores the 3D atom coordinates and bond table. We use it as an
    # intermediate because Open Babel reads SDF reliably and can add
    # the AutoDock-specific information from it.
    sdf = out_path.with_suffix(".sdf")
    Chem.MolToMolFile(mol, str(sdf))

    # --- Step 6: Convert SDF -> PDBQT via Open Babel ---
    # `obabel <input.sdf> -O <output.pdbqt> --gen3d` does three things:
    #   (a) Re-generates 3D coordinates if needed (here they already exist).
    #   (b) Assigns Gasteiger partial charges to every atom (required by Vina
    #       for its electrostatics scoring term).
    #   (c) Detects rotatable bonds and writes the torsion tree into the
    #       PDBQT header (ROOT/BRANCH/ENDBRANCH records), which tells Vina
    #       which bonds can rotate during docking.
    # check=True raises CalledProcessError if obabel returns a non-zero
    # exit code; capture_output=True suppresses Open Babel's verbose
    # stdout/stderr from appearing in the terminal.
    subprocess.run(
        ["obabel", str(sdf), "-O", str(out_path), "--gen3d"],
        check=True,
        capture_output=True,
    )

    # --- Step 7: Remove the intermediate SDF file ---
    # It is no longer needed now that the PDBQT exists; missing_ok=True
    # means no error if it was already deleted or never created.
    sdf.unlink(missing_ok=True)
    return True


# ---------------------------------------------------------------------------
# Function 2: run one AutoDock Vina docking job
# ---------------------------------------------------------------------------

def run_vina(receptor: str, ligand_pdbqt: Path, out_pdbqt: Path,
              center, size, exhaustiveness: int = 8) -> float | None:
    """Run AutoDock Vina for one ligand against one receptor and return
    the best pose's binding-affinity score.

    What Vina does internally
    --------------------------
    1. Reads the receptor PDBQT (protein with charges + atom types).
    2. Reads the ligand PDBQT (small molecule with charges + torsion tree).
    3. Defines a 3D search box (center + size) -- Vina only explores
       ligand positions inside this box.
    4. Uses a stochastic global search (iterated local search with
       Metropolis-style moves) to sample ligand poses inside the box,
       guided by the scoring function.
    5. Refines the best poses with local optimization.
    6. Prints a ranked table of the top poses to stdout.

    The Vina stdout table looks like:
    -------+--------+----------+----------
       #   |   mode |  affinity| dist from best mode
            |        | (kcal/mol)| rmsd l.b.| rmsd u.b.
    -------+--------+----------+----------
       1    |   1    |   -7.234 |   0.000  |   0.000
       2    |   2    |   -7.011 |   1.423  |   2.117
       ...
    We extract the score from pose "1" (the best pose).

    Parameters
    ----------
    receptor : str
        Path to the receptor PDBQT (e.g. "data/receptor.pdbqt").
        Prepared with Meeko (mk_prepare_receptor.py) for FVIIa.
    ligand_pdbqt : Path
        Path to the ligand PDBQT produced by smiles_to_pdbqt().
    out_pdbqt : Path
        Where Vina writes the docked pose(s) in PDBQT format.
        These can be visualized later in PyMOL or UCSF Chimera.
    center : tuple(float, float, float)
        (center_x, center_y, center_z) -- the centroid of the FVIIa
        binding pocket in Angstroms, from config.yaml.
        Currently: (29.90, 14.01, 7.89) derived from residue-based
        pocket calculation.
    size : tuple(float, float, float)
        (size_x, size_y, size_z) -- dimensions of the search box in
        Angstroms. Currently: (25.6, 24.9, 23.1).
        Rule of thumb: the box should be large enough to enclose the
        entire binding pocket + ~4 Angstrom buffer on each side.
    exhaustiveness : int
        Controls how thoroughly Vina searches the pose space. Default
        is 8 (Vina's own default). Higher values (16, 32) improve
        reliability but linearly increase runtime.

    Returns
    -------
    float or None
        The binding affinity of the best pose in kcal/mol (negative
        value; more negative = stronger predicted binding).
        Returns None if Vina's stdout contained no parseable result
        (e.g. Vina crashed, ligand is outside the box, or the box is
        too small for the molecule).
    """
    # Unpack center and size tuples into individual x/y/z values
    # because Vina takes them as separate command-line arguments.
    cx, cy, cz = center
    sx, sy, sz = size

    # Build the full Vina command line. Every parameter maps directly
    # to a config.yaml entry so the same command works for any target
    # just by editing config.yaml.
    cmd = [
        "vina",
        "--receptor",    receptor,          # FVIIa PDBQT
        "--ligand",      str(ligand_pdbqt), # small molecule PDBQT
        "--out",         str(out_pdbqt),    # save docked poses here
        "--center_x",    str(cx),           # pocket centroid X
        "--center_y",    str(cy),           # pocket centroid Y
        "--center_z",    str(cz),           # pocket centroid Z
        "--size_x",      str(sx),           # search box X dimension
        "--size_y",      str(sy),           # search box Y dimension
        "--size_z",      str(sz),           # search box Z dimension
        "--exhaustiveness", str(exhaustiveness),
    ]

    # Run Vina as a subprocess, capturing all output.
    # text=True decodes stdout/stderr from bytes to str automatically.
    # We do NOT use check=True here -- a Vina failure for one molecule
    # should return None rather than raising an exception and aborting
    # the entire docking run.
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse Vina's stdout to extract the best pose score.
    # The results table starts after a header block; each data row
    # begins with the pose number. We look for the line that starts
    # with "1 " (pose 1 = best score), split it, and return column [1]
    # (the affinity value, e.g. "-7.234").
    for line in result.stdout.splitlines():
        if line.strip().startswith("1 "):
            return float(line.split()[1])

    # Vina produced no usable score -- log is handled by the caller.
    return None


# ---------------------------------------------------------------------------
# Function 3: main Stage 5 pipeline
# ---------------------------------------------------------------------------

def run(cfg: dict):
    """Dock the top QSAR-ranked molecules and write the final shortlist.

    Reads Stage 4's output (data/intermediate/stage4/qsar_ranked.parquet),
    docks the top `docking.top_n` molecules with AutoDock Vina, sorts
    the results by Vina score (ascending = best first), and writes the
    top `docking.final_shortlist` molecules to
    data/output/final_shortlist.csv.

    Parameters
    ----------
    cfg : dict
        The full parsed config.yaml. The "docking:" section is used for
        all Vina parameters; "processing.output_dir" for file paths.

    Returns
    -------
    pandas.DataFrame
        Final shortlist with columns:
          id          -- molecule identifier (PubChem CID)
          smiles      -- canonical SMILES
          qsar_score  -- Stage 4 QSAR probability (0-1, higher = more active)
          vina_score  -- AutoDock Vina binding affinity (kcal/mol, more negative = better)
        Sorted by vina_score ascending, capped at docking.final_shortlist rows.
    """
    # Load the "docking:" section of config.yaml into a short alias.
    dcfg = cfg["docking"]

    # -----------------------------------------------------------------
    # Step 1: Load Stage 4's QSAR-ranked output and select top N
    # -----------------------------------------------------------------
    # Stage 4 already sorted the molecules by qsar_score descending and
    # capped at qsar.top_n=500. Here we take only the top
    # docking.top_n=100 of those for actual docking, because docking
    # all 500 would take much longer with limited additional benefit
    # (the highest-confidence QSAR hits are most likely to be real
    # binders and should be prioritized for expensive docking).
    infile = Path(cfg["processing"]["output_dir"]) / "stage4" / "qsar_ranked.parquet"
    df = pd.read_parquet(infile).head(dcfg["top_n"])
    log.info("Stage 5 docking: %d molecules to dock (top_n=%d)", len(df), dcfg["top_n"])

    # -----------------------------------------------------------------
    # Step 2: Prepare output directories
    # -----------------------------------------------------------------
    # ligands/ -- stores the per-molecule PDBQT files used as Vina input
    # results/ -- stores Vina's output PDBQT files (docked pose coordinates)
    # Both directories are created if they don't already exist.
    lig_dir = Path("data/intermediate/docking/ligands")
    out_dir = Path("data/intermediate/docking/results")
    lig_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------
    # Step 3: Read binding pocket box from config
    # -----------------------------------------------------------------
    # These six values define WHERE in the receptor Vina should look for
    # binding poses. They must correspond to the actual FVIIa binding
    # pocket -- wrong values (e.g. the default 0,0,0) would cause Vina
    # to dock into empty space and return meaningless scores.
    # Values set in config.yaml via compute_pocket_box_blind.sh.
    center = (dcfg["center_x"], dcfg["center_y"], dcfg["center_z"])  # (29.90, 14.01, 7.89)
    size   = (dcfg["size_x"],   dcfg["size_y"],   dcfg["size_z"])    # (25.6, 24.9, 23.1)

    # -----------------------------------------------------------------
    # Step 4: Dock each molecule
    # -----------------------------------------------------------------
    scores = []
    n = len(df)

    for i, (_, row) in enumerate(df.iterrows(), start=1):
        try:
            # -- 4a. Build ligand PDBQT from SMILES --
            # The PDBQT file name is based on the molecule's ID so each
            # molecule has its own file and runs can be inspected/resumed.
            lig = lig_dir / f"{row['id']}.pdbqt"
            if not smiles_to_pdbqt(row["smiles"], lig):
                # smiles_to_pdbqt returns False for embedding failures or
                # unparseable SMILES -- skip and continue with next molecule.
                log.warning("Skipping %s: could not generate 3D structure", row["id"])
                continue

            # -- 4b. Run AutoDock Vina --
            # out_pdbqt stores the docked pose coordinates and can be
            # opened in PyMOL/Chimera to visualize how the molecule sits
            # in the FVIIa binding pocket.
            out = out_dir / f"{row['id']}_out.pdbqt"
            score = run_vina(
                dcfg["receptor"], lig, out,
                center, size,
                dcfg["exhaustiveness"],
            )

            if score is None:
                # Vina returned no score -- molecule may be outside the
                # search box or Vina encountered an internal error.
                log.warning("Skipping %s: Vina returned no score", row["id"])
                continue

            scores.append({
                "id":         row["id"],
                "smiles":     row["smiles"],
                "qsar_score": row["qsar_score"],
                "vina_score": score,
            })
            log.info(
                "Docked %d/%d  id=%s  qsar=%.3f  vina=%.2f kcal/mol",
                i, n, row["id"], row["qsar_score"], score,
            )

        except Exception as exc:
            # Catch-all so that one unexpected failure (e.g. obabel
            # error on an unusual structure) does not abort the entire
            # run.
            log.warning("Skipping %s: unexpected error -- %s", row["id"], exc)
            continue

    # -----------------------------------------------------------------
    # Step 5: Assemble results DataFrame
    # -----------------------------------------------------------------
    result = pd.DataFrame(scores)

    if result.empty:
        log.warning("Stage 5: no molecules were successfully docked.")
        return result

    # dropna removes any rows where vina_score is NaN (should not happen
    # given the None check above, but defensive).
    result = result.dropna(subset=["vina_score"])

    # -----------------------------------------------------------------
    # Step 6: Sort by Vina score (ascending = best binding first)
    # -----------------------------------------------------------------
    # Vina scores are in kcal/mol and are NEGATIVE for predicted binders.
    # Sorting ascending puts the most negative (strongest predicted
    # binder) at the top -- this is the final ranking.
    result = result.sort_values("vina_score")

    # -----------------------------------------------------------------
    # Step 7: Truncate to the final shortlist and save
    # -----------------------------------------------------------------
    # docking.final_shortlist (default 34) is the number of molecules
    # the downstream MD simulation steps can realistically handle.
    # These are the pipeline's final hits for further experimental
    # validation or MD simulation.
    final = result.head(dcfg["final_shortlist"])

    out_path = Path("data/output/final_shortlist.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(out_path, index=False)

    log.info(
        "Stage 5 complete: %d docked / %d attempted -> final shortlist %d molecules -> %s",
        len(result), n, len(final), out_path,
    )
    return final


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Run Stage 5 from the command line:
    #   python stage5_docking.py --config config.yaml
    #
    # Output files:
    #   data/intermediate/docking/ligands/<id>.pdbqt   -- ligand input PDBQTs
    #   data/intermediate/docking/results/<id>_out.pdbqt -- Vina pose outputs
    #   data/output/final_shortlist.csv                 -- final ranked hits
    parser = argparse.ArgumentParser(
        description="Stage 5: dock top QSAR hits against receptor with AutoDock Vina."
    )
    parser.add_argument(
        "--config", default="config.yaml",
        help="Path to config.yaml (default: config.yaml)"
    )
    args = parser.parse_args()
    run(load_config(args.config))
