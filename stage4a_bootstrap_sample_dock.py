"""Bootstrap step (between Stage 3 and Stage 4): generate labeled
training data for the QSAR model by docking a random sample of
ADMET-passing molecules.

WHY THIS EXISTS
---------------
Stage 4 (stage4_qsar.py) needs a pre-trained QSAR model
(qsar.model_path). Training that model (train_qsar.py) needs known
"active" and "inactive" molecules with SMILES. If you don't already
have experimentally-validated actives/inactives for your target, this
script generates a reasonable substitute by:

  1. Taking a random sample of `bootstrap.sample_size` molecules from
     Stage 3's output (data/intermediate/stage3/admet_passed.parquet).
  2. Docking each one with AutoDock Vina against your receptor
     (re-using smiles_to_pdbqt/run_vina from stage5_docking.py).
  3. Saving every molecule's Vina binding-affinity score.

The companion script label_qsar_set.py then turns these scores into
actives.csv / inactives.csv (by percentile thresholds), which
train_qsar.py consumes to produce qsar.model_path. Stage 4 can then
score the REST of the library quickly (without docking everything).

REQUIREMENTS:
  - `vina` and `obabel` on PATH.
  - docking.receptor pointing to a prepared receptor PDBQT.
  - docking.center_x/y/z + size_x/y/z set to your real binding pocket.

This step is slow (one Vina run per sampled molecule), which is why
the sample size is kept small (default 300) relative to the full
library (thousands of molecules).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from stage5_docking import run_vina, smiles_to_pdbqt
from utils import load_config, log, save_parquet


def run(cfg: dict):
    """Dock a random sample of Stage 3 molecules and record their
    Vina binding-affinity scores.

    Returns
    -------
    pandas.DataFrame
        Columns: id, smiles, vina_score (one row per successfully
        docked molecule). Also written to
        data/intermediate/bootstrap/docked_sample.parquet.
    """
    dcfg = cfg["docking"]
    bcfg = cfg["bootstrap"]

    # Read the ADMET-passing molecules produced by Stage 3.
    infile = Path(cfg["processing"]["output_dir"]) / "stage3" / "admet_passed.parquet"
    df = pd.read_parquet(infile)

    # Take a reproducible random sample (capped at the available rows).
    n = min(bcfg["sample_size"], len(df))
    sample = df.sample(n=n, random_state=bcfg["seed"]).reset_index(drop=True)
    log.info("Bootstrap docking: sampled %d / %d Stage 3 molecules", n, len(df))

    # Directories for intermediate ligand files and Vina output poses.
    lig_dir = Path("data/intermediate/bootstrap/ligands")
    out_dir = Path("data/intermediate/bootstrap/results")
    lig_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Binding pocket search box, read from config.
    center = (dcfg["center_x"], dcfg["center_y"], dcfg["center_z"])
    size = (dcfg["size_x"], dcfg["size_y"], dcfg["size_z"])

    scores = []
    for i, row in sample.iterrows():
        try:
            # Step 1: build a 3D ligand PDBQT for this molecule.
            lig = lig_dir / f"{row['id']}.pdbqt"
            if not smiles_to_pdbqt(row["smiles"], lig):
                log.warning("Skipping %s: could not generate 3D structure", row["id"])
                continue

            # Step 2: dock the ligand into the receptor's binding pocket.
            out = out_dir / f"{row['id']}_out.pdbqt"
            score = run_vina(dcfg["receptor"], lig, out, center, size, dcfg["exhaustiveness"])
            if score is None:
                log.warning("Skipping %s: Vina returned no score", row["id"])
                continue

            scores.append({"id": row["id"], "smiles": row["smiles"], "vina_score": score})
            log.info("Docked %d/%d  id=%s  vina_score=%.2f", i + 1, n, row["id"], score)
        except Exception as exc:
            # A single problematic molecule (e.g. obabel choking on an
            # unusual structure) shouldn't abort a 300-molecule run.
            log.warning("Skipping %s: %s", row["id"], exc)
            continue

    result = pd.DataFrame(scores)

    out_path = Path(cfg["processing"]["output_dir"]) / "bootstrap" / "docked_sample.parquet"
    save_parquet(result, out_path)
    log.info("Bootstrap docking complete: %d / %d molecules docked successfully", len(result), n)
    return result


if __name__ == "__main__":
    # Command-line entry point:
    #   python bootstrap_sample_dock.py --config config.yaml
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    run(load_config(args.config))
