"""Bootstrap step (between bootstrap_sample_dock.py and train_qsar.py):
turn docking scores from the bootstrap sample into labeled
actives.csv / inactives.csv files for train_qsar.py.

WHY THIS EXISTS
---------------
bootstrap_sample_dock.py produces a Vina binding-affinity score
(more negative = stronger predicted binding) for a random sample of
ADMET-passing molecules. This script ranks that sample by score and
splits it into two groups by percentile:

  - "actives"   = the best-scoring `bootstrap.active_pct` percent
                  (most negative vina_score) -> label 1
  - "inactives" = the worst-scoring `bootstrap.inactive_pct` percent
                  (least negative / weakest binding) -> label 0

These are written to data/input/actives.csv and
data/input/inactives.csv (each with an "id" and "smiles" column,
matching the format train_qsar.py expects).

NOTE: This produces PSEUDO-labels based on a single docking score per
molecule (Vina scores have significant noise/uncertainty). The
resulting QSAR model is a heuristic pre-filter to reduce how many
molecules need full docking in Stage 5 -- it is not a substitute for
experimentally-validated activity data.

Explains the role of this script in the pipeline: it sits between the bootstrap docking step and QSAR training, splits the docked sample by percentile of vina_score, and explicitly flags that these are pseudo-labels (noisy, score-derived — not real assay data).

"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from utils import load_config, log

##- pandas for reading/writing the parquet and CSV files.
##- load_config, log from utils.py — same config-loading and logging infrastructure used by every stage.


def run(cfg: dict):
    """Split the bootstrap-docked sample into actives/inactives CSVs.

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        The actives and inactives DataFrames that were written out.
    """
    bcfg = cfg["bootstrap"]

    infile = Path(cfg["processing"]["output_dir"]) / "bootstrap" / "docked_sample.parquet"
    df = pd.read_parquet(infile)

    # Sort by Vina score ascending: most negative (best binding affinity) first.
    df = df.sort_values("vina_score").reset_index(drop=True)
    n = len(df)

    # Number of rows in the "active" and "inactive" groups, based on
    # the configured percentiles. At least 1 row each (if n > 0).
    n_active = max(1, round(n * bcfg["active_pct"] / 100))
    n_inactive = max(1, round(n * bcfg["inactive_pct"] / 100))

    # Best-scoring (most negative vina_score) rows -> presumed actives.
    actives = df.head(n_active)

    # Worst-scoring (least negative / weakest binding) rows -> presumed
    # inactives. Guaranteed not to overlap with `actives` as long as
    # n_active + n_inactive <= n (true for the default 10% / 50% split).
    inactives = df.tail(n_inactive)

    out_dir = Path("data/input")
    out_dir.mkdir(parents=True, exist_ok=True)

    actives_path = out_dir / "actives.csv"
    inactives_path = out_dir / "inactives.csv"
    actives[["id", "smiles"]].to_csv(actives_path, index=False)
    inactives[["id", "smiles"]].to_csv(inactives_path, index=False)

    log.info(
        "Wrote %d actives (vina_score <= %.2f) -> %s",
        len(actives), actives["vina_score"].max(), actives_path,
    )
    log.info(
        "Wrote %d inactives (vina_score >= %.2f) -> %s",
        len(inactives), inactives["vina_score"].min(), inactives_path,
    )
    return actives, inactives


if __name__ == "__main__":
    # Command-line entry point:
    #   python label_qsar_set.py --config config.yaml
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    run(load_config(args.config))
