"""Stage 3: Rule-based ADMET filtering (Discovery Studio-style).

This stage applies a SECOND, tighter set of physicochemical filters on
top of the survivors from Stage 2 (data/intermediate/stage2/filtered.parquet).

Unlike Stage 2 (which computes descriptors fresh with RDKit and exits
early on the first failed rule), this stage:
  - Re-uses the descriptor columns already computed in Stage 2
    (MW, LogP, HBD, HBA, TPSA, RotBonds, AromaticRings, QED).
  - Requires a molecule to satisfy ALL conditions simultaneously
    (no early exit) before being kept.
  - Adds range checks (min AND max) for LogP and TPSA, and a minimum
    QED (drug-likeness) score, as rough ADMET (Absorption,
    Distribution, Metabolism, Excretion, Toxicity) proxies.

Output: data/intermediate/stage3/admet_passed.parquet, consumed by
Stage 4 (QSAR scoring).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from utils import load_config, log, save_parquet


def admet_pass(row, acfg) -> bool:
    """Check a single molecule's descriptors against ADMET thresholds.

    Parameters
    ----------
    row : pandas.Series
        One row from the Stage 2 output, containing pre-computed
        descriptor columns (LogP, TPSA, MW, HBD, HBA, RotBonds, QED,
        AromaticRings).
    acfg : dict
        The "admet" section of config.yaml, supplying the thresholds.

    Returns
    -------
    bool
        True only if EVERY check below passes (all() requires every
        element to be True). Unlike Stage 2, there is no early exit:
        all 8 conditions are evaluated for every row.
    """
    checks = [
        # LogP (lipophilicity) must fall within a band, not just below
        # a ceiling -- molecules that are too hydrophilic (low LogP)
        # can also have poor membrane permeability / absorption.
        acfg["min_logp"] <= row["LogP"] <= acfg["max_logp"],

        # TPSA (topological polar surface area) range check.
        # Used here as a rough oral-absorption proxy: too low TPSA can
        # indicate poor solubility, too high can indicate poor
        # membrane permeability.
        acfg["min_tpsa"] <= row["TPSA"] <= acfg["max_tpsa"],

        # Molecular weight ceiling (often stricter than Stage 2's,
        # e.g. 450 vs 500) for better ADMET properties.
        row["MW"] <= acfg["max_mw"],

        # Hydrogen bond donor count ceiling.
        row["HBD"] <= acfg["max_hbd"],

        # Hydrogen bond acceptor count ceiling.
        row["HBA"] <= acfg["max_hba"],

        # Rotatable bond count ceiling (flexibility / bioavailability).
        row["RotBonds"] <= acfg["max_rotatable_bonds"],

        # QED (Quantitative Estimate of Drug-likeness, 0-1 scale)
        # must meet a MINIMUM score -- filters out molecules RDKit's
        # composite drug-likeness model scores poorly.
        row["QED"] >= acfg["min_qed"],

        # Aromatic ring count ceiling -- too many aromatic rings is
        # associated with poor solubility / higher attrition risk.
        row["AromaticRings"] <= acfg["max_aromatic_rings"],
    ]
    # all() => every condition in `checks` must be True for the
    # molecule to pass. A single False anywhere causes rejection.
    return all(checks)


def run(cfg: dict):
    """Run Stage 3 ADMET filtering over the Stage 2 output.

    Parameters
    ----------
    cfg : dict
        The full parsed config.yaml.

    Returns
    -------
    pandas.DataFrame
        The subset of Stage 2's molecules that also pass every ADMET
        rule, including their original descriptor columns.
    """
    # Stage 3 reads directly from Stage 2's output -- no re-parsing of
    # SMILES or recomputation of descriptors is needed.
    infile = Path(cfg["processing"]["output_dir"]) / "stage2" / "filtered.parquet"
    df = pd.read_parquet(infile)

    # Load the ADMET threshold dict from config.yaml's "admet:" section.
    acfg = cfg["admet"]

    # Apply admet_pass() to every row (axis=1 => row-wise) to build a
    # boolean Series the same length as df. This is single-threaded
    # (no multiprocessing) since each check is cheap arithmetic on
    # numbers that are already computed.
    mask = df.apply(lambda r: admet_pass(r, acfg), axis=1)

    # Keep only the rows where mask is True. .copy() avoids
    # SettingWithCopyWarning if `passed` is modified later.
    passed = df[mask].copy()

    # Write the surviving molecules (with all their descriptor columns
    # intact) to Stage 3's output file.
    out = Path(cfg["processing"]["output_dir"]) / "stage3" / "admet_passed.parquet"
    save_parquet(passed, out)

    # Report how many molecules passed out of how many were checked.
    log.info("Stage 3 ADMET: %d / %d passed", len(passed), len(df))
    return passed


if __name__ == "__main__":
    # Command-line entry point:
    #   python stage3_admet.py --config config.yaml
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    run(load_config(args.config))
