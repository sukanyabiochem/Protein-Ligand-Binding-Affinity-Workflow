"""Stage 6: QSAR model scoring with Morgan fingerprints.

(Note: internally labeled "Stage 6" from the original pipeline naming,
but in THIS pipeline it is the 4th script run, i.e. "Stage 4".)

This stage takes the ADMET-passing molecules from Stage 3
(data/intermediate/stage3/admet_passed.parquet), converts each one to
a Morgan (ECFP-style) fingerprint, and scores it with a pre-trained
QSAR machine-learning model (loaded via joblib).

Molecules are then:
  1. Filtered to keep only those with qsar_score >= qsar.min_probability
  2. Sorted by score, descending
  3. Truncated to the top qsar.top_n molecules

Output: data/intermediate/stage4/qsar_ranked.parquet, consumed by
Stage 5 (docking).

REQUIREMENTS TO RUN THIS STAGE:
  - config.yaml must have a "qsar:" section with model_path,
    min_probability, and top_n.
  - A trained model file must exist at qsar.model_path (a
    scikit-learn-compatible estimator saved with joblib.dump(),
    accepting a 2048-bit Morgan fingerprint as input). This is
    normally produced by train_qsar.py (run separately, once, before
    this stage).
"""
from __future__ import annotations

import argparse
from multiprocessing import Pool
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import DataStructs

from utils import load_config, log, mol_from_smiles, morgan_fp, save_parquet


def _fp_row(row):
    """Compute a Morgan fingerprint feature vector for one molecule.

    Parameters
    ----------
    row : dict
        A single record with "id" and "smiles" keys (one row from the
        Stage 3 output).

    Returns
    -------
    dict or None
        {"id": ..., "smiles": ..., "fp": numpy array of 0/1, length 2048}
        or None if the SMILES could not be parsed.
    """
    # Step 1: re-parse the canonical SMILES into an RDKit Mol object.
    mol = mol_from_smiles(row["smiles"])
    if mol is None:
        # Should be rare at this point -- the molecule already passed
        # Stages 1-3 -- but guard against edge cases anyway.
        return None

    # Step 2: compute the Morgan (ECFP-like) circular fingerprint.
    # utils.morgan_fp() defaults to radius=2, n_bits=2048, returning an
    # RDKit ExplicitBitVect (a compact bit-vector object).
    fp = morgan_fp(mol)

    # Step 3: convert the RDKit bit-vector into a numpy array so it can
    # be fed into a scikit-learn model. ConvertToNumpyArray resizes
    # `arr` in-place to match the fingerprint length (2048), filling it
    # with 0s and 1s.
    arr = np.zeros((1,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)

    return {"id": row["id"], "smiles": row["smiles"], "fp": arr}


def run(cfg: dict):
    """Run Stage 4 QSAR scoring over the Stage 3 output.

    Parameters
    ----------
    cfg : dict
        The full parsed config.yaml. Must contain a "qsar" section
        (model_path, min_probability, top_n).

    Returns
    -------
    pandas.DataFrame
        The top-ranked molecules (id, smiles, qsar_score) that meet
        or exceed the configured probability threshold, sorted
        descending by score and capped at top_n rows.
    """
    # Stage 4 reads directly from Stage 3's output.
    infile = Path(cfg["processing"]["output_dir"]) / "stage3" / "admet_passed.parquet"
    df = pd.read_parquet(infile)

    # Load the "qsar:" config section (model_path, min_probability, top_n).
    qcfg = cfg["qsar"]

    # Load the pre-trained QSAR model from disk. This must be a
    # scikit-learn-compatible object (has .predict, optionally
    # .predict_proba) saved with joblib.dump() by train_qsar.py.
    model = joblib.load(qcfg["model_path"])

    log.info("Stage 6 QSAR: scoring %d molecules", len(df))

    # ------------------------------------------------------------
    # Step 1: compute Morgan fingerprints for every molecule, in
    # parallel across n_workers processes.
    # ------------------------------------------------------------
    n_workers = cfg["processing"]["n_workers"]
    with Pool(n_workers) as pool:
        fps = pool.map(_fp_row, df.to_dict("records"))

    # Drop any molecules whose SMILES failed to parse (returns None).
    valid = [f for f in fps if f is not None]

    # Stack the individual 2048-length fingerprint vectors into a single
    # 2D feature matrix X of shape (n_molecules, 2048) for the model.
    X = np.vstack([f["fp"] for f in valid])

    # ------------------------------------------------------------
    # Step 2: score every molecule with the trained model.
    # ------------------------------------------------------------
    if hasattr(model, "predict_proba"):
        # Classifier: take the probability of the POSITIVE class
        # (index 1, assumed to mean "active" / "desired hit").
        scores = model.predict_proba(X)[:, 1]
    else:
        # Regressor (or any model without predict_proba): use the
        # raw prediction directly as the score.
        scores = model.predict(X)

    # ------------------------------------------------------------
    # Step 3: assemble results, filter by threshold, rank, and cap.
    # ------------------------------------------------------------
    result = pd.DataFrame({
        "id": [f["id"] for f in valid],
        "smiles": [f["smiles"] for f in valid],
        "qsar_score": scores,
    })

    # Keep only molecules whose score meets the configured minimum
    # probability (e.g. 0.6 = "likely active").
    result = result[result["qsar_score"] >= qcfg["min_probability"]]

    # Sort by score, best first, and keep only the top N molecules
    # (e.g. to cap the list before the expensive docking stage).
    result = result.sort_values("qsar_score", ascending=False).head(qcfg["top_n"])

    # Write the final ranked shortlist to Stage 4's output file.
    out = Path(cfg["processing"]["output_dir"]) / "stage4" / "qsar_ranked.parquet"
    save_parquet(result, out)

    log.info("Stage 6 complete: %d molecules above threshold", len(result))
    return result


if __name__ == "__main__":
    # Command-line entry point:
    #   python stage4_qsar.py --config config.yaml
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    run(load_config(args.config))
