"""Train a QSAR classifier from known (or bootstrap-labeled) actives
and inactives -- MODIFIED version with overfitting-reduction tweaks.

INPUT:
  --actives    CSV with a "smiles" column -- molecules labeled 1 (active)
  --inactives  CSV with a "smiles" column -- molecules labeled 0 (inactive)

  These can come from:
    - Real experimental/literature data for your target, OR
    - bootstrap_sample_dock.py + label_qsar_set.py (docking-derived
      pseudo-labels), which write data/input/actives.csv and
      data/input/inactives.csv automatically.

OUTPUT:
  --output     path to save the trained model (joblib), default
                models/qsar_model.joblib -- this is the file
                stage4_qsar.py loads via qsar.model_path.

MODEL:
  Each molecule's SMILES is converted to a 2048-bit Morgan fingerprint
  (same representation used by stage4_qsar.py). A RandomForestClassifier
  is trained to distinguish actives (1) from inactives (0), with a
  held-out test split and repeated stratified cross-validation reported
  as ROC-AUC.

CHANGES vs train_qsar.py (to reduce overfitting on a small,
imbalanced, high-dimensional dataset -- e.g. 29 actives / 143
inactives, 2048 binary features):

  1. max_depth: 20 -> 8
     Depth-20 trees can carve out a leaf for almost every individual
     training molecule (memorization). Shallower trees are forced to
     find patterns shared by multiple molecules, which generalizes
     better to the held-out test set.

  2. min_samples_leaf: (default 1) -> 3
     Prevents leaves from being created for 1-2 samples, which at this
     sample size is almost always fitting noise rather than signal.

  3. class_weight: (default None) -> "balanced"
     The actives:inactives ratio is roughly 1:5 (from the 10%/50%
     bootstrap split). Without class weighting, a model can reach a
     deceptively high accuracy/AUC just by leaning toward the majority
     (inactive) class. "balanced" reweights each class inversely
     proportional to its frequency.

  4. cross_val_score(cv=5) -> RepeatedStratifiedKFold(n_splits=5, n_repeats=5)
     With ~138 training molecules and a 1:5 class ratio, each fold's
     validation set has only ~28 molecules (~5 actives). A single
     5-fold split can land on an unrepresentatively "easy" or "hard"
     split by chance -- which is the likely reason the original run
     showed CV ROC-AUC=0.907 but Test ROC-AUC=0.668. Repeating the
     5-fold split 5 times (25 train/validate estimates total) gives a
     mean that is much more stable and realistic.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import DataStructs
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score, train_test_split

from utils import mol_from_smiles, morgan_fp


def smiles_to_fps(smiles_list, labels):
    """Convert a list of SMILES + labels into a feature matrix X and
    label vector y, skipping any SMILES that fail to parse.

    Parameters
    ----------
    smiles_list : list[str]
        SMILES strings for both actives and inactives, concatenated.
    labels : list[int]
        Matching labels (1 = active, 0 = inactive).

    Returns
    -------
    X : numpy.ndarray, shape (n_valid, 2048)
        One row per molecule: its 2048-bit Morgan fingerprint as 0/1 ints.
    y : numpy.ndarray, shape (n_valid,)
        The corresponding labels.
    """
    X, y = [], []
    for smi, lab in zip(smiles_list, labels):
        mol = mol_from_smiles(smi)
        if mol is None:
            # Skip invalid SMILES rather than failing the whole run.
            continue
        fp = morgan_fp(mol)
        arr = np.zeros((1,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        X.append(arr)
        y.append(lab)
    return np.vstack(X), np.array(y)


def main(active_file: str, inactive_file: str, model_out: str):
    """Load actives/inactives, train a RandomForest QSAR model, report
    cross-validated and held-out ROC-AUC, and save the model.
    """
    actives = pd.read_csv(active_file)
    inactives = pd.read_csv(inactive_file)

    # Build the combined SMILES list and label vector:
    # label 1 for every active, label 0 for every inactive.
    smiles = list(actives["smiles"]) + list(inactives["smiles"])
    labels = [1] * len(actives) + [0] * len(inactives)

    # Featurize: SMILES -> 2048-bit Morgan fingerprints.
    X, y = smiles_to_fps(smiles, labels)

    # Hold out 20% of the data for a final, unbiased performance check.
    # stratify=y keeps the same actives:inactives ratio in both splits,
    # which matters when the classes are imbalanced (~1:5 here).
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # RandomForest is a robust default for fingerprint-based QSAR:
    # handles high-dimensional sparse binary features well, and
    # predict_proba() gives a probability score (used by stage4_qsar.py).
    #
    # Overfitting-reduction tweaks (see module docstring):
    #   - max_depth=8 (was 20): shallower trees, less memorization
    #   - min_samples_leaf=3: leaves must represent >=3 molecules
    #   - class_weight="balanced": correct for the ~1:5 active:inactive ratio
    model = RandomForestClassifier(
        n_estimators=500,
        max_depth=8,
        min_samples_leaf=3,
        class_weight="balanced",
        n_jobs=-1,
        random_state=42,
    )

    # Repeated stratified 5-fold CV (5 repeats = 25 train/validate splits)
    # on the training set. This gives a much more stable estimate of
    # generalization than a single 5-fold split, especially with only
    # ~138 samples and a 1:5 class ratio (each fold's validation set has
    # just ~28 molecules / ~5 actives -- a single split can be lucky or
    # unlucky by chance).
    cv_splitter = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=42)
    cv = cross_val_score(model, X_train, y_train, cv=cv_splitter, scoring="roc_auc")
    print(f"CV ROC-AUC (5x5 repeated): {cv.mean():.3f} +/- {cv.std():.3f}")

    # Fit on the full training set, then evaluate on the held-out test set.
    model.fit(X_train, y_train)
    auc = roc_auc_score(y_test, model.predict_proba(X_test)[:, 1])
    print(f"Test ROC-AUC: {auc:.3f}")

    # Save the trained model so stage4_qsar.py can load it via
    # joblib.load(qcfg["model_path"]).
    Path(model_out).parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, model_out)
    print(f"Model saved -> {model_out}")


if __name__ == "__main__":
    # Command-line entry point:
    #   python train_qsar_mod.py --actives data/input/actives.csv \
    #                             --inactives data/input/inactives.csv \
    #                             --output models/qsar_model.joblib
    parser = argparse.ArgumentParser()
    parser.add_argument("--actives", required=True, help="CSV with 'smiles' column (known/labeled actives)")
    parser.add_argument("--inactives", required=True, help="CSV with 'smiles' column (decoys/inactives)")
    parser.add_argument("--output", default="models/qsar_model.joblib")
    args = parser.parse_args()
    main(args.actives, args.inactives, args.output)
