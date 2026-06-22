"""Visualize the QSAR training/testing dataset and the trained model's
performance.

This script re-loads the same actives.csv / inactives.csv used by
train_qsar_mod.py, recreates the IDENTICAL train/test split (same
test_size, random_state, and stratify=y), loads the saved model
(models/qsar_model.joblib), and produces four diagnostic plots:

  1. roc_curves.png
     ROC curves for repeated 5-fold cross-validation on the training
     set (5 folds x 5 repeats = 25 splits, light blue, plus the mean)
     overlaid with the held-out test set ROC curve (red). This matches
     train_qsar_mod.py's "CV ROC-AUC (5x5 repeated)" reporting, so the
     plotted CV mean should match the printed 0.843 +/- 0.085. Lets
     you visually compare the CV estimate against the actual test
     performance -- e.g. confirm that the large CV/Test gap seen with
     the original hyperparameters (CV 0.907 vs Test 0.668) is gone
     with the train_qsar_mod.py settings (CV ~0.84 vs Test 0.885).

  2. chemical_space_pca.png
     A 2D PCA projection of the 2048-bit Morgan fingerprints for every
     molecule in the dataset, colored by label (active=red,
     inactive=blue) and marked by split (circle=train, triangle=test).
     Shows whether "active" and "inactive" pseudo-labels occupy
     distinguishable regions of chemical (fingerprint) space.

  3. probability_distributions.png
     Histograms of the model's predicted P(active) for the training
     set and the test set, split by true label. A good model pushes
     true actives toward 1.0 and true inactives toward 0.0. The
     dashed line marks qsar.min_probability=0.6 (the cutoff
     stage4_qsar.py uses to keep/discard molecules).

  4. confusion_matrix.png
     Confusion matrix on the held-out test set at the
     qsar.min_probability threshold (default 0.6), showing counts of
     true/false positives and negatives.

USAGE:
    python plot_qsar_results.py \
        --actives data/input/actives.csv \
        --inactives data/input/inactives.csv \
        --model models/qsar_model.joblib \
        --outdir results/qsar_plots \
        --threshold 0.6

REQUIREMENTS:
    matplotlib (if missing: conda install -c conda-forge matplotlib)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import matplotlib

# Use a non-interactive backend -- this script only saves PNG files,
# it never needs to open a window (works on headless servers too).
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from sklearn.base import clone
from sklearn.decomposition import PCA
from sklearn.metrics import ConfusionMatrixDisplay, auc, confusion_matrix, roc_curve
from sklearn.model_selection import RepeatedStratifiedKFold, train_test_split

# Re-use the exact same featurization (SMILES -> 2048-bit Morgan
# fingerprint matrix) that produced the saved model, so the plots are
# computed on identical features.
from train_qsar_mod import smiles_to_fps
from utils import log

import pandas as pd


def load_dataset(active_file: str, inactive_file: str):
    """Load actives/inactives CSVs and featurize them.

    Returns
    -------
    X : numpy.ndarray, shape (n_valid, 2048)
    y : numpy.ndarray, shape (n_valid,)
        1 = active, 0 = inactive.
    """
    actives = pd.read_csv(active_file)
    inactives = pd.read_csv(inactive_file)
    smiles = list(actives["smiles"]) + list(inactives["smiles"])
    labels = [1] * len(actives) + [0] * len(inactives)
    return smiles_to_fps(smiles, labels)


def plot_roc_cv_and_test(model, X_train, y_train, X_test, y_test, outdir: Path):
    """Plot 5-fold CV ROC curves (training set) overlaid with the
    held-out test set ROC curve from the already-fitted model.
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    # --- Repeated 5-fold CV ROC curves on the training set ---
    # Same scheme as train_qsar_mod.py's reported "CV ROC-AUC (5x5
    # repeated)": 5-fold stratified CV repeated 5 times with different
    # shuffles = 25 train/validate splits total. We train fresh clones
    # of the model (same hyperparameters, but unfitted) on each split's
    # training portion and evaluate on its validation portion, so the
    # mean AUC plotted here matches the number printed during training.
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=42)
    mean_fpr = np.linspace(0, 1, 100)
    interpolated_tprs = []
    fold_aucs = []

    for i, (train_idx, val_idx) in enumerate(cv.split(X_train, y_train)):
        fold_model = clone(model)
        fold_model.fit(X_train[train_idx], y_train[train_idx])

        fpr, tpr, _ = roc_curve(y_train[val_idx], fold_model.predict_proba(X_train[val_idx])[:, 1])
        fold_auc = auc(fpr, tpr)
        fold_aucs.append(fold_auc)
        interpolated_tprs.append(np.interp(mean_fpr, fpr, tpr))

        # With 25 folds, a low alpha keeps individual curves visible
        # without the plot becoming a solid blue smear.
        ax.plot(fpr, tpr, color="tab:blue", alpha=0.10, lw=1,
                label="CV folds (5x5 repeated)" if i == 0 else None)

    mean_tpr = np.mean(interpolated_tprs, axis=0)
    ax.plot(mean_fpr, mean_tpr, color="tab:blue", lw=2,
            label=f"CV mean (AUC={np.mean(fold_aucs):.3f} +/- {np.std(fold_aucs):.3f})")

    # --- Test set ROC curve using the final, already-fitted model ---
    fpr, tpr, _ = roc_curve(y_test, model.predict_proba(X_test)[:, 1])
    test_auc = auc(fpr, tpr)
    ax.plot(fpr, tpr, color="tab:red", lw=2, label=f"Held-out test (AUC={test_auc:.3f})")

    # Reference diagonal: a random classifier (AUC=0.5).
    ax.plot([0, 1], [0, 1], "k--", lw=1, label="Random (AUC=0.5)")

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("QSAR model ROC curves: CV folds vs held-out test set")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()

    out = outdir / "roc_curves.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved %s", out)


def plot_chemical_space(X_train, y_train, X_test, y_test, outdir: Path):
    """PCA scatter of the Morgan fingerprints, colored by label and
    marked by train/test split.
    """
    X_all = np.vstack([X_train, X_test])

    # Reduce 2048 binary fingerprint bits to 2 principal components
    # purely for visualization -- this is NOT how the model itself
    # makes predictions (it uses the full 2048-dim space).
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_all)

    n_train = len(X_train)
    coords_train, coords_test = coords[:n_train], coords[n_train:]

    fig, ax = plt.subplots(figsize=(7, 6))
    for label, color, name in [(1, "tab:red", "active"), (0, "tab:blue", "inactive")]:
        mask = y_train == label
        ax.scatter(coords_train[mask, 0], coords_train[mask, 1],
                   c=color, marker="o", alpha=0.6,
                   label=f"train {name} (n={mask.sum()})")

        mask = y_test == label
        ax.scatter(coords_test[mask, 0], coords_test[mask, 1],
                   c=color, marker="^", alpha=0.9, edgecolors="black", linewidths=0.5,
                   label=f"test {name} (n={mask.sum()})")

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)")
    ax.set_title("Chemical space (PCA of 2048-bit Morgan fingerprints)")
    ax.legend(fontsize=8)
    fig.tight_layout()

    out = outdir / "chemical_space_pca.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved %s", out)


def plot_probability_distributions(model, X_train, y_train, X_test, y_test,
                                     outdir: Path, threshold: float):
    """Histograms of predicted P(active) for train and test sets,
    split by true label.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    bins = np.linspace(0, 1, 21)

    for ax, X, y, title in [
        (axes[0], X_train, y_train, "Training set"),
        (axes[1], X_test, y_test, "Held-out test set"),
    ]:
        probs = model.predict_proba(X)[:, 1]
        ax.hist(probs[y == 0], bins=bins, alpha=0.6, color="tab:blue", label="inactive (y=0)")
        ax.hist(probs[y == 1], bins=bins, alpha=0.6, color="tab:red", label="active (y=1)")
        ax.axvline(threshold, color="black", linestyle="--", lw=1,
                   label=f"qsar.min_probability={threshold}")
        ax.set_xlabel("Predicted P(active)")
        ax.set_title(title)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("Number of molecules")
    fig.suptitle("Predicted-probability distributions by true class")
    fig.tight_layout()

    out = outdir / "probability_distributions.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved %s", out)


def plot_confusion_matrix(model, X_test, y_test, outdir: Path, threshold: float):
    """Confusion matrix on the held-out test set at the given
    probability threshold (default: qsar.min_probability=0.6).
    """
    probs = model.predict_proba(X_test)[:, 1]
    preds = (probs >= threshold).astype(int)
    cm = confusion_matrix(y_test, preds, labels=[0, 1])

    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    disp = ConfusionMatrixDisplay(cm, display_labels=["inactive", "active"])
    disp.plot(ax=ax, cmap="Blues", colorbar=False)
    ax.set_title(f"Test set confusion matrix (threshold={threshold})")
    fig.tight_layout()

    out = outdir / "confusion_matrix.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved %s", out)


def main(active_file: str, inactive_file: str, model_path: str, outdir: str, threshold: float):
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    X, y = load_dataset(active_file, inactive_file)

    # Recreate the EXACT train/test split used by train_qsar_mod.py
    # (same test_size, random_state, and stratify=y) so these plots
    # describe the same train/test molecules the model was evaluated on.
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    model = joblib.load(model_path)

    plot_roc_cv_and_test(model, X_train, y_train, X_test, y_test, out_path)
    plot_chemical_space(X_train, y_train, X_test, y_test, out_path)
    plot_probability_distributions(model, X_train, y_train, X_test, y_test, out_path, threshold)
    plot_confusion_matrix(model, X_test, y_test, out_path, threshold)

    log.info("All plots saved to %s", out_path)


if __name__ == "__main__":
    # Command-line entry point:
    #   python plot_qsar_results.py --actives data/input/actives.csv \
    #                                --inactives data/input/inactives.csv \
    #                                --model models/qsar_model.joblib \
    #                                --outdir results/qsar_plots
    parser = argparse.ArgumentParser()
    parser.add_argument("--actives", default="data/input/actives.csv")
    parser.add_argument("--inactives", default="data/input/inactives.csv")
    parser.add_argument("--model", default="models/qsar_model.joblib")
    parser.add_argument("--outdir", default="results/qsar_plots")
    parser.add_argument("--threshold", type=float, default=0.6,
                         help="Probability cutoff for the confusion matrix "
                              "(should match qsar.min_probability in config.yaml)")
    args = parser.parse_args()
    main(args.actives, args.inactives, args.model, args.outdir, args.threshold)
