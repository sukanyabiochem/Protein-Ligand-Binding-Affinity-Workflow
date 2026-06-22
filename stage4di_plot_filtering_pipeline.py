"""Visualize the virtual screening pipeline filtering funnel and the
chemical diversity of Stage 4 output molecules.

This script reads the parquet files produced by each stage and generates
four diagnostic plots:

  1. pipeline_funnel.png
     Horizontal bar chart showing how many molecules survived each stage:
       Stage 1 (standardize/dedup) -> Stage 2 (Lipinski/PAINS/BRENK)
       -> Stage 3 (ADMET) -> Stage 4 (QSAR top-500).
     Gives a visual sense of how much the library was trimmed at each step.

  2. qsar_score_distribution.png
     Histogram of QSAR predicted-probability scores for ALL Stage 3
     molecules (5,825), scored by re-running the trained model.
     The vertical dashed line marks qsar.min_probability=0.6 (the cutoff
     used in Stage 4), and the molecules that passed (score >= 0.6) are
     shown in a distinct color.
     Helps validate that the threshold captures a meaningful tail of
     the score distribution rather than cutting arbitrarily.

  3. chemical_space_pca.png
     2D PCA of Morgan fingerprints for all Stage 3 molecules (5,825,
     shown in light grey) with the Stage 4 selected molecules (500,
     shown colored by qsar_score) overlaid.
     Shows WHICH REGION of chemical space Stage 4 selected -- if the
     colored dots cluster tightly, the QSAR model found a specific
     structural motif; if they spread out, it selected diverse hits.

  4. tanimoto_dendrogram.png
     Hierarchical clustering dendrogram of the 500 Stage 4 molecules,
     using Tanimoto distance computed from Morgan fingerprints.
     Tanimoto distance = 1 - Tanimoto_similarity (ranges 0-1; 0 = same
     molecule, 1 = no shared fingerprint bits).
     The UPGMA (average-linkage) method is used because it is compatible
     with any distance metric (unlike Ward linkage, which requires
     Euclidean distance).
     Reveals clusters of structurally similar molecules within the 500
     shortlist -- useful for identifying chemical series and ensuring
     structural diversity before docking.

USAGE:
    python plot_filtering_pipeline.py --config config.yaml --outdir results/filtering_plots

REQUIREMENTS:
    matplotlib  -- conda install -c conda-forge matplotlib
    scipy       -- conda install -c conda-forge scipy  (likely already present)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")   # non-interactive backend -- saves PNG without opening a window
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from rdkit import DataStructs
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA

from utils import load_config, log, mol_from_smiles, morgan_fp


# ---------------------------------------------------------------------------
# Helper: build fingerprint matrix from a list of SMILES strings
# ---------------------------------------------------------------------------

def smiles_to_fp_matrix(smiles_list: list[str]):
    """Convert SMILES to a (n_valid, 2048) numpy matrix of Morgan fingerprints
    and return the matrix plus a boolean mask of which rows are valid.

    Molecules that fail to parse (mol_from_smiles returns None) are
    silently skipped; the mask lets the caller align results back to the
    original list if needed.
    """
    fps_rdkit = []   # RDKit ExplicitBitVect objects (needed for Tanimoto)
    fps_numpy = []   # numpy 0/1 arrays (needed for sklearn PCA / model)
    valid_mask = []  # True for rows that were successfully computed

    for smi in smiles_list:
        mol = mol_from_smiles(smi)
        if mol is None:
            valid_mask.append(False)
            continue
        fp = morgan_fp(mol)           # RDKit ExplicitBitVect, 2048 bits
        arr = np.zeros(2048, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fps_rdkit.append(fp)
        fps_numpy.append(arr)
        valid_mask.append(True)

    X = np.vstack(fps_numpy) if fps_numpy else np.empty((0, 2048), dtype=np.uint8)
    return X, fps_rdkit, valid_mask


# ---------------------------------------------------------------------------
# Plot 1: Pipeline filtering funnel
# ---------------------------------------------------------------------------

def plot_pipeline_funnel(stage_counts: dict, outdir: Path):
    """Horizontal bar chart of molecule counts at each pipeline stage.

    Parameters
    ----------
    stage_counts : dict
        Ordered dict: {stage_label: count}, e.g.
        {"Stage 1\\nStandardize": 9999, "Stage 2\\nLipinski+PAINS": 6853, ...}
    outdir : Path
        Directory to save the PNG.
    """
    labels = list(stage_counts.keys())
    counts = list(stage_counts.values())
    max_count = max(counts)

    # Color bars from light to dark blue to represent progressive filtering.
    # Darker = fewer molecules = later / more selective stage.
    cmap = plt.cm.Blues
    colors = [cmap(0.35 + 0.55 * (i / (len(labels) - 1))) for i in range(len(labels))]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    bars = ax.barh(labels, counts, color=colors, edgecolor="white", height=0.55)

    # Annotate each bar with its count and the percentage remaining from Stage 1.
    start = counts[0]
    for bar, count in zip(bars, counts):
        pct = count / start * 100
        ax.text(
            bar.get_width() + max_count * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{count:,}  ({pct:.1f}%)",
            va="center", ha="left", fontsize=9,
        )

    ax.set_xlabel("Number of molecules")
    ax.set_xlim(0, max_count * 1.25)
    ax.set_title("Virtual screening pipeline: molecules surviving each stage")
    ax.invert_yaxis()   # Stage 1 at top, last stage at bottom
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    out = outdir / "pipeline_funnel.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Plot 2: QSAR score distribution for all Stage 3 molecules
# ---------------------------------------------------------------------------

def plot_qsar_score_distribution(scores_all: np.ndarray, threshold: float, outdir: Path):
    """Histogram of QSAR predicted probabilities for all Stage 3 molecules,
    split into 'failed' (score < threshold) and 'passed' (>= threshold).

    Parameters
    ----------
    scores_all : np.ndarray
        Predicted P(active) for every Stage 3 molecule (shape: n_stage3,).
    threshold : float
        qsar.min_probability from config.yaml (e.g. 0.6).
    outdir : Path
    """
    bins = np.linspace(0, 1, 41)   # 40 bins across the full 0-1 range

    passed  = scores_all[scores_all >= threshold]
    failed  = scores_all[scores_all <  threshold]

    fig, ax = plt.subplots(figsize=(8, 4.5))

    ax.hist(failed,  bins=bins, color="#aec6e8", edgecolor="white", label=f"Below threshold (n={len(failed):,})")
    ax.hist(passed,  bins=bins, color="#1f6eb5", edgecolor="white", label=f"Above threshold -- passed Stage 4 (n={len(passed):,})")
    ax.axvline(threshold, color="crimson", linestyle="--", lw=1.5,
               label=f"qsar.min_probability = {threshold}")

    ax.set_xlabel("QSAR predicted P(active)")
    ax.set_ylabel("Number of molecules")
    ax.set_title("QSAR score distribution: all Stage 3 molecules (n={:,})".format(len(scores_all)))
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    out = outdir / "qsar_score_distribution.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Plot 3: Chemical space PCA (Stage 3 vs Stage 4)
# ---------------------------------------------------------------------------

def plot_chemical_space_pca(X_all: np.ndarray, selected_mask: np.ndarray,
                             scores_selected: np.ndarray, outdir: Path):
    """PCA scatter of all Stage 3 molecules, highlighting Stage 4 selections.

    Parameters
    ----------
    X_all : np.ndarray, shape (n_stage3, 2048)
        Morgan fingerprint matrix for all Stage 3 molecules.
    selected_mask : np.ndarray of bool, shape (n_stage3,)
        True for molecules that are in Stage 4's output.
    scores_selected : np.ndarray
        QSAR scores for the selected molecules (used for color mapping).
    outdir : Path
    """
    # Fit PCA on all Stage 3 molecules so the "map" represents the full
    # chemical space, then project both the background and selected onto it.
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_all)

    fig, ax = plt.subplots(figsize=(8, 7))

    # Background: all Stage 3 molecules that were NOT selected (grey)
    not_selected = ~selected_mask
    ax.scatter(
        coords[not_selected, 0], coords[not_selected, 1],
        c="#cccccc", s=5, alpha=0.4, linewidths=0,
        label=f"Stage 3 (not selected, n={not_selected.sum():,})",
        rasterized=True,   # rasterize dense grey scatter for smaller file size
    )

    # Foreground: Stage 4 selected molecules, colored by qsar_score
    sc = ax.scatter(
        coords[selected_mask, 0], coords[selected_mask, 1],
        c=scores_selected, cmap="plasma", s=18, alpha=0.9,
        linewidths=0.3, edgecolors="black",
        label=f"Stage 4 selected (n={selected_mask.sum():,})",
        vmin=0.6, vmax=1.0,   # color scale starts at the 0.6 threshold
    )

    cbar = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("QSAR score (P active)", fontsize=9)

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)")
    ax.set_title("Chemical space (PCA of Morgan fingerprints): Stage 3 vs Stage 4")
    ax.legend(fontsize=8, markerscale=2)
    fig.tight_layout()

    out = outdir / "chemical_space_pca.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Plot 4: Tanimoto dendrogram of Stage 4 molecules
# ---------------------------------------------------------------------------

def plot_tanimoto_dendrogram(fps_rdkit: list, ids: list[str], outdir: Path):
    """Hierarchical clustering dendrogram of Stage 4 molecules using
    Tanimoto distance from Morgan fingerprints.

    Tanimoto similarity between two bit-vectors A and B:
        Tanimoto(A, B) = |A AND B| / |A OR B|
    (number of shared ON bits / number of bits ON in either vector)

    Tanimoto DISTANCE = 1 - Tanimoto_similarity (0 = identical, 1 = no overlap)

    UPGMA (average linkage) merges clusters by the average distance
    between all pairs across the two clusters. It is suitable for any
    distance metric (unlike Ward linkage which requires Euclidean).

    Parameters
    ----------
    fps_rdkit : list of RDKit ExplicitBitVect
        Morgan fingerprints for the Stage 4 molecules.
    ids : list[str]
        Molecule IDs -- used as labels only if n_molecules is small.
    outdir : Path
    """
    n = len(fps_rdkit)
    log.info("Computing %d x %d Tanimoto distance matrix ...", n, n)

    # Build the full n x n Tanimoto SIMILARITY matrix using RDKit's
    # BulkTanimotoSimilarity, which is heavily optimized in C++ and
    # much faster than a Python nested loop.
    sim_matrix = np.zeros((n, n), dtype=np.float32)
    for i, fp in enumerate(fps_rdkit):
        # BulkTanimotoSimilarity returns a list of similarities between
        # fp and every fingerprint in fps_rdkit simultaneously.
        sim_matrix[i] = DataStructs.BulkTanimotoSimilarity(fp, fps_rdkit)

    # Convert similarity to distance: distance = 1 - similarity.
    # The diagonal (self-similarity = 1) becomes 0 (distance to itself = 0).
    dist_matrix = (1.0 - sim_matrix).astype(np.float64)

    # squareform() converts a square symmetric matrix to the condensed
    # upper-triangle form that scipy's linkage() expects.
    condensed = squareform(dist_matrix, checks=False)

    # Hierarchical clustering with UPGMA (average linkage).
    # Z is an (n-1) x 4 array encoding the merge history:
    # each row = [cluster_i, cluster_j, merge_distance, cluster_size]
    log.info("Running hierarchical clustering (UPGMA) ...")
    Z = linkage(condensed, method="average")

    # Plot the dendrogram.
    # With 500 molecules all leaves would overlap; truncate_mode='lastp'
    # shows only the last p=40 merges (top of the tree = most distinct
    # clusters), which gives a readable summary of major structural groups.
    fig, ax = plt.subplots(figsize=(14, 5))
    dendrogram(
        Z,
        ax=ax,
        truncate_mode="lastp",   # show only the top p clusters
        p=40,                    # number of terminal nodes to display
        show_leaf_counts=True,   # show (n) = number of molecules per cluster
        leaf_rotation=90,
        leaf_font_size=7,
        color_threshold=0.6,     # color branches whose merge distance > 0.6
        above_threshold_color="grey",
    )

    ax.set_xlabel("Cluster (number of molecules in parentheses)")
    ax.set_ylabel("Tanimoto distance (1 - similarity)")
    ax.set_title(
        f"Hierarchical clustering of Stage 4 molecules (n={n}) "
        f"by Morgan fingerprint Tanimoto distance\n"
        f"[truncated to top 40 clusters; branch colors: Tanimoto distance < 0.6]"
    )
    ax.axhline(0.6, color="crimson", linestyle="--", lw=1,
               label="Tanimoto distance = 0.6 (Tc similarity = 0.4, rough diversity threshold)")
    ax.legend(fontsize=7)
    fig.tight_layout()

    out = outdir / "tanimoto_dendrogram.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def main(config_path: str, outdir: str):
    cfg = load_config(config_path)
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    inter = Path(cfg["processing"]["output_dir"])

    # ------------------------------------------------------------------
    # 1. Read all stage parquet files and count rows dynamically
    # ------------------------------------------------------------------
    log.info("Reading stage parquet files ...")
    df_s1 = pd.read_parquet(inter / "stage1" / "all_standardized.parquet")
    df_s2 = pd.read_parquet(inter / "stage2" / "filtered.parquet")
    df_s3 = pd.read_parquet(inter / "stage3" / "admet_passed.parquet")
    df_s4 = pd.read_parquet(inter / "stage4" / "qsar_ranked.parquet")

    stage_counts = {
        "Stage 1\nStandardize & dedup":   len(df_s1),
        "Stage 2\nLipinski + PAINS/BRENK": len(df_s2),
        "Stage 3\nADMET filters":          len(df_s3),
        "Stage 4\nQSAR scoring (top-N)":   len(df_s4),
    }
    for label, count in stage_counts.items():
        log.info("  %s: %d molecules", label.replace("\n", " "), count)

    # ------------------------------------------------------------------
    # 2. Plot 1: pipeline funnel
    # ------------------------------------------------------------------
    plot_pipeline_funnel(stage_counts, out_path)

    # ------------------------------------------------------------------
    # 3. Score ALL Stage 3 molecules with the QSAR model so we can show
    #    the full score distribution (Stage 4 only saves the top 500).
    # ------------------------------------------------------------------
    log.info("Loading QSAR model and scoring all %d Stage 3 molecules ...", len(df_s3))
    model = joblib.load(cfg["qsar"]["model_path"])
    threshold = cfg["qsar"]["min_probability"]

    X_s3, fps_s3_rdkit, mask_s3 = smiles_to_fp_matrix(df_s3["smiles"].tolist())
    scores_s3 = model.predict_proba(X_s3)[:, 1]   # P(active) for every valid molecule

    # ------------------------------------------------------------------
    # 4. Plot 2: QSAR score distribution
    # ------------------------------------------------------------------
    plot_qsar_score_distribution(scores_s3, threshold, out_path)

    # ------------------------------------------------------------------
    # 5. Build a boolean mask: which Stage 3 molecules ended up in Stage 4?
    #    Stage 4 IDs are the reference; we mark corresponding Stage 3 rows.
    # ------------------------------------------------------------------
    s4_ids = set(df_s4["id"].astype(str))
    # df_s3 (valid rows only, aligned to X_s3/scores_s3 via mask_s3)
    df_s3_valid = df_s3[mask_s3].reset_index(drop=True)
    selected_mask = np.array([str(mid) in s4_ids for mid in df_s3_valid["id"]])

    # QSAR scores for the selected molecules (for PCA color mapping)
    id_to_score = dict(zip(df_s4["id"].astype(str), df_s4["qsar_score"]))
    scores_selected = np.array([
        id_to_score.get(str(mid), 0.6)
        for mid in df_s3_valid.loc[selected_mask, "id"]
    ])

    # ------------------------------------------------------------------
    # 6. Plot 3: chemical space PCA
    # ------------------------------------------------------------------
    plot_chemical_space_pca(X_s3, selected_mask, scores_selected, out_path)

    # ------------------------------------------------------------------
    # 7. Compute Morgan fingerprints for Stage 4 molecules only
    #    (needed for Tanimoto dendrogram -- we use rdkit ExplicitBitVect)
    # ------------------------------------------------------------------
    log.info("Computing fingerprints for Stage 4 dendrogram ...")
    _, fps_s4_rdkit, mask_s4 = smiles_to_fp_matrix(df_s4["smiles"].tolist())
    ids_s4 = df_s4["id"].astype(str).tolist()

    # ------------------------------------------------------------------
    # 8. Plot 4: Tanimoto dendrogram of Stage 4 molecules
    # ------------------------------------------------------------------
    plot_tanimoto_dendrogram(fps_s4_rdkit, ids_s4, out_path)

    log.info("All plots saved to %s", out_path)
    log.info(
        "Summary: %d -> %d -> %d -> %d molecules (%.1f%% of original library)",
        len(df_s1), len(df_s2), len(df_s3), len(df_s4),
        len(df_s4) / len(df_s1) * 100,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot the virtual screening filtering funnel and Stage 4 chemical diversity."
    )
    parser.add_argument("--config",  default="config.yaml")
    parser.add_argument("--outdir",  default="results/filtering_plots")
    args = parser.parse_args()
    main(args.config, args.outdir)
