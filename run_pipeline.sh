#!/bin/bash
# run_pipeline.sh -- Run the entire virtual screening pipeline with one command.
#
# Usage:
#   bash run_pipeline.sh                    # uses config.yaml (default)
#   bash run_pipeline.sh my_config.yaml     # uses a custom config file

set -e  # stop immediately if any stage fails

CONFIG="${1:-config.yaml}"

echo "============================================"
echo "  Virtual Screening Pipeline"
echo "  Config: $CONFIG"
echo "============================================"

echo ""
echo "[Stage 1] Standardize & deduplicate SMILES"
python stage1_standardize.py --config "$CONFIG"

echo ""
echo "[Stage 2] Lipinski + PAINS/BRENK filters"
python stage2_filters.py --config "$CONFIG"

echo ""
echo "[Stage 3] ADMET property filters"
python stage3_admet.py --config "$CONFIG"

echo ""
echo "[Stage 4a] Bootstrap sample docking"
python stage4a_bootstrap_sample_dock.py --config "$CONFIG"

echo ""
echo "[Stage 4b] Label actives / inactives"
python stage4b_label_qsar_set.py --config "$CONFIG"

echo ""
echo "[Stage 4c] Train QSAR model"
python stage4c_train_qsar.py \
    --actives data/input/actives.csv \
    --inactives data/input/inactives.csv \
    --output models/qsar_model.joblib

echo ""
echo "[Stage 4ci] Generate QSAR diagnostic plots"
python stage4ci_plot_qsar_results.py \
    --actives data/input/actives.csv \
    --inactives data/input/inactives.csv \
    --model models/qsar_model.joblib \
    --outdir results/qsar_plots

echo ""
echo "[Stage 4d] QSAR scoring of full library"
python stage4d_qsar.py --config "$CONFIG"

echo ""
echo "[Stage 4di] Generate pipeline filtering plots"
python stage4di_plot_filtering_pipeline.py \
    --config "$CONFIG" \
    --outdir results/filtering_plots

echo ""
echo "[Stage 5] AutoDock Vina docking"
python stage5_docking.py --config "$CONFIG"

echo ""
echo "============================================"
echo "  Pipeline complete!"
echo "  Final results: data/output/final_shortlist.csv"
echo "  QSAR plots:    results/qsar_plots/"
echo "  Filter plots:  results/filtering_plots/"
echo "============================================"
