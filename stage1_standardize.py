"""Stage 1: Canonicalize SMILES and remove duplicates."""
from __future__ import annotations

import argparse
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

from utils import load_config, log, read_input_chunks, save_parquet, standardize_smiles


def _std_row(row):
    canon = standardize_smiles(row["smiles"])
    if canon is None:
        return None
    return {"id": row["id"], "smiles": canon, "original_id": row["id"]}


def run(cfg: dict):
    out_dir = Path(cfg["processing"]["output_dir"]) / "stage1"
    out_dir.mkdir(parents=True, exist_ok=True)
    chunk_size = cfg["processing"]["chunk_size"]
    n_workers = cfg["processing"]["n_workers"]

    all_parts = []
    for i, chunk in enumerate(read_input_chunks(cfg, chunk_size)):
        log.info("Stage 1 chunk %d: %d molecules", i, len(chunk))
        with Pool(n_workers) as pool:
            results = pool.map(_std_row, chunk.to_dict("records"))
        valid = [r for r in results if r is not None]
        df = pd.DataFrame(valid)
        df = df.drop_duplicates(subset=["smiles"], keep="first")
        path = out_dir / f"chunk_{i:04d}.parquet"
        save_parquet(df, path)
        all_parts.append(df)

    combined = pd.concat(all_parts, ignore_index=True)
    combined = combined.drop_duplicates(subset=["smiles"], keep="first")
    save_parquet(combined, out_dir / "all_standardized.parquet")
    log.info("Stage 1 complete: %d unique valid molecules", len(combined))
    return combined


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    run(load_config(args.config))
