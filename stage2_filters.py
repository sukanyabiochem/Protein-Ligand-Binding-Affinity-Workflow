"""Stage 2-3: Drug-likeness and structural alert filters.

This stage takes the canonicalized, deduplicated molecules produced by
Stage 1 (data/intermediate/stage1/all_standardized.parquet) and removes
any molecule that:
  1. Violates Lipinski's Rule of Five / Veber's rules (drug-likeness), or
  2. Contains a PAINS (Pan-Assay INterference compoundS) or BRENK
     (unstable / toxic / reactive) structural alert.

Surviving molecules, along with their computed 2D descriptors
(MW, LogP, HBD, HBA, TPSA, etc.), are written to
data/intermediate/stage2/filtered.parquet for use by later stages
(ADMET rules, QSAR scoring, docking).
"""
from __future__ import annotations

import argparse
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import FilterCatalog, Lipinski, rdMolDescriptors

from utils import calc_descriptors, load_config, log, mol_from_smiles, save_parquet


def _build_catalog():
    """Build an RDKit FilterCatalog containing the PAINS and BRENK
    substructure libraries.

    - PAINS: substructures known to cause false positives in many types
      of biological assays (e.g. reactive groups that hit many targets
      non-specifically).
    - BRENK: substructures associated with reactivity, instability,
      or toxicity (e.g. Michael acceptors, aldehydes, etc.).

    Both catalogs are combined into a single FilterCatalog: a molecule
    is "flagged" if it matches ANY pattern from EITHER catalog.
    """
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)


# Module-level placeholder for the FilterCatalog object.
# Each multiprocessing worker process gets its own copy of this module,
# so this starts as None in every worker until _init_worker() runs.
CATALOG = None


def _init_worker():
    """Pool initializer: runs once per worker process when the Pool starts.

    Building the FilterCatalog is relatively expensive, so we build it
    ONCE per worker process here (rather than once per molecule inside
    _filter_row), and store it in the worker's global CATALOG variable
    so every call to _filter_row in that process can reuse it.
    """
    global CATALOG
    CATALOG = _build_catalog()


def _filter_row(row, fcfg):
    """Apply all Stage 2-3 checks to a single molecule.

    Parameters
    ----------
    row : dict
        A single record with at least "id" and "smiles" keys
        (one row from the Stage 1 output).
    fcfg : dict
        The "filters" section of config.yaml, containing the
        Lipinski/Veber thresholds and PAINS/BRENK toggle flags.

    Returns
    -------
    dict or None
        A dict of computed descriptors (+ id/smiles) if the molecule
        PASSES every check, otherwise None (molecule is dropped).
    """
    # Step 1: Re-parse the canonical SMILES into an RDKit Mol object.
    # (mol_from_smiles also calls Chem.SanitizeMol internally.)
    mol = mol_from_smiles(row["smiles"])
    if mol is None:
        # Should be rare since Stage 1 already validated the SMILES,
        # but guard against any edge cases.
        return None

    # Step 2: Compute the full 2D descriptor set for this molecule.
    # calc_descriptors() (in utils.py) returns a dict with keys:
    #   MW, LogP, HBD, HBA, TPSA, RotBonds, HeavyAtoms,
    #   AromaticRings, QED, FractionCSP3, NumRings
    desc = calc_descriptors(mol)

    # Carry the identifier and canonical SMILES through to the output
    # so downstream stages can trace each row back to its molecule.
    desc["id"] = row["id"]
    desc["smiles"] = row["smiles"]

    # ------------------------------------------------------------
    # Step 3: Lipinski's Rule of Five + Veber's rules.
    # Each check below is an "early exit" — the FIRST rule a molecule
    # fails causes it to be rejected immediately (return None), so a
    # molecule does not need to fail every rule to be dropped.
    # ------------------------------------------------------------

    # Molecular weight must not exceed the configured maximum
    # (Lipinski: MW <= 500 typically).
    if desc["MW"] > fcfg["max_mw"]:
        return None

    # Calculated LogP (lipophilicity) must not exceed the maximum
    # (Lipinski: LogP <= 5 typically).
    if desc["LogP"] > fcfg["max_logp"]:
        return None

    # Number of hydrogen bond donors (e.g. OH, NH groups) must not
    # exceed the maximum (Lipinski: HBD <= 5 typically).
    if desc["HBD"] > fcfg["max_hbd"]:
        return None

    # Number of hydrogen bond acceptors (e.g. O, N atoms) must not
    # exceed the maximum (Lipinski: HBA <= 10 typically).
    if desc["HBA"] > fcfg["max_hba"]:
        return None

    # Number of rotatable bonds must not exceed the maximum
    # (Veber rule: RotBonds <= 10, relates to oral bioavailability).
    if desc["RotBonds"] > fcfg["max_rotatable_bonds"]:
        return None

    # Topological Polar Surface Area (TPSA) must fall within the
    # configured [min_tpsa, max_tpsa] range
    # (Veber rule: TPSA <= 140 A^2 typically; a minimum is also
    # enforced here to filter out overly non-polar molecules).
    if desc["TPSA"] < fcfg["min_tpsa"] or desc["TPSA"] > fcfg["max_tpsa"]:
        return None

    # Total heavy (non-hydrogen) atom count must not exceed the
    # configured maximum (a general size/complexity cap).
    if desc["HeavyAtoms"] > fcfg["max_heavy_atoms"]:
        return None

    # ------------------------------------------------------------
    # Step 4: Structural alert screening (PAINS / BRENK).
    # Only run if at least one of the two toggles is enabled in config.
    # NOTE: because _build_catalog() merges PAINS and BRENK into ONE
    # catalog, this check cannot distinguish "PAINS only" vs
    # "BRENK only" — enabling either flag screens against BOTH.
    # ------------------------------------------------------------
    if fcfg.get("exclude_pains") or fcfg.get("exclude_brenk"):
        # HasMatch() returns True if the molecule contains ANY
        # substructure pattern from the PAINS/BRENK catalog.
        if CATALOG.HasMatch(mol):
            return None

    # Molecule passed every check -> return its descriptors so it can
    # be written to the Stage 2 output file.
    return desc


def run(cfg: dict, input_path: str | None = None):
    """Run Stage 2-3 over the full Stage 1 output.

    Parameters
    ----------
    cfg : dict
        The full parsed config.yaml.
    input_path : str or None
        Optional override for the directory containing
        "stage1/all_standardized.parquet". If not given, defaults to
        cfg["processing"]["output_dir"].

    Returns
    -------
    pandas.DataFrame
        The molecules (with descriptors) that passed all filters.
    """
    # Load the Lipinski/Veber thresholds + PAINS/BRENK toggles.
    fcfg = cfg["filters"]

    # Stage 2 always reads Stage 1's combined, deduplicated output.
    infile = Path(input_path or cfg["processing"]["output_dir"]) / "stage1" / "all_standardized.parquet"
    df = pd.read_parquet(infile)
    log.info("Stage 2-3: filtering %d molecules", len(df))

    n_workers = cfg["processing"]["n_workers"]

    # Convert the DataFrame to a list of plain dicts so each row can be
    # sent to a worker process individually.
    rows = df.to_dict("records")

    # Create a worker pool. `initializer=_init_worker` makes each worker
    # build its own FilterCatalog exactly once (see _init_worker above).
    with Pool(n_workers, initializer=_init_worker) as pool:
        # starmap is used (instead of map) because _filter_row needs
        # TWO arguments per call: (row, fcfg). We zip each row with the
        # same fcfg dict.
        results = pool.starmap(_filter_row, [(r, fcfg) for r in rows])

    # Drop the None results (molecules that failed any filter) and
    # build a DataFrame from the surviving descriptor dicts.
    passed = pd.DataFrame([r for r in results if r is not None])

    # Write the surviving molecules + descriptors to Stage 2's output.
    out = Path(cfg["processing"]["output_dir"]) / "stage2" / "filtered.parquet"
    save_parquet(passed, out)

    # Report how many molecules passed, out of how many were checked.
    log.info("Stage 2-3 complete: %d / %d passed (%.1f%%)",
             len(passed), len(df), 100 * len(passed) / max(len(df), 1))
    return passed


if __name__ == "__main__":
    # Command-line entry point:
    #   python stage2_filters.py --config config.yaml [--input PATH]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--input", default=None)
    args = parser.parse_args()
    run(load_config(args.config), args.input)
