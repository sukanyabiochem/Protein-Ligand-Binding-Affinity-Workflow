"""Shared utilities for the screening pipeline."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterator

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, QED, rdFingerprintGenerator, rdMolDescriptors

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def load_config(path: str = "config.yaml") -> dict:
    import yaml
    with open(path) as f:
        return yaml.safe_load(f)


def standardize_smiles(smiles: str) -> str | None:
    """Canonicalize SMILES; return None if invalid."""
    if not smiles or not isinstance(smiles, str):
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


def mol_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def calc_descriptors(mol) -> dict:
    """RDKit 2D descriptors used for ADMET rules and QSAR."""
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "RotBonds": Lipinski.NumRotatableBonds(mol),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "QED": QED.qed(mol),
        "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
        "NumRings": rdMolDescriptors.CalcNumRings(mol),
    }


# Module-level cached generator -- avoids rebuilding the generator object
# on every call (5,000+ times in Stage 4) and silences the deprecation
# warning from the old AllChem.GetMorganFingerprintAsBitVect API.
_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def morgan_fp(mol, radius: int = 2, n_bits: int = 2048):
    if radius == 2 and n_bits == 2048:
        return _morgan_gen.GetFingerprint(mol)
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=n_bits
    ).GetFingerprint(mol)


def read_input_chunks(cfg: dict, chunk_size: int) -> Iterator[pd.DataFrame]:
    """Read SMILES library in chunks."""
    inp = cfg["input"]
    path = Path(inp["file"])
    fmt = inp.get("format", path.suffix.lstrip(".")).lower()
    id_col = inp.get("id_column", "id")
    smi_col = inp.get("smiles_column", "smiles")

    if fmt in ("csv", "tsv"):
        sep = "\t" if fmt == "tsv" else ","
        for chunk in pd.read_csv(path, sep=sep, chunksize=chunk_size):
            chunk = chunk.rename(columns={id_col: "id", smi_col: "smiles"})
            yield chunk[["id", "smiles"]]
    elif fmt in ("smi", "smiles"):
        rows, batch = [], []
        with open(path) as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                smi = parts[0]
                mid = parts[1] if len(parts) > 1 else f"mol_{i}"
                batch.append({"id": mid, "smiles": smi})
                if len(batch) >= chunk_size:
                    yield pd.DataFrame(batch)
                    batch = []
        if batch:
            yield pd.DataFrame(batch)
    elif fmt == "sdf":
        suppl = Chem.SDMolSupplier(str(path), removeHs=True)
        batch = []
        for i, mol in enumerate(suppl):
            if mol is None:
                continue
            smi = Chem.MolToSmiles(mol)
            mid = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i}"
            batch.append({"id": mid, "smiles": smi})
            if len(batch) >= chunk_size:
                yield pd.DataFrame(batch)
                batch = []
        if batch:
            yield pd.DataFrame(batch)
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def save_parquet(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    log.info("Saved %d rows --> %s", len(df), path)
