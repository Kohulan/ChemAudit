"""
One-time script to generate Morgan fingerprint matrix from ChEMBL approved drugs CSV.

Reads backend/data/chembl_approved_drugs.csv, computes Morgan fingerprints
(radius=2, 2048 bits) using GetMorganGenerator, converts to NumPy arrays,
and saves the matrix to backend/data/chembl_approved_fps.npz.

Usage:
    python backend/scripts/generate_chembl_fps.py

Output:
    backend/data/chembl_approved_fps.npz — compressed NumPy archive with key "fps"
    containing a 2D uint8 array of shape (n_valid_mols, 2048).
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Resolve paths relative to repo root (script can be run from anywhere)
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent  # backend/
_CSV_PATH = _REPO_ROOT / "data" / "chembl_approved_drugs.csv"
_NPZ_PATH = _REPO_ROOT / "data" / "chembl_approved_fps.npz"

# Morgan generator: radius=2, 2048 bits
_GEN = GetMorganGenerator(radius=2, fpSize=2048)


def generate_fps(csv_path: Path, npz_path: Path) -> int:
    """
    Generate fingerprint matrix from SMILES CSV and save as .npz.

    Args:
        csv_path: Path to input CSV with 'smiles' column.
        npz_path: Path to output .npz file.

    Returns:
        Number of valid molecules processed.
    """
    if not csv_path.exists():
        logger.error("CSV not found: %s", csv_path)
        sys.exit(1)

    smiles_list: list[str] = []
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            smiles_list.append(row["smiles"].strip())

    logger.info("Loaded %d SMILES from %s", len(smiles_list), csv_path)

    fp_arrays: list[np.ndarray] = []
    skipped = 0

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None or mol.GetNumAtoms() == 0:
            logger.warning("Skipping invalid SMILES: %s", smi)
            skipped += 1
            continue

        fp = _GEN.GetFingerprint(mol)
        arr = np.zeros(2048, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fp_arrays.append(arr)

    if not fp_arrays:
        logger.error("No valid molecules found — aborting.")
        sys.exit(1)

    matrix = np.array(fp_arrays, dtype=np.uint8)
    np.savez_compressed(npz_path, fps=matrix)
    logger.info(
        "Saved %d fingerprints (%d skipped) to %s",
        len(fp_arrays),
        skipped,
        npz_path,
    )
    return len(fp_arrays)


if __name__ == "__main__":
    generate_fps(_CSV_PATH, _NPZ_PATH)
