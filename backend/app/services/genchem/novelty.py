"""
Novelty checker for generative chemistry molecules.

Computes Tanimoto similarity between a query molecule and a reference set of
ChEMBL approved drug fingerprints. A molecule is considered novel if its
maximum Tanimoto similarity to the reference set is below the threshold.

The ChEMBL fingerprint matrix is loaded lazily from a pre-generated .npz file
(produced by backend/scripts/generate_chembl_fps.py). The file must exist at
backend/data/chembl_approved_fps.npz.

Reference set: Morgan fingerprints (radius=2, 2048 bits) of ~50 representative
approved drugs. Can be expanded to the full ~2,500 ChEMBL approved set by
re-running generate_chembl_fps.py with a larger input CSV.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

logger = logging.getLogger(__name__)

# Morgan generator singleton: radius=2, 2048 bits (consistent with generate_chembl_fps.py)
_GEN = GetMorganGenerator(radius=2, fpSize=2048)

# Lazy-loaded ChEMBL fingerprint reference set
_CHEMBL_FPS: Optional[list[DataStructs.ExplicitBitVect]] = None
_CHEMBL_LOADED: bool = False

# Path to the pre-generated .npz fingerprint file
_CHEMBL_NPZ_PATH = (
    Path(__file__).parent.parent.parent.parent / "data" / "chembl_approved_fps.npz"
)


def _numpy_row_to_bitvect(arr: np.ndarray) -> DataStructs.ExplicitBitVect:
    """
    Convert a 1D numpy uint8 binary array to an RDKit ExplicitBitVect.

    Args:
        arr: 1D uint8 array of length fpSize.

    Returns:
        RDKit ExplicitBitVect with bits set where arr[i] == 1.
    """
    bv = DataStructs.ExplicitBitVect(len(arr))
    for i, val in enumerate(arr):
        if val:
            bv.SetBit(i)
    return bv


def _load_chembl_fps() -> list[DataStructs.ExplicitBitVect]:
    """
    Load the ChEMBL approved drug fingerprints from the pre-generated .npz file.

    Returns a cached list of ExplicitBitVect objects on subsequent calls.
    Returns an empty list (with a warning) if the .npz file is not found.

    Returns:
        List of RDKit ExplicitBitVect objects, one per reference molecule.
    """
    global _CHEMBL_FPS, _CHEMBL_LOADED

    if _CHEMBL_LOADED:
        return _CHEMBL_FPS or []

    _CHEMBL_LOADED = True

    if not _CHEMBL_NPZ_PATH.exists():
        logger.warning(
            "ChEMBL fingerprint file not found: %s — novelty check will treat all molecules "
            "as novel. Run backend/scripts/generate_chembl_fps.py to generate it.",
            _CHEMBL_NPZ_PATH,
        )
        _CHEMBL_FPS = []
        return []

    try:
        data = np.load(_CHEMBL_NPZ_PATH)
        matrix = data["fps"]  # shape: (n_mols, fpSize)
        _CHEMBL_FPS = [_numpy_row_to_bitvect(row) for row in matrix]
        logger.info(
            "Loaded %d ChEMBL reference fingerprints from %s",
            len(_CHEMBL_FPS),
            _CHEMBL_NPZ_PATH,
        )
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to load ChEMBL fps: %s — falling back to empty reference.", exc)
        _CHEMBL_FPS = []

    return _CHEMBL_FPS or []


def check_novelty(
    mol: Chem.Mol,
    threshold: float = 0.85,
    user_fps: Optional[list[DataStructs.ExplicitBitVect]] = None,
) -> tuple[bool, Optional[float]]:
    """
    Check whether a molecule is novel relative to the ChEMBL approved drug set.

    A molecule is considered novel if its maximum Tanimoto similarity to any
    reference fingerprint is at or below the given threshold.

    Args:
        mol: RDKit molecule to check.
        threshold: Maximum allowed Tanimoto similarity (default 0.85).
            Molecules with max_sim > threshold are classified as non-novel.
        user_fps: Optional custom reference fingerprint list. If provided,
            overrides the built-in ChEMBL reference set (useful for testing).

    Returns:
        Tuple of (is_novel: bool, max_similarity: float | None).
        max_similarity is None if the reference set is empty.
    """
    query_fp = _GEN.GetFingerprint(mol)

    ref_fps = user_fps if user_fps is not None else _load_chembl_fps()

    if not ref_fps:
        return (True, None)

    sims = DataStructs.BulkTanimotoSimilarity(query_fp, ref_fps)
    max_sim = max(sims)

    return (max_sim <= threshold, max_sim)
