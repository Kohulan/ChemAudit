"""
MCS (Maximum Common Substructure) Pairwise Comparison Service

Computes maximum common substructure between two molecules using rdFMCS.FindMCS
with a hardcoded timeout of 10 seconds (Phase 13 decision D-09). Returns MCS
SMARTS, atom/bond counts, Tanimoto similarity, and 8 property deltas.
"""

from __future__ import annotations

import logging
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdFMCS
from rdkit.Chem.QED import qed

logger = logging.getLogger(__name__)

# Morgan fingerprint generator (module-level singleton for reuse)
_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Property calculators: (name, callable)
_PROPERTY_CALCULATORS: list[tuple[str, Any]] = [
    ("MW", Descriptors.MolWt),
    ("LogP", Descriptors.MolLogP),
    ("TPSA", Descriptors.TPSA),
    ("QED", qed),
    ("HBA", Descriptors.NumHAcceptors),
    ("HBD", Descriptors.NumHDonors),
    ("RotBonds", Descriptors.NumRotatableBonds),
    ("RingCount", Descriptors.RingCount),
]


def compute_mcs_comparison(smiles_a: str, smiles_b: str) -> dict:
    """
    Compute MCS comparison between two molecules.

    Finds the maximum common substructure, computes Tanimoto similarity
    using Morgan fingerprints, and calculates property deltas for 8
    physicochemical properties.

    Args:
        smiles_a: SMILES string for molecule A.
        smiles_b: SMILES string for molecule B.

    Returns:
        Dict with keys:
        - ``mcs_smarts``: MCS as SMARTS string
        - ``num_atoms``: number of atoms in MCS
        - ``num_bonds``: number of bonds in MCS
        - ``timed_out``: whether FindMCS timed out
        - ``tanimoto``: Tanimoto similarity (0-1)
        - ``property_deltas``: list of 8 property delta dicts
        - ``smiles_a``: input SMILES A
        - ``smiles_b``: input SMILES B

    Raises:
        ValueError: If either SMILES cannot be parsed.
    """
    mol_a = Chem.MolFromSmiles(smiles_a)
    if mol_a is None:
        raise ValueError(f"Cannot parse SMILES A: {smiles_a}")

    mol_b = Chem.MolFromSmiles(smiles_b)
    if mol_b is None:
        raise ValueError(f"Cannot parse SMILES B: {smiles_b}")

    # Step 1: Find MCS with hardcoded timeout=10 (D-09, never configurable)
    mcs = rdFMCS.FindMCS([mol_a, mol_b], timeout=10)

    # Step 2: Compute Tanimoto similarity using Morgan fingerprints
    fp_a = _MORGAN_GEN.GetFingerprint(mol_a)
    fp_b = _MORGAN_GEN.GetFingerprint(mol_b)
    tanimoto = DataStructs.TanimotoSimilarity(fp_a, fp_b)

    # Step 3: Compute property deltas for 8 properties
    property_deltas: list[dict] = []
    for prop_name, calc in _PROPERTY_CALCULATORS:
        val_a = float(calc(mol_a))
        val_b = float(calc(mol_b))
        property_deltas.append(
            {
                "property": prop_name,
                "mol_a": val_a,
                "mol_b": val_b,
                "delta": val_a - val_b,
            }
        )

    return {
        "mcs_smarts": mcs.smartsString,
        "num_atoms": mcs.numAtoms,
        "num_bonds": mcs.numBonds,
        "timed_out": bool(mcs.canceled),
        "tanimoto": float(tanimoto),
        "property_deltas": property_deltas,
        "smiles_a": smiles_a,
        "smiles_b": smiles_b,
    }
