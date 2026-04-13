"""
Registration Hash Computation Service

Computes RDKit RegistrationHash for compound uniqueness with
``enable_tautomer_hash_v2=True`` (Phase 13 decision D-15).
Stores RDKit version in every result. Groups molecules by hash
to detect collisions (same hash = same compound or tautomers).
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import RegistrationHash, rdMolDescriptors

logger = logging.getLogger(__name__)


def compute_registration_hashes(results: list[dict[str, Any]]) -> dict:
    """
    Compute registration hashes for a batch of molecules.

    For each result with a valid SMILES, computes the RegistrationHash
    using ``enable_tautomer_hash_v2=True``. Groups molecules by hash
    to detect collisions.

    Args:
        results: Batch result dicts, each containing at least ``smiles``
            and ``index`` keys.

    Returns:
        Dict with keys:
        - ``per_molecule``: list of dicts with ``index``, ``smiles``,
          ``hash``, ``canonical_smiles``, ``formula``
        - ``unique_count``: number of unique hashes
        - ``total_count``: total processed molecules
        - ``collision_groups``: list of dicts for hash groups with count > 1
        - ``rdkit_version``: version string from rdBase.rdkitVersion
        - ``tautomer_hash_v2``: always True
    """
    per_molecule: list[dict] = []
    hash_to_indices: dict[str, list[int]] = defaultdict(list)

    for result in results:
        smiles = result.get("smiles", "")
        index = result.get("index", 0)

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.debug(
                "Registration hash: skipping invalid SMILES at index %d: %s",
                index,
                smiles,
            )
            continue

        try:
            # Compute registration hash with tautomer hash v2 enabled (D-15)
            layers = RegistrationHash.GetMolLayers(mol, enable_tautomer_hash_v2=True)
            mol_hash = RegistrationHash.GetMolHash(layers)

            canonical_smiles = Chem.MolToSmiles(mol)
            formula = rdMolDescriptors.CalcMolFormula(mol)

            per_molecule.append(
                {
                    "index": index,
                    "smiles": smiles,
                    "hash": mol_hash,
                    "canonical_smiles": canonical_smiles,
                    "formula": formula,
                }
            )

            hash_to_indices[mol_hash].append(index)

        except Exception:
            logger.warning(
                "Registration hash computation failed for index %d: %s",
                index,
                smiles,
                exc_info=True,
            )

    # Build collision groups (hashes shared by >1 molecule)
    collision_groups: list[dict] = []
    for mol_hash, indices in hash_to_indices.items():
        if len(indices) > 1:
            collision_groups.append(
                {
                    "hash": mol_hash,
                    "molecule_indices": sorted(indices),
                    "count": len(indices),
                }
            )

    unique_hashes = len(hash_to_indices)

    return {
        "per_molecule": per_molecule,
        "unique_count": unique_hashes,
        "total_count": len(per_molecule),
        "collision_groups": collision_groups,
        "rdkit_version": rdBase.rdkitVersion,
        "tautomer_hash_v2": True,
    }
