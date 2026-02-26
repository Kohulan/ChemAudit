"""
Scaffold Analysis Service

Computes Murcko scaffold decomposition, generic scaffold grouping, Shannon
entropy diversity metrics, and R-group decomposition for batch results.

Designed for use as a user-triggered expensive analytics computation via the
``run_expensive_analytics`` Celery task.
"""

import logging
import math
from typing import Any

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _extract_scaffold_smiles(mol: Chem.Mol) -> tuple[str, str]:
    """
    Extract Murcko scaffold and generic scaffold SMILES from a molecule.

    Uses the double-GetScaffoldForMol pattern to correctly remove exocyclic
    substituents that MakeScaffoldGeneric converts to single-bonded atoms.

    See: https://github.com/rdkit/rdkit/discussions/6844

    Args:
        mol: Valid RDKit molecule.

    Returns:
        Tuple of (scaffold_smiles, generic_scaffold_smiles).
        Both are empty strings if the molecule has no ring system.
    """
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return "", ""

    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold is None:
            return "", ""
        scaffold_smiles = Chem.MolToSmiles(scaffold)
    except Exception:
        logger.debug("Murcko scaffold extraction failed for a molecule — skipping.")
        return "", ""

    try:
        generic_scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        # IMPORTANT: second GetScaffoldForMol call removes exocyclic substituents
        # that MakeScaffoldGeneric converted from double- to single-bonded atoms.
        generic_scaffold = MurckoScaffold.GetScaffoldForMol(generic_scaffold)
        generic_scaffold_smiles = Chem.MolToSmiles(generic_scaffold)
    except Exception:
        logger.debug("Generic scaffold extraction failed — falling back to Murcko scaffold.")
        generic_scaffold_smiles = scaffold_smiles

    return scaffold_smiles, generic_scaffold_smiles


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_scaffold_analysis(results: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Compute Murcko scaffold diversity analysis for a batch of validation results.

    Steps:
    1. Parse each successful result's SMILES into an RDKit molecule.
    2. Extract Murcko and generic scaffolds using the double-call pattern.
    3. Group molecules by scaffold SMILES (acyclic molecules share the "" group).
    4. Compute Shannon entropy diversity metric from scaffold frequency distribution.
    5. Build a capped frequency distribution (top 50 scaffolds + "Other" bucket).

    Args:
        results: List of batch validation result dicts. Each dict must have at
            least ``"smiles"`` and optionally ``"status"`` fields. Only results
            with status ``"valid"`` or a truthy ``"mol"`` are processed.

    Returns:
        Dict matching the ``ScaffoldResult`` Pydantic schema:
        {
            "scaffolds": [{"scaffold_smiles", "generic_scaffold_smiles",
                           "molecule_indices", "count"}, ...],
            "unique_scaffold_count": int,
            "shannon_entropy": float,
            "frequency_distribution": {scaffold_smiles: count, ...},
        }
    """
    if not results:
        return {
            "scaffolds": [],
            "unique_scaffold_count": 0,
            "shannon_entropy": 0.0,
            "frequency_distribution": {},
        }

    # scaffold_smiles -> {generic, indices}
    scaffold_map: dict[str, dict[str, Any]] = {}

    for idx, result in enumerate(results):
        smiles = result.get("smiles") or result.get("original_smiles") or ""
        if not smiles:
            logger.debug("compute_scaffold_analysis: result %d missing SMILES — skipping.", idx)
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.debug(
                "compute_scaffold_analysis: result %d SMILES unparseable — skipping.", idx
            )
            continue

        scaffold_smi, generic_smi = _extract_scaffold_smiles(mol)

        if scaffold_smi not in scaffold_map:
            scaffold_map[scaffold_smi] = {
                "generic_scaffold_smiles": generic_smi,
                "molecule_indices": [],
            }
        scaffold_map[scaffold_smi]["molecule_indices"].append(idx)

    # Build sorted scaffold list (descending by count)
    scaffolds = []
    for smi, data in scaffold_map.items():
        scaffolds.append(
            {
                "scaffold_smiles": smi,
                "generic_scaffold_smiles": data["generic_scaffold_smiles"],
                "molecule_indices": data["molecule_indices"],
                "count": len(data["molecule_indices"]),
            }
        )
    scaffolds.sort(key=lambda g: g["count"], reverse=True)

    unique_scaffold_count = len(scaffolds)

    # Shannon entropy over scaffold frequency distribution
    total_n = sum(g["count"] for g in scaffolds)
    if total_n == 0 or unique_scaffold_count <= 1:
        shannon_entropy = 0.0
    else:
        shannon_entropy = 0.0
        for group in scaffolds:
            p_i = group["count"] / total_n
            if p_i > 0:
                shannon_entropy -= p_i * math.log2(p_i)

    # Frequency distribution — top 50 scaffolds, remainder in "Other"
    MAX_DISPLAY = 50
    frequency_distribution: dict[str, int] = {}
    other_count = 0

    for i, group in enumerate(scaffolds):
        if i < MAX_DISPLAY:
            frequency_distribution[group["scaffold_smiles"]] = group["count"]
        else:
            other_count += group["count"]

    if other_count > 0:
        frequency_distribution["Other"] = other_count

    return {
        "scaffolds": scaffolds,
        "unique_scaffold_count": unique_scaffold_count,
        "shannon_entropy": round(shannon_entropy, 6),
        "frequency_distribution": frequency_distribution,
    }


def compute_rgroup_decomposition(
    results: list[dict[str, Any]],
    core_smarts: str | None = None,
) -> dict[str, Any]:
    """
    Perform R-group decomposition of batch molecules around a user-specified core.

    Steps:
    1. Validate the core SMARTS string; return an error dict if invalid.
    2. Parse all successful molecules from the result list.
    3. Run ``RGroupDecompose`` from ``rdkit.Chem.rdRGroupDecomposition``.
    4. Build per-molecule decomposition dicts with R-group SMILES.
    5. Count unmatched molecules.

    Args:
        results: List of batch validation result dicts (same format as
            ``compute_scaffold_analysis``).
        core_smarts: SMARTS pattern defining the core scaffold. Required;
            if None or empty, all molecules are reported as unmatched.

    Returns:
        Dict matching the ``RGroupResult`` Pydantic schema:
        {
            "core_smarts": str,
            "decomposition": [{"molecule_index", "core", "rgroups"}, ...],
            "unmatched_count": int,
        }
        May contain an ``"error"`` key when SMARTS is invalid or decomposition
        fails catastrophically.
    """
    total = len(results)

    if not core_smarts:
        return {
            "core_smarts": core_smarts or "",
            "decomposition": [],
            "unmatched_count": total,
            "error": "No core SMARTS provided",
        }

    # Validate SMARTS
    core = Chem.MolFromSmarts(core_smarts)
    if core is None:
        return {
            "core_smarts": core_smarts,
            "decomposition": [],
            "unmatched_count": total,
            "error": "Invalid SMARTS pattern",
        }

    # Parse molecules and track original indices
    mols: list[Chem.Mol] = []
    mol_indices: list[int] = []

    for idx, result in enumerate(results):
        smiles = result.get("smiles") or result.get("original_smiles") or ""
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        mols.append(mol)
        mol_indices.append(idx)

    if not mols:
        return {
            "core_smarts": core_smarts,
            "decomposition": [],
            "unmatched_count": total,
        }

    try:
        from rdkit.Chem import rdRGroupDecomposition

        matched_rows, unmatched_indices = rdRGroupDecomposition.RGroupDecompose(
            [core], mols, asSmiles=True
        )
    except Exception as exc:
        logger.exception("compute_rgroup_decomposition: RGroupDecompose raised an exception.")
        return {
            "core_smarts": core_smarts,
            "decomposition": [],
            "unmatched_count": total,
            "error": f"RGroupDecompose failed: {exc}",
        }

    decomposition: list[dict[str, Any]] = []

    for row_idx, row in enumerate(matched_rows):
        original_idx = mol_indices[row_idx]
        # row is a dict like {"Core": "c1ccccc1", "R1": "[Cc1ccccc1:1]CC", ...}
        rgroups = {k: v for k, v in row.items() if k.startswith("R")}
        decomposition.append(
            {
                "molecule_index": original_idx,
                "core": core_smarts,
                "rgroups": rgroups,
            }
        )

    matched_count = len(decomposition)
    unmatched_count = total - matched_count

    return {
        "core_smarts": core_smarts,
        "decomposition": decomposition,
        "unmatched_count": unmatched_count,
    }
