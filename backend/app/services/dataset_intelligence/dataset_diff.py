"""
Dataset Diff Engine Service Module (DSET-03).

Compares two datasets by InChIKey set operations to identify added,
removed, modified, and unchanged molecules. Property-level changes
are tracked for modified molecules. Handles column mismatches between
datasets by comparing only shared columns (Pitfall 7).
"""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)


def _compute_property_changes(
    primary_mol: dict,
    comparison_mol: dict,
    compare_columns: list[str],
) -> list[dict]:
    """Compare property values between two molecules for shared columns.

    Skips 'smiles' and 'inchikey' columns (structural identity already
    matched by InChIKey).

    Args:
        primary_mol: Primary molecule dict with 'properties' key.
        comparison_mol: Comparison molecule dict with 'properties' key.
        compare_columns: List of column names to compare.

    Returns:
        List of change dicts: ``[{"column": str, "old_value": Any,
        "new_value": Any}]`` for columns that differ.
    """
    changes: list[dict] = []
    primary_props = primary_mol.get("properties", {})
    comparison_props = comparison_mol.get("properties", {})

    for col in compare_columns:
        if col in ("smiles", "inchikey"):
            continue
        old_val = primary_props.get(col)
        new_val = comparison_props.get(col)
        if old_val != new_val:
            changes.append({
                "column": col,
                "old_value": old_val,
                "new_value": new_val,
            })

    return changes


def compute_dataset_diff(
    primary_molecules: list[dict],
    comparison_molecules: list[dict],
) -> dict:
    """Compute InChIKey-based diff between two datasets.

    Identifies added, removed, modified, and unchanged molecules.
    Only compares property columns present in both datasets (per Pitfall 7).

    Args:
        primary_molecules: List of molecule dicts from the primary dataset.
            Each dict: ``{"index": int, "smiles": str, "inchikey": str | None,
            "properties": dict}``
        comparison_molecules: List of molecule dicts from the comparison dataset.

    Returns:
        Dict with keys:
            added (list[dict]): Molecules in comparison but not primary.
            removed (list[dict]): Molecules in primary but not comparison.
            modified (list[dict]): Molecules in both with property changes.
            added_count (int), removed_count (int), modified_count (int),
            unchanged_count (int), unique_columns_primary (int),
            unique_columns_comparison (int).
    """
    # Build lookup dicts by InChIKey, skipping None
    primary_by_ik: dict[str, dict] = {}
    for mol in primary_molecules:
        ik = mol.get("inchikey")
        if ik is not None:
            primary_by_ik[ik] = mol

    comparison_by_ik: dict[str, dict] = {}
    for mol in comparison_molecules:
        ik = mol.get("inchikey")
        if ik is not None:
            comparison_by_ik[ik] = mol

    primary_keys = set(primary_by_ik.keys())
    comparison_keys = set(comparison_by_ik.keys())

    added_keys = comparison_keys - primary_keys
    removed_keys = primary_keys - comparison_keys
    common_keys = primary_keys & comparison_keys

    # Determine shared columns for comparison (exclude smiles/inchikey)
    primary_columns: set[str] = set()
    for mol in primary_by_ik.values():
        primary_columns.update(mol.get("properties", {}).keys())

    comparison_columns: set[str] = set()
    for mol in comparison_by_ik.values():
        comparison_columns.update(mol.get("properties", {}).keys())

    shared_columns = primary_columns & comparison_columns
    shared_columns.discard("smiles")
    shared_columns.discard("inchikey")

    unique_primary = primary_columns - comparison_columns
    unique_comparison = comparison_columns - primary_columns

    compare_columns = sorted(shared_columns)

    # Build result lists
    added: list[dict] = []
    for ik in sorted(added_keys):
        mol = comparison_by_ik[ik]
        added.append({
            "inchikey": ik,
            "smiles": mol.get("smiles", ""),
            "row_index": mol.get("index", 0),
            "properties": mol.get("properties", {}),
        })

    removed: list[dict] = []
    for ik in sorted(removed_keys):
        mol = primary_by_ik[ik]
        removed.append({
            "inchikey": ik,
            "smiles": mol.get("smiles", ""),
            "row_index": mol.get("index", 0),
            "properties": mol.get("properties", {}),
        })

    modified: list[dict] = []
    unchanged_count = 0
    for ik in sorted(common_keys):
        primary_mol = primary_by_ik[ik]
        comparison_mol = comparison_by_ik[ik]
        changes = _compute_property_changes(primary_mol, comparison_mol, compare_columns)
        if changes:
            modified.append({
                "inchikey": ik,
                "smiles": primary_mol.get("smiles", ""),
                "row_index": primary_mol.get("index", 0),
                "changes": changes,
                "properties": comparison_mol.get("properties", {}),
            })
        else:
            unchanged_count += 1

    return {
        "added": added,
        "removed": removed,
        "modified": modified,
        "added_count": len(added),
        "removed_count": len(removed),
        "modified_count": len(modified),
        "unchanged_count": unchanged_count,
        "unique_columns_primary": len(unique_primary),
        "unique_columns_comparison": len(unique_comparison),
    }
