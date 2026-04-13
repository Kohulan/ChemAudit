"""
Contradictory Label Detection Service Module (DSET-02).

Detects molecules with the same InChIKey but divergent activity values,
indicating potential data quality issues (e.g., mislabeled compounds,
assay errors, or SAR cliffs).

Implements fold-difference calculation with division-by-zero safety
(Pitfall 5: min value clamped to 1e-10).
"""

from __future__ import annotations

import logging
from collections import defaultdict

logger = logging.getLogger(__name__)

# Priority-1 activity column names (exact match, case-insensitive)
_EXACT_ACTIVITY_NAMES = {
    "activity", "pic50", "ic50", "ki", "ec50", "kd", "potency",
}

# Priority-2 activity column substrings (case-insensitive)
_ACTIVITY_SUBSTRINGS = {"activity", "ic50"}


def detect_contradictory_labels(
    molecules: list[dict],
    activity_column: str,
    fold_threshold: float = 10.0,
) -> list[dict]:
    """Detect molecules with the same InChIKey but divergent activity values.

    Groups molecules by InChIKey and computes fold-difference between max and
    min activity values within each group. Returns contradictions where the
    fold-difference exceeds the threshold.

    Args:
        molecules: List of molecule dicts, each with keys:
            ``{"index": int, "smiles": str, "inchikey": str | None,
              "properties": dict}``
        activity_column: Name of the activity column in properties.
        fold_threshold: Minimum fold-difference to flag as contradictory.
            Default 10.0 (per D-11).

    Returns:
        List of contradiction dicts sorted by fold_difference descending:
        ``[{"inchikey": str, "entries": [{"row_index": int, "smiles": str,
           "activity": float}], "fold_difference": float, "entry_count": int,
           "smiles": str}]``
    """
    # Group molecules by InChIKey
    ik_groups: dict[str, list[dict]] = defaultdict(list)
    for mol_dict in molecules:
        ik = mol_dict.get("inchikey")
        if ik is None:
            continue

        # Extract activity value
        props = mol_dict.get("properties", {})
        raw_value = props.get(activity_column)
        if raw_value is None:
            continue

        try:
            value = float(raw_value)
        except (ValueError, TypeError):
            continue

        # Skip zero and negative values (Pitfall 5: division safety)
        if value <= 0:
            continue

        ik_groups[ik].append({
            "row_index": mol_dict.get("index", 0),
            "smiles": mol_dict.get("smiles", ""),
            "activity": value,
        })

    # Find contradictions
    contradictions: list[dict] = []
    for ik, entries in ik_groups.items():
        if len(entries) < 2:
            continue

        activities = [e["activity"] for e in entries]
        max_val = max(activities)
        min_val = min(activities)

        # Division-by-zero safety: clamp denominator to 1e-10
        fold_diff = max_val / max(min_val, 1e-10)

        if fold_diff > fold_threshold:
            contradictions.append({
                "inchikey": ik,
                "entries": sorted(entries, key=lambda e: e["activity"]),
                "fold_difference": round(fold_diff, 1),
                "entry_count": len(entries),
                "smiles": entries[0]["smiles"],
            })

    # Sort by fold_difference descending
    contradictions.sort(key=lambda c: c["fold_difference"], reverse=True)

    return contradictions


def detect_numeric_columns(
    column_names: list[str],
    sample_values: dict[str, list],
) -> list[dict]:
    """Detect and prioritize numeric columns from dataset column names.

    Checks if >80% of non-empty sample values are numeric and assigns
    priority based on activity-like naming heuristics:
        1 = exact match (activity, pIC50, IC50, Ki, EC50, Kd, potency)
        2 = substring match containing 'activity' or 'IC50'
        3 = other numeric column

    Args:
        column_names: List of column names to evaluate.
        sample_values: Dict mapping column names to sample value lists.

    Returns:
        List of ``{"name": str, "priority": int}`` sorted by priority asc,
        then alphabetical.
    """
    result: list[dict] = []

    for col_name in column_names:
        values = sample_values.get(col_name, [])
        if not values:
            continue

        # Check if >80% of non-empty values are numeric
        non_empty = [v for v in values if v is not None and str(v).strip() != ""]
        if not non_empty:
            continue

        numeric_count = 0
        for v in non_empty:
            try:
                float(v)
                numeric_count += 1
            except (ValueError, TypeError):
                pass

        if numeric_count / len(non_empty) <= 0.8:
            continue

        # Assign priority
        col_lower = col_name.lower()
        if col_lower in _EXACT_ACTIVITY_NAMES:
            priority = 1
        elif any(sub in col_lower for sub in _ACTIVITY_SUBSTRINGS):
            priority = 2
        else:
            priority = 3

        result.append({"name": col_name, "priority": priority})

    # Sort by priority ascending, then alphabetical
    result.sort(key=lambda r: (r["priority"], r["name"]))

    return result
