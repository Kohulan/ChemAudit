"""
Analytics Statistics Module

Computes property distribution statistics, pairwise Pearson correlations,
IQR-based outlier detection, and a composite batch quality score from batch results.

The `compute_all_statistics` function is the primary entry point and is called
by `run_cheap_analytics` in analytics_tasks.py after batch aggregation completes.
"""

from __future__ import annotations

import logging
import math
from typing import Callable

import numpy as np

from app.schemas.analytics import (
    OutlierInfo,
    PropertyCorrelation,
    PropertyStats,
    QualityScore,
    StatisticsResult,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Property extraction constants
# ---------------------------------------------------------------------------

# Lambdas that extract a numeric value from a single molecule result dict.
# Each lambda must return a float/int or None (None means "not available").
PROPERTY_EXTRACTORS: dict[str, Callable[[dict], float | None]] = {
    "validation_score": lambda r: (r.get("validation") or {}).get("overall_score"),
    "qed_score": lambda r: ((r.get("scoring") or {}).get("druglikeness") or {}).get(
        "qed_score"
    ),
    "sa_score": lambda r: ((r.get("scoring") or {}).get("admet") or {}).get("sa_score"),
    "ml_readiness_score": lambda r: (
        (r.get("scoring") or {}).get("ml_readiness") or {}
    ).get("score"),
    "fsp3": lambda r: ((r.get("scoring") or {}).get("admet") or {}).get("fsp3"),
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _extract_property_values(
    results: list[dict],
    prop_name: str,
    extractor: Callable[[dict], float | None],
) -> list[tuple[int, float]]:
    """
    Extract numeric property values from a list of molecule result dicts.

    Iterates over all results (regardless of success/error status) and attempts
    to extract the property value using the provided extractor callable.

    Args:
        results: List of molecule result dicts from batch processing.
        prop_name: Human-readable property name (used only for logging).
        extractor: Callable that maps a result dict to a float value or None.

    Returns:
        List of (molecule_index, value) tuples for results where the value
        was successfully extracted and is not None.
    """
    extracted: list[tuple[int, float]] = []
    for idx, result in enumerate(results):
        try:
            value = extractor(result)
            if value is not None:
                extracted.append((idx, float(value)))
        except Exception:
            logger.debug("Failed to extract %s from result at index %d", prop_name, idx)
    return extracted


# ---------------------------------------------------------------------------
# Public computation functions
# ---------------------------------------------------------------------------


def compute_property_stats(results: list[dict]) -> list[PropertyStats]:
    """
    Compute descriptive statistics for each numeric molecular property.

    For each property in PROPERTY_EXTRACTORS, extracts available values and
    computes mean, median, std, quartiles, IQR, min, and max using numpy.
    Properties with fewer than 2 available values are skipped.

    Args:
        results: List of molecule result dicts from batch processing.

    Returns:
        List of PropertyStats models, one per property with sufficient data.
    """
    stats_list: list[PropertyStats] = []

    for prop_name, extractor in PROPERTY_EXTRACTORS.items():
        pairs = _extract_property_values(results, prop_name, extractor)
        if len(pairs) < 2:
            continue

        vals = np.array([v for _, v in pairs], dtype=float)
        mean = float(np.mean(vals))
        median = float(np.median(vals))
        std = float(np.std(vals, ddof=1))
        q1, q3 = float(np.percentile(vals, 25)), float(np.percentile(vals, 75))
        iqr = q3 - q1
        min_val = float(np.min(vals))
        max_val = float(np.max(vals))

        stats_list.append(
            PropertyStats(
                property_name=prop_name,
                mean=mean,
                median=median,
                std=std,
                q1=q1,
                q3=q3,
                iqr=iqr,
                min=min_val,
                max=max_val,
                count=len(vals),
            )
        )

    return stats_list


def compute_correlations(results: list[dict]) -> list[PropertyCorrelation]:
    """
    Compute pairwise Pearson correlations for all numeric molecular properties.

    Only includes properties with at least 10 available values. Handles edge
    cases where one or both arrays are constant (corrcoef returns NaN) by
    substituting 0.0. Skips pairs where the absolute correlation rounds to 0.

    Args:
        results: List of molecule result dicts from batch processing.

    Returns:
        List of PropertyCorrelation models for all pairs with non-zero correlation.
    """
    # Gather aligned arrays for each property
    prop_arrays: dict[str, np.ndarray] = {}
    prop_indices: dict[str, set[int]] = {}

    for prop_name, extractor in PROPERTY_EXTRACTORS.items():
        pairs = _extract_property_values(results, prop_name, extractor)
        if len(pairs) < 10:
            continue
        prop_indices[prop_name] = {idx for idx, _ in pairs}
        # Build full-length array with NaN for missing indices, but for
        # pairwise correlation we need the same set of molecules — use
        # a dict keyed by molecule index for alignment.
        prop_arrays[prop_name] = np.array([v for _, v in pairs], dtype=float)

    # For pairwise correlation we need the same molecules in both arrays.
    # Build a lookup: {prop_name: {mol_idx: value}}
    prop_value_maps: dict[str, dict[int, float]] = {}
    for prop_name, extractor in PROPERTY_EXTRACTORS.items():
        pairs = _extract_property_values(results, prop_name, extractor)
        if len(pairs) < 10:
            continue
        prop_value_maps[prop_name] = {idx: v for idx, v in pairs}

    prop_names = list(prop_value_maps.keys())
    correlations: list[PropertyCorrelation] = []

    for i, name_a in enumerate(prop_names):
        for name_b in prop_names[i + 1 :]:
            map_a = prop_value_maps[name_a]
            map_b = prop_value_maps[name_b]

            # Use only indices present in both properties
            common_indices = sorted(set(map_a.keys()) & set(map_b.keys()))
            if len(common_indices) < 2:
                continue

            arr_a = np.array([map_a[i] for i in common_indices], dtype=float)
            arr_b = np.array([map_b[i] for i in common_indices], dtype=float)

            corr_matrix = np.corrcoef(arr_a, arr_b)
            pearson_r = float(corr_matrix[0, 1])

            # Replace NaN (constant arrays) with 0.0
            if math.isnan(pearson_r):
                pearson_r = 0.0

            # Skip zero correlations for brevity
            if pearson_r == 0.0:
                continue

            correlations.append(
                PropertyCorrelation(
                    property_a=name_a,
                    property_b=name_b,
                    pearson_r=round(pearson_r, 4),
                )
            )

    return correlations


def compute_outliers(results: list[dict]) -> list[OutlierInfo]:
    """
    Detect outlier molecules using the IQR fence method.

    For each property in PROPERTY_EXTRACTORS, computes IQR fences:
        lower = Q1 - 1.5 * IQR
        upper = Q3 + 1.5 * IQR

    Molecules whose property value falls outside these fences are flagged.
    Properties with fewer than 4 available values are skipped.

    Args:
        results: List of molecule result dicts from batch processing.

    Returns:
        List of OutlierInfo models for each (molecule, property) pair flagged.
    """
    outliers: list[OutlierInfo] = []

    for prop_name, extractor in PROPERTY_EXTRACTORS.items():
        pairs = _extract_property_values(results, prop_name, extractor)
        if len(pairs) < 4:
            continue

        vals = np.array([v for _, v in pairs], dtype=float)
        q1, q3 = float(np.percentile(vals, 25)), float(np.percentile(vals, 75))
        iqr = q3 - q1
        lower_fence = q1 - 1.5 * iqr
        upper_fence = q3 + 1.5 * iqr

        for mol_idx, value in pairs:
            if value < lower_fence or value > upper_fence:
                outliers.append(
                    OutlierInfo(
                        molecule_index=mol_idx,
                        property_name=prop_name,
                        value=round(value, 4),
                        lower_fence=round(lower_fence, 4),
                        upper_fence=round(upper_fence, 4),
                    )
                )

    return outliers


def compute_quality_score(results: list[dict]) -> QualityScore:
    """
    Compute a composite 0-100 batch quality score with three weighted components.

    Components and weights:
    - Validity (40%): Fraction of results with status="success" * 100.
    - Diversity (35%): Shannon entropy of Murcko scaffold frequency distribution,
      normalized by log2(unique_scaffold_count). Scaffolds are extracted inline
      using RDKit (does not depend on the scaffold_analysis service).
    - Drug-likeness (25%): Fraction of molecules with lipinski_passed=True * 100.
      Falls back to 50% (neutral) if no Lipinski data is available.

    Composite: score = validity_pct*0.40 + diversity_pct*0.35 + druglikeness_pct*0.25

    Args:
        results: List of molecule result dicts from batch processing.

    Returns:
        QualityScore model with rounded component and composite values.
    """
    total = len(results)
    if total == 0:
        return QualityScore(
            score=0.0,
            validity_pct=0.0,
            diversity_pct=0.0,
            druglikeness_pct=0.0,
        )

    # --- Validity (40%) ---
    successful = sum(1 for r in results if r.get("status") == "success")
    validity_pct = (successful / total) * 100.0

    # --- Diversity (35%): inline Murcko scaffold extraction ---
    diversity_pct = 0.0
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        scaffold_counts: dict[str, int] = {}
        for result in results:
            if result.get("status") != "success":
                continue
            smiles = result.get("smiles") or result.get("standardized_smiles")
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
            if scaffold_mol is None:
                continue
            scaffold_smi = Chem.MolToSmiles(scaffold_mol)
            scaffold_counts[scaffold_smi] = scaffold_counts.get(scaffold_smi, 0) + 1

        unique_scaffolds = len(scaffold_counts)
        if unique_scaffolds > 1:
            total_scaffold_mols = sum(scaffold_counts.values())
            entropy = 0.0
            for count in scaffold_counts.values():
                p = count / total_scaffold_mols
                if p > 0:
                    entropy -= p * math.log2(p)
            max_entropy = math.log2(unique_scaffolds)
            diversity_pct = (entropy / max_entropy) * 100.0 if max_entropy > 0 else 0.0
        else:
            diversity_pct = 0.0

    except Exception:
        logger.warning(
            "compute_quality_score: scaffold diversity computation failed; defaulting to 0",
            exc_info=True,
        )
        diversity_pct = 0.0

    # --- Drug-likeness (25%): Lipinski pass rate ---
    lipinski_passes = 0
    lipinski_total = 0
    for result in results:
        if result.get("status") != "success":
            continue
        druglikeness = ((result.get("scoring") or {}).get("druglikeness") or {})
        lipinski_passed = druglikeness.get("lipinski_passed")
        if lipinski_passed is not None:
            lipinski_total += 1
            if lipinski_passed:
                lipinski_passes += 1

    if lipinski_total > 0:
        druglikeness_pct = (lipinski_passes / lipinski_total) * 100.0
    else:
        druglikeness_pct = 50.0  # Neutral fallback when no Lipinski data available

    # --- Composite score ---
    score = validity_pct * 0.40 + diversity_pct * 0.35 + druglikeness_pct * 0.25

    return QualityScore(
        score=round(score, 1),
        validity_pct=round(validity_pct, 1),
        diversity_pct=round(diversity_pct, 1),
        druglikeness_pct=round(druglikeness_pct, 1),
    )


def compute_all_statistics(results: list[dict]) -> StatisticsResult:
    """
    Compute all batch statistics and return a StatisticsResult model.

    Aggregates:
    - Property distribution statistics (mean, median, std, quartiles, IQR, min, max)
    - Pairwise Pearson correlations for numeric properties
    - IQR-based outlier flags
    - Composite batch quality score (40% validity / 35% diversity / 25% drug-likeness)

    Called by `run_cheap_analytics` in analytics_tasks.py immediately after batch
    aggregation completes (no user trigger required).

    Args:
        results: List of molecule result dicts from batch processing.

    Returns:
        StatisticsResult Pydantic model combining all four sub-results.
    """
    property_stats = compute_property_stats(results)
    correlations = compute_correlations(results)
    outliers = compute_outliers(results)
    quality_score = compute_quality_score(results)

    return StatisticsResult(
        property_stats=property_stats,
        correlations=correlations,
        outliers=outliers,
        quality_score=quality_score,
    )
