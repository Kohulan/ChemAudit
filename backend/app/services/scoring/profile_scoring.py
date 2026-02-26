"""
Profile-based desirability scoring.

Computes per-property desirability (0-1) and a weighted geometric mean
profile score (0-100) following QED-style methodology (Bickerton et al., 2012).
"""

from math import prod
from typing import Any, Optional


def desirability(
    value: Optional[float],
    min_val: float,
    max_val: float,
) -> Optional[float]:
    """Compute 0-1 desirability for a property value given an acceptable range.

    Returns 1.0 if value is within [min_val, max_val].
    Linear falloff outside the range, reaching 0.0 at 100% beyond the range width.
    """
    if value is None:
        return None
    if min_val <= value <= max_val:
        return 1.0
    range_width = max(max_val - min_val, 1e-6)
    if value < min_val:
        distance = (min_val - value) / range_width
    else:
        distance = (value - max_val) / range_width
    return max(0.0, 1.0 - distance)


def profile_score(
    properties: dict[str, Any],
    thresholds: dict[str, dict[str, Any]],
    weights: dict[str, float],
) -> Optional[float]:
    """Compute weighted geometric mean of desirabilities scaled to 0-100.

    Args:
        properties: Property name -> value (e.g. {"mw": 342.1, "logp": 2.3})
        thresholds: Property name -> {"min": float, "max": float}
        weights: Property name -> weight (default 1.0)

    Returns:
        Score 0-100, or None if no evaluable properties.
    """
    scores: list[float] = []
    total_weight = 0.0
    for key, thresh in thresholds.items():
        val = properties.get(key)
        if val is None:
            continue
        t_min = thresh.get("min", float("-inf"))
        t_max = thresh.get("max", float("inf"))
        d = desirability(val, t_min, t_max)
        if d is None:
            continue
        w = weights.get(key, 1.0)
        scores.append(d**w)
        total_weight += w

    if not scores or total_weight == 0:
        return None

    geometric_mean = prod(scores) ** (1.0 / total_weight)
    return round(geometric_mean * 100, 1)


def compute_profile_result(
    properties: dict[str, Any],
    profile_id: int,
    profile_name: str,
    thresholds: dict[str, dict[str, Any]],
    weights: dict[str, float],
) -> dict[str, Any]:
    """Build the full profile scoring result dict for a single molecule.

    Returns dict with score, profile metadata, and per-property breakdown.
    """
    score = profile_score(properties, thresholds, weights)
    prop_details: dict[str, Any] = {}
    for key, thresh in thresholds.items():
        val = properties.get(key)
        t_min = thresh.get("min", float("-inf"))
        t_max = thresh.get("max", float("inf"))
        d = desirability(val, t_min, t_max) if val is not None else None
        in_range = (t_min <= val <= t_max) if val is not None else None
        prop_details[key] = {
            "value": val,
            "min": t_min if t_min != float("-inf") else None,
            "max": t_max if t_max != float("inf") else None,
            "in_range": in_range,
            "desirability": round(d, 3) if d is not None else None,
        }

    return {
        "profile_id": profile_id,
        "profile_name": profile_name,
        "score": score,
        "properties": prop_details,
    }
