"""
Multi-Parameter Optimization (MPO) Scoring

Implements two MPO frameworks:
1. CNS MPO preset — Wager et al. 2010 (Pfizer), piecewise linear desirability
2. Custom MPO — configurable desirability functions (sigmoid/ramp/step)

The CNS MPO implementation uses ONLY the 4 RDKit-computable properties
(cLogP, TPSA, MW, HBD). The full Wager 2010 paper includes 6 components
(adding pKa and CLogD), but those require additional tools not in scope.
Our max_score = 4.0 (not 6.0) to reflect the 4-component implementation.

IMPORTANT: CNS MPO uses PIECEWISE LINEAR desirability functions, not sigmoid.
This matches Wager 2010 exactly. See Pitfall 4 in RESEARCH.md.

References:
- Wager et al. ACS Chem Neurosci 2010;1:435-449 (CNS MPO)
- Bickerton et al. Nat Chem 2012;4:90-98 (general MPO)
"""

import math
from typing import Any

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Lipinski, rdMolDescriptors

# ---------------------------------------------------------------------------
# CNS MPO piecewise linear desirability functions — Wager 2010
# ---------------------------------------------------------------------------


def _cns_mpo_clogp(clogp: float) -> float:
    """
    CNS MPO desirability for cLogP (Wager 2010).

    Piecewise linear: <=1 → 1.0; >=5 → 0.0; between → linear decay.
    """
    if clogp <= 1:
        return 1.0
    if clogp >= 5:
        return 0.0
    return 1.0 - (clogp - 1) / 4.0


def _cns_mpo_tpsa(tpsa: float) -> float:
    """
    CNS MPO desirability for TPSA (Wager 2010).

    Piecewise: [40, 90] → 1.0; <20 or >120 → 0.0; linear ramps in between.
    """
    if 40 <= tpsa <= 90:
        return 1.0
    if tpsa < 20 or tpsa > 120:
        return 0.0
    if tpsa < 40:
        return (tpsa - 20) / 20.0
    # tpsa > 90 and tpsa <= 120
    return 1.0 - (tpsa - 90) / 30.0


def _cns_mpo_mw(mw: float) -> float:
    """
    CNS MPO desirability for MW (Wager 2010).

    Piecewise linear: <=360 → 1.0; >=500 → 0.0; between → linear decay.
    """
    if mw <= 360:
        return 1.0
    if mw >= 500:
        return 0.0
    return 1.0 - (mw - 360) / 140.0


def _cns_mpo_hbd(hbd: int) -> float:
    """
    CNS MPO desirability for HBD (Wager 2010).

    Returns max(0, 1 - hbd * 0.25): 0→1.0, 1→0.75, 2→0.5, 3→0.25, 4+→0.0.
    """
    return max(0.0, 1.0 - hbd * 0.25)


# ---------------------------------------------------------------------------
# CNS MPO main function — PROF-07
# ---------------------------------------------------------------------------


def compute_cns_mpo(mol: Chem.Mol) -> dict:
    """
    Compute CNS MPO score using Wager 2010 piecewise linear functions.

    Computes 4 desirability components (cLogP, TPSA, MW, HBD) and sums them.
    Max score = 4.0 (reflecting 4-component RDKit-available implementation).

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict with keys: score, max_score, components (clogp, tpsa, mw, hbd).
    """
    clogp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    mw = Descriptors.MolWt(mol)
    hbd = int(Lipinski.NumHDonors(mol))

    d_clogp = _cns_mpo_clogp(clogp)
    d_tpsa = _cns_mpo_tpsa(tpsa)
    d_mw = _cns_mpo_mw(mw)
    d_hbd = _cns_mpo_hbd(hbd)

    score = d_clogp + d_tpsa + d_mw + d_hbd

    return {
        "score": round(score, 2),
        "max_score": 4.0,
        "components": {
            "clogp": round(d_clogp, 3),
            "tpsa": round(d_tpsa, 3),
            "mw": round(d_mw, 3),
            "hbd": round(d_hbd, 3),
        },
    }


# ---------------------------------------------------------------------------
# Custom MPO framework — PROF-07
# ---------------------------------------------------------------------------


def _desirability_value(value: float, low: float, high: float, shape: str) -> float:
    """
    Compute a single desirability value using the specified curve shape.

    Shapes:
    - sigmoid: S-shaped transition from low to high boundary
    - ramp:    linear ramp from 0 at low to 1 at high
    - step:    0 below high threshold, 1 at or above

    Args:
        value: Actual property value.
        low:   Lower boundary of the favorable range.
        high:  Upper boundary of the favorable range.
        shape: One of "sigmoid", "ramp", "step".

    Returns:
        Desirability value in [0, 1].
    """
    if shape == "sigmoid":
        if high == low:
            return 1.0 if value >= high else 0.0
        x = (value - low) / (high - low)
        clamped = max(0.0, min(1.0, x))
        return 1.0 / (1.0 + math.exp(-12.0 * (clamped - 0.5)))
    elif shape == "ramp":
        if high == low:
            return 1.0 if value >= high else 0.0
        return max(0.0, min(1.0, (value - low) / (high - low)))
    elif shape == "step":
        return 1.0 if value >= high else 0.0
    else:
        # Unknown shape — default to ramp
        if high == low:
            return 1.0 if value >= high else 0.0
        return max(0.0, min(1.0, (value - low) / (high - low)))


def _get_descriptor_value(mol: Chem.Mol, property_name: str) -> float:
    """
    Get a named molecular descriptor value using RDKit.

    Supports: MW, LogP, TPSA, HBD, HBA, RotBonds, Fsp3, QED, SA, PFI.

    Args:
        mol: RDKit molecule object.
        property_name: Descriptor name (case-sensitive).

    Returns:
        Float descriptor value, or 0.0 if unknown.
    """
    dispatch = {
        "MW": lambda m: Descriptors.MolWt(m),
        "LogP": lambda m: Descriptors.MolLogP(m),
        "TPSA": lambda m: Descriptors.TPSA(m),
        "HBD": lambda m: float(Lipinski.NumHDonors(m)),
        "HBA": lambda m: float(Lipinski.NumHAcceptors(m)),
        "RotBonds": lambda m: float(Lipinski.NumRotatableBonds(m)),
        "Fsp3": lambda m: rdMolDescriptors.CalcFractionCSP3(m),
        "QED": lambda m: QED.qed(m),
        "SA": _get_sa_score,
        "PFI": lambda m: Descriptors.MolLogP(m) + Descriptors.NumAromaticRings(m),
        "NumRings": lambda m: float(rdMolDescriptors.CalcNumRings(m)),
        "HeavyAtoms": lambda m: float(m.GetNumHeavyAtoms()),
    }
    func = dispatch.get(property_name)
    if func is None:
        return 0.0
    try:
        return float(func(mol))
    except Exception:
        return 0.0


def _get_sa_score(mol: Chem.Mol) -> float:
    """Get SA_Score, with graceful fallback."""
    try:
        from rdkit.Contrib.SA_Score import sascorer  # type: ignore

        return sascorer.calculateScore(mol)
    except ImportError:
        pass
    try:
        import os
        import sys

        from rdkit import RDConfig

        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer  # type: ignore

        return sascorer.calculateScore(mol)
    except Exception:
        return 5.0  # Neutral fallback


def compute_custom_mpo(mol: Chem.Mol, profile: list[dict]) -> dict:
    """
    Compute a custom MPO score using user-defined desirability functions.

    Each entry in profile defines one property:
    - property: Name of the descriptor (MW, LogP, TPSA, etc.)
    - low:      Lower boundary of favorable range
    - high:     Upper boundary of favorable range
    - weight:   Relative importance weight (default 1.0)
    - shape:    Curve type: "sigmoid" | "ramp" | "step" (default "sigmoid")

    The normalized score = sum(desirability * weight) / sum(weight).

    Args:
        mol: RDKit molecule object.
        profile: List of property dicts (see MPOProperty schema).

    Returns:
        Dict with keys: score, max_score, normalized, components.
    """
    weighted_sum = 0.0
    total_weight = 0.0
    components = []

    for entry in profile:
        prop_name = entry.get("property", "")
        low = float(entry.get("low", 0))
        high = float(entry.get("high", 1))
        weight = float(entry.get("weight", 1.0))
        shape = entry.get("shape", "sigmoid")

        value = _get_descriptor_value(mol, prop_name)
        desirability = _desirability_value(value, low, high, shape)
        weighted_contribution = desirability * weight

        weighted_sum += weighted_contribution
        total_weight += weight

        components.append(
            {
                "property": prop_name,
                "value": round(value, 4),
                "desirability": round(desirability, 4),
                "weight": weight,
                "contribution": round(weighted_contribution, 4),
            }
        )

    normalized = weighted_sum / total_weight if total_weight > 0 else 0.0

    return {
        "score": round(weighted_sum, 4),
        "max_score": round(total_weight, 4),
        "normalized": round(normalized, 4),
        "components": components,
    }


# ---------------------------------------------------------------------------
# Oral Drug MPO preset — D-12
# ---------------------------------------------------------------------------

ORAL_DRUG_MPO_PRESET: list[dict[str, Any]] = [
    {"property": "MW", "low": 200, "high": 500, "weight": 1.0, "shape": "ramp"},
    {"property": "LogP", "low": -0.4, "high": 5.6, "weight": 1.0, "shape": "ramp"},
    {"property": "TPSA", "low": 20, "high": 130, "weight": 1.0, "shape": "ramp"},
    {"property": "HBD", "low": 0, "high": 5, "weight": 1.0, "shape": "ramp"},
    {"property": "RotBonds", "low": 0, "high": 10, "weight": 1.0, "shape": "ramp"},
]
