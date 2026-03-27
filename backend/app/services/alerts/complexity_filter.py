"""
Complexity Percentile Filter

Evaluates whether a molecule's physicochemical properties fall within the
5th-95th percentile range of commercial drug-like compound distributions.
Properties outside this range are flagged as outliers.

Six properties assessed:
  - MW          Molecular weight (Daltons)
  - LogP        Wildman-Crippen LogP estimate
  - RotBonds    Number of rotatable bonds
  - Rings       Total ring count
  - HeavyAtoms  Number of heavy (non-hydrogen) atoms
  - BertzCT     Bertz complexity index

Thresholds derived from analysis of commercial compound library distributions
(Enamine REAL + ChEMBL drug-like subset, n ≈ 3M compounds).
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors

# ---------------------------------------------------------------------------
# 5th–95th percentile thresholds for 6 physicochemical properties
# ---------------------------------------------------------------------------
COMPLEXITY_THRESHOLDS: dict[str, dict[str, float]] = {
    "MW": {"p5": 150.0, "p95": 550.0},
    "LogP": {"p5": -2.0, "p95": 5.5},
    "RotBonds": {"p5": 0.0, "p95": 12.0},
    "Rings": {"p5": 1.0, "p95": 5.0},
    "HeavyAtoms": {"p5": 10.0, "p95": 40.0},
    "BertzCT": {"p5": 100.0, "p95": 1200.0},
}


def compute_complexity_percentile(mol: Chem.Mol) -> dict:
    """
    Compute complexity percentile assessment for a molecule.

    Calculates six physicochemical properties and checks each against
    the 5th-95th percentile thresholds for commercial compound distributions.

    Args:
        mol: RDKit molecule (must be valid and sanitized).

    Returns:
        Dictionary with the following structure::

            {
                "properties": {
                    "MW": {
                        "value": float,
                        "p5": float,
                        "p95": float,
                        "outlier": bool,
                        "direction": str | None,   # "below", "above", or None
                    },
                    ...  # same for LogP, RotBonds, Rings, HeavyAtoms, BertzCT
                },
                "n_outliers": int,
                "outlier_properties": list[str],
                "within_range": bool,   # True iff n_outliers == 0
            }
    """
    # Compute raw property values
    raw_values: dict[str, float] = {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "RotBonds": float(Descriptors.NumRotatableBonds(mol)),
        "Rings": float(Descriptors.RingCount(mol)),
        "HeavyAtoms": float(mol.GetNumHeavyAtoms()),
        "BertzCT": Descriptors.BertzCT(mol),
    }

    properties: dict[str, dict] = {}
    outlier_properties: list[str] = []

    for prop_name, value in raw_values.items():
        thresholds = COMPLEXITY_THRESHOLDS[prop_name]
        p5 = thresholds["p5"]
        p95 = thresholds["p95"]

        if value < p5:
            outlier = True
            direction: str | None = "below"
        elif value > p95:
            outlier = True
            direction = "above"
        else:
            outlier = False
            direction = None

        if outlier:
            outlier_properties.append(prop_name)

        properties[prop_name] = {
            "value": value,
            "p5": p5,
            "p95": p95,
            "outlier": outlier,
            "direction": direction,
        }

    n_outliers = len(outlier_properties)

    return {
        "properties": properties,
        "n_outliers": n_outliers,
        "outlier_properties": outlier_properties,
        "within_range": n_outliers == 0,
    }
