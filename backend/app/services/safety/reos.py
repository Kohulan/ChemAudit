"""
REOS (Rapid Elimination of Swill) Filter (SAFE-04)

Evaluates molecules against 7 physicochemical property ranges to rapidly
identify undesirable compounds. Originally published by Walters et al. (1999).
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# REOS property ranges: property -> (lower_bound, upper_bound) inclusive
REOS_RANGES: dict[str, tuple[float, float]] = {
    "MW": (200.0, 500.0),
    "LogP": (-5.0, 5.0),
    "HBD": (0.0, 5.0),
    "HBA": (0.0, 10.0),
    "RotBonds": (0.0, 8.0),
    "TPSA": (0.0, 150.0),
    "Rings": (0.0, 4.0),
}


def compute_reos(mol: Chem.Mol) -> dict:
    """Evaluate a molecule against the REOS 7-property filter.

    Checks each property against its acceptable range. A molecule passes REOS
    only if all 7 properties are within their defined ranges.

    Args:
        mol: RDKit molecule to evaluate.

    Returns:
        Dict with keys:
          - passed (bool): True if no violations
          - violations (list[dict]): Each with property, value, range, exceeded
          - n_violations (int): Number of property violations
          - descriptors (dict): Computed values for all 7 properties
    """
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = float(rdMolDescriptors.CalcNumHBD(mol))
    hba = float(rdMolDescriptors.CalcNumHBA(mol))
    rot_bonds = float(rdMolDescriptors.CalcNumRotatableBonds(mol))
    tpsa = Descriptors.TPSA(mol)
    rings = float(rdMolDescriptors.CalcNumRings(mol))

    descriptors = {
        "MW": mw,
        "LogP": logp,
        "HBD": hbd,
        "HBA": hba,
        "RotBonds": rot_bonds,
        "TPSA": tpsa,
        "Rings": rings,
    }

    violations: list[dict] = []

    for prop, (lo, hi) in REOS_RANGES.items():
        value = descriptors[prop]
        if not (lo <= value <= hi):
            violations.append(
                {
                    "property": prop,
                    "value": value,
                    "range": [lo, hi],
                    "exceeded": value > hi,
                }
            )

    return {
        "passed": len(violations) == 0,
        "violations": violations,
        "n_violations": len(violations),
        "descriptors": descriptors,
    }
