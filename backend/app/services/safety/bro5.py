"""
Beyond-Rule-of-5 (bRo5) Check (SAFE-03)

Evaluates molecules with MW > 500 (macrocycles, PROTACs, natural products)
against relaxed physicochemical thresholds. Returns N/A for standard-sized
molecules where Ro5 applies.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def compute_bro5(mol: Chem.Mol) -> dict:
    """Evaluate a molecule against beyond-Rule-of-5 (bRo5) criteria.

    bRo5 applies only to molecules with MW > 500. For smaller molecules,
    standard Ro5 applies and this check returns not applicable.

    bRo5 thresholds (Doak et al., 2014):
      - MW > 1000
      - LogP > 10
      - HBD > 6
      - TPSA > 250
      - RotBonds > 20

    Args:
        mol: RDKit molecule to evaluate.

    Returns:
        Dict with keys:
          - applicable (bool): True if MW > 500
          - passed (bool): True if no bRo5 violations (or not applicable)
          - message (str): Present when not applicable
          - violations (list[dict]): Each with property, value, threshold, direction
          - values (dict): Computed property values (empty when not applicable)
    """
    mw = Descriptors.MolWt(mol)

    if mw <= 500:
        return {
            "applicable": False,
            "passed": True,
            "message": "MW <= 500 -- standard Ro5 applies",
            "violations": [],
            "values": {},
        }

    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    tpsa = Descriptors.TPSA(mol)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    values = {
        "MW": mw,
        "LogP": logp,
        "HBD": hbd,
        "TPSA": tpsa,
        "RotBonds": rot_bonds,
    }

    violations: list[dict] = []

    checks = [
        ("MW", mw, 1000.0),
        ("LogP", logp, 10.0),
        ("HBD", float(hbd), 6.0),
        ("TPSA", tpsa, 250.0),
        ("RotBonds", float(rot_bonds), 20.0),
    ]

    for prop, value, threshold in checks:
        if value > threshold:
            violations.append(
                {
                    "property": prop,
                    "value": value,
                    "threshold": threshold,
                    "direction": ">",
                }
            )

    return {
        "applicable": True,
        "passed": len(violations) == 0,
        "violations": violations,
        "values": values,
    }
