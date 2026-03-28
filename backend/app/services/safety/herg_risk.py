"""
hERG Liability Assessment (SAFE-02)

Rule-based 4-factor amphiphile hERG risk assessment using physicochemical
descriptors. Returns a risk score from 0-4 and human-readable flag descriptions.
"""
from __future__ import annotations

from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

# Module-level compiled singleton for basic nitrogen SMARTS
_BASIC_N_SMARTS: Optional[Chem.Mol] = None


def _get_basic_n_pattern() -> Chem.Mol:
    """Return lazily compiled basic nitrogen SMARTS pattern."""
    global _BASIC_N_SMARTS
    if _BASIC_N_SMARTS is None:
        _BASIC_N_SMARTS = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")
    return _BASIC_N_SMARTS


def compute_herg_risk(mol: Chem.Mol) -> dict:
    """Compute hERG liability risk using the 4-factor amphiphile rule.

    Evaluates four physicochemical factors known to correlate with hERG channel
    blockade: logP, molecular weight in the amphiphilic range, low TPSA, and
    presence of a basic nitrogen.

    Args:
        mol: RDKit molecule to assess.

    Returns:
        Dict with keys:
          - herg_risk (str): "low", "moderate", or "high"
          - risk_score (int): 0-4
          - max_score (int): always 4
          - flags (list[str]): human-readable descriptions of flagged factors
          - descriptors (dict): logp, mw, tpsa, has_basic_nitrogen values
    """
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    has_basic_n = mol.HasSubstructMatch(_get_basic_n_pattern())

    risk_score = 0
    flags: list[str] = []

    if logp > 3.7:
        risk_score += 1
        flags.append(f"LogP={logp:.2f} > 3.7")

    if 250 < mw < 500:
        risk_score += 1
        flags.append(f"MW={mw:.0f} in amphiphilic range 250-500")

    if tpsa < 75:
        risk_score += 1
        flags.append(f"TPSA={tpsa:.0f} < 75 (lipophilic)")

    if has_basic_n:
        risk_score += 1
        flags.append("Contains basic nitrogen")

    if risk_score <= 1:
        herg_risk = "low"
    elif risk_score == 2:
        herg_risk = "moderate"
    else:
        herg_risk = "high"

    return {
        "herg_risk": herg_risk,
        "risk_score": risk_score,
        "max_score": 4,
        "flags": flags,
        "descriptors": {
            "logp": logp,
            "mw": mw,
            "tpsa": tpsa,
            "has_basic_nitrogen": has_basic_n,
        },
    }
