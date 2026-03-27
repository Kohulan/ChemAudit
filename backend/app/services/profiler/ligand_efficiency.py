"""
Extended Ligand Efficiency Metrics

Computes 5 ligand efficiency metrics from measured activity data:
- LE   (Ligand Efficiency): 1.4 * pIC50 / HA
- LLE  (Lipophilic Ligand Efficiency): pIC50 - LogP
- LELP (Ligand Efficiency Lipophilicity Price): LogP / LE
- BEI  (Binding Efficiency Index): pIC50 * 1000 / MW
- SEI  (Surface Efficiency Index): pIC50 * 100 / TPSA

Supports activity types: IC50_nM, IC50_uM, Ki_nM, pIC50, pKd.

References:
- LE: Hopkins et al. Drug Discov Today 2004;9:430-431
- LLE: Leeson & Springthorpe Nat Rev Drug Discov 2007;6:881-890
- LELP: Keseru & Makara Drug Discov Today 2009;14:741-748
- BEI/SEI: Abad-Zapatero & Metz Drug Discov Today 2005;10:464-469
"""

import math

from rdkit import Chem
from rdkit.Chem import Descriptors


def compute_ligand_efficiency(
    mol: Chem.Mol,
    activity_value: float,
    activity_type: str,
) -> dict:
    """
    Compute extended ligand efficiency metrics from supplied activity data.

    Activity type conversions:
    - IC50_nM: pIC50 = -log10(value * 1e-9)
    - IC50_uM: pIC50 = -log10(value * 1e-6)
    - Ki_nM:   pIC50 = -log10(value * 1e-9)  (treated as equivalent)
    - pIC50:   pIC50 = value (used directly)
    - pKd:     pIC50 = value (used directly as equivalent)

    Args:
        mol: RDKit molecule object.
        activity_value: Numeric activity measurement.
        activity_type: One of IC50_nM, IC50_uM, Ki_nM, pIC50, pKd.

    Returns:
        Dict with keys: pIC50, LE, LLE, LELP, BEI, SEI.
        On unknown activity_type: {"error": "Unknown activity_type: {type}"}.
    """
    # Convert activity to pIC50
    if activity_type == "pIC50":
        pic50 = activity_value
    elif activity_type == "pKd":
        pic50 = activity_value
    elif activity_type in ("IC50_nM", "Ki_nM"):
        pic50 = -math.log10(activity_value * 1e-9)
    elif activity_type == "IC50_uM":
        pic50 = -math.log10(activity_value * 1e-6)
    else:
        return {"error": f"Unknown activity_type: {activity_type}"}

    ha = mol.GetNumHeavyAtoms()
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)

    # Core metrics
    le = 1.4 * pic50 / ha if ha > 0 else 0.0
    lle = pic50 - logp
    lelp = logp / le if le != 0 else float("inf")
    bei = pic50 * 1000 / mw if mw > 0 else 0.0
    sei = pic50 * 100 / tpsa if tpsa > 0 else 0.0

    return {
        "pIC50": round(pic50, 2),
        "LE": round(le, 3),
        "LLE": round(lle, 2),
        "LELP": round(lelp, 2),
        "BEI": round(bei, 2),
        "SEI": round(sei, 2),
    }
