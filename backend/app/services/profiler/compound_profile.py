"""
Compound Profiler — Core Property Computation

Implements 6 profiling metrics:
- PFI (Property Forecast Index) — Young et al., Drug Discov Today 2011
- #Stars (QikProp-equivalent 95th-percentile outlier count)
- Abbott Bioavailability Score — Martin, J Med Chem 2005
- Consensus LogP (Wildman-Crippen + XLOGP3 approximation)
- Skin Permeation (Potts-Guy log Kp, 1992)
- 3D Shape Descriptors (PMI/NPR/PBF) with graceful failure

All functions are module-level (not class-based) following the existing
scoring service pattern. No SYBA or SCScore imports — those belong in
sa_comparison.py (Plan 02).
"""

import os
import sys
from typing import Optional

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# 95th-percentile ranges for 8 properties (QikProp-equivalent, Oprea 2000 / AstraZeneca data)
DRUG_95TH_RANGES: dict[str, tuple[float, float]] = {
    "MW": (130, 725),
    "LogP": (-2, 6.5),
    "HBA": (2, 20),
    "HBD": (0, 6),
    "TPSA": (7, 200),
    "RotBonds": (0, 15),
    "Fsp3": (0, 1),
    "NumRings": (0, 6),
}


# ---------------------------------------------------------------------------
# PFI — PROF-01
# ---------------------------------------------------------------------------


def compute_pfi(mol: Chem.Mol) -> dict:
    """
    Compute the Property Forecast Index (PFI).

    PFI = cLogP + Number of aromatic rings.
    Thresholds: <5 low, 5-7 moderate, >=7 high (Young et al., 2011).

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict with keys: pfi, clogp, aromatic_rings, risk.
    """
    clogp = Descriptors.MolLogP(mol)
    aromatic_rings = Descriptors.NumAromaticRings(mol)
    pfi = clogp + aromatic_rings

    if pfi < 5:
        risk = "low"
    elif pfi < 7:
        risk = "moderate"
    else:
        risk = "high"

    return {
        "pfi": round(pfi, 2),
        "clogp": round(clogp, 2),
        "aromatic_rings": int(aromatic_rings),
        "risk": risk,
    }


# ---------------------------------------------------------------------------
# Stars — PROF-02
# ---------------------------------------------------------------------------


def compute_stars(mol: Chem.Mol) -> dict:
    """
    Compute the #stars count (QikProp-equivalent 95th-percentile outlier count).

    Checks 8 molecular properties against the 95th-percentile ranges of oral drugs.
    Returns all 8 properties in the details list regardless of violation status.

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict with keys: stars (int, count of violations), details (list of 8 property dicts).
        Each detail dict has: property, value, range_low, range_high, violated.
    """
    props = {
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "HBA": int(Lipinski.NumHAcceptors(mol)),
        "HBD": int(Lipinski.NumHDonors(mol)),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "RotBonds": int(Lipinski.NumRotatableBonds(mol)),
        "Fsp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "NumRings": int(rdMolDescriptors.CalcNumRings(mol)),
    }

    details = []
    stars_count = 0

    for prop_name, (range_low, range_high) in DRUG_95TH_RANGES.items():
        value = props[prop_name]
        violated = not (range_low <= value <= range_high)
        if violated:
            stars_count += 1
        details.append(
            {
                "property": prop_name,
                "value": value,
                "range_low": range_low,
                "range_high": range_high,
                "violated": violated,
            }
        )

    return {
        "stars": stars_count,
        "details": details,
    }


# ---------------------------------------------------------------------------
# Abbott Bioavailability Score — PROF-03
# ---------------------------------------------------------------------------


def compute_abbott_score(mol: Chem.Mol) -> dict:
    """
    Compute the Abbott Bioavailability Score (4-class probability).

    Based on Martin, J Med Chem 2005. Uses TPSA, Lipinski violations, and
    formal charge to assign one of four oral bioavailability probabilities:
    0.11 (11%), 0.17 (17%), 0.56 (56%), 0.85 (85%).

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict with keys: abbott_score, probability_pct, tpsa, lipinski_violations.
    """
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    tpsa = Descriptors.TPSA(mol)

    # Count Lipinski violations
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1

    # Formal charge — approximate from GetFormalCharge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

    # 4-class assignment per Martin 2005
    if tpsa > 150:
        score = 0.11
    elif violations > 1:
        score = 0.17
    elif total_charge == 0 and violations == 0:
        score = 0.85
    else:
        score = 0.56

    return {
        "abbott_score": round(score, 2),
        "probability_pct": int(round(score * 100)),
        "tpsa": round(tpsa, 1),
        "lipinski_violations": violations,
    }


# ---------------------------------------------------------------------------
# Consensus LogP — PROF-04
# ---------------------------------------------------------------------------


def compute_consensus_logp(mol: Chem.Mol) -> dict:
    """
    Compute consensus LogP from two methods.

    Uses Wildman-Crippen (RDKit built-in) and an XLOGP3 approximation
    (0.92 * Wildman-Crippen + 0.12, per empirical correction).

    The xlogp3_is_approximation flag is always True — per Pitfall 7,
    this value is an approximation and must be disclosed to users.

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict with keys: consensus_logp, wildman_crippen, xlogp3_approx,
        xlogp3_is_approximation.
    """
    wc = Descriptors.MolLogP(mol)
    xlogp3_approx = 0.92 * wc + 0.12
    consensus = (wc + xlogp3_approx) / 2

    return {
        "consensus_logp": round(consensus, 2),
        "wildman_crippen": round(wc, 2),
        "xlogp3_approx": round(xlogp3_approx, 2),
        "xlogp3_is_approximation": True,
    }


# ---------------------------------------------------------------------------
# Skin Permeation — PROF-05
# ---------------------------------------------------------------------------


def compute_skin_permeation(mol: Chem.Mol) -> dict:
    """
    Compute skin permeation coefficient (log Kp).

    Uses the Potts-Guy 1992 equation:
        log Kp = -2.71 + 0.71 * logP - 0.0061 * MW

    Classification: <-6 low, -6 to -4 moderate, >=-4 high.

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict with keys: log_kp, classification.
    """
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    log_kp = -2.71 + 0.71 * logp - 0.0061 * mw

    if log_kp < -6:
        classification = "low"
    elif log_kp < -4:
        classification = "moderate"
    else:
        classification = "high"

    return {
        "log_kp": round(log_kp, 2),
        "classification": classification,
    }


# ---------------------------------------------------------------------------
# 3D Shape Descriptors — PROF-06
# ---------------------------------------------------------------------------


def _classify_shape(npr1: float, npr2: float) -> str:
    """
    Classify molecular shape based on NPR values.

    Mapping per Sauer & Schwarz (2003):
    - sphere: npr1 < 0.3 and npr2 > 0.7
    - rod:    npr1 > 0.5 and npr2 < 0.6
    - disc:   otherwise
    """
    if npr1 < 0.3 and npr2 > 0.7:
        return "sphere"
    if npr1 > 0.5 and npr2 < 0.6:
        return "rod"
    return "disc"


def compute_3d_shape(mol: Chem.Mol) -> Optional[dict]:
    """
    Compute 3D shape descriptors using ETKDGv3 conformer generation.

    Uses a 2-attempt cascade:
    1. ETKDGv3 (best quality)
    2. ETKDG (fallback)

    Applies MMFF94 force field optimization after successful embedding.
    Returns None on all failures — NEVER raises an exception.

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict with keys: pmi1, pmi2, pmi3, npr1, npr2, pbf, shape_class
        or None if conformer generation fails.
    """
    try:
        from rdkit.Chem import Descriptors3D

        mol3d = Chem.AddHs(mol)

        # Attempt 1: ETKDGv3
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        embed_result = AllChem.EmbedMolecule(mol3d, params)

        if embed_result == -1:
            # Attempt 2: ETKDG fallback
            mol3d = Chem.AddHs(mol)
            params2 = AllChem.EmbedParameters()
            params2.randomSeed = 42
            embed_result = AllChem.EmbedMolecule(mol3d, AllChem.ETKDGv3())
            if embed_result == -1:
                return None

        # MMFF optimization
        try:
            AllChem.MMFFOptimizeMolecule(mol3d, maxIters=200)
        except Exception:
            pass  # Continue with unoptimized conformer

        # Compute 3D descriptors
        pmi1 = Descriptors3D.PMI1(mol3d)
        pmi2 = Descriptors3D.PMI2(mol3d)
        pmi3 = Descriptors3D.PMI3(mol3d)
        npr1 = Descriptors3D.NPR1(mol3d)
        npr2 = Descriptors3D.NPR2(mol3d)
        pbf = Descriptors3D.PBF(mol3d)

        return {
            "pmi1": round(pmi1, 4),
            "pmi2": round(pmi2, 4),
            "pmi3": round(pmi3, 4),
            "npr1": round(npr1, 4),
            "npr2": round(npr2, 4),
            "pbf": round(pbf, 4),
            "shape_class": _classify_shape(npr1, npr2),
        }

    except Exception:
        return None


# ---------------------------------------------------------------------------
# Full Profile (composite, lazy 3D) — D-26
# ---------------------------------------------------------------------------


def compute_full_profile(mol: Chem.Mol) -> dict:
    """
    Compute a composite property profile for a molecule.

    Calls all core profiling functions EXCEPT 3D shape (lazy, per D-26).
    Also includes the drug-likeness rules grid from the existing scoring layer.

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict containing: pfi, stars, abbott, consensus_logp, skin_permeation,
        and druglikeness (from existing calculate_druglikeness).
    """
    from app.services.scoring.druglikeness import calculate_druglikeness

    pfi_result = compute_pfi(mol)
    stars_result = compute_stars(mol)
    abbott_result = compute_abbott_score(mol)
    consensus_logp_result = compute_consensus_logp(mol)
    skin_result = compute_skin_permeation(mol)

    # Include existing drug-likeness rules (Lipinski, QED, Veber, etc.) per D-03
    dl_result = calculate_druglikeness(mol)

    return {
        "pfi": pfi_result,
        "stars": stars_result,
        "abbott": abbott_result,
        "consensus_logp": consensus_logp_result,
        "skin_permeation": skin_result,
        "druglikeness": {
            "lipinski": {
                "passed": dl_result.lipinski.passed,
                "violations": dl_result.lipinski.violations,
            },
            "qed": round(dl_result.qed.score, 3),
            "veber": {"passed": dl_result.veber.passed},
        },
    }
