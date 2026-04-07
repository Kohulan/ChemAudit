"""
Composite scorer for generative chemistry (REINVENT-compatible).

Computes a 0-1 composite score for a SMILES string based on:
  - Validity (1.0 if parseable)
  - QED (drug-likeness, 0-1)
  - Alert-free binary (1.0 if no structural alerts, 0.0 otherwise)
  - SA Score normalized (1.0 when sa=1, 0.0 when sa=10)

Weights are preset-specific per D-15, enabling differentiated scoring
across drug_like, lead_like, fragment_like, and permissive configurations.
"""

from __future__ import annotations

import logging
import os
import sys
from typing import Optional

from rdkit import Chem, RDConfig
from rdkit.Chem import QED

from app.services.alerts.kazius_rules import screen_kazius
from app.services.alerts.nibr_filters import screen_nibr
from app.services.genchem.filter_config import FilterConfig
from app.services.scoring.safety_filters import _scorer

# SA Score via RDKit Contrib (Phase 7 pattern)
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore  # noqa: E402

logger = logging.getLogger(__name__)


def _has_any_alerts(mol: Chem.Mol, config: FilterConfig) -> bool:
    """
    Check if a molecule has any structural alerts according to the config flags.

    Mirrors the alert stage logic in filter_pipeline.py.

    Args:
        mol: RDKit molecule object.
        config: FilterConfig specifying which alert catalogs to use.

    Returns:
        True if any alert is found, False otherwise.
    """
    if config.use_pains:
        if _scorer.get_pains_alerts(mol):
            return True
    if config.use_brenk:
        if _scorer.get_brenk_alerts(mol):
            return True
    if config.use_kazius:
        if screen_kazius(mol):
            return True
    if config.use_nibr:
        if screen_nibr(mol):
            return True
    return False


def score_for_generative(smiles: str, config: FilterConfig) -> Optional[float]:
    """
    Compute a composite 0-1 score for use in generative model reward functions.

    Returns None if the SMILES is invalid or yields an empty molecule, ensuring
    the calling code can distinguish "invalid structure" from "low-scoring structure".

    Per D-14/Pitfall 4: always return None (not 0.0) for invalid SMILES so that
    REINVENT and similar tools can filter out non-parseable candidates cleanly.

    Args:
        smiles: SMILES string to score.
        config: FilterConfig with weight vector and alert catalog selection.

    Returns:
        Float 0.0-1.0 for valid SMILES, or None for invalid/empty structures.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return None

    # Validity component — always 1.0 since parsing succeeded
    validity = 1.0

    # QED — drug-likeness (0-1 range, higher = more drug-like)
    qed_score = QED.qed(mol)

    # Alert-free component — binary: 1.0 if no alerts, 0.0 if any match
    alert_free = 0.0 if _has_any_alerts(mol, config) else 1.0

    # SA Score normalized to [0, 1]: 1.0 when sa=1 (easy), 0.0 when sa=10 (hard)
    sa_raw = sascorer.calculateScore(mol)
    sa_normalized = max(0.0, (10.0 - sa_raw) / 9.0)

    # Composite score via preset weight vector (D-15)
    composite = (
        config.weight_validity * validity
        + config.weight_qed * qed_score
        + config.weight_alert_free * alert_free
        + config.weight_sa * sa_normalized
    )

    return round(composite, 4)
