"""
FilterConfig dataclass and preset configurations for the structure filter pipeline.

Defines property thresholds, alert catalog selection, novelty settings, and
composite scorer weight vectors for 4 preset configurations per D-12 and D-15.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class FilterConfig:
    """
    Configuration for the structure filter pipeline.

    Specifies property thresholds (D-12), alert catalog selection (D-12),
    novelty settings (D-18, D-19), and composite scorer weight vectors (D-15).
    """

    # Property thresholds (per D-12)
    min_mw: float = 200.0
    max_mw: float = 500.0
    min_logp: float = -1.0
    max_logp: float = 5.0
    max_tpsa: float = 140.0
    max_rot_bonds: int = 10
    max_rings: int | None = None  # None = no ring limit
    max_sa_score: float = 5.0
    # Alert catalog selection (per D-12)
    use_pains: bool = True
    use_brenk: bool = True
    use_kazius: bool = True
    use_nibr: bool = False
    # Novelty (per D-18, D-19)
    enable_novelty: bool = False
    novelty_threshold: float = 0.85
    # Composite scorer weights (per D-15)
    weight_validity: float = 0.3
    weight_qed: float = 0.3
    weight_alert_free: float = 0.2
    weight_sa: float = 0.2


# Per D-15: each preset has its own weight vector.
# - drug_like: balanced across all 4 components (default from feature guide)
# - lead_like: higher QED weight (lead optimization cares more about drug-likeness)
# - fragment_like: higher alert_free weight (fragments must be clean building blocks)
# - permissive: higher validity weight, lower alert_free (alerts disabled anyway)
PRESETS: dict[str, FilterConfig] = {
    "drug_like": FilterConfig(
        # D-12 thresholds: MW 200-500, LogP -1 to 5, TPSA <=140, RotBonds <=10, SA <=5.0
        # Alerts: PAINS + Brenk + Kazius
        # D-15 weights: balanced default
        weight_validity=0.3,
        weight_qed=0.3,
        weight_alert_free=0.2,
        weight_sa=0.2,
    ),
    "lead_like": FilterConfig(
        max_mw=350.0,
        max_logp=3.5,
        max_rot_bonds=7,
        max_sa_score=4.0,
        use_nibr=True,
        # D-15 weights: emphasize QED for lead optimization
        weight_validity=0.2,
        weight_qed=0.4,
        weight_alert_free=0.2,
        weight_sa=0.2,
    ),
    "fragment_like": FilterConfig(
        min_mw=100.0,
        max_mw=300.0,
        max_logp=3.0,
        max_tpsa=100.0,
        max_rot_bonds=3,
        max_rings=3,
        max_sa_score=3.0,
        use_brenk=False,
        use_kazius=False,
        # D-15 weights: emphasize alert-free (clean building blocks) and SA (easy to synthesize)
        weight_validity=0.2,
        weight_qed=0.2,
        weight_alert_free=0.3,
        weight_sa=0.3,
    ),
    "permissive": FilterConfig(
        min_mw=100.0,
        max_mw=800.0,
        min_logp=-5.0,
        max_logp=8.0,
        max_tpsa=200.0,
        max_rot_bonds=15,
        max_sa_score=7.0,
        use_pains=False,
        use_brenk=False,
        use_kazius=False,
        # D-15 weights: emphasize validity (main filter), lower alert_free (alerts off anyway)
        weight_validity=0.4,
        weight_qed=0.3,
        weight_alert_free=0.1,
        weight_sa=0.2,
    ),
}
