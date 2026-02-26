"""
Scoring Services

Provides comprehensive molecular scoring including:
- ML-readiness scoring
- NP-likeness scoring (with fragment breakdown)
- Scaffold extraction
- Drug-likeness scoring (Lipinski, QED, Veber, Ro3, Ghose, Egan, Muegge)
- Consensus drug-likeness scoring (5 rule sets)
- Lead-likeness assessment
- Safety filters (PAINS, Brenk, NIH, ZINC, ChEMBL)
- ADMET predictions (SAscore, ESOL, Fsp3, CNS MPO, Pfizer/GSK rules)
- Aggregator likelihood prediction
- Salt/counterion inventory
- Ligand efficiency
- Property breakdown (TPSA, LogP per-atom, Bertz detail, Fsp3 detail)
- Bioavailability radar (6-axis normalized assessment)
- BOILED-Egg classification (GI/BBB permeation)
- Property radar comparison
"""

from app.services.scoring.admet import ADMETScorer, calculate_admet
from app.services.scoring.aggregator import (
    AggregatorScorer,
    calculate_aggregator_likelihood,
)
from app.services.scoring.bioavailability_radar import (
    BioavailabilityRadarScorer,
    calculate_bioavailability_radar,
    calculate_boiled_egg,
    calculate_radar_comparison,
)
from app.services.scoring.druglikeness import (
    DrugLikenessScorer,
    calculate_consensus,
    calculate_druglikeness,
    calculate_lead_likeness,
)
from app.services.scoring.ml_readiness import MLReadinessScorer, calculate_ml_readiness
from app.services.scoring.np_likeness import (
    NPLikenessScorer,
    calculate_np_breakdown,
    calculate_np_likeness,
)
from app.services.scoring.property_breakdown import (
    PropertyBreakdownScorer,
    calculate_bertz_detail,
    calculate_fsp3_detail,
    calculate_logp_breakdown,
    calculate_tpsa_breakdown,
)
from app.services.scoring.safety_filters import (
    SafetyFilterScorer,
    calculate_safety_filters,
    get_pains_alerts,
    is_pains_clean,
)
from app.services.scoring.salt_inventory import (
    LigandEfficiencyScorer,
    SaltInventoryScorer,
    calculate_ligand_efficiency,
    calculate_salt_inventory,
)
from app.services.scoring.scaffold import extract_scaffold

__all__ = [
    # ML Readiness
    "MLReadinessScorer",
    "calculate_ml_readiness",
    # NP Likeness
    "NPLikenessScorer",
    "calculate_np_likeness",
    "calculate_np_breakdown",
    # Scaffold
    "extract_scaffold",
    # Drug-likeness
    "DrugLikenessScorer",
    "calculate_druglikeness",
    "calculate_consensus",
    "calculate_lead_likeness",
    # Safety Filters
    "SafetyFilterScorer",
    "calculate_safety_filters",
    "get_pains_alerts",
    "is_pains_clean",
    # ADMET
    "ADMETScorer",
    "calculate_admet",
    # Aggregator
    "AggregatorScorer",
    "calculate_aggregator_likelihood",
    # Salt Inventory
    "SaltInventoryScorer",
    "calculate_salt_inventory",
    # Ligand Efficiency
    "LigandEfficiencyScorer",
    "calculate_ligand_efficiency",
    # Property Breakdown
    "PropertyBreakdownScorer",
    "calculate_tpsa_breakdown",
    "calculate_logp_breakdown",
    "calculate_bertz_detail",
    "calculate_fsp3_detail",
    # Bioavailability Radar
    "BioavailabilityRadarScorer",
    "calculate_bioavailability_radar",
    "calculate_boiled_egg",
    "calculate_radar_comparison",
]
