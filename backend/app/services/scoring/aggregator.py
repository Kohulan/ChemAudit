"""
Aggregator Likelihood Scorer

Predicts the likelihood of a compound forming colloidal aggregates in assays.
Aggregators can cause false positives in HTS campaigns.

Based on:
- Shoichet Lab aggregator advisor (aggregator.orchard.ucsf.edu)
- Irwin et al. (2015) - "An Aggregation Advisor for Ligand Discovery"
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


@dataclass
class AggregatorLikelihoodResult:
    """Aggregator likelihood prediction result."""

    likelihood: str  # "low", "moderate", "high"
    risk_score: float  # 0-1 scale
    logp: float
    tpsa: float
    mw: float
    aromatic_rings: int
    risk_factors: List[str] = field(default_factory=list)
    interpretation: str = ""
    confidence: float = 0.0  # 0-1 scale
    evidence: List[Dict[str, Any]] = field(default_factory=list)


# Known aggregator SMILES patterns (representative examples)
# These are simplified SMARTS patterns for common aggregator scaffolds
AGGREGATOR_PATTERNS = [
    # Rhodanines
    ("[#6]1=[#6]C(=O)N([#6,#1])C1=S", "Rhodanine"),
    # Quinones
    ("O=C1C=CC(=O)C=C1", "Quinone"),
    # Catechols (can aggregate at high concentrations)
    ("c1cc(O)c(O)cc1", "Catechol"),
    # Curcumin-like
    ("CC(=O)C=Cc1ccc(O)c(OC)c1", "Curcumin-like"),
    # Phenol-sulfonamides with extended conjugation
    ("c1ccc(S(=O)(=O)N)cc1", "Sulfonamide"),
    # Flavonoids (some can aggregate)
    ("O=C1CC(c2ccccc2)Oc2ccccc12", "Flavone scaffold"),
    # Long chain fatty acids/lipids
    ("CCCCCCCCCCCC", "Long aliphatic chain"),
    # Acridines
    ("c1ccc2nc3ccccc3cc2c1", "Acridine"),
    # Phenothiazines
    ("c1ccc2c(c1)Nc1ccccc1S2", "Phenothiazine"),
    # Extended conjugated polyenes
    ("C=CC=CC=CC=C", "Extended conjugation"),
    # Azo compounds
    ("N=Nc1ccccc1", "Azo dye"),
]


class AggregatorScorer:
    """
    Predicts likelihood of compound aggregation.

    Uses multiple indicators:
    1. Lipophilicity (LogP) - high LogP increases aggregation risk
    2. Molecular size - larger molecules may aggregate more
    3. Aromatic content - highly aromatic compounds may stack
    4. Structural patterns - known aggregator scaffolds
    5. TPSA - low TPSA (high hydrophobicity) increases risk
    """

    # Thresholds based on literature
    LOGP_HIGH_RISK = 4.0  # LogP > 4 increases aggregation risk
    LOGP_MODERATE_RISK = 3.0
    TPSA_LOW_RISK = 40  # TPSA < 40 increases risk
    MW_HIGH = 500  # Large molecules may aggregate
    AROMATIC_RINGS_HIGH = 4  # Many aromatic rings increase stacking

    def __init__(self):
        """Initialize aggregator patterns."""
        self._patterns = []
        for smarts, name in AGGREGATOR_PATTERNS:
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    self._patterns.append((pattern, name))
            except Exception:
                pass

    # Total number of independent indicators checked
    TOTAL_INDICATORS = 6

    def score(self, mol: Chem.Mol) -> AggregatorLikelihoodResult:
        """
        Predict aggregator likelihood for a molecule.

        Args:
            mol: RDKit molecule object

        Returns:
            AggregatorLikelihoodResult with prediction, confidence, and evidence
        """
        risk_factors = []
        risk_score = 0.0
        evidence: List[Dict[str, Any]] = []
        triggered_count = 0

        # Calculate descriptors
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        mw = Descriptors.MolWt(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)

        # LogP risk assessment
        logp_triggered = logp > self.LOGP_MODERATE_RISK
        evidence.append({
            "indicator": "lipophilicity",
            "value": f"{logp:.2f}",
            "threshold": f">{self.LOGP_MODERATE_RISK}",
            "triggered": logp_triggered,
        })
        if logp > self.LOGP_HIGH_RISK:
            risk_score += 0.3
            risk_factors.append(f"High lipophilicity (LogP={logp:.1f}>4)")
            triggered_count += 1
        elif logp > self.LOGP_MODERATE_RISK:
            risk_score += 0.15
            risk_factors.append(f"Moderate lipophilicity (LogP={logp:.1f}>3)")
            triggered_count += 1

        # TPSA risk assessment (low TPSA = high hydrophobicity)
        tpsa_triggered = tpsa < self.TPSA_LOW_RISK
        evidence.append({
            "indicator": "low_polarity",
            "value": f"{tpsa:.1f}",
            "threshold": f"<{self.TPSA_LOW_RISK}",
            "triggered": tpsa_triggered,
        })
        if tpsa_triggered:
            risk_score += 0.2
            risk_factors.append(f"Low polarity (TPSA={tpsa:.0f}<40)")
            triggered_count += 1

        # Aromatic ring stacking risk
        arom_triggered = aromatic_rings >= 3
        evidence.append({
            "indicator": "aromatic_stacking",
            "value": str(aromatic_rings),
            "threshold": f">={self.AROMATIC_RINGS_HIGH}",
            "triggered": arom_triggered,
        })
        if aromatic_rings >= self.AROMATIC_RINGS_HIGH:
            risk_score += 0.2
            risk_factors.append(f"Many aromatic rings ({aromatic_rings})")
            triggered_count += 1
        elif aromatic_rings >= 3:
            risk_score += 0.1
            risk_factors.append(f"Multiple aromatic rings ({aromatic_rings})")
            triggered_count += 1

        # Molecular size
        mw_triggered = mw > self.MW_HIGH
        evidence.append({
            "indicator": "molecular_size",
            "value": f"{mw:.1f}",
            "threshold": f">{self.MW_HIGH}",
            "triggered": mw_triggered,
        })
        if mw_triggered:
            risk_score += 0.1
            risk_factors.append(f"Large molecule (MW={mw:.0f}>500)")
            triggered_count += 1

        # Check for known aggregator scaffolds
        matched_patterns = self._check_aggregator_patterns(mol)
        pattern_triggered = len(matched_patterns) > 0
        evidence.append({
            "indicator": "aggregator_scaffolds",
            "value": ", ".join(matched_patterns) if matched_patterns else "none",
            "threshold": "any known scaffold",
            "triggered": pattern_triggered,
        })
        if matched_patterns:
            risk_score += 0.2 * min(len(matched_patterns), 2)
            for pattern_name in matched_patterns[:3]:
                risk_factors.append(f"Contains {pattern_name} scaffold")
            triggered_count += 1

        # Calculate heavy atom count and check for highly conjugated systems
        num_atoms = mol.GetNumHeavyAtoms()
        num_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_fraction = num_aromatic_atoms / num_atoms if num_atoms > 0 else 0

        conj_triggered = aromatic_fraction > 0.7 and num_atoms > 20
        evidence.append({
            "indicator": "high_conjugation",
            "value": f"{aromatic_fraction:.2f}",
            "threshold": ">0.7 (with >20 atoms)",
            "triggered": conj_triggered,
        })
        if conj_triggered:
            risk_score += 0.15
            risk_factors.append("Highly conjugated aromatic system")
            triggered_count += 1

        # Cap risk score at 1.0
        risk_score = min(1.0, risk_score)

        # Confidence = fraction of indicators that triggered
        confidence = triggered_count / self.TOTAL_INDICATORS

        # Determine likelihood category
        if risk_score >= 0.6:
            likelihood = "high"
            interpretation = (
                "High aggregation risk. This compound has multiple features "
                "associated with colloidal aggregation. Consider counter-screening "
                "with detergent or dynamic light scattering."
            )
        elif risk_score >= 0.3:
            likelihood = "moderate"
            interpretation = (
                "Moderate aggregation risk. Some features suggest possible "
                "aggregation. Validate hits with orthogonal assays."
            )
        else:
            likelihood = "low"
            interpretation = (
                "Low aggregation risk. Compound properties are within ranges "
                "typically associated with non-aggregating compounds."
            )

        return AggregatorLikelihoodResult(
            likelihood=likelihood,
            risk_score=round(risk_score, 2),
            logp=round(logp, 2),
            tpsa=round(tpsa, 2),
            mw=round(mw, 2),
            aromatic_rings=aromatic_rings,
            risk_factors=risk_factors,
            interpretation=interpretation,
            confidence=round(confidence, 2),
            evidence=evidence,
        )

    def _check_aggregator_patterns(self, mol: Chem.Mol) -> List[str]:
        """Check for known aggregator scaffold patterns."""
        return [
            name for pattern, name in self._patterns if mol.HasSubstructMatch(pattern)
        ]


# Module-level instance
_scorer = AggregatorScorer()


def calculate_aggregator_likelihood(mol: Chem.Mol) -> AggregatorLikelihoodResult:
    """
    Predict aggregator likelihood for a molecule.

    Args:
        mol: RDKit molecule object

    Returns:
        AggregatorLikelihoodResult with prediction
    """
    return _scorer.score(mol)
