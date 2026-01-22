"""
ML-Readiness Scorer

Calculates how suitable a molecule is for machine learning applications.
Score breakdown:
- Descriptors (40 points): Successful descriptor calculations
- Fingerprints (40 points): Successful fingerprint generation
- Size/Elements (20 points): Molecular weight and atom count constraints
"""
from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, MACCSkeys, rdMolDescriptors


@dataclass
class MLReadinessBreakdown:
    """Breakdown of ML-readiness score components."""

    descriptors_score: float = 0.0
    descriptors_max: float = 40.0
    descriptors_successful: int = 0
    descriptors_total: int = 0
    failed_descriptors: list = field(default_factory=list)

    fingerprints_score: float = 0.0
    fingerprints_max: float = 40.0
    fingerprints_successful: list = field(default_factory=list)
    fingerprints_failed: list = field(default_factory=list)

    size_score: float = 0.0
    size_max: float = 20.0
    molecular_weight: Optional[float] = None
    num_atoms: Optional[int] = None
    size_category: str = "unknown"


@dataclass
class MLReadinessResult:
    """Result of ML-readiness scoring."""

    score: int
    breakdown: MLReadinessBreakdown
    interpretation: str
    failed_descriptors: list


class MLReadinessScorer:
    """
    Scores molecules for ML-readiness based on descriptor calculability,
    fingerprint generation success, and size constraints.
    """

    # Fingerprint types with their point allocations
    FINGERPRINT_TYPES = [
        ("morgan", 15),
        ("maccs", 15),
        ("atompair", 10),
    ]

    # Size thresholds
    OPTIMAL_MW_RANGE = (100, 900)
    OPTIMAL_ATOM_RANGE = (3, 100)
    ACCEPTABLE_MW_RANGE = (50, 1200)
    ACCEPTABLE_ATOM_RANGE = (1, 150)

    def score(self, mol: Chem.Mol) -> MLReadinessResult:
        """
        Calculate ML-readiness score for a molecule.

        Args:
            mol: RDKit molecule object

        Returns:
            MLReadinessResult with score, breakdown, and interpretation
        """
        breakdown = MLReadinessBreakdown()

        # Calculate descriptor score (40 points max)
        self._score_descriptors(mol, breakdown)

        # Calculate fingerprint score (40 points max)
        self._score_fingerprints(mol, breakdown)

        # Calculate size score (20 points max)
        self._score_size(mol, breakdown)

        # Calculate total score
        total_score = int(
            breakdown.descriptors_score +
            breakdown.fingerprints_score +
            breakdown.size_score
        )
        total_score = max(0, min(100, total_score))

        # Generate interpretation
        interpretation = self._get_interpretation(total_score, breakdown)

        return MLReadinessResult(
            score=total_score,
            breakdown=breakdown,
            interpretation=interpretation,
            failed_descriptors=breakdown.failed_descriptors,
        )

    def _score_descriptors(
        self, mol: Chem.Mol, breakdown: MLReadinessBreakdown
    ) -> None:
        """Score descriptor calculability (40 points max)."""
        try:
            # Use CalcMolDescriptors with missingVal=None to track failures
            # The silent=True prevents RDKit warnings from cluttering output
            descriptors = Descriptors.CalcMolDescriptors(
                mol, missingVal=None, silent=True
            )

            total = len(descriptors)
            successful = 0
            failed = []

            for name, value in descriptors.items():
                if value is not None:
                    successful += 1
                else:
                    failed.append(name)

            breakdown.descriptors_total = total
            breakdown.descriptors_successful = successful
            breakdown.failed_descriptors = failed

            if total > 0:
                breakdown.descriptors_score = (
                    breakdown.descriptors_max * (successful / total)
                )

        except Exception as e:
            # If descriptor calculation completely fails
            breakdown.descriptors_score = 0.0
            breakdown.failed_descriptors = [f"CalcMolDescriptors error: {str(e)}"]

    def _score_fingerprints(
        self, mol: Chem.Mol, breakdown: MLReadinessBreakdown
    ) -> None:
        """Score fingerprint generation success (40 points max)."""
        total_score = 0.0
        successful = []
        failed = []

        for fp_name, points in self.FINGERPRINT_TYPES:
            try:
                if fp_name == "morgan":
                    # Morgan fingerprint (radius=2, 2048 bits)
                    fp = AllChem.GetMorganFingerprintAsBitVect(
                        mol, radius=2, nBits=2048
                    )
                elif fp_name == "maccs":
                    # MACCS keys
                    fp = MACCSkeys.GenMACCSKeys(mol)
                elif fp_name == "atompair":
                    # Atom pair fingerprint
                    fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                        mol, nBits=2048
                    )
                else:
                    continue

                # If we get here without exception, fingerprint was successful
                if fp is not None:
                    total_score += points
                    successful.append(fp_name)
                else:
                    failed.append(fp_name)

            except Exception:
                failed.append(fp_name)

        breakdown.fingerprints_score = total_score
        breakdown.fingerprints_successful = successful
        breakdown.fingerprints_failed = failed

    def _score_size(self, mol: Chem.Mol, breakdown: MLReadinessBreakdown) -> None:
        """Score molecule size constraints (20 points max)."""
        try:
            mw = Descriptors.MolWt(mol)
            num_atoms = mol.GetNumAtoms()

            breakdown.molecular_weight = mw
            breakdown.num_atoms = num_atoms

            # Check optimal range
            mw_optimal = (
                self.OPTIMAL_MW_RANGE[0] <= mw <= self.OPTIMAL_MW_RANGE[1]
            )
            atoms_optimal = (
                self.OPTIMAL_ATOM_RANGE[0] <= num_atoms <= self.OPTIMAL_ATOM_RANGE[1]
            )

            # Check acceptable range
            mw_acceptable = (
                self.ACCEPTABLE_MW_RANGE[0] <= mw <= self.ACCEPTABLE_MW_RANGE[1]
            )
            atoms_acceptable = (
                self.ACCEPTABLE_ATOM_RANGE[0] <= num_atoms <= self.ACCEPTABLE_ATOM_RANGE[1]
            )

            if mw_optimal and atoms_optimal:
                breakdown.size_score = 20.0
                breakdown.size_category = "optimal"
            elif mw_acceptable and atoms_acceptable:
                breakdown.size_score = 10.0
                breakdown.size_category = "acceptable"
            else:
                breakdown.size_score = 0.0
                breakdown.size_category = "out_of_range"

        except Exception:
            breakdown.size_score = 0.0
            breakdown.size_category = "error"

    def _get_interpretation(
        self, score: int, breakdown: MLReadinessBreakdown
    ) -> str:
        """Generate human-readable interpretation of the score."""
        parts = []

        # Overall assessment
        if score >= 80:
            parts.append("Excellent ML-readiness.")
        elif score >= 60:
            parts.append("Good ML-readiness with minor limitations.")
        elif score >= 40:
            parts.append("Moderate ML-readiness; some calculations may fail.")
        else:
            parts.append("Poor ML-readiness; significant computation issues likely.")

        # Specific issues
        if breakdown.failed_descriptors:
            count = len(breakdown.failed_descriptors)
            parts.append(f"{count} descriptors could not be calculated.")

        if breakdown.fingerprints_failed:
            parts.append(
                f"Failed fingerprints: {', '.join(breakdown.fingerprints_failed)}."
            )

        if breakdown.size_category == "out_of_range":
            parts.append("Molecule size outside typical ML training ranges.")
        elif breakdown.size_category == "acceptable":
            parts.append("Molecule size is within acceptable but not optimal range.")

        return " ".join(parts)


# Module-level convenience function
_scorer = MLReadinessScorer()


def calculate_ml_readiness(mol: Chem.Mol) -> MLReadinessResult:
    """
    Calculate ML-readiness score for a molecule.

    Args:
        mol: RDKit molecule object

    Returns:
        MLReadinessResult with score (0-100), breakdown, and interpretation
    """
    return _scorer.score(mol)
