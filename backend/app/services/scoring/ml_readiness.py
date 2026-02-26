"""
ML-Readiness Scorer

Calculates how suitable a molecule is for machine learning applications
using a scientifically meaningful 4-dimension scoring system (0-100):

  1. Structural Quality   (20 pts) — binary pass/fail structural checks
  2. Property Profile      (35 pts) — desirability-scored physicochemical properties
  3. Complexity & Feasibility (25 pts) — QED, SA Score, Fsp3, stereocenters
  4. Representation Quality (20 pts) — descriptor/fingerprint completeness + informativeness

Citations:
  - Bickerton et al. (2012) Nature Chemistry 4:90-98  (QED, desirability)
  - Lipinski et al. (1997) Adv. Drug Deliv. Rev. 23:3-25
  - Ertl & Schuffenhauer (2009) J. Cheminf. 1:8  (SA Score)
  - Lovering et al. (2009) J. Med. Chem. 52:6752-6756  (Fsp3)
  - Rogers & Hahn (2010) J. Chem. Inf. Model. 50:742-754  (Morgan FP)
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem, Descriptors, MACCSkeys, rdFingerprintGenerator, rdMolDescriptors

from app.services.scoring.profile_scoring import desirability

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class MLDimensionItem:
    """A single scored item within a dimension."""

    name: str
    score: float
    max_score: float
    passed: Optional[bool] = None
    subtitle: Optional[str] = None
    tooltip: Optional[str] = None


@dataclass
class MLDimension:
    """A scored dimension with its constituent items."""

    name: str
    score: float
    max_score: float
    items: List[MLDimensionItem] = field(default_factory=list)
    details: Dict = field(default_factory=dict)


@dataclass
class MLReadinessBreakdown:
    """Breakdown of the ML-readiness score into 4 dimensions."""

    structural_quality: MLDimension = field(default_factory=lambda: MLDimension(
        name="Structural Quality", score=0, max_score=20,
    ))
    property_profile: MLDimension = field(default_factory=lambda: MLDimension(
        name="Property Profile", score=0, max_score=35,
    ))
    complexity_feasibility: MLDimension = field(default_factory=lambda: MLDimension(
        name="Complexity & Feasibility", score=0, max_score=25,
    ))
    representation_quality: MLDimension = field(default_factory=lambda: MLDimension(
        name="Representation Quality", score=0, max_score=20,
    ))


@dataclass
class MLReadinessResult:
    """Result of ML-readiness scoring."""

    score: int
    label: str
    color: str
    breakdown: MLReadinessBreakdown
    caveats: List[str] = field(default_factory=list)
    supplementary: Dict = field(default_factory=dict)
    interpretation: str = ""


# ---------------------------------------------------------------------------
# Interpretation helpers
# ---------------------------------------------------------------------------

_THRESHOLDS = [
    (85, "Excellent", "green"),
    (70, "Good", "teal"),
    (50, "Moderate", "amber"),
    (30, "Limited", "orange"),
    (0, "Poor", "red"),
]


def _label_and_color(score: int) -> tuple:
    for threshold, label, color in _THRESHOLDS:
        if score >= threshold:
            return label, color
    return "Poor", "red"


# ---------------------------------------------------------------------------
# Scorer
# ---------------------------------------------------------------------------

class MLReadinessScorer:
    """Scores molecules for ML-readiness across 4 scientific dimensions."""

    # Fingerprint types with point allocations for Dimension 4
    FINGERPRINT_TYPES = [
        ("morgan", 8),
        ("morgan_features", 8),
        ("maccs", 8),
        ("atompair", 4),
        ("topological_torsion", 4),
        ("rdkit_fp", 4),
        ("avalon", 4),
    ]

    # Property profile ranges (Dimension 2)
    PROPERTY_RANGES = [
        ("mw", 200, 500, 6),
        ("logp", 0.5, 4.0, 6),
        ("tpsa", 40, 120, 5),
        ("hbd", 0, 3, 4),
        ("hba", 2, 8, 4),
        ("rotatable_bonds", 1, 8, 5),
        ("aromatic_rings", 1, 3, 5),
    ]

    def __init__(self):
        """Initialize fingerprint generators."""
        self._morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        self._morgan_feat_gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=2,
            fpSize=2048,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
        )
        self._atompair_gen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
        self._torsion_gen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
        self._rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)

    def score(self, mol: Chem.Mol) -> MLReadinessResult:
        """Calculate ML-readiness score for a molecule."""
        breakdown = MLReadinessBreakdown()
        caveats: List[str] = []

        self._score_structural_quality(mol, breakdown, caveats)
        self._score_property_profile(mol, breakdown)
        self._score_complexity_feasibility(mol, breakdown)
        self._score_representation_quality(mol, breakdown)

        total = int(
            breakdown.structural_quality.score
            + breakdown.property_profile.score
            + breakdown.complexity_feasibility.score
            + breakdown.representation_quality.score
        )
        total = max(0, min(100, total))

        label, color = _label_and_color(total)
        interpretation = self._get_interpretation(total, breakdown, caveats)

        # Supplementary (informational, not scored)
        supplementary = self._get_supplementary(mol)

        return MLReadinessResult(
            score=total,
            label=label,
            color=color,
            breakdown=breakdown,
            caveats=caveats,
            supplementary=supplementary,
            interpretation=interpretation,
        )

    # ------------------------------------------------------------------
    # Dimension 1: Structural Quality (20 pts)
    # ------------------------------------------------------------------

    def _score_structural_quality(
        self,
        mol: Chem.Mol,
        breakdown: MLReadinessBreakdown,
        caveats: List[str],
    ) -> None:
        """Binary checks for structural soundness."""
        from app.services.validation.checks.deep_complexity import (
            ChargedSpeciesCheck,
            PolymerDetectionCheck,
        )
        from app.services.validation.checks.deep_composition import (
            InorganicFilterCheck,
            IsotopeLabelDetectionCheck,
            MixtureDetectionCheck,
            RadicalDetectionCheck,
            TrivialMoleculeCheck,
        )

        dim = breakdown.structural_quality
        items: List[MLDimensionItem] = []

        # 1. Single component (5 pts)
        mixture_result = MixtureDetectionCheck().run(mol)
        passed = mixture_result.passed
        items.append(MLDimensionItem(
            "Single component", 5.0 if passed else 0.0, 5.0, passed,
            tooltip="Checks for multi-component mixtures (salts, solvates). "
            "Mixtures cause descriptor failures in ML pipelines.",
        ))
        if not passed:
            caveats.append("Multi-component mixture detected")

        # 2. Standard organic elements (5 pts)
        inorganic_result = InorganicFilterCheck().run(mol)
        passed = inorganic_result.passed
        items.append(MLDimensionItem(
            "Standard organic elements", 5.0 if passed else 0.0, 5.0, passed,
            tooltip="Checks for inorganic atoms or metals. Non-organic elements "
            "have poor descriptor coverage in standard ML featurizers.",
        ))
        if not passed:
            is_inorganic = inorganic_result.details.get("is_inorganic", False)
            if is_inorganic:
                caveats.append("Inorganic compound detected")
            else:
                caveats.append("Organometallic or non-standard elements detected")

        # 3. No radicals (3 pts)
        radical_result = RadicalDetectionCheck().run(mol)
        passed = radical_result.passed
        items.append(MLDimensionItem(
            "No radicals", 3.0 if passed else 0.0, 3.0, passed,
            tooltip="Checks for unpaired electrons. Radical species are poorly "
            "represented in most molecular descriptor sets.",
        ))
        if not passed:
            caveats.append("Radical electrons detected")

        # 4. Reasonable charge |net| <= 2 (3 pts)
        charge_result = ChargedSpeciesCheck().run(mol)
        net_charge = abs(charge_result.details.get("net_charge", 0))
        passed = net_charge <= 2
        items.append(MLDimensionItem(
            "Reasonable charge", 3.0 if passed else 0.0, 3.0, passed,
            tooltip="Checks net formal charge \u2264 \u00b12. Highly charged species "
            "behave anomalously in property predictions.",
        ))

        # 5. No dummy atoms (4 pts)
        polymer_result = PolymerDetectionCheck().run(mol)
        has_dummy = polymer_result.details.get("has_dummy_atoms", False)
        passed = not has_dummy
        items.append(MLDimensionItem(
            "No dummy atoms", 4.0 if passed else 0.0, 4.0, passed,
            tooltip="Checks for placeholder atoms (R-groups, *). Dummy atoms "
            "break fingerprint and descriptor calculations.",
        ))

        # Caveat-only checks (not scored)
        isotope_result = IsotopeLabelDetectionCheck().run(mol)
        if not isotope_result.passed:
            caveats.append("Isotope-labeled atoms present")

        trivial_result = TrivialMoleculeCheck().run(mol)
        if not trivial_result.passed:
            caveats.append("Trivially small molecule")

        dim.items = items
        dim.score = round(sum(i.score for i in items), 1)

    # ------------------------------------------------------------------
    # Dimension 2: Property Profile (35 pts)
    # ------------------------------------------------------------------

    def _score_property_profile(
        self, mol: Chem.Mol, breakdown: MLReadinessBreakdown
    ) -> None:
        """Desirability-scored physicochemical properties."""
        dim = breakdown.property_profile
        items: List[MLDimensionItem] = []

        # Calculate properties
        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        except Exception:
            dim.score = 0.0
            return

        prop_values = {
            "mw": mw,
            "logp": logp,
            "tpsa": tpsa,
            "hbd": float(hbd),
            "hba": float(hba),
            "rotatable_bonds": float(rot_bonds),
            "aromatic_rings": float(aromatic_rings),
        }

        prop_config = {
            "mw": {
                "label": "MW fitness",
                "tooltip": "Desirability score for molecular weight. Ideal range "
                "based on drug-like chemical space (Bickerton 2012).",
                "fmt": lambda v, lo, hi: f"{v:.1f} Da (ideal: {lo}\u2013{hi})",
            },
            "logp": {
                "label": "LogP fitness",
                "tooltip": "Desirability score for partition coefficient. "
                "Measures lipophilicity suitability for ML training sets.",
                "fmt": lambda v, lo, hi: f"{v:.2f} (ideal: {lo}\u2013{hi})",
            },
            "tpsa": {
                "label": "TPSA fitness",
                "tooltip": "Desirability score for topological polar surface area. "
                "Relates to membrane permeability.",
                "fmt": lambda v, lo, hi: f"{v:.1f} \u00c5\u00b2 (ideal: {lo}\u2013{hi})",
            },
            "hbd": {
                "label": "HBD fitness",
                "tooltip": "Desirability score for hydrogen bond donor count.",
                "fmt": lambda v, lo, hi: f"{int(v)} donors (ideal: {lo}\u2013{hi})",
            },
            "hba": {
                "label": "HBA fitness",
                "tooltip": "Desirability score for hydrogen bond acceptor count.",
                "fmt": lambda v, lo, hi: f"{int(v)} acceptors (ideal: {lo}\u2013{hi})",
            },
            "rotatable_bonds": {
                "label": "RotBonds fitness",
                "tooltip": "Desirability score for rotatable bond count. "
                "Relates to conformational flexibility.",
                "fmt": lambda v, lo, hi: f"{int(v)} bonds (ideal: {lo}\u2013{hi})",
            },
            "aromatic_rings": {
                "label": "AromRings fitness",
                "tooltip": "Desirability score for aromatic ring count. "
                "Relates to solubility and planarity.",
                "fmt": lambda v, lo, hi: f"{int(v)} rings (ideal: {lo}\u2013{hi})",
            },
        }

        for prop_name, min_val, max_val, max_pts in self.PROPERTY_RANGES:
            value = prop_values[prop_name]
            d = desirability(value, min_val, max_val) or 0.0
            pts = round(d * max_pts, 2)
            cfg = prop_config[prop_name]
            items.append(MLDimensionItem(
                cfg["label"], pts, float(max_pts),
                subtitle=cfg["fmt"](value, min_val, max_val),
                tooltip=cfg["tooltip"],
            ))

        dim.items = items
        dim.score = round(sum(i.score for i in items), 1)
        dim.details = prop_values

    # ------------------------------------------------------------------
    # Dimension 3: Complexity & Feasibility (25 pts)
    # ------------------------------------------------------------------

    def _score_complexity_feasibility(
        self, mol: Chem.Mol, breakdown: MLReadinessBreakdown
    ) -> None:
        """QED, SA Score, Fsp3, and stereocenter manageability."""
        from app.services.validation.checks.deep_stereo_tautomer import (
            StereoisomerEnumerationCheck,
        )

        dim = breakdown.complexity_feasibility
        items: List[MLDimensionItem] = []

        # QED (8 pts)
        try:
            from rdkit.Chem import QED
            qed_value = QED.qed(mol)
        except Exception:
            qed_value = 0.0
        qed_pts = round(qed_value * 8, 2)
        items.append(MLDimensionItem(
            "QED contribution", qed_pts, 8.0,
            subtitle=f"QED = {qed_value:.3f}",
            tooltip="Quantitative Estimate of Drug-likeness (Bickerton 2012). "
            "Integrates MW, LogP, HBD, HBA, PSA, RotBonds, AromRings, and alerts "
            "into a single 0–1 score. Linearly scaled to points.",
        ))

        # SA Score (8 pts) — inverse mapping
        try:
            import os
            import sys

            from rdkit.Chem import RDConfig
            sa_model_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
            if sa_model_path not in sys.path:
                sys.path.insert(0, sa_model_path)
            from sascorer import calculateScore
            sa_score = calculateScore(mol)
        except Exception:
            sa_score = 5.0  # Neutral fallback
        sa_pts = self._sa_score_to_points(sa_score)
        items.append(MLDimensionItem(
            "Synthesizability", sa_pts, 8.0,
            subtitle=f"SA = {sa_score:.1f} (lower = easier)",
            tooltip="Synthetic Accessibility Score (Ertl 2009). Range 1 (easy) to 10 "
            "(hard). Inversely mapped: SA ≤ 3 → full points, SA ≥ 7 → zero.",
        ))

        # Fsp3 (4 pts)
        try:
            fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
        except Exception:
            fsp3 = 0.0
        fsp3_d = desirability(fsp3, 0.2, 0.6) or 0.0
        fsp3_pts = round(fsp3_d * 4, 2)
        items.append(MLDimensionItem(
            "Fsp3 fitness", fsp3_pts, 4.0,
            subtitle=f"Fsp3 = {fsp3:.3f} (ideal: 0.2\u20130.6)",
            tooltip="Fraction of sp3-hybridized carbons (Lovering 2009). Higher Fsp3 "
            "correlates with clinical success. Desirability-scored against 0.2–0.6 range.",
        ))

        # Stereocenters (5 pts)
        try:
            stereo_result = StereoisomerEnumerationCheck().run(mol)
            total_centers = stereo_result.details.get("total_centers", 0)
            undefined_count = stereo_result.details.get("undefined_count", 0)
        except Exception:
            total_centers = 0
            undefined_count = 0

        stereo_pts = self._stereocenter_points(total_centers, undefined_count)
        if undefined_count > 0:
            stereo_subtitle = f"{total_centers} centers ({undefined_count} undefined)"
        elif total_centers == 0:
            stereo_subtitle = "No stereocenters"
        else:
            stereo_subtitle = f"{total_centers} centers (\u22644 ideal)"
        items.append(MLDimensionItem(
            "Stereo. manageability", stereo_pts, 5.0,
            subtitle=stereo_subtitle,
            tooltip="Scores how manageable the stereocenter count is for ML. "
            "0–4 centers: full points. 5–8: partial. >8: zero. "
            ">50% undefined centers halves the score.",
        ))

        dim.items = items
        dim.score = round(sum(i.score for i in items), 1)
        dim.details = {
            "qed": round(qed_value, 3),
            "sa_score": round(sa_score, 2),
            "fsp3": round(fsp3, 3),
            "total_stereocenters": total_centers,
            "undefined_stereocenters": undefined_count,
        }

    @staticmethod
    def _sa_score_to_points(sa: float) -> float:
        """Convert SA Score (1-10) to points (0-8). Lower SA = more points."""
        if sa <= 3:
            return 8.0
        elif sa <= 5:
            # Linear interpolation: 3→8, 5→4
            return round(8.0 - (sa - 3) * 2.0, 2)
        elif sa <= 7:
            # Linear interpolation: 5→4, 7→0
            return round(4.0 - (sa - 5) * 2.0, 2)
        else:
            return 0.0

    @staticmethod
    def _stereocenter_points(total: int, undefined: int) -> float:
        """Score stereocenter manageability (0-5 pts)."""
        if total <= 4:
            base = 5.0
        elif total <= 8:
            # Linear: 5→5, 8→2
            base = round(5.0 - (total - 4) * 0.75, 2)
        else:
            base = 0.0

        # Penalty: >50% undefined → halve the score
        if total > 0 and undefined / total > 0.5:
            base *= 0.5

        return round(base, 2)

    # ------------------------------------------------------------------
    # Dimension 4: Representation Quality (20 pts)
    # ------------------------------------------------------------------

    def _score_representation_quality(
        self, mol: Chem.Mol, breakdown: MLReadinessBreakdown
    ) -> None:
        """Descriptor completeness, fingerprint generation + informativeness, conformer."""
        dim = breakdown.representation_quality
        items: List[MLDimensionItem] = []

        # 1. Descriptor completeness (5 pts)
        desc_pts, desc_total, desc_success = self._score_descriptors(mol)
        items.append(MLDimensionItem(
            "Descriptor completeness", desc_pts, 5.0,
            subtitle=f"{desc_success}/{desc_total} computed successfully",
            tooltip="Fraction of RDKit molecular descriptors that compute without "
            "error. Includes 217 standard + 192 AUTOCORR2D + 42 MQN = 451 total.",
        ))

        # 2. Fingerprint generation (5 pts)
        fp_pts, fp_successful, fp_failed, fps = self._score_fingerprints(mol)
        items.append(MLDimensionItem(
            "Fingerprint generation", fp_pts, 5.0,
            subtitle=f"{len(fp_successful)}/7 types generated",
            tooltip="Tests generation of 7 fingerprint types: Morgan, Morgan features, "
            "MACCS, atom pair, topological torsion, RDKit, and Avalon.",
        ))

        # 3. Fingerprint informativeness (5 pts)
        info_pts = self._score_fingerprint_informativeness(fps)
        items.append(MLDimensionItem(
            "Fingerprint informativeness", info_pts, 5.0,
            tooltip="Average bit density across generated fingerprints. "
            "Ideal range: 1–30%. Too sparse (<1%) or saturated (>30%) "
            "fingerprints carry less discriminative information for ML models.",
        ))

        # 4. Conformer generation (5 pts)
        conf_pts = self._score_conformer_generation(mol)
        items.append(MLDimensionItem(
            "Conformer generation", conf_pts, 5.0,
            subtitle="ETKDGv3" if conf_pts == 5.0
            else "ETKDGv3 (fallback)" if conf_pts > 0 else "Failed",
            tooltip="Tests 3D conformer embedding via ETKDGv3 (Riniker 2015). "
            "Falls back to random coordinates if standard embedding fails. "
            "Conformer feasibility indicates suitability for 3D-QSAR and docking.",
        ))

        dim.items = items
        dim.score = round(sum(i.score for i in items), 1)
        dim.details = {
            "descriptors_successful": desc_success,
            "descriptors_total": desc_total,
            "fingerprints_successful": fp_successful,
            "fingerprints_failed": fp_failed,
        }

    def _score_descriptors(self, mol: Chem.Mol) -> tuple:
        """Score descriptor completeness (5 pts max). Returns (pts, total, successful).

        Tests 451 descriptors: 217 standard (CalcMolDescriptors) +
        192 AUTOCORR2D + 42 MQN.
        """
        import math

        total = 0
        successful = 0

        # Standard descriptors (217)
        try:
            descriptors = Descriptors.CalcMolDescriptors(mol, missingVal=None, silent=True)
            total += len(descriptors)
            successful += sum(
                1 for v in descriptors.values()
                if v is not None and not (isinstance(v, float) and math.isnan(v))
            )
        except Exception:
            pass

        # AUTOCORR2D (192 values)
        try:
            autocorr = rdMolDescriptors.CalcAUTOCORR2D(mol)
            total += len(autocorr)
            successful += sum(
                1 for v in autocorr
                if v is not None and not (isinstance(v, float) and math.isnan(v))
            )
        except Exception:
            total += 192  # Count as failures

        # MQN (42 values)
        try:
            mqn = rdMolDescriptors.MQNs_(mol)
            total += len(mqn)
            successful += sum(
                1 for v in mqn
                if v is not None and not (isinstance(v, float) and math.isnan(v))
            )
        except Exception:
            total += 42  # Count as failures

        if total > 0:
            pts = round(5.0 * (successful / total), 2)
        else:
            pts = 0.0
        return pts, total, successful

    def _score_fingerprints(self, mol: Chem.Mol) -> tuple:
        """Score FP generation (5 pts max). Returns (pts, successful_list, failed_list, fp_dict)."""
        total_weight = sum(w for _, w in self.FINGERPRINT_TYPES)
        earned = 0.0
        successful = []
        failed = []
        fps = {}

        for fp_name, weight in self.FINGERPRINT_TYPES:
            try:
                fp = self._generate_fingerprint(mol, fp_name)
                if fp is not None:
                    earned += weight
                    successful.append(fp_name)
                    fps[fp_name] = fp
                else:
                    failed.append(fp_name)
            except Exception:
                failed.append(fp_name)

        pts = round(5.0 * (earned / total_weight), 2) if total_weight > 0 else 0.0
        return pts, successful, failed, fps

    def _generate_fingerprint(self, mol: Chem.Mol, fp_name: str):
        """Generate a specific fingerprint type."""
        if fp_name == "morgan":
            return self._morgan_gen.GetFingerprint(mol)
        elif fp_name == "morgan_features":
            return self._morgan_feat_gen.GetFingerprint(mol)
        elif fp_name == "maccs":
            return MACCSkeys.GenMACCSKeys(mol)
        elif fp_name == "atompair":
            return self._atompair_gen.GetFingerprint(mol)
        elif fp_name == "topological_torsion":
            return self._torsion_gen.GetFingerprint(mol)
        elif fp_name == "rdkit_fp":
            return self._rdkit_gen.GetFingerprint(mol)
        elif fp_name == "avalon":
            return pyAvalonTools.GetAvalonFP(mol)
        return None

    def _score_fingerprint_informativeness(self, fps: dict) -> float:
        """Score average bit density across fingerprints (5 pts). Ideal: 1-30%."""
        if not fps:
            return 0.0

        densities = []
        for fp in fps.values():
            try:
                num_bits = fp.GetNumBits()
                on_bits = fp.GetNumOnBits()
                if num_bits > 0:
                    densities.append(on_bits / num_bits)
            except AttributeError:
                # Some FP types (e.g. Avalon) return different object types
                continue

        if not densities:
            return 0.0

        avg_density = sum(densities) / len(densities)

        # Ideal range 1-30%
        if 0.01 <= avg_density <= 0.30:
            return 5.0
        elif avg_density < 0.01:
            # Too sparse
            return round(5.0 * (avg_density / 0.01), 2)
        elif avg_density <= 0.45:
            # Slightly saturated — linear falloff 30%→45%
            return round(5.0 * (1.0 - (avg_density - 0.30) / 0.15), 2)
        else:
            return 0.0

    def _score_conformer_generation(self, mol: Chem.Mol) -> float:
        """Score 3D conformer feasibility (5 pts). ETKDGv3 + randomCoords fallback."""
        try:
            work_mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.maxIterations = 500

            result = AllChem.EmbedMolecule(work_mol, params)
            if result == 0:
                return 5.0

            # Fallback: random coords
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(work_mol, params)
            if result == 0:
                return 3.0  # Partial credit for needing fallback

            return 0.0
        except Exception:
            return 0.0

    # ------------------------------------------------------------------
    # Supplementary (informational, not scored)
    # ------------------------------------------------------------------

    def _get_supplementary(self, mol: Chem.Mol) -> Dict:
        """Get supplementary info (Lipinski violations, Veber pass/fail)."""
        supplementary: Dict = {}
        try:
            from app.services.scoring.druglikeness import calculate_druglikeness
            dl = calculate_druglikeness(mol, include_extended=False)
            supplementary["lipinski_violations"] = dl.lipinski.violations
            supplementary["lipinski_passed"] = dl.lipinski.passed
            supplementary["veber_passed"] = dl.veber.passed
        except Exception:
            pass
        return supplementary

    # ------------------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------------------

    def _get_interpretation(
        self,
        score: int,
        breakdown: MLReadinessBreakdown,
        caveats: List[str],
    ) -> str:
        """Generate human-readable interpretation with practical guidance."""
        label, _ = _label_and_color(score)

        # Practical meaning per tier
        tier_guidance = {
            "Excellent": "suitable for most ML workflows without modification",
            "Good": "suitable for ML with minor preprocessing",
            "Moderate": "usable but may need filtering or augmentation in ML pipelines",
            "Limited": "likely to cause issues in ML models — consider alternatives",
            "Poor": "not recommended for ML use without significant preprocessing",
        }
        guidance = tier_guidance.get(label, "")
        parts = [f"{label} ML-readiness ({score}/100) — {guidance}."]

        dims = [
            breakdown.structural_quality,
            breakdown.property_profile,
            breakdown.complexity_feasibility,
            breakdown.representation_quality,
        ]

        # Classify dimensions by strength
        strong = []
        weak = []
        for dim in dims:
            pct = dim.score / dim.max_score if dim.max_score > 0 else 0
            if pct >= 0.8:
                strong.append(dim.name)
            elif pct < 0.5:
                weak.append(f"{dim.name} ({dim.score:.0f}/{dim.max_score:.0f})")

        if strong:
            parts.append(f"Strongest: {', '.join(strong)}.")
        if weak:
            parts.append(f"Needs improvement: {', '.join(weak)}.")

        if caveats:
            parts.append(f"Caveats: {'; '.join(caveats)}.")

        return " ".join(parts)


# Module-level convenience
_scorer = MLReadinessScorer()


def calculate_ml_readiness(mol: Chem.Mol) -> MLReadinessResult:
    """Calculate ML-readiness score for a molecule.

    Args:
        mol: RDKit molecule object

    Returns:
        MLReadinessResult with score (0-100), 4-dimension breakdown,
        caveats, supplementary info, and interpretation.
    """
    return _scorer.score(mol)
