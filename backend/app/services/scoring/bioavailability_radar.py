"""
Bioavailability Radar, BOILED-Egg, and Property Comparison Scorer

Provides:
- Bioavailability radar: 6-axis normalized assessment (LIPO, SIZE, POLAR, INSOLU, INSATU, FLEX)
- BOILED-Egg classification: GI absorption / BBB permeation via elliptical model
- Property radar comparison: multi-molecule radar profiles with reference
"""

from dataclasses import dataclass, field
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors


@dataclass
class RadarAxis:
    """A single axis of the bioavailability radar."""

    name: str  # "LIPO", "SIZE", "POLAR", "INSOLU", "INSATU", "FLEX"
    actual_value: float
    normalized: float  # 0-1
    optimal_min: float
    optimal_max: float
    in_range: bool
    property_name: str  # "WLOGP", "MW", "TPSA", "LogS_ESOL", "Fsp3", "RotBonds"
    unit: str


@dataclass
class BioavailabilityRadarResult:
    """Bioavailability radar result with 6 axes."""

    axes: List[RadarAxis] = field(default_factory=list)
    overall_in_range_count: int = 0
    interpretation: str = ""


@dataclass
class EllipseParams:
    """Ellipse parameters for BOILED-Egg model."""

    cx: float
    cy: float
    a: float  # semi-axis x (TPSA axis)
    b: float  # semi-axis y (WLOGP axis)


@dataclass
class BoiledEggResult:
    """BOILED-Egg classification result."""

    wlogp: float
    tpsa: float
    gi_absorbed: bool
    bbb_permeant: bool
    region: str  # "yolk", "white", "grey"
    gi_ellipse: Optional[EllipseParams] = None
    bbb_ellipse: Optional[EllipseParams] = None
    interpretation: str = ""


@dataclass
class RadarProfile:
    """A single molecule's radar profile."""

    smiles: str
    axes: List[RadarAxis] = field(default_factory=list)
    is_reference: bool = False


@dataclass
class RadarComparisonResult:
    """Multi-molecule radar comparison."""

    profiles: List[RadarProfile] = field(default_factory=list)
    reference: Optional[RadarProfile] = None


# ESOL model coefficients (Delaney 2004) — duplicated from admet.py to avoid import cycle
_ESOL_INTERCEPT = 0.16
_ESOL_LOGP_COEF = -0.63
_ESOL_MW_COEF = -0.0062
_ESOL_ROTB_COEF = 0.066
_ESOL_AP_COEF = -0.74


def _calculate_esol_logs(mol: Chem.Mol) -> float:
    """Calculate ESOL LogS for the bioavailability radar INSOLU axis."""
    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    num_atoms = mol.GetNumHeavyAtoms()
    num_aromatic = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    ap = num_aromatic / num_atoms if num_atoms > 0 else 0
    return (
        _ESOL_INTERCEPT
        + _ESOL_LOGP_COEF * logp
        + _ESOL_MW_COEF * mw
        + _ESOL_ROTB_COEF * rotatable_bonds
        + _ESOL_AP_COEF * ap
    )


def _normalize_axis(value: float, opt_min: float, opt_max: float) -> float:
    """
    Normalize an axis value to 0-1 scale.

    If value is within [opt_min, opt_max], returns 1.0.
    If outside, decreases linearly toward 0 based on distance from range.
    """
    axis_range = opt_max - opt_min
    if axis_range <= 0:
        return 1.0 if opt_min <= value <= opt_max else 0.0

    if opt_min <= value <= opt_max:
        return 1.0
    elif value < opt_min:
        return max(0.0, 1.0 - (opt_min - value) / axis_range)
    else:
        return max(0.0, 1.0 - (value - opt_max) / axis_range)


class BioavailabilityRadarScorer:
    """
    Provides bioavailability radar, BOILED-Egg, and property comparison.
    """

    # Bioavailability radar axis definitions: (name, property_name, opt_min, opt_max, unit)
    AXES = [
        ("LIPO", "WLOGP", -0.7, 5.0, ""),
        ("SIZE", "MW", 150.0, 500.0, "Da"),
        ("POLAR", "TPSA", 20.0, 130.0, "A^2"),
        ("INSOLU", "LogS_ESOL", -6.0, 0.0, "log mol/L"),
        ("INSATU", "Fsp3", 0.25, 1.0, ""),
        ("FLEX", "RotBonds", 0.0, 9.0, ""),
    ]

    # BOILED-Egg ellipse parameters (Daina & Zoete 2016)
    # GI ellipse (white): center and semi-axes on (TPSA, WLOGP) plane
    GI_ELLIPSE = EllipseParams(cx=71.051, cy=2.292, a=142.081, b=4.580)
    # BBB ellipse (yolk): tighter region
    BBB_ELLIPSE = EllipseParams(cx=45.542, cy=2.137, a=79.579, b=3.177)

    def calculate_radar(self, mol: Chem.Mol) -> BioavailabilityRadarResult:
        """
        Calculate bioavailability radar with 6 normalized axes.

        Args:
            mol: RDKit molecule object

        Returns:
            BioavailabilityRadarResult with 6 axes
        """
        # Compute descriptors
        logp = Crippen.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        tpsa = Descriptors.TPSA(mol)
        log_s = _calculate_esol_logs(mol)
        fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
        rotbonds = Lipinski.NumRotatableBonds(mol)

        values = {
            "WLOGP": logp,
            "MW": mw,
            "TPSA": tpsa,
            "LogS_ESOL": log_s,
            "Fsp3": fsp3,
            "RotBonds": float(rotbonds),
        }

        axes = []
        for name, prop_name, opt_min, opt_max, unit in self.AXES:
            actual = values[prop_name]
            normalized = _normalize_axis(actual, opt_min, opt_max)
            in_range = opt_min <= actual <= opt_max

            axes.append(
                RadarAxis(
                    name=name,
                    actual_value=round(actual, 4),
                    normalized=round(normalized, 4),
                    optimal_min=opt_min,
                    optimal_max=opt_max,
                    in_range=in_range,
                    property_name=prop_name,
                    unit=unit,
                )
            )

        in_range_count = sum(1 for a in axes if a.in_range)

        if in_range_count == 6:
            interpretation = "All properties within optimal ranges — excellent bioavailability profile."
        elif in_range_count >= 4:
            interpretation = (
                f"{in_range_count}/6 properties in range — good bioavailability profile."
            )
        elif in_range_count >= 2:
            interpretation = (
                f"{in_range_count}/6 properties in range — moderate bioavailability profile."
            )
        else:
            interpretation = (
                f"{in_range_count}/6 properties in range — poor bioavailability profile."
            )

        return BioavailabilityRadarResult(
            axes=axes,
            overall_in_range_count=in_range_count,
            interpretation=interpretation,
        )

    def calculate_boiled_egg(self, mol: Chem.Mol) -> BoiledEggResult:
        """
        Calculate BOILED-Egg classification for GI absorption / BBB permeation.

        Uses elliptical model from Daina & Zoete (2016).
        TPSA is x-axis, WLOGP is y-axis.

        Args:
            mol: RDKit molecule object

        Returns:
            BoiledEggResult with classification
        """
        wlogp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)

        # Point-in-ellipse test: ((x-cx)/a)^2 + ((y-cy)/b)^2 <= 1
        gi_test = (
            ((tpsa - self.GI_ELLIPSE.cx) / self.GI_ELLIPSE.a) ** 2
            + ((wlogp - self.GI_ELLIPSE.cy) / self.GI_ELLIPSE.b) ** 2
        ) <= 1

        bbb_test = (
            ((tpsa - self.BBB_ELLIPSE.cx) / self.BBB_ELLIPSE.a) ** 2
            + ((wlogp - self.BBB_ELLIPSE.cy) / self.BBB_ELLIPSE.b) ** 2
        ) <= 1

        if bbb_test:
            region = "yolk"
            interpretation = (
                "BBB permeant (yolk region): predicted to cross the blood-brain barrier. "
                "Also predicted to be GI-absorbed."
            )
        elif gi_test:
            region = "white"
            interpretation = (
                "GI absorbed (white region): predicted to be passively absorbed "
                "from the gastrointestinal tract. Not predicted to cross the BBB."
            )
        else:
            region = "grey"
            interpretation = (
                "Grey region: predicted to have low passive GI absorption and "
                "not to cross the BBB."
            )

        return BoiledEggResult(
            wlogp=round(wlogp, 4),
            tpsa=round(tpsa, 4),
            gi_absorbed=gi_test,
            bbb_permeant=bbb_test,
            region=region,
            gi_ellipse=self.GI_ELLIPSE,
            bbb_ellipse=self.BBB_ELLIPSE,
            interpretation=interpretation,
        )

    def calculate_comparison(
        self, smiles_list: List[str]
    ) -> RadarComparisonResult:
        """
        Compare multiple molecules using radar profiles.

        Args:
            smiles_list: List of SMILES strings to compare

        Returns:
            RadarComparisonResult with profiles and reference
        """
        profiles = []

        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            radar = self.calculate_radar(mol)
            profiles.append(
                RadarProfile(
                    smiles=smi,
                    axes=radar.axes,
                    is_reference=False,
                )
            )

        # Build reference profile from average drug-like values
        ref_values = {
            "WLOGP": 2.5,
            "MW": 350.0,
            "TPSA": 75.0,
            "LogS_ESOL": -3.0,
            "Fsp3": 0.35,
            "RotBonds": 5.0,
        }

        ref_axes = []
        for name, prop_name, opt_min, opt_max, unit in self.AXES:
            actual = ref_values[prop_name]
            normalized = _normalize_axis(actual, opt_min, opt_max)
            in_range = opt_min <= actual <= opt_max

            ref_axes.append(
                RadarAxis(
                    name=name,
                    actual_value=round(actual, 4),
                    normalized=round(normalized, 4),
                    optimal_min=opt_min,
                    optimal_max=opt_max,
                    in_range=in_range,
                    property_name=prop_name,
                    unit=unit,
                )
            )

        reference = RadarProfile(
            smiles="reference",
            axes=ref_axes,
            is_reference=True,
        )

        return RadarComparisonResult(
            profiles=profiles,
            reference=reference,
        )


# Module-level instance and convenience functions
_scorer = BioavailabilityRadarScorer()


def calculate_bioavailability_radar(mol: Chem.Mol) -> BioavailabilityRadarResult:
    """Calculate bioavailability radar with 6 normalized axes."""
    return _scorer.calculate_radar(mol)


def calculate_boiled_egg(mol: Chem.Mol) -> BoiledEggResult:
    """Calculate BOILED-Egg GI/BBB classification."""
    return _scorer.calculate_boiled_egg(mol)


def calculate_radar_comparison(smiles_list: List[str]) -> RadarComparisonResult:
    """Compare multiple molecules using radar profiles."""
    return _scorer.calculate_comparison(smiles_list)
