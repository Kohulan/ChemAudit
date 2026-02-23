"""
Salt/Counterion Inventory Scorer

Identifies and classifies salt fragments, counterions, and solvents
in a molecule without requiring prior standardization.

Also provides ligand efficiency calculations.
"""

from dataclasses import dataclass, field
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, QED

from app.services.standardization.fragment_dict import classify_fragment


@dataclass
class SaltFragment:
    """A single fragment from a salt-form molecule."""

    smiles: str
    name: str
    category: str  # "counterion", "salt", "solvent", "drug", "unknown"
    mw: float
    heavy_atom_count: int


@dataclass
class SaltInventoryResult:
    """Salt/counterion inventory result."""

    has_salts: bool
    parent_smiles: str
    fragments: List[SaltFragment] = field(default_factory=list)
    total_fragments: int = 0
    interpretation: str = ""


@dataclass
class LigandEfficiencyResult:
    """Ligand efficiency calculation result."""

    le: Optional[float]
    heavy_atom_count: int
    activity_value: Optional[float] = None
    activity_type: Optional[str] = None
    proxy_used: bool = False
    interpretation: str = ""


class SaltInventoryScorer:
    """
    Identifies and classifies fragments in salt-form molecules.

    Works standalone without prior standardization by splitting on
    dot-separated SMILES and classifying each fragment.
    """

    def calculate(
        self,
        mol: Chem.Mol,
        provenance_fragments: Optional[List[str]] = None,
    ) -> SaltInventoryResult:
        """
        Calculate salt inventory for a molecule.

        Args:
            mol: RDKit molecule object
            provenance_fragments: Optional list of fragment SMILES from
                standardization provenance. If provided, these are used
                instead of splitting the input molecule.

        Returns:
            SaltInventoryResult with fragment classification
        """
        smiles = Chem.MolToSmiles(mol)

        if provenance_fragments is not None:
            fragment_smiles_list = provenance_fragments
        elif "." in smiles:
            fragment_smiles_list = smiles.split(".")
        else:
            # Single component — no salts
            return SaltInventoryResult(
                has_salts=False,
                parent_smiles=smiles,
                total_fragments=1,
                interpretation="Single-component molecule with no salt forms detected.",
            )

        if not fragment_smiles_list:
            return SaltInventoryResult(
                has_salts=False,
                parent_smiles=smiles,
                total_fragments=0,
                interpretation="No fragments found.",
            )

        # Classify each fragment
        classified = []
        for frag_smi in fragment_smiles_list:
            frag_mol = Chem.MolFromSmiles(frag_smi)
            if frag_mol is None:
                classified.append(
                    SaltFragment(
                        smiles=frag_smi,
                        name="unknown",
                        category="unknown",
                        mw=0.0,
                        heavy_atom_count=0,
                    )
                )
                continue

            info = classify_fragment(frag_smi)
            ha_count = frag_mol.GetNumHeavyAtoms()
            mw = Descriptors.MolWt(frag_mol)

            classified.append(
                SaltFragment(
                    smiles=info["smiles"],
                    name=info["name"] or "unknown",
                    category=info["role"],
                    mw=round(mw, 2),
                    heavy_atom_count=ha_count,
                )
            )

        # Identify the parent (largest fragment by heavy atom count)
        if classified:
            parent_idx = max(range(len(classified)), key=lambda i: classified[i].heavy_atom_count)
            classified[parent_idx] = SaltFragment(
                smiles=classified[parent_idx].smiles,
                name=classified[parent_idx].name
                if classified[parent_idx].category != "unknown"
                else "parent compound",
                category="drug",
                mw=classified[parent_idx].mw,
                heavy_atom_count=classified[parent_idx].heavy_atom_count,
            )
            parent_smiles = classified[parent_idx].smiles
        else:
            parent_smiles = smiles

        # Non-drug fragments
        non_drug = [f for f in classified if f.category != "drug"]
        has_salts = len(non_drug) > 0

        # Build interpretation
        if has_salts:
            salt_names = [f.name for f in non_drug if f.name != "unknown"]
            if salt_names:
                interpretation = (
                    f"Salt-form molecule with {len(non_drug)} non-drug "
                    f"fragment(s): {', '.join(salt_names)}."
                )
            else:
                interpretation = (
                    f"Salt-form molecule with {len(non_drug)} unidentified "
                    f"fragment(s)."
                )
        else:
            interpretation = "Single-component molecule with no salt forms detected."

        return SaltInventoryResult(
            has_salts=has_salts,
            parent_smiles=parent_smiles,
            fragments=classified,
            total_fragments=len(classified),
            interpretation=interpretation,
        )


class LigandEfficiencyScorer:
    """
    Calculates ligand efficiency metrics.

    Standard LE = 1.37 * pIC50 / heavy_atom_count
    When no activity is provided, uses a BEI proxy based on QED.
    """

    def calculate(
        self,
        mol: Chem.Mol,
        activity_value: Optional[float] = None,
        activity_type: str = "pIC50",
    ) -> LigandEfficiencyResult:
        """
        Calculate ligand efficiency.

        Args:
            mol: RDKit molecule object
            activity_value: Activity value (e.g., pIC50, pKi)
            activity_type: Type of activity measurement

        Returns:
            LigandEfficiencyResult with efficiency metrics
        """
        heavy_atom_count = mol.GetNumHeavyAtoms()

        if heavy_atom_count == 0:
            return LigandEfficiencyResult(
                le=None,
                heavy_atom_count=0,
                interpretation="Cannot compute LE for molecule with 0 heavy atoms.",
            )

        if activity_value is not None:
            # Standard LE formula: deltaG / N_heavy ~ 1.37 * pActivity / N_heavy
            le = 1.37 * activity_value / heavy_atom_count
            interpretation = (
                f"LE = {le:.3f} kcal/mol/atom from {activity_type}={activity_value}. "
            )
            if le >= 0.3:
                interpretation += "Excellent ligand efficiency."
            elif le >= 0.2:
                interpretation += "Good ligand efficiency."
            else:
                interpretation += "Poor ligand efficiency."

            return LigandEfficiencyResult(
                le=round(le, 4),
                heavy_atom_count=heavy_atom_count,
                activity_value=activity_value,
                activity_type=activity_type,
                proxy_used=False,
                interpretation=interpretation,
            )
        else:
            # BEI proxy: QED * 1000 / MW
            mw = Descriptors.MolWt(mol)
            if mw <= 0:
                return LigandEfficiencyResult(
                    le=None,
                    heavy_atom_count=heavy_atom_count,
                    proxy_used=True,
                    interpretation="Cannot compute proxy LE: MW is zero.",
                )

            try:
                qed_score = QED.qed(mol)
            except Exception:
                qed_score = 0.0

            bei_proxy = qed_score * 1000 / mw
            interpretation = (
                f"BEI proxy = {bei_proxy:.3f} (QED-based, no activity data). "
            )
            if bei_proxy >= 1.5:
                interpretation += "High efficiency proxy."
            elif bei_proxy >= 0.8:
                interpretation += "Moderate efficiency proxy."
            else:
                interpretation += "Low efficiency proxy."

            return LigandEfficiencyResult(
                le=round(bei_proxy, 4),
                heavy_atom_count=heavy_atom_count,
                proxy_used=True,
                interpretation=interpretation,
            )


# Module-level instances
_salt_scorer = SaltInventoryScorer()
_le_scorer = LigandEfficiencyScorer()


def calculate_salt_inventory(
    mol: Chem.Mol,
    provenance_fragments: Optional[List[str]] = None,
) -> SaltInventoryResult:
    """
    Calculate salt/counterion inventory for a molecule.

    Args:
        mol: RDKit molecule object
        provenance_fragments: Optional fragment SMILES from standardization

    Returns:
        SaltInventoryResult with fragment classification
    """
    return _salt_scorer.calculate(mol, provenance_fragments)


def calculate_ligand_efficiency(
    mol: Chem.Mol,
    activity_value: Optional[float] = None,
    activity_type: str = "pIC50",
) -> LigandEfficiencyResult:
    """
    Calculate ligand efficiency.

    Args:
        mol: RDKit molecule object
        activity_value: Optional activity value (pIC50, pKi, etc.)
        activity_type: Type of activity measurement

    Returns:
        LigandEfficiencyResult with efficiency metrics
    """
    return _le_scorer.calculate(mol, activity_value, activity_type)
