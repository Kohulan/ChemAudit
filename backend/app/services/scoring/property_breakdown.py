"""
Property Breakdown Scorer

Provides per-atom and per-functional-group property decomposition:
- TPSA breakdown: per-atom topological polar surface area contributions
- LogP breakdown: per-atom Crippen LogP contributions (with H folding)
- Bertz complexity detail: graph-level metrics
- Fsp3 detail: per-carbon hybridization enumeration
"""

from dataclasses import dataclass, field
from typing import List

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, GraphDescriptors, rdMolDescriptors


@dataclass
class AtomContribution:
    """Per-atom property contribution."""

    atom_index: int
    symbol: str
    contribution: float


@dataclass
class FunctionalGroupContribution:
    """Per-functional-group property contribution."""

    group_name: str
    contribution: float
    atom_indices: List[int] = field(default_factory=list)


@dataclass
class TPSABreakdownResult:
    """TPSA per-atom breakdown result."""

    total: float
    atom_contributions: List[AtomContribution] = field(default_factory=list)
    functional_group_summary: List[FunctionalGroupContribution] = field(
        default_factory=list
    )


@dataclass
class LogPBreakdownResult:
    """LogP per-atom breakdown result."""

    total: float
    atom_contributions: List[AtomContribution] = field(default_factory=list)
    functional_group_summary: List[FunctionalGroupContribution] = field(
        default_factory=list
    )


@dataclass
class BertzDetailResult:
    """Bertz complexity detail result."""

    bertz_ct: float
    num_bonds: int
    num_atoms: int
    num_rings: int
    num_aromatic_rings: int
    ring_complexity: float
    interpretation: str = ""


@dataclass
class CarbonHybridization:
    """Per-carbon hybridization data."""

    atom_index: int
    symbol: str = "C"
    hybridization: str = "sp3"  # "sp", "sp2", "sp3", "other"


@dataclass
class Fsp3DetailResult:
    """Fsp3 per-carbon detail result."""

    fsp3: float
    total_carbons: int
    sp3_count: int
    sp2_count: int
    sp_count: int
    per_carbon: List[CarbonHybridization] = field(default_factory=list)
    interpretation: str = ""


def _classify_tpsa_group(atom, mol) -> str:
    """Classify an atom contributing to TPSA into a functional group."""
    symbol = atom.GetSymbol()
    neighbors = atom.GetNeighbors()
    num_hs = atom.GetTotalNumHs()

    if symbol == "N":
        # Check for amide: N bonded to C which is bonded to O via double bond
        for nbr in neighbors:
            if nbr.GetSymbol() == "C":
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetSymbol() == "O":
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            if num_hs >= 1:
                                return "amide NH"
                            return "amide N"
        if atom.GetIsAromatic():
            return "aromatic N"
        if num_hs == 2:
            return "primary amine"
        if num_hs == 1:
            return "secondary amine"
        if num_hs == 0:
            return "tertiary amine"
        return "nitrogen"
    elif symbol == "O":
        if atom.GetIsAromatic():
            return "aromatic O"
        # Check if in a ring
        if atom.IsInRing():
            return "ring oxygen"
        # Check for carboxyl: O bonded to C which also has =O
        for nbr in neighbors:
            if nbr.GetSymbol() == "C":
                has_double_o = False
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() != atom.GetIdx() and nbr2.GetSymbol() == "O":
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            has_double_o = True
                if has_double_o and num_hs >= 1:
                    return "carboxyl OH"
                if has_double_o and num_hs == 0:
                    return "ester oxygen"
        # Check carbonyl
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond and bond.GetBondTypeAsDouble() == 2.0 and nbr.GetSymbol() == "C":
                return "carbonyl"
        if num_hs >= 1:
            return "hydroxyl"
        return "ether oxygen"
    elif symbol == "S":
        if atom.GetIsAromatic():
            return "aromatic S"
        # Check for sulfonamide/sulfone
        double_o_count = 0
        for nbr in neighbors:
            if nbr.GetSymbol() == "O":
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond and bond.GetBondTypeAsDouble() == 2.0:
                    double_o_count += 1
        if double_o_count >= 2:
            return "sulfone"
        if double_o_count == 1:
            return "sulfoxide"
        if num_hs >= 1:
            return "thiol"
        return "thioether"
    elif symbol == "P":
        return "phosphorus"
    else:
        return f"{symbol} (polar)"


def _classify_logp_group(atom, mol) -> str:
    """Classify an atom for LogP functional group summary."""
    symbol = atom.GetSymbol()
    neighbors = atom.GetNeighbors()

    if symbol == "C":
        if atom.GetIsAromatic():
            return "aromatic C"
        # Check if bonded to any heteroatom
        has_hetero = any(
            n.GetSymbol() not in ("C", "H") for n in neighbors
        )
        if has_hetero:
            return "C (hetero-bonded)"
        return "aliphatic C"
    elif symbol == "N":
        if atom.GetIsAromatic():
            return "aromatic N"
        return "nitrogen"
    elif symbol == "O":
        if atom.GetIsAromatic():
            return "aromatic O"
        num_hs = atom.GetTotalNumHs()
        if num_hs >= 1:
            return "hydroxyl/OH"
        # Check for carbonyl
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond and bond.GetBondTypeAsDouble() == 2.0:
                return "carbonyl O"
        return "ether O"
    elif symbol == "S":
        return "sulfur"
    elif symbol == "F":
        return "fluorine"
    elif symbol == "Cl":
        return "chlorine"
    elif symbol == "Br":
        return "bromine"
    elif symbol == "I":
        return "iodine"
    else:
        return f"{symbol}"


class PropertyBreakdownScorer:
    """
    Provides per-atom and per-group property decomposition.
    """

    def calculate_tpsa_breakdown(self, mol: Chem.Mol) -> TPSABreakdownResult:
        """
        Calculate TPSA per-atom breakdown.

        Uses rdMolDescriptors._CalcTPSAContribs for per-atom TPSA values.
        Does NOT call AddHs — TPSA contribs work on heavy atoms.

        Args:
            mol: RDKit molecule object

        Returns:
            TPSABreakdownResult with per-atom and per-group contributions
        """
        contribs = rdMolDescriptors._CalcTPSAContribs(mol)

        atom_contributions = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            atom_contributions.append(
                AtomContribution(
                    atom_index=i,
                    symbol=atom.GetSymbol(),
                    contribution=round(contribs[i], 4),
                )
            )

        # Build functional group summary for atoms with nonzero contribution
        group_map: dict = {}
        for i, ac in enumerate(atom_contributions):
            if ac.contribution > 0:
                atom = mol.GetAtomWithIdx(i)
                group_name = _classify_tpsa_group(atom, mol)
                if group_name not in group_map:
                    group_map[group_name] = {"contribution": 0.0, "atom_indices": []}
                group_map[group_name]["contribution"] += ac.contribution
                group_map[group_name]["atom_indices"].append(i)

        functional_group_summary = [
            FunctionalGroupContribution(
                group_name=name,
                contribution=round(data["contribution"], 4),
                atom_indices=data["atom_indices"],
            )
            for name, data in sorted(
                group_map.items(), key=lambda x: x[1]["contribution"], reverse=True
            )
        ]

        total = round(sum(c.contribution for c in atom_contributions), 4)

        return TPSABreakdownResult(
            total=total,
            atom_contributions=atom_contributions,
            functional_group_summary=functional_group_summary,
        )

    def calculate_logp_breakdown(self, mol: Chem.Mol) -> LogPBreakdownResult:
        """
        Calculate LogP per-atom breakdown using Crippen contributions.

        CRITICAL: Must call Chem.AddHs(mol) before getting atom contribs,
        then fold H contributions back to parent heavy atoms.

        Args:
            mol: RDKit molecule object

        Returns:
            LogPBreakdownResult with per-atom and per-group contributions
        """
        mol_h = Chem.AddHs(mol)
        contribs = Crippen._GetAtomContribs(mol_h)

        num_heavy = mol.GetNumAtoms()

        # Initialize heavy atom contributions
        heavy_contribs = [0.0] * num_heavy

        for i in range(mol_h.GetNumAtoms()):
            atom = mol_h.GetAtomWithIdx(i)
            logp_contrib = contribs[i][0]  # (logP, MR) tuple

            if atom.GetAtomicNum() == 1:
                # Hydrogen: fold into its neighbor (the heavy atom)
                neighbors = atom.GetNeighbors()
                if neighbors:
                    nbr_idx = neighbors[0].GetIdx()
                    if nbr_idx < num_heavy:
                        heavy_contribs[nbr_idx] += logp_contrib
            else:
                if i < num_heavy:
                    heavy_contribs[i] += logp_contrib

        atom_contributions = []
        for i in range(num_heavy):
            atom = mol.GetAtomWithIdx(i)
            atom_contributions.append(
                AtomContribution(
                    atom_index=i,
                    symbol=atom.GetSymbol(),
                    contribution=round(heavy_contribs[i], 4),
                )
            )

        # Build functional group summary
        group_map: dict = {}
        for i, ac in enumerate(atom_contributions):
            atom = mol.GetAtomWithIdx(i)
            group_name = _classify_logp_group(atom, mol)
            if group_name not in group_map:
                group_map[group_name] = {"contribution": 0.0, "atom_indices": []}
            group_map[group_name]["contribution"] += ac.contribution
            group_map[group_name]["atom_indices"].append(i)

        functional_group_summary = [
            FunctionalGroupContribution(
                group_name=name,
                contribution=round(data["contribution"], 4),
                atom_indices=data["atom_indices"],
            )
            for name, data in sorted(
                group_map.items(), key=lambda x: x[1]["contribution"], reverse=True
            )
        ]

        total = round(sum(c.contribution for c in atom_contributions), 4)

        return LogPBreakdownResult(
            total=total,
            atom_contributions=atom_contributions,
            functional_group_summary=functional_group_summary,
        )

    def calculate_bertz_detail(self, mol: Chem.Mol) -> BertzDetailResult:
        """
        Calculate Bertz complexity with graph-level detail.

        Args:
            mol: RDKit molecule object

        Returns:
            BertzDetailResult with complexity metrics
        """
        bertz_ct = GraphDescriptors.BertzCT(mol)
        num_bonds = mol.GetNumBonds()
        num_atoms = mol.GetNumAtoms()
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        num_aromatic_rings = Descriptors.NumAromaticRings(mol)

        # Ring complexity: fraction of bonds in rings
        ring_bond_count = sum(
            1 for bond in mol.GetBonds() if bond.IsInRing()
        )
        ring_complexity = ring_bond_count / num_bonds if num_bonds > 0 else 0.0

        # Interpretation based on Bertz CT ranges
        if bertz_ct < 200:
            interpretation = "Simple molecule"
        elif bertz_ct < 500:
            interpretation = "Moderate complexity"
        elif bertz_ct < 1000:
            interpretation = "Complex molecule"
        else:
            interpretation = "Highly complex molecule"

        return BertzDetailResult(
            bertz_ct=round(bertz_ct, 2),
            num_bonds=num_bonds,
            num_atoms=num_atoms,
            num_rings=num_rings,
            num_aromatic_rings=num_aromatic_rings,
            ring_complexity=round(ring_complexity, 4),
            interpretation=interpretation,
        )

    def calculate_fsp3_detail(self, mol: Chem.Mol) -> Fsp3DetailResult:
        """
        Calculate Fsp3 with per-carbon hybridization enumeration.

        Args:
            mol: RDKit molecule object

        Returns:
            Fsp3DetailResult with per-carbon hybridization data
        """
        fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)

        per_carbon = []
        sp3_count = 0
        sp2_count = 0
        sp_count = 0
        total_carbons = 0

        hyb_map = {
            Chem.rdchem.HybridizationType.SP3: "sp3",
            Chem.rdchem.HybridizationType.SP2: "sp2",
            Chem.rdchem.HybridizationType.SP: "sp",
        }

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon
                total_carbons += 1
                hyb = atom.GetHybridization()
                hyb_str = hyb_map.get(hyb, "other")

                if hyb_str == "sp3":
                    sp3_count += 1
                elif hyb_str == "sp2":
                    sp2_count += 1
                elif hyb_str == "sp":
                    sp_count += 1

                per_carbon.append(
                    CarbonHybridization(
                        atom_index=atom.GetIdx(),
                        symbol="C",
                        hybridization=hyb_str,
                    )
                )

        # Interpretation
        if fsp3 < 0.25:
            interpretation = "Flat molecule (low 3D character)"
        elif fsp3 < 0.5:
            interpretation = "Moderate 3D character"
        else:
            interpretation = "High 3D character"

        return Fsp3DetailResult(
            fsp3=round(fsp3, 4),
            total_carbons=total_carbons,
            sp3_count=sp3_count,
            sp2_count=sp2_count,
            sp_count=sp_count,
            per_carbon=per_carbon,
            interpretation=interpretation,
        )


# Module-level instance and convenience functions
_scorer = PropertyBreakdownScorer()


def calculate_tpsa_breakdown(mol: Chem.Mol) -> TPSABreakdownResult:
    """Calculate TPSA per-atom breakdown."""
    return _scorer.calculate_tpsa_breakdown(mol)


def calculate_logp_breakdown(mol: Chem.Mol) -> LogPBreakdownResult:
    """Calculate LogP per-atom breakdown with H folding."""
    return _scorer.calculate_logp_breakdown(mol)


def calculate_bertz_detail(mol: Chem.Mol) -> BertzDetailResult:
    """Calculate Bertz complexity with graph-level detail."""
    return _scorer.calculate_bertz_detail(mol)


def calculate_fsp3_detail(mol: Chem.Mol) -> Fsp3DetailResult:
    """Calculate Fsp3 with per-carbon hybridization."""
    return _scorer.calculate_fsp3_detail(mol)
