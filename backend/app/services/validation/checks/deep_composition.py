"""
Deep Composition Validation Checks

Milestone 1.2 checks covering chemical composition guards:
mixture detection, solvent contamination, inorganic/organometallic filtering,
radical detection, isotope label detection, and trivial molecule flagging.

Checks:
    - mixture_detection (DVAL-06)
    - solvent_contamination (DVAL-07)
    - inorganic_filter (DVAL-08)
    - radical_detection (DVAL-09)
    - isotope_label_detection (DVAL-10)
    - trivial_molecule (DVAL-11)
"""

from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors

from app.schemas.common import Severity

from ..registry import CheckRegistry
from .base import BaseCheck, CheckResult

# ---------------------------------------------------------------------------
# Common solvents (SMILES-based canonical matching + substructure)
# ---------------------------------------------------------------------------
COMMON_SOLVENTS: Dict[str, str] = {
    "water": "O",
    "DMSO": "CS(=O)C",
    "DMF": "CN(C)C=O",
    "acetonitrile": "CC#N",
    "methanol": "CO",
    "ethanol": "CCO",
    "acetone": "CC(=O)C",
    "THF": "C1CCCO1",
    "DCM": "ClCCl",
    "chloroform": "ClC(Cl)Cl",
    "toluene": "Cc1ccccc1",
    "hexane": "CCCCCC",
    "diethyl_ether": "CCOCC",
    "ethyl_acetate": "CCOC(=O)C",
    "isopropanol": "CC(O)C",
}

# Precompute canonical SMILES and mol objects for fast matching
_SOLVENT_CANON_SMILES: Dict[str, str] = {}
_SOLVENT_MOLS: Dict[str, Optional[Chem.Mol]] = {}
for _name, _smi in COMMON_SOLVENTS.items():
    _m = Chem.MolFromSmiles(_smi)
    if _m:
        _SOLVENT_CANON_SMILES[_name] = Chem.MolToSmiles(_m)
        _SOLVENT_MOLS[_name] = _m
    else:
        _SOLVENT_CANON_SMILES[_name] = _smi
        _SOLVENT_MOLS[_name] = None

# MolVS patterns categorised as solvents by name
_MOLVS_SOLVENT_KEYWORDS = {"water", "hydroxide", "ammonia", "ammonium"}

# ---------------------------------------------------------------------------
# Metal atomic numbers for organometallic detection
# ---------------------------------------------------------------------------
METAL_ATOMIC_NUMS = {
    # Alkali metals
    3,   # Li
    11,  # Na
    19,  # K
    37,  # Rb
    55,  # Cs
    # Alkaline earth metals
    4,   # Be
    12,  # Mg
    20,  # Ca
    38,  # Sr
    56,  # Ba
    # Transition metals
    21,  # Sc
    22,  # Ti
    23,  # V
    24,  # Cr
    25,  # Mn
    26,  # Fe
    27,  # Co
    28,  # Ni
    29,  # Cu
    30,  # Zn
    39,  # Y
    40,  # Zr
    41,  # Nb
    42,  # Mo
    43,  # Tc
    44,  # Ru
    45,  # Rh
    46,  # Pd
    47,  # Ag
    48,  # Cd
    72,  # Hf
    73,  # Ta
    74,  # W
    75,  # Re
    76,  # Os
    77,  # Ir
    78,  # Pt
    79,  # Au
    80,  # Hg
    # Post-transition metals
    13,  # Al
    31,  # Ga
    49,  # In
    50,  # Sn
    81,  # Tl
    82,  # Pb
    83,  # Bi
    # Lanthanides
    57,  # La
    58,  # Ce
    59,  # Pr
    60,  # Nd
    61,  # Pm
    62,  # Sm
    63,  # Eu
    64,  # Gd
    65,  # Tb
    66,  # Dy
    67,  # Ho
    68,  # Er
    69,  # Tm
    70,  # Yb
    71,  # Lu
    # Actinides (common)
    89,  # Ac
    90,  # Th
    92,  # U
    94,  # Pu
}

# ---------------------------------------------------------------------------
# Isotope common name mapping: (atomic_num, isotope_mass) → common_name
# ---------------------------------------------------------------------------
ISOTOPE_COMMON_NAMES: Dict[Tuple[int, int], str] = {
    (1, 2): "deuterium",
    (1, 3): "tritium",
    (6, 13): "carbon-13",
    (6, 14): "carbon-14",
    (7, 15): "nitrogen-15",
    (8, 17): "oxygen-17",
    (8, 18): "oxygen-18",
    (9, 18): "fluorine-18",
    (15, 32): "phosphorus-32",
    (16, 35): "sulfur-35",
}


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _fragment_mol_to_smiles(frag: Chem.Mol) -> str:
    """Return canonical SMILES for a fragment, or empty string on failure."""
    try:
        return Chem.MolToSmiles(frag)
    except Exception:
        return ""


def _fragment_mw(frag: Chem.Mol) -> float:
    """Return molecular weight of fragment."""
    try:
        return round(Descriptors.MolWt(frag), 4)
    except Exception:
        return 0.0


def _classify_fragment(
    frag_mol: Chem.Mol,
    frag_smiles: str,
    frag_mw: float,
    is_largest: bool,
) -> Tuple[str, Optional[str]]:
    """
    Classify a fragment as 'drug', 'salt', 'solvent', or 'unknown'.

    Returns (classification, pattern_name) where pattern_name is the matched
    MolVS pattern name (or None if no match).
    """
    try:
        from molvs.fragment import REMOVE_FRAGMENTS  # type: ignore[import]

        frag_heavy_count = frag_mol.GetNumHeavyAtoms()
        for fp in REMOVE_FRAGMENTS:
            try:
                pattern_mol = fp.smarts
                # Only match if the fragment has the same heavy atom count as the pattern
                # to avoid false positives from broad SMARTS (e.g. [#7] matching all amines)
                pattern_heavy = pattern_mol.GetNumHeavyAtoms() if pattern_mol else -1
                if pattern_heavy != frag_heavy_count:
                    continue
                matches = frag_mol.HasSubstructMatch(pattern_mol)
            except Exception:
                try:
                    smarts_mol = Chem.MolFromSmarts(fp.smarts_str)
                    if smarts_mol is None:
                        continue
                    if smarts_mol.GetNumHeavyAtoms() != frag_heavy_count:
                        continue
                    matches = frag_mol.HasSubstructMatch(smarts_mol)
                except Exception:
                    matches = False

            if matches:
                name_lower = fp.name.lower()
                if any(kw in name_lower for kw in _MOLVS_SOLVENT_KEYWORDS):
                    return "solvent", fp.name
                return "salt", fp.name
    except ImportError:
        pass

    # Heuristic fallback
    has_carbon = any(a.GetAtomicNum() == 6 for a in frag_mol.GetAtoms())

    # Largest carbon-containing fragment defaults to 'drug'
    if is_largest and has_carbon:
        return "drug", None
    if has_carbon and frag_mw > 100:
        return "drug", None
    if not has_carbon or frag_mw < 50:
        return "salt", None
    return "unknown", None


# ---------------------------------------------------------------------------
# DVAL-06: Mixture Detection
# ---------------------------------------------------------------------------

@CheckRegistry.register("mixture_detection")
class MixtureDetectionCheck(BaseCheck):
    """
    Detect multi-fragment mixtures and classify each fragment.

    DVAL-06: Input containing a dot-separated SMILES (mixture) is flagged
    with each fragment classified as drug, salt, solvent, or unknown,
    accompanied by SMILES, MW, heavy atom count, and matched pattern name.
    """

    name = "mixture_detection"
    description = "Detect multi-fragment mixtures and classify each component"
    category = "chemical_composition"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check for multiple disconnected fragments and classify each.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult — passes for single fragments, WARNING for mixtures
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check mixture composition of None molecule",
            )

        try:
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
            num_frags = len(frags)

            if num_frags == 1:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="Input is a single connected molecule",
                    details={"num_fragments": 1},
                )

            # Sort fragments: largest (by heavy atom count) first
            frags_sorted = sorted(frags, key=lambda f: f.GetNumHeavyAtoms(), reverse=True)

            fragment_infos: List[Dict[str, Any]] = []
            classification_counts: Dict[str, int] = {}

            for i, frag in enumerate(frags_sorted):
                smiles = _fragment_mol_to_smiles(frag)
                mw = _fragment_mw(frag)
                heavy_count = frag.GetNumHeavyAtoms()
                is_largest = i == 0

                classification, pattern_name = _classify_fragment(frag, smiles, mw, is_largest)

                fragment_infos.append({
                    "smiles": smiles,
                    "molecular_weight": mw,
                    "heavy_atom_count": heavy_count,
                    "classification": classification,
                    "pattern_name": pattern_name,
                })
                classification_counts[classification] = classification_counts.get(classification, 0) + 1

            # Build human-readable summary
            summary_parts = [f"{v} {k}" for k, v in classification_counts.items()]
            summary = ", ".join(summary_parts)

            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.WARNING,
                message=f"Input contains {num_frags} fragments: {summary}.",
                affected_atoms=[],
                details={
                    "num_fragments": num_frags,
                    "fragments": fragment_infos,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error checking mixture composition: {str(e)}",
            )


# ---------------------------------------------------------------------------
# DVAL-07: Solvent Contamination
# ---------------------------------------------------------------------------

@CheckRegistry.register("solvent_contamination")
class SolventContaminationCheck(BaseCheck):
    """
    Detect common lab solvents contaminating the molecule input.

    DVAL-07: A molecule containing a known solvent fragment (DMSO, DMF,
    water, etc.) returns a solvent contamination warning with the solvent name.
    """

    name = "solvent_contamination"
    description = "Detect common laboratory solvents in the input"
    category = "chemical_composition"

    def _is_solvent_fragment(self, frag_mol: Chem.Mol) -> Optional[Dict[str, Any]]:
        """
        Check if a fragment matches a known solvent.

        Returns a dict with solvent info if matched, else None.
        """
        frag_smiles = _fragment_mol_to_smiles(frag_mol)

        # Check against COMMON_SOLVENTS via canonical SMILES comparison
        for name, canon_smi in _SOLVENT_CANON_SMILES.items():
            if frag_smiles == canon_smi:
                return {
                    "name": name,
                    "smiles": frag_smiles,
                    "molecular_weight": _fragment_mw(frag_mol),
                }

        # Substructure match against known solvents
        for name, sol_mol in _SOLVENT_MOLS.items():
            if sol_mol is None:
                continue
            try:
                if (frag_mol.GetNumHeavyAtoms() == sol_mol.GetNumHeavyAtoms()
                        and frag_mol.HasSubstructMatch(sol_mol)
                        and sol_mol.HasSubstructMatch(frag_mol)):
                    return {
                        "name": name,
                        "smiles": frag_smiles,
                        "molecular_weight": _fragment_mw(frag_mol),
                    }
            except Exception:
                continue

        # Check MolVS REMOVE_FRAGMENTS for water/ammonia
        # Only match if the fragment has the same heavy atom count as the pattern
        # to avoid false positives from broad SMARTS like [#8] or [#7]
        try:
            from molvs.fragment import REMOVE_FRAGMENTS  # type: ignore[import]

            frag_heavy_count = frag_mol.GetNumHeavyAtoms()
            for fp in REMOVE_FRAGMENTS:
                name_lower = fp.name.lower()
                if not any(kw in name_lower for kw in _MOLVS_SOLVENT_KEYWORDS):
                    continue
                try:
                    pattern_mol = fp.smarts
                    pattern_heavy = pattern_mol.GetNumHeavyAtoms() if pattern_mol else 0
                    # Only match if fragment has the same heavy atom count as pattern
                    if pattern_heavy != frag_heavy_count:
                        continue
                    if frag_mol.HasSubstructMatch(pattern_mol):
                        return {
                            "name": fp.name,
                            "smiles": frag_smiles,
                            "molecular_weight": _fragment_mw(frag_mol),
                        }
                except Exception:
                    pass
        except ImportError:
            pass

        return None

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check all fragments for common laboratory solvents.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult — passes if no solvents detected, WARNING if found
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check solvent contamination of None molecule",
            )

        try:
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
            num_frags = len(frags)
            solvents_found: List[Dict[str, Any]] = []

            for frag in frags:
                match = self._is_solvent_fragment(frag)
                if match:
                    # Avoid duplicate entries by name
                    if not any(s["name"] == match["name"] for s in solvents_found):
                        solvents_found.append(match)

            if not solvents_found:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="No solvent contamination detected",
                    details={"solvents_found": [], "is_pure_solvent": False},
                )

            is_pure_solvent = num_frags == 1
            solvent_names = ", ".join(s["name"] for s in solvents_found)

            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.WARNING,
                message=f"Detected solvent contaminant(s): {solvent_names}.",
                affected_atoms=[],
                details={
                    "solvents_found": solvents_found,
                    "is_pure_solvent": is_pure_solvent,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error checking solvent contamination: {str(e)}",
            )


# ---------------------------------------------------------------------------
# DVAL-08: Inorganic / Organometallic Filter
# ---------------------------------------------------------------------------

@CheckRegistry.register("inorganic_filter")
class InorganicFilterCheck(BaseCheck):
    """
    Flag inorganic molecules and organometallics.

    DVAL-08: Molecules with no carbon atoms (inorganic) or that contain
    metal atoms alongside carbon (organometallic) are flagged with atom details.
    """

    name = "inorganic_filter"
    description = "Detect inorganic or organometallic molecules"
    category = "chemical_composition"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check for absence of carbon (inorganic) or presence of metals (organometallic).

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult — WARNING for organometallic, ERROR for purely inorganic
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check inorganic/organometallic nature of None molecule",
            )

        try:
            has_carbon = False
            metal_atoms: List[Dict[str, Any]] = []
            element_counts: Dict[str, int] = {}

            for atom in mol.GetAtoms():
                atomic_num = atom.GetAtomicNum()
                symbol = atom.GetSymbol()

                if atomic_num == 6:
                    has_carbon = True

                element_counts[symbol] = element_counts.get(symbol, 0) + 1

                if atomic_num in METAL_ATOMIC_NUMS:
                    metal_atoms.append({
                        "atom_idx": atom.GetIdx(),
                        "symbol": symbol,
                        "atomic_num": atomic_num,
                    })

            is_inorganic = not has_carbon
            is_organometallic = has_carbon and len(metal_atoms) > 0

            if not is_inorganic and not is_organometallic:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="Molecule is organic (contains carbon, no metals)",
                    details={
                        "has_carbon": has_carbon,
                        "is_inorganic": False,
                        "is_organometallic": False,
                        "metal_atoms": [],
                        "element_counts": element_counts,
                    },
                )

            affected_atoms = [m["atom_idx"] for m in metal_atoms]
            element_summary = ", ".join(
                f"{sym}({cnt})" for sym, cnt in element_counts.items()
            )

            if is_organometallic:
                metal_symbols = ", ".join(
                    sorted({m["symbol"] for m in metal_atoms})
                )
                return CheckResult(
                    check_name=self.name,
                    passed=False,
                    severity=Severity.WARNING,
                    message=(
                        f"Molecule is organometallic: contains metal(s) {metal_symbols} "
                        f"alongside carbon. Elements: {element_summary}."
                    ),
                    affected_atoms=affected_atoms,
                    details={
                        "has_carbon": True,
                        "is_inorganic": False,
                        "is_organometallic": True,
                        "metal_atoms": metal_atoms,
                        "element_counts": element_counts,
                    },
                )

            # Purely inorganic
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Molecule is inorganic: no carbon atoms. Elements: {element_summary}.",
                affected_atoms=affected_atoms,
                details={
                    "has_carbon": False,
                    "is_inorganic": True,
                    "is_organometallic": False,
                    "metal_atoms": metal_atoms,
                    "element_counts": element_counts,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error checking inorganic/organometallic nature: {str(e)}",
            )


# ---------------------------------------------------------------------------
# DVAL-09: Radical Detection
# ---------------------------------------------------------------------------

@CheckRegistry.register("radical_detection")
class RadicalDetectionCheck(BaseCheck):
    """
    Detect atoms bearing radical electrons.

    DVAL-09: A molecule with radical electrons is flagged with the affected
    atom indices and number of radical electrons per atom.
    """

    name = "radical_detection"
    description = "Detect atoms with radical (unpaired) electrons"
    category = "chemical_composition"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Iterate all atoms and flag those with non-zero radical electron counts.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult — passes if no radicals, WARNING if found
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check radicals of None molecule",
            )

        try:
            radical_atoms: List[Dict[str, Any]] = []

            for atom in mol.GetAtoms():
                num_radical = atom.GetNumRadicalElectrons()
                if num_radical > 0:
                    radical_atoms.append({
                        "atom_idx": atom.GetIdx(),
                        "symbol": atom.GetSymbol(),
                        "num_radical_electrons": num_radical,
                    })

            if not radical_atoms:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="No radical electrons detected",
                    details={"radical_atoms": [], "total_radical_electrons": 0},
                )

            total_radical = sum(r["num_radical_electrons"] for r in radical_atoms)
            affected_atoms = [r["atom_idx"] for r in radical_atoms]

            atom_descriptions = ", ".join(
                f"{r['symbol']}(idx={r['atom_idx']}, radicals={r['num_radical_electrons']})"
                for r in radical_atoms
            )

            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.WARNING,
                message=f"Detected {total_radical} radical electron(s) on {atom_descriptions}.",
                affected_atoms=affected_atoms,
                details={
                    "radical_atoms": radical_atoms,
                    "total_radical_electrons": total_radical,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error checking radical electrons: {str(e)}",
            )


# ---------------------------------------------------------------------------
# DVAL-10: Isotope Label Detection
# ---------------------------------------------------------------------------

@CheckRegistry.register("isotope_label_detection")
class IsotopeLabelDetectionCheck(BaseCheck):
    """
    Detect isotope-labeled atoms.

    DVAL-10: A molecule with isotope labels (e.g., deuterium, 13C) is flagged
    with isotope details per atom, including common isotope names.
    """

    name = "isotope_label_detection"
    description = "Detect isotope-labeled atoms (deuterium, 13C, etc.)"
    category = "chemical_composition"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Find all atoms with non-default isotope labels.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult — passes if no labels (INFO), INFO severity if labels found
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check isotope labels of None molecule",
            )

        try:
            labeled_atoms: List[Dict[str, Any]] = []

            for atom in mol.GetAtoms():
                isotope = atom.GetIsotope()
                if isotope > 0:
                    atomic_num = atom.GetAtomicNum()
                    symbol = atom.GetSymbol()
                    common_name = ISOTOPE_COMMON_NAMES.get((atomic_num, isotope))

                    labeled_atoms.append({
                        "atom_idx": atom.GetIdx(),
                        "symbol": symbol,
                        "isotope_mass": isotope,
                        "common_name": common_name,
                    })

            if not labeled_atoms:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="No isotope labels detected",
                    details={"labeled_atoms": [], "total_labeled": 0},
                )

            total_labeled = len(labeled_atoms)
            affected_atoms = [la["atom_idx"] for la in labeled_atoms]

            summary_parts = []
            for la in labeled_atoms:
                label_str = la["common_name"] if la["common_name"] else f"{la['isotope_mass']}{la['symbol']}"
                summary_parts.append(label_str)

            summary = ", ".join(summary_parts)

            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.INFO,
                message=f"Found {total_labeled} isotope-labeled atom(s): {summary}.",
                affected_atoms=affected_atoms,
                details={
                    "labeled_atoms": labeled_atoms,
                    "total_labeled": total_labeled,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error checking isotope labels: {str(e)}",
            )


# ---------------------------------------------------------------------------
# DVAL-11: Trivial Molecule Check
# ---------------------------------------------------------------------------

@CheckRegistry.register("trivial_molecule")
class TrivialMoleculeCheck(BaseCheck):
    """
    Flag trivially small molecules as too small for meaningful validation.

    DVAL-11: A molecule with heavy atom count <= 3 is flagged as too small
    for reliable chemical validation and ML scoring.
    """

    name = "trivial_molecule"
    description = "Flag molecules that are too small for meaningful validation"
    category = "chemical_composition"

    HEAVY_ATOM_THRESHOLD = 3  # Aligned with ChEMBL fragment filter

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check that the molecule has enough heavy atoms for meaningful validation.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult — passes if > HEAVY_ATOM_THRESHOLD, ERROR if too small
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check trivial molecule status of None molecule",
            )

        try:
            heavy_atom_count = mol.GetNumHeavyAtoms()
            num_bonds = mol.GetNumBonds()
            is_single_atom = num_bonds == 0 and heavy_atom_count <= 1

            if heavy_atom_count > self.HEAVY_ATOM_THRESHOLD:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message=(
                        f"Molecule has {heavy_atom_count} heavy atom(s) — "
                        f"sufficient for validation"
                    ),
                    details={
                        "heavy_atom_count": heavy_atom_count,
                        "num_bonds": num_bonds,
                        "is_single_atom": is_single_atom,
                        "threshold": self.HEAVY_ATOM_THRESHOLD,
                    },
                )

            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=(
                    f"Molecule has only {heavy_atom_count} heavy atom(s) — "
                    f"too small for meaningful validation."
                ),
                affected_atoms=[],
                details={
                    "heavy_atom_count": heavy_atom_count,
                    "num_bonds": num_bonds,
                    "is_single_atom": is_single_atom,
                    "threshold": self.HEAVY_ATOM_THRESHOLD,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error checking trivial molecule: {str(e)}",
            )
