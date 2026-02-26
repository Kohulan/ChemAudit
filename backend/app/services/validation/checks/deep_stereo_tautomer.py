"""
Deep Stereo/Tautomer Validation Checks (M1.1)

Implements DVAL-01 through DVAL-05:
- DVAL-01/02: Stereoisomer enumeration (undefined centers + enumerated SMILES)
- DVAL-03:    Tautomer detection (canonical form detection + count)
- DVAL-04:    Aromatic system validation (unusual ring sizes + charged aromatics)
- DVAL-05:    Coordinate dimension detection (2D/3D/no_coordinates/degenerate)
"""

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.MolStandardize import rdMolStandardize

from app.schemas.common import Severity

from ..registry import CheckRegistry
from .base import BaseCheck, CheckResult


def _get_largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    """
    Return the largest fragment of a molecule by heavy atom count.

    Args:
        mol: RDKit molecule object (possibly multi-fragment)

    Returns:
        Largest fragment as a new Mol object
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not frags:
        return mol
    return max(frags, key=lambda f: f.GetNumAtoms())


@CheckRegistry.register("stereoisomer_enumeration")
class StereoisomerEnumerationCheck(BaseCheck):
    """
    Enumerate all possible stereoisomers for undefined stereocenters.

    Covers DVAL-01 (undefined stereocenters) and DVAL-02 (stereoisomer enumeration)
    as a single operation: detects undefined tetrahedral chiral centers and returns
    all enumerated stereoisomer SMILES up to a configurable cap.
    """

    name = "stereoisomer_enumeration"
    description = "Enumerate stereoisomers for undefined stereocenters"
    category = "stereo_tautomer"

    ENUMERATION_CAP = 128  # 2^7 — return count-only above this

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Find undefined stereocenters and enumerate possible stereoisomers.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with WARNING if undefined centers found, INFO if all defined,
            affected_atoms listing undefined center indices, and details with
            undefined_count, total_centers, atom_indices, stereoisomer_smiles,
            enumeration_cap, and cap_exceeded flag.
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot enumerate stereoisomers for None molecule",
                affected_atoms=[],
                details={},
            )

        try:
            # Use largest fragment to avoid multi-fragment confusion
            work_mol = _get_largest_fragment(mol)

            # Find all chiral centers including unassigned ones
            chiral_centers = Chem.FindMolChiralCenters(
                work_mol, includeUnassigned=True, useLegacyImplementation=False
            )

            # Filter for undefined centers (marked with '?')
            undefined_centers = [
                (idx, stereo) for idx, stereo in chiral_centers if stereo == "?"
            ]

            total_centers = len(chiral_centers)
            undefined_count = len(undefined_centers)
            atom_indices = [idx for idx, _ in undefined_centers]

            if undefined_count == 0:
                # Enumerate the single stereoisomer (the molecule as-is)
                opts = StereoEnumerationOptions(unique=True, onlyStereoGroups=False)
                isomers = list(EnumerateStereoisomers(work_mol, options=opts))
                stereoisomer_smiles = [Chem.MolToSmiles(i) for i in isomers]

                msg = (
                    f"All {total_centers} stereocenter(s) properly defined."
                    if total_centers > 0
                    else "No chiral centers present."
                )
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message=msg,
                    affected_atoms=[],
                    details={
                        "undefined_count": 0,
                        "total_centers": total_centers,
                        "atom_indices": [],
                        "stereoisomer_smiles": stereoisomer_smiles,
                        "enumeration_cap": self.ENUMERATION_CAP,
                        "cap_exceeded": False,
                    },
                )

            # Enumerate stereoisomers for undefined centers
            opts = StereoEnumerationOptions(unique=True, onlyStereoGroups=False)
            isomers = list(EnumerateStereoisomers(work_mol, options=opts))
            isomer_count = len(isomers)

            cap_exceeded = isomer_count > self.ENUMERATION_CAP
            if cap_exceeded:
                stereoisomer_smiles: list[str] = []
            else:
                stereoisomer_smiles = [Chem.MolToSmiles(i) for i in isomers]

            indices_str = ", ".join(str(i) for i in atom_indices)
            msg = (
                f"Found {undefined_count} undefined stereocenter(s) at atom "
                f"indices {indices_str}. {isomer_count} stereoisomer(s) enumerated."
            )

            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.WARNING,
                message=msg,
                affected_atoms=atom_indices,
                details={
                    "undefined_count": undefined_count,
                    "total_centers": total_centers,
                    "atom_indices": atom_indices,
                    "stereoisomer_smiles": stereoisomer_smiles,
                    "enumeration_cap": self.ENUMERATION_CAP,
                    "cap_exceeded": cap_exceeded,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error enumerating stereoisomers: {str(e)}",
                affected_atoms=[],
                details={},
            )


@CheckRegistry.register("tautomer_detection")
class TautomerDetectionCheck(BaseCheck):
    """
    Detect tautomeric forms and identify whether the input is the canonical tautomer.

    Covers DVAL-03. Uses RDKit's MolStandardize TautomerEnumerator (not MolVS).
    For multi-fragment molecules, operates on the largest fragment only.
    """

    name = "tautomer_detection"
    description = "Detect tautomers and identify canonical tautomeric form"
    category = "stereo_tautomer"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Enumerate tautomers and check if input is the canonical tautomer.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with INFO severity (having tautomers is informational),
            details including tautomer_count, canonical_smiles, is_canonical_form,
            and tautomer_smiles list.
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot detect tautomers for None molecule",
                affected_atoms=[],
                details={},
            )

        try:
            # Use largest fragment to avoid multi-fragment enumeration confusion
            work_mol = _get_largest_fragment(mol)

            enumerator = rdMolStandardize.TautomerEnumerator()

            # Get canonical tautomer (Canonicalize, NOT GetCanonicalTautomer)
            canonical = enumerator.Canonicalize(work_mol)
            canonical_smiles = Chem.MolToSmiles(canonical)

            # Enumerate all tautomers
            tautomers = enumerator.Enumerate(work_mol)
            tautomer_count = len(tautomers)
            tautomer_smiles = [Chem.MolToSmiles(t) for t in tautomers]

            # Determine if input is the canonical form
            input_smiles = Chem.MolToSmiles(work_mol)
            is_canonical_form = input_smiles == canonical_smiles

            canonical_str = "Is" if is_canonical_form else "Is not"
            msg = (
                f"Molecule has {tautomer_count} tautomeric form(s). "
                f"{canonical_str} the canonical tautomer."
            )

            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.INFO,
                message=msg,
                affected_atoms=[],
                details={
                    "tautomer_count": tautomer_count,
                    "canonical_smiles": canonical_smiles,
                    "is_canonical_form": is_canonical_form,
                    "tautomer_smiles": tautomer_smiles,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error detecting tautomers: {str(e)}",
                affected_atoms=[],
                details={},
            )


@CheckRegistry.register("aromatic_system_validation")
class AromaticSystemValidationCheck(BaseCheck):
    """
    Validate aromatic systems for unusual ring sizes and charged aromatics.

    Covers DVAL-04. Extends the existing AromaticityCheck (which handles
    Kekulization failures) with ring-size and formal-charge analysis:
    - Unusual aromatic ring sizes: aromatic rings that are not 5 or 6 membered
    - Charged aromatics: aromatic atoms with non-zero formal charge
    """

    name = "aromatic_system_validation"
    description = "Detect unusual aromatic ring sizes and charged aromatic atoms"
    category = "stereo_tautomer"

    STANDARD_AROMATIC_SIZES = {5, 6}

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check for unusual aromatic ring sizes and charged aromatic atoms.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with WARNING if unusual rings or charged aromatics found,
            affected_atoms listing flagged atom indices, and details with
            unusual_ring_sizes and charged_aromatics lists.
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot validate aromatic system of None molecule",
                affected_atoms=[],
                details={},
            )

        try:
            ri = mol.GetRingInfo()
            unusual_rings: list[dict] = []
            charged_aromatics: list[dict] = []
            affected_atom_set: set[int] = set()

            # Check 1: Unusual aromatic ring sizes (not 5 or 6)
            for ring_atoms in ri.AtomRings():
                ring_size = len(ring_atoms)
                # Check if all atoms in ring are aromatic
                all_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms)
                if all_aromatic and ring_size not in self.STANDARD_AROMATIC_SIZES:
                    unusual_rings.append({
                        "ring_size": ring_size,
                        "atom_indices": list(ring_atoms),
                    })
                    affected_atom_set.update(ring_atoms)

            # Check 2: Charged aromatic atoms
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic() and atom.GetFormalCharge() != 0:
                    charged_aromatics.append({
                        "atom_idx": atom.GetIdx(),
                        "symbol": atom.GetSymbol(),
                        "charge": atom.GetFormalCharge(),
                    })
                    affected_atom_set.add(atom.GetIdx())

            affected_atoms = sorted(affected_atom_set)

            if not unusual_rings and not charged_aromatics:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="Aromatic systems have standard ring sizes and no charged atoms.",
                    affected_atoms=[],
                    details={
                        "unusual_ring_sizes": [],
                        "charged_aromatics": [],
                    },
                )

            # Build descriptive message
            parts: list[str] = []
            if unusual_rings:
                sizes = [r["ring_size"] for r in unusual_rings]
                parts.append(f"Found {len(unusual_rings)} unusual aromatic ring(s) with size(s): {sizes}.")
            if charged_aromatics:
                atoms_desc = [
                    f"{c['symbol']}({'+' if c['charge'] > 0 else ''}{c['charge']}) at idx {c['atom_idx']}"
                    for c in charged_aromatics
                ]
                parts.append(f"Found {len(charged_aromatics)} charged aromatic atom(s): {', '.join(atoms_desc)}.")
            msg = " ".join(parts)

            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.WARNING,
                message=msg,
                affected_atoms=affected_atoms,
                details={
                    "unusual_ring_sizes": unusual_rings,
                    "charged_aromatics": charged_aromatics,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error validating aromatic system: {str(e)}",
                affected_atoms=[],
                details={},
            )


@CheckRegistry.register("coordinate_dimension")
class CoordinateDimensionCheck(BaseCheck):
    """
    Detect the coordinate dimensionality of a molecule.

    Covers DVAL-05. Identifies whether a molecule has:
    - 2D coordinates (conformer present, all z == 0.0)
    - 3D coordinates (conformer present, at least one non-zero z)
    - No coordinates (no conformer attached)
    - Degenerate coordinates (all atoms at the same position)

    This is an informational check — always passes.
    """

    name = "coordinate_dimension"
    description = "Detect coordinate dimensionality: 2D, 3D, no_coordinates, or degenerate"
    category = "stereo_tautomer"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Determine the coordinate dimension of the molecule.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with INFO severity and details containing
            dimension (2d/3d/no_coordinates/degenerate) and num_conformers.
            Always passes — this is a purely informational check.
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot detect coordinate dimension of None molecule",
                affected_atoms=[],
                details={},
            )

        try:
            num_conformers = mol.GetNumConformers()

            if num_conformers == 0:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="Molecule has no coordinate conformer (SMILES input).",
                    affected_atoms=[],
                    details={
                        "dimension": "no_coordinates",
                        "num_conformers": 0,
                    },
                )

            conf = mol.GetConformer(0)
            num_atoms = mol.GetNumAtoms()

            if num_atoms == 0:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="Molecule has no atoms; no coordinate dimension to detect.",
                    affected_atoms=[],
                    details={
                        "dimension": "no_coordinates",
                        "num_conformers": num_conformers,
                    },
                )

            # Collect all positions
            positions = [conf.GetAtomPosition(i) for i in range(num_atoms)]

            # Check for degenerate coordinates: all atoms at same (x, y, z)
            xs = {p.x for p in positions}
            ys = {p.y for p in positions}
            zs = {p.z for p in positions}

            if len(xs) == 1 and len(ys) == 1 and len(zs) == 1 and num_atoms > 1:
                return CheckResult(
                    check_name=self.name,
                    passed=True,
                    severity=Severity.INFO,
                    message="Molecule has degenerate coordinates (all atoms at identical position).",
                    affected_atoms=[],
                    details={
                        "dimension": "degenerate",
                        "num_conformers": num_conformers,
                    },
                )

            # Check z-coordinates
            has_nonzero_z = any(p.z != 0.0 for p in positions)

            if has_nonzero_z:
                dimension = "3d"
                msg = "Molecule has 3D coordinates."
            else:
                dimension = "2d"
                msg = "Molecule has 2D coordinates (all z-coordinates are zero)."

            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.INFO,
                message=msg,
                affected_atoms=[],
                details={
                    "dimension": dimension,
                    "num_conformers": num_conformers,
                },
            )

        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error detecting coordinate dimension: {str(e)}",
                affected_atoms=[],
                details={},
            )
