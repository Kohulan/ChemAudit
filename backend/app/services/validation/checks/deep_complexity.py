"""
Deep Validation Checks: Structural Complexity (M1.3)

Covers DVAL-12 through DVAL-17:
- DVAL-12: Hypervalent atom detection
- DVAL-13: Polymer/high-MW detection
- DVAL-14: Ring strain (3/4-membered rings)
- DVAL-15: Macrocycle detection (>12 atoms)
- DVAL-16: Charged species and zwitterion detection
- DVAL-17: Explicit hydrogen audit
"""

from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from app.schemas.common import Severity

from ..registry import CheckRegistry
from .base import BaseCheck, CheckResult


@CheckRegistry.register("hypervalent_atoms")
class HypervalentAtomCheck(BaseCheck):
    """
    Detect atoms with valence exceeding their allowed maximum.

    Covers DVAL-12. Uses RDKit's periodic table to determine allowed valences
    per element. Flags atoms where actual total valence exceeds the maximum
    allowed value (excluding elements with unrestricted valence such as
    transition metals where -1 is in the allowed list).

    Severity: WARNING — hypervalent atoms are chemically suspect but may
    represent valid exotic species (e.g., phosphorus, sulfur in expanded valence).
    """

    name = "hypervalent_atoms"
    description = "Detect atoms exceeding their allowed valence"
    category = "structural_complexity"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check for hypervalent atoms in the molecule.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with WARNING severity if hypervalent atoms found
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check hypervalent atoms of None molecule",
            )

        pt = Chem.GetPeriodicTable()
        hypervalent_atoms = []

        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            # Skip dummy atoms (R groups, attachment points)
            if atomic_num == 0:
                continue

            try:
                symbol = atom.GetSymbol()
                allowed_valences = list(pt.GetValenceList(symbol))

                # -1 in allowed valences means unrestricted (e.g., transition metals)
                if -1 in allowed_valences:
                    continue

                if not allowed_valences:
                    continue

                actual_valence = atom.GetTotalValence()
                max_allowed = max(allowed_valences)

                if actual_valence > max_allowed:
                    hypervalent_atoms.append(
                        {
                            "atom_idx": atom.GetIdx(),
                            "symbol": symbol,
                            "actual_valence": actual_valence,
                            "allowed_valences": allowed_valences,
                        }
                    )
            except Exception:
                # Some exotic atoms may fail GetValenceList; skip silently
                continue

        if not hypervalent_atoms:
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.INFO,
                message="No hypervalent atoms detected",
                details={"hypervalent_atoms": []},
            )

        n = len(hypervalent_atoms)
        summary_parts = [
            f"{a['symbol']}({a['actual_valence']})" for a in hypervalent_atoms[:5]
        ]
        summary = ", ".join(summary_parts)
        if n > 5:
            summary += f", ... (+{n - 5} more)"

        return CheckResult(
            check_name=self.name,
            passed=False,
            severity=Severity.WARNING,
            message=f"Found {n} hypervalent atom(s): {summary}.",
            affected_atoms=[a["atom_idx"] for a in hypervalent_atoms],
            details={"hypervalent_atoms": hypervalent_atoms},
        )


@CheckRegistry.register("polymer_detection")
class PolymerDetectionCheck(BaseCheck):
    """
    Detect possible polymer or macromolecule inputs.

    Covers DVAL-13. Uses three heuristics:
    1. SGroup markers (SRU / COP type) from V3000 format
    2. Molecular weight > 1500 Da heuristic
    3. Dummy atoms (atomic_num == 0) indicating attachment points

    Severity: INFO — polymers and high-MW natural products are not errors,
    but users should be aware that validation may be less meaningful.
    """

    name = "polymer_detection"
    description = "Detect possible polymer or macromolecule structures"
    category = "structural_complexity"

    # MW threshold above which a molecule is considered a possible polymer
    MW_THRESHOLD = 1500.0

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check for polymer indicators.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with INFO severity if polymer indicators found
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check polymer indicators of None molecule",
            )

        # --- Path 1: SGroup markers ---
        has_sgroup_markers = False
        sgroup_types: List[str] = []
        try:
            sgroups = Chem.GetMolSubstanceGroups(mol)
            for sg in sgroups:
                sg_type = sg.GetProp("TYPE") if sg.HasProp("TYPE") else ""
                sgroup_types.append(sg_type)
                # SRU = Structural Repeat Unit, COP = Copolymer
                if "SRU" in sg_type.upper() or "COP" in sg_type.upper():
                    has_sgroup_markers = True
        except Exception:
            pass

        # --- Path 2: Molecular weight heuristic ---
        try:
            molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
        except Exception:
            molecular_weight = 0.0
        exceeds_mw_threshold = molecular_weight > self.MW_THRESHOLD

        # --- Path 3: Dummy atoms (attachment points) ---
        dummy_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
        has_dummy_atoms = len(dummy_atoms) > 0

        details: Dict[str, Any] = {
            "has_sgroup_markers": has_sgroup_markers,
            "sgroup_types": sgroup_types,
            "molecular_weight": round(molecular_weight, 4),
            "exceeds_mw_threshold": exceeds_mw_threshold,
            "has_dummy_atoms": has_dummy_atoms,
            "dummy_atom_count": len(dummy_atoms),
        }

        is_flagged = has_sgroup_markers or exceeds_mw_threshold or has_dummy_atoms

        if not is_flagged:
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.INFO,
                message="No polymer indicators detected",
                details=details,
            )

        reasons = []
        if has_sgroup_markers:
            reasons.append(f"SGroup markers present (types: {', '.join(sgroup_types)})")
        if exceeds_mw_threshold:
            reasons.append(f"MW {molecular_weight:.1f} Da exceeds {self.MW_THRESHOLD} Da threshold")
        if has_dummy_atoms:
            reasons.append(f"{len(dummy_atoms)} dummy atom(s) indicating attachment points")

        reason_str = "; ".join(reasons)

        return CheckResult(
            check_name=self.name,
            passed=False,
            severity=Severity.INFO,
            message=f"Possible polymer detected: {reason_str}.",
            affected_atoms=dummy_atoms,
            details=details,
        )


@CheckRegistry.register("ring_strain")
class RingStrainCheck(BaseCheck):
    """
    Detect strained small rings (3 or 4-membered).

    Covers DVAL-14. Uses ring size heuristic only (no force field calculation).
    3-membered rings (cyclopropane) and 4-membered rings (cyclobutane) have
    significant angle strain affecting stability and reactivity.

    Severity: WARNING — ring strain can affect compound stability and
    ML model applicability.
    """

    name = "ring_strain"
    description = "Detect 3 or 4-membered strained rings"
    category = "structural_complexity"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check for 3- and 4-membered rings.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with WARNING severity if strained rings found
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check ring strain of None molecule",
            )

        try:
            ring_info = mol.GetRingInfo()
            atom_rings = ring_info.AtomRings()
        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error retrieving ring information: {str(e)}",
            )

        strained_rings = []
        for ring in atom_rings:
            ring_size = len(ring)
            if ring_size in (3, 4):
                strained_rings.append(
                    {
                        "ring_size": ring_size,
                        "atom_indices": list(ring),
                    }
                )

        details = {
            "strained_rings": strained_rings,
            "total_strained_rings": len(strained_rings),
        }

        if not strained_rings:
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.INFO,
                message="No strained rings detected",
                details=details,
            )

        # Collect all affected atom indices (union across all strained rings)
        affected = set()
        for ring in strained_rings:
            affected.update(ring["atom_indices"])

        n = len(strained_rings)
        sizes = [str(r["ring_size"]) for r in strained_rings]
        sizes_str = ", ".join(sizes)

        return CheckResult(
            check_name=self.name,
            passed=False,
            severity=Severity.WARNING,
            message=f"Found {n} strained ring(s): {sizes_str}-membered ring(s).",
            affected_atoms=sorted(affected),
            details=details,
        )


@CheckRegistry.register("macrocycle_detection")
class MacrocycleDetectionCheck(BaseCheck):
    """
    Detect macrocyclic rings (> 12 atoms).

    Covers DVAL-15. Uses SSSR (Smallest Set of Smallest Rings) from RDKit.
    Note: SSSR has a limitation for fused bicyclic systems — the peripheral
    macrocyclic path may not be captured. This is documented in the output.

    Severity: INFO — macrocycles are valid drug candidates (cyclic peptides,
    macrolide antibiotics). This is informational context for the user.
    """

    name = "macrocycle_detection"
    description = "Detect macrocyclic rings with more than 12 atoms"
    category = "structural_complexity"

    # Threshold: rings with more than this many atoms are macrocycles
    MACROCYCLE_THRESHOLD = 12

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check for macrocyclic rings (> 12 atoms).

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with INFO severity if macrocycles found
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check macrocycles of None molecule",
            )

        try:
            ring_info = mol.GetRingInfo()
            atom_rings = ring_info.AtomRings()
        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Error retrieving ring information: {str(e)}",
            )

        macrocycles = []
        for ring in atom_rings:
            ring_size = len(ring)
            if ring_size > self.MACROCYCLE_THRESHOLD:
                macrocycles.append(
                    {
                        "ring_size": ring_size,
                        "atom_indices": list(ring),
                    }
                )

        details = {
            "macrocycles": macrocycles,
            "total_macrocycles": len(macrocycles),
            "sssr_note": (
                "Uses SSSR; fused bicyclic peripheral paths may not be detected"
            ),
        }

        if not macrocycles:
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.INFO,
                message="No macrocyclic rings detected",
                details=details,
            )

        # Collect all affected atom indices (union across all macrocycles)
        affected = set()
        for ring in macrocycles:
            affected.update(ring["atom_indices"])

        n = len(macrocycles)

        return CheckResult(
            check_name=self.name,
            passed=False,
            severity=Severity.INFO,
            message=f"Found {n} macrocyclic ring(s) with > 12 atoms.",
            affected_atoms=sorted(affected),
            details=details,
        )


@CheckRegistry.register("charged_species")
class ChargedSpeciesCheck(BaseCheck):
    """
    Detect charged atoms and identify zwitterions.

    Covers DVAL-16. Collects all formal charges on atoms, calculates net
    charge, and identifies zwitterions (net charge 0 but both positive and
    negative atoms present).

    Severity: INFO — charged species are common and valid drug structures.
    This is contextual information for the user.
    """

    name = "charged_species"
    description = "Detect charged atoms and zwitterionic structures"
    category = "structural_complexity"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Check for charged atoms and zwitterions.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with INFO severity if charged species found
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot check charged species of None molecule",
            )

        positive_atoms = []
        negative_atoms = []

        for atom in mol.GetAtoms():
            charge = atom.GetFormalCharge()
            if charge > 0:
                positive_atoms.append(
                    {
                        "atom_idx": atom.GetIdx(),
                        "symbol": atom.GetSymbol(),
                        "charge": charge,
                    }
                )
            elif charge < 0:
                negative_atoms.append(
                    {
                        "atom_idx": atom.GetIdx(),
                        "symbol": atom.GetSymbol(),
                        "charge": charge,
                    }
                )

        net_charge = sum(a["charge"] for a in positive_atoms) + sum(
            a["charge"] for a in negative_atoms
        )

        has_positive = len(positive_atoms) > 0
        has_negative = len(negative_atoms) > 0
        is_zwitterion = net_charge == 0 and has_positive and has_negative

        total_charged_atoms = len(positive_atoms) + len(negative_atoms)

        details = {
            "net_charge": net_charge,
            "positive_atoms": positive_atoms,
            "negative_atoms": negative_atoms,
            "is_zwitterion": is_zwitterion,
            "total_charged_atoms": total_charged_atoms,
        }

        if total_charged_atoms == 0:
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.INFO,
                message="No charged atoms detected. Net charge: 0.",
                details=details,
            )

        zwitterion_str = "Is a zwitterion." if is_zwitterion else "Is not a zwitterion."
        affected_atoms = (
            [a["atom_idx"] for a in positive_atoms]
            + [a["atom_idx"] for a in negative_atoms]
        )

        return CheckResult(
            check_name=self.name,
            passed=False,
            severity=Severity.INFO,
            message=(
                f"Net charge {net_charge}. {zwitterion_str} "
                f"{total_charged_atoms} charged atom(s)."
            ),
            affected_atoms=affected_atoms,
            details=details,
        )


@CheckRegistry.register("explicit_hydrogen_audit")
class ExplicitHydrogenAuditCheck(BaseCheck):
    """
    Audit explicit hydrogen specification in the molecule.

    Covers DVAL-17. Reports atoms with explicit hydrogen counts set
    (GetNumExplicitHs() > 0) and detects whether the molecule has been
    processed with AddHs() (separate H atom objects bonded to heavy atoms).

    Per RDKit pitfall #5: use GetNumExplicitHs() which returns H count
    stored on the heavy atom, not by counting H atoms in mol.GetAtoms().
    However, we also count actual H atom objects for completeness.

    Severity: INFO — explicit H specification is informational context.
    """

    name = "explicit_hydrogen_audit"
    description = "Audit explicit hydrogen specification patterns"
    category = "structural_complexity"

    def run(self, mol: Chem.Mol) -> CheckResult:
        """
        Audit explicit hydrogen patterns.

        Args:
            mol: RDKit molecule object

        Returns:
            CheckResult with INFO severity reporting explicit H patterns
        """
        if mol is None:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message="Cannot audit explicit hydrogens of None molecule",
            )

        # Collect atoms with explicit H count stored on the heavy atom
        atoms_with_explicit_h = []
        total_explicit_h = 0

        for atom in mol.GetAtoms():
            # Skip actual H atom objects — we handle those separately
            if atom.GetAtomicNum() == 1:
                continue
            explicit_h = atom.GetNumExplicitHs()
            if explicit_h > 0:
                atoms_with_explicit_h.append(
                    {
                        "atom_idx": atom.GetIdx(),
                        "symbol": atom.GetSymbol(),
                        "explicit_h_count": explicit_h,
                    }
                )
                total_explicit_h += explicit_h

        # Detect actual H atom objects (from AddHs() processing)
        h_atom_objects = [
            atom for atom in mol.GetAtoms()
            if atom.GetAtomicNum() == 1
        ]
        has_h_atom_objects = len(h_atom_objects) > 0
        h_atom_object_count = len(h_atom_objects)

        details = {
            "atoms_with_explicit_h": atoms_with_explicit_h,
            "total_explicit_h": total_explicit_h,
            "has_h_atom_objects": has_h_atom_objects,
            "h_atom_object_count": h_atom_object_count,
        }

        n = len(atoms_with_explicit_h)

        if n == 0 and not has_h_atom_objects:
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.INFO,
                message="No explicit hydrogen specification detected",
                details=details,
            )

        parts = []
        if n > 0:
            parts.append(f"Found {n} atom(s) with explicit hydrogen specification.")
        if has_h_atom_objects:
            parts.append(
                f"Molecule contains {h_atom_object_count} hydrogen atom object(s) "
                "(likely processed with AddHs())."
            )

        return CheckResult(
            check_name=self.name,
            passed=False,
            severity=Severity.INFO,
            message=" ".join(parts),
            affected_atoms=[a["atom_idx"] for a in atoms_with_explicit_h],
            details=details,
        )
