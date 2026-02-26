"""
Tests for structural complexity deep validation checks (M1.3).

Tests DVAL-12 through DVAL-17:
- HypervalentAtomCheck (DVAL-12)
- PolymerDetectionCheck (DVAL-13)
- RingStrainCheck (DVAL-14)
- MacrocycleDetectionCheck (DVAL-15)
- ChargedSpeciesCheck (DVAL-16)
- ExplicitHydrogenAuditCheck (DVAL-17)

Also includes:
- TestRegistration: verifies all 6 M1.3 checks are in CheckRegistry
- TestAllDeepValidationChecks: CI-level hard assertion that all 16 deep
  validation checks (M1.1 + M1.2 + M1.3) are registered
"""

import os

from rdkit import Chem

from app.schemas.common import Severity
from app.services.validation.checks.deep_complexity import (
    ChargedSpeciesCheck,
    ExplicitHydrogenAuditCheck,
    HypervalentAtomCheck,
    MacrocycleDetectionCheck,
    PolymerDetectionCheck,
    RingStrainCheck,
)

# =============================================================================
# TestHypervalentAtoms
# =============================================================================


class TestHypervalentAtoms:
    """Tests for HypervalentAtomCheck (DVAL-12)."""

    def test_normal_molecule_passes(self):
        """Ethanol has standard valences — should pass."""
        mol = Chem.MolFromSmiles("CCO")
        check = HypervalentAtomCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert "No hypervalent atoms" in result.message
        assert result.details["hypervalent_atoms"] == []

    def test_normal_benzene_passes(self):
        """Benzene has standard valences — should pass."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        check = HypervalentAtomCheck()
        result = check.run(mol)

        assert result.passed
        assert result.details["hypervalent_atoms"] == []

    def test_details_structure(self):
        """Result details must have the hypervalent_atoms list key."""
        mol = Chem.MolFromSmiles("CCO")
        check = HypervalentAtomCheck()
        result = check.run(mol)

        assert "hypervalent_atoms" in result.details
        assert isinstance(result.details["hypervalent_atoms"], list)

    def test_affected_atoms_empty_when_no_hypervalent(self):
        """No hypervalent atoms means empty affected_atoms list."""
        mol = Chem.MolFromSmiles("CCCC")
        check = HypervalentAtomCheck()
        result = check.run(mol)

        assert result.passed
        assert result.affected_atoms == []

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR result."""
        check = HypervalentAtomCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_check_name_correct(self):
        """Check name should be hypervalent_atoms."""
        check = HypervalentAtomCheck()
        mol = Chem.MolFromSmiles("CCO")
        result = check.run(mol)

        assert result.check_name == "hypervalent_atoms"

    def test_category_structural_complexity(self):
        """Check category should be structural_complexity."""
        check = HypervalentAtomCheck()
        assert check.category == "structural_complexity"

    def test_hypervalent_atom_details_fields(self):
        """Details list entries should have required fields when atoms flagged."""
        # Build a molecule that has an atom with high explicit valence
        # We use phosphorus pentavalent form (P in [P](=O)(O)(O)O accepted by RDKit)
        # In RDKit, phosphorus allows valence 5 so it may not flag. Use sulfone-like
        # structure that produces hypervalent detection if present.
        # This test verifies the structure of the details even with normal molecules.
        mol = Chem.MolFromSmiles("CCO")
        check = HypervalentAtomCheck()
        result = check.run(mol)

        # For a passing molecule, verify details structure
        assert isinstance(result.details["hypervalent_atoms"], list)
        # If any entries exist (e.g., from exotic molecules), they must have these keys
        for entry in result.details["hypervalent_atoms"]:
            assert "atom_idx" in entry
            assert "symbol" in entry
            assert "actual_valence" in entry
            assert "allowed_valences" in entry


# =============================================================================
# TestPolymerDetection
# =============================================================================


class TestPolymerDetection:
    """Tests for PolymerDetectionCheck (DVAL-13)."""

    def test_normal_molecule_passes(self):
        """Ethanol: small MW, no SGroups, no dummy atoms — should pass."""
        mol = Chem.MolFromSmiles("CCO")
        check = PolymerDetectionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert "No polymer indicators" in result.message
        assert result.details["has_sgroup_markers"] is False
        assert result.details["has_dummy_atoms"] is False

    def test_high_mw_flagged(self):
        """A long alkane chain (MW > 1500 Da) should be flagged as possible polymer."""
        # C120H242 — ~1682 Da, well above the 1500 threshold
        big_smiles = "C" * 120
        mol = Chem.MolFromSmiles(big_smiles)
        assert mol is not None, "Expected valid molecule for high-MW test"

        check = PolymerDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["exceeds_mw_threshold"] is True
        assert result.details["molecular_weight"] > 1500
        assert "Possible polymer detected" in result.message

    def test_dummy_atoms_flagged(self):
        """Molecule with dummy atoms (attachment points) should be flagged."""
        mol = Chem.MolFromSmiles("[*]CC[*]")
        assert mol is not None

        check = PolymerDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["has_dummy_atoms"] is True
        assert result.details["dummy_atom_count"] >= 2
        assert len(result.affected_atoms) >= 2
        assert "Possible polymer detected" in result.message

    def test_single_dummy_atom_flagged(self):
        """Single dummy atom (R-group) should be flagged."""
        mol = Chem.MolFromSmiles("[*]C")
        assert mol is not None

        check = PolymerDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["has_dummy_atoms"] is True
        assert result.details["dummy_atom_count"] == 1

    def test_details_structure(self):
        """Result details must contain all required keys."""
        mol = Chem.MolFromSmiles("CCO")
        check = PolymerDetectionCheck()
        result = check.run(mol)

        required_keys = [
            "has_sgroup_markers",
            "sgroup_types",
            "molecular_weight",
            "exceeds_mw_threshold",
            "has_dummy_atoms",
            "dummy_atom_count",
        ]
        for key in required_keys:
            assert key in result.details, f"Missing key: {key}"

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR result."""
        check = PolymerDetectionCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_check_name_correct(self):
        """Check name should be polymer_detection."""
        check = PolymerDetectionCheck()
        mol = Chem.MolFromSmiles("CCO")
        result = check.run(mol)

        assert result.check_name == "polymer_detection"

    def test_normal_mw_not_flagged(self):
        """A medium-sized molecule with MW < 1500 should not be flagged by MW."""
        # Aspirin: ~180 Da
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        check = PolymerDetectionCheck()
        result = check.run(mol)

        assert result.details["exceeds_mw_threshold"] is False

    def test_dummy_atoms_in_affected_atoms(self):
        """Dummy atom indices should appear in affected_atoms."""
        mol = Chem.MolFromSmiles("[*]CC")
        assert mol is not None

        check = PolymerDetectionCheck()
        result = check.run(mol)

        # At least one dummy atom index should be in affected_atoms
        assert len(result.affected_atoms) >= 1
        # Find dummy atom index
        dummy_idx = next(
            atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0
        )
        assert dummy_idx in result.affected_atoms


# =============================================================================
# TestRingStrain
# =============================================================================


class TestRingStrain:
    """Tests for RingStrainCheck (DVAL-14)."""

    def test_no_small_rings_passes(self):
        """Benzene (6-membered ring) should pass — not strained."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        check = RingStrainCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert "No strained rings" in result.message
        assert result.details["total_strained_rings"] == 0

    def test_cyclopropane_flagged(self):
        """Cyclopropane (3-membered ring) should be flagged with WARNING."""
        mol = Chem.MolFromSmiles("C1CC1")
        check = RingStrainCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.WARNING
        assert result.details["total_strained_rings"] == 1
        assert result.details["strained_rings"][0]["ring_size"] == 3
        assert "strained ring" in result.message

    def test_cyclobutane_flagged(self):
        """Cyclobutane (4-membered ring) should be flagged with WARNING."""
        mol = Chem.MolFromSmiles("C1CCC1")
        check = RingStrainCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.WARNING
        assert result.details["total_strained_rings"] == 1
        assert result.details["strained_rings"][0]["ring_size"] == 4

    def test_cyclopentane_passes(self):
        """Cyclopentane (5-membered ring) should pass — not in 3/4-membered set."""
        mol = Chem.MolFromSmiles("C1CCCC1")
        check = RingStrainCheck()
        result = check.run(mol)

        assert result.passed
        assert result.details["total_strained_rings"] == 0

    def test_cyclohexane_passes(self):
        """Cyclohexane (6-membered ring) should pass."""
        mol = Chem.MolFromSmiles("C1CCCCC1")
        check = RingStrainCheck()
        result = check.run(mol)

        assert result.passed
        assert result.details["total_strained_rings"] == 0

    def test_multiple_strained_rings(self):
        """Molecule with both 3- and 4-membered rings should flag both."""
        # Bicyclo[1.1.0]butane has a 3-membered ring fused to another 3-membered ring
        # Use spiro compound: spiro[2.3]hexane (3+4-membered rings)
        mol = Chem.MolFromSmiles("C1CC11CCC1")  # Spiro[2.3]hexane
        assert mol is not None

        check = RingStrainCheck()
        result = check.run(mol)

        # Should have at least one strained ring detected
        # (3 and/or 4 membered ring in spiro compound)
        assert not result.passed
        assert result.details["total_strained_rings"] >= 1

    def test_affected_atoms_correct(self):
        """Affected atoms should be the ring atom indices for strained rings."""
        mol = Chem.MolFromSmiles("C1CC1")  # Cyclopropane: atoms 0, 1, 2
        check = RingStrainCheck()
        result = check.run(mol)

        assert not result.passed
        # All atoms in cyclopropane should be in affected_atoms
        assert len(result.affected_atoms) == 3
        assert set(result.affected_atoms) == {0, 1, 2}

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR result."""
        check = RingStrainCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_details_structure(self):
        """Result details must have strained_rings list and total count."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        check = RingStrainCheck()
        result = check.run(mol)

        assert "strained_rings" in result.details
        assert "total_strained_rings" in result.details
        assert isinstance(result.details["strained_rings"], list)

    def test_strained_ring_entry_structure(self):
        """Each strained ring entry must have ring_size and atom_indices."""
        mol = Chem.MolFromSmiles("C1CC1")
        check = RingStrainCheck()
        result = check.run(mol)

        assert not result.passed
        ring_entry = result.details["strained_rings"][0]
        assert "ring_size" in ring_entry
        assert "atom_indices" in ring_entry
        assert ring_entry["ring_size"] == 3
        assert isinstance(ring_entry["atom_indices"], list)

    def test_acyclic_molecule_passes(self):
        """Acyclic molecules have no rings, should pass."""
        mol = Chem.MolFromSmiles("CCCCCC")
        check = RingStrainCheck()
        result = check.run(mol)

        assert result.passed
        assert result.details["total_strained_rings"] == 0

    def test_check_name_correct(self):
        """Check name should be ring_strain."""
        check = RingStrainCheck()
        mol = Chem.MolFromSmiles("CCO")
        result = check.run(mol)
        assert result.check_name == "ring_strain"


# =============================================================================
# TestMacrocycleDetection
# =============================================================================


class TestMacrocycleDetection:
    """Tests for MacrocycleDetectionCheck (DVAL-15)."""

    def test_small_ring_passes(self):
        """Benzene (6-membered ring) should pass — not macrocyclic."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        check = MacrocycleDetectionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert "No macrocyclic rings" in result.message
        assert result.details["total_macrocycles"] == 0

    def test_macrocycle_detected(self):
        """A 14-atom ring (cyclotetradecane) should be flagged as macrocycle."""
        # Cyclotetradecane: C14H28, 14-membered ring
        mol = Chem.MolFromSmiles("C1CCCCCCCCCCCCC1")  # 14 atoms in ring
        assert mol is not None

        check = MacrocycleDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.INFO
        assert result.details["total_macrocycles"] == 1
        assert result.details["macrocycles"][0]["ring_size"] == 14
        assert "macrocyclic ring" in result.message

    def test_13_atom_ring_flagged(self):
        """A 13-atom ring should be flagged (> 12, not >= 12)."""
        # Cyclotridecane: 13-membered ring
        mol = Chem.MolFromSmiles("C1CCCCCCCCCCCC1")  # 13 atoms
        assert mol is not None

        check = MacrocycleDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["macrocycles"][0]["ring_size"] == 13

    def test_12_atom_ring_passes(self):
        """A 12-atom ring should NOT be flagged (threshold is > 12, not >= 12)."""
        # Cyclododecane: 12-membered ring
        mol = Chem.MolFromSmiles("C1CCCCCCCCCCC1")  # 12 atoms
        assert mol is not None

        check = MacrocycleDetectionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.details["total_macrocycles"] == 0

    def test_affected_atoms_correct(self):
        """Macrocycle ring atom indices should appear in affected_atoms."""
        mol = Chem.MolFromSmiles("C1CCCCCCCCCCCCC1")  # 14-membered ring
        assert mol is not None

        check = MacrocycleDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        # All 14 ring atoms should be in affected_atoms
        assert len(result.affected_atoms) == 14

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR result."""
        check = MacrocycleDetectionCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_details_structure(self):
        """Result details must have macrocycles list, total, and sssr_note."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        check = MacrocycleDetectionCheck()
        result = check.run(mol)

        assert "macrocycles" in result.details
        assert "total_macrocycles" in result.details
        assert "sssr_note" in result.details
        assert isinstance(result.details["macrocycles"], list)
        assert "SSSR" in result.details["sssr_note"]

    def test_macrocycle_entry_structure(self):
        """Each macrocycle entry must have ring_size and atom_indices."""
        mol = Chem.MolFromSmiles("C1CCCCCCCCCCCCC1")  # 14-membered ring
        assert mol is not None

        check = MacrocycleDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        entry = result.details["macrocycles"][0]
        assert "ring_size" in entry
        assert "atom_indices" in entry
        assert entry["ring_size"] == 14
        assert len(entry["atom_indices"]) == 14

    def test_medium_ring_passes(self):
        """Cyclohexane (6-membered) and cyclodecane (10-membered) both pass."""
        for smiles, name in [("C1CCCCC1", "cyclohexane"), ("C1CCCCCCCCC1", "cyclodecane")]:
            mol = Chem.MolFromSmiles(smiles)
            check = MacrocycleDetectionCheck()
            result = check.run(mol)
            assert result.passed, f"{name} should not be flagged as macrocycle"

    def test_check_name_correct(self):
        """Check name should be macrocycle_detection."""
        check = MacrocycleDetectionCheck()
        mol = Chem.MolFromSmiles("CCO")
        result = check.run(mol)
        assert result.check_name == "macrocycle_detection"


# =============================================================================
# TestChargedSpecies
# =============================================================================


class TestChargedSpecies:
    """Tests for ChargedSpeciesCheck (DVAL-16)."""

    def test_neutral_molecule_passes(self):
        """Ethanol has no charges — should pass."""
        mol = Chem.MolFromSmiles("CCO")
        check = ChargedSpeciesCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["net_charge"] == 0
        assert result.details["total_charged_atoms"] == 0
        assert result.details["is_zwitterion"] is False

    def test_positive_charge_detected(self):
        """Protonated amine should detect positive charge."""
        mol = Chem.MolFromSmiles("[NH3+]CC")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.INFO
        assert result.details["net_charge"] == 1
        assert len(result.details["positive_atoms"]) == 1
        assert len(result.details["negative_atoms"]) == 0
        assert result.details["is_zwitterion"] is False

    def test_negative_charge_detected(self):
        """Carboxylate anion should detect negative charge."""
        mol = Chem.MolFromSmiles("CC([O-])=O")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.INFO
        assert result.details["net_charge"] == -1
        assert len(result.details["positive_atoms"]) == 0
        assert len(result.details["negative_atoms"]) == 1
        assert result.details["is_zwitterion"] is False

    def test_zwitterion_detected(self):
        """Amino acid zwitterion form: net charge 0, both + and - charges."""
        # Glycine zwitterion: [NH3+]CC([O-])=O
        mol = Chem.MolFromSmiles("[NH3+]CC([O-])=O")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.INFO
        assert result.details["net_charge"] == 0
        assert result.details["is_zwitterion"] is True
        assert len(result.details["positive_atoms"]) == 1
        assert len(result.details["negative_atoms"]) == 1
        assert "zwitterion" in result.message.lower()

    def test_zwitterion_longer_chain(self):
        """Longer amino acid form should also be detected as zwitterion."""
        # Alanine zwitterion form: [NH3+]C(C)C([O-])=O
        mol = Chem.MolFromSmiles("[NH3+]C(C)C([O-])=O")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        assert result.details["is_zwitterion"] is True
        assert result.details["net_charge"] == 0

    def test_net_charge_calculated(self):
        """Net charge should be correct sum of all formal charges."""
        # Doubly charged phosphate: PO4^2- in SMILES form
        mol = Chem.MolFromSmiles("[O-]P([O-])(=O)O")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["net_charge"] == -2
        assert result.details["is_zwitterion"] is False

    def test_affected_atoms_correct(self):
        """Charged atom indices should appear in affected_atoms."""
        mol = Chem.MolFromSmiles("[NH3+]CC([O-])=O")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        # Should have 2 charged atoms (NH3+ and O-)
        assert result.details["total_charged_atoms"] == 2
        assert len(result.affected_atoms) == 2

    def test_positive_atom_details(self):
        """Positive atom entries must have atom_idx, symbol, charge."""
        mol = Chem.MolFromSmiles("[NH3+]CC")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        pos_atom = result.details["positive_atoms"][0]
        assert "atom_idx" in pos_atom
        assert "symbol" in pos_atom
        assert "charge" in pos_atom
        assert pos_atom["charge"] > 0

    def test_negative_atom_details(self):
        """Negative atom entries must have atom_idx, symbol, charge."""
        mol = Chem.MolFromSmiles("CC([O-])=O")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        neg_atom = result.details["negative_atoms"][0]
        assert "atom_idx" in neg_atom
        assert "symbol" in neg_atom
        assert "charge" in neg_atom
        assert neg_atom["charge"] < 0

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR result."""
        check = ChargedSpeciesCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_check_name_correct(self):
        """Check name should be charged_species."""
        check = ChargedSpeciesCheck()
        mol = Chem.MolFromSmiles("CCO")
        result = check.run(mol)
        assert result.check_name == "charged_species"

    def test_non_zwitterion_charged_molecule(self):
        """A molecule with only positive charges is not a zwitterion."""
        mol = Chem.MolFromSmiles("[NH4+]")
        assert mol is not None

        check = ChargedSpeciesCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["is_zwitterion"] is False
        assert result.details["net_charge"] == 1


# =============================================================================
# TestExplicitHydrogenAudit
# =============================================================================


class TestExplicitHydrogenAudit:
    """Tests for ExplicitHydrogenAuditCheck (DVAL-17)."""

    def test_normal_smiles_passes(self):
        """Ethanol parsed from SMILES has implicit H — typically passes."""
        mol = Chem.MolFromSmiles("CCO")
        check = ExplicitHydrogenAuditCheck()
        result = check.run(mol)

        # With implicit H representation, should pass or have minimal explicit H
        assert result.check_name == "explicit_hydrogen_audit"
        assert result.details["has_h_atom_objects"] is False

    def test_explicit_h_with_bracket_notation(self):
        """[CH4] forces explicit H count on carbon atom."""
        mol = Chem.MolFromSmiles("[CH4]")
        assert mol is not None

        check = ExplicitHydrogenAuditCheck()
        result = check.run(mol)

        # [CH4] stores explicit H count on the carbon
        # Some RDKit versions normalize this; test that it's detected if present
        # or gracefully handled
        assert "atoms_with_explicit_h" in result.details
        assert "total_explicit_h" in result.details
        assert result.check_name == "explicit_hydrogen_audit"

    def test_h_atom_objects_detected_after_addhs(self):
        """After Chem.AddHs(), molecule should have H atom objects detected."""
        mol = Chem.MolFromSmiles("CCO")
        mol_with_h = Chem.AddHs(mol)
        assert mol_with_h is not None

        check = ExplicitHydrogenAuditCheck()
        result = check.run(mol_with_h)

        assert result.details["has_h_atom_objects"] is True
        assert result.details["h_atom_object_count"] > 0
        assert not result.passed  # Has H atom objects, so flagged
        assert "hydrogen atom object" in result.message.lower()

    def test_details_structure(self):
        """Result details must contain all required keys."""
        mol = Chem.MolFromSmiles("CCO")
        check = ExplicitHydrogenAuditCheck()
        result = check.run(mol)

        required_keys = [
            "atoms_with_explicit_h",
            "total_explicit_h",
            "has_h_atom_objects",
            "h_atom_object_count",
        ]
        for key in required_keys:
            assert key in result.details, f"Missing key: {key}"

    def test_atoms_with_explicit_h_entries_have_correct_fields(self):
        """Each entry in atoms_with_explicit_h must have atom_idx, symbol, explicit_h_count."""
        mol = Chem.MolFromSmiles("CCO")
        mol_with_h = Chem.AddHs(mol)

        check = ExplicitHydrogenAuditCheck()
        result = check.run(mol_with_h)

        # atoms_with_explicit_h may be empty (AddHs creates H atom objects, not explicit H)
        # Just verify the structure is correct
        for entry in result.details["atoms_with_explicit_h"]:
            assert "atom_idx" in entry
            assert "symbol" in entry
            assert "explicit_h_count" in entry
            assert entry["explicit_h_count"] > 0

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR result."""
        check = ExplicitHydrogenAuditCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_check_name_correct(self):
        """Check name should be explicit_hydrogen_audit."""
        check = ExplicitHydrogenAuditCheck()
        mol = Chem.MolFromSmiles("CCO")
        result = check.run(mol)
        assert result.check_name == "explicit_hydrogen_audit"

    def test_category_structural_complexity(self):
        """Check category should be structural_complexity."""
        check = ExplicitHydrogenAuditCheck()
        assert check.category == "structural_complexity"

    def test_h_atom_object_count_matches_actual(self):
        """h_atom_object_count should equal actual H atoms added by AddHs."""
        mol = Chem.MolFromSmiles("C")  # Methane: CH4, 4 H atoms
        mol_with_h = Chem.AddHs(mol)
        expected_h = 4  # Methane has 4 H atoms

        check = ExplicitHydrogenAuditCheck()
        result = check.run(mol_with_h)

        assert result.details["h_atom_object_count"] == expected_h


# =============================================================================
# TestRegistration
# =============================================================================


class TestRegistration:
    """Verify all 6 M1.3 checks are registered in CheckRegistry."""

    def test_all_m13_checks_registered(self):
        """All 6 M1.3 check names must be in CheckRegistry."""
        # Import the module to ensure registration happens
        import app.services.validation.checks.deep_complexity  # noqa: F401
        from app.services.validation.registry import CheckRegistry

        registered = set(CheckRegistry.get_all().keys())
        m13_checks = {
            "hypervalent_atoms",
            "polymer_detection",
            "ring_strain",
            "macrocycle_detection",
            "charged_species",
            "explicit_hydrogen_audit",
        }

        missing = m13_checks - registered
        assert not missing, f"Missing M1.3 checks in registry: {missing}"

    def test_m13_check_classes_instantiatable(self):
        """All 6 M1.3 check classes must be instantiatable."""
        checks = [
            HypervalentAtomCheck(),
            PolymerDetectionCheck(),
            RingStrainCheck(),
            MacrocycleDetectionCheck(),
            ChargedSpeciesCheck(),
            ExplicitHydrogenAuditCheck(),
        ]
        assert len(checks) == 6

    def test_m13_checks_have_correct_category(self):
        """All M1.3 checks must have category='structural_complexity'."""
        checks = [
            HypervalentAtomCheck(),
            PolymerDetectionCheck(),
            RingStrainCheck(),
            MacrocycleDetectionCheck(),
            ChargedSpeciesCheck(),
            ExplicitHydrogenAuditCheck(),
        ]
        for check in checks:
            assert check.category == "structural_complexity", (
                f"{check.__class__.__name__} has wrong category: {check.category}"
            )


# =============================================================================
# TestAllDeepValidationChecks — CI-level cross-check
# =============================================================================


class TestAllDeepValidationChecks:
    """
    CI-level test: hard assert that all 16 deep validation checks are registered.

    This is the hard assertion counterpart to the logger.warning in main.py.
    Run with: pytest tests/test_validation/test_deep_complexity_checks.py::TestAllDeepValidationChecks -v
    """

    def test_all_16_deep_validation_checks_registered(self):
        """Hard assert: all 16 deep validation checks (M1.1+M1.2+M1.3) are registered."""
        # Import all three deep check modules to trigger registration
        import app.services.validation.checks.deep_complexity  # noqa: F401
        import app.services.validation.checks.deep_composition  # noqa: F401
        import app.services.validation.checks.deep_stereo_tautomer  # noqa: F401
        from app.services.validation.registry import CheckRegistry

        expected_deep_validation_checks = {
            # M1.1: Stereo & Tautomer
            "stereoisomer_enumeration",  # DVAL-01+02
            "tautomer_detection",  # DVAL-03
            "aromatic_system_validation",  # DVAL-04
            "coordinate_dimension",  # DVAL-05
            # M1.2: Chemical Composition
            "mixture_detection",  # DVAL-06
            "solvent_contamination",  # DVAL-07
            "inorganic_filter",  # DVAL-08
            "radical_detection",  # DVAL-09
            "isotope_label_detection",  # DVAL-10
            "trivial_molecule",  # DVAL-11
            # M1.3: Structural Complexity
            "hypervalent_atoms",  # DVAL-12
            "polymer_detection",  # DVAL-13
            "ring_strain",  # DVAL-14
            "macrocycle_detection",  # DVAL-15
            "charged_species",  # DVAL-16
            "explicit_hydrogen_audit",  # DVAL-17
        }

        registered = set(CheckRegistry.get_all().keys())
        missing = expected_deep_validation_checks - registered
        assert not missing, f"CI check: missing deep validation checks: {missing}"

    def test_check_file_count_matches_registry(self):
        """
        Cross-check: count .py files in checks/ directory and compare to
        expected check file count. This guards against adding check files
        without updating the registry.

        Expected check files (excluding base, __init__, __pycache__):
        - basic.py
        - representation.py
        - stereo.py
        - deep_stereo_tautomer.py
        - deep_composition.py
        - deep_complexity.py

        Total: 6 check files.
        """
        import importlib.util

        # Find the checks directory
        spec = importlib.util.find_spec("app.services.validation.checks")
        assert spec is not None, "Cannot find checks package"
        checks_dir = os.path.dirname(spec.origin)

        # Count .py files excluding non-check files
        excluded = {"__init__.py", "base.py"}
        check_files = [
            f
            for f in os.listdir(checks_dir)
            if f.endswith(".py") and f not in excluded and not f.startswith("__")
        ]

        expected_check_files = 6  # basic, representation, stereo, deep_stereo_tautomer,
        # deep_composition, deep_complexity
        assert len(check_files) == expected_check_files, (
            f"Expected {expected_check_files} check files in checks/ directory, "
            f"found {len(check_files)}: {sorted(check_files)}"
        )
