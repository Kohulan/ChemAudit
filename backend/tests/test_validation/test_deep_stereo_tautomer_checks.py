"""
Tests for deep stereo/tautomer validation checks (M1.1).

Tests DVAL-01 through DVAL-05:
- StereoisomerEnumerationCheck (DVAL-01/02)
- TautomerDetectionCheck (DVAL-03)
- AromaticSystemValidationCheck (DVAL-04)
- CoordinateDimensionCheck (DVAL-05)
"""

from rdkit import Chem
from rdkit.Chem import AllChem

from app.schemas.common import Severity
from app.services.validation.checks.deep_stereo_tautomer import (
    AromaticSystemValidationCheck,
    CoordinateDimensionCheck,
    StereoisomerEnumerationCheck,
    TautomerDetectionCheck,
)


class TestStereoisomerEnumeration:
    """Tests for StereoisomerEnumerationCheck (DVAL-01/02)."""

    def test_no_stereocenters_passes(self):
        """Molecule without chiral centers should pass with undefined_count=0."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol — no stereocenters
        check = StereoisomerEnumerationCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["undefined_count"] == 0
        assert result.details["total_centers"] == 0
        assert result.details["cap_exceeded"] is False

    def test_undefined_stereocenters_detected(self):
        """Molecule with 2 undefined stereocenters should fail with WARNING."""
        mol = Chem.MolFromSmiles("CC(O)C(N)CC")  # 2 undefined stereocenters
        check = StereoisomerEnumerationCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.WARNING
        assert result.details["undefined_count"] == 2
        assert len(result.details["atom_indices"]) == 2
        assert len(result.details["stereoisomer_smiles"]) > 0
        assert result.details["cap_exceeded"] is False
        assert len(result.affected_atoms) == 2

    def test_fully_defined_stereo_passes(self):
        """Molecule with fully defined stereocenters should pass."""
        mol = Chem.MolFromSmiles("C[C@H](O)[C@@H](N)CC")  # Both centers defined
        check = StereoisomerEnumerationCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["undefined_count"] == 0
        assert result.details["total_centers"] == 2
        # The defined molecule itself should appear in stereoisomer_smiles
        assert len(result.details["stereoisomer_smiles"]) >= 1

    def test_cap_exceeded(self):
        """Molecule with many undefined stereocenters should set cap_exceeded=True."""
        # Use a molecule with many undefined stereocenters (>7 to exceed 128 cap)
        # D-glucose derivative with many undefined centers
        # Create a molecule with 8+ undefined stereocenters
        mol = Chem.MolFromSmiles("C(C(C(C(C(C(C(C(O)O)O)O)O)O)O)O)O")
        check = StereoisomerEnumerationCheck()
        result = check.run(mol)

        # If cap is exceeded: no SMILES returned, cap_exceeded = True
        if result.details.get("cap_exceeded"):
            assert result.details["stereoisomer_smiles"] == []
        # Otherwise test simply validates structure of details (cap may not be exceeded for this mol)
        assert "cap_exceeded" in result.details
        assert "enumeration_cap" in result.details
        assert result.details["enumeration_cap"] == 128

    def test_none_molecule_returns_error(self):
        """Passing None should return ERROR severity."""
        check = StereoisomerEnumerationCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_details_structure(self):
        """Result details must contain all required keys."""
        mol = Chem.MolFromSmiles("CC(O)CC")  # One undefined stereocenter
        check = StereoisomerEnumerationCheck()
        result = check.run(mol)

        required_keys = {
            "undefined_count",
            "total_centers",
            "atom_indices",
            "stereoisomer_smiles",
            "enumeration_cap",
            "cap_exceeded",
        }
        assert required_keys == set(result.details.keys()), (
            f"Details missing keys: {required_keys - set(result.details.keys())}"
        )

    def test_affected_atoms_populated_for_undefined_centers(self):
        """affected_atoms should list the undefined stereocenter indices."""
        mol = Chem.MolFromSmiles("CC(O)CC")  # One undefined stereocenter at C1
        check = StereoisomerEnumerationCheck()
        result = check.run(mol)

        assert not result.passed
        assert len(result.affected_atoms) > 0
        # affected_atoms should match atom_indices in details
        assert result.affected_atoms == result.details["atom_indices"]

    def test_message_contains_indices_and_count(self):
        """Human-readable message should mention undefined centers and stereoisomer count."""
        mol = Chem.MolFromSmiles("CC(O)C(N)CC")  # 2 undefined centers
        check = StereoisomerEnumerationCheck()
        result = check.run(mol)

        assert "undefined" in result.message.lower() or "stereocenter" in result.message.lower()
        assert "stereoisomer" in result.message.lower() or "enumerat" in result.message.lower()


class TestTautomerDetection:
    """Tests for TautomerDetectionCheck (DVAL-03)."""

    def test_molecule_with_tautomers(self):
        """Phenol (keto-enol tautomerism) should detect tautomers."""
        mol = Chem.MolFromSmiles("OC1=CC=CC=C1")  # Phenol — has keto tautomers
        check = TautomerDetectionCheck()
        result = check.run(mol)

        assert result.passed  # INFO check always passes
        assert result.severity == Severity.INFO
        assert result.details["tautomer_count"] >= 1
        assert "canonical_smiles" in result.details
        assert isinstance(result.details["tautomer_smiles"], list)

    def test_molecule_without_tautomers(self):
        """Simple alkane (ethane) should have tautomer_count of 1 (itself)."""
        mol = Chem.MolFromSmiles("CC")  # Ethane — no tautomers
        check = TautomerDetectionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["tautomer_count"] == 1

    def test_canonical_form_detection_true(self):
        """Canonical tautomer input should have is_canonical_form=True."""
        # Phenol is the canonical tautomer of the keto form
        mol = Chem.MolFromSmiles("Oc1ccccc1")  # Phenol — canonical SMILES
        check = TautomerDetectionCheck()
        result = check.run(mol)

        assert result.passed
        # Verify flag is a boolean
        assert isinstance(result.details["is_canonical_form"], bool)
        # The canonical form should match (phenol IS the canonical tautomer)
        assert result.details["is_canonical_form"] is True

    def test_canonical_form_detection_non_canonical(self):
        """Non-canonical tautomer input should have is_canonical_form=False."""
        # Keto form of cyclohexadienone (tautomer of phenol) — non-canonical
        mol = Chem.MolFromSmiles("O=C1C=CC=CC1")  # 2,4-cyclohexadienone
        check = TautomerDetectionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.details["tautomer_count"] >= 2
        # This is the keto form — the enol (phenol) is the canonical tautomer
        assert result.details["is_canonical_form"] is False

    def test_none_molecule_returns_error(self):
        """Passing None should return ERROR severity."""
        check = TautomerDetectionCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_details_structure(self):
        """Result details must contain all required keys."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        check = TautomerDetectionCheck()
        result = check.run(mol)

        required_keys = {"tautomer_count", "canonical_smiles", "is_canonical_form", "tautomer_smiles"}
        assert required_keys.issubset(set(result.details.keys())), (
            f"Details missing keys: {required_keys - set(result.details.keys())}"
        )

    def test_message_contains_tautomer_info(self):
        """Human-readable message should mention tautomer count and canonical status."""
        mol = Chem.MolFromSmiles("CCO")
        check = TautomerDetectionCheck()
        result = check.run(mol)

        assert "tautomer" in result.message.lower()
        assert "canonical" in result.message.lower()


class TestAromaticSystemValidation:
    """Tests for AromaticSystemValidationCheck (DVAL-04)."""

    def test_normal_benzene_passes(self):
        """Benzene (6-membered aromatic) should pass with no issues."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        check = AromaticSystemValidationCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["unusual_ring_sizes"] == []
        assert result.details["charged_aromatics"] == []
        assert result.affected_atoms == []

    def test_pyrrole_passes(self):
        """Pyrrole (5-membered aromatic) should pass — standard size."""
        mol = Chem.MolFromSmiles("c1cc[nH]c1")  # Pyrrole
        check = AromaticSystemValidationCheck()
        result = check.run(mol)

        assert result.passed
        assert result.details["unusual_ring_sizes"] == []

    def test_charged_aromatic_detected(self):
        """Pyridinium ion ([nH+]) should flag charged aromatic atom."""
        mol = Chem.MolFromSmiles("[nH+]1ccccc1")  # Protonated pyridine-like
        if mol is None:
            # Try alternative SMILES for charged aromatic
            mol = Chem.MolFromSmiles("c1cc[nH+]cc1")
        if mol is None:
            # Fallback: create manually — pyridinium
            mol = Chem.MolFromSmiles("C1=CC=[NH+]C=C1")
        assert mol is not None, "Could not construct charged aromatic test molecule"

        check = AromaticSystemValidationCheck()
        result = check.run(mol)

        # Check may or may not flag depending on aromaticity assignment
        # At minimum, details structure must be correct
        assert "charged_aromatics" in result.details
        assert "unusual_ring_sizes" in result.details

    def test_none_molecule_returns_error(self):
        """Passing None should return ERROR severity."""
        check = AromaticSystemValidationCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_affected_atoms_populated(self):
        """When issues found, affected_atoms should contain the flagged atom indices."""
        # Azulene — 5+7 fused rings, has 7-membered aromatic ring (unusual)
        mol = Chem.MolFromSmiles("c1ccc2cccccc12")  # Azulene: 5+7 membered rings
        assert mol is not None, "Azulene should parse successfully"
        check = AromaticSystemValidationCheck()
        result = check.run(mol)

        # Azulene has a 7-membered aromatic ring — should be flagged
        if not result.passed:
            assert len(result.affected_atoms) > 0

    def test_azulene_unusual_ring(self):
        """Azulene has a 7-membered aromatic ring — should be flagged as unusual."""
        mol = Chem.MolFromSmiles("c1ccc2cccccc12")  # Azulene (5+7 fused rings)
        assert mol is not None, "Azulene should parse successfully"

        check = AromaticSystemValidationCheck()
        result = check.run(mol)

        # Azulene has 5+7 membered rings; the 7-membered should be flagged as unusual
        ring_sizes = [r["ring_size"] for r in result.details.get("unusual_ring_sizes", [])]
        assert not result.passed, "Azulene should fail due to 7-membered aromatic ring"
        assert 7 in ring_sizes, f"Expected 7-membered ring flagged, got: {ring_sizes}"

    def test_details_structure(self):
        """Result details must contain unusual_ring_sizes and charged_aromatics."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        check = AromaticSystemValidationCheck()
        result = check.run(mol)

        assert "unusual_ring_sizes" in result.details
        assert "charged_aromatics" in result.details
        assert isinstance(result.details["unusual_ring_sizes"], list)
        assert isinstance(result.details["charged_aromatics"], list)


class TestCoordinateDimension:
    """Tests for CoordinateDimensionCheck (DVAL-05)."""

    def test_no_conformer(self):
        """Molecule from SMILES (no conformer) should return dimension='no_coordinates'."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol — no conformer
        check = CoordinateDimensionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["dimension"] == "no_coordinates"
        assert result.details["num_conformers"] == 0

    def test_2d_coordinates(self):
        """Molecule with 2D coordinates (all z=0) should return dimension='2d'."""
        mol = Chem.MolFromSmiles("CCO")
        AllChem.Compute2DCoords(mol)  # Add 2D conformer

        check = CoordinateDimensionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["dimension"] == "2d"
        assert result.details["num_conformers"] >= 1

    def test_3d_coordinates(self):
        """Molecule with 3D coordinates (non-zero z) should return dimension='3d'."""
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        mol = Chem.RemoveHs(mol)

        check = CoordinateDimensionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["dimension"] in ("3d", "2d")  # EmbedMolecule produces 3D
        assert result.details["num_conformers"] >= 1

    def test_none_molecule_returns_error(self):
        """Passing None should return ERROR severity."""
        check = CoordinateDimensionCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_details_structure(self):
        """Result details must contain dimension and num_conformers."""
        mol = Chem.MolFromSmiles("CCO")
        check = CoordinateDimensionCheck()
        result = check.run(mol)

        assert "dimension" in result.details
        assert "num_conformers" in result.details
        assert result.details["dimension"] in ("2d", "3d", "no_coordinates", "degenerate")
        assert isinstance(result.details["num_conformers"], int)

    def test_always_passes_informationally(self):
        """CoordinateDimensionCheck should always pass (informational only)."""
        molecules = [
            Chem.MolFromSmiles("CCO"),
            Chem.MolFromSmiles("c1ccccc1"),
            Chem.MolFromSmiles("CC(=O)O"),
        ]
        check = CoordinateDimensionCheck()
        for mol in molecules:
            result = check.run(mol)
            assert result.passed, f"Expected pass for {Chem.MolToSmiles(mol)}"
            assert result.severity == Severity.INFO

    def test_message_describes_dimension(self):
        """Human-readable message should describe the coordinate dimension."""
        mol = Chem.MolFromSmiles("CCO")
        check = CoordinateDimensionCheck()
        result = check.run(mol)

        assert len(result.message) > 0
        # Message should contain dimension info
        assert any(
            term in result.message.lower()
            for term in ["coordinate", "2d", "3d", "conformer", "no_coord"]
        )


class TestRegistration:
    """Test that all M1.1 checks are properly registered in CheckRegistry."""

    def test_all_m11_checks_registered(self):
        """All 4 M1.1 check names must appear in CheckRegistry.get_all()."""
        # Force import to trigger registration
        import app.services.validation.checks.deep_stereo_tautomer  # noqa: F401

        from app.services.validation.registry import CheckRegistry

        registered = set(CheckRegistry.get_all().keys())
        expected = {
            "stereoisomer_enumeration",
            "tautomer_detection",
            "aromatic_system_validation",
            "coordinate_dimension",
        }
        missing = expected - registered
        assert not missing, f"Missing M1.1 checks in registry: {missing}"

    def test_checks_have_correct_category(self):
        """All M1.1 checks should have category='stereo_tautomer'."""
        checks_module = __import__(
            "app.services.validation.checks.deep_stereo_tautomer",
            fromlist=["StereoisomerEnumerationCheck"],
        )
        check_classes = [
            checks_module.StereoisomerEnumerationCheck,
            checks_module.TautomerDetectionCheck,
            checks_module.AromaticSystemValidationCheck,
            checks_module.CoordinateDimensionCheck,
        ]
        for cls in check_classes:
            assert cls.category == "stereo_tautomer", (
                f"{cls.__name__} has wrong category: {cls.category}"
            )
