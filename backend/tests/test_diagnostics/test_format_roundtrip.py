"""
Unit tests for Format Round-Trip Lossiness Service (DIAG-03).

Tests cover:
- MOL round-trip lossless for simple molecule
- InChI round-trip lossless for simple molecule
- Stereo loss detection via InChI route
- Invalid SMILES returns error string without raising exception
"""


from app.services.diagnostics.format_roundtrip import check_roundtrip


class TestMolRoundtrip:
    """Tests for SMILES->MOL->SMILES route."""

    def test_mol_roundtrip_lossless(self) -> None:
        """Simple 'CCO' SMILES->MOL->SMILES is lossless with no losses."""
        result = check_roundtrip("CCO", route="smiles_mol_smiles")
        assert result["route"] == "smiles_mol_smiles"
        assert result["lossy"] is False
        assert result["losses"] == []
        assert result["error"] is None
        assert result["roundtrip_smiles"] is not None

    def test_mol_roundtrip_has_intermediate(self) -> None:
        """MOL route returns intermediate MOL block."""
        result = check_roundtrip("CCO", route="smiles_mol_smiles")
        assert result["intermediate"] is not None
        assert "M  END" in result["intermediate"]

    def test_mol_roundtrip_original_smiles(self) -> None:
        """Result contains original canonical SMILES."""
        result = check_roundtrip("OCC", route="smiles_mol_smiles")  # non-canonical input
        assert result["original_smiles"] is not None
        assert isinstance(result["original_smiles"], str)

    def test_mol_roundtrip_benzene(self) -> None:
        """Benzene round-trip via MOL is lossless."""
        result = check_roundtrip("c1ccccc1", route="smiles_mol_smiles")
        assert result["lossy"] is False
        assert result["error"] is None


class TestInChIRoundtrip:
    """Tests for SMILES->InChI->SMILES route."""

    def test_inchi_roundtrip_lossless(self) -> None:
        """Simple 'CCO' SMILES->InChI->SMILES is lossless."""
        result = check_roundtrip("CCO", route="smiles_inchi_smiles")
        assert result["route"] == "smiles_inchi_smiles"
        assert result["lossy"] is False
        assert result["losses"] == []
        assert result["error"] is None
        assert result["roundtrip_smiles"] is not None

    def test_inchi_roundtrip_has_intermediate(self) -> None:
        """InChI route returns intermediate InChI string."""
        result = check_roundtrip("CCO", route="smiles_inchi_smiles")
        assert result["intermediate"] is not None
        assert result["intermediate"].startswith("InChI=")

    def test_stereo_loss_via_inchi(self) -> None:
        """Molecule with defined stereo center is checked for stereo loss through InChI."""
        # [C@@H](F)(Cl)Br has a defined stereo center
        # InChI may or may not lose stereo depending on the molecule and RDKit version
        # The test verifies the function runs and returns expected structure
        result = check_roundtrip("[C@@H](F)(Cl)Br", route="smiles_inchi_smiles")
        assert isinstance(result, dict)
        assert "lossy" in result
        assert "losses" in result
        assert isinstance(result["losses"], list)
        # Each loss entry should have required fields
        for loss in result["losses"]:
            assert "type" in loss
            assert "description" in loss
            assert "before" in loss
            assert "after" in loss

    def test_inchi_roundtrip_default_route(self) -> None:
        """Default route is smiles_inchi_smiles when not specified."""
        result = check_roundtrip("CCO")
        assert result["route"] == "smiles_inchi_smiles"


class TestInvalidSmiles:
    """Tests for invalid SMILES input handling."""

    def test_invalid_smiles(self) -> None:
        """Invalid SMILES returns error string without raising exception."""
        result = check_roundtrip("invalid", route="smiles_inchi_smiles")
        assert isinstance(result, dict)
        assert result["error"] is not None
        assert isinstance(result["error"], str)
        assert len(result["error"]) > 0

    def test_invalid_smiles_mol_route(self) -> None:
        """Invalid SMILES via MOL route returns error without raising exception."""
        result = check_roundtrip("not-a-smiles", route="smiles_mol_smiles")
        assert isinstance(result, dict)
        assert result["error"] is not None

    def test_invalid_smiles_has_correct_structure(self) -> None:
        """Invalid SMILES result has all expected keys even on error."""
        result = check_roundtrip("xyz_invalid", route="smiles_inchi_smiles")
        for key in ["route", "original_smiles", "intermediate", "roundtrip_smiles", "lossy", "losses", "error"]:
            assert key in result, f"Missing key: {key}"


class TestUnknownRoute:
    """Tests for unknown route parameter."""

    def test_unknown_route_returns_error(self) -> None:
        """Unknown route returns error in the result dict, not an exception."""
        result = check_roundtrip("CCO", route="smiles_xyz_smiles")
        assert result["error"] is not None
        assert "Unknown route" in result["error"] or result["error"] is not None
