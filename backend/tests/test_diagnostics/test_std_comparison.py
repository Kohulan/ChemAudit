"""
Unit tests for Cross-Pipeline Standardization Comparison Service (DIAG-04).

Tests cover:
- All-agree result for simple molecule (ethanol)
- Disagreement for salt molecule (RDKit removes salt, Minimal keeps it)
- Exactly 3 pipelines present in result
- ValueError raised for invalid SMILES
"""

import pytest

from app.services.diagnostics.std_comparison import compare_pipelines


class TestAllAgree:
    """Tests for molecules where all pipelines produce the same result."""

    def test_all_agree(self) -> None:
        """Simple 'CCO' molecule gives all_agree=True and disagreements=0."""
        result = compare_pipelines("CCO")
        assert result["all_agree"] is True
        assert result["disagreements"] == 0

    def test_all_agree_result_structure(self) -> None:
        """Result contains all required keys."""
        result = compare_pipelines("CCO")
        for key in [
            "pipelines",
            "disagreements",
            "structural_disagreements",
            "all_agree",
            "property_comparison",
        ]:
            assert key in result, f"Missing key: {key}"


class TestThreePipelines:
    """Tests verifying exactly 3 pipelines are present."""

    def test_three_pipelines_present(self) -> None:
        """Result has exactly 3 pipelines in the list."""
        result = compare_pipelines("CCO")
        assert len(result["pipelines"]) == 3

    def test_pipeline_names(self) -> None:
        """Each pipeline result has a 'name' field."""
        result = compare_pipelines("CCO")
        pipeline_names = [p["name"] for p in result["pipelines"]]
        assert len(pipeline_names) == 3
        assert all(isinstance(n, str) and len(n) > 0 for n in pipeline_names)

    def test_pipeline_properties(self) -> None:
        """Non-failed pipelines have expected properties."""
        result = compare_pipelines("CCO")
        for pipeline in result["pipelines"]:
            if "error" not in pipeline:
                assert "smiles" in pipeline
                assert "inchikey" in pipeline
                assert "mw" in pipeline


class TestCrossPipelineDisagree:
    """Tests for molecules where pipelines may disagree."""

    def test_cross_pipeline_disagree(self) -> None:
        """Salt molecule may cause disagreements between pipelines."""
        # Sodium acetate — RDKit MolStandardize removes salt/counterion;
        # Minimal pipeline keeps it
        result = compare_pipelines("CC(=O)[O-].[Na+]")
        # Should return valid result (some pipelines may disagree)
        assert isinstance(result, dict)
        assert "pipelines" in result
        assert len(result["pipelines"]) == 3
        assert "disagreements" in result
        assert isinstance(result["disagreements"], int)
        # For a salt, pipelines are likely to disagree
        assert result["disagreements"] >= 0  # At least 0, likely > 0

    def test_salt_disagree_disagreements_positive(self) -> None:
        """Salt molecule should produce at least some disagreements."""
        result = compare_pipelines("CC(=O)[O-].[Na+]")
        # RDKit LargestFragmentChooser removes Na+, Minimal keeps it
        # This should result in at least 1 disagreement
        assert result["disagreements"] > 0

    def test_property_comparison_has_all_properties(self) -> None:
        """Property comparison list covers all 6 expected properties."""
        result = compare_pipelines("CCO")
        prop_names = [pc["property"] for pc in result["property_comparison"]]
        for expected_prop in ["smiles", "inchikey", "mw", "formula", "charge", "stereo_count"]:
            assert expected_prop in prop_names, f"Missing property: {expected_prop}"

    def test_property_comparison_structure(self) -> None:
        """Each property comparison row has correct structure."""
        result = compare_pipelines("CCO")
        for pc in result["property_comparison"]:
            assert "property" in pc
            assert "values" in pc
            assert "agrees" in pc
            assert "structural" in pc
            assert isinstance(pc["values"], list)
            assert len(pc["values"]) == 3  # One per pipeline


class TestInvalidSmiles:
    """Tests for invalid SMILES raising ValueError."""

    def test_invalid_smiles_raises(self) -> None:
        """Invalid SMILES raises ValueError (not HTTPException)."""
        with pytest.raises(ValueError):
            compare_pipelines("invalid_xyz")

    def test_truly_invalid_smiles_raises(self) -> None:
        """A clearly invalid SMILES string raises ValueError."""
        with pytest.raises(ValueError):
            compare_pipelines("not_a_smiles_XYZ")


class TestMcsHighlights:
    """MCS-based structural diff highlights on each pipeline result."""

    def test_all_agree_highlights_empty(self) -> None:
        """When pipelines agree, no atoms/bonds are flagged as divergent."""
        result = compare_pipelines("CCO")
        assert result["all_agree"] is True
        for pipeline in result["pipelines"]:
            assert pipeline["highlight_atoms"] == []
            assert pipeline["highlight_bonds"] == []

    def test_disagree_yields_highlights(self) -> None:
        """Salt SMILES: RDKit strips the counterion; at least one pipeline
        has non-empty highlights (atoms the others retained)."""
        result = compare_pipelines("CC(=O)[O-].[Na+]")
        assert result["structural_disagreements"] >= 1
        any_highlighted = any(
            len(p["highlight_atoms"]) > 0 or len(p["highlight_bonds"]) > 0
            for p in result["pipelines"]
            if "error" not in p
        )
        assert any_highlighted, "Expected at least one pipeline to flag divergent atoms"

    def test_highlight_indices_in_range(self) -> None:
        """Every highlighted index is a valid atom/bond index of the pipeline's SMILES."""
        from rdkit import Chem

        result = compare_pipelines("CC(=O)[O-].[Na+]")
        for pipeline in result["pipelines"]:
            if "error" in pipeline:
                continue
            mol = Chem.MolFromSmiles(pipeline["smiles"])
            assert mol is not None
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()
            for idx in pipeline["highlight_atoms"]:
                assert 0 <= idx < num_atoms, (
                    f"Atom index {idx} out of range (0, {num_atoms}) for "
                    f"{pipeline['name']} SMILES {pipeline['smiles']!r}"
                )
            for idx in pipeline["highlight_bonds"]:
                assert 0 <= idx < num_bonds, (
                    f"Bond index {idx} out of range (0, {num_bonds}) for "
                    f"{pipeline['name']} SMILES {pipeline['smiles']!r}"
                )

    def test_highlight_fields_always_present(self) -> None:
        """Every pipeline result exposes highlight_atoms and highlight_bonds."""
        result = compare_pipelines("CCO")
        for pipeline in result["pipelines"]:
            assert "highlight_atoms" in pipeline
            assert "highlight_bonds" in pipeline
            assert isinstance(pipeline["highlight_atoms"], list)
            assert isinstance(pipeline["highlight_bonds"], list)
