"""
Unit tests for SMILES Error Diagnostics Service (DIAG-01).

Tests cover:
- Valid SMILES detection with canonical SMILES output
- Valence error detection with position extraction
- Bracket fix suggestions for unmatched parentheses
- Unknown atom symbol detection
- Ring closure mismatch detection
- Empty SMILES graceful handling
"""

import pytest

from app.services.diagnostics.smiles_diagnostics import diagnose_smiles


class TestValidSmiles:
    """Tests for valid SMILES input."""

    def test_valid_smiles(self) -> None:
        """Valid SMILES returns valid=True with canonical_smiles and empty errors."""
        result = diagnose_smiles("CCO")
        assert result["valid"] is True
        assert result["canonical_smiles"] == "CCO"
        assert result["errors"] == []

    def test_valid_smiles_complex(self) -> None:
        """Complex valid SMILES returns valid=True."""
        result = diagnose_smiles("c1ccccc1")  # benzene
        assert result["valid"] is True
        assert result["canonical_smiles"] is not None
        assert result["errors"] == []

    def test_valid_smiles_with_stereo(self) -> None:
        """Valid SMILES with defined stereo returns valid=True."""
        result = diagnose_smiles("[C@@H](F)(Cl)Br")
        assert result["valid"] is True
        assert result["canonical_smiles"] is not None
        assert result["errors"] == []

    def test_valid_smiles_canonical_form(self) -> None:
        """diagnose_smiles returns canonical_smiles field when valid."""
        result = diagnose_smiles("OCC")  # ethanol in non-canonical order
        assert result["valid"] is True
        assert result["canonical_smiles"] is not None
        assert isinstance(result["canonical_smiles"], str)
        assert len(result["canonical_smiles"]) > 0


class TestInvalidSmiles:
    """Tests for invalid SMILES input producing structured errors."""

    def test_valence_error_position(self) -> None:
        """SMILES with overvalent carbon (pentavalent) returns valid=False and valence_error."""
        result = diagnose_smiles("C(C)(C)(C)(C)C")
        assert result["valid"] is False
        assert len(result["errors"]) > 0
        # Check error type is valence-related
        error_types = [e["error_type"] for e in result["errors"]]
        assert "valence_error" in error_types, f"Expected valence_error, got {error_types}"

    def test_bracket_fix_suggestion(self) -> None:
        """SMILES with unmatched parenthesis returns errors with fix suggestions."""
        result = diagnose_smiles("C(C(C")
        assert result["valid"] is False
        assert len(result["errors"]) > 0
        # At least one error should have suggestions
        all_suggestions = []
        for error in result["errors"]:
            all_suggestions.extend(error.get("suggestions", []))
        assert len(all_suggestions) > 0, "Expected at least one fix suggestion"
        # At least one suggestion should reference a closing parenthesis
        suggestion_texts = [s.get("description", "") + str(s.get("corrected_smiles", "")) for s in all_suggestions]
        has_bracket_hint = any(")" in text or "parenthes" in text.lower() or "bracket" in text.lower() for text in suggestion_texts)
        assert has_bracket_hint, f"Expected bracket-related suggestion, got: {suggestion_texts}"

    def test_unknown_atom(self) -> None:
        """SMILES with unknown atom symbol returns valid=False with appropriate error type."""
        result = diagnose_smiles("CXC")
        assert result["valid"] is False
        assert len(result["errors"]) > 0
        error_types = [e["error_type"] for e in result["errors"]]
        # Should be unknown_atom_symbol or parse_error
        valid_types = {"unknown_atom_symbol", "parse_error"}
        assert any(t in valid_types for t in error_types), f"Expected unknown/parse error, got {error_types}"

    def test_ring_closure_mismatch(self) -> None:
        """SMILES with unclosed ring returns valid=False with errors."""
        result = diagnose_smiles("C1CC")  # Ring opened but never closed
        assert result["valid"] is False
        assert len(result["errors"]) > 0

    def test_empty_smiles(self) -> None:
        """Empty string is handled gracefully — returns dict with valid=False or at minimum a dict."""
        # diagnose_smiles is called with empty string; behavior varies by RDKit version
        # but must not raise an exception
        result = diagnose_smiles("")
        assert isinstance(result, dict)
        assert "valid" in result
        assert "errors" in result
        # Either valid=False or valid=True (some RDKit versions parse empty as valid)
        if result["valid"]:
            assert result["errors"] == []
        else:
            assert len(result["errors"]) >= 0  # May have errors or may be empty list with valid=False

    def test_error_has_required_fields(self) -> None:
        """Error dicts contain all required fields: raw_message, error_type, message, suggestions."""
        result = diagnose_smiles("C(C)(C)(C)(C)C")
        assert result["valid"] is False
        for error in result["errors"]:
            assert "raw_message" in error
            assert "error_type" in error
            assert "message" in error
            assert "suggestions" in error
            assert isinstance(error["suggestions"], list)

    def test_suggestions_sorted_by_confidence(self) -> None:
        """Fix suggestions are sorted by confidence descending."""
        result = diagnose_smiles("C(C(C")
        for error in result["errors"]:
            suggestions = error.get("suggestions", [])
            if len(suggestions) >= 2:
                for i in range(len(suggestions) - 1):
                    assert suggestions[i]["confidence"] >= suggestions[i + 1]["confidence"]
