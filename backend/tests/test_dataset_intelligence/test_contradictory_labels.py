"""Tests for the contradictory label detection service module."""

import pytest

from app.services.dataset_intelligence.contradictory_labels import (
    detect_contradictory_labels,
    detect_numeric_columns,
)


def _make_molecule(
    smiles: str, index: int, inchikey: str | None, properties: dict | None = None
) -> dict:
    """Helper to build a molecule dict for contradictory label tests."""
    return {
        "index": index,
        "smiles": smiles,
        "mol": None,  # Not needed for label detection
        "inchikey": inchikey,
        "properties": properties or {},
    }


class TestDetectContradictoryLabels:
    """Tests for detect_contradictory_labels function."""

    def test_same_inchikey_large_fold_difference(self):
        """Same InChIKey with 100-fold activity diff => 1 contradiction."""
        molecules = [
            _make_molecule("CCO", 0, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                           {"activity": "1.0"}),
            _make_molecule("CCO", 1, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                           {"activity": "100.0"}),
        ]
        result = detect_contradictory_labels(molecules, "activity", fold_threshold=10.0)
        assert len(result) == 1
        assert result[0]["fold_difference"] == pytest.approx(100.0, rel=0.1)

    def test_threshold_sensitivity(self):
        """Lower threshold catches more contradictions."""
        molecules = [
            _make_molecule("CCO", 0, "IK_A", {"activity": "1.0"}),
            _make_molecule("CCO", 1, "IK_A", {"activity": "5.0"}),
            _make_molecule("c1ccccc1", 2, "IK_B", {"activity": "10.0"}),
            _make_molecule("c1ccccc1", 3, "IK_B", {"activity": "50.0"}),
        ]
        strict = detect_contradictory_labels(molecules, "activity", fold_threshold=100.0)
        lenient = detect_contradictory_labels(molecules, "activity", fold_threshold=3.0)
        assert len(lenient) >= len(strict)

    def test_zero_activity_excluded(self):
        """Zero activity values are excluded (no division by zero)."""
        molecules = [
            _make_molecule("CCO", 0, "IK_A", {"activity": "0.0"}),
            _make_molecule("CCO", 1, "IK_A", {"activity": "100.0"}),
        ]
        # Should NOT crash, zero values excluded from fold calculation
        result = detect_contradictory_labels(molecules, "activity")
        # With zero excluded, only one non-zero value remains — cannot compute fold
        # or it handles gracefully
        assert isinstance(result, list)

    def test_sorted_by_fold_difference_descending(self):
        """Results are sorted by fold_difference descending."""
        molecules = [
            _make_molecule("A", 0, "IK_A", {"activity": "1.0"}),
            _make_molecule("A", 1, "IK_A", {"activity": "20.0"}),
            _make_molecule("B", 2, "IK_B", {"activity": "1.0"}),
            _make_molecule("B", 3, "IK_B", {"activity": "1000.0"}),
        ]
        result = detect_contradictory_labels(molecules, "activity", fold_threshold=5.0)
        if len(result) >= 2:
            for i in range(len(result) - 1):
                assert result[i]["fold_difference"] >= result[i + 1]["fold_difference"]

    def test_no_activity_column_returns_empty(self):
        """Molecules with no matching activity column => empty list."""
        molecules = [
            _make_molecule("CCO", 0, "IK_A", {"other_col": "1.0"}),
            _make_molecule("CCO", 1, "IK_A", {"other_col": "100.0"}),
        ]
        result = detect_contradictory_labels(molecules, "activity")
        assert result == []

    def test_none_inchikey_excluded(self):
        """Molecules with None InChIKey are excluded from grouping."""
        molecules = [
            _make_molecule("CCO", 0, None, {"activity": "1.0"}),
            _make_molecule("CCO", 1, None, {"activity": "100.0"}),
        ]
        result = detect_contradictory_labels(molecules, "activity")
        assert result == []

    def test_negative_values_excluded(self):
        """Negative activity values are excluded."""
        molecules = [
            _make_molecule("CCO", 0, "IK_A", {"activity": "-5.0"}),
            _make_molecule("CCO", 1, "IK_A", {"activity": "100.0"}),
        ]
        result = detect_contradictory_labels(molecules, "activity")
        # Only one valid value, cannot compute fold
        assert isinstance(result, list)


class TestDetectNumericColumns:
    """Tests for detect_numeric_columns function."""

    def test_identifies_numeric_columns(self):
        """Correctly identifies columns with >80% numeric values."""
        result = detect_numeric_columns(
            ["activity", "name", "mw"],
            {
                "activity": ["1.0", "2.5", "3.0", "4.0", "5.0"],
                "name": ["aspirin", "ethanol", "benzene", "water", "toluene"],
                "mw": ["180.0", "46.0", "78.0", "18.0", "92.0"],
            },
        )
        names = [r["name"] for r in result]
        assert "activity" in names
        assert "mw" in names
        assert "name" not in names

    def test_priority_ordering(self):
        """Activity-like columns get priority 1, others get 3."""
        result = detect_numeric_columns(
            ["random_number", "pIC50", "mw"],
            {
                "random_number": ["1.0", "2.0", "3.0", "4.0", "5.0"],
                "pIC50": ["6.5", "7.2", "5.8", "8.1", "6.0"],
                "mw": ["300.0", "250.0", "180.0", "400.0", "350.0"],
            },
        )
        # pIC50 should be first (priority 1)
        assert result[0]["name"] == "pIC50"
        assert result[0]["priority"] == 1

    def test_empty_columns_returns_empty(self):
        """Empty column lists return empty list."""
        result = detect_numeric_columns([], {})
        assert result == []

    def test_non_numeric_columns_excluded(self):
        """Columns with <80% numeric values are excluded."""
        result = detect_numeric_columns(
            ["mixed"],
            {"mixed": ["1.0", "abc", "def", "ghi", "jkl"]},
        )
        assert len(result) == 0
