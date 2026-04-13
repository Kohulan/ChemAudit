"""Tests for the dataset diff service module."""


from app.services.dataset_intelligence.dataset_diff import (
    _compute_property_changes,
    compute_dataset_diff,
)


def _make_molecule(
    smiles: str, index: int, inchikey: str | None, properties: dict | None = None
) -> dict:
    """Helper to build a molecule dict for diff tests."""
    return {
        "index": index,
        "smiles": smiles,
        "mol": None,
        "inchikey": inchikey,
        "properties": properties or {},
    }


class TestComputeDatasetDiff:
    """Tests for compute_dataset_diff function."""

    def test_basic_diff_counts(self):
        """Primary=[A, B], comparison=[B, C, D] => added=2, removed=1."""
        primary = [
            _make_molecule("A", 0, "IK_A", {"mw": "100"}),
            _make_molecule("B", 1, "IK_B", {"mw": "200"}),
        ]
        comparison = [
            _make_molecule("B", 0, "IK_B", {"mw": "200"}),
            _make_molecule("C", 1, "IK_C", {"mw": "300"}),
            _make_molecule("D", 2, "IK_D", {"mw": "400"}),
        ]
        result = compute_dataset_diff(primary, comparison)
        assert result["added_count"] == 2  # C and D
        assert result["removed_count"] == 1  # A

    def test_modified_detection(self):
        """Same InChIKey with different MW shows up as modified."""
        primary = [
            _make_molecule("A", 0, "IK_A", {"mw": "100"}),
        ]
        comparison = [
            _make_molecule("A", 0, "IK_A", {"mw": "150"}),
        ]
        result = compute_dataset_diff(primary, comparison)
        assert result["modified_count"] == 1
        assert len(result["modified"]) == 1
        assert len(result["modified"][0]["changes"]) > 0
        assert result["modified"][0]["changes"][0]["column"] == "mw"

    def test_identical_datasets(self):
        """Identical datasets => all counts 0 except unchanged."""
        data = [
            _make_molecule("A", 0, "IK_A", {"mw": "100"}),
            _make_molecule("B", 1, "IK_B", {"mw": "200"}),
        ]
        result = compute_dataset_diff(data, data)
        assert result["added_count"] == 0
        assert result["removed_count"] == 0
        assert result["modified_count"] == 0
        assert result["unchanged_count"] == 2

    def test_none_inchikey_excluded(self):
        """Molecules with None InChIKey are excluded from diff."""
        primary = [
            _make_molecule("A", 0, None, {"mw": "100"}),
        ]
        comparison = [
            _make_molecule("B", 0, None, {"mw": "200"}),
        ]
        result = compute_dataset_diff(primary, comparison)
        # Both excluded — nothing to diff
        assert result["added_count"] == 0
        assert result["removed_count"] == 0
        assert result["unchanged_count"] == 0

    def test_column_mismatch_tracking(self):
        """Primary has 'extra1', comparison has 'extra2' => unique columns counted."""
        primary = [
            _make_molecule("A", 0, "IK_A", {"mw": "100", "extra1": "x"}),
        ]
        comparison = [
            _make_molecule("A", 0, "IK_A", {"mw": "100", "extra2": "y"}),
        ]
        result = compute_dataset_diff(primary, comparison)
        assert result["unique_columns_primary"] >= 1
        assert result["unique_columns_comparison"] >= 1
        # Only shared columns compared — "mw" is shared, unchanged
        assert result["modified_count"] == 0
        assert result["unchanged_count"] == 1


class TestComputePropertyChanges:
    """Tests for _compute_property_changes helper."""

    def test_detects_changes(self):
        """Detects property value changes between two molecules."""
        primary = {"properties": {"mw": "100", "logp": "2.0"}}
        comparison = {"properties": {"mw": "150", "logp": "2.0"}}
        changes = _compute_property_changes(primary, comparison, ["mw", "logp"])
        assert len(changes) == 1
        assert changes[0]["column"] == "mw"
        assert changes[0]["old_value"] == "100"
        assert changes[0]["new_value"] == "150"

    def test_no_changes(self):
        """Returns empty list when properties are identical."""
        primary = {"properties": {"mw": "100"}}
        comparison = {"properties": {"mw": "100"}}
        changes = _compute_property_changes(primary, comparison, ["mw"])
        assert changes == []
