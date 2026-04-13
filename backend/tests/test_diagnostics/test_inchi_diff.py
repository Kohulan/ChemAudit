"""
Unit tests for InChI Layer Diff Service (DIAG-02).

Tests cover:
- Identical InChI strings returning identical=True
- Stereo layer differences between enantiomers
- Formula layer differences between different molecules
- Missing layer handling (one InChI has it, other does not)
- Non-standard InChI version prefix parsing
"""


from app.services.diagnostics.inchi_diff import diff_inchi_layers, parse_inchi_layers

# Standard test InChI strings
ETHANOL_INCHI = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
METHANOL_INCHI = "InChI=1S/CH4O/c1-2/h2H,1H3"

# Stereo isomers of lactic acid
L_LACTIC_INCHI = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1"
D_LACTIC_INCHI = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m1/s1"


class TestIdenticalInChI:
    """Tests for identical InChI comparison."""

    def test_identical_inchi(self) -> None:
        """Same InChI for both returns identical=True with all layer_rows matching."""
        result = diff_inchi_layers(ETHANOL_INCHI, ETHANOL_INCHI)
        assert result["identical"] is True
        assert all(row["match"] for row in result["layer_rows"])

    def test_identical_inchi_layer_rows_present(self) -> None:
        """Identical InChI result has layer_rows with expected layers."""
        result = diff_inchi_layers(ETHANOL_INCHI, ETHANOL_INCHI)
        layer_names = [row["layer"] for row in result["layer_rows"]]
        assert "formula" in layer_names
        assert "connections" in layer_names
        assert "hydrogens" in layer_names


class TestStereoDiff:
    """Tests for stereo layer differences."""

    def test_stereo_layer_diff(self) -> None:
        """Two enantiomers differing in stereo layer return identical=False."""
        result = diff_inchi_layers(L_LACTIC_INCHI, D_LACTIC_INCHI)
        assert result["identical"] is False
        # At least one layer_row should have match=False
        mismatched = [row for row in result["layer_rows"] if not row["match"]]
        assert len(mismatched) > 0, "Expected at least one mismatched layer for stereoisomers"

    def test_stereo_layer_diff_has_stereo_rows(self) -> None:
        """Stereo layer rows are included in the comparison table."""
        result = diff_inchi_layers(L_LACTIC_INCHI, D_LACTIC_INCHI)
        layer_names = [row["layer"] for row in result["layer_rows"]]
        # L/D lactic acid differ in stereo_parity (m layer)
        assert "stereo_parity" in layer_names


class TestFormulaDiff:
    """Tests for formula layer differences."""

    def test_formula_diff(self) -> None:
        """Two different molecules return identical=False with formula match=False."""
        result = diff_inchi_layers(ETHANOL_INCHI, METHANOL_INCHI)
        assert result["identical"] is False
        formula_rows = [row for row in result["layer_rows"] if row["layer"] == "formula"]
        assert len(formula_rows) == 1
        assert formula_rows[0]["match"] is False

    def test_formula_values_populated(self) -> None:
        """Formula layer row has non-None values for both InChI strings."""
        result = diff_inchi_layers(ETHANOL_INCHI, METHANOL_INCHI)
        formula_rows = [row for row in result["layer_rows"] if row["layer"] == "formula"]
        assert len(formula_rows) == 1
        assert formula_rows[0]["value_a"] is not None
        assert formula_rows[0]["value_b"] is not None
        assert formula_rows[0]["value_a"] != formula_rows[0]["value_b"]


class TestMissingLayer:
    """Tests for InChI strings where one has a layer the other lacks."""

    def test_missing_layer(self) -> None:
        """One InChI with stereo layer, other without returns row with value_b None."""
        # Simple InChI without stereo layer
        simple_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)"
        # InChI with stereo layer
        stereo_inchi = L_LACTIC_INCHI

        result = diff_inchi_layers(simple_inchi, stereo_inchi)
        # Result should parse and contain rows for all layers present in either
        assert isinstance(result["layer_rows"], list)
        assert len(result["layer_rows"]) > 0
        # Should be non-identical since stereo differs
        assert result["identical"] is False

    def test_layers_a_and_b_populated(self) -> None:
        """layers_a and layers_b dicts are populated with parsed values."""
        result = diff_inchi_layers(ETHANOL_INCHI, METHANOL_INCHI)
        assert isinstance(result["layers_a"], dict)
        assert isinstance(result["layers_b"], dict)
        assert "formula" in result["layers_a"]
        assert "formula" in result["layers_b"]
        assert result["layers_a"]["formula"] == "C2H6O"
        assert result["layers_b"]["formula"] == "CH4O"


class TestNonStandardInChI:
    """Tests for non-standard InChI format handling."""

    def test_non_standard_inchi(self) -> None:
        """InChI starting with 'InChI=1/' (not 1S) parses without error."""
        non_standard = "InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = diff_inchi_layers(non_standard, ETHANOL_INCHI)
        # Should parse without raising an exception
        assert isinstance(result, dict)
        assert "identical" in result
        assert "layer_rows" in result
        assert "layers_a" in result

    def test_parse_inchi_layers_non_standard_version(self) -> None:
        """parse_inchi_layers extracts version from non-standard InChI prefix."""
        non_standard = "InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3"
        layers = parse_inchi_layers(non_standard)
        assert layers.get("version") == "1"
        assert layers.get("formula") == "C2H6O"
