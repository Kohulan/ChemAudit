"""
Tests for batch file parsing (SDF and CSV).
"""

import pytest

from app.services.batch.file_parser import (
    _is_headerless_delimited_file,
    _looks_like_smiles_value,
    detect_csv_columns,
    parse_csv,
    parse_sdf,
    validate_file_content_type,
)


class TestParseSDF:
    """Tests for SDF file parsing."""

    def test_parse_valid_sdf(self):
        """Test parsing a valid SDF with multiple molecules."""
        # Create a simple SDF content with 2 molecules
        sdf_content = b"""
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
> <_Name>
Ethanol

$$$$

     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
> <_Name>
Methanol

$$$$
"""
        molecules = parse_sdf(sdf_content)

        assert len(molecules) == 2
        assert molecules[0].name == "Ethanol"
        assert molecules[0].smiles != ""
        assert molecules[0].parse_error is None
        assert molecules[1].name == "Methanol"
        assert molecules[1].smiles != ""

    def test_parse_sdf_with_invalid_entry(self):
        """Test that invalid molecules don't crash the entire batch."""
        # SDF with one valid and one that will fail to parse
        sdf_content = b"""
     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
> <_Name>
Methanol

$$$$
INVALID MOLECULE DATA
$$$$
"""
        molecules = parse_sdf(sdf_content)

        # Should get at least the valid molecule
        assert len(molecules) >= 1
        valid_mols = [m for m in molecules if m.parse_error is None]
        assert len(valid_mols) >= 1
        assert valid_mols[0].name == "Methanol"

    def test_parse_sdf_generates_index_names(self):
        """Test that molecules without names get index-based names."""
        sdf_content = b"""
     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
"""
        molecules = parse_sdf(sdf_content)

        assert len(molecules) == 1
        assert molecules[0].name == "mol_0"
        assert molecules[0].index == 0


class TestParseCSV:
    """Tests for CSV file parsing."""

    def test_parse_valid_csv_with_smiles_column(self):
        """Test parsing CSV with standard SMILES column."""
        csv_content = b"""SMILES,Name,MolWeight
CCO,Ethanol,46.07
C,Methane,16.04
CC(=O)O,AceticAcid,60.05
"""
        molecules = parse_csv(csv_content)

        assert len(molecules) == 3
        assert molecules[0].smiles == "CCO"
        assert molecules[0].name == "Ethanol"
        # Note: CSV parser only reads SMILES and Name columns for efficiency
        # Properties are not populated from CSV files
        assert molecules[0].properties == {}

    def test_parse_csv_with_custom_smiles_column(self):
        """Test parsing CSV with non-standard SMILES column name."""
        csv_content = b"""compound_smiles,compound_id
CCO,ETH001
C,MTH001
"""
        molecules = parse_csv(csv_content, smiles_column="compound_smiles")

        assert len(molecules) == 2
        assert molecules[0].smiles == "CCO"

    def test_parse_csv_case_insensitive_column(self):
        """Test that SMILES column detection is case-insensitive."""
        csv_content = b"""smiles,name
CCO,Ethanol
"""
        molecules = parse_csv(csv_content, smiles_column="SMILES")

        assert len(molecules) == 1
        assert molecules[0].smiles == "CCO"

    def test_parse_csv_missing_smiles_column_raises(self):
        """Test that missing SMILES column raises ValueError."""
        csv_content = b"""Name,MolWeight
Ethanol,46.07
"""
        with pytest.raises(ValueError) as exc_info:
            parse_csv(csv_content, smiles_column="SMILES")

        assert "SMILES column 'SMILES' not found" in str(exc_info.value)

    def test_parse_csv_handles_empty_smiles(self):
        """Test that empty SMILES values create error entries."""
        csv_content = b"""SMILES,Name
CCO,Ethanol
,Empty
"""
        molecules = parse_csv(csv_content)

        assert len(molecules) == 2
        assert molecules[0].smiles == "CCO"
        assert molecules[0].parse_error is None
        # Empty SMILES should have a parse error
        assert molecules[1].parse_error is not None
        assert "Empty" in molecules[1].parse_error or "invalid" in molecules[1].parse_error.lower()

    def test_parse_csv_auto_detect_name_column(self):
        """Test auto-detection of name column."""
        csv_content = b"""SMILES,ID
CCO,ETH001
"""
        molecules = parse_csv(csv_content)

        assert molecules[0].name == "ETH001"

    def test_parse_tsv_content(self):
        """Test parsing tab-separated content (TSV format)."""
        tsv_content = b"""SMILES\tName\tMolWeight
CCO\tEthanol\t46.07
C\tMethane\t16.04
CC(=O)O\tAceticAcid\t60.05
"""
        molecules = parse_csv(tsv_content)

        assert len(molecules) == 3
        assert molecules[0].smiles == "CCO"
        assert molecules[0].name == "Ethanol"
        assert molecules[1].smiles == "C"
        assert molecules[1].name == "Methane"

    def test_parse_mixed_delimiter_preference(self):
        """Test that pandas correctly handles tab-separated data."""
        # Tab-separated with no commas in data
        tsv_content = b"""SMILES\tID
CCO\tETH001
C\tMTH001
"""
        molecules = parse_csv(tsv_content)

        assert len(molecules) == 2
        assert molecules[0].smiles == "CCO"
        assert molecules[0].name == "ETH001"


class TestDetectCSVColumns:
    """Tests for CSV column detection."""

    def test_detect_columns_returns_all_columns(self):
        """Test that all columns are returned."""
        csv_content = b"""SMILES,Name,MW,LogP
CCO,Ethanol,46.07,-0.31
"""
        result = detect_csv_columns(csv_content)

        assert "columns" in result
        assert "SMILES" in result["columns"]
        assert "Name" in result["columns"]
        assert "MW" in result["columns"]

    def test_detect_columns_suggests_smiles(self):
        """Test that SMILES column is suggested."""
        csv_content = b"""compound_smiles,name
CCO,Ethanol
"""
        result = detect_csv_columns(csv_content)

        assert result["suggested_smiles"] == "compound_smiles"

    def test_detect_columns_row_count(self):
        """Test row count estimation."""
        csv_content = b"""SMILES,Name
CCO,Ethanol
C,Methane
CC,Ethane
"""
        result = detect_csv_columns(csv_content)

        assert result["row_count_estimate"] == 3


class TestSingleColumnDelimitedFiles:
    """Tests for single-column and headerless delimited text files.

    Covers the bug where a TSV/CSV/TXT containing only SMILES (one per line,
    with or without a header row) was incorrectly rejected by the batch
    uploader with 'CSV file must use comma or tab delimiters'.
    """

    def test_single_column_tsv_with_header(self):
        """A TSV with just a SMILES header and one SMILES per line parses."""
        content = b"SMILES\nCCO\nCC(=O)O\nc1ccccc1\n"
        molecules = parse_csv(content)

        assert len(molecules) == 3
        assert [m.smiles for m in molecules] == ["CCO", "CC(=O)O", "c1ccccc1"]
        assert all(m.parse_error is None for m in molecules)

    def test_single_column_csv_with_header(self):
        """A CSV with just a SMILES header and one SMILES per line parses."""
        content = b"SMILES\nCCO\nCC(=O)O\nc1ccccc1\n"
        molecules = parse_csv(content)

        assert len(molecules) == 3
        assert [m.smiles for m in molecules] == ["CCO", "CC(=O)O", "c1ccccc1"]

    def test_headerless_single_column(self):
        """A file with no header, just SMILES per line, synthesizes SMILES column."""
        content = b"CCO\nCC(=O)O\nc1ccccc1\nCCCCCC\n"
        molecules = parse_csv(content)

        assert len(molecules) == 4
        assert [m.smiles for m in molecules] == [
            "CCO",
            "CC(=O)O",
            "c1ccccc1",
            "CCCCCC",
        ]
        assert all(m.parse_error is None for m in molecules)

    def test_headerless_multi_column_tab(self):
        """Headerless TSV with SMILES + extra column uses SMILES col_0."""
        content = b"CCO\tethanol\nCC(=O)O\tacetic\nc1ccccc1\tbenzene\n"
        molecules = parse_csv(content)

        assert len(molecules) == 3
        assert molecules[0].smiles == "CCO"
        assert molecules[1].smiles == "CC(=O)O"
        assert molecules[2].smiles == "c1ccccc1"

    def test_crlf_line_endings(self):
        """A file with Windows CRLF line endings parses correctly."""
        content = b"SMILES\r\nCCO\r\nCC(=O)O\r\n"
        molecules = parse_csv(content)

        assert len(molecules) == 2
        assert molecules[0].smiles == "CCO"
        assert molecules[1].smiles == "CC(=O)O"

    def test_headerless_crlf(self):
        """Headerless file with CRLF endings still detected as headerless."""
        content = b"CCO\r\nCC(=O)O\r\nc1ccccc1\r\n"
        molecules = parse_csv(content)

        assert len(molecules) == 3
        assert molecules[0].smiles == "CCO"

    def test_detect_columns_single_column_with_header(self):
        """Column detection on single-column headered TSV."""
        content = b"SMILES\nCCO\nCC(=O)O\n"
        result = detect_csv_columns(content)

        assert result["columns"] == ["SMILES"]
        assert result["suggested_smiles"] == "SMILES"
        assert result["column_samples"].get("SMILES") == "CCO"
        assert result["row_count_estimate"] == 2

    def test_detect_columns_headerless(self):
        """Column detection synthesizes SMILES for a headerless file."""
        content = b"CCO\nCC(=O)O\nc1ccccc1\n"
        result = detect_csv_columns(content)

        assert "SMILES" in result["columns"]
        assert result["suggested_smiles"] == "SMILES"
        # Row count should NOT subtract a non-existent header
        assert result["row_count_estimate"] == 3

    def test_headered_file_with_valid_smiles_column_is_not_headerless(self):
        """Standard CSV with 'SMILES' header must not be misclassified."""
        content = b"SMILES,Name\nCCO,ethanol\nCC(=O)O,acetic\n"
        molecules = parse_csv(content)

        assert len(molecules) == 2
        assert molecules[0].smiles == "CCO"
        assert molecules[0].name == "ethanol"

    def test_headerless_detection_helper_positive(self):
        """_is_headerless_delimited_file returns True for clear headerless data."""
        assert _is_headerless_delimited_file(b"CCO\nCC(=O)O\nc1ccccc1\n", ",")

    def test_headerless_detection_helper_rejects_header(self):
        """_is_headerless_delimited_file returns False when first row is a header."""
        assert not _is_headerless_delimited_file(b"SMILES\nCCO\nCC(=O)O\n", ",")
        assert not _is_headerless_delimited_file(b"SMILES,Name\nCCO,ethanol\n", ",")

    def test_headerless_detection_helper_requires_two_valid_rows(self):
        """First row parseable alone isn't enough — need two to avoid false positives."""
        # Second line is intentionally invalid — should NOT be classified headerless
        assert not _is_headerless_delimited_file(b"CCO\nNOT-A-SMILES-!@#\n", ",")

    def test_looks_like_smiles_rejects_header_keywords(self):
        """Common column names must not be detected as SMILES even if RDKit parses them."""
        # 'N' parses as nitrogen but is a common column-name letter
        assert not _looks_like_smiles_value("smiles")
        assert not _looks_like_smiles_value("SMILES")
        assert not _looks_like_smiles_value("Name")
        assert not _looks_like_smiles_value("mol")
        assert not _looks_like_smiles_value("")

    def test_looks_like_smiles_accepts_real_smiles(self):
        """Valid SMILES strings are detected."""
        assert _looks_like_smiles_value("CCO")
        assert _looks_like_smiles_value("CC(=O)O")
        assert _looks_like_smiles_value("c1ccccc1")

    def test_looks_like_smiles_rejects_whitespace(self):
        """SMILES never contain whitespace; reject anything that does."""
        assert not _looks_like_smiles_value("C C O")
        assert not _looks_like_smiles_value("CCO ethanol")

    def test_headerless_user_selects_first_row_cell_as_smiles_column(self):
        """Frontend shows first row as column names for headerless files. Users
        pick e.g. ``CCO`` as the SMILES column; backend must accept that.
        """
        content = b"CCO\nCC(=O)O\nc1ccccc1\n"
        # Simulate UI sending the first-row value as the column name
        molecules = parse_csv(content, smiles_column="CCO")
        assert len(molecules) == 3
        assert molecules[0].smiles == "CCO"
        assert molecules[1].smiles == "CC(=O)O"

    def test_headerless_multi_column_user_picks_by_position(self):
        """Headerless multi-column file: user's column selection is mapped by
        position in the first data row, not by synthesized name."""
        content = b"CCO\tethanol\nCC(=O)O\tacetic_acid\nc1ccccc1\tbenzene\n"
        molecules = parse_csv(content, smiles_column="CCO", name_column="ethanol")
        assert len(molecules) == 3
        assert molecules[0].smiles == "CCO"
        assert molecules[0].name == "ethanol"
        assert molecules[1].name == "acetic_acid"

    def test_headerless_user_picks_synthesized_name(self):
        """The synthesized ``SMILES`` name must also be accepted directly."""
        content = b"CCO\nCC(=O)O\n"
        molecules = parse_csv(content, smiles_column="SMILES")
        assert len(molecules) == 2
        assert molecules[0].smiles == "CCO"


class TestValidateFileContentType:
    """Tests for the pre-parse content validator."""

    def test_accepts_single_column_tsv(self):
        """A TSV with only SMILES per line (no tabs) must be accepted."""
        content = b"SMILES\nCCO\nCC(=O)O\n"
        is_valid, err = validate_file_content_type(content, "csv", "data.tsv")
        assert is_valid, f"Expected valid, got error: {err}"

    def test_accepts_headerless_tsv(self):
        """A headerless TSV with raw SMILES must be accepted."""
        content = b"CCO\nCC(=O)O\nc1ccccc1\n"
        is_valid, err = validate_file_content_type(content, "csv", "data.tsv")
        assert is_valid, f"Expected valid, got error: {err}"

    def test_accepts_standard_csv(self):
        """Standard multi-column CSV still accepted (regression)."""
        content = b"SMILES,Name\nCCO,ethanol\n"
        is_valid, err = validate_file_content_type(content, "csv", "data.csv")
        assert is_valid

    def test_accepts_crlf(self):
        """CRLF-terminated single-column file accepted."""
        content = b"SMILES\r\nCCO\r\nCC(=O)O\r\n"
        is_valid, err = validate_file_content_type(content, "csv", "data.tsv")
        assert is_valid

    def test_rejects_empty_file(self):
        """Empty file is rejected."""
        is_valid, err = validate_file_content_type(b"", "csv", "empty.csv")
        assert not is_valid

    def test_rejects_whitespace_only(self):
        """File with only whitespace is rejected."""
        is_valid, err = validate_file_content_type(b"\n\n   \n", "csv", "blank.csv")
        assert not is_valid

    def test_rejects_binary_content(self):
        """Binary content (non-printable bytes) is rejected."""
        content = b"SMILES\nCCO\n" + (b"\x01\x02\x03\x04" * 30)
        is_valid, err = validate_file_content_type(content, "csv", "evil.csv")
        assert not is_valid

    def test_rejects_executable_disguise(self):
        """Windows EXE pretending to be CSV is rejected."""
        is_valid, err = validate_file_content_type(b"MZ\x90\x00" + b"payload", "csv", "fake.csv")
        assert not is_valid

    def test_rejects_script_injection(self):
        """Script injection attempt is rejected."""
        is_valid, err = validate_file_content_type(
            b"SMILES\n<script>alert(1)</script>\n", "csv", "xss.csv"
        )
        assert not is_valid
