"""
Unit tests for File Pre-Validator Service (DIAG-05).

Tests cover:
- Valid SDF: well-formed block, valid=True, no issues
- Missing M END: SDF block without terminator, valid=False
- Malformed count line: bad atom/bond count field
- Valid CSV: well-formed CSV with SMILES column, valid=True
- CSV encoding warning: Latin-1 bytes trigger encoding_fallback warning
- CSV missing SMILES column: header without SMILES/smi/canonical_smiles
"""


from app.services.diagnostics.file_prevalidator import prevalidate_csv, prevalidate_sdf

# Minimal valid SDF block — includes header, blank, comment, counts line, atoms, bonds, M END, $$$$
VALID_SDF_BLOCK = (
    b"\n"
    b"  Mrv2211 03272602002D\n"
    b"\n"
    b"  3  2  0  0  0  0            999 V2000\n"
    b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    b"    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    b"    3.0800    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    b"  1  2  1  0  0  0  0\n"
    b"  2  3  1  0  0  0  0\n"
    b"M  END\n"
    b"$$$$\n"
)

# SDF block missing M  END terminator
MISSING_M_END_SDF = (
    b"\n"
    b"  Mrv2211 03272602002D\n"
    b"\n"
    b"  3  2  0  0  0  0            999 V2000\n"
    b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    b"    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    b"    3.0800    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    b"  1  2  1  0  0  0  0\n"
    b"  2  3  1  0  0  0  0\n"
    b"$$$$\n"
)

# SDF block with malformed counts line (only 1 field on line 4)
MALFORMED_COUNT_LINE_SDF = (
    b"\n"
    b"  Mrv2211 03272602002D\n"
    b"\n"
    b"BAD\n"  # Only one field, should have at least 2
    b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    b"M  END\n"
    b"$$$$\n"
)

# Valid CSV with SMILES column
VALID_CSV = b"smiles,name\nCCO,ethanol\nCCN,ethylamine\n"

# CSV with standard header but Latin-1 encoded content
LATIN1_CSV = b"smiles,name\nCCO,caf\xe9\n"

# CSV without SMILES header
NO_SMILES_COLUMN_CSV = b"compound_id,formula\n1,C2H6O\n2,C2H7N\n"


class TestValidSDF:
    """Tests for well-formed SDF content."""

    def test_valid_sdf(self) -> None:
        """Well-formed SDF content returns valid=True and issue_count=0."""
        result = prevalidate_sdf(VALID_SDF_BLOCK)
        assert result["valid"] is True
        assert result["issue_count"] == 0
        assert result["issues"] == []
        assert result["file_type"] == "sdf"

    def test_valid_sdf_total_blocks(self) -> None:
        """Valid SDF returns correct total_blocks count."""
        result = prevalidate_sdf(VALID_SDF_BLOCK)
        assert result["total_blocks"] == 1

    def test_valid_sdf_multiple_blocks(self) -> None:
        """Multiple valid SDF blocks are counted correctly."""
        two_blocks = VALID_SDF_BLOCK + VALID_SDF_BLOCK
        result = prevalidate_sdf(two_blocks)
        assert result["total_blocks"] == 2
        assert result["valid"] is True

    def test_valid_sdf_has_required_keys(self) -> None:
        """SDF result contains all required response keys."""
        result = prevalidate_sdf(VALID_SDF_BLOCK)
        for key in ["file_type", "total_blocks", "issue_count", "issues", "valid"]:
            assert key in result, f"Missing key: {key}"


class TestMissingMEnd:
    """Tests for SDF blocks missing M  END terminator."""

    def test_missing_m_end(self) -> None:
        """SDF block without 'M  END' returns valid=False with missing_m_end issue."""
        result = prevalidate_sdf(MISSING_M_END_SDF)
        assert result["valid"] is False
        issue_types = [i["issue_type"] for i in result["issues"]]
        assert "missing_m_end" in issue_types

    def test_missing_m_end_severity_error(self) -> None:
        """missing_m_end issue has severity=error."""
        result = prevalidate_sdf(MISSING_M_END_SDF)
        missing_end_issues = [i for i in result["issues"] if i["issue_type"] == "missing_m_end"]
        assert len(missing_end_issues) > 0
        assert missing_end_issues[0]["severity"] == "error"

    def test_missing_m_end_has_block_number(self) -> None:
        """missing_m_end issue has block number set."""
        result = prevalidate_sdf(MISSING_M_END_SDF)
        missing_end_issues = [i for i in result["issues"] if i["issue_type"] == "missing_m_end"]
        assert missing_end_issues[0]["block"] == 1


class TestMalformedCountLine:
    """Tests for SDF blocks with malformed counts line."""

    def test_malformed_count_line(self) -> None:
        """SDF with bad counts line reports malformed_count_line issue."""
        result = prevalidate_sdf(MALFORMED_COUNT_LINE_SDF)
        issue_types = [i["issue_type"] for i in result["issues"]]
        assert "malformed_count_line" in issue_types

    def test_malformed_count_line_severity_warning(self) -> None:
        """malformed_count_line issue has severity=warning (not error)."""
        result = prevalidate_sdf(MALFORMED_COUNT_LINE_SDF)
        malformed_issues = [i for i in result["issues"] if i["issue_type"] == "malformed_count_line"]
        assert len(malformed_issues) > 0
        assert malformed_issues[0]["severity"] == "warning"


class TestValidCSV:
    """Tests for well-formed CSV content."""

    def test_valid_csv(self) -> None:
        """Well-formed CSV with SMILES column returns valid=True."""
        result = prevalidate_csv(VALID_CSV)
        assert result["valid"] is True
        assert result["file_type"] == "csv"

    def test_valid_csv_total_rows(self) -> None:
        """Valid CSV returns correct total_rows (excluding header)."""
        result = prevalidate_csv(VALID_CSV)
        assert result["total_rows"] == 2

    def test_valid_csv_encoding(self) -> None:
        """Valid UTF-8 CSV reports encoding='utf-8'."""
        result = prevalidate_csv(VALID_CSV)
        assert result["encoding"] == "utf-8"

    def test_valid_csv_has_required_keys(self) -> None:
        """CSV result contains all required response keys."""
        result = prevalidate_csv(VALID_CSV)
        for key in ["file_type", "total_rows", "encoding", "issue_count", "issues", "valid"]:
            assert key in result, f"Missing key: {key}"


class TestCSVEncodingWarning:
    """Tests for CSV files with non-UTF-8 encoding."""

    def test_csv_encoding_warning(self) -> None:
        """CSV with Latin-1 bytes produces an encoding_fallback warning."""
        result = prevalidate_csv(LATIN1_CSV)
        issue_types = [i["issue_type"] for i in result["issues"]]
        assert "encoding_fallback" in issue_types

    def test_csv_encoding_fallback_severity(self) -> None:
        """encoding_fallback issue has severity=warning."""
        result = prevalidate_csv(LATIN1_CSV)
        encoding_issues = [i for i in result["issues"] if i["issue_type"] == "encoding_fallback"]
        assert len(encoding_issues) > 0
        assert encoding_issues[0]["severity"] == "warning"

    def test_csv_encoding_fallback_still_valid(self) -> None:
        """CSV with encoding fallback is still valid (warning only, not error)."""
        result = prevalidate_csv(LATIN1_CSV)
        # Warnings don't affect validity — only errors do
        assert result["valid"] is True

    def test_csv_encoding_reports_latin1(self) -> None:
        """CSV falling back to Latin-1 reports encoding='latin-1'."""
        result = prevalidate_csv(LATIN1_CSV)
        assert result["encoding"] == "latin-1"


class TestCSVMissingSmiles:
    """Tests for CSV files lacking a SMILES column."""

    def test_csv_missing_smiles_column(self) -> None:
        """CSV without SMILES header produces missing_smiles_column issue."""
        result = prevalidate_csv(NO_SMILES_COLUMN_CSV)
        issue_types = [i["issue_type"] for i in result["issues"]]
        assert "missing_smiles_column" in issue_types

    def test_csv_missing_smiles_column_severity(self) -> None:
        """missing_smiles_column issue has severity=warning."""
        result = prevalidate_csv(NO_SMILES_COLUMN_CSV)
        smiles_issues = [i for i in result["issues"] if i["issue_type"] == "missing_smiles_column"]
        assert len(smiles_issues) > 0
        assert smiles_issues[0]["severity"] == "warning"
