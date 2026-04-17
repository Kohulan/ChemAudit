"""
File Parser Module

Parses SDF and CSV files into molecule data for batch processing.
Handles errors per-molecule without crashing the entire batch.
Includes security validations for production use.
"""

import io
import logging
import re
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Dict, Iterator, List, Optional, Tuple

import pandas as pd
from rdkit import Chem, RDLogger

logger = logging.getLogger(__name__)


@contextmanager
def _silence_rdkit() -> Iterator[None]:
    """Suppress RDKit log output for the duration of the block.

    Used when probing whether a string parses as SMILES — failures are
    expected here and the noisy default logging is not useful.
    """
    RDLogger.DisableLog("rdApp.*")
    try:
        yield
    finally:
        RDLogger.EnableLog("rdApp.*")


# =============================================================================
# File Content Type Validation
# =============================================================================

# Suspicious patterns that should never appear in chemical data files
SUSPICIOUS_PATTERNS = [
    rb"<script",  # JavaScript injection
    rb"javascript:",  # JavaScript protocol
    rb"<iframe",  # Iframe injection
    rb"onerror\s*=",  # Event handlers
    rb"onclick\s*=",  # Event handlers
    rb"onload\s*=",  # Event handlers
    rb"<object",  # Object tags
    rb"<embed",  # Embed tags
    rb"data:text/html",  # Data URI HTML
    rb"\x00",  # Null bytes (binary injection)
]

# SDF file markers
SDF_MARKERS = [
    b"M  END",  # End of molecule connection table
    b"$$$$",  # End of molecule record
    b"V2000",  # V2000 molfile format
    b"V3000",  # V3000 molfile format
]

# CSV-like content patterns
CSV_PATTERNS = [
    rb"^[^,\t]+[,\t][^,\t]+",  # Has delimiters
]


def validate_file_content_type(
    content: bytes,
    expected_type: str,
    filename: str = "",
) -> Tuple[bool, Optional[str]]:
    """
    Validate that file content matches the expected type.

    Performs magic byte / content analysis to detect:
    - Mismatched file types (e.g., executable disguised as CSV)
    - Malicious content injection attempts
    - Binary files masquerading as text

    Args:
        content: Raw file bytes
        expected_type: Expected file type ('sdf' or 'csv')
        filename: Original filename for logging

    Returns:
        Tuple of (is_valid, error_message)
        - (True, None) if content matches expected type
        - (False, error_message) if validation fails
    """
    if not content:
        return False, "Empty file"

    # Check for suspicious patterns in first 100KB
    sample = content[:102400]

    for pattern in SUSPICIOUS_PATTERNS:
        if re.search(pattern, sample, re.IGNORECASE):
            logger.warning(f"Suspicious pattern detected in file '{filename}': {pattern}")
            return False, "File contains potentially malicious content"

    # Check file magic bytes / signatures
    # Common executable/binary signatures that should never be in chemical files
    dangerous_signatures = [
        (b"MZ", "Windows executable"),  # PE/EXE
        (b"\x7fELF", "Linux executable"),  # ELF
        (b"PK\x03\x04", "ZIP archive"),  # ZIP (could be docx, xlsx, etc)
        (b"%PDF", "PDF document"),  # PDF
        (b"\xd0\xcf\x11\xe0", "OLE document"),  # DOC, XLS (old Office)
        (b"Rar!", "RAR archive"),  # RAR
        (b"\x1f\x8b", "GZIP compressed"),  # GZIP
        (b"BZh", "BZIP2 compressed"),  # BZIP2
    ]

    for sig, desc in dangerous_signatures:
        if content.startswith(sig):
            logger.warning(f"Invalid file type for '{filename}': detected {desc}")
            return False, f"Invalid file type: {desc} files are not allowed"

    # Validate based on expected type
    if expected_type.lower() == "sdf":
        return _validate_sdf_content(content, filename)
    elif expected_type.lower() == "csv":
        return _validate_csv_content(content, filename)
    else:
        return False, f"Unknown file type: {expected_type}"


def _validate_sdf_content(content: bytes, filename: str) -> Tuple[bool, Optional[str]]:
    """
    Validate that content looks like a valid SDF file.

    Args:
        content: Raw file bytes
        filename: Original filename

    Returns:
        Tuple of (is_valid, error_message)
    """
    # Check for SDF markers in the file
    found_markers = []
    for marker in SDF_MARKERS:
        if marker in content:
            found_markers.append(marker)

    # A valid SDF should have at least M END or $$$$
    if not any(marker in found_markers for marker in [b"M  END", b"$$$$"]):
        logger.warning(f"File '{filename}' does not contain SDF markers")
        return (
            False,
            "File does not appear to be a valid SDF file (missing molecule delimiters)",
        )

    # Check that content is primarily text-based
    # SDF files should be mostly printable ASCII
    sample = content[:10000]
    non_printable = sum(
        1
        for b in sample
        if b < 32 and b not in (9, 10, 13)  # Allow tab, newline, carriage return
    )
    if non_printable > len(sample) * 0.05:  # More than 5% non-printable
        logger.warning(f"File '{filename}' contains excessive non-printable characters")
        return False, "File contains too many non-printable characters for an SDF file"

    return True, None


def _validate_csv_content(content: bytes, filename: str) -> Tuple[bool, Optional[str]]:
    """
    Validate that content looks like a valid delimited text file.

    Accepts CSV, TSV, and plain-text files — including single-column files
    that contain no explicit delimiters in the header (e.g. a TSV with one
    SMILES per line). Rejects binary files, unreadable encodings, and
    content too small to be usable.

    Args:
        content: Raw file bytes
        filename: Original filename

    Returns:
        Tuple of (is_valid, error_message)
    """
    # Check that content is primarily text-based
    sample = content[:10000]

    # Count non-printable characters (excluding common whitespace)
    non_printable = sum(
        1
        for b in sample
        if b < 32 and b not in (9, 10, 13)  # Allow tab, newline, carriage return
    )
    if non_printable > len(sample) * 0.01:  # More than 1% non-printable
        logger.warning(f"File '{filename}' contains non-printable characters")
        return False, "File contains too many non-printable characters for a CSV file"

    # Try to decode as UTF-8 or Latin-1
    try:
        text = sample.decode("utf-8")
    except UnicodeDecodeError:
        try:
            text = sample.decode("latin-1")
        except UnicodeDecodeError:
            return (
                False,
                "File is not valid text (could not decode as UTF-8 or Latin-1)",
            )

    # Require at least one non-empty line. splitlines() handles all common
    # line endings (LF, CRLF, CR) so old-Mac-style files are not rejected here.
    non_empty_lines = [line for line in text.splitlines() if line.strip()]
    if not non_empty_lines:
        return False, "File appears to be empty"

    return True, None


def detect_suspicious_content(content: bytes) -> List[str]:
    """
    Scan content for suspicious patterns and return list of findings.

    This is useful for logging/auditing rather than blocking.

    Args:
        content: Raw file bytes

    Returns:
        List of suspicious pattern descriptions found
    """
    findings = []
    sample = content[:102400]

    patterns_with_names = [
        (rb"<script", "JavaScript tag"),
        (rb"javascript:", "JavaScript protocol"),
        (rb"<iframe", "Iframe tag"),
        (rb"onerror\s*=", "onerror event handler"),
        (rb"onclick\s*=", "onclick event handler"),
        (rb"onload\s*=", "onload event handler"),
        (rb"<object", "Object tag"),
        (rb"<embed", "Embed tag"),
        (rb"data:text/html", "Data URI with HTML"),
        (rb"\x00", "Null byte"),
    ]

    for pattern, name in patterns_with_names:
        if re.search(pattern, sample, re.IGNORECASE):
            findings.append(name)

    return findings


# Security limits
MAX_COLUMN_NAME_LENGTH = 256
MAX_PROPERTY_VALUE_LENGTH = 10000
MAX_SMILES_LENGTH = 5000
FORBIDDEN_COLUMN_PATTERNS = [
    r'[<>"\'\\/;`]',  # Characters that could be used in injection attacks
]


@dataclass
class MoleculeData:
    """Data extracted from a molecule in a batch file."""

    smiles: str
    name: Optional[str] = None
    index: int = 0
    properties: Dict[str, Any] = field(default_factory=dict)
    parse_error: Optional[str] = None


def _sanitize_string(value: str, max_length: int = MAX_PROPERTY_VALUE_LENGTH) -> str:
    """Sanitize a string value by truncating and removing control characters."""
    if not isinstance(value, str):
        value = str(value)
    # Remove control characters except newlines and tabs
    value = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]", "", value)
    # Truncate to max length
    return value[:max_length]


def _validate_column_name(name: str) -> bool:
    """Validate that a column name is safe to use."""
    if not name or len(name) > MAX_COLUMN_NAME_LENGTH:
        return False
    for pattern in FORBIDDEN_COLUMN_PATTERNS:
        if re.search(pattern, name):
            return False
    return True


def _detect_delimiter(content: bytes) -> str:
    """
    Detect the delimiter used in a delimited text file.

    Examines the first non-empty line to determine if the file uses
    tabs or commas as delimiters. Single-column files (no delimiter in
    the first line) fall through to the comma default; the parser
    treats them as a one-column file either way.

    Args:
        content: Raw file bytes

    Returns:
        Delimiter character ('\\t' for tab, ',' for comma)
    """
    # Get the first line (decode robustly, splitlines handles LF/CRLF/CR)
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        try:
            text = content.decode("latin-1")
        except UnicodeDecodeError:
            return ","  # Default to comma

    first_line = ""
    for line in text.splitlines():
        if line.strip():
            first_line = line
            break

    tab_count = first_line.count("\t")
    comma_count = first_line.count(",")

    # Prefer tab when at least as common as comma — tabs rarely appear in SMILES
    # or chemical names, so a tab in the header is a strong signal for TSV.
    if tab_count > 0 and tab_count >= comma_count:
        return "\t"
    return ","


# Column-name keywords that must never be treated as SMILES data, even when
# RDKit would technically parse them (e.g. "mol", "N", "I", "S"). Used by
# the headerless-file detector to avoid misclassifying common headers.
_HEADER_KEYWORDS = frozenset(
    {
        "smiles",
        "smi",
        "canonical_smiles",
        "isomeric_smiles",
        "input_smiles",
        "mol_smiles",
        "structure",
        "mol",
        "molecule",
        "name",
        "id",
        "compound",
        "compound_id",
        "compound_name",
        "title",
        "label",
        "identifier",
        "inchi",
        "inchikey",
        "inchi_key",
        "cas",
        "cas_number",
        "cas_rn",
        "activity",
        "value",
        "target",
        "assay",
        "ic50",
        "ki",
        "mw",
        "molweight",
        "logp",
        "tpsa",
    }
)


def _looks_like_smiles_value(value: str) -> bool:
    """Return True if ``value`` parses as a valid SMILES via RDKit.

    Used to detect whether the first row of a delimited file is data
    (headerless) or a column header. A value is considered SMILES-like
    only when:

    - It is non-empty, bounded, and contains no whitespace
      (SMILES strings never contain whitespace).
    - It is not a common column-name keyword (see ``_HEADER_KEYWORDS``) —
      even if RDKit would parse it (e.g. ``N`` as nitrogen).
    - RDKit produces a sanitized molecule with at least one atom.
    """
    if not value or not isinstance(value, str):
        return False
    cleaned = value.strip()
    if not cleaned or len(cleaned) > MAX_SMILES_LENGTH:
        return False
    if any(ch.isspace() for ch in cleaned):
        return False
    if cleaned.lower() in _HEADER_KEYWORDS:
        return False
    try:
        with _silence_rdkit():
            mol = Chem.MolFromSmiles(cleaned)
    except Exception:
        return False
    return mol is not None and mol.GetNumAtoms() > 0


def _is_headerless_delimited_file(content: bytes, delimiter: str) -> bool:
    """Detect whether a delimited text file has no header row.

    Returns True only when the file's first data column is unambiguously
    SMILES data rather than a column name. To avoid misclassifying a
    chemistry-like header (e.g. a single row named ``CO``), both the first
    AND second non-empty lines' leading cells must parse as valid SMILES.
    Files with a single non-empty line fall back to the first-line check.

    Args:
        content: Raw file bytes.
        delimiter: Detected column delimiter.

    Returns:
        True if the file should be read with ``header=None``.
    """
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        try:
            text = content.decode("latin-1")
        except UnicodeDecodeError:
            return False

    non_empty = [line for line in text.splitlines() if line.strip()]
    if not non_empty:
        return False

    first_cell = non_empty[0].split(delimiter)[0].strip()
    if not _looks_like_smiles_value(first_cell):
        return False

    if len(non_empty) == 1:
        return True

    second_cell = non_empty[1].split(delimiter)[0].strip()
    return _looks_like_smiles_value(second_cell)


@dataclass
class _DelimitedFileInfo:
    """Shape information about a delimited text file used by the parser helpers."""

    delimiter: str
    is_headerless: bool
    columns: List[str]
    read_kwargs: Dict[str, Any]


def _inspect_delimited_file(content: bytes) -> _DelimitedFileInfo:
    """Inspect a delimited text file to determine delimiter, header mode, columns.

    Produces the ``read_kwargs`` a subsequent ``pandas.read_csv`` call must use
    to parse the file consistently with the detected shape. For headerless
    files, column names are synthesized so downstream consumers always see a
    ``SMILES`` column in position 0.
    """
    delimiter = _detect_delimiter(content)
    is_headerless = _is_headerless_delimited_file(content, delimiter)

    base_kwargs: Dict[str, Any] = {
        "sep": delimiter,
        "dtype": str,
        "na_filter": False,
    }

    if is_headerless:
        probe = pd.read_csv(
            io.BytesIO(content),
            nrows=1,
            header=None,
            **base_kwargs,
        )
        n_cols = len(probe.columns)
        synthetic = ["SMILES"] + [f"col_{i}" for i in range(1, n_cols)]
        read_kwargs = {**base_kwargs, "header": None, "names": synthetic}
        return _DelimitedFileInfo(
            delimiter=delimiter,
            is_headerless=True,
            columns=synthetic,
            read_kwargs=read_kwargs,
        )

    header_probe = pd.read_csv(io.BytesIO(content), nrows=0, **base_kwargs)
    return _DelimitedFileInfo(
        delimiter=delimiter,
        is_headerless=False,
        columns=header_probe.columns.tolist(),
        read_kwargs=base_kwargs,
    )


def _match_column_case_insensitive(columns: List[str], target: str) -> Optional[str]:
    """Find a column by case-insensitive exact-name match, or None."""
    lowered = target.lower()
    for col in columns:
        if str(col).lower() == lowered:
            return col
    return None


_NAME_COLUMN_CANDIDATES: Tuple[str, ...] = (
    "Name",
    "name",
    "NAME",
    "ID",
    "id",
    "Compound",
    "compound",
    "Molecule",
    "molecule",
    "Title",
    "title",
)


def _auto_detect_name_column(columns: List[str]) -> Optional[str]:
    """Pick a likely name/ID column from a list of columns, or None."""
    for candidate in _NAME_COLUMN_CANDIDATES:
        if candidate in columns:
            return candidate
    return None


def _peek_first_row_cells(content: bytes, delimiter: str) -> List[str]:
    """Return the cells of the first non-empty line, split by ``delimiter``.

    Used to map user-supplied column selections (which the frontend derives
    from the first row of a headerless file) to synthesized column positions.
    """
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        try:
            text = content.decode("latin-1")
        except UnicodeDecodeError:
            return []

    for line in text.splitlines():
        if line.strip():
            return [cell.strip() for cell in line.split(delimiter)]
    return []


def _resolve_column_for_headerless(
    requested: str,
    synthesized_columns: List[str],
    first_row_cells: List[str],
) -> Optional[str]:
    """Map a user-supplied column name to a synthesized column for headerless files.

    The frontend's local column detector treats the first data row as
    column headers, so it shows the user cells like ``CCO`` or ``ethanol``
    as column names. The backend synthesizes ``SMILES`` / ``col_N``. When
    the file is headerless, we accept the user's selection if it matches
    either the synthesized name or the first-row cell at that same position.
    """
    if not requested:
        return None
    req = requested.strip().lower()
    for idx, synth in enumerate(synthesized_columns):
        if req == synth.lower():
            return synth
        if idx < len(first_row_cells) and req == first_row_cells[idx].lower():
            return synth
    return None


def _validate_smiles_format(smiles: str) -> Tuple[bool, Optional[str]]:
    """
    Basic validation of SMILES string format before RDKit parsing.

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not smiles or not isinstance(smiles, str):
        return False, "Empty or invalid SMILES"

    smiles = smiles.strip()

    if len(smiles) > MAX_SMILES_LENGTH:
        return False, f"SMILES too long (>{MAX_SMILES_LENGTH} chars)"

    # Check for obviously invalid characters (basic sanity check)
    # SMILES can contain: letters, numbers, @, #, =, -, +, [, ], (, ), /, \, ., %
    if re.search(r'[<>"\';`]', smiles):
        return False, "SMILES contains invalid characters"

    return True, None


def parse_sdf(file_content: bytes, max_file_size_mb: int = 500) -> List[MoleculeData]:
    """
    Parse an SDF file into a list of molecule data.

    Args:
        file_content: Raw bytes of the SDF file
        max_file_size_mb: Maximum file size in MB

    Returns:
        List of MoleculeData objects (includes entries with parse_error for invalid molecules)

    Raises:
        ValueError: If file exceeds size limit

    Note:
        Individual parse errors are captured per-molecule rather than failing the entire batch.
    """
    # Security: Check file size
    file_size_mb = len(file_content) / (1024 * 1024)
    if file_size_mb > max_file_size_mb:
        raise ValueError(
            f"File too large: {file_size_mb:.1f}MB exceeds limit of {max_file_size_mb}MB"
        )

    molecules = []
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(file_content))

    for idx, mol in enumerate(supplier):
        if mol is None:
            # Failed to parse this molecule
            molecules.append(
                MoleculeData(
                    smiles="",
                    name=f"mol_{idx}",
                    index=idx,
                    parse_error=f"Failed to parse molecule at index {idx}",
                )
            )
            continue

        try:
            # Get SMILES
            smiles = Chem.MolToSmiles(mol)

            # Extract name from _Name property or generate index-based
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{idx}"
            if not name.strip():
                name = f"mol_{idx}"
            name = _sanitize_string(name, max_length=256)

            # Extract all other properties with sanitization
            properties = {}
            for prop_name in mol.GetPropsAsDict():
                if prop_name != "_Name":
                    try:
                        prop_value = mol.GetProp(prop_name)
                        # Sanitize property name and value
                        safe_name = _sanitize_string(prop_name, max_length=MAX_COLUMN_NAME_LENGTH)
                        if _validate_column_name(safe_name):
                            properties[safe_name] = _sanitize_string(prop_value)
                    except Exception:
                        pass

            molecules.append(
                MoleculeData(
                    smiles=smiles,
                    name=name,
                    index=idx,
                    properties=properties,
                )
            )
        except Exception as e:
            molecules.append(
                MoleculeData(
                    smiles="",
                    name=f"mol_{idx}",
                    index=idx,
                    parse_error=f"Error extracting data: {str(e)[:200]}",
                )
            )

    return molecules


def parse_csv(
    file_content: bytes,
    smiles_column: str = "SMILES",
    name_column: Optional[str] = None,
    max_file_size_mb: int = 500,
) -> List[MoleculeData]:
    """
    Parse a CSV/TSV/TXT file into a list of molecule data.

    Automatically detects:
    - Whether the file uses comma or tab delimiters.
    - Whether the file has a header row, by probing the first two rows with
      RDKit. Headerless files have a ``SMILES`` column name synthesized
      so downstream callers need no special handling.

    Args:
        file_content: Raw bytes of the delimited text file
        smiles_column: Name of the column containing SMILES strings
        name_column: Optional name of the column containing molecule names
        max_file_size_mb: Maximum file size in MB

    Returns:
        List of MoleculeData objects

    Raises:
        ValueError: If the SMILES column is not found or file exceeds size limit

    Note:
        Individual parse errors are captured per-molecule.
    """
    # Security: Check file size
    file_size_mb = len(file_content) / (1024 * 1024)
    if file_size_mb > max_file_size_mb:
        raise ValueError(
            f"File too large: {file_size_mb:.1f}MB exceeds limit of {max_file_size_mb}MB"
        )

    try:
        info = _inspect_delimited_file(file_content)
    except Exception as exc:
        raise ValueError(f"Failed to parse CSV file: {str(exc)[:200]}")

    all_columns = info.columns

    # Security: Validate column count
    if len(all_columns) > 1000:
        raise ValueError("Too many columns in CSV (maximum 1000)")

    smiles_col_actual = _match_column_case_insensitive(all_columns, smiles_column)
    if smiles_col_actual is None and info.is_headerless:
        # The frontend's local column detector reads the first data row of a
        # headerless file as column names, so the user's ``smiles_column`` can
        # arrive as e.g. ``"CCO"`` instead of the synthesized ``"SMILES"``.
        # Map the user's selection to the synthesized column at the same
        # position in the first data row.
        first_row_cells = _peek_first_row_cells(file_content, info.delimiter)
        smiles_col_actual = _resolve_column_for_headerless(
            smiles_column, all_columns, first_row_cells
        )
    if smiles_col_actual is None:
        raise ValueError(f"SMILES column '{smiles_column}' not found in CSV")

    # Determine name column. Honour the user's explicit ``name_column`` first;
    # accept first-row-cell aliasing for headerless files the same way we do
    # for the SMILES column. For headered files without an explicit selection,
    # fall back to the keyword-based auto-detection.
    name_col_actual: Optional[str] = None
    if name_column:
        name_col_actual = _match_column_case_insensitive(all_columns, name_column)
        if name_col_actual is None and info.is_headerless:
            first_row_cells = _peek_first_row_cells(file_content, info.delimiter)
            name_col_actual = _resolve_column_for_headerless(
                name_column, all_columns, first_row_cells
            )
    elif not info.is_headerless:
        name_col_actual = _auto_detect_name_column(all_columns)

    # Read only the columns we need (much faster for wide files)
    cols_to_read = [smiles_col_actual]
    if name_col_actual:
        cols_to_read.append(name_col_actual)

    try:
        df = pd.read_csv(
            io.BytesIO(file_content),
            usecols=cols_to_read,
            **info.read_kwargs,
        )
    except Exception as exc:
        raise ValueError(f"Failed to parse CSV file: {str(exc)[:200]}")

    molecules: List[MoleculeData] = []
    for idx, row in df.iterrows():
        smiles = str(row[smiles_col_actual]).strip()

        is_valid, error_msg = _validate_smiles_format(smiles)
        if not is_valid:
            molecules.append(
                MoleculeData(
                    smiles="",
                    name=f"row_{idx}",
                    index=idx,
                    parse_error=error_msg or "Invalid SMILES",
                )
            )
            continue

        if name_col_actual and row[name_col_actual]:
            name = _sanitize_string(str(row[name_col_actual]).strip(), max_length=256)
        else:
            name = f"row_{idx}"

        molecules.append(
            MoleculeData(
                smiles=smiles,
                name=name,
                index=idx,
                properties={},
            )
        )

    return molecules


def _count_file_lines(content: bytes) -> int:
    """Count physical lines in a byte buffer, tolerant of CRLF/LF/CR endings.

    ``content.count(b"\\n")`` alone misses old-Mac CR-only files. Using the
    larger of LF-count and CR-count gives the correct line total for all
    three common line-ending styles.
    """
    lf = content.count(b"\n")
    cr = content.count(b"\r")
    total = max(lf, cr)
    # Files that end without a trailing terminator still contain that last line
    if content and not content.endswith((b"\n", b"\r")):
        total += 1
    return max(total, 0)


def detect_csv_columns(
    file_content: bytes,
    max_file_size_mb: int = 500,
) -> Dict[str, Any]:
    """
    Detect column names in a delimited text file for user selection.

    Automatically detects:
    - Whether the file uses comma or tab delimiters.
    - Whether the file has a header row. For headerless files a synthetic
      ``SMILES`` column is surfaced so the user sees a consistent choice.

    Args:
        file_content: Raw bytes of the delimited text file
        max_file_size_mb: Maximum file size in MB

    Returns:
        Dictionary with column info and suggestions

    Raises:
        ValueError: If file cannot be read or exceeds size limit
    """
    # Security: Check file size
    file_size_mb = len(file_content) / (1024 * 1024)
    if file_size_mb > max_file_size_mb:
        raise ValueError(
            f"File too large: {file_size_mb:.1f}MB exceeds limit of {max_file_size_mb}MB"
        )

    try:
        info = _inspect_delimited_file(file_content)
    except ValueError:
        raise
    except Exception as exc:
        logger.error(f"Error detecting CSV columns: {exc}")
        raise ValueError(f"Failed to read CSV file: {str(exc)[:100]}")

    try:
        df_preview = pd.read_csv(
            io.BytesIO(file_content),
            nrows=5,
            **info.read_kwargs,
        )
    except ValueError:
        raise
    except Exception as exc:
        logger.error(f"Error reading CSV preview: {exc}")
        raise ValueError(f"Failed to read CSV file: {str(exc)[:100]}")

    columns = [col for col in info.columns if _validate_column_name(col)]
    if not columns:
        raise ValueError("No valid columns found in CSV")

    # Suggest SMILES column. For headerless files the synthesized ``SMILES``
    # column is always position 0 and always the right answer.
    if info.is_headerless:
        suggested_smiles: Optional[str] = "SMILES" if "SMILES" in columns else None
    else:
        suggested_smiles = None
        smiles_keywords = [
            "smiles",
            "smi",
            "canonical_smiles",
            "isomeric_smiles",
            "structure",
        ]
        for col in columns:
            col_lower = str(col).lower()
            for keyword in smiles_keywords:
                if keyword in col_lower:
                    suggested_smiles = col
                    break
            if suggested_smiles:
                break

    # Suggest name/ID column only for headered files; synthesized ``col_N``
    # labels from headerless parsing carry no semantic meaning.
    suggested_name: Optional[str] = None
    if not info.is_headerless:
        name_keywords = ["name", "id", "compound", "molecule", "title", "identifier"]
        for col in columns:
            if col == suggested_smiles:
                continue
            col_lower = str(col).lower()
            for keyword in name_keywords:
                if keyword in col_lower:
                    suggested_name = col
                    break
            if suggested_name:
                break

    # First non-empty sample value per column (capped to first 20 columns)
    column_samples: Dict[str, str] = {}
    for col in columns[:20]:
        for val in df_preview[col]:
            if val and str(val).strip():
                column_samples[col] = _sanitize_string(str(val), max_length=100)
                break

    # Row count estimate. Subtract the header row only when the file has one.
    header_offset = 0 if info.is_headerless else 1
    if file_size_mb < 1:
        line_count = _count_file_lines(file_content)
        row_count = max(1, line_count - header_offset)
    else:
        try:
            sample_size = min(102400, len(file_content))
            sample_bytes = file_content[:sample_size]
            lines_in_sample = _count_file_lines(sample_bytes)
            if lines_in_sample > 1:
                avg_line_length = sample_size / lines_in_sample
                row_count = max(1, int(len(file_content) / avg_line_length) - header_offset)
            else:
                row_count = max(1, int(file_size_mb * 5000))
        except Exception:
            row_count = max(1, int(file_size_mb * 5000))

    return {
        "columns": columns,
        "suggested_smiles": suggested_smiles,
        "suggested_name": suggested_name,
        "column_samples": column_samples,
        "row_count_estimate": row_count,
        "file_size_mb": round(file_size_mb, 2),
    }
