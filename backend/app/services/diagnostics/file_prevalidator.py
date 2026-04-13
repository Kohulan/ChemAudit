"""
File Pre-Validator Service (DIAG-05).

Provides structural pre-validation for SDF and CSV files before molecule parsing.
Reports issues with absolute file line numbers, block integrity, and encoding.
"""

import importlib.util
import re
import sys
from pathlib import Path

# Import SDF_MARKERS and SUSPICIOUS_PATTERNS directly from file_parser module.
# Using importlib bypasses app.services.batch.__init__.py which imports Celery tasks
# and triggers the full app config/Redis stack — not appropriate for a lightweight service.
_fp_path = Path(__file__).parent.parent / "batch" / "file_parser.py"
if _fp_path.exists() and "app.services.batch.file_parser" not in sys.modules:
    _fp_spec = importlib.util.spec_from_file_location("app.services.batch.file_parser", _fp_path)
    _fp_mod = importlib.util.module_from_spec(_fp_spec)
    sys.modules["app.services.batch.file_parser"] = _fp_mod
    _fp_spec.loader.exec_module(_fp_mod)

from app.services.batch.file_parser import SUSPICIOUS_PATTERNS  # noqa: E402

# Common SMILES column name variants (case-sensitive check against normalized header)
_SMILES_COLUMN_NAMES = {"smiles", "canonical_smiles", "smi"}


def _scan_suspicious(content: bytes) -> list[dict]:
    """Scan raw bytes for security-suspicious patterns.

    Args:
        content: Raw file bytes.

    Returns:
        List of issue dicts with severity=error for each suspicious pattern found.
    """
    issues = []
    sample = content[:102400]
    for pattern in SUSPICIOUS_PATTERNS:
        if re.search(pattern, sample, re.IGNORECASE):
            pattern_str = pattern.decode("utf-8", errors="replace") if isinstance(pattern, bytes) else str(pattern)
            issues.append({
                "block": None,
                "line": 1,
                "issue_type": "suspicious_content",
                "severity": "error",
                "description": f"Security: suspicious pattern detected — {pattern_str}",
            })
    return issues


def prevalidate_sdf(content: bytes) -> dict:
    """Pre-validate an SDF file for structural integrity issues.

    Checks performed:
    - Security: suspicious patterns (XSS, injection)
    - Encoding: UTF-8 decode with replacement character detection
    - Per-block: M END presence (absolute line numbers, 1-indexed)
    - Per-block: counts line format (atoms/bonds line)

    Args:
        content: Raw SDF file bytes.

    Returns:
        Dict with keys:
            file_type (str): Always "sdf".
            total_blocks (int): Number of molecule blocks found.
            issue_count (int): Total number of issues detected.
            issues (list[dict]): Each issue has block, line, issue_type, severity, description.
            valid (bool): True if no error-severity issues found.
    """
    issues: list[dict] = []

    # Security check
    suspicious = _scan_suspicious(content)
    issues.extend(suspicious)

    # Decode with replacement to detect encoding problems
    text = content.decode("utf-8", errors="replace")
    replacement_count = text.count("\ufffd")
    if replacement_count > 0:
        issues.append({
            "block": None,
            "line": 1,
            "issue_type": "encoding_error",
            "severity": "warning",
            "description": (
                f"File contains {replacement_count} character(s) that could not be decoded as UTF-8; "
                "consider re-saving as UTF-8"
            ),
        })

    # Split on $$$$ record separator — last element may be empty trailing string
    raw_blocks = text.split("$$$$")

    # Count the separator lines to track absolute line offsets per block
    # We re-split on newlines for per-block offset tracking
    # Build a map from block index to starting line offset (1-indexed)
    block_start_lines: list[int] = []
    cumulative_line = 1
    for block_text in raw_blocks[:-1]:  # Exclude trailing empty element
        block_start_lines.append(cumulative_line)
        # Number of lines this block occupies plus the $$$$ separator line
        block_line_count = block_text.count("\n") + 1  # +1 for the $$$$ line itself
        cumulative_line += block_line_count

    # Validate each block
    total_blocks = len(raw_blocks) - 1  # Last element is trailing content after last $$$$
    if total_blocks <= 0:
        # No $$$$ found — treat the whole file as one block
        total_blocks = 1
        raw_blocks = [text, ""]
        block_start_lines = [1]

    for block_idx, block_text in enumerate(raw_blocks[:total_blocks]):
        block_num = block_idx + 1
        file_line_offset = block_start_lines[block_idx] if block_idx < len(block_start_lines) else 1

        block_lines = block_text.split("\n")

        # Check for M  END presence
        m_end_present = any("M  END" in line for line in block_lines)
        if not m_end_present:
            issues.append({
                "block": block_num,
                "line": file_line_offset,
                "issue_type": "missing_m_end",
                "severity": "error",
                "description": f"Block {block_num}: missing 'M  END' terminator (expected near line {file_line_offset})",
            })

        # Check counts line (line index 3, i.e., 4th line of the block, 1-indexed = file_line_offset + 3)
        if len(block_lines) >= 4:
            counts_line = block_lines[3]
            parts = counts_line.split()
            if len(parts) < 2:
                issues.append({
                    "block": block_num,
                    "line": file_line_offset + 3,
                    "issue_type": "malformed_count_line",
                    "severity": "warning",
                    "description": (
                        f"Block {block_num}: malformed counts line at file line {file_line_offset + 3} — "
                        "expected at least 2 fields (atom count, bond count)"
                    ),
                })

    error_count = sum(1 for i in issues if i["severity"] == "error")

    return {
        "file_type": "sdf",
        "total_blocks": total_blocks,
        "issue_count": len(issues),
        "issues": issues,
        "valid": error_count == 0,
    }


def prevalidate_csv(content: bytes) -> dict:
    """Pre-validate a CSV file for structural issues.

    Checks performed:
    - Encoding: UTF-8 decode; fallback to Latin-1 with warning
    - Security: suspicious patterns
    - Header: SMILES column detection
    - Rows: empty row detection

    Args:
        content: Raw CSV file bytes.

    Returns:
        Dict with keys:
            file_type (str): Always "csv".
            total_rows (int): Number of data rows (excluding header).
            encoding (str): Detected encoding used ("utf-8" or "latin-1").
            issue_count (int): Total number of issues detected.
            issues (list[dict]): Each issue has block, line, issue_type, severity, description.
            valid (bool): True if no error-severity issues found.
    """
    issues: list[dict] = []
    encoding = "utf-8"

    # Try UTF-8 decode first, then fall back to Latin-1
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        try:
            text = content.decode("latin-1")
            encoding = "latin-1"
            issues.append({
                "block": None,
                "line": 1,
                "issue_type": "encoding_fallback",
                "severity": "warning",
                "description": (
                    "File could not be decoded as UTF-8; fell back to Latin-1 encoding. "
                    "Consider re-saving as UTF-8 for maximum compatibility."
                ),
            })
        except UnicodeDecodeError:
            return {
                "file_type": "csv",
                "total_rows": 0,
                "encoding": "unknown",
                "issue_count": 1,
                "issues": [{
                    "block": None,
                    "line": 1,
                    "issue_type": "encoding_error",
                    "severity": "error",
                    "description": "File could not be decoded as UTF-8 or Latin-1",
                }],
                "valid": False,
            }

    # Security check on raw bytes
    suspicious = _scan_suspicious(content)
    issues.extend(suspicious)

    # Parse lines
    lines = text.splitlines()
    if not lines:
        issues.append({
            "block": None,
            "line": 1,
            "issue_type": "empty_file",
            "severity": "error",
            "description": "File is empty",
        })
        return {
            "file_type": "csv",
            "total_rows": 0,
            "encoding": encoding,
            "issue_count": len(issues),
            "issues": issues,
            "valid": False,
        }

    header_line = lines[0]
    data_lines = lines[1:]
    total_rows = len(data_lines)

    # Check for SMILES column in header (case-insensitive)
    # Detect delimiter from header
    delimiter = "\t" if header_line.count("\t") >= header_line.count(",") else ","
    header_columns = [col.strip() for col in header_line.split(delimiter)]

    # --- Duplicate column detection (case-insensitive) ---
    seen_columns: dict[str, list[int]] = {}
    for pos, col_name in enumerate(header_columns, start=1):
        key = col_name.strip().lower()
        if key:
            seen_columns.setdefault(key, []).append(pos)

    for col_key, positions in seen_columns.items():
        if len(positions) > 1:
            issues.append({
                "block": None,
                "line": 1,
                "issue_type": "duplicate_columns",
                "severity": "warning",
                "description": (
                    f"Column '{col_key}' appears {len(positions)} times "
                    f"(positions {', '.join(str(p) for p in positions)}). "
                    f"Duplicate headers may cause data loss during processing."
                ),
            })

    header_lower = {col.lower() for col in header_columns}

    if not (_SMILES_COLUMN_NAMES & header_lower):
        issues.append({
            "block": None,
            "line": 1,
            "issue_type": "missing_smiles_column",
            "severity": "warning",
            "description": (
                "No SMILES column detected in header. "
                f"Expected one of: smiles, SMILES, Smiles, canonical_smiles, smi. "
                f"Found columns: {', '.join(header_columns[:10])}"
            ),
        })

    # Check for empty rows
    empty_row_count = sum(1 for line in data_lines if not line.strip())
    if empty_row_count > 0:
        issues.append({
            "block": None,
            "line": None,
            "issue_type": "empty_rows",
            "severity": "warning",
            "description": f"{empty_row_count} empty row(s) detected in the data section",
        })

    error_count = sum(1 for i in issues if i["severity"] == "error")

    return {
        "file_type": "csv",
        "total_rows": total_rows,
        "encoding": encoding,
        "issue_count": len(issues),
        "issues": issues,
        "valid": error_count == 0,
    }
