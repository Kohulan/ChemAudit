"""
SMILES Error Diagnostics Service (DIAG-01).

Provides position-specific SMILES error extraction and fix suggestions using
a dual strategy: RDKit DetectChemistryProblems for parseable-but-problematic
SMILES, and log capture via rdBase.LogToPythonLogger for unparseable SMILES.
"""

import logging
import re
from typing import Optional

from rdkit import Chem, rdBase


# Module logger for capturing RDKit error messages
_rdkit_log_capture_logger = logging.getLogger("rdkit.log_capture")

# Mapping from error type to human-readable description
_ERROR_DESCRIPTIONS = {
    "unmatched_bracket": "Unmatched parenthesis in SMILES string",
    "valence_error": "Atom has an impossible valence for its element",
    "ring_closure_mismatch": "Ring closure digit appears unpaired or mismatched",
    "unknown_atom_symbol": "Unrecognized atom symbol in SMILES string",
    "invalid_charge": "Invalid charge specification on atom",
    "parse_error": "SMILES string could not be parsed by RDKit",
}


def _extract_position(error_msg: str) -> Optional[int]:
    """Extract character position from an RDKit error message.

    Args:
        error_msg: Raw RDKit error message string.

    Returns:
        Zero-indexed character position if found, else None.
    """
    match = re.search(r"at position (\d+)", error_msg, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None


def _classify_error(error_msg: str) -> str:
    """Classify an RDKit error message into a structured error type.

    Args:
        error_msg: Raw RDKit error message string.

    Returns:
        Error type string: one of unmatched_bracket, valence_error,
        ring_closure_mismatch, unknown_atom_symbol, invalid_charge, parse_error.
    """
    msg_lower = error_msg.lower()
    if "bracket" in msg_lower or "parenthes" in msg_lower or "unmatched" in msg_lower:
        return "unmatched_bracket"
    if "valence" in msg_lower:
        return "valence_error"
    if "ring" in msg_lower or "closure" in msg_lower:
        return "ring_closure_mismatch"
    if "atom" in msg_lower and ("symbol" in msg_lower or "type" in msg_lower or "unknown" in msg_lower):
        return "unknown_atom_symbol"
    if "charge" in msg_lower:
        return "invalid_charge"
    return "parse_error"


def _human_message(error_type: str) -> str:
    """Map an error type to a human-readable description.

    Args:
        error_type: Structured error type string.

    Returns:
        Human-readable description of the error.
    """
    return _ERROR_DESCRIPTIONS.get(error_type, "An error occurred while parsing the SMILES string.")


def _build_suggestions(
    error_type: str, smiles: str, position: Optional[int]
) -> list[dict]:
    """Build fix suggestions for a SMILES error.

    Args:
        error_type: Structured error type string.
        smiles: The original SMILES string.
        position: Zero-indexed position of the error, if known.

    Returns:
        List of suggestion dicts with keys: description, corrected_smiles, confidence.
        Sorted by confidence descending. Maximum 3 items.
    """
    suggestions: list[dict] = []

    if error_type == "unmatched_bracket":
        # Check round parentheses
        open_count = smiles.count("(")
        close_count = smiles.count(")")
        if open_count > close_count:
            diff = open_count - close_count
            corrected = smiles + ")" * diff
            suggestions.append({
                "description": f"Append {diff} closing parenthesis character(s) to balance the string",
                "corrected_smiles": corrected,
                "confidence": 0.8,
            })
        elif close_count > open_count:
            diff = close_count - open_count
            # Remove trailing close parentheses
            corrected = smiles.rstrip(")")
            if corrected.count(")") < open_count:
                # Only suggest if removing trailing ones fixes it
                corrected = smiles[:smiles.rfind(")")] + smiles[smiles.rfind(")") + 1:]
            suggestions.append({
                "description": f"Remove {diff} extra closing parenthesis character(s)",
                "corrected_smiles": None,
                "confidence": 0.6,
            })

        # Check square brackets
        sq_open = smiles.count("[")
        sq_close = smiles.count("]")
        if sq_open > sq_close:
            diff = sq_open - sq_close
            corrected = smiles + "]" * diff
            suggestions.append({
                "description": f"Append {diff} closing square bracket(s) to balance the atom specification",
                "corrected_smiles": corrected,
                "confidence": 0.75,
            })
        elif sq_close > sq_open:
            suggestions.append({
                "description": "Remove extra closing square bracket(s) from atom specification",
                "corrected_smiles": None,
                "confidence": 0.6,
            })

    elif error_type == "valence_error":
        suggestions.append({
            "description": "Check the valence of the atom at the indicated position; "
                           "consider adding explicit hydrogens with [AtomH] notation",
            "corrected_smiles": None,
            "confidence": 0.4,
        })
        suggestions.append({
            "description": "Verify the bond orders around the atom are chemically valid",
            "corrected_smiles": None,
            "confidence": 0.3,
        })

    elif error_type == "ring_closure_mismatch":
        # Find unpaired ring closure digits
        digits_found: dict[str, int] = {}
        i = 0
        while i < len(smiles):
            if smiles[i] == "%" and i + 2 < len(smiles) and smiles[i + 1:i + 3].isdigit():
                key = smiles[i:i + 3]
                digits_found[key] = digits_found.get(key, 0) + 1
                i += 3
            elif smiles[i].isdigit():
                key = smiles[i]
                digits_found[key] = digits_found.get(key, 0) + 1
                i += 1
            else:
                i += 1
        unpaired = [k for k, v in digits_found.items() if v % 2 != 0]
        if unpaired:
            suggestions.append({
                "description": f"Ring closure digit(s) {unpaired} appear an odd number of times; "
                               "each ring closure digit must appear exactly twice",
                "corrected_smiles": None,
                "confidence": 0.65,
            })

    elif error_type == "unknown_atom_symbol":
        suggestions.append({
            "description": "Check that all atom symbols are valid element symbols or "
                           "organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I)",
            "corrected_smiles": None,
            "confidence": 0.5,
        })

    elif error_type == "invalid_charge":
        suggestions.append({
            "description": "Ensure charge is specified as +, -, ++, --, or a number inside square brackets, "
                           "e.g. [NH4+] for ammonium",
            "corrected_smiles": None,
            "confidence": 0.5,
        })

    else:
        suggestions.append({
            "description": "Verify the SMILES string against standard SMILES notation rules",
            "corrected_smiles": None,
            "confidence": 0.3,
        })

    # Sort by confidence descending and limit to 3
    suggestions.sort(key=lambda s: s["confidence"], reverse=True)
    return suggestions[:3]


class _RDKitLogCapture(logging.Handler):
    """Logging handler that captures RDKit log messages into a list."""

    def __init__(self) -> None:
        super().__init__()
        self.records: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record.getMessage())


def diagnose_smiles(smiles: str) -> dict:
    """Diagnose a SMILES string for errors, providing position and fix suggestions.

    Uses a dual strategy:
    - For parseable-but-problematic SMILES: DetectChemistryProblems gives typed warnings/errors.
    - For unparseable SMILES: captures RDKit C++ error log output to extract position info.

    Args:
        smiles: Raw SMILES string to diagnose.

    Returns:
        Dict with keys:
            valid (bool): Whether the SMILES is valid.
            canonical_smiles (str | None): Canonical SMILES if valid, else None.
            warnings (list[dict]): List of warnings on valid SMILES.
            errors (list[dict]): List of errors on invalid SMILES, each with:
                raw_message, position (int | None), error_type, message, suggestions.
    """
    # --- Step 1: Attempt to parse without sanitization ---
    handler = _RDKitLogCapture()
    rdkit_logger = logging.getLogger("rdkit")

    rdBase.LogToPythonLogger()
    rdkit_logger.addHandler(handler)

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    finally:
        rdkit_logger.removeHandler(handler)

    if mol is None:
        # Completely unparseable — use captured log messages
        captured = handler.records if handler.records else []

        # If no log messages were captured, try again with sanitize=True to trigger them
        if not captured:
            handler2 = _RDKitLogCapture()
            rdBase.LogToPythonLogger()
            rdkit_logger.addHandler(handler2)
            try:
                Chem.MolFromSmiles(smiles)
            finally:
                rdkit_logger.removeHandler(handler2)
            captured = handler2.records

        errors = []
        if captured:
            for raw_msg in captured:
                error_type = _classify_error(raw_msg)
                position = _extract_position(raw_msg)
                errors.append({
                    "raw_message": raw_msg,
                    "position": position,
                    "error_type": error_type,
                    "message": _human_message(error_type),
                    "suggestions": _build_suggestions(error_type, smiles, position),
                })
        else:
            # Fallback: no messages captured — generic error
            error_type = "parse_error"
            errors.append({
                "raw_message": "SMILES could not be parsed",
                "position": None,
                "error_type": error_type,
                "message": _human_message(error_type),
                "suggestions": _build_suggestions(error_type, smiles, None),
            })

        return {
            "valid": False,
            "canonical_smiles": None,
            "warnings": [],
            "errors": errors,
        }

    # --- Step 2: Parseable — detect chemistry problems ---
    problems = Chem.DetectChemistryProblems(mol)

    if problems:
        # Has chemistry problems — treat as invalid
        errors = []
        for prob in problems:
            raw_msg = prob.Message()
            error_type = _classify_error(raw_msg)
            position = _extract_position(raw_msg)
            errors.append({
                "raw_message": raw_msg,
                "position": position,
                "error_type": error_type,
                "message": _human_message(error_type),
                "suggestions": _build_suggestions(error_type, smiles, position),
            })
        return {
            "valid": False,
            "canonical_smiles": None,
            "warnings": [],
            "errors": errors,
        }

    # --- Step 3: No problems — fully sanitize and return canonical SMILES ---
    try:
        Chem.SanitizeMol(mol)
        canonical = Chem.MolToSmiles(mol, canonical=True)
    except Exception as exc:
        error_type = "parse_error"
        raw_msg = str(exc)
        return {
            "valid": False,
            "canonical_smiles": None,
            "warnings": [],
            "errors": [{
                "raw_message": raw_msg,
                "position": _extract_position(raw_msg),
                "error_type": error_type,
                "message": _human_message(error_type),
                "suggestions": _build_suggestions(error_type, smiles, None),
            }],
        }

    # Collect DetectChemistryProblems warnings on the sanitized mol (informational)
    sanitized_problems = Chem.DetectChemistryProblems(mol)
    warnings = []
    for prob in sanitized_problems:
        atom_idx = None
        try:
            atom_idx = prob.GetAtomIdx()
        except Exception:
            pass
        warnings.append({
            "atom_index": atom_idx,
            "type": prob.GetType(),
            "message": prob.Message(),
        })

    return {
        "valid": True,
        "canonical_smiles": canonical,
        "warnings": warnings,
        "errors": [],
    }
