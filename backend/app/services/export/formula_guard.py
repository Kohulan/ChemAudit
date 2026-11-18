"""
Spreadsheet formula injection guard for exported user text.

Molecule names, SMILES and dataset properties come from uploaded files. A cell
that starts with = + - @ (or tab / carriage return, or an array formula "{=")
is run as a formula when the CSV or XLSX is opened in Excel, LibreOffice or
Google Sheets. Such text is prefixed with a single quote (OWASP CSV injection
guidance). Values that parse as numbers are left alone, so -1.5 stays -1.5.
"""

from typing import Any

import pandas as pd

_FORMULA_STARTS = ("=", "+", "-", "@", "\t", "\r", "{=")

# Options for pd.ExcelWriter(engine="xlsxwriter"): never turn strings into
# formulas or hyperlinks.
XLSX_ENGINE_KWARGS = {"options": {"strings_to_formulas": False, "strings_to_urls": False}}


def neutralize_cell(value: Any) -> Any:
    """Return ``value`` with a leading quote if it is text a spreadsheet would run."""
    if not isinstance(value, str) or not value.startswith(_FORMULA_STARTS):
        return value
    try:
        float(value)
        return value
    except ValueError:
        return "'" + value


def neutralize_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Copy of ``df`` with every text cell and column name neutralised."""
    out = df.copy()
    text_cols = out.select_dtypes(include=["object", "string"]).columns
    for col in text_cols:
        out[col] = out[col].map(neutralize_cell)
    out.columns = [neutralize_cell(c) for c in out.columns]
    return out


def neutralize_row(row: dict) -> dict:
    """Copy of a csv.DictWriter row with every value neutralised (keys unchanged)."""
    return {k: neutralize_cell(v) for k, v in row.items()}
