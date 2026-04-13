"""
CSV Exporter

Exports batch results to CSV format with full audit data.
"""

from io import BytesIO, StringIO
from typing import Any, Dict, List

import pandas as pd

from .audit_columns import extract_flat_row, get_flat_headers, get_identity_row
from .base import BaseExporter, ExporterFactory, ExportFormat


class CSVExporter(BaseExporter):
    """Export batch results to CSV format with full audit columns."""

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """Export results to CSV format with all audit data.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing UTF-8 encoded CSV data
        """
        rows = []
        for idx, result in enumerate(results):
            row = get_identity_row(idx, result)
            row.update(extract_flat_row(result))
            rows.append(row)

        columns = ["index", "name", "input_smiles"] + get_flat_headers()
        df = pd.DataFrame(rows, columns=columns)

        string_buffer = StringIO()
        df.to_csv(string_buffer, index=False)

        bytes_buffer = BytesIO()
        bytes_buffer.write(string_buffer.getvalue().encode("utf-8"))
        bytes_buffer.seek(0)
        return bytes_buffer

    @property
    def media_type(self) -> str:
        """MIME type for CSV."""
        return "text/csv"

    @property
    def file_extension(self) -> str:
        """File extension for CSV."""
        return "csv"


# Register with factory
ExporterFactory.register(ExportFormat.CSV, CSVExporter)
