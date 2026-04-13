"""
JSON Exporter

Exports batch results to JSON format with structured audit sections.
"""

import json
from datetime import datetime, timezone
from io import BytesIO
from typing import Any, Dict, List

try:
    import orjson

    HAS_ORJSON = True
except ImportError:
    HAS_ORJSON = False

from .audit_columns import extract_nested, get_identity_row
from .base import BaseExporter, ExporterFactory, ExportFormat


class JSONExporter(BaseExporter):
    """Export batch results to JSON format with structured audit sections."""

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """Export results to structured JSON with audit sections per molecule.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing structured JSON data
        """
        structured_results = []
        for idx, result in enumerate(results):
            entry = dict(get_identity_row(idx, result))
            entry["status"] = result.get("status", "error")
            entry.update(extract_nested(result))
            structured_results.append(entry)

        export_data = {
            "metadata": {
                "export_date": datetime.now(timezone.utc).isoformat(),
                "total_count": len(results),
                "format_version": "2.0",
                "tool": "ChemAudit",
            },
            "results": structured_results,
        }

        if HAS_ORJSON:
            json_bytes = orjson.dumps(
                export_data, option=orjson.OPT_INDENT_2 | orjson.OPT_APPEND_NEWLINE
            )
        else:
            json_str = json.dumps(export_data, indent=2, ensure_ascii=False)
            json_bytes = json_str.encode("utf-8")

        bytes_buffer = BytesIO(json_bytes)
        bytes_buffer.seek(0)
        return bytes_buffer

    @property
    def media_type(self) -> str:
        """MIME type for JSON."""
        return "application/json"

    @property
    def file_extension(self) -> str:
        """File extension for JSON."""
        return "json"


# Register with factory
ExporterFactory.register(ExportFormat.JSON, JSONExporter)
