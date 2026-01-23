"""
Excel Exporter

Exports batch results to Excel format with conditional formatting using XlsxWriter.
"""
from io import BytesIO
from typing import List, Dict, Any
import pandas as pd

from .base import BaseExporter, ExportFormat, ExporterFactory


class ExcelExporter(BaseExporter):
    """Export batch results to Excel format with formatting."""

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """
        Export results to Excel format with conditional formatting.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing Excel data
        """
        # Extract relevant fields from results
        rows = []
        successful_count = 0
        total_score = 0
        total_ml_score = 0
        alert_distribution: Dict[str, int] = {}

        for idx, result in enumerate(results):
            # Get validation data
            validation = result.get("validation", {})
            scoring = result.get("scoring", {})
            alerts = result.get("alerts", {})
            status = result.get("status", "error")

            if status == "success":
                successful_count += 1

            # Get scores
            overall_score = validation.get("overall_score", 0)
            ml_score = scoring.get("ml_readiness_score", 0) if scoring else 0

            if status == "success":
                total_score += overall_score
                total_ml_score += ml_score

            # Count alerts
            alerts_count = 0
            if alerts:
                for catalog in ["pains", "brenk", "nih", "glaxo"]:
                    catalog_data = alerts.get(catalog, {})
                    if isinstance(catalog_data, dict):
                        matches = catalog_data.get("matches", [])
                        count = len(matches)
                        alerts_count += count
                        if count > 0:
                            alert_distribution[catalog] = alert_distribution.get(catalog, 0) + count

            # Collect issues for summary
            issues = validation.get("issues", [])
            issues_summary = "; ".join(
                [f"{issue.get('check_name', 'unknown')}: {issue.get('message', '')}"
                 for issue in issues[:3]]  # Limit to first 3 issues
            )

            row = {
                "index": idx + 1,
                "name": result.get("name", ""),
                "input_smiles": result.get("smiles", ""),
                "canonical_smiles": validation.get("canonical_smiles", ""),
                "inchikey": validation.get("inchikey", ""),
                "overall_score": overall_score,
                "ml_readiness_score": ml_score,
                "np_likeness_score": scoring.get("np_likeness_score", 0) if scoring else 0,
                "alerts_count": alerts_count,
                "issues_summary": issues_summary,
                "standardized_smiles": result.get("standardized_smiles", ""),
            }
            rows.append(row)

        # Create DataFrame
        df = pd.DataFrame(rows)

        # Create BytesIO buffer
        bytes_buffer = BytesIO()

        # Write to Excel with xlsxwriter engine
        with pd.ExcelWriter(bytes_buffer, engine="xlsxwriter") as writer:
            # Write main results sheet
            df.to_excel(writer, sheet_name="Results", index=False)

            # Get workbook and worksheet objects
            workbook = writer.book
            worksheet = writer.sheets["Results"]

            # Define formats for conditional coloring
            green_format = workbook.add_format({"bg_color": "#C6EFCE"})
            yellow_format = workbook.add_format({"bg_color": "#FFEB9C"})
            red_format = workbook.add_format({"bg_color": "#FFC7CE"})

            # Apply conditional formatting to score columns (F, G, H - zero-indexed: 5, 6, 7)
            # overall_score column (F, index 5)
            worksheet.conditional_format(
                1,
                5,
                len(df),
                5,
                {
                    "type": "cell",
                    "criteria": ">=",
                    "value": 80,
                    "format": green_format,
                },
            )
            worksheet.conditional_format(
                1,
                5,
                len(df),
                5,
                {
                    "type": "cell",
                    "criteria": "between",
                    "minimum": 50,
                    "maximum": 79,
                    "format": yellow_format,
                },
            )
            worksheet.conditional_format(
                1,
                5,
                len(df),
                5,
                {
                    "type": "cell",
                    "criteria": "<",
                    "value": 50,
                    "format": red_format,
                },
            )

            # ml_readiness_score column (G, index 6)
            worksheet.conditional_format(
                1,
                6,
                len(df),
                6,
                {
                    "type": "cell",
                    "criteria": ">=",
                    "value": 80,
                    "format": green_format,
                },
            )
            worksheet.conditional_format(
                1,
                6,
                len(df),
                6,
                {
                    "type": "cell",
                    "criteria": "between",
                    "minimum": 50,
                    "maximum": 79,
                    "format": yellow_format,
                },
            )
            worksheet.conditional_format(
                1,
                6,
                len(df),
                6,
                {
                    "type": "cell",
                    "criteria": "<",
                    "value": 50,
                    "format": red_format,
                },
            )

            # Freeze first row (header)
            worksheet.freeze_panes(1, 0)

            # Auto-fit column widths (approximate)
            for idx, col in enumerate(df.columns):
                # Calculate max width
                max_width = max(
                    df[col].astype(str).map(len).max(),
                    len(str(col))
                )
                # Add some padding
                worksheet.set_column(idx, idx, min(max_width + 2, 50))

            # Create summary sheet
            total_count = len(results)
            avg_score = total_score / successful_count if successful_count > 0 else 0
            avg_ml_score = total_ml_score / successful_count if successful_count > 0 else 0

            summary_data = {
                "Metric": [
                    "Total Molecules",
                    "Successful Validations",
                    "Failed Validations",
                    "Average Overall Score",
                    "Average ML-Readiness Score",
                ],
                "Value": [
                    total_count,
                    successful_count,
                    total_count - successful_count,
                    f"{avg_score:.2f}",
                    f"{avg_ml_score:.2f}",
                ],
            }

            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name="Summary", index=False)

            # Add alert distribution to summary sheet
            if alert_distribution:
                summary_worksheet = writer.sheets["Summary"]
                row_offset = len(summary_data["Metric"]) + 3

                summary_worksheet.write(row_offset, 0, "Alert Distribution")
                row_offset += 1
                summary_worksheet.write(row_offset, 0, "Catalog")
                summary_worksheet.write(row_offset, 1, "Count")
                row_offset += 1

                for catalog, count in alert_distribution.items():
                    summary_worksheet.write(row_offset, 0, catalog.upper())
                    summary_worksheet.write(row_offset, 1, count)
                    row_offset += 1

        # Seek to beginning
        bytes_buffer.seek(0)

        return bytes_buffer

    @property
    def media_type(self) -> str:
        """MIME type for Excel."""
        return "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"

    @property
    def file_extension(self) -> str:
        """File extension for Excel."""
        return "xlsx"


# Register with factory
ExporterFactory.register(ExportFormat.EXCEL, ExcelExporter)
