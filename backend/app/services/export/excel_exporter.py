"""
Excel Exporter

Exports batch results to Excel format with conditional formatting using XlsxWriter.
Optionally embeds 2D chemical structure images rendered via RDKit.
"""

import logging
from io import BytesIO
from typing import Any, Dict, List, Optional

import pandas as pd

from .base import (
    BaseExporter,
    ExporterFactory,
    ExportFormat,
    count_alerts,
    count_alerts_by_catalog,
)

logger = logging.getLogger(__name__)

# Image dimensions (pixels) and corresponding Excel cell sizing
STRUCTURE_IMG_WIDTH = 150
STRUCTURE_IMG_HEIGHT = 150
# Excel row height in points (1 pt ≈ 1.33 px). Add padding for cell borders.
STRUCTURE_ROW_HEIGHT = 118
# Excel column width in character units (1 char ≈ 7.5 px at default font)
STRUCTURE_COL_WIDTH = 22


class ExcelExporter(BaseExporter):
    """Export batch results to Excel format with formatting."""

    def __init__(self, include_images: bool = False):
        self._include_images = include_images

    @staticmethod
    def _mol_to_png_bytes(smiles: str) -> Optional[BytesIO]:
        """Render a SMILES string as a PNG image and return as BytesIO buffer."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            img = Draw.MolToImage(mol, size=(STRUCTURE_IMG_WIDTH, STRUCTURE_IMG_HEIGHT))
            buf = BytesIO()
            img.save(buf, format="PNG")
            buf.seek(0)
            return buf
        except Exception:
            logger.debug("Failed to render structure for SMILES: %s", smiles, exc_info=True)
            return None

    def _embed_structure_images(self, worksheet: Any, results: List[Dict[str, Any]]) -> None:
        """Render and embed 2D structure PNGs into the worksheet's structure column."""
        img_col = 1  # Column B (after index)

        for row_idx, result in enumerate(results):
            excel_row = row_idx + 1  # +1 for header row
            worksheet.set_row(excel_row, STRUCTURE_ROW_HEIGHT)

            # Prefer canonical SMILES, fall back to input
            validation = result.get("validation") or {}
            smiles = validation.get("canonical_smiles") or result.get("smiles", "")
            if not smiles:
                continue

            img_buf = self._mol_to_png_bytes(smiles)
            if img_buf is None:
                continue

            worksheet.insert_image(
                excel_row,
                img_col,
                f"structure_{row_idx}.png",
                {
                    "image_data": img_buf,
                    "x_offset": 4,
                    "y_offset": 4,
                    "positioning": 1,  # Move and size with cells
                },
            )

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
        qed_scores = []
        sa_scores = []
        lipinski_passes = 0
        lipinski_total = 0
        safety_passes = 0
        safety_total = 0
        alert_distribution: Dict[str, int] = {}

        for idx, result in enumerate(results):
            validation = result.get("validation") or {}
            scoring = result.get("scoring") or {}
            alerts = result.get("alerts") or {}
            status = result.get("status", "error")

            if status == "success":
                successful_count += 1

            overall_score = validation.get("overall_score", 0)
            ml_score = scoring.get("ml_readiness_score", 0) if scoring else 0

            qed_score = None
            sa_score = None
            lipinski_passed = None
            safety_passed = None

            if scoring:
                druglikeness = scoring.get("druglikeness") or {}
                if druglikeness and "error" not in druglikeness:
                    qed_score = druglikeness.get("qed_score")
                    if qed_score is not None:
                        qed_scores.append(qed_score)
                    lipinski_passed = druglikeness.get("lipinski_passed")
                    if lipinski_passed is not None:
                        lipinski_total += 1
                        if lipinski_passed:
                            lipinski_passes += 1

                admet = scoring.get("admet") or {}
                if admet and "error" not in admet:
                    sa_score = admet.get("sa_score")
                    if sa_score is not None:
                        sa_scores.append(sa_score)

                safety_filters = scoring.get("safety_filters") or {}
                if safety_filters and "error" not in safety_filters:
                    safety_passed = safety_filters.get("all_passed")
                    if safety_passed is not None:
                        safety_total += 1
                        if safety_passed:
                            safety_passes += 1

            if status == "success":
                total_score += overall_score
                total_ml_score += ml_score

            # Count alerts and update distribution
            alerts_count = count_alerts(alerts)
            for catalog, cnt in count_alerts_by_catalog(alerts).items():
                alert_distribution[catalog] = alert_distribution.get(catalog, 0) + cnt

            # Collect issues for summary (limit to first 3)
            issues = validation.get("issues", [])
            issues_summary = "; ".join(
                f"{issue.get('check_name', 'unknown')}: {issue.get('message', '')}"
                for issue in issues[:3]
            )

            row = {"index": idx + 1}
            if self._include_images:
                row["structure"] = ""
            row.update(
                {
                    "name": result.get("name", ""),
                    "input_smiles": result.get("smiles", ""),
                    "canonical_smiles": validation.get("canonical_smiles", ""),
                    "inchikey": validation.get("inchikey", ""),
                    "overall_score": overall_score,
                    "ml_readiness_score": ml_score,
                    "qed_score": qed_score if qed_score is not None else "",
                    "sa_score": sa_score if sa_score is not None else "",
                    "np_likeness_score": scoring.get("np_likeness_score", 0)
                    if scoring
                    else 0,
                    "alerts_count": alerts_count,
                    "issues_summary": issues_summary,
                    "standardized_smiles": result.get("standardized_smiles", ""),
                }
            )
            rows.append(row)

        df = pd.DataFrame(rows)
        bytes_buffer = BytesIO()

        with pd.ExcelWriter(bytes_buffer, engine="xlsxwriter") as writer:
            df.to_excel(writer, sheet_name="Results", index=False)
            workbook = writer.book
            worksheet = writer.sheets["Results"]

            # Column offset: 1 extra column when structure images are present
            img_offset = 1 if self._include_images else 0

            # Define formats for conditional coloring
            green_format = workbook.add_format({"bg_color": "#C6EFCE"})
            yellow_format = workbook.add_format({"bg_color": "#FFEB9C"})
            red_format = workbook.add_format({"bg_color": "#FFC7CE"})

            # Score columns are at indices 5 and 6 without images, shifted by 1 with images
            score_col_indices = [5 + img_offset, 6 + img_offset]
            for col_idx in score_col_indices:
                worksheet.conditional_format(
                    1,
                    col_idx,
                    len(df),
                    col_idx,
                    {
                        "type": "cell",
                        "criteria": ">=",
                        "value": 80,
                        "format": green_format,
                    },
                )
                worksheet.conditional_format(
                    1,
                    col_idx,
                    len(df),
                    col_idx,
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
                    col_idx,
                    len(df),
                    col_idx,
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
                if self._include_images and col == "structure":
                    # Fixed width for image column
                    worksheet.set_column(idx, idx, STRUCTURE_COL_WIDTH)
                    continue
                # Calculate max width
                max_width = max(df[col].astype(str).map(len).max(), len(str(col)))
                # Add some padding
                worksheet.set_column(idx, idx, min(max_width + 2, 50))

            if self._include_images:
                self._embed_structure_images(worksheet, results)

            # Create summary sheet
            total_count = len(results)
            avg_score = total_score / successful_count if successful_count > 0 else 0
            avg_ml_score = (
                total_ml_score / successful_count if successful_count > 0 else 0
            )
            avg_qed = sum(qed_scores) / len(qed_scores) if qed_scores else 0
            avg_sa = sum(sa_scores) / len(sa_scores) if sa_scores else 0
            lipinski_pass_rate = (
                (lipinski_passes / lipinski_total) * 100 if lipinski_total > 0 else 0
            )
            safety_pass_rate = (
                (safety_passes / safety_total) * 100 if safety_total > 0 else 0
            )

            summary_data = {
                "Metric": [
                    "Total Molecules",
                    "Successful Validations",
                    "Failed Validations",
                    "Average Overall Score",
                    "Average ML-Readiness Score",
                    "Average QED Score",
                    "Average SA Score",
                    "Lipinski Pass Rate (%)",
                    "Safety Pass Rate (%)",
                ],
                "Value": [
                    total_count,
                    successful_count,
                    total_count - successful_count,
                    f"{avg_score:.2f}",
                    f"{avg_ml_score:.2f}",
                    f"{avg_qed:.2f}" if qed_scores else "-",
                    f"{avg_sa:.1f}" if sa_scores else "-",
                    f"{lipinski_pass_rate:.1f}" if lipinski_total > 0 else "-",
                    f"{safety_pass_rate:.1f}" if safety_total > 0 else "-",
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
