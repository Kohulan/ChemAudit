"""
Excel Exporter

Exports batch results to Excel format with conditional formatting using XlsxWriter.
Supports two layout modes:
  - "single" (default): one "Results" sheet with all audit columns prefixed by section.
  - "multi": one sheet per audit section plus a "Summary" sheet.

Optionally embeds 2D chemical structure images rendered via RDKit.
"""

import logging
from collections import OrderedDict
from io import BytesIO
from typing import Any, Dict, List, Optional

import pandas as pd

from .audit_columns import (
    AUDIT_SECTIONS,
    extract_by_section,
    extract_flat_row,
    get_identity_row,
)
from .base import (
    BaseExporter,
    ExporterFactory,
    ExportFormat,
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

# Headers whose columns receive green/yellow/red conditional formatting.
# Match is performed as a substring of the column header.
_SCORE_HEADER_SUBSTRINGS = ("Score (0-100)", "QED (0-1)")


def _is_score_column(header: str) -> bool:
    """Return True when *header* should receive score conditional formatting."""
    return any(sub in header for sub in _SCORE_HEADER_SUBSTRINGS)


class ExcelExporter(BaseExporter):
    """Export batch results to Excel format with formatting.

    Args:
        include_images: Embed 2D structure PNGs in the first data sheet.
        sheet_layout: ``"single"`` (default) creates one Results sheet with all
            prefixed audit columns; ``"multi"`` creates one sheet per audit section.
    """

    def __init__(self, include_images: bool = False, sheet_layout: str = "single"):
        self._include_images = include_images
        self._sheet_layout = sheet_layout

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

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

    def _embed_structure_images(
        self, worksheet: Any, results: List[Dict[str, Any]], *, img_col: int
    ) -> None:
        """Render and embed 2D structure PNGs into the worksheet's structure column.

        Args:
            worksheet: XlsxWriter Worksheet instance.
            results: List of batch result dictionaries.
            img_col: 0-based column index where the "structure" placeholder lives.
        """

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

    @staticmethod
    def _apply_score_formatting(
        workbook: Any,
        worksheet: Any,
        headers: List[str],
        n_data_rows: int,
    ) -> None:
        """Apply green/yellow/red conditional formatting to score columns.

        Args:
            workbook: XlsxWriter Workbook instance.
            worksheet: XlsxWriter Worksheet instance.
            headers: Ordered list of column header strings (matches DataFrame columns).
            n_data_rows: Number of data rows (excluding header).
        """
        green_format = workbook.add_format({"bg_color": "#C6EFCE"})
        yellow_format = workbook.add_format({"bg_color": "#FFEB9C"})
        red_format = workbook.add_format({"bg_color": "#FFC7CE"})

        for col_idx, header in enumerate(headers):
            if not _is_score_column(str(header)):
                continue
            worksheet.conditional_format(
                1,
                col_idx,
                n_data_rows,
                col_idx,
                {"type": "cell", "criteria": ">=", "value": 80, "format": green_format},
            )
            worksheet.conditional_format(
                1,
                col_idx,
                n_data_rows,
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
                n_data_rows,
                col_idx,
                {"type": "cell", "criteria": "<", "value": 50, "format": red_format},
            )

    @staticmethod
    def _autofit_columns(worksheet: Any, df: "pd.DataFrame", include_images: bool) -> None:
        """Set approximate column widths based on content.

        Args:
            worksheet: XlsxWriter Worksheet instance.
            df: DataFrame whose columns are already written to *worksheet*.
            include_images: When True, the "structure" column gets a fixed width.
        """
        for idx, col in enumerate(df.columns):
            if include_images and col == "structure":
                worksheet.set_column(idx, idx, STRUCTURE_COL_WIDTH)
                continue
            max_width = max(df[col].astype(str).map(len).max(), len(str(col)))
            worksheet.set_column(idx, idx, min(max_width + 2, 50))

    # ------------------------------------------------------------------
    # Summary stats collector
    # ------------------------------------------------------------------

    @staticmethod
    def _collect_summary_stats(
        results: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Walk *results* once and return summary statistics.

        Returns a dict with keys used by :meth:`_write_summary_sheet`.
        """
        successful_count = 0
        total_score = 0
        total_ml_score = 0
        qed_scores: List[float] = []
        sa_scores: List[float] = []
        lipinski_passes = 0
        lipinski_total = 0
        safety_passes = 0
        safety_total = 0
        alert_distribution: Dict[str, int] = {}

        for result in results:
            validation = result.get("validation") or {}
            scoring = result.get("scoring") or {}
            alerts = result.get("alerts") or {}
            status = result.get("status", "error")

            overall_score = validation.get("overall_score", 0)
            ml_readiness = scoring.get("ml_readiness") or {}
            ml_score = ml_readiness.get("score", 0) if isinstance(ml_readiness, dict) else 0

            if status == "success":
                successful_count += 1
                total_score += overall_score
                total_ml_score += ml_score

            if scoring:
                druglikeness = scoring.get("druglikeness") or {}
                if druglikeness and "error" not in druglikeness:
                    qed_score = druglikeness.get("qed_score")
                    if qed_score is not None:
                        qed_scores.append(qed_score)
                    lipinski_passed = druglikeness.get("lipinski_passed")
                    if lipinski_passed is None:
                        lipinski_sub = druglikeness.get("lipinski") or {}
                        lipinski_passed = lipinski_sub.get("passed")
                    if lipinski_passed is not None:
                        lipinski_total += 1
                        if lipinski_passed:
                            lipinski_passes += 1

                admet = scoring.get("admet") or {}
                if admet and "error" not in admet:
                    sa_sub = admet.get("synthetic_accessibility") or {}
                    sa_score = (
                        sa_sub.get("score") if isinstance(sa_sub, dict) else admet.get("sa_score")
                    )
                    if sa_score is not None:
                        sa_scores.append(sa_score)

                safety_filters = scoring.get("safety_filters") or {}
                if safety_filters and "error" not in safety_filters:
                    safety_passed = safety_filters.get("all_passed")
                    if safety_passed is not None:
                        safety_total += 1
                        if safety_passed:
                            safety_passes += 1

            for catalog, cnt in count_alerts_by_catalog(alerts).items():
                alert_distribution[catalog] = alert_distribution.get(catalog, 0) + cnt

        return {
            "total_count": len(results),
            "successful_count": successful_count,
            "total_score": total_score,
            "total_ml_score": total_ml_score,
            "qed_scores": qed_scores,
            "sa_scores": sa_scores,
            "lipinski_passes": lipinski_passes,
            "lipinski_total": lipinski_total,
            "safety_passes": safety_passes,
            "safety_total": safety_total,
            "alert_distribution": alert_distribution,
        }

    @staticmethod
    def _write_summary_sheet(writer: Any, stats: Dict[str, Any]) -> None:
        """Write the Summary sheet using pre-collected *stats*.

        Args:
            writer: pandas ExcelWriter (xlsxwriter engine).
            stats: Dict returned by :meth:`_collect_summary_stats`.
        """
        total_count = stats["total_count"]
        successful_count = stats["successful_count"]
        total_score = stats["total_score"]
        total_ml_score = stats["total_ml_score"]
        qed_scores = stats["qed_scores"]
        sa_scores = stats["sa_scores"]
        lipinski_passes = stats["lipinski_passes"]
        lipinski_total = stats["lipinski_total"]
        safety_passes = stats["safety_passes"]
        safety_total = stats["safety_total"]
        alert_distribution = stats["alert_distribution"]

        avg_score = total_score / successful_count if successful_count > 0 else 0
        avg_ml_score = total_ml_score / successful_count if successful_count > 0 else 0
        avg_qed = sum(qed_scores) / len(qed_scores) if qed_scores else 0
        avg_sa = sum(sa_scores) / len(sa_scores) if sa_scores else 0
        lipinski_pass_rate = (lipinski_passes / lipinski_total) * 100 if lipinski_total > 0 else 0
        safety_pass_rate = (safety_passes / safety_total) * 100 if safety_total > 0 else 0

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

    # ------------------------------------------------------------------
    # Layout builders
    # ------------------------------------------------------------------

    def _build_single_sheet(
        self,
        writer: Any,
        results: List[Dict[str, Any]],
    ) -> None:
        """Write a single "Results" sheet with all audit columns prefixed by section."""
        rows = []
        for idx, result in enumerate(results):
            identity = get_identity_row(idx, result)
            if self._include_images:
                identity["structure"] = ""  # placeholder column for images
            flat = extract_flat_row(result)
            row: Dict[str, Any] = dict(identity)
            row.update(flat)
            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_excel(writer, sheet_name="Results", index=False)

        workbook = writer.book
        worksheet = writer.sheets["Results"]

        # Conditional formatting on score columns
        self._apply_score_formatting(workbook, worksheet, list(df.columns), len(df))

        # Freeze first row (header)
        worksheet.freeze_panes(1, 0)

        # Auto-fit column widths
        self._autofit_columns(worksheet, df, self._include_images)

        # Embed structure images if requested
        if self._include_images:
            img_col = list(df.columns).index("structure")
            self._embed_structure_images(worksheet, results, img_col=img_col)

    def _build_multi_sheet(
        self,
        writer: Any,
        results: List[Dict[str, Any]],
    ) -> None:
        """Write one sheet per audit section, each with identity + section columns."""
        # Pre-compute per-result: identity rows and section extractions (once each)
        identities = []
        all_by_section = []
        for idx, result in enumerate(results):
            identities.append(get_identity_row(idx, result))
            all_by_section.append(extract_by_section(result))

        first_sheet = True
        for section in AUDIT_SECTIONS:
            rows = []
            for idx, result in enumerate(results):
                identity = OrderedDict(identities[idx])
                if first_sheet and self._include_images:
                    identity["structure"] = ""

                section_data = all_by_section[idx].get(section.name, {})

                row: Dict[str, Any] = dict(identity)
                row.update(section_data)
                rows.append(row)

            df = pd.DataFrame(rows)
            sheet_name = section.name  # e.g. "Validation", "Deep Validation", …
            df.to_excel(writer, sheet_name=sheet_name, index=False)

            workbook = writer.book
            worksheet = writer.sheets[sheet_name]

            # Conditional formatting (no section prefix in multi-sheet headers)
            self._apply_score_formatting(workbook, worksheet, list(df.columns), len(df))

            # Freeze header row
            worksheet.freeze_panes(1, 0)

            # Auto-fit columns
            self._autofit_columns(worksheet, df, self._include_images and first_sheet)

            # Embed structure images on the first sheet only
            if first_sheet and self._include_images:
                img_col = list(df.columns).index("structure")
                self._embed_structure_images(worksheet, results, img_col=img_col)

            first_sheet = False

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """Export results to Excel format with conditional formatting.

        Args:
            results: List of batch result dictionaries.

        Returns:
            BytesIO buffer containing Excel data.
        """
        stats = self._collect_summary_stats(results)
        bytes_buffer = BytesIO()

        with pd.ExcelWriter(bytes_buffer, engine="xlsxwriter") as writer:
            if self._sheet_layout == "multi":
                self._build_multi_sheet(writer, results)
            else:
                self._build_single_sheet(writer, results)

            self._write_summary_sheet(writer, stats)

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
