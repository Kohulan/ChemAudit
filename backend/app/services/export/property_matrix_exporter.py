"""
Property Matrix Exporter

Exports all computed molecular properties as both a flat CSV and a multi-sheet Excel file.
The two files are packaged into a single zip archive.
"""

from io import BytesIO, StringIO
from typing import Any, Dict, List
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd

from .base import (
    BaseExporter,
    ExporterFactory,
    ExportFormat,
    count_alerts,
    count_alerts_by_catalog,
    extract_alert_names,
)


class PropertyMatrixExporter(BaseExporter):
    """Export all computed properties as a flat CSV + multi-sheet Excel zip.

    The zip contains:
    - properties_flat.csv: All properties as columns
    - properties.xlsx: 4-sheet Excel (Descriptors, Scores, Alerts, Properties)
    """

    def _extract_properties(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Extract all property fields from a batch result.

        Args:
            result: Single batch result dictionary

        Returns:
            Dictionary of property name -> value
        """
        validation = result.get("validation") or {}
        scoring = result.get("scoring") or {}
        alerts = result.get("alerts") or {}
        properties = result.get("properties") or {}

        # Druglikeness sub-dict
        druglikeness = scoring.get("druglikeness") or {}
        admet = scoring.get("admet") or {}

        # Alert analysis
        alert_count = count_alerts(alerts)
        alert_names_list = extract_alert_names(alerts)
        alert_catalogs = count_alerts_by_catalog(alerts)

        return {
            # Identifiers
            "smiles": result.get("smiles", ""),
            "name": result.get("name", ""),
            "inchikey": validation.get("inchikey", ""),
            "canonical_smiles": validation.get("canonical_smiles", ""),
            # Descriptors
            "mw": properties.get("mw") or druglikeness.get("mw", ""),
            "logp": properties.get("logp") or druglikeness.get("logp", ""),
            "tpsa": properties.get("tpsa") or druglikeness.get("tpsa", ""),
            "hbd": properties.get("hbd") or druglikeness.get("hbd", ""),
            "hba": properties.get("hba") or druglikeness.get("hba", ""),
            "rotatable_bonds": properties.get("rotatable_bonds")
            or druglikeness.get("rotatable_bonds", ""),
            "aromatic_rings": properties.get("aromatic_rings")
            or druglikeness.get("aromatic_rings", ""),
            "fsp3": properties.get("fsp3") or druglikeness.get("fsp3", ""),
            "heavy_atom_count": properties.get("heavy_atom_count", ""),
            "ring_count": properties.get("ring_count", ""),
            # Scores
            "overall_score": validation.get("overall_score", ""),
            "ml_readiness_score": scoring.get("ml_readiness_score", ""),
            "qed_score": druglikeness.get("qed_score", ""),
            "sa_score": admet.get("sa_score", ""),
            "np_likeness_score": scoring.get("np_likeness_score", ""),
            "lipinski_passed": druglikeness.get("lipinski_passed", ""),
            # Alerts
            "alert_count": alert_count,
            "alert_names": "; ".join(alert_names_list),
            "pains_count": alert_catalogs.get("pains", 0),
            "brenk_count": alert_catalogs.get("brenk", 0),
            "nih_count": alert_catalogs.get("nih", 0),
            "glaxo_count": alert_catalogs.get("glaxo", 0),
            # Status
            "status": result.get("status", ""),
        }

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """Export property matrix as a zip with flat CSV and multi-sheet Excel.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing zip with 2 files
        """
        # Extract all properties
        rows = [self._extract_properties(r) for r in results]

        if not rows:
            rows = []

        all_df = pd.DataFrame(rows)

        zip_buffer = BytesIO()
        with ZipFile(zip_buffer, "w", ZIP_DEFLATED) as zf:
            # Flat CSV
            csv_buf = StringIO()
            all_df.to_csv(csv_buf, index=False)
            zf.writestr("properties_flat.csv", csv_buf.getvalue())

            # Multi-sheet Excel
            excel_buf = BytesIO()
            with pd.ExcelWriter(excel_buf, engine="xlsxwriter") as writer:
                # Descriptors sheet
                desc_cols = [
                    "smiles", "name", "inchikey", "mw", "logp", "tpsa",
                    "hbd", "hba", "rotatable_bonds", "aromatic_rings",
                    "fsp3", "heavy_atom_count", "ring_count",
                ]
                desc_df = all_df[[c for c in desc_cols if c in all_df.columns]]
                desc_df.to_excel(writer, sheet_name="Descriptors", index=False)

                # Scores sheet
                score_cols = [
                    "smiles", "name", "overall_score", "ml_readiness_score",
                    "qed_score", "sa_score", "np_likeness_score", "lipinski_passed",
                ]
                score_df = all_df[[c for c in score_cols if c in all_df.columns]]
                score_df.to_excel(writer, sheet_name="Scores", index=False)

                # Alerts sheet
                alert_cols = [
                    "smiles", "name", "alert_count", "alert_names",
                    "pains_count", "brenk_count", "nih_count", "glaxo_count",
                ]
                alert_df = all_df[[c for c in alert_cols if c in all_df.columns]]
                alert_df.to_excel(writer, sheet_name="Alerts", index=False)

                # Properties sheet (all remaining)
                all_df.to_excel(writer, sheet_name="Properties", index=False)

            # CRITICAL: seek(0) after ExcelWriter context manager closes
            excel_buf.seek(0)
            zf.writestr("properties.xlsx", excel_buf.read())

        zip_buffer.seek(0)
        return zip_buffer

    @property
    def media_type(self) -> str:
        """MIME type for zip."""
        return "application/zip"

    @property
    def file_extension(self) -> str:
        """File extension for zip."""
        return "zip"


# Register with factory
ExporterFactory.register(ExportFormat.PROPERTY_MATRIX, PropertyMatrixExporter)
