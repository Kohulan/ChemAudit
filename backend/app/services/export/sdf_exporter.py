"""
SDF Exporter

Exports batch results to SDF format using RDKit SDWriter.
"""

import logging
from io import BytesIO, StringIO
from typing import Any, Dict, List

from rdkit import Chem

from .audit_columns import extract_flat_row
from .base import BaseExporter, ExporterFactory, ExportFormat, extract_alert_names

logger = logging.getLogger(__name__)


class SDFExporter(BaseExporter):
    """Export batch results to SDF format with properties."""

    def __init__(self, include_audit: bool = False) -> None:
        """
        Initialise the SDF exporter.

        Args:
            include_audit: When True, attach all audit data as SDF properties
                           (without the section-prefix, e.g. ``Parsability (Pass/Fail)``).
        """
        self._include_audit = include_audit

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """
        Export results to SDF format.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing SDF data
        """
        # Use StringIO for text-based SDF output
        string_buffer = StringIO()

        # Create SDWriter
        writer = Chem.SDWriter(string_buffer)

        skipped_count = 0

        for idx, result in enumerate(results):
            # Prefer standardized SMILES, fall back to canonical, then original input
            validation = result.get("validation") or {}
            smiles = (
                result.get("standardized_smiles")
                or validation.get("canonical_smiles")
                or result.get("smiles")
            )

            if not smiles:
                logger.warning("Skipping molecule at index %d: no valid SMILES", idx)
                skipped_count += 1
                continue

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning("Skipping molecule at index %d: invalid SMILES '%s'", idx, smiles)
                skipped_count += 1
                continue

            mol.SetProp("_Name", result.get("name", f"mol_{idx + 1}"))

            overall_score = validation.get("overall_score", 0)
            mol.SetProp("overall_score", str(overall_score))

            scoring = result.get("scoring") or {}
            ml_readiness = scoring.get("ml_readiness") or {}
            ml_score = ml_readiness.get("score", 0) if isinstance(ml_readiness, dict) else 0
            mol.SetProp("ml_readiness_score", str(ml_score))

            inchikey = validation.get("inchikey", "")
            if inchikey:
                mol.SetProp("inchikey", inchikey)

            alert_names = extract_alert_names(result.get("alerts") or {})
            if alert_names:
                mol.SetProp("alerts", ", ".join(alert_names))

            if self._include_audit:
                flat = extract_flat_row(result)
                for prefixed_key, value in flat.items():
                    # Strip "[Section] " prefix: "[Validation] Parsability" -> "Parsability"
                    short_key = (
                        prefixed_key.split("] ", 1)[1] if "] " in prefixed_key else prefixed_key
                    )
                    str_value = str(value)
                    if str_value not in ("N/A", ""):
                        mol.SetProp(short_key, str_value)

            try:
                writer.write(mol)
            except Exception as e:
                logger.warning("Failed to write molecule at index %d: %s", idx, e)
                skipped_count += 1

        writer.close()

        if skipped_count > 0:
            logger.info("Skipped %d molecules with invalid SMILES", skipped_count)

        # Convert to BytesIO
        bytes_buffer = BytesIO()
        bytes_buffer.write(string_buffer.getvalue().encode("utf-8"))
        bytes_buffer.seek(0)

        return bytes_buffer

    @property
    def media_type(self) -> str:
        """MIME type for SDF."""
        return "chemical/x-mdl-sdfile"

    @property
    def file_extension(self) -> str:
        """File extension for SDF."""
        return "sdf"


# Register with factory
ExporterFactory.register(ExportFormat.SDF, SDFExporter)
