"""
Deduplication Exporter

Exports batch results with deduplication grouping: a summary CSV (one row per group)
and a full annotated CSV (every molecule with group assignment).
"""

from collections import defaultdict
from io import BytesIO, StringIO
from typing import Any, Dict, List
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd
from rdkit import Chem

from .base import BaseExporter, ExporterFactory, ExportFormat


class DedupExporter(BaseExporter):
    """Export deduplicated batch results as a zip with summary and annotated CSVs.

    Summary CSV: one row per unique group (representative, group size, member indices).
    Annotated CSV: every molecule with group_id, is_representative, and dedup_level columns.
    """

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """Export deduplication grouping as a zip archive.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing zip with 2 CSV files
        """
        # Build canonical SMILES groups
        groups: Dict[str, list[int]] = defaultdict(list)
        canonical_map: Dict[int, str] = {}

        for idx, result in enumerate(results):
            smiles = result.get("smiles", "")
            canonical = None

            if result.get("status") == "success":
                # Try to get canonical SMILES from validation result
                validation = result.get("validation") or {}
                canonical = validation.get("canonical_smiles")

            if not canonical and smiles:
                # Fall back to RDKit canonicalization
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    canonical = Chem.MolToSmiles(mol)

            if not canonical:
                # Unique group per molecule with unparseable SMILES
                canonical = f"__unparseable_{idx}"

            canonical_map[idx] = canonical
            groups[canonical].append(idx)

        # Assign group IDs (sorted by first appearance)
        group_id_map: Dict[str, int] = {}
        for group_idx, (canonical, indices) in enumerate(groups.items()):
            group_id_map[canonical] = group_idx + 1

        # Build summary rows
        summary_rows = []
        for canonical, indices in groups.items():
            representative_idx = min(indices)
            representative_result = results[representative_idx]
            representative_smiles = representative_result.get("smiles", "")

            summary_rows.append({
                "group_id": group_id_map[canonical],
                "representative_smiles": representative_smiles,
                "canonical_smiles": canonical if not canonical.startswith("__unparseable") else "",
                "group_size": len(indices),
                "duplicate_indices": ",".join(str(i) for i in indices),
                "dedup_level": "exact",
            })

        # Build annotated rows
        annotated_rows = []
        for idx, result in enumerate(results):
            canonical = canonical_map[idx]
            gid = group_id_map[canonical]
            representative_idx = min(groups[canonical])

            row = {
                "index": idx,
                "smiles": result.get("smiles", ""),
                "name": result.get("name", ""),
                "status": result.get("status", ""),
                "is_representative": idx == representative_idx,
                "group_id": gid,
                "dedup_level": "exact",
            }

            # Include validation score if available
            validation = result.get("validation") or {}
            row["overall_score"] = validation.get("overall_score", "")
            row["inchikey"] = validation.get("inchikey", "")

            annotated_rows.append(row)

        # Write to zip
        zip_buffer = BytesIO()
        with ZipFile(zip_buffer, "w", ZIP_DEFLATED) as zf:
            # Summary CSV
            summary_buf = StringIO()
            summary_df = pd.DataFrame(summary_rows)
            summary_df.to_csv(summary_buf, index=False)
            zf.writestr("dedup_summary.csv", summary_buf.getvalue())

            # Annotated CSV
            annotated_buf = StringIO()
            annotated_df = pd.DataFrame(annotated_rows)
            annotated_df.to_csv(annotated_buf, index=False)
            zf.writestr("dedup_annotated.csv", annotated_buf.getvalue())

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
ExporterFactory.register(ExportFormat.DEDUP, DedupExporter)
