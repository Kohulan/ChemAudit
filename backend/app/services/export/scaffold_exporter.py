"""
Scaffold-Grouped Exporter

Exports batch results as a CSV with Murcko scaffold SMILES and scaffold group assignment.
"""

from io import BytesIO, StringIO
from typing import Any, Dict, List

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from .base import BaseExporter, ExporterFactory, ExportFormat


class ScaffoldExporter(BaseExporter):
    """Export batch results as a CSV grouped by Murcko scaffold.

    Each molecule gets scaffold_smiles (Murcko framework) and scaffold_group (integer ID).
    Acyclic molecules get empty scaffold_smiles and scaffold_group=0.
    """

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """Export scaffold-grouped results as CSV.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing CSV
        """
        # Compute scaffolds for all molecules
        scaffold_smiles_list: list[str] = []
        for result in results:
            smiles = result.get("smiles", "")
            scaffold = ""
            if result.get("status") == "success" and smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    try:
                        core = MurckoScaffold.GetScaffoldForMol(mol)
                        scaffold = Chem.MolToSmiles(core)
                        # Acyclic molecules produce empty scaffold
                        if not scaffold or scaffold == "":
                            scaffold = ""
                    except Exception:
                        scaffold = ""
            scaffold_smiles_list.append(scaffold)

        # Assign scaffold groups: sort unique scaffolds alphabetically, assign IDs
        unique_scaffolds = sorted(set(s for s in scaffold_smiles_list if s))
        scaffold_to_group: Dict[str, int] = {s: i + 1 for i, s in enumerate(unique_scaffolds)}
        # Empty scaffold (acyclic) gets group 0
        scaffold_to_group[""] = 0

        # Build rows
        rows = []
        for idx, result in enumerate(results):
            validation = result.get("validation") or {}
            scoring = result.get("scoring") or {}
            scaffold = scaffold_smiles_list[idx]

            row = {
                "index": idx,
                "smiles": result.get("smiles", ""),
                "name": result.get("name", ""),
                "status": result.get("status", ""),
                "scaffold_smiles": scaffold,
                "scaffold_group": scaffold_to_group.get(scaffold, 0),
                "overall_score": validation.get("overall_score", ""),
                "ml_readiness_score": scoring.get("ml_readiness_score", "") if scoring else "",
                "inchikey": validation.get("inchikey", ""),
            }
            rows.append(row)

        # Write CSV
        csv_buf = StringIO()
        df = pd.DataFrame(rows)
        df.to_csv(csv_buf, index=False)

        bytes_buffer = BytesIO()
        bytes_buffer.write(csv_buf.getvalue().encode("utf-8"))
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
ExporterFactory.register(ExportFormat.SCAFFOLD, ScaffoldExporter)
