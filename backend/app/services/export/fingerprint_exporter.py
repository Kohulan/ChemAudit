"""
Fingerprint Exporter

Exports batch results as fingerprint matrices (Morgan, MACCS, RDKit) in CSV, npy, and npz formats.
All files are packaged into a single zip archive.
"""

from io import BytesIO, StringIO
from typing import Any, Dict, List
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

from .base import BaseExporter, ExporterFactory, ExportFormat


class FingerprintExporter(BaseExporter):
    """Export batch results as fingerprint matrices in multiple formats.

    Produces a zip containing 9 files: 3 FP types (morgan, maccs, rdkit) x 3 formats (csv, npy, npz).
    Molecules with errors or invalid SMILES are skipped.
    """

    FP_CONFIGS = {
        "morgan": {"bits": 2048},
        "maccs": {"bits": 167},
        "rdkit": {"bits": 2048},
    }

    def _compute_fingerprint(
        self, mol: Chem.Mol, fp_type: str
    ) -> list[int] | None:
        """Compute fingerprint for a molecule.

        Args:
            mol: RDKit Mol object
            fp_type: One of 'morgan', 'maccs', 'rdkit'

        Returns:
            List of bit values or None if computation fails
        """
        try:
            if fp_type == "morgan":
                gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
                fp = gen.GetFingerprint(mol)
                return list(fp.ToList())
            elif fp_type == "maccs":
                fp = MACCSkeys.GenMACCSKeys(mol)
                return list(fp.ToList())
            elif fp_type == "rdkit":
                gen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
                fp = gen.GetFingerprint(mol)
                return list(fp.ToList())
        except Exception:
            return None
        return None

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """Export fingerprint matrices as a zip archive.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing zip with 9 files
        """
        # Filter valid molecules
        valid_entries: list[tuple[str, Chem.Mol]] = []
        for result in results:
            if result.get("status") != "success":
                continue
            smiles = result.get("smiles", "")
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            name = result.get("name", smiles)
            valid_entries.append((name, mol))

        zip_buffer = BytesIO()
        with ZipFile(zip_buffer, "w", ZIP_DEFLATED) as zf:
            for fp_type, config in self.FP_CONFIGS.items():
                nbits = config["bits"]
                names: list[str] = []
                fp_matrix: list[list[int]] = []

                for name, mol in valid_entries:
                    bits = self._compute_fingerprint(mol, fp_type)
                    if bits is not None:
                        names.append(name)
                        fp_matrix.append(bits)

                if not fp_matrix:
                    # Write empty files for this FP type
                    arr = np.zeros((0, nbits), dtype=np.int8)
                else:
                    arr = np.array(fp_matrix, dtype=np.int8)

                # CSV with name column + bit columns
                csv_buf = StringIO()
                col_names = ["name"] + [f"bit_{i}" for i in range(arr.shape[1] if arr.shape[0] > 0 else nbits)]
                df = pd.DataFrame(
                    arr, columns=[f"bit_{i}" for i in range(arr.shape[1] if arr.shape[0] > 0 else nbits)]
                )
                df.insert(0, "name", names if names else [])
                df.to_csv(csv_buf, index=False)
                zf.writestr(f"{fp_type}.csv", csv_buf.getvalue())

                # NPY
                npy_buf = BytesIO()
                np.save(npy_buf, arr)
                npy_buf.seek(0)
                zf.writestr(f"{fp_type}.npy", npy_buf.read())

                # NPZ (compressed)
                npz_buf = BytesIO()
                np.savez_compressed(npz_buf, fingerprints=arr)
                npz_buf.seek(0)
                zf.writestr(f"{fp_type}.npz", npz_buf.read())

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
ExporterFactory.register(ExportFormat.FINGERPRINT, FingerprintExporter)
