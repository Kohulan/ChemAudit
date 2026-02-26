"""
Tests for Advanced Exporters

Tests for fingerprint, deduplication, scaffold, and property matrix exporters.
"""

import zipfile
from io import BytesIO

import numpy as np
import pandas as pd

from app.services.export.base import ExporterFactory, ExportFormat
from app.services.export.dedup_exporter import DedupExporter
from app.services.export.fingerprint_exporter import FingerprintExporter
from app.services.export.property_matrix_exporter import PropertyMatrixExporter
from app.services.export.scaffold_exporter import ScaffoldExporter

# Test data: 3 valid + 1 error molecule
SAMPLE_RESULTS = [
    {
        "smiles": "CCO",
        "name": "Ethanol",
        "index": 0,
        "status": "success",
        "validation": {
            "canonical_smiles": "CCO",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "overall_score": 95,
            "issues": [],
        },
        "alerts": {"pains": {"matches": []}, "brenk": {"matches": []}},
        "scoring": {"ml_readiness_score": 88, "np_likeness_score": -0.5},
        "standardized_smiles": "CCO",
        "properties": {"mw": 46.07, "logp": -0.31, "tpsa": 20.23, "hbd": 1, "hba": 1},
    },
    {
        "smiles": "c1ccccc1",
        "name": "Benzene",
        "index": 1,
        "status": "success",
        "validation": {
            "canonical_smiles": "c1ccccc1",
            "inchikey": "UHOVQNZJYSORNB-UHFFFAOYSA-N",
            "overall_score": 85,
            "issues": [],
        },
        "alerts": {
            "pains": {"matches": [{"pattern_name": "anil_di_alk_D(258)", "smarts": "[Br]"}]},
            "brenk": {"matches": []},
        },
        "scoring": {"ml_readiness_score": 75, "np_likeness_score": 1.2},
        "standardized_smiles": "c1ccccc1",
        "properties": {"mw": 78.11, "logp": 1.56, "tpsa": 0.0, "hbd": 0, "hba": 0},
    },
    {
        "smiles": "CCO",
        "name": "Ethanol-dup",
        "index": 2,
        "status": "success",
        "validation": {
            "canonical_smiles": "CCO",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "overall_score": 90,
            "issues": [],
        },
        "alerts": {"pains": {"matches": []}, "brenk": {"matches": []}},
        "scoring": {"ml_readiness_score": 80, "np_likeness_score": -0.3},
        "standardized_smiles": "CCO",
        "properties": {"mw": 46.07, "logp": -0.31, "tpsa": 20.23, "hbd": 1, "hba": 1},
    },
    {
        "smiles": "INVALID",
        "name": "Bad Molecule",
        "index": 3,
        "status": "error",
        "error": "Invalid SMILES",
        "validation": {},
        "alerts": {},
        "scoring": {},
    },
]


class TestFingerprintExporter:
    """Test FingerprintExporter."""

    def test_fingerprint_exporter_produces_zip(self):
        """Feed 3 valid + 1 error molecule, verify zip contains 9 files."""
        exporter = FingerprintExporter()
        buf = exporter.export(SAMPLE_RESULTS)

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            names = zf.namelist()
            assert len(names) == 9
            expected = [
                "morgan.csv", "morgan.npy", "morgan.npz",
                "maccs.csv", "maccs.npy", "maccs.npz",
                "rdkit.csv", "rdkit.npy", "rdkit.npz",
            ]
            for name in expected:
                assert name in names, f"Missing {name} in zip"

    def test_fingerprint_csv_columns(self):
        """Verify CSV has name + bit columns and correct row count."""
        exporter = FingerprintExporter()
        buf = exporter.export(SAMPLE_RESULTS)

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            morgan_csv = zf.read("morgan.csv").decode("utf-8")
            df = pd.read_csv(BytesIO(morgan_csv.encode("utf-8")))
            # 3 valid molecules (INVALID skipped)
            assert len(df) == 3
            assert "name" in df.columns
            # Morgan FP has 2048 bits
            assert len(df.columns) == 2049  # name + 2048 bits

    def test_fingerprint_numpy_shape(self):
        """Verify numpy arrays have correct shape."""
        exporter = FingerprintExporter()
        buf = exporter.export(SAMPLE_RESULTS)

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            # Morgan: 3 mols x 2048 bits
            morgan_npy = np.load(BytesIO(zf.read("morgan.npy")))
            assert morgan_npy.shape == (3, 2048)

            # MACCS: 3 mols x 167 bits
            maccs_npy = np.load(BytesIO(zf.read("maccs.npy")))
            assert maccs_npy.shape == (3, 167)

            # RDKit: 3 mols x 2048 bits
            rdkit_npy = np.load(BytesIO(zf.read("rdkit.npy")))
            assert rdkit_npy.shape == (3, 2048)

    def test_fingerprint_npz_loadable(self):
        """Verify npz files can be loaded."""
        exporter = FingerprintExporter()
        buf = exporter.export(SAMPLE_RESULTS)

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            data = np.load(BytesIO(zf.read("morgan.npz")))
            assert "fingerprints" in data.files
            assert data["fingerprints"].shape == (3, 2048)

    def test_fingerprint_exporter_empty_results(self):
        """Feed empty list, verify graceful handling."""
        exporter = FingerprintExporter()
        buf = exporter.export([])

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            names = zf.namelist()
            assert len(names) == 9

            morgan_npy = np.load(BytesIO(zf.read("morgan.npy")))
            assert morgan_npy.shape[0] == 0

    def test_fingerprint_media_type(self):
        """Verify media type and extension."""
        exporter = FingerprintExporter()
        assert exporter.media_type == "application/zip"
        assert exporter.file_extension == "zip"

    def test_fingerprint_factory_registration(self):
        """Verify fingerprint exporter registered in factory."""
        exporter = ExporterFactory.create(ExportFormat.FINGERPRINT)
        assert isinstance(exporter, FingerprintExporter)


class TestDedupExporter:
    """Test DedupExporter."""

    def test_dedup_exporter_with_duplicates(self):
        """Feed results with 2 identical SMILES, verify grouping."""
        exporter = DedupExporter()
        buf = exporter.export(SAMPLE_RESULTS)

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            # Summary CSV
            summary_csv = zf.read("dedup_summary.csv").decode("utf-8")
            summary_df = pd.read_csv(BytesIO(summary_csv.encode("utf-8")))

            # CCO appears twice (indices 0 and 2), benzene once, INVALID once
            # So we should have 3 groups
            assert len(summary_df) == 3

            # Find the CCO group
            cco_groups = summary_df[summary_df["representative_smiles"] == "CCO"]
            assert len(cco_groups) == 1
            assert cco_groups.iloc[0]["group_size"] == 2

            # Annotated CSV
            annotated_csv = zf.read("dedup_annotated.csv").decode("utf-8")
            annotated_df = pd.read_csv(BytesIO(annotated_csv.encode("utf-8")))

            # All 4 molecules present
            assert len(annotated_df) == 4

            # CCO molecules should have same group_id
            cco_rows = annotated_df[annotated_df["smiles"] == "CCO"]
            assert len(cco_rows) == 2
            group_ids = cco_rows["group_id"].unique()
            assert len(group_ids) == 1  # Same group

            # First occurrence should be representative
            reps = cco_rows[cco_rows["is_representative"] == True]  # noqa: E712
            assert len(reps) == 1
            assert reps.iloc[0]["index"] == 0  # min index

    def test_dedup_exporter_no_duplicates(self):
        """Feed 3 unique molecules, verify 3 groups of size 1."""
        unique_results = [
            {
                "smiles": "CCO", "name": "Ethanol", "index": 0, "status": "success",
                "validation": {"canonical_smiles": "CCO", "overall_score": 95},
                "alerts": {}, "scoring": {},
            },
            {
                "smiles": "c1ccccc1", "name": "Benzene", "index": 1, "status": "success",
                "validation": {"canonical_smiles": "c1ccccc1", "overall_score": 85},
                "alerts": {}, "scoring": {},
            },
            {
                "smiles": "CC(=O)O", "name": "Acetic acid", "index": 2, "status": "success",
                "validation": {"canonical_smiles": "CC(=O)O", "overall_score": 90},
                "alerts": {}, "scoring": {},
            },
        ]
        exporter = DedupExporter()
        buf = exporter.export(unique_results)

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            summary_csv = zf.read("dedup_summary.csv").decode("utf-8")
            summary_df = pd.read_csv(BytesIO(summary_csv.encode("utf-8")))
            assert len(summary_df) == 3
            assert (summary_df["group_size"] == 1).all()

    def test_dedup_media_type(self):
        """Verify media type and extension."""
        exporter = DedupExporter()
        assert exporter.media_type == "application/zip"
        assert exporter.file_extension == "zip"

    def test_dedup_factory_registration(self):
        """Verify dedup exporter registered in factory."""
        exporter = ExporterFactory.create(ExportFormat.DEDUP)
        assert isinstance(exporter, DedupExporter)


class TestScaffoldExporter:
    """Test ScaffoldExporter."""

    def test_scaffold_exporter_groups(self):
        """Feed molecules with shared/different scaffolds, verify assignment."""
        exporter = ScaffoldExporter()
        buf = exporter.export(SAMPLE_RESULTS)

        buf.seek(0)
        content = buf.read().decode("utf-8")
        df = pd.read_csv(BytesIO(content.encode("utf-8")))

        assert "scaffold_smiles" in df.columns
        assert "scaffold_group" in df.columns

        # Benzene has a scaffold, ethanol is acyclic
        benzene_row = df[df["name"] == "Benzene"]
        assert len(benzene_row) == 1
        assert benzene_row.iloc[0]["scaffold_smiles"] != ""
        assert benzene_row.iloc[0]["scaffold_group"] > 0

    def test_scaffold_exporter_acyclic(self):
        """Feed acyclic molecule, verify empty scaffold and group 0."""
        acyclic_results = [
            {
                "smiles": "CCCCC", "name": "Pentane", "index": 0, "status": "success",
                "validation": {"canonical_smiles": "CCCCC", "overall_score": 80},
                "alerts": {}, "scoring": {},
            },
        ]
        exporter = ScaffoldExporter()
        buf = exporter.export(acyclic_results)

        buf.seek(0)
        content = buf.read().decode("utf-8")
        df = pd.read_csv(BytesIO(content.encode("utf-8")))

        assert len(df) == 1
        assert df.iloc[0]["scaffold_smiles"] == "" or pd.isna(df.iloc[0]["scaffold_smiles"])
        assert df.iloc[0]["scaffold_group"] == 0

    def test_scaffold_media_type(self):
        """Verify media type and extension."""
        exporter = ScaffoldExporter()
        assert exporter.media_type == "text/csv"
        assert exporter.file_extension == "csv"

    def test_scaffold_factory_registration(self):
        """Verify scaffold exporter registered in factory."""
        exporter = ExporterFactory.create(ExportFormat.SCAFFOLD)
        assert isinstance(exporter, ScaffoldExporter)


class TestPropertyMatrixExporter:
    """Test PropertyMatrixExporter."""

    def test_property_matrix_csv_and_excel(self):
        """Feed 3 molecules, verify zip has 2 files with correct structure."""
        exporter = PropertyMatrixExporter()
        buf = exporter.export(SAMPLE_RESULTS[:3])  # Exclude error molecule

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            names = zf.namelist()
            assert "properties_flat.csv" in names
            assert "properties.xlsx" in names

            # Check CSV
            csv_content = zf.read("properties_flat.csv").decode("utf-8")
            csv_df = pd.read_csv(BytesIO(csv_content.encode("utf-8")))
            assert len(csv_df) == 3
            assert "smiles" in csv_df.columns
            assert "mw" in csv_df.columns
            assert "overall_score" in csv_df.columns
            assert "alert_count" in csv_df.columns

            # Check Excel has 4 sheets
            excel_buf = BytesIO(zf.read("properties.xlsx"))
            xls = pd.ExcelFile(excel_buf)
            assert "Descriptors" in xls.sheet_names
            assert "Scores" in xls.sheet_names
            assert "Alerts" in xls.sheet_names
            assert "Properties" in xls.sheet_names

    def test_property_matrix_empty_results(self):
        """Feed empty list, verify graceful handling."""
        exporter = PropertyMatrixExporter()
        buf = exporter.export([])

        buf.seek(0)
        with zipfile.ZipFile(buf, "r") as zf:
            names = zf.namelist()
            assert "properties_flat.csv" in names
            assert "properties.xlsx" in names

    def test_property_matrix_media_type(self):
        """Verify media type and extension."""
        exporter = PropertyMatrixExporter()
        assert exporter.media_type == "application/zip"
        assert exporter.file_extension == "zip"

    def test_property_matrix_factory_registration(self):
        """Verify property matrix exporter registered in factory."""
        exporter = ExporterFactory.create(ExportFormat.PROPERTY_MATRIX)
        assert isinstance(exporter, PropertyMatrixExporter)
