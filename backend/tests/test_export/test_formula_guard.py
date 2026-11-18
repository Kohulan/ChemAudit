"""
Exports must not turn user text into spreadsheet formulas or forged SDF records.

Regression tests for security scan findings F12, F17, F19, F21, F22 and F28.
"""

import copy
import csv
import importlib
import io
import zipfile

import openpyxl
import pandas as pd
import pytest

from app.services.export import CSVExporter, ExcelExporter, SDFExporter
from app.services.export.formula_guard import neutralize_cell, neutralize_frame
from tests.test_export.test_exporters import SAMPLE_RESULTS

EVIL = '=HYPERLINK("http://x","y")'


@pytest.mark.parametrize("value", ["=1+1", "+1A", "@SUM(A1)", "\tx", "\rx", "-A1", "{=A1}"])
def test_neutralize_cell_prefixes_formula_starts(value):
    assert neutralize_cell(value) == "'" + value


@pytest.mark.parametrize("value", ["-1.5", "+3", "1e-5", "CCO", "", "(-)-menthol", 5, None, 2.5])
def test_neutralize_cell_leaves_data_alone(value):
    assert neutralize_cell(value) == value


def test_neutralize_frame_covers_values_and_headers():
    df = pd.DataFrame({"=evil": ["=1+1", "ok"], "n": [-1, 2]})
    out = neutralize_frame(df)
    assert list(out.columns) == ["'=evil", "n"]
    assert out["'=evil"].tolist() == ["'=1+1", "ok"]
    assert out["n"].tolist() == [-1, 2]


def _evil_results():
    results = copy.deepcopy(SAMPLE_RESULTS)
    results[0]["name"] = EVIL
    return results


def _xlsx_has_formula(data: bytes) -> bool:
    wb = openpyxl.load_workbook(io.BytesIO(data))
    return any(
        c.data_type == "f" or (isinstance(c.value, str) and c.value.startswith("="))
        for ws in wb.worksheets
        for row in ws.iter_rows()
        for c in row
    )


def _csv_has_formula_cell(text: str) -> bool:
    return any(cell.startswith("=") for row in csv.reader(io.StringIO(text)) for cell in row)


def test_csv_export_neutralizes_name():
    text = CSVExporter().export(_evil_results()).getvalue().decode()
    assert "'=HYPERLINK" in text
    assert not _csv_has_formula_cell(text)


@pytest.mark.parametrize("layout", ["single", "multi"])
def test_excel_export_has_no_formula_cells(layout):
    data = ExcelExporter(sheet_layout=layout).export(_evil_results()).getvalue()
    assert not _xlsx_has_formula(data)


@pytest.mark.parametrize(
    "module,cls",
    [
        ("property_matrix_exporter", "PropertyMatrixExporter"),
        ("scaffold_exporter", "ScaffoldExporter"),
        ("dedup_exporter", "DedupExporter"),
        ("fingerprint_exporter", "FingerprintExporter"),
    ],
)
def test_zip_and_csv_exporters_have_no_formula_cells(module, cls):
    exporter = getattr(importlib.import_module(f"app.services.export.{module}"), cls)()
    data = exporter.export(_evil_results()).getvalue()
    if not zipfile.is_zipfile(io.BytesIO(data)):
        assert not _csv_has_formula_cell(data.decode())
        return
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        for name in zf.namelist():
            blob = zf.read(name)
            if name.endswith(".xlsx"):
                assert not _xlsx_has_formula(blob), name
            elif name.endswith(".csv"):
                assert not _csv_has_formula_cell(blob.decode()), name


def test_sdf_export_name_newlines_cannot_forge_records():
    results = copy.deepcopy(SAMPLE_RESULTS)
    results[0]["name"] = "a\nM  END\n> <overall_score>\n100\n\n$$$$\nforged"
    clean = SDFExporter().export(copy.deepcopy(SAMPLE_RESULTS)).getvalue().decode()
    text = SDFExporter().export(results).getvalue().decode()
    # A record ends at a line that is exactly "$$$$".
    assert text.splitlines().count("$$$$") == clean.splitlines().count("$$$$")


async def _body(response) -> str:
    if hasattr(response, "body_iterator"):
        chunks = [c async for c in response.body_iterator]
        return "".join(c if isinstance(c, str) else c.decode() for c in chunks)
    return response.body.decode()


async def test_qsar_csv_neutralizes_cells():
    from app.services.qsar_ready.builders import _build_csv_response

    resp = _build_csv_response("job", [{"original_smiles": "=1+1", "status": "ok"}])
    assert not _csv_has_formula_cell(await _body(resp))


async def test_structure_filter_csv_neutralizes_cells():
    from app.api.routes.structure_filter import _build_full_csv_response

    resp = _build_full_csv_response("job", [{"smiles": "=1+1", "status": "ok"}])
    assert not _csv_has_formula_cell(await _body(resp))


async def test_curated_dataset_csv_neutralizes_cells_and_headers(client):
    import json
    from unittest.mock import MagicMock, patch

    rows = [{"smiles": "CCO", "=cmd|' /C calc'!A0": "=1+1", "activity": "-1.5"}]
    redis = MagicMock()
    redis.get.return_value = json.dumps({"curated_csv_rows": rows})
    job = "123e4567-e89b-12d3-a456-426614174000"
    with patch("app.api.routes.dataset_intelligence._get_redis", return_value=redis):
        r = await client.get(f"/api/v1/dataset/{job}/download/csv")
    assert r.status_code == 200
    assert not _csv_has_formula_cell(r.text)
    assert ",-1.5" in r.text
