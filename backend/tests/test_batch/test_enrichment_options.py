"""Tests for batch enrichment options in Celery task (Task 8)."""

from app.services.batch.tasks import _process_single_molecule


def test_process_molecule_with_profiling():
    mol_data = {"smiles": "c1ccccc1", "name": "benzene", "index": 0}
    opts = {"include_profiling": True}
    result = _process_single_molecule(mol_data, safety_options=opts)
    assert result["status"] == "success"
    assert result.get("profiling") is not None
    assert "pfi" in result["profiling"]


def test_process_molecule_with_safety_assessment():
    mol_data = {"smiles": "c1ccccc1", "name": "benzene", "index": 0}
    opts = {"include_safety_assessment": True}
    result = _process_single_molecule(mol_data, safety_options=opts)
    assert result["status"] == "success"
    assert result.get("safety_assessment") is not None


def test_process_molecule_without_enrichment():
    mol_data = {"smiles": "c1ccccc1", "name": "benzene", "index": 0}
    result = _process_single_molecule(mol_data)
    assert result.get("profiling") is None
    assert result.get("safety_assessment") is None
