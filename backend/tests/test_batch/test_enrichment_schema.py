"""Tests for batch enrichment schema fields (Task 7)."""

from app.schemas.batch import BatchResultItem


def test_batch_result_item_has_profiling_field():
    item = BatchResultItem(
        smiles="CCO", index=0, status="success", profiling={"pfi": {"pfi": 3.2}}
    )
    assert item.profiling is not None


def test_batch_result_item_has_safety_assessment_field():
    item = BatchResultItem(
        smiles="CCO", index=0, status="success", safety_assessment={"cyp_result": {}}
    )
    assert item.safety_assessment is not None


def test_batch_result_item_enrichment_fields_optional():
    item = BatchResultItem(smiles="CCO", index=0, status="success")
    assert item.profiling is None
    assert item.safety_assessment is None
