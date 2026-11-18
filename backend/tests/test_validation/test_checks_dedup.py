"""Duplicate check names must not run a check more than once (security scan finding F14)."""

import pytest
from pydantic import ValidationError
from rdkit import Chem

from app.schemas.validation import ValidationRequest
from app.services.validation.engine import ValidationEngine


def test_engine_runs_each_check_once():
    results, _ = ValidationEngine().validate(Chem.MolFromSmiles("CCO"), ["valence"] * 5)
    assert [r.check_name for r in results] == ["valence"]


def test_request_deduplicates_checks():
    assert ValidationRequest(molecule="CCO", checks=["valence", "valence"]).checks == ["valence"]


def test_request_rejects_oversized_checks_list():
    with pytest.raises(ValidationError):
        ValidationRequest(molecule="CCO", checks=["valence"] * 51)
