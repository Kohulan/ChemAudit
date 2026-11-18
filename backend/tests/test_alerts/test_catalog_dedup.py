"""Duplicate catalog names must not repeat full screening (security scan finding F13)."""

import pytest
from pydantic import ValidationError
from rdkit import Chem

from app.schemas.alerts import AlertScreenRequest
from app.services.alerts.alert_manager import AlertManager


def test_catalogs_deduplicated_case_insensitively():
    req = AlertScreenRequest(molecule="CCO", catalogs=["PAINS", "pains", "PAINS"])
    assert req.catalogs == ["PAINS"]


def test_all_collapses_the_list():
    assert AlertScreenRequest(molecule="CCO", catalogs=["PAINS", "ALL"]).catalogs == ["ALL"]


def test_oversized_catalog_list_rejected():
    with pytest.raises(ValidationError):
        AlertScreenRequest(molecule="CCO", catalogs=["PAINS"] * 100)


def test_screen_deduplicates_catalogs():
    mol = Chem.MolFromSmiles("O=C(C=Cc1ccccc1)c1ccccc1")
    result = AlertManager().screen(mol, catalogs=["BRENK", "brenk", "BRENK"])
    assert result.screened_catalogs == ["BRENK"]
