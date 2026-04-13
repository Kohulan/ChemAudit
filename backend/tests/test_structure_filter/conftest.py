"""Fixtures for GenChem filter pipeline and scorer tests."""

import pytest
from rdkit import Chem


@pytest.fixture
def valid_smiles_list():
    """A list of valid, drug-like SMILES."""
    return ["CCO", "c1ccccc1", "CC(=O)O", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"]


@pytest.fixture
def mixed_smiles_list():
    """Mixed list including invalid SMILES, empty string, and duplicate."""
    return ["CCO", "invalid###", "c1ccccc1", "", "CCO"]  # last is duplicate of first


@pytest.fixture
def large_mw_smiles():
    """
    A glycoside-like SMILES with MW ~560 and LogP ~1.0 that passes all structural
    alert screens (no PAINS/Brenk/Kazius matches) but exceeds the drug_like MW
    threshold of 500. SA score ~4.6 (passes permissive max_sa=7.0).
    Used to test property-stage rejection specifically.
    """
    return "OC(CC(=O)c1ccc(OC2OC(C)C(O)C(N(C)C)C2)cc1)c1ccc(OC2OC(CO)C(O)C(N)C2)cc1"


@pytest.fixture
def ethanol_mol():
    """Ethanol as an RDKit molecule."""
    return Chem.MolFromSmiles("CCO")
