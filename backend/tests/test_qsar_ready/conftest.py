"""
Fixtures for QSAR-Ready pipeline tests.
"""

import pytest

from app.services.qsar_ready.pipeline import QSARReadyConfig


@pytest.fixture
def aspirin_smiles() -> str:
    """Aspirin SMILES."""
    return "CC(=O)OC1=CC=CC=C1C(=O)O"


@pytest.fixture
def leucine_smiles() -> str:
    """L-leucine SMILES with stereo."""
    return "CC(C)C[C@H](N)C(=O)O"


@pytest.fixture
def r_leucine_smiles() -> str:
    """R-leucine (D-leucine) SMILES with stereo."""
    return "CC(C)C[C@@H](N)C(=O)O"


@pytest.fixture
def sodium_acetate_smiles() -> str:
    """Sodium acetate SMILES (contains metal)."""
    return "CC(=O)O[Na]"


@pytest.fixture
def deuterium_benzene_smiles() -> str:
    """Deuterium-labelled benzene SMILES."""
    return "[2H]c1ccccc1"


@pytest.fixture
def default_config() -> QSARReadyConfig:
    """Default QSAR-Ready configuration."""
    return QSARReadyConfig()


@pytest.fixture
def qsar_2d_config() -> QSARReadyConfig:
    """QSAR-2D preset configuration."""
    return QSARReadyConfig.qsar_2d()


@pytest.fixture
def qsar_3d_config() -> QSARReadyConfig:
    """QSAR-3D preset configuration."""
    return QSARReadyConfig.qsar_3d()


@pytest.fixture
def minimal_config() -> QSARReadyConfig:
    """Minimal preset configuration."""
    return QSARReadyConfig.minimal()
