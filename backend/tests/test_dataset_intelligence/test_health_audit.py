"""Tests for the dataset health audit service module."""

from unittest.mock import patch

import pytest

from app.services.dataset_intelligence.health_audit import (
    DatasetHealthResult,
    DRUG_SPACE_REFERENCE,
    compute_health_score,
)
from rdkit import Chem
from rdkit.Chem import inchi as rdkit_inchi


def _make_molecule(smiles: str, index: int, properties: dict | None = None) -> dict:
    """Helper to build a molecule dict for compute_health_score."""
    mol = Chem.MolFromSmiles(smiles)
    inchikey = None
    if mol is not None:
        inchikey = rdkit_inchi.MolToInchiKey(mol)
    return {
        "index": index,
        "smiles": smiles,
        "mol": mol,
        "inchikey": inchikey,
        "properties": properties or {},
    }


class TestComputeHealthScore:
    """Tests for compute_health_score function."""

    def test_parsability_score_with_invalid(self):
        """3 valid + 1 invalid SMILES => parsability_score=0.75, overall in 0-100."""
        molecules = [
            _make_molecule("CCO", 0),
            _make_molecule("c1ccccc1", 1),
            _make_molecule("CC(=O)O", 2),
            _make_molecule("INVALID", 3),
        ]
        result = compute_health_score(molecules)
        assert isinstance(result, DatasetHealthResult)
        assert result.parsability_score == pytest.approx(0.75)
        assert 0 <= result.overall_score <= 100

    def test_all_valid_stereo_score(self):
        """All valid SMILES with no stereo issues => stereo_score=1.0."""
        molecules = [
            _make_molecule("CCO", 0),
            _make_molecule("c1ccccc1", 1),
            _make_molecule("CC(=O)O", 2),
        ]
        result = compute_health_score(molecules)
        assert result.parsability_score == pytest.approx(1.0)
        assert result.stereo_score == pytest.approx(1.0)

    def test_duplicate_inchikeys_uniqueness(self):
        """Duplicate InChIKeys (submit 'CCO' twice) => uniqueness_score < 1.0."""
        molecules = [
            _make_molecule("CCO", 0),
            _make_molecule("CCO", 1),
            _make_molecule("c1ccccc1", 2),
        ]
        result = compute_health_score(molecules)
        assert result.uniqueness_score < 1.0
        assert result.duplicate_count > 0
        assert len(result.dedup_groups) > 0

    def test_property_distributions_keys(self):
        """property_distributions has 'mw', 'logp', 'tpsa' with bins and counts."""
        molecules = [
            _make_molecule("CCO", 0),
            _make_molecule("c1ccccc1", 1),
        ]
        result = compute_health_score(molecules)
        assert "mw" in result.property_distributions
        assert "logp" in result.property_distributions
        assert "tpsa" in result.property_distributions
        for key in ["mw", "logp", "tpsa"]:
            dist = result.property_distributions[key]
            assert "bins" in dist
            assert "counts" in dist

    def test_std_sample_size_capped_at_500(self):
        """std_sample_size <= 500 when given many molecules."""
        # Use 600 copies of CCO
        molecules = [_make_molecule("CCO", i) for i in range(600)]
        # Mock compare_pipelines to avoid heavy computation
        with patch(
            "app.services.dataset_intelligence.health_audit.compare_pipelines"
        ) as mock_cp:
            mock_cp.return_value = {"all_agree": True}
            result = compute_health_score(molecules)
        assert result.std_sample_size <= 500

    def test_default_weights_sum_to_one(self):
        """Default weights sum to 1.0."""
        molecules = [_make_molecule("CCO", 0)]
        result = compute_health_score(molecules)
        assert sum(result.weights.values()) == pytest.approx(1.0)

    def test_overall_score_is_weighted_composite(self):
        """Overall score is a weighted composite of sub-scores * 100."""
        molecules = [
            _make_molecule("CCO", 0),
            _make_molecule("c1ccccc1", 1),
        ]
        result = compute_health_score(molecules)
        expected = (
            result.parsability_score * result.weights["parsability"]
            + result.stereo_score * result.weights["stereo"]
            + result.uniqueness_score * result.weights["uniqueness"]
            + result.alert_score * result.weights["alerts"]
            + result.std_consistency_score * result.weights["std_consistency"]
        ) * 100
        assert result.overall_score == pytest.approx(expected, abs=0.01)

    def test_molecule_count_tracking(self):
        """Molecule count and parse failure tracking."""
        molecules = [
            _make_molecule("CCO", 0),
            _make_molecule("INVALID", 1),
        ]
        result = compute_health_score(molecules)
        assert result.molecule_count == 2
        assert result.parse_failures == 1


class TestDrugSpaceReference:
    """Tests for the DRUG_SPACE_REFERENCE constant."""

    def test_has_required_keys(self):
        """DRUG_SPACE_REFERENCE has mw, logp, tpsa keys."""
        assert "mw" in DRUG_SPACE_REFERENCE
        assert "logp" in DRUG_SPACE_REFERENCE
        assert "tpsa" in DRUG_SPACE_REFERENCE

    def test_has_bin_edges_and_counts(self):
        """Each entry has bin_edges and reference_counts."""
        for key in ["mw", "logp", "tpsa"]:
            assert "bin_edges" in DRUG_SPACE_REFERENCE[key]
            assert "reference_counts" in DRUG_SPACE_REFERENCE[key]
