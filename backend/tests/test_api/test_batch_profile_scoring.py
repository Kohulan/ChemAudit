"""Test that profile scoring is applied when profile data is in safety_options."""
from app.services.batch.tasks import _process_single_molecule


def test_profile_score_in_result_when_profile_provided():
    """When safety_options includes profile data, result should have scoring.profile."""
    mol_data = {"smiles": "c1ccc(CC(=O)O)cc1", "name": "test", "index": 0}
    safety_options = {
        "profile_id": 1,
        "profile_name": "Drug-like (Lipinski)",
        "profile_thresholds": {
            "mw": {"min": 0, "max": 500},
            "logp": {"min": -5, "max": 5},
            "hbd": {"min": 0, "max": 5},
            "hba": {"min": 0, "max": 10},
        },
        "profile_weights": {"mw": 1.0, "logp": 1.0, "hbd": 1.0, "hba": 1.0},
    }
    result = _process_single_molecule(mol_data, safety_options=safety_options)
    profile = result["scoring"]["profile"]
    assert profile["profile_id"] == 1
    assert profile["profile_name"] == "Drug-like (Lipinski)"
    assert 0 <= profile["score"] <= 100
    assert "mw" in profile["properties"]
    assert profile["properties"]["mw"]["in_range"] is True


def test_no_profile_score_when_profile_not_provided():
    """When no profile data in safety_options, result should not have scoring.profile."""
    mol_data = {"smiles": "c1ccc(CC(=O)O)cc1", "name": "test", "index": 0}
    result = _process_single_molecule(mol_data)
    assert "profile" not in result.get("scoring", {})
