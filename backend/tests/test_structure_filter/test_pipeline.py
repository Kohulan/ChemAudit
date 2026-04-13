"""
Unit tests for the GenChem filter pipeline.

Tests cover:
- PRESETS configuration (4 presets with correct D-12 thresholds)
- Preset weight vector differentiation (D-15)
- filter_batch() stage count (6 stages)
- Per-stage rejection behavior (parse, property, dedup)
- Fragment-like and permissive preset behavior
"""

from __future__ import annotations

from app.services.genchem.filter_config import PRESETS
from app.services.genchem.filter_pipeline import filter_batch

# ---------------------------------------------------------------------------
# Preset configuration tests
# ---------------------------------------------------------------------------


def test_presets_defined():
    """PRESETS dict must have exactly 4 keys."""
    assert set(PRESETS.keys()) == {"drug_like", "lead_like", "fragment_like", "permissive"}


def test_drug_like_preset_values():
    """Verify drug_like preset matches D-12 specification."""
    dl = PRESETS["drug_like"]
    assert dl.min_mw == 200.0
    assert dl.max_mw == 500.0
    assert dl.max_tpsa == 140.0
    assert dl.max_rot_bonds == 10
    assert dl.max_sa_score == 5.0
    assert dl.use_pains is True
    assert dl.use_brenk is True
    assert dl.use_kazius is True
    assert dl.use_nibr is False


def test_lead_like_preset_values():
    """Verify lead_like preset matches D-12 specification."""
    ll = PRESETS["lead_like"]
    assert ll.max_mw == 350.0
    assert ll.max_logp == 3.5
    assert ll.max_rot_bonds == 7
    assert ll.use_nibr is True


def test_fragment_like_preset_values():
    """Verify fragment_like preset matches D-12 specification."""
    fl = PRESETS["fragment_like"]
    assert fl.max_mw == 300.0
    assert fl.max_rings == 3
    assert fl.max_sa_score == 3.0
    assert fl.use_pains is True
    assert fl.use_brenk is False
    assert fl.use_kazius is False


def test_permissive_preset_values():
    """Verify permissive preset matches D-12 specification."""
    p = PRESETS["permissive"]
    assert p.max_mw == 800.0
    assert p.max_sa_score == 7.0
    assert p.use_pains is False
    assert p.use_brenk is False
    assert p.use_kazius is False
    assert p.use_nibr is False


def test_preset_weights_differ():
    """
    Per D-15: each preset must have a unique (validity, qed, alert_free, sa) weight vector.
    No two presets may share the same weight combination.
    """
    weight_vectors = {
        name: (
            cfg.weight_validity,
            cfg.weight_qed,
            cfg.weight_alert_free,
            cfg.weight_sa,
        )
        for name, cfg in PRESETS.items()
    }
    # All 4 weight vectors must be distinct
    assert len(set(weight_vectors.values())) == 4, (
        f"D-15 violation: some presets share weight vectors: {weight_vectors}"
    )


def test_drug_like_weight_vector():
    """Verify drug_like weight vector per D-15."""
    dl = PRESETS["drug_like"]
    assert (dl.weight_validity, dl.weight_qed, dl.weight_alert_free, dl.weight_sa) == (
        0.3,
        0.3,
        0.2,
        0.2,
    )


def test_lead_like_weight_vector():
    """Verify lead_like weight vector per D-15 (emphasizes QED)."""
    ll = PRESETS["lead_like"]
    assert (ll.weight_validity, ll.weight_qed, ll.weight_alert_free, ll.weight_sa) == (
        0.2,
        0.4,
        0.2,
        0.2,
    )


def test_fragment_like_weight_vector():
    """Verify fragment_like weight vector per D-15 (emphasizes alert_free and SA)."""
    fl = PRESETS["fragment_like"]
    assert (fl.weight_validity, fl.weight_qed, fl.weight_alert_free, fl.weight_sa) == (
        0.2,
        0.2,
        0.3,
        0.3,
    )


def test_permissive_weight_vector():
    """Verify permissive weight vector per D-15 (emphasizes validity)."""
    p = PRESETS["permissive"]
    assert (p.weight_validity, p.weight_qed, p.weight_alert_free, p.weight_sa) == (
        0.4,
        0.3,
        0.1,
        0.2,
    )


def test_fragment_like_alerts_pains_only():
    """Fragment-like config only enables PAINS (Brenk/Kazius/NIBR are disabled)."""
    fl = PRESETS["fragment_like"]
    assert fl.use_pains is True
    assert fl.use_brenk is False
    assert fl.use_kazius is False
    assert fl.use_nibr is False


# ---------------------------------------------------------------------------
# filter_batch() tests
# ---------------------------------------------------------------------------


def test_filter_single_stages():
    """filter_batch with 1 valid SMILES returns 6 StageResult entries."""
    # Use permissive preset so CCO isn't rejected by property thresholds
    result = filter_batch(["CCO"], PRESETS["permissive"])
    assert result.input_count == 1
    assert len(result.stages) == 6


def test_filter_batch_counts(mixed_smiles_list):
    """filter_batch with mixed input tracks input_count correctly."""
    result = filter_batch(mixed_smiles_list, PRESETS["drug_like"])
    assert result.input_count == 5


def test_parse_stage_rejects_invalid(mixed_smiles_list):
    """'invalid###' must be rejected at the parse stage."""
    result = filter_batch(mixed_smiles_list, PRESETS["drug_like"])
    invalid_mol = next(m for m in result.molecules if m.smiles == "invalid###")
    assert invalid_mol.status == "rejected"
    assert invalid_mol.failed_at == "parse"


def test_empty_mol_rejected(mixed_smiles_list):
    """Empty string SMILES must be rejected at the parse stage (Pitfall 1)."""
    result = filter_batch(mixed_smiles_list, PRESETS["drug_like"])
    empty_mol = next(m for m in result.molecules if m.smiles == "")
    assert empty_mol.status == "rejected"
    assert empty_mol.failed_at == "parse"


def test_drug_like_mw_reject(large_mw_smiles):
    """Large MW molecule must be rejected at the property stage with drug_like preset."""
    result = filter_batch([large_mw_smiles], PRESETS["drug_like"])
    assert result.molecules[0].status == "rejected"
    assert result.molecules[0].failed_at == "property"


def test_permissive_passes_large_mw(large_mw_smiles):
    """Large MW molecule should NOT be rejected at property stage with permissive preset."""
    result = filter_batch([large_mw_smiles], PRESETS["permissive"])
    # The molecule should pass property stage (may still fail SA or alerts)
    property_stage = next(s for s in result.stages if s.stage_name == "property")
    assert property_stage.rejected_count == 0, (
        f"Large MW molecule rejected at property with permissive preset: "
        f"{result.molecules[0].rejection_reason}"
    )


def test_dedup_stage():
    """Second occurrence of identical SMILES must be marked as duplicate."""
    # Use ibuprofen (MW~206, passes all drug_like filters) to ensure molecules reach dedup stage
    ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    result = filter_batch([ibuprofen, ibuprofen], PRESETS["drug_like"])
    # Find the molecule with duplicate status
    duplicates = [m for m in result.molecules if m.status == "duplicate"]
    assert len(duplicates) >= 1
    assert duplicates[0].failed_at == "dedup"
    assert "Duplicate of row" in (duplicates[0].rejection_reason or "")


def test_dedup_marks_second_occurrence():
    """Dedup stage: first occurrence passes, second is marked as duplicate of row 1."""
    # Use ibuprofen (MW~206, passes all drug_like filters) to ensure molecules reach dedup stage
    ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    result = filter_batch([ibuprofen, ibuprofen], PRESETS["drug_like"])
    # First molecule (index 0) should pass
    assert result.molecules[0].status == "passed"
    # Second molecule (index 1) should be duplicate
    assert result.molecules[1].status == "duplicate"
    assert result.molecules[1].failed_at == "dedup"
    # 1-indexed row reference: first occurrence is row 1
    assert "row 1" in (result.molecules[1].rejection_reason or "")


def test_all_pass_ibuprofen():
    """Ibuprofen (MW~206, no alerts) should pass all stages with drug_like preset."""
    ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    result = filter_batch([ibuprofen], PRESETS["drug_like"])
    assert result.output_count == 1
    assert result.molecules[0].status == "passed"


def test_output_count_reflects_passed():
    """output_count must equal the number of molecules with status='passed'."""
    ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    result = filter_batch([ibuprofen, "invalid###", "c1ccccc1"], PRESETS["drug_like"])
    passed_count = sum(1 for m in result.molecules if m.status == "passed")
    assert result.output_count == passed_count


def test_stage_names_present():
    """All 6 expected stage names must be present in result.stages."""
    result = filter_batch(["CCO"], PRESETS["permissive"])
    stage_names = {s.stage_name for s in result.stages}
    expected = {"parse", "valence", "alerts", "property", "sa", "dedup"}
    assert expected == stage_names


def test_stage_indices_sequential():
    """Stage indices must be 1 through 6 in order."""
    result = filter_batch(["CCO"], PRESETS["permissive"])
    indices = [s.stage_index for s in result.stages]
    assert indices == list(range(1, 7))


def test_valence_stage_no_rejections():
    """Valence stage must always have 0 rejections (handled by RDKit during parse)."""
    result = filter_batch(["CCO", "c1ccccc1", "CC(=O)O"], PRESETS["permissive"])
    valence_stage = next(s for s in result.stages if s.stage_name == "valence")
    assert valence_stage.rejected_count == 0


def test_novelty_stage_not_in_default():
    """Novelty stage must NOT appear in results when enable_novelty=False (default)."""
    result = filter_batch(["CCO"], PRESETS["permissive"])
    novelty_stages = [s for s in result.stages if s.stage_name == "novelty"]
    assert len(novelty_stages) == 0
