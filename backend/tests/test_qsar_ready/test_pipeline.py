"""
Unit tests for the QSAR-Ready pipeline service.

Tests cover:
- All 10 pipeline steps
- 3 preset configurations (QSAR-2D, QSAR-3D, Minimal)
- InChIKey tracking (original vs standardized, inchikey_changed flag)
- Error handling (invalid SMILES, null molecule guards)
- Batch processing with deduplication
"""

import pytest

from app.services.qsar_ready.pipeline import (
    QSARReadyConfig,
    QSARReadyResult,
    QSARStepResult,
    qsar_ready_batch,
    qsar_ready_single,
)


class TestAspirinAllSteps:
    """Test aspirin through the default pipeline."""

    def test_aspirin_all_steps(self, aspirin_smiles, default_config):
        """qsar_ready_single returns status='ok' with 10 step results for aspirin."""
        result = qsar_ready_single(aspirin_smiles, default_config)

        assert result.status == "ok"
        assert len(result.steps) == 10
        assert result.curated_smiles is not None
        assert result.curated_smiles != ""

    def test_aspirin_step_names(self, aspirin_smiles, default_config):
        """All 10 steps have correct step_name and step_index."""
        result = qsar_ready_single(aspirin_smiles, default_config)

        expected_names = [
            "parse",
            "metals",
            "desalt",
            "normalize",
            "neutralize",
            "tautomer",
            "stereo",
            "isotope",
            "filter",
            "canonical",
        ]
        for i, (step, expected_name) in enumerate(zip(result.steps, expected_names)):
            assert step.step_name == expected_name, (
                f"Step {i+1} name mismatch: got {step.step_name}, expected {expected_name}"
            )
            assert step.step_index == i + 1, (
                f"Step {i+1} index mismatch: got {step.step_index}"
            )

    def test_aspirin_step_has_required_fields(self, aspirin_smiles, default_config):
        """Each step result has required fields."""
        result = qsar_ready_single(aspirin_smiles, default_config)

        for step in result.steps:
            assert hasattr(step, "step_name")
            assert hasattr(step, "step_index")
            assert hasattr(step, "enabled")
            assert hasattr(step, "status")
            assert step.status in ("applied", "no_change", "skipped", "error")


class TestInChIKeyTracking:
    """Test InChIKey computation and change detection."""

    def test_inchikey_unchanged(self, aspirin_smiles, default_config):
        """Aspirin has inchikey_changed=False through default config."""
        result = qsar_ready_single(aspirin_smiles, default_config)

        assert result.original_inchikey is not None
        assert result.standardized_inchikey is not None
        assert result.inchikey_changed is False
        assert result.original_inchikey == result.standardized_inchikey

    def test_inchikey_capture_order(self, sodium_acetate_smiles, default_config):
        """Sodium acetate has inchikey_changed=True (metal disconnect changes structure)."""
        result = qsar_ready_single(sodium_acetate_smiles, default_config)

        # Metal disconnect + desalt changes the molecule: [Na+] removed, leaving acetic acid
        assert result.original_inchikey is not None
        assert result.standardized_inchikey is not None
        assert result.inchikey_changed is True
        assert result.original_inchikey != result.standardized_inchikey

    def test_inchikeys_are_strings(self, aspirin_smiles, default_config):
        """InChIKeys are proper non-empty strings."""
        result = qsar_ready_single(aspirin_smiles, default_config)

        assert isinstance(result.original_inchikey, str)
        assert len(result.original_inchikey) == 27  # Standard InChIKey length
        assert isinstance(result.standardized_inchikey, str)
        assert len(result.standardized_inchikey) == 27


class TestInvalidSMILES:
    """Test handling of invalid input."""

    def test_invalid_smiles(self, default_config):
        """Invalid SMILES returns status='rejected' with step 1 error."""
        result = qsar_ready_single("not_valid_smiles", default_config)

        assert result.status == "rejected"
        assert len(result.steps) >= 1
        assert result.steps[0].status == "error"
        assert result.steps[0].step_name == "parse"

    def test_invalid_smiles_no_curated(self, default_config):
        """Invalid SMILES has no curated_smiles."""
        result = qsar_ready_single("not_valid_smiles", default_config)

        assert result.curated_smiles is None

    def test_empty_string_smiles(self, default_config):
        """Empty string SMILES returns rejected."""
        result = qsar_ready_single("", default_config)

        assert result.status == "rejected"
        assert result.steps[0].status == "error"


class TestStepSkipping:
    """Test step enable/disable functionality."""

    def test_step_skipped(self, aspirin_smiles, default_config):
        """With enable_metals=False, metals step has status='skipped' and enabled=False."""
        config = QSARReadyConfig(enable_metals=False)
        result = qsar_ready_single(aspirin_smiles, config)

        metals_step = result.steps[1]  # metals is step index 2 (0-indexed: 1)
        assert metals_step.step_name == "metals"
        assert metals_step.enabled is False
        assert metals_step.status == "skipped"

    def test_step_skipped_still_10_steps(self, aspirin_smiles):
        """Even with steps disabled, result still has 10 steps."""
        config = QSARReadyConfig(
            enable_metals=False, enable_neutralize=False, enable_tautomer=False
        )
        result = qsar_ready_single(aspirin_smiles, config)

        assert len(result.steps) == 10


class TestQSAR2DPreset:
    """Test QSAR-2D preset configuration."""

    def test_qsar_2d_config(self):
        """QSAR-2D preset has enable_stereo_strip=True, enable_isotope_strip=True."""
        config = QSARReadyConfig.qsar_2d()

        assert config.enable_stereo_strip is True
        assert config.enable_isotope_strip is True
        assert config.enable_metals is True
        assert config.enable_desalt is True
        assert config.enable_normalize is True
        assert config.enable_neutralize is True
        assert config.enable_tautomer is True

    def test_qsar_2d_strips_stereo(self, leucine_smiles, qsar_2d_config):
        """L-leucine curated_smiles has no @ character after QSAR-2D processing."""
        result = qsar_ready_single(leucine_smiles, qsar_2d_config)

        assert result.status == "ok"
        assert result.curated_smiles is not None
        assert "@" not in result.curated_smiles, (
            f"Stereo not stripped: {result.curated_smiles}"
        )


class TestQSAR3DPreset:
    """Test QSAR-3D preset configuration."""

    def test_qsar_3d_config(self):
        """QSAR-3D preset has enable_stereo_strip=False."""
        config = QSARReadyConfig.qsar_3d()

        assert config.enable_stereo_strip is False
        assert config.enable_isotope_strip is True

    def test_qsar_3d_preserves_stereo(self, leucine_smiles, qsar_3d_config):
        """L-leucine curated_smiles preserves @ character after QSAR-3D processing."""
        result = qsar_ready_single(leucine_smiles, qsar_3d_config)

        assert result.status == "ok"
        assert result.curated_smiles is not None
        assert "@" in result.curated_smiles, (
            f"Stereo was stripped unexpectedly: {result.curated_smiles}"
        )


class TestMinimalPreset:
    """Test Minimal preset configuration."""

    def test_minimal_config(self):
        """Minimal preset has expected disabled steps."""
        config = QSARReadyConfig.minimal()

        assert config.enable_metals is False
        assert config.enable_neutralize is False
        assert config.enable_tautomer is False
        assert config.enable_stereo_strip is False
        assert config.enable_isotope_strip is False

    def test_minimal_preset_steps_skipped(self, aspirin_smiles, minimal_config):
        """Minimal preset: metals/neutralize/tautomer/stereo/isotope steps are skipped."""
        result = qsar_ready_single(aspirin_smiles, minimal_config)

        assert result.status == "ok"

        step_map = {s.step_name: s for s in result.steps}

        for skipped_step in ["metals", "neutralize", "tautomer", "stereo", "isotope"]:
            assert step_map[skipped_step].status == "skipped", (
                f"Expected {skipped_step} to be skipped in minimal preset"
            )
            assert step_map[skipped_step].enabled is False


class TestIsotopeStrip:
    """Test isotope stripping step."""

    def test_isotope_strip(self, deuterium_benzene_smiles, qsar_2d_config):
        """Deuterium benzene with enable_isotope_strip=True produces SMILES without isotope labels."""
        result = qsar_ready_single(deuterium_benzene_smiles, qsar_2d_config)

        assert result.status == "ok"
        assert result.curated_smiles is not None
        # Deuterium [2H] should be stripped
        assert "[2H]" not in result.curated_smiles, (
            f"Isotope not stripped: {result.curated_smiles}"
        )

    def test_isotope_step_applied(self, deuterium_benzene_smiles, qsar_2d_config):
        """Isotope step shows 'applied' status for deuterium compound."""
        result = qsar_ready_single(deuterium_benzene_smiles, qsar_2d_config)

        isotope_step = result.steps[7]  # isotope is step index 8 (0-indexed: 7)
        assert isotope_step.step_name == "isotope"
        assert isotope_step.status == "applied"


class TestCompositionFilter:
    """Test composition filter step (step 9)."""

    def test_composition_filter_rejects_single_heavy_atom(self, default_config):
        """Sodium ion '[Na+]' after desalt (1 heavy atom) returns rejected with 'heavy atom' reason."""
        # Use a molecule that after desalt becomes too small
        # '[Na+]' alone has only 1 heavy atom, below min_heavy_atoms=3
        result = qsar_ready_single("[Na+]", default_config)

        assert result.status == "rejected"
        assert result.rejection_reason is not None
        assert "heavy atom" in result.rejection_reason.lower()

    def test_composition_filter_passes_aspirin(self, aspirin_smiles, default_config):
        """Aspirin passes all composition filters."""
        result = qsar_ready_single(aspirin_smiles, default_config)

        # Filter step should not reject aspirin
        assert result.status == "ok"
        filter_step = result.steps[8]  # filter is step index 9 (0-indexed: 8)
        assert filter_step.step_name == "filter"
        assert filter_step.status in ("applied", "no_change")


class TestNullMolGuard:
    """Test that steps handle None molecules gracefully."""

    def test_null_mol_guard_no_crash(self, default_config):
        """Pipeline handles problematic molecules gracefully without crashing."""
        # Use a molecule that might create issues
        # Test that even edge-case SMILES don't cause unhandled exceptions
        result = qsar_ready_single("C", default_config)  # Methane

        # Should return some result (might be rejected due to filter, but not crash)
        assert result is not None
        assert isinstance(result, QSARReadyResult)
        assert result.status in ("ok", "rejected", "error")


class TestBatchDeduplication:
    """Test batch processing with deduplication."""

    def test_batch_deduplication_helper(self, aspirin_smiles, default_config):
        """qsar_ready_batch with duplicate aspirin produces 1 ok + 1 duplicate."""
        results = qsar_ready_batch([aspirin_smiles, aspirin_smiles], default_config)

        assert len(results) == 2

        statuses = [r.status for r in results]
        assert "ok" in statuses
        assert "duplicate" in statuses

    def test_batch_deduplication_first_ok(self, aspirin_smiles, default_config):
        """First occurrence of duplicate molecule has status='ok'."""
        results = qsar_ready_batch([aspirin_smiles, aspirin_smiles], default_config)

        assert results[0].status == "ok"
        assert results[1].status == "duplicate"

    def test_batch_deduplication_reason(self, aspirin_smiles, default_config):
        """Duplicate entry has rejection_reason indicating which molecule it duplicates."""
        results = qsar_ready_batch([aspirin_smiles, aspirin_smiles], default_config)

        duplicate = results[1]
        assert duplicate.rejection_reason is not None
        assert "0" in duplicate.rejection_reason or "duplicate" in duplicate.rejection_reason.lower()

    def test_batch_multiple_molecules(self, aspirin_smiles, leucine_smiles, default_config):
        """Batch with different molecules returns all with status='ok'."""
        results = qsar_ready_batch([aspirin_smiles, leucine_smiles], default_config)

        assert len(results) == 2
        for result in results:
            assert result.status == "ok"

    def test_batch_invalid_smiles(self, aspirin_smiles, default_config):
        """Batch with invalid SMILES handles it gracefully."""
        results = qsar_ready_batch([aspirin_smiles, "not_valid", aspirin_smiles], default_config)

        assert len(results) == 3
        assert results[0].status == "ok"
        assert results[1].status == "rejected"
        # Third is duplicate of first
        assert results[2].status == "duplicate"


class TestEnantiomerDedup:
    """Test enantiomer deduplication via QSAR-2D stereo stripping."""

    def test_enantiomer_dedup(self, leucine_smiles, r_leucine_smiles, qsar_2d_config):
        """R-leucine and S-leucine are deduplicated when stereo is stripped."""
        results = qsar_ready_batch(
            [leucine_smiles, r_leucine_smiles], qsar_2d_config
        )

        assert len(results) == 2
        # After stereo stripping, both should have same InChIKey
        assert results[0].standardized_inchikey == results[1].standardized_inchikey
        # First should be ok, second should be duplicate
        statuses = {results[0].status, results[1].status}
        assert statuses == {"ok", "duplicate"}


class TestResultDataclass:
    """Test QSARReadyResult dataclass structure."""

    def test_result_has_original_smiles(self, aspirin_smiles, default_config):
        """Result preserves the original SMILES input."""
        result = qsar_ready_single(aspirin_smiles, default_config)

        assert result.original_smiles == aspirin_smiles

    def test_step_result_has_before_after_smiles(self, aspirin_smiles, default_config):
        """Enabled step results have before_smiles and after_smiles."""
        result = qsar_ready_single(aspirin_smiles, default_config)

        # Step 1 (parse) captures before/after for all enabled steps
        for step in result.steps:
            if step.enabled and step.status in ("applied", "no_change"):
                # Filter step doesn't capture before/after the same way
                if step.step_name not in ("filter", "canonical"):
                    assert step.before_smiles is not None or step.status == "no_change"
