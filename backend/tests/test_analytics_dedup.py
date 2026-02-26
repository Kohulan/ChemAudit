"""
Tests for the multi-level deduplication service.

Covers all 4 dedup levels (exact, tautomeric, stereo_insensitive, salt_form),
the all-levels aggregator, and edge cases (empty batch, error results,
unparseable SMILES, representative-index selection).
"""

from __future__ import annotations

from app.services.analytics.deduplication import (
    compute_all_dedup_levels,
    compute_exact_dedup,
    compute_saltform_dedup,
    compute_stereo_dedup,
    compute_tautomer_dedup,
)

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _make_result(
    smiles: str,
    index: int,
    status: str = "success",
    standardized_smiles: str | None = None,
) -> dict:
    """
    Build a minimal batch-result dict that matches the format expected by all
    dedup functions.

    Args:
        smiles: Input SMILES string (may be invalid for edge-case tests).
        index: Original submission index for this molecule.
        status: ``"success"`` | ``"error"`` | other status values.
        standardized_smiles: Optional pre-standardized SMILES stored in the
            ``standardization`` sub-dict (used by salt-form fast path).

    Returns:
        Dict with ``smiles``, ``index``, ``status``, and ``standardization`` keys.
    """
    std: dict = {}
    if standardized_smiles is not None:
        std["standardized_smiles"] = standardized_smiles
    return {"smiles": smiles, "index": index, "status": status, "standardization": std}


# ---------------------------------------------------------------------------
# Exact dedup tests
# ---------------------------------------------------------------------------


class TestExactDedup:
    """Tests for compute_exact_dedup."""

    def test_exact_dedup_finds_duplicates(self):
        """Two identical SMILES strings are grouped as duplicates."""
        results = [
            _make_result("CCO", 0),
            _make_result("CCO", 1),
        ]
        groups = compute_exact_dedup(results)
        assert len(groups) == 1
        g = groups[0]
        assert g["level"] == "exact"
        assert g["count"] == 2
        assert g["representative_index"] == 0
        assert g["duplicate_indices"] == [1]

    def test_exact_dedup_no_duplicates(self):
        """Three structurally distinct molecules produce no dedup groups."""
        results = [
            _make_result("CCO", 0),
            _make_result("CCC", 1),
            _make_result("c1ccccc1", 2),
        ]
        groups = compute_exact_dedup(results)
        assert groups == []

    def test_exact_dedup_canonical_normalization(self):
        """'C(O)C' and 'CCO' are the same molecule after canonicalization."""
        results = [
            _make_result("C(O)C", 0),
            _make_result("CCO", 1),
        ]
        groups = compute_exact_dedup(results)
        assert len(groups) == 1
        g = groups[0]
        assert g["count"] == 2
        # Both should collapse to canonical 'CCO'
        assert g["group_key"] == "CCO"

    def test_exact_dedup_skips_errors(self):
        """Results with status='error' are excluded from dedup."""
        results = [
            _make_result("CCO", 0),
            _make_result("CCO", 1, status="error"),
            _make_result("CCO", 2),
        ]
        groups = compute_exact_dedup(results)
        # Only indices 0 and 2 are valid — still a duplicate pair.
        assert len(groups) == 1
        assert 1 not in groups[0]["duplicate_indices"]
        assert 1 != groups[0]["representative_index"]

    def test_exact_dedup_representative_is_min_index(self):
        """In a group of 3 duplicates with indices [5, 2, 8], representative is 2."""
        results = [
            _make_result("CCO", 5),
            _make_result("CCO", 2),
            _make_result("CCO", 8),
        ]
        groups = compute_exact_dedup(results)
        assert len(groups) == 1
        g = groups[0]
        assert g["representative_index"] == 2
        assert sorted(g["duplicate_indices"]) == [5, 8]
        assert g["count"] == 3


# ---------------------------------------------------------------------------
# Tautomeric dedup tests
# ---------------------------------------------------------------------------


class TestTautomerDedup:
    """Tests for compute_tautomer_dedup."""

    def test_tautomer_dedup_keto_enol_pair(self):
        """Keto form CC(=O)C and enol form CC(O)=C are tautomers — grouped."""
        results = [
            _make_result("CC(=O)C", 0),   # acetone (keto)
            _make_result("CC(O)=C", 1),   # enol form
        ]
        groups = compute_tautomer_dedup(results)
        # Both should canonicalize to the same tautomer
        assert len(groups) == 1
        g = groups[0]
        assert g["level"] == "tautomeric"
        assert g["count"] == 2
        assert g["representative_index"] == 0

    def test_tautomer_dedup_distinct_molecules(self):
        """Two structurally unrelated molecules are not tautomers."""
        results = [
            _make_result("CCO", 0),
            _make_result("c1ccccc1", 1),
        ]
        groups = compute_tautomer_dedup(results)
        assert groups == []


# ---------------------------------------------------------------------------
# Stereo-insensitive dedup tests
# ---------------------------------------------------------------------------


class TestStereoDedup:
    """Tests for compute_stereo_dedup."""

    def test_stereo_dedup_enantiomers_grouped(self):
        """C[C@@H](O)F and C[C@H](O)F are enantiomers — grouped as stereo duplicates."""
        results = [
            _make_result("C[C@@H](O)F", 0),
            _make_result("C[C@H](O)F", 1),
        ]
        groups = compute_stereo_dedup(results)
        assert len(groups) == 1
        g = groups[0]
        assert g["level"] == "stereo_insensitive"
        assert g["count"] == 2
        assert g["representative_index"] == 0

    def test_stereo_dedup_different_molecules_not_grouped(self):
        """CCO and CCCO are constitutionally distinct — NOT grouped."""
        results = [
            _make_result("CCO", 0),
            _make_result("CCCO", 1),
        ]
        groups = compute_stereo_dedup(results)
        assert groups == []

    def test_stereo_dedup_strips_stereo_layers(self):
        """Stereo-stripped InChI key must not contain /t, /m, /s sub-layers."""
        results = [
            _make_result("C[C@@H](O)F", 0),
            _make_result("C[C@H](O)F", 1),
        ]
        groups = compute_stereo_dedup(results)
        assert len(groups) == 1
        # The group key should be the stereo-stripped InChI (no /t /m /s)
        key = groups[0]["group_key"]
        for layer in ["/t", "/m", "/s"]:
            assert layer not in key, f"Stereo layer {layer!r} still present in key: {key}"


# ---------------------------------------------------------------------------
# Salt-form dedup tests
# ---------------------------------------------------------------------------


class TestSaltFormDedup:
    """Tests for compute_saltform_dedup."""

    def test_saltform_dedup_acid_salt_grouped(self):
        """'CCN.Cl' and 'CCN' have the same parent compound (ethylamine)."""
        results = [
            _make_result("CCN.Cl", 0),
            _make_result("CCN", 1),
        ]
        groups = compute_saltform_dedup(results)
        assert len(groups) == 1
        g = groups[0]
        assert g["level"] == "salt_form"
        assert g["count"] == 2
        assert g["representative_index"] == 0
        # Parent of both is ethylamine
        assert g["group_key"] == "CCN"

    def test_saltform_dedup_different_parents_not_grouped(self):
        """Different parent compounds are not grouped."""
        results = [
            _make_result("CCN", 0),
            _make_result("CCCN", 1),
        ]
        groups = compute_saltform_dedup(results)
        assert groups == []

    def test_saltform_dedup_fast_path_standardized_smiles(self):
        """If standardized_smiles is present in the result, it is used as fast path."""
        # Supply already-desalted SMILES in standardization block
        results = [
            _make_result("CCN.[Na+]", 0, standardized_smiles="CCN"),
            _make_result("CCN", 1),
        ]
        groups = compute_saltform_dedup(results)
        # Both should map to CCN parent
        assert len(groups) == 1
        assert groups[0]["representative_index"] == 0


# ---------------------------------------------------------------------------
# Aggregator tests
# ---------------------------------------------------------------------------


class TestComputeAllDedupLevels:
    """Tests for compute_all_dedup_levels."""

    def test_returns_deduplication_result_schema(self):
        """Return value is a DeduplicationResult with all expected keys."""
        from app.schemas.analytics import DeduplicationResult

        results = [_make_result("CCO", 0), _make_result("CCC", 1)]
        out = compute_all_dedup_levels(results)
        assert isinstance(out, DeduplicationResult)
        assert hasattr(out, "exact")
        assert hasattr(out, "tautomeric")
        assert hasattr(out, "stereo_insensitive")
        assert hasattr(out, "salt_form")
        assert hasattr(out, "total_unique")
        assert set(out.total_unique.keys()) == {
            "exact",
            "tautomeric",
            "stereo_insensitive",
            "salt_form",
        }

    def test_compute_all_dedup_levels_mixed_batch(self):
        """Mixed batch: exact pair + stereo pair + unique → all levels report correct groups."""
        results = [
            _make_result("CCO", 0),          # exact dup of index 2
            _make_result("CCC", 1),          # unique
            _make_result("CCO", 2),          # exact dup of index 0
            _make_result("C[C@@H](O)F", 3),  # stereo dup of index 4
            _make_result("C[C@H](O)F", 4),   # stereo dup of index 3
        ]
        out = compute_all_dedup_levels(results)

        # Exact: 1 group (CCO × 2).  The return value is a DeduplicationResult
        # Pydantic model, so list elements are DeduplicationGroup objects —
        # access via attribute, not subscript.
        assert len(out.exact) == 1
        assert out.exact[0].representative_index == 0

        # total_unique for exact: 5 success - 1 removed = 4
        assert out.total_unique["exact"] == 4

        # Stereo: 2 groups — CCO pair (no stereo → same InChI) + enantiomers pair.
        # Both the CCO molecules and the enantiomer pair each form one stereo group,
        # so total stereo groups == 2, total_unique = 5 - 2 = 3.
        assert len(out.stereo_insensitive) == 2
        stereo_reps = {g.representative_index for g in out.stereo_insensitive}
        assert 0 in stereo_reps  # CCO group representative
        assert 3 in stereo_reps  # enantiomers group representative
        assert out.total_unique["stereo_insensitive"] == 3

    def test_compute_all_dedup_levels_total_unique_no_dups(self):
        """No duplicates → total_unique equals total successful molecule count."""
        results = [
            _make_result("CCO", 0),
            _make_result("CCC", 1),
            _make_result("c1ccccc1", 2),
        ]
        out = compute_all_dedup_levels(results)
        for level, count in out.total_unique.items():
            assert count == 3, f"Expected 3 unique at {level!r}, got {count}"


# ---------------------------------------------------------------------------
# Edge-case tests
# ---------------------------------------------------------------------------


class TestEdgeCases:
    """Edge cases: empty batch, all errors, unparseable SMILES."""

    def test_dedup_empty_batch(self):
        """Empty results list → all levels return empty lists, total_unique all zero."""
        out = compute_all_dedup_levels([])
        assert out.exact == []
        assert out.tautomeric == []
        assert out.stereo_insensitive == []
        assert out.salt_form == []
        for level, count in out.total_unique.items():
            assert count == 0, f"Expected 0 at {level!r}, got {count}"

    def test_dedup_unparseable_smiles_no_crash(self):
        """Results with invalid SMILES are skipped gracefully — no exception raised.

        Unparseable molecules are excluded from grouping but still count toward
        ``total_success`` (their status is "success"), so ``total_unique`` equals
        the total successful count (3) with no groups formed.
        """
        results = [
            _make_result("NOT_A_SMILES", 0),
            _make_result("ALSO_INVALID$$", 1),
            _make_result("CCO", 2),
        ]
        # Should complete without error; no duplicate groups since unparseable
        # molecules are skipped and CCO appears only once.
        out = compute_all_dedup_levels(results)
        assert out.exact == []
        assert out.tautomeric == []
        assert out.stereo_insensitive == []
        assert out.salt_form == []
        # total_unique = total_success (3) - 0 duplicates removed = 3
        for count in out.total_unique.values():
            assert count == 3

    def test_dedup_all_errors_skipped(self):
        """All error results → no duplicates detected, total_unique all zero."""
        results = [
            _make_result("CCO", 0, status="error"),
            _make_result("CCO", 1, status="error"),
        ]
        out = compute_all_dedup_levels(results)
        assert out.exact == []
        for count in out.total_unique.values():
            assert count == 0

    def test_dedup_single_molecule_no_groups(self):
        """A batch with one successful molecule produces no dedup groups."""
        results = [_make_result("CCO", 0)]
        out = compute_all_dedup_levels(results)
        assert out.exact == []
        assert out.tautomeric == []
        assert out.stereo_insensitive == []
        assert out.salt_form == []
