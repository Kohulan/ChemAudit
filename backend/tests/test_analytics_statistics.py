"""
Tests for backend/app/services/analytics/statistics.py

Covers:
- Property distribution statistics
- Pairwise Pearson correlations
- IQR-based outlier detection
- Composite batch quality score
- Edge cases (empty batch, insufficient data)
- Integration via compute_all_statistics
"""

from __future__ import annotations

import pytest

from app.services.analytics.statistics import (
    PROPERTY_EXTRACTORS,
    compute_all_statistics,
    compute_correlations,
    compute_outliers,
    compute_property_stats,
    compute_quality_score,
)

# ---------------------------------------------------------------------------
# Helper: build a synthetic result dict
# ---------------------------------------------------------------------------


def _make_result_with_scores(
    smiles: str,
    index: int,
    validation_score: float | None = None,
    qed: float | None = None,
    lipinski_passed: bool | None = None,
    sa_score: float | None = None,
    fsp3: float | None = None,
    ml_readiness_score: float | None = None,
    status: str = "success",
) -> dict:
    """
    Build a minimal molecule result dict with the nested scoring structure
    used by the batch processing pipeline.

    All score arguments are optional; omitting them means the field will
    not appear in the result (simulating a run where that scorer was skipped).
    """
    result: dict = {
        "index": index,
        "smiles": smiles,
        "status": status,
        "validation": {},
        "scoring": {
            "druglikeness": {},
            "admet": {},
            "ml_readiness": {},
        },
    }

    if validation_score is not None:
        result["validation"]["overall_score"] = validation_score

    if qed is not None:
        result["scoring"]["druglikeness"]["qed_score"] = qed

    if lipinski_passed is not None:
        result["scoring"]["druglikeness"]["lipinski_passed"] = lipinski_passed

    if sa_score is not None:
        result["scoring"]["admet"]["sa_score"] = sa_score

    if fsp3 is not None:
        result["scoring"]["admet"]["fsp3"] = fsp3

    if ml_readiness_score is not None:
        result["scoring"]["ml_readiness"]["score"] = ml_readiness_score

    return result


def _make_batch(
    n: int,
    base_qed: float = 0.7,
    base_validation: float = 80.0,
    lipinski_passed: bool = True,
) -> list[dict]:
    """Build a batch of n results with consistent scoring for baseline tests."""
    # Use a selection of simple SMILES that RDKit can parse
    smiles_pool = [
        "CCO",
        "CC(=O)O",
        "c1ccccc1",
        "CC1=CC=CC=C1",
        "CC(C)CC1=CC=CC=C1",
        "CC(C)C1=CC=CC=C1",
        "CCCCO",
        "CCCCN",
        "CC1CCCCC1",
        "c1ccc2ccccc2c1",
    ]
    results = []
    for i in range(n):
        smi = smiles_pool[i % len(smiles_pool)]
        results.append(
            _make_result_with_scores(
                smiles=smi,
                index=i,
                validation_score=base_validation + (i % 10),
                qed=round(base_qed + (i % 5) * 0.02, 3),
                lipinski_passed=lipinski_passed,
                sa_score=2.0 + (i % 3) * 0.1,
                fsp3=0.2 + (i % 4) * 0.05,
                ml_readiness_score=75.0 + (i % 5),
            )
        )
    return results


# ===========================================================================
# Property statistics tests
# ===========================================================================


class TestPropertyStatsBasic:
    """test_property_stats_basic: known QED scores produce expected stats."""

    def test_mean_median_std_within_range(self):
        """Batch of 10 molecules with known QED scores gives expected stats."""
        qed_values = [0.5, 0.6, 0.7, 0.8, 0.9, 0.5, 0.6, 0.7, 0.8, 0.9]
        results = [
            _make_result_with_scores(smiles="CCO", index=i, qed=qed_values[i])
            for i in range(10)
        ]

        stats_list = compute_property_stats(results)
        qed_stat = next((s for s in stats_list if s.property_name == "qed_score"), None)

        assert qed_stat is not None, "qed_score should appear in property stats"
        assert abs(qed_stat.mean - 0.7) < 0.01, f"Expected mean ~0.7, got {qed_stat.mean}"
        assert abs(qed_stat.median - 0.7) < 0.01, f"Expected median ~0.7, got {qed_stat.median}"
        assert qed_stat.std > 0, "Std should be positive for varying values"
        assert qed_stat.min == pytest.approx(0.5, abs=0.001)
        assert qed_stat.max == pytest.approx(0.9, abs=0.001)
        assert qed_stat.count == 10


class TestPropertyStatsAllProperties:
    """test_property_stats_all_properties: all 5 properties produce stats."""

    def test_all_five_properties_present(self):
        """All 5 properties in PROPERTY_EXTRACTORS should appear when data available."""
        results = _make_batch(10)

        stats_list = compute_property_stats(results)
        returned_names = {s.property_name for s in stats_list}

        # All extractors should have produced stats
        assert returned_names == set(PROPERTY_EXTRACTORS.keys()), (
            f"Expected {set(PROPERTY_EXTRACTORS.keys())}, got {returned_names}"
        )


class TestPropertyStatsSkipInsufficient:
    """test_property_stats_skip_insufficient: single value for property is skipped."""

    def test_single_value_property_skipped(self):
        """Only 1 value for qed_score; that property should be omitted."""
        # Only one result has a QED score
        results = [
            _make_result_with_scores(smiles="CCO", index=0, qed=0.8, validation_score=80.0),
            _make_result_with_scores(
                smiles="CCO", index=1, qed=None, validation_score=70.0
            ),
            _make_result_with_scores(
                smiles="CCO", index=2, qed=None, validation_score=75.0
            ),
        ]

        stats_list = compute_property_stats(results)
        names = {s.property_name for s in stats_list}
        assert "qed_score" not in names, "qed_score with only 1 value should be skipped"
        # validation_score has 3 values, should be present
        assert "validation_score" in names


# ===========================================================================
# Correlation tests
# ===========================================================================


class TestCorrelationsPositive:
    """test_correlations_positive: correlated properties yield positive pearson_r."""

    def test_positively_correlated_scores(self):
        """Build batch where validation_score and ml_readiness_score are linearly correlated."""
        results = []
        # Create 15 results where both scores increase together
        for i in range(15):
            results.append(
                _make_result_with_scores(
                    smiles="CCO",
                    index=i,
                    validation_score=50.0 + i * 3.0,
                    ml_readiness_score=50.0 + i * 2.5,
                    qed=0.5 + i * 0.01,
                    sa_score=2.0 + i * 0.05,
                    fsp3=0.2 + i * 0.01,
                )
            )

        correlations = compute_correlations(results)
        pair = next(
            (
                c
                for c in correlations
                if set([c.property_a, c.property_b])
                == {"validation_score", "ml_readiness_score"}
            ),
            None,
        )
        assert pair is not None, "validation_score/ml_readiness_score pair not found"
        assert pair.pearson_r > 0.9, (
            f"Expected strong positive correlation, got {pair.pearson_r}"
        )


class TestCorrelationsNanHandling:
    """test_correlations_nan_handling: constant property yields near-zero correlation, excluded."""

    def test_constant_qed_excluded_from_correlations(self):
        """All QED values constant (0.7) => corrcoef returns near-zero or NaN => pair excluded.

        The service converts NaN to 0.0 and skips |pearson_r| < 1e-10 pairs (floating-point
        noise from constant arrays). This ensures no NaN propagates into results.
        """
        import math as _math

        # 15 results where QED is constant (all 0.7) and other scores vary
        results = [
            _make_result_with_scores(
                smiles="CCO",
                index=i,
                validation_score=float(50 + i * 2),
                qed=0.7,  # constant — produces near-zero corrcoef
                sa_score=2.0 + i * 0.1,
                fsp3=0.2 + i * 0.02,
                ml_readiness_score=50.0 + i * 1.5,
            )
            for i in range(15)
        ]

        correlations = compute_correlations(results)

        # Verify no NaN values in any returned correlation
        for corr in correlations:
            assert not _math.isnan(corr.pearson_r), f"NaN pearson_r found: {corr}"

        # Constant QED pairs should be excluded (near-zero pearson_r filtered out)
        qed_pairs = [
            c for c in correlations
            if c.property_a == "qed_score" or c.property_b == "qed_score"
        ]
        assert len(qed_pairs) == 0, (
            f"Constant QED pairs should be excluded from correlations, got: {qed_pairs}"
        )


class TestCorrelationsInsufficientData:
    """test_correlations_insufficient_data: <10 points for a property excluded."""

    def test_property_with_8_values_excluded(self):
        """Property with only 8 data points should be excluded from correlations."""
        results = []
        # 15 results total — QED only present in first 8
        for i in range(15):
            qed = 0.5 + i * 0.02 if i < 8 else None
            results.append(
                _make_result_with_scores(
                    smiles="CCO",
                    index=i,
                    validation_score=float(60 + i),
                    qed=qed,
                    ml_readiness_score=60.0 + i * 0.5,
                    sa_score=2.0,
                    fsp3=0.3,
                )
            )

        correlations = compute_correlations(results)
        for corr in correlations:
            assert "qed_score" not in (corr.property_a, corr.property_b), (
                "qed_score with 8 data points should be excluded from correlations"
            )


# ===========================================================================
# Outlier detection tests
# ===========================================================================


class TestOutlierDetectionBasic:
    """test_outlier_detection_basic: one extreme value is flagged."""

    def test_single_extreme_outlier_flagged(self):
        """Batch with one very low QED (0.01) when others are 0.5–0.9 should be flagged."""
        results = []
        # 9 normal molecules
        for i in range(9):
            results.append(
                _make_result_with_scores(
                    smiles="CCO",
                    index=i,
                    qed=0.5 + i * 0.04,
                )
            )
        # 1 extreme outlier
        results.append(
            _make_result_with_scores(smiles="CCO", index=9, qed=0.01)
        )

        outliers = compute_outliers(results)
        qed_outliers = [o for o in outliers if o.property_name == "qed_score"]
        assert len(qed_outliers) >= 1, "Expected at least 1 QED outlier"
        assert any(o.molecule_index == 9 for o in qed_outliers), (
            "The extreme QED=0.01 molecule should be flagged"
        )


class TestOutlierDetectionNoOutliers:
    """test_outlier_detection_no_outliers: normally distributed values yield no outliers."""

    def test_uniform_values_no_outliers(self):
        """Closely clustered QED values should not produce outliers."""
        # Narrow range: 0.60 to 0.75 (no extreme values)
        results = [
            _make_result_with_scores(smiles="CCO", index=i, qed=0.60 + i * 0.01)
            for i in range(16)
        ]

        outliers = compute_outliers(results)
        qed_outliers = [o for o in outliers if o.property_name == "qed_score"]
        assert len(qed_outliers) == 0, (
            f"Expected no QED outliers in uniform data, got {qed_outliers}"
        )


class TestOutlierDetectionInsufficient:
    """test_outlier_detection_insufficient: fewer than 4 values → no outliers computed."""

    def test_three_values_no_outliers_computed(self):
        """Only 3 QED values available; outlier detection should be skipped for this property."""
        results = [
            _make_result_with_scores(smiles="CCO", index=0, qed=0.1),
            _make_result_with_scores(smiles="CCO", index=1, qed=0.5),
            _make_result_with_scores(smiles="CCO", index=2, qed=0.9),
        ]

        outliers = compute_outliers(results)
        qed_outliers = [o for o in outliers if o.property_name == "qed_score"]
        assert len(qed_outliers) == 0, "Fewer than 4 values should skip outlier detection"


# ===========================================================================
# Quality score tests
# ===========================================================================


class TestQualityScorePerfectBatch:
    """test_quality_score_perfect_batch: all valid, diverse, Lipinski → score near 100."""

    def test_perfect_batch_high_score(self):
        """All molecules valid, diverse scaffolds, all pass Lipinski → score should be high."""
        # Use structurally diverse SMILES with different scaffolds
        diverse_smiles = [
            "CCO",           # ethanol
            "c1ccccc1",      # benzene
            "C1CCCCC1",      # cyclohexane
            "c1ccncc1",      # pyridine
            "CC(=O)O",       # acetic acid
            "CCCCN",         # butylamine
            "c1ccc2ccccc2c1",  # naphthalene
            "CC1CCCCC1",     # methylcyclohexane
        ]
        results = [
            _make_result_with_scores(
                smiles=smi,
                index=i,
                validation_score=100.0,
                qed=0.9,
                lipinski_passed=True,
                sa_score=2.0,
                fsp3=0.3,
                ml_readiness_score=100.0,
                status="success",
            )
            for i, smi in enumerate(diverse_smiles)
        ]

        qs = compute_quality_score(results)

        # Validity: 100% → contributes 40
        assert qs.validity_pct == pytest.approx(100.0)
        # Drug-likeness: 100% → contributes 25
        assert qs.druglikeness_pct == pytest.approx(100.0)
        # Overall score should be high (>= 65, validity + drug-likeness alone = 65)
        assert qs.score >= 65.0, f"Expected score >= 65, got {qs.score}"


class TestQualityScorePoorBatch:
    """test_quality_score_poor_batch: many errors, same scaffold, low Lipinski → score < 50."""

    def test_poor_batch_low_score(self):
        """50% errors, single scaffold, low Lipinski pass rate → score well below 50."""
        results = []
        # 5 error molecules
        for i in range(5):
            results.append(
                _make_result_with_scores(
                    smiles="CCO", index=i, status="error", lipinski_passed=None
                )
            )
        # 5 success molecules, all same scaffold, all fail Lipinski
        for i in range(5, 10):
            results.append(
                _make_result_with_scores(
                    smiles="c1ccccc1",
                    index=i,
                    status="success",
                    qed=0.3,
                    lipinski_passed=False,
                )
            )

        qs = compute_quality_score(results)
        assert qs.score < 50.0, f"Expected poor quality score < 50, got {qs.score}"


class TestQualityScoreValidityWeight:
    """test_quality_score_validity_weight: 50% errors → validity = 50, contributes 20."""

    def test_validity_component_weight(self):
        """50% error rate → validity_pct=50 → validity contributes exactly 0.4*50=20."""
        results = []
        # 5 success, 5 error
        for i in range(5):
            results.append(
                _make_result_with_scores(
                    smiles="c1ccccc1", index=i, status="success", lipinski_passed=None
                )
            )
        for i in range(5, 10):
            results.append(
                _make_result_with_scores(
                    smiles="CCO", index=i, status="error", lipinski_passed=None
                )
            )

        qs = compute_quality_score(results)
        assert qs.validity_pct == pytest.approx(50.0), (
            f"Expected validity_pct=50.0, got {qs.validity_pct}"
        )
        # Validity contribution to composite should be 20
        validity_contribution = qs.validity_pct * 0.40
        assert validity_contribution == pytest.approx(20.0, abs=0.1), (
            f"Expected validity contribution 20.0, got {validity_contribution}"
        )


class TestQualityScoreDiversitySingleScaffold:
    """test_quality_score_diversity_single_scaffold: all same scaffold → diversity = 0."""

    def test_single_scaffold_zero_diversity(self):
        """All molecules with identical SMILES → single scaffold → diversity_pct = 0."""
        results = [
            _make_result_with_scores(
                smiles="c1ccccc1",
                index=i,
                status="success",
                qed=0.7,
                lipinski_passed=True,
            )
            for i in range(8)
        ]

        qs = compute_quality_score(results)
        assert qs.diversity_pct == pytest.approx(0.0), (
            f"Expected diversity_pct=0.0 for single scaffold, got {qs.diversity_pct}"
        )
        # Diversity contribution should be 0 to composite
        diversity_contribution = qs.diversity_pct * 0.35
        assert diversity_contribution == pytest.approx(0.0), (
            f"Expected diversity contribution 0.0, got {diversity_contribution}"
        )


class TestQualityScoreNoLipinskiData:
    """test_quality_score_no_lipinski_data: no Lipinski data → default 50%."""

    def test_fallback_to_50_percent(self):
        """When no Lipinski data present, drug-likeness defaults to 50%."""
        results = [
            _make_result_with_scores(
                smiles="CCO",
                index=i,
                status="success",
                qed=None,
                lipinski_passed=None,  # No Lipinski data
                validation_score=80.0,
            )
            for i in range(5)
        ]

        qs = compute_quality_score(results)
        assert qs.druglikeness_pct == pytest.approx(50.0), (
            f"Expected druglikeness_pct=50.0 (fallback), got {qs.druglikeness_pct}"
        )


# ===========================================================================
# Integration and edge case tests
# ===========================================================================


class TestComputeAllStatistics:
    """test_compute_all_statistics: full integration returns dict with all 4 keys."""

    def test_returns_all_four_keys(self):
        """compute_all_statistics should return StatisticsResult with all 4 sub-results."""
        results = _make_batch(15)

        stats_result = compute_all_statistics(results)

        assert hasattr(stats_result, "property_stats"), "Missing property_stats"
        assert hasattr(stats_result, "correlations"), "Missing correlations"
        assert hasattr(stats_result, "outliers"), "Missing outliers (list)"
        assert hasattr(stats_result, "quality_score"), "Missing quality_score"

    def test_model_dump_produces_serializable_dict(self):
        """model_dump() (used by analytics_tasks.py) should produce a plain dict."""
        results = _make_batch(10)
        stats_result = compute_all_statistics(results)
        dumped = stats_result.model_dump()

        assert isinstance(dumped, dict), "model_dump() should return a dict"
        assert "property_stats" in dumped
        assert "correlations" in dumped
        assert "outliers" in dumped
        assert "quality_score" in dumped


class TestEmptyBatch:
    """test_empty_batch: empty results → valid empty response (no crash)."""

    def test_empty_property_stats(self):
        """Empty batch produces empty property_stats list."""
        stats_list = compute_property_stats([])
        assert stats_list == []

    def test_empty_correlations(self):
        """Empty batch produces empty correlations list."""
        correlations = compute_correlations([])
        assert correlations == []

    def test_empty_outliers(self):
        """Empty batch produces empty outliers list."""
        outliers = compute_outliers([])
        assert outliers == []

    def test_empty_quality_score(self):
        """Empty batch produces quality score with all zeros."""
        qs = compute_quality_score([])
        assert qs.score == pytest.approx(0.0)
        assert qs.validity_pct == pytest.approx(0.0)
        assert qs.diversity_pct == pytest.approx(0.0)
        assert qs.druglikeness_pct == pytest.approx(0.0)

    def test_compute_all_statistics_empty(self):
        """compute_all_statistics on empty batch returns valid StatisticsResult."""
        result = compute_all_statistics([])
        assert result.property_stats == []
        assert result.correlations == []
        assert result.outliers == []
        assert result.quality_score.score == pytest.approx(0.0)
