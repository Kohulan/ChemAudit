"""Test desirability and profile scoring functions."""
import pytest

from app.services.scoring.profile_scoring import desirability, profile_score


class TestDesirability:
    def test_value_in_range_returns_one(self):
        assert desirability(300, 0, 500) == 1.0

    def test_value_at_boundary_returns_one(self):
        assert desirability(500, 0, 500) == 1.0
        assert desirability(0, 0, 500) == 1.0

    def test_value_above_max_linear_falloff(self):
        # 50% beyond range width -> d = 0.5
        result = desirability(750, 0, 500)
        assert result == pytest.approx(0.5, abs=0.01)

    def test_value_below_min_linear_falloff(self):
        # 100 below min of 200, range width 300 -> distance = 100/300 ~ 0.33
        result = desirability(100, 200, 500)
        assert result == pytest.approx(0.667, abs=0.01)

    def test_value_far_outside_returns_zero(self):
        # > 100% beyond range -> clamped to 0
        assert desirability(1100, 0, 500) == 0.0

    def test_none_value_returns_none(self):
        assert desirability(None, 0, 500) is None


class TestProfileScore:
    def test_all_in_range_returns_100(self):
        properties = {"mw": 300, "logp": 2.5}
        thresholds = {"mw": {"min": 0, "max": 500}, "logp": {"min": -5, "max": 5}}
        weights = {"mw": 1.0, "logp": 1.0}
        assert profile_score(properties, thresholds, weights) == 100.0

    def test_one_out_of_range_reduces_score(self):
        properties = {"mw": 300, "logp": 7.5}
        thresholds = {"mw": {"min": 0, "max": 500}, "logp": {"min": -5, "max": 5}}
        weights = {"mw": 1.0, "logp": 1.0}
        score = profile_score(properties, thresholds, weights)
        assert 0 < score < 100

    def test_missing_property_skipped(self):
        properties = {"mw": 300}
        thresholds = {"mw": {"min": 0, "max": 500}, "logp": {"min": -5, "max": 5}}
        weights = {"mw": 1.0, "logp": 1.0}
        # Only mw evaluated, in range -> 100
        assert profile_score(properties, thresholds, weights) == 100.0

    def test_empty_thresholds_returns_none(self):
        assert profile_score({"mw": 300}, {}, {}) is None

    def test_weights_affect_score(self):
        properties = {"mw": 600, "logp": 2.5}
        thresholds = {"mw": {"min": 0, "max": 500}, "logp": {"min": -5, "max": 5}}
        # Heavy weight on mw (which is out of range) -> lower score
        score_heavy = profile_score(properties, thresholds, {"mw": 3.0, "logp": 1.0})
        score_light = profile_score(properties, thresholds, {"mw": 1.0, "logp": 1.0})
        assert score_heavy < score_light
