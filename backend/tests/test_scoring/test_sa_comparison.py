"""
Tests for SA comparison service: SA Score + SCScore + SYBA side-by-side.

Tests cover:
- SA Score computation and classification (always available via RDKit Contrib)
- SCScore graceful fallback when not installed or weight file fails
- SYBA graceful fallback when subprocess fails or times out
- GPL-3.0 isolation: no module-level syba import in sa_comparison.py
"""

import ast
from pathlib import Path
from unittest.mock import patch

import pytest
from rdkit import Chem

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
SA_COMPARISON_MODULE = Path(__file__).parent.parent.parent / "app" / "services" / "profiler" / "sa_comparison.py"


class TestSAComparison:
    """Tests for SA comparison with graceful fallbacks."""

    @pytest.fixture
    def aspirin_mol(self):
        """Aspirin RDKit molecule."""
        return Chem.MolFromSmiles(ASPIRIN_SMILES)

    # ------------------------------------------------------------------
    # SA Score tests (always available — RDKit Contrib bundled)
    # ------------------------------------------------------------------

    def test_sa_score_aspirin(self, aspirin_mol):
        """SA Score for aspirin should be approximately 1.58 (tolerance +-0.1)."""
        from app.services.profiler.sa_comparison import _compute_sa_score

        result = _compute_sa_score(aspirin_mol)

        assert result["available"] is True
        assert abs(result["score"] - 1.58) < 0.1, f"Expected ~1.58, got {result['score']}"
        assert result["scale"] == "1-10"
        assert "classification" in result

    def test_sa_score_classification_easy(self, aspirin_mol):
        """SA score < 3 returns 'easy' classification."""
        from app.services.profiler.sa_comparison import _compute_sa_score

        result = _compute_sa_score(aspirin_mol)
        # Aspirin SA score is ~1.58 which is < 3 -> "easy"
        assert result["classification"] == "easy"

    def test_sa_score_classification_moderate(self):
        """SA score between 3 and 5 returns 'moderate' classification."""
        from app.services.profiler.sa_comparison import _compute_sa_score

        # Mock sascorer to return a moderate score (lazy loader seam)
        with patch("app.services.profiler.sa_comparison._get_sascorer") as mock_get:
            mock_get.return_value.calculateScore.return_value = 4.0
            mol = Chem.MolFromSmiles("C")
            result = _compute_sa_score(mol)

        assert result["classification"] == "moderate"

    def test_sa_score_classification_difficult(self):
        """SA score >= 5 returns 'difficult' classification."""
        from app.services.profiler.sa_comparison import _compute_sa_score

        # Mock sascorer to return a difficult score (lazy loader seam)
        with patch("app.services.profiler.sa_comparison._get_sascorer") as mock_get:
            mock_get.return_value.calculateScore.return_value = 6.0
            mol = Chem.MolFromSmiles("C")
            result = _compute_sa_score(mol)

        assert result["classification"] == "difficult"

    # ------------------------------------------------------------------
    # SCScore tests (optional dependency — graceful fallback)
    # ------------------------------------------------------------------

    def test_scscore_available_for_aspirin(self, aspirin_mol):
        """With vendored weights present, SCScore is available and matches the model."""
        from app.services.profiler.sa_comparison import _compute_scscore

        result = _compute_scscore(ASPIRIN_SMILES)

        assert result["available"] is True
        assert result["scale"] == "1-5"
        assert result["score"] == 1.59  # reproduces upstream standalone model

    def test_scscore_unavailable_when_weights_missing(self, aspirin_mol):
        """When weights can't be loaded, scscore result must have available=False."""
        from app.services.profiler import sa_comparison

        with patch.object(sa_comparison, "_load_scscore_weights", return_value=None):
            result = sa_comparison._compute_scscore(ASPIRIN_SMILES)

        assert result["available"] is False
        assert "error" in result

    def test_scscore_unavailable_message(self, aspirin_mol):
        """SCScore unavailable result must include descriptive error message."""
        from app.services.profiler import sa_comparison

        with patch.object(sa_comparison, "_load_scscore_weights", return_value=None):
            result = sa_comparison._compute_scscore(ASPIRIN_SMILES)

        assert "scscore" in result["error"].lower() or "not available" in result["error"].lower()

    # ------------------------------------------------------------------
    # SYBA score -> classification mapping (worker mechanism is patched out;
    # the persistent-worker IPC/lifecycle is covered by test_syba_worker.py)
    # ------------------------------------------------------------------

    _SYBA_PREDICT = "app.services.profiler.sa_comparison._syba_via_subprocess"

    def test_syba_unavailable_returns_available_false(self):
        """When the worker can't produce a score, syba result is available=False."""
        from app.services.profiler.sa_comparison import _compute_syba

        with patch(self._SYBA_PREDICT, return_value=None):
            result = _compute_syba(ASPIRIN_SMILES)

        assert result["available"] is False
        assert "error" in result

    def test_syba_success(self):
        """A valid score yields available=True with the score echoed back."""
        from app.services.profiler.sa_comparison import _compute_syba

        with patch(self._SYBA_PREDICT, return_value=42.5):
            result = _compute_syba(ASPIRIN_SMILES)

        assert result["available"] is True
        assert result["score"] == 42.5
        assert result["classification"] == "moderate"  # 0 < 42.5 <= 50

    def test_syba_classification_easy(self):
        from app.services.profiler.sa_comparison import _compute_syba

        with patch(self._SYBA_PREDICT, return_value=75.0):
            result = _compute_syba(ASPIRIN_SMILES)
        assert result["classification"] == "easy"

    def test_syba_classification_moderate(self):
        from app.services.profiler.sa_comparison import _compute_syba

        with patch(self._SYBA_PREDICT, return_value=25.0):
            result = _compute_syba(ASPIRIN_SMILES)
        assert result["classification"] == "moderate"

    def test_syba_classification_difficult(self):
        from app.services.profiler.sa_comparison import _compute_syba

        with patch(self._SYBA_PREDICT, return_value=-50.0):
            result = _compute_syba(ASPIRIN_SMILES)
        assert result["classification"] == "difficult"

    # ------------------------------------------------------------------
    # Full comparison function tests
    # ------------------------------------------------------------------

    def test_full_comparison_returns_all_keys(self, aspirin_mol):
        """compute_sa_comparison must return all 4 keys: sa_score, scscore, syba, rascore."""
        from app.services.profiler.sa_comparison import compute_sa_comparison

        with patch(
            "app.services.profiler.sa_comparison._syba_via_subprocess", return_value=None
        ):
            result = compute_sa_comparison(aspirin_mol, ASPIRIN_SMILES)

        assert "sa_score" in result, "Missing 'sa_score' key"
        assert "scscore" in result, "Missing 'scscore' key"
        assert "syba" in result, "Missing 'syba' key"
        assert "rascore" in result, "Missing 'rascore' key"

    def test_rascore_always_unavailable(self, aspirin_mol):
        """RAscore slot must always have available=False (out of scope for v3.0)."""
        from app.services.profiler.sa_comparison import compute_sa_comparison

        with patch(
            "app.services.profiler.sa_comparison._syba_via_subprocess", return_value=None
        ):
            result = compute_sa_comparison(aspirin_mol, ASPIRIN_SMILES)

        assert result["rascore"]["available"] is False

    def test_sa_score_always_in_result(self, aspirin_mol):
        """SA Score must always be present and available=True even when other scorers fail."""
        from app.services.profiler.sa_comparison import compute_sa_comparison

        with patch(
            "app.services.profiler.sa_comparison._syba_via_subprocess", return_value=None
        ):
            result = compute_sa_comparison(aspirin_mol, ASPIRIN_SMILES)

        assert result["sa_score"]["available"] is True
        assert "score" in result["sa_score"]

    # ------------------------------------------------------------------
    # GPL-3.0 isolation test
    # ------------------------------------------------------------------

    def test_no_syba_module_import(self):
        """sa_comparison.py must NOT have top-level 'import syba' or 'from syba' statements.

        Uses AST parsing to verify GPL-3.0 isolation — SYBA must only be called
        via subprocess, never imported at module level in this Apache-2.0 codebase.
        """
        assert SA_COMPARISON_MODULE.exists(), f"Module not found: {SA_COMPARISON_MODULE}"

        source = SA_COMPARISON_MODULE.read_text()
        tree = ast.parse(source)

        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert not alias.name.startswith("syba"), (
                        f"Found forbidden top-level 'import syba' at line {node.lineno}"
                    )
            elif isinstance(node, ast.ImportFrom):
                if node.module and node.module.startswith("syba"):
                    pytest.fail(
                        f"Found forbidden top-level 'from syba import ...' at line {node.lineno}"
                    )
