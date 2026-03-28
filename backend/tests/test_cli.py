"""Tests for ChemAudit CLI tool."""

import json
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from app.cli import app

runner = CliRunner()

MOCK_VALIDATION_RESPONSE = {
    "status": "completed",
    "overall_score": 85,
    "molecule_info": {
        "canonical_smiles": "CCO",
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    },
    "issues": [],
}

MOCK_SCORING_RESPONSE = {
    "ml_readiness": {"score": 75, "label": "good"},
    "druglikeness": {
        "lipinski": {"passed": True, "violations": 0},
        "qed": {"score": 0.6},
    },
    "np_likeness": {"score": -1.5, "interpretation": "synthetic-like"},
    "admet": {
        "synthetic_accessibility": {"score": 1.5, "classification": "easy"},
    },
}

MOCK_STANDARDIZATION_RESPONSE = {
    "result": {
        "original_smiles": "CCO",
        "standardized_smiles": "CCO",
        "success": True,
        "steps_applied": [],
    }
}

MOCK_PROFILE_RESPONSE = {
    "properties": {"MW": 46.07, "LogP": -0.31},
    "druglikeness": {"lipinski": {"passed": True}},
    "admet": {"synthetic_accessibility": {"score": 1.5}},
}


class TestValidateCommand:
    """Tests for the validate subcommand."""

    @patch("app.cli._http_post")
    def test_validate_smiles_json_output(self, mock_http):
        """Test validate with --smiles and --format json."""
        mock_http.return_value = MOCK_VALIDATION_RESPONSE
        result = runner.invoke(app, ["validate", "--smiles", "CCO", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["status"] == "completed"
        assert data["overall_score"] == 85

    @patch("app.cli._http_post")
    def test_validate_smiles_table_output(self, mock_http):
        """Test validate with --format table produces Rich table."""
        mock_http.return_value = MOCK_VALIDATION_RESPONSE
        result = runner.invoke(app, ["validate", "--smiles", "CCO", "--format", "table"])
        assert result.exit_code == 0
        assert "Validation Result" in result.output

    @patch("app.cli._http_post")
    def test_validate_http_error(self, mock_http):
        """Test validate handles HTTP errors gracefully."""
        import httpx

        mock_response = MagicMock()
        mock_response.status_code = 500
        mock_response.request = MagicMock()
        mock_http.side_effect = httpx.HTTPStatusError(
            "Server error", request=mock_response.request, response=mock_response
        )
        result = runner.invoke(app, ["validate", "--smiles", "CCO", "--format", "json"])
        assert result.exit_code == 1


class TestScoreCommand:
    """Tests for the score subcommand."""

    @patch("app.cli._http_post")
    def test_score_smiles(self, mock_http):
        """Test score with --smiles and --format json."""
        mock_http.return_value = MOCK_SCORING_RESPONSE
        result = runner.invoke(app, ["score", "--smiles", "CCO", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "ml_readiness" in data or "druglikeness" in data


class TestStandardizeCommand:
    """Tests for the standardize subcommand."""

    @patch("app.cli._http_post")
    def test_standardize_smiles(self, mock_http):
        """Test standardize with --smiles and --format json."""
        mock_http.return_value = MOCK_STANDARDIZATION_RESPONSE
        result = runner.invoke(app, ["standardize", "--smiles", "CCO", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["result"]["original_smiles"] == "CCO"


class TestProfileCommand:
    """Tests for the profile subcommand."""

    @patch("app.cli._http_post")
    def test_profile_smiles(self, mock_http):
        """Test profile with --smiles and --format json."""
        mock_http.return_value = MOCK_PROFILE_RESPONSE
        result = runner.invoke(app, ["profile", "--smiles", "CCO", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "properties" in data or "druglikeness" in data


class TestStdinInput:
    """Tests for stdin pipe input."""

    @patch("app.cli._http_post")
    def test_stdin_input(self, mock_http):
        """Test piping SMILES via stdin."""
        mock_http.return_value = MOCK_VALIDATION_RESPONSE
        result = runner.invoke(app, ["validate", "--format", "json"], input="CCO\n")
        assert result.exit_code == 0
        # Verify the mock was called with payload containing the SMILES
        mock_http.assert_called_once()
        call_args = mock_http.call_args
        assert call_args[0][2]["molecule"] == "CCO"


class TestNoInput:
    """Tests for missing input handling."""

    def test_no_input_shows_error(self):
        """Test that no input produces an error."""
        # Provide empty string to force stdin to be non-TTY but empty
        result = runner.invoke(app, ["validate", "--format", "json"], input="")
        # Should fail because no SMILES provided
        assert result.exit_code != 0 or "Error" in result.output or "Provide" in result.output


class TestServerOption:
    """Tests for custom server URL."""

    @patch("app.cli._http_post")
    def test_server_option(self, mock_http):
        """Test --server option uses custom URL."""
        mock_http.return_value = MOCK_VALIDATION_RESPONSE
        result = runner.invoke(
            app,
            ["validate", "--smiles", "CCO", "--server", "http://custom:9000", "--format", "json"],
        )
        assert result.exit_code == 0
        mock_http.assert_called_once()
        call_args = mock_http.call_args
        assert call_args[0][0] == "http://custom:9000"


class TestHelpOutput:
    """Tests for help and no-args behavior."""

    def test_no_args_shows_help(self):
        """Test that running without arguments shows help text."""
        result = runner.invoke(app, [])
        # Typer no_args_is_help=True may exit with 0 or 2 depending on version
        assert result.exit_code in (0, 2)
        assert "Usage" in result.output or "chemaudit" in result.output.lower()

    def test_validate_help(self):
        """Test validate --help shows usage information."""
        result = runner.invoke(app, ["validate", "--help"])
        assert result.exit_code == 0
        assert "smiles" in result.output.lower()

    def test_score_help(self):
        """Test score --help shows usage information."""
        result = runner.invoke(app, ["score", "--help"])
        assert result.exit_code == 0
        assert "smiles" in result.output.lower()
