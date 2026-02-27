"""
Tests for Phase 1: Startup validation of insecure default secrets.

Verifies that Settings rejects insecure defaults in production mode
and only warns in debug mode.
"""

import os
from unittest.mock import patch

import pytest
from pydantic import ValidationError


class TestSettingsSecurityValidation:
    """Test that insecure default secrets are rejected in production."""

    def test_rejects_default_secret_key_in_production(self):
        """SECRET_KEY='CHANGE_ME_IN_PRODUCTION' must fail when DEBUG=False."""
        env = {
            "DEBUG": "false",
            "SECRET_KEY": "CHANGE_ME_IN_PRODUCTION",
            "API_KEY_ADMIN_SECRET": "a-real-secret-here",
            "CSRF_SECRET_KEY": "another-real-secret",
        }
        with patch.dict(os.environ, env, clear=False):
            from importlib import reload

            import app.core.config as config_module

            with pytest.raises(ValidationError, match="SECRET_KEY"):
                reload(config_module)
            # Restore working settings for other tests
            os.environ["DEBUG"] = "true"
            reload(config_module)

    def test_rejects_default_admin_secret_in_production(self):
        """API_KEY_ADMIN_SECRET='CHANGE_ME_IN_PRODUCTION' must fail when DEBUG=False."""
        env = {
            "DEBUG": "false",
            "SECRET_KEY": "a-real-secret-here",
            "API_KEY_ADMIN_SECRET": "CHANGE_ME_IN_PRODUCTION",
            "CSRF_SECRET_KEY": "another-real-secret",
        }
        with patch.dict(os.environ, env, clear=False):
            from importlib import reload

            import app.core.config as config_module

            with pytest.raises(ValidationError, match="API_KEY_ADMIN_SECRET"):
                reload(config_module)
            os.environ["DEBUG"] = "true"
            reload(config_module)

    def test_rejects_default_csrf_secret_in_production(self):
        """CSRF_SECRET_KEY='CHANGE_ME_IN_PRODUCTION' must fail when DEBUG=False."""
        env = {
            "DEBUG": "false",
            "SECRET_KEY": "a-real-secret-here",
            "API_KEY_ADMIN_SECRET": "another-real-secret",
            "CSRF_SECRET_KEY": "CHANGE_ME_IN_PRODUCTION",
        }
        with patch.dict(os.environ, env, clear=False):
            from importlib import reload

            import app.core.config as config_module

            with pytest.raises(ValidationError, match="CSRF_SECRET_KEY"):
                reload(config_module)
            os.environ["DEBUG"] = "true"
            reload(config_module)

    def test_allows_defaults_in_debug_mode(self):
        """Insecure defaults should only warn (not fail) when DEBUG=True."""
        env = {
            "DEBUG": "true",
            "SECRET_KEY": "CHANGE_ME_IN_PRODUCTION",
            "API_KEY_ADMIN_SECRET": "CHANGE_ME_IN_PRODUCTION",
            "CSRF_SECRET_KEY": "CHANGE_ME_IN_PRODUCTION",
        }
        with patch.dict(os.environ, env, clear=False):
            from importlib import reload

            import app.core.config as config_module

            # Should NOT raise — just warn
            reload(config_module)
            assert config_module.settings.DEBUG is True

    def test_allows_strong_secrets_in_production(self):
        """Real secrets should pass validation in production mode."""
        env = {
            "DEBUG": "false",
            "SECRET_KEY": "super-secret-production-key-12345",
            "API_KEY_ADMIN_SECRET": "admin-secret-production-67890",
            "CSRF_SECRET_KEY": "csrf-secret-production-abcde",
        }
        with patch.dict(os.environ, env, clear=False):
            from importlib import reload

            import app.core.config as config_module

            reload(config_module)
            assert config_module.settings.DEBUG is False
            # Restore debug mode for other tests
            os.environ["DEBUG"] = "true"
            reload(config_module)
