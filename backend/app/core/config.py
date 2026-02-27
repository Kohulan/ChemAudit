"""
Configuration management for ChemAudit backend.

Uses pydantic-settings for environment variable management.
"""

import logging
from typing import List

from pydantic import ConfigDict, computed_field, model_validator
from pydantic_settings import BaseSettings

_config_logger = logging.getLogger(__name__)

# Insecure default values that must be replaced in production
_INSECURE_DEFAULTS = frozenset({"CHANGE_ME_IN_PRODUCTION", ""})


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    # App metadata
    APP_NAME: str = "ChemAudit"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = False

    # Database — placeholder; override via DATABASE_URL env var in all environments
    DATABASE_URL: str = "postgresql+asyncpg://user:pass@localhost:5432/chemaudit"

    # Redis
    REDIS_URL: str = "redis://localhost:6379"

    # ==========================================================================
    # Security Settings
    # ==========================================================================
    # Secret key for signing tokens (REQUIRED in production)
    SECRET_KEY: str = "CHANGE_ME_IN_PRODUCTION"

    # Admin secret for API key management endpoints (REQUIRED)
    API_KEY_ADMIN_SECRET: str = "CHANGE_ME_IN_PRODUCTION"

    # CSRF protection secret
    CSRF_SECRET_KEY: str = "CHANGE_ME_IN_PRODUCTION"

    # API Key expiration settings
    API_KEY_DEFAULT_EXPIRY_DAYS: int = 90
    API_KEY_MAX_EXPIRY_DAYS: int = 365

    # ==========================================================================
    # CORS Settings
    # ==========================================================================
    # Store as string, parse to list via computed_field
    CORS_ORIGINS_STR: str = "http://localhost:3002,http://127.0.0.1:3002"

    @computed_field
    @property
    def CORS_ORIGINS(self) -> List[str]:  # noqa: N802
        """Parse CORS origins from comma-separated string."""
        return [
            origin.strip()
            for origin in self.CORS_ORIGINS_STR.split(",")
            if origin.strip()
        ]

    # Explicit CORS methods (not wildcard for security)
    CORS_ALLOW_METHODS: List[str] = ["GET", "POST", "PUT", "DELETE", "OPTIONS", "PATCH"]

    # Explicit CORS headers (not wildcard for security)
    CORS_ALLOW_HEADERS: List[str] = [
        "Accept",
        "Accept-Language",
        "Content-Language",
        "Content-Type",
        "Authorization",
        "X-API-Key",
        "X-CSRF-Token",
        "X-Requested-With",
    ]

    # ==========================================================================
    # Rate Limiting
    # ==========================================================================
    RATE_LIMIT_ENABLED: bool = True

    # IP banning for repeated rate limit violations
    RATE_LIMIT_BAN_THRESHOLD: int = 100  # violations before ban
    RATE_LIMIT_BAN_DURATION_MINUTES: int = 60  # ban duration

    # ==========================================================================
    # Caching
    # ==========================================================================
    VALIDATION_CACHE_TTL: int = 3600  # 1 hour in seconds
    VALIDATION_CACHE_ENABLED: bool = True
    # Batch result TTL — 24h to support analytics (INFRA-01)
    BATCH_RESULT_TTL: int = 86400

    # ==========================================================================
    # Metrics
    # ==========================================================================
    ENABLE_METRICS: bool = False  # Enable Prometheus metrics endpoint

    # ==========================================================================
    # Validation limits
    # ==========================================================================
    MAX_MOLECULE_LENGTH: int = 10000
    MAX_BATCH_SIZE: int = 10000
    MAX_FILE_SIZE_MB: int = 500

    # Deployment profile
    DEPLOYMENT_PROFILE: str = "medium"

    # ==========================================================================
    # SMTP Settings (for batch completion email notifications)
    # ==========================================================================
    SMTP_HOST: str = "localhost"
    SMTP_PORT: int = 587
    SMTP_USER: str = ""
    SMTP_PASSWORD: str = ""
    SMTP_FROM: str = "noreply@chemaudit.local"
    SMTP_TLS: bool = True
    NOTIFICATION_EMAIL: str = ""  # Default recipient for batch completion emails (empty = disabled)

    # ==========================================================================
    # Webhook Settings (for batch completion notifications)
    # ==========================================================================
    WEBHOOK_URL: str = ""  # Endpoint to POST on batch completion (empty = disabled)
    WEBHOOK_SECRET: str = ""  # HMAC signing secret for webhook payloads

    # ==========================================================================
    # OPSIN Settings (IUPAC name to SMILES conversion)
    # ==========================================================================
    OPSIN_JAR_PATH: str = "/app/data/opsin.jar"

    # ==========================================================================
    # Application Base URL (for webhook payloads and email links)
    # ==========================================================================
    BASE_URL: str = "http://localhost:3002"

    # ==========================================================================
    # External API endpoints
    # ==========================================================================
    COCONUT_API_URL: str = "https://coconut.naturalproducts.net/api"
    COCONUT_API_TOKEN: str = ""  # Required for COCONUT API v2
    PUBCHEM_API_URL: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    CHEMBL_API_URL: str = "https://www.ebi.ac.uk/chembl/api/data"
    EXTERNAL_API_TIMEOUT: int = 30  # seconds

    model_config = ConfigDict(
        env_file=".env",
        extra="ignore",
    )

    @model_validator(mode="after")
    def _check_insecure_defaults(self) -> "Settings":
        """Reject insecure default secrets when DEBUG is False."""
        secret_fields = ("SECRET_KEY", "API_KEY_ADMIN_SECRET", "CSRF_SECRET_KEY")
        for field_name in secret_fields:
            value = getattr(self, field_name)
            if value in _INSECURE_DEFAULTS:
                if not self.DEBUG:
                    raise ValueError(
                        f"{field_name} still has an insecure default value. "
                        f"Set a strong secret via environment variable before "
                        f"running in production (DEBUG=False)."
                    )
                _config_logger.warning(
                    "%s has an insecure default value. "
                    "This is acceptable for local development only.",
                    field_name,
                )
        return self


settings = Settings()
