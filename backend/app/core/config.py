"""
Configuration management for ChemStructVal backend.

Uses pydantic-settings for environment variable management.
"""
from pydantic import ConfigDict
from pydantic_settings import BaseSettings
from typing import List


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    # App metadata
    APP_NAME: str = "ChemStructVal"
    APP_VERSION: str = "0.1.0"
    DEBUG: bool = False

    # Database
    DATABASE_URL: str = "postgresql+asyncpg://chemstructval:chemstructval@localhost:5432/chemstructval"

    # Redis
    REDIS_URL: str = "redis://localhost:6379"

    # CORS
    CORS_ORIGINS: List[str] = ["http://localhost:3000", "http://127.0.0.1:3000"]

    # Rate Limiting
    RATE_LIMIT_ENABLED: bool = True

    # Validation limits
    MAX_MOLECULE_LENGTH: int = 10000
    MAX_BATCH_SIZE: int = 10000

    model_config = ConfigDict(
        env_file=".env",
        extra="ignore",
    )


settings = Settings()
