"""
BatchPermalink ORM Model

Shareable permalink for batch validation results.
"""

from datetime import datetime

from sqlalchemy import DateTime, Integer, String, Text, func
from sqlalchemy.orm import Mapped, mapped_column

from app.db import Base


class BatchPermalink(Base):
    """Shareable permalink with batch snapshot data."""

    __tablename__ = "batch_permalinks"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    short_id: Mapped[str] = mapped_column(
        String(16), unique=True, index=True, nullable=False
    )
    job_id: Mapped[str] = mapped_column(String(50), nullable=False)
    snapshot_data: Mapped[str | None] = mapped_column(Text, nullable=True)
    settings: Mapped[str | None] = mapped_column(Text, nullable=True)
    created_at: Mapped[datetime] = mapped_column(
        DateTime, server_default=func.now()
    )
    expires_at: Mapped[datetime | None] = mapped_column(DateTime, nullable=True)
