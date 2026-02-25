"""
ValidationAuditEntry ORM Model

Append-only audit trail for all validation events.
"""

from datetime import datetime

from sqlalchemy import DateTime, Integer, String, Text, func
from sqlalchemy.orm import Mapped, mapped_column

from app.db import Base


class ValidationAuditEntry(Base):
    """Append-only validation audit trail entry."""

    __tablename__ = "validation_audit"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    smiles: Mapped[str] = mapped_column(Text, nullable=False)
    inchikey: Mapped[str | None] = mapped_column(String(27), nullable=True)
    outcome: Mapped[str] = mapped_column(String(10), nullable=False)
    score: Mapped[int | None] = mapped_column(Integer, nullable=True)
    job_id: Mapped[str | None] = mapped_column(String(50), nullable=True)
    molecule_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    pass_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    fail_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    api_key_hash: Mapped[str | None] = mapped_column(String(64), nullable=True)
    session_id: Mapped[str | None] = mapped_column(
        String(36), nullable=True, index=True
    )
    source: Mapped[str] = mapped_column(String(10), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime, server_default=func.now(), index=True
    )
