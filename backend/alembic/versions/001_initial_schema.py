"""Initial schema: scoring_profiles, bookmarks, batch_permalinks, validation_audit

Revision ID: 001
Revises:
Create Date: 2026-02-24
"""

from typing import Sequence, Union

import sqlalchemy as sa

from alembic import op

revision: str = "001"
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Create all Phase 6 tables."""
    # scoring_profiles
    op.create_table(
        "scoring_profiles",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column("name", sa.String(100), nullable=False),
        sa.Column("description", sa.Text, nullable=True),
        sa.Column("thresholds", sa.Text, nullable=False, server_default="{}"),
        sa.Column("weights", sa.Text, nullable=False, server_default="{}"),
        sa.Column("is_preset", sa.Boolean, server_default="false"),
        sa.Column("is_active", sa.Boolean, server_default="true"),
        sa.Column("created_at", sa.DateTime, server_default=sa.func.now()),
        sa.Column("updated_at", sa.DateTime, server_default=sa.func.now()),
    )

    # bookmarks
    op.create_table(
        "bookmarks",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column("smiles", sa.Text, nullable=False),
        sa.Column("name", sa.String(200), nullable=True),
        sa.Column("inchikey", sa.String(27), nullable=True),
        sa.Column("tags", sa.Text, nullable=True),
        sa.Column("notes", sa.Text, nullable=True),
        sa.Column("source", sa.String(50), nullable=True),
        sa.Column("job_id", sa.String(50), nullable=True),
        sa.Column("api_key_hash", sa.String(64), nullable=True),
        sa.Column("created_at", sa.DateTime, server_default=sa.func.now()),
    )
    op.create_index("ix_bookmarks_api_key_hash", "bookmarks", ["api_key_hash"])

    # batch_permalinks
    op.create_table(
        "batch_permalinks",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column("short_id", sa.String(16), unique=True, nullable=False),
        sa.Column("job_id", sa.String(50), nullable=False),
        sa.Column("snapshot_data", sa.Text, nullable=True),
        sa.Column("settings", sa.Text, nullable=True),
        sa.Column("created_at", sa.DateTime, server_default=sa.func.now()),
        sa.Column("expires_at", sa.DateTime, nullable=True),
    )
    op.create_index("ix_batch_permalinks_short_id", "batch_permalinks", ["short_id"])

    # validation_audit
    op.create_table(
        "validation_audit",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column("smiles", sa.Text, nullable=False),
        sa.Column("inchikey", sa.String(27), nullable=True),
        sa.Column("outcome", sa.String(10), nullable=False),
        sa.Column("score", sa.Integer, nullable=True),
        sa.Column("job_id", sa.String(50), nullable=True),
        sa.Column("molecule_count", sa.Integer, nullable=True),
        sa.Column("pass_count", sa.Integer, nullable=True),
        sa.Column("fail_count", sa.Integer, nullable=True),
        sa.Column("api_key_hash", sa.String(64), nullable=True),
        sa.Column("source", sa.String(10), nullable=False),
        sa.Column("created_at", sa.DateTime, server_default=sa.func.now()),
    )
    op.create_index("ix_validation_audit_created_at", "validation_audit", ["created_at"])


def downgrade() -> None:
    """Drop all Phase 6 tables."""
    op.drop_table("validation_audit")
    op.drop_table("batch_permalinks")
    op.drop_table("bookmarks")
    op.drop_table("scoring_profiles")
