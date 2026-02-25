"""Add session_id columns and enable row-level security.

Revision ID: 002
Revises: 001
Create Date: 2026-02-25
"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op

revision: str = "002"
down_revision: Union[str, None] = "001"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # --- 1. Add session_id columns ---
    op.add_column("bookmarks", sa.Column("session_id", sa.String(36), nullable=True))
    op.create_index("ix_bookmarks_session_id", "bookmarks", ["session_id"])

    op.add_column("validation_audit", sa.Column("session_id", sa.String(36), nullable=True))
    op.create_index("ix_validation_audit_session_id", "validation_audit", ["session_id"])

    # --- 2. Purge unscoped legacy data ---
    op.execute("DELETE FROM bookmarks WHERE session_id IS NULL AND api_key_hash IS NULL")
    op.execute(
        "DELETE FROM validation_audit WHERE session_id IS NULL AND api_key_hash IS NULL"
    )

    # --- 3. Enable Row-Level Security ---
    op.execute("ALTER TABLE bookmarks ENABLE ROW LEVEL SECURITY")
    op.execute("ALTER TABLE bookmarks FORCE ROW LEVEL SECURITY")
    op.execute("""
        CREATE POLICY bookmark_isolation ON bookmarks
        USING (
            session_id = NULLIF(current_setting('app.session_id', true), '')
            OR api_key_hash = NULLIF(current_setting('app.api_key_hash', true), '')
        )
    """)

    op.execute("ALTER TABLE validation_audit ENABLE ROW LEVEL SECURITY")
    op.execute("ALTER TABLE validation_audit FORCE ROW LEVEL SECURITY")
    op.execute("""
        CREATE POLICY audit_isolation ON validation_audit
        USING (
            session_id = NULLIF(current_setting('app.session_id', true), '')
            OR api_key_hash = NULLIF(current_setting('app.api_key_hash', true), '')
        )
    """)


def downgrade() -> None:
    # Remove RLS
    op.execute("DROP POLICY IF EXISTS audit_isolation ON validation_audit")
    op.execute("ALTER TABLE validation_audit DISABLE ROW LEVEL SECURITY")

    op.execute("DROP POLICY IF EXISTS bookmark_isolation ON bookmarks")
    op.execute("ALTER TABLE bookmarks DISABLE ROW LEVEL SECURITY")

    # Remove columns
    op.drop_index("ix_validation_audit_session_id", "validation_audit")
    op.drop_column("validation_audit", "session_id")
    op.drop_index("ix_bookmarks_session_id", "bookmarks")
    op.drop_column("bookmarks", "session_id")
