"""Add owner columns to scoring_profiles.

Revision ID: 003
Revises: 002
Create Date: 2026-11-18

Idempotent on purpose: run_migrations.py stamps an untracked database at 002 when
bookmarks.session_id exists, and a database created by a newer create_all already
has these columns. Existing custom profiles keep NULL owners and become invisible.
"""

from typing import Sequence, Union

from alembic import op

revision: str = "003"
down_revision: Union[str, None] = "002"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.execute("ALTER TABLE scoring_profiles ADD COLUMN IF NOT EXISTS session_id VARCHAR(36)")
    op.execute("ALTER TABLE scoring_profiles ADD COLUMN IF NOT EXISTS api_key_hash VARCHAR(64)")
    op.execute(
        "CREATE INDEX IF NOT EXISTS ix_scoring_profiles_session_id "
        "ON scoring_profiles (session_id)"
    )
    op.execute(
        "CREATE INDEX IF NOT EXISTS ix_scoring_profiles_api_key_hash "
        "ON scoring_profiles (api_key_hash)"
    )


def downgrade() -> None:
    op.execute("DROP INDEX IF EXISTS ix_scoring_profiles_api_key_hash")
    op.execute("DROP INDEX IF EXISTS ix_scoring_profiles_session_id")
    op.execute("ALTER TABLE scoring_profiles DROP COLUMN IF EXISTS api_key_hash")
    op.execute("ALTER TABLE scoring_profiles DROP COLUMN IF EXISTS session_id")
