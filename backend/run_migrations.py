"""
Startup migration runner.

Ensures Alembic migrations are applied before the app starts.
Handles three scenarios:
  1. Fresh DB (no tables) — skip, create_all will handle it, then stamp at head
  2. Existing DB without alembic tracking — stamp current state, then upgrade
  3. Existing DB with alembic tracking — just upgrade
"""

import asyncio
import logging
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("migrations")

CREATE_ALEMBIC_TABLE_SQL = (
    "CREATE TABLE IF NOT EXISTS alembic_version "
    "(version_num VARCHAR(32) NOT NULL, "
    "CONSTRAINT alembic_version_pkc PRIMARY KEY (version_num))"
)


def _inspect_db(sync_conn) -> dict:
    """Inspect database state synchronously (called via run_sync)."""
    from sqlalchemy import inspect as sa_inspect

    inspector = sa_inspect(sync_conn)
    tables = inspector.get_table_names()
    has_tables = "bookmarks" in tables
    has_session_id = False
    if has_tables:
        columns = [c["name"] for c in inspector.get_columns("bookmarks")]
        has_session_id = "session_id" in columns
    return {
        "has_alembic": "alembic_version" in tables,
        "has_tables": has_tables,
        "has_session_id": has_session_id,
    }


async def _stamp_alembic(engine, revision: str) -> None:
    """Create the alembic_version table and stamp it at the given revision."""
    from sqlalchemy import text

    async with engine.begin() as conn:
        await conn.execute(text(CREATE_ALEMBIC_TABLE_SQL))
        await conn.execute(text("DELETE FROM alembic_version"))
        await conn.execute(
            text("INSERT INTO alembic_version (version_num) VALUES (:rev)"),
            {"rev": revision},
        )
    logger.info("Stamped alembic_version at revision %s", revision)


async def _inspect_and_stamp() -> bool:
    """Inspect DB and stamp alembic if needed. Returns whether upgrade is needed."""
    from sqlalchemy.ext.asyncio import create_async_engine

    db_url = os.environ.get("DATABASE_URL", "")
    if not db_url:
        logger.warning("DATABASE_URL not set, skipping migrations")
        return False

    engine = create_async_engine(db_url)
    try:
        async with engine.connect() as conn:
            info = await conn.run_sync(_inspect_db)

        if not info["has_tables"]:
            logger.info("No tables found — app startup will create them via create_all")
            return False

        if not info["has_alembic"]:
            logger.info("No alembic_version table — stamping current schema state")
            revision = "002" if info["has_session_id"] else "001"
            await _stamp_alembic(engine, revision)
    finally:
        await engine.dispose()

    return True


def main() -> None:
    """Run database migrations: inspect, stamp if needed, then upgrade."""
    needs_upgrade = asyncio.run(_inspect_and_stamp())
    if not needs_upgrade:
        return

    from alembic.config import Config

    from alembic import command

    logger.info("Running alembic upgrade head...")
    cfg = Config("alembic.ini")
    command.upgrade(cfg, "head")
    logger.info("Migrations complete")


if __name__ == "__main__":
    main()
