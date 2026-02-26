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


def _inspect_db(sync_conn) -> dict:
    """Inspect database state synchronously (called via run_sync)."""
    from sqlalchemy import inspect

    insp = inspect(sync_conn)
    tables = insp.get_table_names()
    has_tables = "bookmarks" in tables
    has_alembic = "alembic_version" in tables
    has_session_id = False
    if has_tables:
        columns = [c["name"] for c in insp.get_columns("bookmarks")]
        has_session_id = "session_id" in columns
    return {
        "has_alembic": has_alembic,
        "has_tables": has_tables,
        "has_session_id": has_session_id,
    }


async def _inspect_and_stamp():
    """Inspect DB and stamp alembic if needed. Returns whether upgrade is needed."""
    from sqlalchemy import text
    from sqlalchemy.ext.asyncio import create_async_engine

    db_url = os.environ.get("DATABASE_URL", "")
    if not db_url:
        logger.warning("DATABASE_URL not set, skipping migrations")
        return False

    engine = create_async_engine(db_url)

    try:
        async with engine.connect() as conn:
            info = await conn.run_sync(_inspect_db)
    finally:
        await engine.dispose()

    has_alembic = info["has_alembic"]
    has_tables = info["has_tables"]
    has_session_id = info["has_session_id"]

    if not has_tables:
        logger.info("No tables found — app startup will create them via create_all")
        return False

    if not has_alembic:
        logger.info("No alembic_version table — stamping current schema state")
        stamp_rev = "002" if has_session_id else "001"

        engine = create_async_engine(db_url)
        try:
            async with engine.begin() as conn:
                await conn.execute(text(
                    "CREATE TABLE IF NOT EXISTS alembic_version "
                    "(version_num VARCHAR(32) NOT NULL, "
                    "CONSTRAINT alembic_version_pkc PRIMARY KEY (version_num))"
                ))
                await conn.execute(text("DELETE FROM alembic_version"))
                await conn.execute(text(
                    f"INSERT INTO alembic_version (version_num) VALUES ('{stamp_rev}')"
                ))
            logger.info("Stamped alembic_version at revision %s", stamp_rev)
        finally:
            await engine.dispose()

    return True  # Upgrade needed


def main():
    # Phase 1: Async inspect + stamp
    needs_upgrade = asyncio.run(_inspect_and_stamp())

    if not needs_upgrade:
        return

    # Phase 2: Sync alembic upgrade (has its own asyncio.run inside env.py)
    logger.info("Running alembic upgrade head...")
    from alembic.config import Config
    from alembic import command

    cfg = Config("alembic.ini")
    command.upgrade(cfg, "head")
    logger.info("Migrations complete")


if __name__ == "__main__":
    main()
