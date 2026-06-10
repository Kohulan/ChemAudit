"""
Database ORM Foundation

SQLAlchemy async engine, session factory, and declarative base for ChemAudit.
"""

from sqlalchemy import create_engine
from sqlalchemy.ext.asyncio import AsyncAttrs, async_sessionmaker, create_async_engine
from sqlalchemy.orm import DeclarativeBase, sessionmaker

from app.core.config import settings

# Async engine — used by the FastAPI request path (asyncpg / aiosqlite).
engine = create_async_engine(settings.DATABASE_URL, echo=False, pool_pre_ping=True)
async_session = async_sessionmaker(engine, expire_on_commit=False)


def _sync_db_url(url: str) -> str:
    """Map an async SQLAlchemy URL to its synchronous-driver equivalent.

    Celery tasks run in a synchronous context (the default prefork pool has no
    event loop), so they use a real sync engine instead of bridging async code
    with asyncio.run() — which is fragile under gevent/eventlet/asyncio pools.
    """
    return url.replace("+asyncpg", "+psycopg").replace("+aiosqlite", "+pysqlite")


# Sync engine + session factory — used by Celery tasks (audit writes, purges).
sync_engine = create_engine(_sync_db_url(settings.DATABASE_URL), echo=False, pool_pre_ping=True)
sync_session = sessionmaker(sync_engine, expire_on_commit=False)


class Base(AsyncAttrs, DeclarativeBase):
    """Base class for all ORM models."""

    pass


async def get_db():
    """FastAPI dependency that yields an async database session."""
    async with async_session() as session:
        yield session
