"""
Database ORM Foundation

SQLAlchemy async engine, session factory, and declarative base for ChemAudit.
"""

from sqlalchemy.ext.asyncio import AsyncAttrs, async_sessionmaker, create_async_engine
from sqlalchemy.orm import DeclarativeBase

from app.core.config import settings

engine = create_async_engine(settings.DATABASE_URL, echo=False, pool_pre_ping=True)
async_session = async_sessionmaker(engine, expire_on_commit=False)


class Base(AsyncAttrs, DeclarativeBase):
    """Base class for all ORM models."""

    pass


async def get_db():
    """FastAPI dependency that yields an async database session."""
    async with async_session() as session:
        yield session
