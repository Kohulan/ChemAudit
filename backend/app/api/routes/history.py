"""
Audit Trail History API Routes

Paginated read access to the append-only validation audit trail.
No UPDATE or DELETE endpoints — audit trail is immutable.
"""

from datetime import datetime
from typing import Optional

from fastapi import APIRouter, Depends, Query
from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession

from app.db import get_db
from app.db.models.audit import ValidationAuditEntry
from app.schemas.history import AuditEntryResponse, AuditHistoryResponse

router = APIRouter()


@router.get("/history", response_model=AuditHistoryResponse)
async def get_history(
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=50, ge=1, le=200),
    date_from: Optional[datetime] = Query(None, description="Filter from date (ISO8601)"),
    date_to: Optional[datetime] = Query(None, description="Filter to date (ISO8601)"),
    outcome: Optional[str] = Query(None, description="Filter by outcome (pass, warn, fail)"),
    source: Optional[str] = Query(None, description="Filter by source (single, batch)"),
    smiles_search: Optional[str] = Query(None, description="SMILES substring search"),
    api_key_hash: Optional[str] = Query(None, description="Filter by API key hash"),
    db: AsyncSession = Depends(get_db),
):
    """
    Get paginated audit trail with optional filters.

    Results are ordered by created_at DESC (newest first).
    """
    # Build base query
    query = select(ValidationAuditEntry)
    count_query = select(func.count(ValidationAuditEntry.id))

    # Apply filters
    if date_from is not None:
        query = query.where(ValidationAuditEntry.created_at >= date_from)
        count_query = count_query.where(ValidationAuditEntry.created_at >= date_from)
    if date_to is not None:
        query = query.where(ValidationAuditEntry.created_at <= date_to)
        count_query = count_query.where(ValidationAuditEntry.created_at <= date_to)
    if outcome is not None:
        query = query.where(ValidationAuditEntry.outcome == outcome)
        count_query = count_query.where(ValidationAuditEntry.outcome == outcome)
    if source is not None:
        query = query.where(ValidationAuditEntry.source == source)
        count_query = count_query.where(ValidationAuditEntry.source == source)
    if smiles_search is not None:
        query = query.where(ValidationAuditEntry.smiles.contains(smiles_search))
        count_query = count_query.where(
            ValidationAuditEntry.smiles.contains(smiles_search)
        )
    if api_key_hash is not None:
        query = query.where(ValidationAuditEntry.api_key_hash == api_key_hash)
        count_query = count_query.where(
            ValidationAuditEntry.api_key_hash == api_key_hash
        )

    # Get total count
    total_result = await db.execute(count_query)
    total = total_result.scalar() or 0

    # Paginate
    offset = (page - 1) * page_size
    query = query.order_by(ValidationAuditEntry.created_at.desc())
    query = query.offset(offset).limit(page_size)

    result = await db.execute(query)
    entries = result.scalars().all()

    return AuditHistoryResponse(
        entries=[
            AuditEntryResponse(
                id=e.id,
                smiles=e.smiles,
                inchikey=e.inchikey,
                outcome=e.outcome,
                score=e.score,
                job_id=e.job_id,
                molecule_count=e.molecule_count,
                pass_count=e.pass_count,
                fail_count=e.fail_count,
                source=e.source,
                created_at=e.created_at,
            )
            for e in entries
        ],
        total=total,
        page=page,
        page_size=page_size,
    )


@router.get("/history/stats")
async def get_history_stats(db: AsyncSession = Depends(get_db)):
    """
    Get summary statistics for the audit trail.

    Returns total validations and outcome distribution.
    """
    # Total count
    total_result = await db.execute(
        select(func.count(ValidationAuditEntry.id))
    )
    total = total_result.scalar() or 0

    # Outcome distribution
    outcome_result = await db.execute(
        select(
            ValidationAuditEntry.outcome,
            func.count(ValidationAuditEntry.id),
        ).group_by(ValidationAuditEntry.outcome)
    )
    outcome_dist = {row[0]: row[1] for row in outcome_result.all()}

    # Source distribution
    source_result = await db.execute(
        select(
            ValidationAuditEntry.source,
            func.count(ValidationAuditEntry.id),
        ).group_by(ValidationAuditEntry.source)
    )
    source_dist = {row[0]: row[1] for row in source_result.all()}

    return {
        "total_validations": total,
        "outcome_distribution": outcome_dist,
        "source_distribution": source_dist,
    }
