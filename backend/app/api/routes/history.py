"""
Audit Trail History API Routes

Paginated read access to the append-only validation audit trail.
No UPDATE or DELETE endpoints — audit trail is immutable.
"""

from datetime import datetime
from typing import Optional

from fastapi import APIRouter, Depends, Query, Request
from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.security import get_api_key, hash_api_key_for_lookup
from app.core.session import get_data_scope, set_rls_context
from app.db import get_db
from app.db.models.audit import ValidationAuditEntry
from app.schemas.history import AuditEntryResponse, AuditHistoryResponse

router = APIRouter()


@router.get("/history", response_model=AuditHistoryResponse)
async def get_history(
    request: Request,
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=50, ge=1, le=200),
    date_from: Optional[datetime] = Query(None, description="Filter from date (ISO8601)"),
    date_to: Optional[datetime] = Query(None, description="Filter to date (ISO8601)"),
    outcome: Optional[str] = Query(None, description="Filter by outcome (pass, warn, fail)"),
    source: Optional[str] = Query(None, description="Filter by source (single, batch)"),
    smiles_search: Optional[str] = Query(None, description="SMILES substring search"),
    api_key: Optional[str] = Depends(get_api_key),
    db: AsyncSession = Depends(get_db),
):
    """
    Get paginated audit trail with optional filters.

    Results are scoped to the caller's session or API key.
    Results are ordered by created_at DESC (newest first).
    """
    # Determine data scope from session cookie or API key
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
    await set_rls_context(db, session_id, api_key_hash)

    # Build base query
    query = select(ValidationAuditEntry)
    count_query = select(func.count(ValidationAuditEntry.id))

    # Apply session/API key scope filter
    if api_key_hash:
        query = query.where(ValidationAuditEntry.api_key_hash == api_key_hash)
        count_query = count_query.where(ValidationAuditEntry.api_key_hash == api_key_hash)
    elif session_id:
        query = query.where(ValidationAuditEntry.session_id == session_id)
        count_query = count_query.where(ValidationAuditEntry.session_id == session_id)
    else:
        return AuditHistoryResponse(entries=[], total=0, page=page, page_size=page_size)

    # Apply user-specified filters
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
        count_query = count_query.where(ValidationAuditEntry.smiles.contains(smiles_search))

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
async def get_history_stats(
    request: Request,
    api_key: Optional[str] = Depends(get_api_key),
    db: AsyncSession = Depends(get_db),
):
    """
    Get summary statistics for the audit trail.

    Results are scoped to the caller's session or API key.
    Returns total validations, outcome distribution, and source distribution.
    """
    # Determine data scope from session cookie or API key
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
    await set_rls_context(db, session_id, api_key_hash)

    # Build scope filter
    if api_key_hash:
        scope_filter = ValidationAuditEntry.api_key_hash == api_key_hash
    elif session_id:
        scope_filter = ValidationAuditEntry.session_id == session_id
    else:
        return {
            "total_validations": 0,
            "outcome_distribution": {},
            "source_distribution": {},
        }

    # Total count
    total_result = await db.execute(select(func.count(ValidationAuditEntry.id)).where(scope_filter))
    total = total_result.scalar() or 0

    # Outcome distribution
    outcome_result = await db.execute(
        select(
            ValidationAuditEntry.outcome,
            func.count(ValidationAuditEntry.id),
        )
        .where(scope_filter)
        .group_by(ValidationAuditEntry.outcome)
    )
    outcome_dist = {row[0]: row[1] for row in outcome_result.all()}

    # Source distribution
    source_result = await db.execute(
        select(
            ValidationAuditEntry.source,
            func.count(ValidationAuditEntry.id),
        )
        .where(scope_filter)
        .group_by(ValidationAuditEntry.source)
    )
    source_dist = {row[0]: row[1] for row in source_result.all()}

    return {
        "total_validations": total,
        "outcome_distribution": outcome_dist,
        "source_distribution": source_dist,
    }
