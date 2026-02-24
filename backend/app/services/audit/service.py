"""
Validation Audit Trail Service

Append-only logging of all validation events. This service only performs
INSERT operations — no UPDATE or DELETE on audit entries.
"""

import logging
from typing import Optional

from sqlalchemy.ext.asyncio import AsyncSession

from app.db.models.audit import ValidationAuditEntry

logger = logging.getLogger(__name__)


async def log_validation_event(
    db: AsyncSession,
    smiles: str,
    inchikey: Optional[str],
    outcome: str,
    score: Optional[int],
    source: str,
    job_id: Optional[str] = None,
    molecule_count: Optional[int] = None,
    pass_count: Optional[int] = None,
    fail_count: Optional[int] = None,
    api_key_hash: Optional[str] = None,
) -> ValidationAuditEntry:
    """
    Log a validation event to the audit trail (append-only).

    Args:
        db: AsyncSession
        smiles: SMILES string (empty for batch entries)
        inchikey: InChIKey if available
        outcome: 'pass', 'warn', or 'fail'
        score: Validation score (0-100)
        source: 'single' or 'batch'
        job_id: Batch job ID (for batch events)
        molecule_count: Total molecules (for batch events)
        pass_count: Passed count (for batch events)
        fail_count: Failed count (for batch events)
        api_key_hash: Hashed API key if present

    Returns:
        Created audit entry
    """
    entry = ValidationAuditEntry(
        smiles=smiles,
        inchikey=inchikey,
        outcome=outcome,
        score=score,
        source=source,
        job_id=job_id,
        molecule_count=molecule_count,
        pass_count=pass_count,
        fail_count=fail_count,
        api_key_hash=api_key_hash,
    )
    db.add(entry)
    await db.commit()
    await db.refresh(entry)
    return entry


async def log_batch_event(
    db: AsyncSession,
    job_id: str,
    molecule_count: int,
    pass_count: int,
    fail_count: int,
    api_key_hash: Optional[str] = None,
) -> ValidationAuditEntry:
    """
    Log a batch completion event to the audit trail.

    Convenience wrapper that computes outcome from pass/fail counts.

    Args:
        db: AsyncSession
        job_id: Batch job ID
        molecule_count: Total molecules processed
        pass_count: Number passing validation
        fail_count: Number failing validation
        api_key_hash: Hashed API key if present

    Returns:
        Created audit entry
    """
    if fail_count == 0:
        outcome = "pass"
    elif fail_count < molecule_count / 2:
        outcome = "warn"
    else:
        outcome = "fail"

    return await log_validation_event(
        db=db,
        smiles="",
        inchikey=None,
        outcome=outcome,
        score=None,
        source="batch",
        job_id=job_id,
        molecule_count=molecule_count,
        pass_count=pass_count,
        fail_count=fail_count,
        api_key_hash=api_key_hash,
    )
