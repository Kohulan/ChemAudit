"""
Audit Trail History Schemas

Pydantic schemas for validation audit trail entries and history API.
"""

from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel


class AuditEntryResponse(BaseModel):
    """Schema for a single audit trail entry."""

    id: int
    smiles: str
    inchikey: Optional[str] = None
    outcome: str
    score: Optional[int] = None
    job_id: Optional[str] = None
    molecule_count: Optional[int] = None
    pass_count: Optional[int] = None
    fail_count: Optional[int] = None
    source: str
    created_at: datetime

    model_config = {"from_attributes": True}


class AuditHistoryResponse(BaseModel):
    """Schema for paginated audit trail response."""

    entries: List[AuditEntryResponse]
    total: int
    page: int
    page_size: int
