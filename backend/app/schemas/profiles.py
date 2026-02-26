"""
Scoring Profile Schemas

Pydantic schemas for custom scoring profiles with property thresholds and weights.
"""

from datetime import datetime
from typing import Any, Dict, Optional

from pydantic import BaseModel, Field


class ThresholdRange(BaseModel):
    """Min/max range for a single property threshold."""

    min: Optional[float] = None
    max: Optional[float] = None


class ScoringProfileBase(BaseModel):
    """Base fields for a scoring profile."""

    name: str = Field(..., min_length=1, max_length=100, description="Profile name")
    description: Optional[str] = Field(None, description="Profile description")
    thresholds: Dict[str, Dict[str, Any]] = Field(
        default_factory=dict,
        description="Property thresholds, e.g. {'mw': {'min': 150, 'max': 500}}",
    )
    weights: Dict[str, float] = Field(
        default_factory=dict,
        description="Property weights, e.g. {'mw': 1.0, 'logp': 0.8}",
    )


class ScoringProfileCreate(ScoringProfileBase):
    """Schema for creating a new scoring profile."""

    pass


class ScoringProfileUpdate(BaseModel):
    """Schema for updating a scoring profile (all fields optional)."""

    name: Optional[str] = Field(None, min_length=1, max_length=100)
    description: Optional[str] = None
    thresholds: Optional[Dict[str, Dict[str, Any]]] = None
    weights: Optional[Dict[str, float]] = None


class ScoringProfileResponse(ScoringProfileBase):
    """Schema for scoring profile response."""

    id: int
    is_preset: bool = False
    is_active: bool = True
    created_at: datetime
    updated_at: datetime

    model_config = {"from_attributes": True}


class ScoringProfileExport(BaseModel):
    """Schema for JSON export/import of a scoring profile."""

    name: str
    description: Optional[str] = None
    thresholds: Dict[str, Dict[str, Any]] = Field(default_factory=dict)
    weights: Dict[str, float] = Field(default_factory=dict)


class DuplicateRequest(BaseModel):
    """Request body for duplicating a profile."""

    name: str = Field(..., min_length=1, max_length=100, description="Name for the copy")
