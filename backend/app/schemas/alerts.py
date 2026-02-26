"""
Alert Screening Schemas

Pydantic schemas for structural alert screening requests and responses.
"""

from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field, field_validator

from app.schemas.validation import MoleculeInfo


class AlertSeverity(str, Enum):
    """Severity level for structural alerts."""

    CRITICAL = "critical"
    WARNING = "warning"
    INFO = "info"


class AlertResultSchema(BaseModel):
    """Schema for a single structural alert match."""

    pattern_name: str = Field(..., description="Name of the matched pattern")
    description: str = Field(..., description="Description of the alert")
    severity: AlertSeverity = Field(
        default=AlertSeverity.WARNING, description="Alert severity level"
    )
    matched_atoms: List[int] = Field(
        default_factory=list, description="Atom indices matched by the pattern"
    )
    catalog_source: str = Field(..., description="Source catalog (PAINS, BRENK, etc.)")
    smarts: Optional[str] = Field(
        None, description="SMARTS pattern that matched (if available)"
    )
    reference: Optional[str] = Field(
        None, description="Citation for the original paper defining this alert"
    )
    scope: Optional[str] = Field(
        None, description="What aspect of compound quality this filter screens for"
    )
    filter_set: Optional[str] = Field(
        None, description="Filter set identifier from RDKit (e.g. ChEMBL23_Dundee)"
    )
    catalog_description: Optional[str] = Field(
        None, description="Human-readable name of the source catalog"
    )
    category: Optional[str] = Field(
        None,
        description="Concern category: Reactive Group, Metabolic Liability, Toxicophore, "
        "Assay Interference, Physicochemical, or Unwanted Functionality",
    )


class AlertScreenRequest(BaseModel):
    """Schema for alert screening request."""

    molecule: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="Molecule string (SMILES, InChI, or MOL block)",
    )
    format: str = Field(
        default="auto",
        pattern="^(auto|smiles|inchi|mol)$",
        description="Input format (auto-detected if not specified)",
    )
    catalogs: List[str] = Field(
        default=["PAINS"],
        description="List of catalogs to screen against (PAINS, BRENK, NIH, ZINC, ALL)",
    )

    @field_validator("molecule")
    @classmethod
    def sanitize_molecule_input(cls, v: str) -> str:
        """Sanitize molecule input to prevent injection attacks."""
        dangerous = ["<", ">", "&", ";", "|", "$", "`"]
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in molecule string")
        return v.strip()

    @field_validator("catalogs")
    @classmethod
    def validate_catalogs(cls, v: List[str]) -> List[str]:
        """Validate and normalize catalog names."""
        valid_catalogs = {
            "PAINS",
            "PAINS_A",
            "PAINS_B",
            "PAINS_C",
            "BRENK",
            "NIH",
            "ZINC",
            "CHEMBL_BMS",
            "CHEMBL_DUNDEE",
            "CHEMBL_GLAXO",
            "CHEMBL_INPHARMATICA",
            "CHEMBL_LINT",
            "CHEMBL_MLSMR",
            "CHEMBL_SURECHEMBL",
            "ALL",
        }
        normalized = []
        for cat in v:
            cat_upper = cat.upper()
            if cat_upper not in valid_catalogs:
                raise ValueError(
                    f"Unknown catalog: {cat}. Valid: {sorted(valid_catalogs)}"
                )
            normalized.append(cat_upper)
        return normalized


class AlertScreenResponse(BaseModel):
    """Schema for alert screening response."""

    status: str = Field(default="completed", description="Screening status")
    molecule_info: MoleculeInfo = Field(
        ..., description="Information about the screened molecule"
    )
    alerts: List[AlertResultSchema] = Field(
        default_factory=list, description="List of matched alerts"
    )
    total_alerts: int = Field(default=0, description="Total number of alerts found")
    screened_catalogs: List[str] = Field(
        default_factory=list, description="List of catalogs screened"
    )
    has_critical: bool = Field(
        default=False, description="True if any critical alerts found"
    )
    has_warning: bool = Field(
        default=False, description="True if any warning alerts found"
    )
    execution_time_ms: int = Field(..., description="Execution time in milliseconds")
    educational_note: str = Field(
        default="Structural alerts are warnings for investigation, not automatic failures. "
        "Many approved drugs contain PAINS patterns (e.g., 87 FDA-approved drugs).",
        description="Educational context about alert interpretation",
    )


class CatalogInfoSchema(BaseModel):
    """Schema for catalog information."""

    name: str = Field(..., description="Catalog name")
    description: str = Field(..., description="Catalog description")
    pattern_count: str = Field(..., description="Approximate pattern count")
    severity: str = Field(..., description="Default severity level for this catalog")
    note: Optional[str] = Field(None, description="Additional notes")
    reference: Optional[str] = Field(None, description="Citation for the original paper")
    scope: Optional[str] = Field(None, description="What the filter set screens for")
    doi: Optional[str] = Field(None, description="DOI link for the reference paper")
    pmid: Optional[str] = Field(None, description="PubMed ID for the reference paper")


class CatalogListResponse(BaseModel):
    """Schema for listing available catalogs."""

    catalogs: Dict[str, CatalogInfoSchema] = Field(
        ..., description="Available catalogs"
    )
    default_catalogs: List[str] = Field(
        default=["PAINS"], description="Default catalogs if none specified"
    )
