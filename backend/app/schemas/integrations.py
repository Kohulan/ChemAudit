"""
Schemas for external integration requests and responses.

Covers COCONUT, PubChem, ChEMBL integrations, Universal Identifier Resolution,
and Cross-Database Comparison.
"""

from typing import Optional

from pydantic import BaseModel, Field


# COCONUT Natural Products Database
class COCONUTRequest(BaseModel):
    """Request for COCONUT lookup."""

    smiles: Optional[str] = None
    inchikey: Optional[str] = None


class COCONUTResult(BaseModel):
    """COCONUT natural product result."""

    found: bool
    coconut_id: Optional[str] = None
    name: Optional[str] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    organism: Optional[str] = None
    organism_type: Optional[str] = None
    nplikeness: Optional[float] = None
    url: Optional[str] = None


# PubChem Cross-Reference
class PubChemRequest(BaseModel):
    """Request for PubChem lookup."""

    smiles: Optional[str] = None
    inchikey: Optional[str] = None


class PubChemResult(BaseModel):
    """PubChem compound result."""

    found: bool
    cid: Optional[int] = None
    iupac_name: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    synonyms: Optional[list[str]] = None
    url: Optional[str] = None


# ChEMBL Bioactivity Data
class ChEMBLRequest(BaseModel):
    """Request for ChEMBL bioactivity lookup."""

    smiles: Optional[str] = None
    inchikey: Optional[str] = None


class BioactivityData(BaseModel):
    """ChEMBL bioactivity record."""

    target_chembl_id: str
    target_name: Optional[str] = None
    target_type: Optional[str] = None
    activity_type: str
    activity_value: Optional[float] = None
    activity_unit: Optional[str] = None
    assay_chembl_id: str
    document_chembl_id: Optional[str] = None


class ChEMBLResult(BaseModel):
    """ChEMBL molecule and bioactivity result."""

    found: bool
    chembl_id: Optional[str] = None
    pref_name: Optional[str] = None
    molecule_type: Optional[str] = None
    max_phase: Optional[int] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    bioactivities: list[BioactivityData] = Field(default_factory=list)
    bioactivity_count: int = 0
    url: Optional[str] = None


# =============================================================================
# Universal Identifier Resolution
# =============================================================================


class CrossReferences(BaseModel):
    """Cross-database references for a resolved compound."""

    pubchem_cid: Optional[int] = None
    chembl_id: Optional[str] = None
    coconut_id: Optional[str] = None
    drugbank_id: Optional[str] = None
    chebi_id: Optional[str] = None
    unii: Optional[str] = None
    cas: Optional[str] = None
    wikipedia_url: Optional[str] = None
    kegg_id: Optional[str] = None


class ResolvedCompound(BaseModel):
    """Result of universal identifier resolution."""

    resolved: bool
    identifier_type_detected: str = "unknown"
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    iupac_name: Optional[str] = None
    resolution_source: str = "none"
    resolution_chain: list[str] = Field(default_factory=list)
    cross_references: CrossReferences = Field(default_factory=CrossReferences)
    confidence: str = "none"


# =============================================================================
# Cross-Database Comparison
# =============================================================================


class DatabaseEntry(BaseModel):
    """A single database's representation of a compound."""

    database: str
    found: bool
    canonical_smiles: Optional[str] = None
    kekulized_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    name: Optional[str] = None
    url: Optional[str] = None


class PropertyComparison(BaseModel):
    """Comparison of a single property across databases."""

    property_name: str
    values: dict[str, Optional[str]] = Field(default_factory=dict)
    status: str = "missing"
    detail: Optional[str] = None


class ConsistencyResult(BaseModel):
    """Result of cross-database comparison."""

    entries: list[DatabaseEntry] = Field(default_factory=list)
    comparisons: list[PropertyComparison] = Field(default_factory=list)
    overall_verdict: str = "no_data"
    summary: str = ""


# =============================================================================
# SureChEMBL Patent Lookup
# =============================================================================


class SureChEMBLRequest(BaseModel):
    """Request for SureChEMBL patent lookup."""

    smiles: Optional[str] = None
    inchikey: Optional[str] = None


class SureChEMBLResult(BaseModel):
    """SureChEMBL patent presence result."""

    found: bool
    schembl_id: Optional[str] = None
    url: Optional[str] = None
    patent_count: Optional[int] = None
    source: Optional[str] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_weight: Optional[float] = None
