"""
Universal Identifier Resolution API Route.

Accepts any chemical identifier and resolves it to a canonical structure.
"""

from typing import Optional

from fastapi import APIRouter, Depends, Request
from pydantic import BaseModel, Field, field_validator

from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.integrations import ResolvedCompound
from app.services.integrations.resolver import resolve_identifier

router = APIRouter()


class ResolveRequest(BaseModel):
    """Request body for identifier resolution."""

    identifier: str = Field(..., min_length=1, description="Any chemical identifier")
    identifier_type: str = Field(
        default="auto",
        description=(
            "Explicit type hint (auto, smiles, inchi, inchikey, pubchem_cid, "
            "chembl_id, cas, drugbank_id, chebi_id, unii, wikipedia, name)"
        ),
    )

    @field_validator("identifier")
    @classmethod
    def strip_identifier(cls, v: str) -> str:
        return v.strip()


@router.post("/resolve", response_model=ResolvedCompound)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def resolve_compound(
    request: Request,
    body: ResolveRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Resolve any chemical identifier to a canonical structure.

    Accepts: compound name, SMILES, InChI, InChIKey, PubChem CID,
    ChEMBL ID, CAS number, DrugBank ID, ChEBI ID, UNII, or Wikipedia URL.

    Returns canonical SMILES, InChIKey, molecular properties,
    resolution provenance chain, and cross-database references.
    """
    return await resolve_identifier(body.identifier, body.identifier_type)
