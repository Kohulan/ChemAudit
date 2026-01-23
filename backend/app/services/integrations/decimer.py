"""
DECIMER OCSR integration client.

DECIMER (Deep Learning for Chemical Image Recognition) is an OCSR tool
that converts chemical structure images to SMILES strings.

This client validates DECIMER output and provides confidence scoring.
"""
import httpx
from typing import Optional
from rdkit import Chem

from app.core.config import settings
from app.schemas.integrations import DECIMERRequest, DECIMERValidation


class DECIMERClient:
    """
    Client for DECIMER OCSR API integration.

    Note: DECIMER API is used for image-to-SMILES conversion.
    This client focuses on validating the OCSR output.
    """

    def __init__(self):
        self.base_url = settings.DECIMER_API_URL
        self.timeout = settings.EXTERNAL_API_TIMEOUT

    async def ocsr_from_image(self, image_data: bytes) -> Optional[str]:
        """
        Convert chemical structure image to SMILES using DECIMER OCSR.

        Args:
            image_data: Image file bytes (PNG, JPG, etc.)

        Returns:
            SMILES string if successful, None if failed
        """
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.post(
                    f"{self.base_url}/ocsr",
                    files={"image": image_data},
                )
                response.raise_for_status()

                data = response.json()
                return data.get("smiles")

        except (httpx.HTTPError, KeyError, ValueError):
            # External API failure - return None gracefully
            return None


def validate_ocsr_result(request: DECIMERRequest) -> DECIMERValidation:
    """
    Validate DECIMER OCSR output SMILES string.

    Checks if the SMILES is valid and provides confidence-adjusted validation.

    Args:
        request: DECIMER validation request with SMILES and optional confidence

    Returns:
        Validation result with canonical SMILES and identifiers
    """
    smiles = request.smiles
    confidence = request.confidence

    # Try to parse the SMILES
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)

        if mol is None:
            return DECIMERValidation(
                smiles=smiles,
                confidence=confidence,
                is_valid=False,
                validation_message="Invalid SMILES string - cannot parse",
            )

        # Try to sanitize
        try:
            Chem.SanitizeMol(mol)
        except (ValueError, RuntimeError) as e:
            return DECIMERValidation(
                smiles=smiles,
                confidence=confidence,
                is_valid=False,
                validation_message=f"Invalid chemistry: {str(e)}",
            )

        # Generate canonical SMILES and identifiers
        canonical_smiles = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol)

        # Confidence-based validation message
        if confidence is not None:
            if confidence >= 0.9:
                message = "High confidence OCSR result - valid structure"
            elif confidence >= 0.7:
                message = "Moderate confidence OCSR result - manual verification recommended"
            else:
                message = "Low confidence OCSR result - manual verification strongly recommended"
        else:
            message = "Valid structure - no confidence score available"

        return DECIMERValidation(
            smiles=smiles,
            confidence=confidence,
            is_valid=True,
            validation_message=message,
            canonical_smiles=canonical_smiles,
            inchi=inchi,
            inchikey=inchikey,
        )

    except Exception as e:
        return DECIMERValidation(
            smiles=smiles,
            confidence=confidence,
            is_valid=False,
            validation_message=f"Validation error: {str(e)}",
        )
