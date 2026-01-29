"""
Custom exceptions and exception handlers for ChemVault.

All application exceptions inherit from ChemVaultException and
return structured JSON responses.
"""

from fastapi import Request, status
from fastapi.responses import JSONResponse

from app.core.config import settings


class ChemVaultException(Exception):
    """Base exception for ChemVault application."""

    def __init__(
        self,
        message: str,
        status_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR,
        details: dict | None = None,
    ):
        self.message = message
        self.status_code = status_code
        self.details = details or {}
        super().__init__(self.message)


class ParseError(ChemVaultException):
    """Exception raised when molecule parsing fails."""

    def __init__(
        self, message: str = "Failed to parse molecule", details: dict | None = None
    ):
        super().__init__(
            message=message,
            status_code=status.HTTP_400_BAD_REQUEST,
            details=details,
        )


class ValidationError(ChemVaultException):
    """Exception raised when molecule validation fails."""

    def __init__(self, message: str = "Validation failed", details: dict | None = None):
        super().__init__(
            message=message,
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            details=details,
        )


class NotFoundError(ChemVaultException):
    """Exception raised when requested resource is not found."""

    def __init__(
        self, message: str = "Resource not found", details: dict | None = None
    ):
        super().__init__(
            message=message,
            status_code=status.HTTP_404_NOT_FOUND,
            details=details,
        )


async def chemvault_exception_handler(
    request: Request, exc: ChemVaultException
) -> JSONResponse:
    """Handle ChemVault exceptions by returning structured JSON."""
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "error": exc.message,
            "details": exc.details,
        },
    )


async def generic_exception_handler(request: Request, exc: Exception) -> JSONResponse:
    """Handle unexpected exceptions by returning generic error."""
    content = {"error": "Internal server error"}
    # Only expose exception details in debug mode to prevent information leakage
    if settings.DEBUG:
        content["details"] = {"message": str(exc)}
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content=content,
    )
