"""
Custom exceptions and exception handlers for ChemStructVal.

All application exceptions inherit from ChemStructValException and
return structured JSON responses.
"""
from fastapi import Request, status
from fastapi.responses import JSONResponse
from typing import Optional


class ChemStructValException(Exception):
    """Base exception for ChemStructVal application."""

    def __init__(
        self,
        message: str,
        status_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR,
        details: Optional[dict] = None,
    ):
        self.message = message
        self.status_code = status_code
        self.details = details or {}
        super().__init__(self.message)


class ParseError(ChemStructValException):
    """Exception raised when molecule parsing fails."""

    def __init__(self, message: str = "Failed to parse molecule", details: Optional[dict] = None):
        super().__init__(
            message=message,
            status_code=status.HTTP_400_BAD_REQUEST,
            details=details,
        )


class ValidationError(ChemStructValException):
    """Exception raised when molecule validation fails."""

    def __init__(
        self, message: str = "Validation failed", details: Optional[dict] = None
    ):
        super().__init__(
            message=message,
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            details=details,
        )


class NotFoundError(ChemStructValException):
    """Exception raised when requested resource is not found."""

    def __init__(self, message: str = "Resource not found", details: Optional[dict] = None):
        super().__init__(
            message=message,
            status_code=status.HTTP_404_NOT_FOUND,
            details=details,
        )


async def chemstructval_exception_handler(
    request: Request, exc: ChemStructValException
) -> JSONResponse:
    """Handle ChemStructVal exceptions by returning structured JSON."""
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "error": exc.message,
            "details": exc.details,
        },
    )


async def generic_exception_handler(request: Request, exc: Exception) -> JSONResponse:
    """Handle unexpected exceptions by returning generic error."""
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={
            "error": "Internal server error",
            "details": {"message": str(exc)},
        },
    )
