"""
Base Exporter and Factory Pattern

Provides abstract base class for exporters and factory for creating exporters.
"""
from abc import ABC, abstractmethod
from enum import Enum
from io import BytesIO
from typing import List, Dict, Any


class ExportFormat(str, Enum):
    """Export format options (extends str for JSON serialization)."""

    CSV = "csv"
    EXCEL = "excel"
    SDF = "sdf"
    JSON = "json"
    PDF = "pdf"


class BaseExporter(ABC):
    """Abstract base class for all exporters."""

    @abstractmethod
    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """
        Export batch results to specific format.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing exported data
        """
        pass

    @property
    @abstractmethod
    def media_type(self) -> str:
        """MIME type for this export format."""
        pass

    @property
    @abstractmethod
    def file_extension(self) -> str:
        """File extension for this export format (without dot)."""
        pass


class ExporterFactory:
    """Factory for creating exporters based on format."""

    _exporters: Dict[ExportFormat, type] = {}

    @classmethod
    def register(cls, format: ExportFormat, exporter_class: type) -> None:
        """Register an exporter class for a format."""
        cls._exporters[format] = exporter_class

    @classmethod
    def create(cls, format: ExportFormat) -> BaseExporter:
        """
        Create exporter instance for specified format.

        Args:
            format: ExportFormat enum value

        Returns:
            BaseExporter instance

        Raises:
            ValueError: If format not supported
        """
        exporter_class = cls._exporters.get(format)
        if not exporter_class:
            raise ValueError(f"Unsupported export format: {format}")
        return exporter_class()
