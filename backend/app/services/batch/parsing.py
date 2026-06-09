"""Shared parsing helpers for the batch API routes."""

import json
from typing import Any


def parse_json_field(value: Any) -> Any:
    """Parse a field that may be stored as a JSON string or already a dict/list."""
    return json.loads(value) if isinstance(value, str) else value
