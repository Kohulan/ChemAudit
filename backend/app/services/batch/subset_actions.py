"""
Batch Subset Actions

Services for re-validating, re-scoring, and exporting subsets of batch results.
"""

import uuid
from io import BytesIO
from typing import List, Optional

from app.services.batch.result_aggregator import result_storage
from app.services.batch.tasks import process_batch_job
from app.services.export.base import ExporterFactory, ExportFormat


def _fetch_subset_smiles(job_id: str, indices: List[int]) -> List[dict]:
    """
    Fetch molecule data for a subset of indices from an existing batch job.

    Args:
        job_id: Original job ID
        indices: Molecule indices to extract

    Returns:
        List of molecule dicts suitable for batch processing

    Raises:
        ValueError: If no results found for the given job/indices
    """
    result_data = result_storage.get_results(
        job_id=job_id,
        page=1,
        page_size=10000,
    )
    results = result_data.get("results", [])

    indices_set = set(indices)
    subset = [r for r in results if r.get("index") in indices_set]

    if not subset:
        raise ValueError(f"No results found for job {job_id} with given indices")

    mol_dicts = [
        {
            "smiles": r.get("smiles", ""),
            "name": r.get("name") or f"mol_{r.get('index', i)}",
            "index": i,
            "properties": {},
            "parse_error": None,
        }
        for i, r in enumerate(subset)
    ]
    return mol_dicts


def revalidate_subset(job_id: str, indices: List[int]) -> str:
    """
    Re-validate a subset of molecules from an existing batch job.

    Creates a new batch job with the selected molecules.

    Args:
        job_id: Original job ID
        indices: Molecule indices to re-validate

    Returns:
        New job_id
    """
    mol_dicts = _fetch_subset_smiles(job_id, indices)
    new_job_id = str(uuid.uuid4())
    process_batch_job(new_job_id, mol_dicts, safety_options={})
    return new_job_id


def rescore_subset(
    job_id: str, indices: List[int], profile_id: Optional[int] = None
) -> str:
    """
    Re-score a subset of molecules, optionally with a custom scoring profile.

    Creates a new batch job for the selected molecules.

    Args:
        job_id: Original job ID
        indices: Molecule indices to re-score
        profile_id: Optional scoring profile ID to apply

    Returns:
        New job_id
    """
    mol_dicts = _fetch_subset_smiles(job_id, indices)
    new_job_id = str(uuid.uuid4())

    # Pass profile_id in safety_options for downstream handling
    safety_options = {}
    if profile_id is not None:
        safety_options["scoring_profile_id"] = profile_id

    process_batch_job(new_job_id, mol_dicts, safety_options=safety_options)
    return new_job_id


def export_subset(job_id: str, indices: List[int], format: str) -> BytesIO:
    """
    Export a subset of molecules in the specified format.

    Args:
        job_id: Original job ID
        indices: Molecule indices to export
        format: Export format string (e.g., "csv", "excel", "sdf")

    Returns:
        BytesIO buffer with exported data
    """
    result_data = result_storage.get_results(
        job_id=job_id,
        page=1,
        page_size=10000,
    )
    results = result_data.get("results", [])

    indices_set = set(indices)
    subset = [r for r in results if r.get("index") in indices_set]

    if not subset:
        raise ValueError(f"No results found for job {job_id} with given indices")

    export_format = ExportFormat(format)
    exporter = ExporterFactory.create(export_format)
    return exporter.export(subset)
