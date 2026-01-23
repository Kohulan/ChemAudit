"""
Export API Routes

Endpoints for exporting batch results in various formats.
"""
from datetime import datetime
from typing import Optional

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import StreamingResponse

from app.services.export.base import ExportFormat, ExporterFactory
from app.services.batch.result_aggregator import result_storage
from app.services.batch.progress_tracker import progress_tracker


router = APIRouter()

# Chunk size for streaming (1MB)
CHUNK_SIZE = 1024 * 1024


@router.get("/batch/{job_id}/export")
async def export_batch_results(
    job_id: str,
    format: ExportFormat = Query(..., description="Export format"),
    score_min: Optional[int] = Query(
        None, ge=0, le=100, description="Minimum validation score"
    ),
    score_max: Optional[int] = Query(
        None, ge=0, le=100, description="Maximum validation score"
    ),
    status: Optional[str] = Query(
        None, description="Filter by status (success, error, warning)"
    ),
):
    """
    Export batch results to specified format.

    - **job_id**: Job identifier from upload response
    - **format**: Export format (csv, excel, sdf, json)
    - **score_min**: Filter results by minimum score
    - **score_max**: Filter results by maximum score
    - **status**: Filter results by status

    Returns file download with appropriate Content-Disposition header.

    Supported formats:
    - CSV: Comma-separated values with all validation data
    - Excel: Formatted spreadsheet with conditional coloring and summary sheet
    - SDF: Structure-data file with properties attached to molecules
    - JSON: Full result objects with metadata
    """
    # Check if job exists
    progress_info = progress_tracker.get_progress(job_id)
    if not progress_info:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    # Get results with filters
    result_data = result_storage.get_results(
        job_id=job_id,
        page=1,
        page_size=10000,  # Get all results (max batch size is 10K)
        status_filter=status,
        min_score=score_min,
        max_score=score_max,
    )

    results = result_data.get("results", [])

    # Handle empty results
    if not results:
        raise HTTPException(
            status_code=404,
            detail="No results found matching the specified filters",
        )

    # Create exporter
    try:
        exporter = ExporterFactory.create(format)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))

    # Export results
    try:
        export_buffer = exporter.export(results)
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to export results: {str(e)}",
        )

    # Generate filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"batch_{job_id[:8]}_{timestamp}.{exporter.file_extension}"

    # Return streaming response
    def iterfile():
        """Stream file in chunks."""
        export_buffer.seek(0)
        while chunk := export_buffer.read(CHUNK_SIZE):
            yield chunk

    return StreamingResponse(
        iterfile(),
        media_type=exporter.media_type,
        headers={
            "Content-Disposition": f'attachment; filename="{filename}"',
        },
    )
