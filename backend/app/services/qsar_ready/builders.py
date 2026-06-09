"""QSAR-ready batch export and summary builders.

Response/summary construction helpers extracted from the QSAR-ready API
routes to keep the route handlers thin.
"""

import csv
import io
import json

from fastapi import (
    HTTPException,
)
from fastapi.responses import StreamingResponse

from app.schemas.qsar_ready import (
    QSARBatchSummary,
)


def _build_csv_response(job_id: str, results: list) -> StreamingResponse:
    """Build a CSV StreamingResponse from results (D-12)."""
    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=[
            "original_smiles",
            "curated_smiles",
            "original_inchikey",
            "curated_inchikey",
            "inchikey_changed",
            "status",
            "rejection_reason",
            "steps_applied",
        ],
    )
    writer.writeheader()
    for r in results:
        steps = r.get("steps", [])
        applied_steps = ",".join(
            s["step_name"] for s in steps if s.get("status") == "applied"
        )
        writer.writerow(
            {
                "original_smiles": r.get("original_smiles", ""),
                "curated_smiles": r.get("curated_smiles", ""),
                "original_inchikey": r.get("original_inchikey", ""),
                "curated_inchikey": r.get("standardized_inchikey", ""),
                "inchikey_changed": str(r.get("inchikey_changed", False)).lower(),
                "status": r.get("status", ""),
                "rejection_reason": r.get("rejection_reason", ""),
                "steps_applied": applied_steps,
            }
        )

    output.seek(0)
    return StreamingResponse(
        iter([output.getvalue()]),
        media_type="text/csv",
        headers={
            "Content-Disposition": f"attachment; filename=qsar_results_{job_id[:8]}.csv"
        },
    )


def _build_sdf_response(job_id: str, results: list) -> StreamingResponse:
    """Build an SDF StreamingResponse from curated_smiles (Pitfall 8 from RESEARCH.md)."""
    from rdkit import Chem

    sdf_blocks = []
    for r in results:
        curated_smiles = r.get("curated_smiles")
        if not curated_smiles:
            continue
        mol = Chem.MolFromSmiles(curated_smiles)
        if mol is None:
            continue
        # Use original InChIKey as molecule title
        title = r.get("original_inchikey") or r.get("original_smiles", "")[:60]
        mol.SetProp("_Name", title)
        mol.SetProp("original_smiles", r.get("original_smiles", ""))
        mol.SetProp("curated_smiles", curated_smiles)
        mol.SetProp("status", r.get("status", ""))
        try:
            block = Chem.MolToMolBlock(mol)
            sdf_blocks.append(block)
            sdf_blocks.append("$$$$\n")
        except Exception:
            continue

    sdf_content = "\n".join(sdf_blocks) if sdf_blocks else ""

    return StreamingResponse(
        iter([sdf_content]),
        media_type="chemical/x-mdl-sdfile",
        headers={
            "Content-Disposition": f"attachment; filename=qsar_results_{job_id[:8]}.sdf"
        },
    )


def _build_json_response(job_id: str, results: list, config_dict: dict) -> StreamingResponse:
    """Build a full-provenance JSON StreamingResponse (D-12)."""
    summary = _compute_summary_dict(results)
    duplicates = [
        {
            "original_smiles": r.get("original_smiles"),
            "standardized_inchikey": r.get("standardized_inchikey"),
            "rejection_reason": r.get("rejection_reason"),
        }
        for r in results
        if r.get("status") == "duplicate"
    ]

    payload = {
        "job_id": job_id,
        "summary": summary,
        "config": config_dict,
        "duplicates": duplicates,
        "results": results,
    }
    json_str = json.dumps(payload, indent=2, default=str)
    return StreamingResponse(
        iter([json_str]),
        media_type="application/json",
        headers={
            "Content-Disposition": f"attachment; filename=qsar_results_{job_id[:8]}.json"
        },
    )


# =============================================================================
# Utility helpers
# =============================================================================


def _validate_uuid(job_id: str) -> None:
    """Validate UUID format; raises HTTPException 400 if invalid."""
    import re

    uuid_re = re.compile(
        r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", re.I
    )
    if not uuid_re.match(job_id):
        raise HTTPException(
            status_code=400,
            detail="Invalid job ID format — must be a UUID",
        )


def _normalize_result(r: dict) -> dict:
    """Normalize a result dict to match QSARReadyResultSchema field names."""
    # Ensure steps have correct structure
    steps = r.get("steps", [])
    normalized_steps = []
    for s in steps:
        normalized_steps.append(
            {
                "step_name": s.get("step_name", ""),
                "step_index": s.get("step_index", 0),
                "enabled": s.get("enabled", True),
                "status": s.get("status", "skipped"),
                "before_smiles": s.get("before_smiles"),
                "after_smiles": s.get("after_smiles"),
                "detail": s.get("detail"),
            }
        )
    return {
        "original_smiles": r.get("original_smiles", ""),
        "original_inchikey": r.get("original_inchikey"),
        "curated_smiles": r.get("curated_smiles"),
        "standardized_inchikey": r.get("standardized_inchikey"),
        "inchikey_changed": r.get("inchikey_changed", False),
        "status": r.get("status", "error"),
        "rejection_reason": r.get("rejection_reason"),
        "steps": normalized_steps,
    }


def _compute_summary(results: list) -> QSARBatchSummary:
    """Compute summary statistics from a list of result dicts."""
    ok = sum(1 for r in results if r.get("status") == "ok")
    rejected = sum(1 for r in results if r.get("status") == "rejected")
    duplicate = sum(1 for r in results if r.get("status") == "duplicate")
    error = sum(1 for r in results if r.get("status") == "error")

    steps_applied_counts: dict = {}
    for r in results:
        for step in r.get("steps", []):
            if step.get("status") == "applied":
                name = step.get("step_name", "")
                steps_applied_counts[name] = steps_applied_counts.get(name, 0) + 1

    return QSARBatchSummary(
        total=len(results),
        ok=ok,
        rejected=rejected,
        duplicate=duplicate,
        error=error,
        steps_applied_counts=steps_applied_counts,
    )


def _compute_summary_dict(results: list) -> dict:
    """Compute summary statistics as a plain dict."""
    summary = _compute_summary(results)
    return summary.model_dump()
