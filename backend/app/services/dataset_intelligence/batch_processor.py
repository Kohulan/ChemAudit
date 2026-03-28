"""
Dataset Intelligence Batch Processor (Phase 12).

Provides a Celery task for async dataset health auditing, contradictory label
detection, and curation report generation.

Progress is published to the Redis channel ``dataset:{job_id}`` -- distinct from
``batch:{job_id}`` (regular batch), ``qsar:{job_id}`` (QSAR pipeline), and
``genchem:{job_id}`` (GenChem filter) channels to avoid channel collision
(Pitfall 4 from RESEARCH.md).

Redis key structure:
- ``dataset:meta:{job_id}``    -- hash with status, progress, current_stage, filename
- ``dataset:results:{job_id}`` -- JSON blob with full audit results (24h TTL)

File handling uses temp file path on disk (NOT base64) per Pitfall 6 -- the
upload route writes to a NamedTemporaryFile and passes the path as a string.
The task deletes the temp file in a finally block after processing.
"""

import hashlib
import json
import logging
import os

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi as rdkit_inchi

from app.celery_app import celery_app

logger = logging.getLogger(__name__)

# Heuristic column names for SMILES detection (case-insensitive)
_SMILES_COLUMN_NAMES = {
    "smiles", "smi", "canonical_smiles", "isomeric_smiles",
    "structure", "input_smiles", "mol_smiles",
}


def _detect_smiles_column(columns: list[str]) -> str | None:
    """Detect the SMILES column name from a list of column names.

    Checks exact case-insensitive match first, then substring match.

    Args:
        columns: List of column names from the dataset.

    Returns:
        Matched column name, or None if no match found.
    """
    for col in columns:
        if col.lower() in _SMILES_COLUMN_NAMES:
            return col
    # Fallback: substring match
    for col in columns:
        col_lower = col.lower()
        if "smiles" in col_lower or "smi" == col_lower:
            return col
    return None


def _parse_csv_full(file_path: str, smiles_column: str | None) -> tuple[list[dict], list[str]]:
    """Parse a CSV file into molecule dicts with all properties.

    Unlike file_parser.parse_csv which only reads the SMILES column,
    this reads ALL columns so property distributions and contradictory
    label detection have access to the full data (Pitfall 9).

    Args:
        file_path: Path to the CSV file on disk.
        smiles_column: Explicit SMILES column name, or None for auto-detect.

    Returns:
        Tuple of (molecule list, column names).
    """
    df = pd.read_csv(file_path, dtype=str, na_filter=False)
    columns = df.columns.tolist()

    # Determine SMILES column
    smi_col = smiles_column
    if smi_col is None:
        smi_col = _detect_smiles_column(columns)
    if smi_col is None:
        raise ValueError(
            "Could not detect SMILES column. "
            "Please specify the smiles_column parameter."
        )

    # Verify column exists (case-insensitive)
    actual_col = None
    for col in columns:
        if col.lower() == smi_col.lower():
            actual_col = col
            break
    if actual_col is None:
        raise ValueError(f"SMILES column '{smi_col}' not found in CSV")

    molecules: list[dict] = []
    for idx, row in df.iterrows():
        smiles = str(row[actual_col]).strip()
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        ik = None
        if mol is not None and mol.GetNumAtoms() > 0:
            try:
                ik = rdkit_inchi.MolToInchiKey(
                    rdkit_inchi.MolFromInchi(rdkit_inchi.MolToInchi(mol))[0]
                )
            except Exception:
                try:
                    ik = rdkit_inchi.InchiToInchiKey(rdkit_inchi.MolToInchi(mol))
                except Exception:
                    pass

        properties = {col: str(row[col]) for col in columns if col != actual_col}
        molecules.append({
            "index": int(idx),
            "smiles": smiles,
            "mol": mol,
            "inchikey": ik,
            "properties": properties,
        })

    return molecules, columns


def _parse_sdf_full(file_path: str) -> tuple[list[dict], list[str]]:
    """Parse an SDF file into molecule dicts with all properties.

    Args:
        file_path: Path to the SDF file on disk.

    Returns:
        Tuple of (molecule list, column names).
    """
    with open(file_path, "rb") as f:
        content = f.read()

    from app.services.batch.file_parser import parse_sdf

    parsed = parse_sdf(content)
    all_prop_names: set[str] = set()

    molecules: list[dict] = []
    for mol_data in parsed:
        smiles = mol_data.smiles
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        ik = None
        if mol is not None and mol.GetNumAtoms() > 0:
            try:
                ik = rdkit_inchi.InchiToInchiKey(rdkit_inchi.MolToInchi(mol))
            except Exception:
                pass

        props = dict(mol_data.properties)
        all_prop_names.update(props.keys())
        molecules.append({
            "index": mol_data.index,
            "smiles": smiles,
            "mol": mol,
            "inchikey": ik,
            "properties": props,
        })

    columns = ["SMILES"] + sorted(all_prop_names)
    return molecules, columns


@celery_app.task(
    bind=True,
    queue="default",
    name="app.services.dataset_intelligence.batch_processor.process_dataset_audit",
)
def process_dataset_audit(
    self,
    file_path: str,
    filename: str,
    file_type: str,
    job_id: str,
    options: dict,
) -> dict:
    """Process a dataset audit job asynchronously.

    Runs the full health audit pipeline in 6 stages:
    1. Parse file (0-10%)
    2. Detect numeric columns (10-15%)
    3. Health audit with 5 sub-scores (15-70%)
    4. Contradictory label detection (70-85%)
    5. Build curation report + curated CSV (85-95%)
    6. Store results in Redis (95-100%)

    Progress is published to the ``dataset:{job_id}`` Redis channel for
    real-time WebSocket updates via /ws/dataset/{job_id}.

    Args:
        file_path: Path to the uploaded temp file on disk (NOT base64).
        filename: Original filename for metadata.
        file_type: File type: 'csv' or 'sdf'.
        job_id: Unique job identifier.
        options: Processing options dict with optional keys:
            smiles_column (str|None), activity_column (str|None).

    Returns:
        Dict with audit results summary.
    """
    import redis as sync_redis

    from app.core.config import settings

    r = sync_redis.from_url(settings.REDIS_URL, decode_responses=True)

    def publish_progress(stage: str, pct: float) -> None:
        """Publish progress to Redis channel and update metadata."""
        msg = json.dumps({
            "job_id": job_id,
            "status": "processing",
            "progress": round(pct, 1),
            "current_stage": stage,
        })
        try:
            r.publish(f"dataset:{job_id}", msg)
            r.hset(f"dataset:meta:{job_id}", mapping={
                "status": "processing",
                "progress": str(round(pct, 1)),
                "current_stage": stage,
            })
        except Exception:
            logger.warning("Failed to publish progress for dataset job %s", job_id)

    try:
        # ---- Stage 1 (0-10%): Parse file ----
        publish_progress("parsing", 0.0)

        smiles_column = options.get("smiles_column")
        activity_column = options.get("activity_column")

        if file_type == "sdf":
            molecules, columns = _parse_sdf_full(file_path)
        else:
            molecules, columns = _parse_csv_full(file_path, smiles_column)

        # Compute file metadata
        with open(file_path, "rb") as f:
            file_bytes = f.read()
        file_size_bytes = len(file_bytes)
        sha256_hash = hashlib.sha256(file_bytes).hexdigest()

        publish_progress("parsing", 10.0)

        # ---- Stage 2 (10-15%): Detect numeric columns ----
        publish_progress("detecting_columns", 10.0)

        from app.services.dataset_intelligence.contradictory_labels import (
            detect_numeric_columns,
        )

        # Build sample values for numeric column detection
        sample_values: dict[str, list] = {}
        sample_size = min(len(molecules), 100)
        for col in columns:
            vals: list = []
            for mol_dict in molecules[:sample_size]:
                v = mol_dict.get("properties", {}).get(col)
                if v is not None:
                    vals.append(v)
            sample_values[col] = vals

        numeric_columns = detect_numeric_columns(columns, sample_values)

        # Auto-detect activity column if not specified
        if activity_column is None:
            for nc in numeric_columns:
                if nc["priority"] <= 2:
                    activity_column = nc["name"]
                    break

        publish_progress("detecting_columns", 15.0)

        # ---- Stage 3 (15-70%): Health audit ----
        from app.services.dataset_intelligence.health_audit import compute_health_score

        def health_progress_cb(stage_name: str, frac: float) -> None:
            """Map health audit internal progress to 15-70% range."""
            mapped_pct = 15.0 + frac * 55.0
            publish_progress(f"health_{stage_name}", mapped_pct)

        health_result = compute_health_score(
            molecules,
            progress_callback=health_progress_cb,
        )

        publish_progress("health_complete", 70.0)

        # ---- Stage 4 (70-85%): Contradictory labels ----
        publish_progress("contradictory_labels", 70.0)

        from app.services.dataset_intelligence.contradictory_labels import (
            detect_contradictory_labels,
        )

        contradictions: list[dict] = []
        if activity_column:
            contradictions = detect_contradictory_labels(
                molecules, activity_column
            )

        publish_progress("contradictory_labels", 85.0)

        # ---- Stage 5 (85-95%): Build curation report + curated CSV ----
        publish_progress("building_report", 85.0)

        from app.services.dataset_intelligence.curation_report import (
            build_curation_report,
            build_curated_csv_rows,
        )

        file_metadata = {
            "filename": filename,
            "format": file_type,
            "row_count": len(molecules),
            "file_size_bytes": file_size_bytes,
            "sha256_hash": sha256_hash,
        }

        curation_report = build_curation_report(
            health_result, file_metadata, contradictions
        )
        curated_rows = build_curated_csv_rows(molecules, health_result)

        publish_progress("building_report", 95.0)

        # ---- Stage 6 (95-100%): Store results in Redis ----
        publish_progress("storing_results", 95.0)

        # Build sub-scores list for the response schema
        sub_scores = [
            {
                "name": "parsability",
                "score": health_result.parsability_score,
                "weight": health_result.weights.get("parsability", 0.25),
                "count": health_result.parse_failures,
                "total": health_result.molecule_count,
            },
            {
                "name": "stereo",
                "score": health_result.stereo_score,
                "weight": health_result.weights.get("stereo", 0.15),
                "count": health_result.stereo_undefined_count,
                "total": health_result.molecule_count,
            },
            {
                "name": "uniqueness",
                "score": health_result.uniqueness_score,
                "weight": health_result.weights.get("uniqueness", 0.20),
                "count": health_result.duplicate_count,
                "total": health_result.molecule_count,
            },
            {
                "name": "alerts",
                "score": health_result.alert_score,
                "weight": health_result.weights.get("alerts", 0.20),
                "count": health_result.alert_hit_count,
                "total": health_result.molecule_count,
            },
            {
                "name": "std_consistency",
                "score": health_result.std_consistency_score,
                "weight": health_result.weights.get("std_consistency", 0.20),
                "count": health_result.std_disagreement_count,
                "total": health_result.std_sample_size,
            },
        ]

        # Strip mol objects before serialization (not JSON-serializable)
        serializable_molecules = []
        for mol_dict in molecules:
            serializable_molecules.append({
                "index": mol_dict["index"],
                "smiles": mol_dict["smiles"],
                "inchikey": mol_dict.get("inchikey"),
                "properties": mol_dict.get("properties", {}),
            })

        result_data = {
            "job_id": job_id,
            "status": "complete",
            "health_audit": {
                "overall_score": health_result.overall_score,
                "sub_scores": sub_scores,
                "weights": health_result.weights,
                "molecule_count": health_result.molecule_count,
                "issues": health_result.issues,
                "property_distributions": health_result.property_distributions,
                "std_pipeline_comparison": health_result.std_pipeline_comparison,
                "std_sample_size": health_result.std_sample_size,
                "dedup_groups": health_result.dedup_groups,
            },
            "contradictions": contradictions,
            "numeric_columns": numeric_columns,
            "curation_report": curation_report,
            "curated_csv_rows": curated_rows,
            "curated_csv_available": bool(curated_rows),
            "molecules": serializable_molecules,
        }

        r.set(
            f"dataset:results:{job_id}",
            json.dumps(result_data),
            ex=settings.BATCH_RESULT_TTL,
        )

        # Update metadata to complete
        r.hset(f"dataset:meta:{job_id}", mapping={
            "status": "complete",
            "progress": "100",
            "current_stage": "",
        })
        r.expire(f"dataset:meta:{job_id}", settings.BATCH_RESULT_TTL)

        # Publish final progress
        r.publish(
            f"dataset:{job_id}",
            json.dumps({
                "job_id": job_id,
                "status": "complete",
                "progress": 100,
                "current_stage": None,
            }),
        )

        publish_progress("complete", 100.0)

        logger.info(
            "Dataset audit complete for job %s: %d molecules, score=%.1f",
            job_id, health_result.molecule_count, health_result.overall_score,
        )

        return {"job_id": job_id, "status": "complete"}

    except Exception as exc:
        logger.error("Dataset audit failed for job %s: %s", job_id, exc, exc_info=True)

        # Update metadata to error
        try:
            r.hset(f"dataset:meta:{job_id}", mapping={
                "status": "error",
                "progress": "0",
                "current_stage": "",
                "error": str(exc)[:500],
            })
            r.expire(f"dataset:meta:{job_id}", 3600)
            r.publish(
                f"dataset:{job_id}",
                json.dumps({
                    "job_id": job_id,
                    "status": "error",
                    "error": str(exc)[:500],
                }),
            )
        except Exception:
            logger.warning("Failed to publish error status for dataset job %s", job_id)

        raise

    finally:
        # Cleanup: delete temp file
        try:
            if os.path.exists(file_path):
                os.unlink(file_path)
                logger.debug("Deleted temp file: %s", file_path)
        except OSError as e:
            logger.warning("Failed to delete temp file %s: %s", file_path, e)
