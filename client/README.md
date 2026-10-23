# ChemAudit Python Client

Official Python client for the ChemAudit API - comprehensive chemical structure validation, standardization, and ML-readiness assessment.

## Installation

The client is not on PyPI yet. Install it from GitHub:

```bash
pip install "git+https://github.com/Kohulan/ChemAudit.git#subdirectory=client"
```

Or install from a local checkout:

```bash
cd client
pip install -e .
```

## Quick Start

```python
from chemaudit import ChemAuditClient

# Create client (an API key is optional)
client = ChemAuditClient(
    base_url="http://localhost:8000",
    api_key="your-api-key",  # Optional
)

# Validate a molecule
result = client.validate("CCO")
print(f"Validation Score: {result.overall_score}")
print(f"Issues: {len(result.issues)}")

# Screen for structural alerts
alerts = client.screen_alerts("c1ccccc1")
print(f"Total Alerts: {alerts.total_alerts}")

# Calculate ML-readiness score
scores = client.score("CCO")
print(f"ML-Readiness: {scores.ml_readiness['score']}")

# Standardize molecule
standard = client.standardize("[Na+].CC(=O)[O-]")
print(f"Standardized: {standard.result.standardized_smiles}")

# Close client when done
client.close()
```

## Features

### Synchronous-Only Design

This client is **intentionally synchronous-only** using `httpx.Client`. This design choice provides:

- **Simpler API** for most use cases (scripts, notebooks, CLI tools)
- **No async runtime complexity** for users
- **Batch operations with polling** that work well synchronously
- Easy integration with existing synchronous code

**For async usage:**
- Wrap calls with `asyncio.to_thread()` in async contexts
- Use `httpx.AsyncClient` directly for full async control

```python
# Async usage example
import asyncio
from chemaudit import ChemAuditClient

async def validate_async():
    client = ChemAuditClient()
    result = await asyncio.to_thread(client.validate, "CCO")
    return result
```

### Single Molecule Operations

```python
# Validate a molecule
result = client.validate(
    molecule="CCO",
    format="smiles",  # auto, smiles, inchi, mol, iupac
)

# Access results
print(f"Overall Score: {result.overall_score}")
for issue in result.issues:
    print(f"- {issue.check_name} ({issue.severity.value}): {issue.message}")

# Screen structural alerts
alerts = client.screen_alerts(
    molecule="c1ccccc1",
    catalogs=["PAINS", "BRENK"],  # Optional
)
print(f"Screened: {alerts.screened_catalogs}, critical: {alerts.has_critical}")
for alert in alerts.alerts:
    print(f"- {alert.pattern_name} [{alert.catalog_source}, {alert.severity}]: "
          f"{alert.description}")

# Calculate scores (each score block is a dict, or None if not computed)
scores = client.score("CCO")
ml_readiness = scores.ml_readiness["score"]
np_likeness = scores.np_likeness["score"]
scaffold = scores.scaffold["scaffold_smiles"]

# Request only some score blocks
scores = client.score("CCO", include=["ml_readiness", "druglikeness"])

# Standardize molecule
result = client.standardize(
    molecule="[Na+].CC(=O)[O-]",
    tautomer=False,  # Preserves E/Z stereochemistry
)
print(f"Standardized: {result.result.standardized_smiles}")
for step in result.result.steps_applied:
    print(f"- {step['step_name']}: {step['description']}")
```

### Batch Processing

```python
# Submit batch job
response = client.submit_batch(
    file_path="molecules.csv",
    smiles_column="smiles",
    include_standardization=False,
)
job_id = response.job_id

# Wait for completion (blocking with polling)
job = client.wait_for_batch(
    job_id,
    poll_interval=10.0,  # Status endpoint allows 10 requests/minute; keep >= 6 s
    timeout=3600.0,      # Max 1 hour
)
print(f"{job.status.value}: {job.processed}/{job.total}")

# Get results (paginated)
results = client.get_batch_results(
    job_id,
    page=1,
    page_size=100,
    score_min=50,  # Filter by score
)

# Iterate through all results
for item in client.iter_batch_results(job_id):
    if item.status == "success":
        score = item.validation["overall_score"]
        print(f"{item.smiles}: {score}")

# Export results (file name comes from the server unless output_path is given)
path = client.export_batch(
    job_id,
    format="excel",  # csv, excel, sdf, json, pdf
    score_min=50,
)
print(f"Exported to: {path}")

# Get statistics
stats = client.get_batch_stats(job_id)
print(f"Average Score: {stats.avg_validation_score}")
print(f"Alert Distribution: {stats.alert_summary}")
```

### Context Manager

```python
# Automatic cleanup
with ChemAuditClient() as client:
    result = client.validate("CCO")
    print(result.overall_score)
# Client automatically closed
```

## API Reference

### ChemAuditClient

**Constructor:**
```python
ChemAuditClient(
    base_url: str = "http://localhost:8000",
    api_key: Optional[str] = None,
    timeout: float = 30.0,
    max_retries: int = 3,
)
```

**Methods:**

- `validate(molecule, format="auto")` - Validate single molecule
- `screen_alerts(molecule, catalogs=None)` - Screen for structural alerts
- `score(molecule, include=None)` - Calculate ML-readiness and other scores
- `standardize(molecule, tautomer=False)` - Standardize molecule
- `submit_batch(file_path, smiles_column="smiles", include_standardization=False)` - Submit batch job
- `get_batch_status(job_id)` - Get batch job status (raises `BatchJobNotFoundError` if the job does not exist)
- `get_batch_results(job_id, page=1, page_size=100, score_min=None, score_max=None, status=None)` - Get batch results
- `get_batch_stats(job_id)` - Get batch statistics
- `cancel_batch(job_id)` - Cancel running batch job
- `wait_for_batch(job_id, poll_interval=10.0, timeout=3600.0)` - Wait for batch completion (the status endpoint allows 10 requests/minute, so keep `poll_interval` at 6 s or more)
- `iter_batch_results(job_id, page_size=100, **filters)` - Iterate all batch results
- `export_batch(job_id, format="csv", output_path=None, **filters)` - Export batch results
- `close()` - Close HTTP client

### Models

**ValidationResult:**
- `overall_score: int` - Validation score (0-100)
- `issues: List[CheckResult]` - Failed checks
- `all_checks: List[CheckResult]` - All check results
- `molecule_info: MoleculeInfo` - Molecule information

`CheckResult.severity` is a `Severity` enum with lowercase values: `critical`, `error`, `warning`, `info`, `pass`.

**AlertScreeningResult:**
- `alerts: List[AlertResult]` - Matches, each with `pattern_name`, `description`, `severity`, `catalog_source`, `smarts`, `matched_atoms`
- `screened_catalogs: List[str]`, `total_alerts: int`, `has_critical: bool`, `has_warning: bool`

**ScoreResult:**
- `ml_readiness`, `np_likeness`, `scaffold`, `druglikeness`, `safety_filters`, `admet`, `aggregator`, `consensus`, `lead_likeness` - Each a dict, or `None` if that block was not computed

**StandardizationResult:**
- `molecule_info: MoleculeInfo`, `execution_time_ms: int`
- `result: StandardizationDetails` - `standardized_smiles`, `steps_applied`, `checker_issues`, `excluded_fragments`, `stereo_comparison`, `mass_change_percent`

**BatchJob:**
- `status: BatchJobStatus` - Enum (`pending`, `processing`, `complete`, `failed`, `cancelled`); use `.value` for the string
- `progress`, `processed`, `total`, `eta_seconds`, `error_message`

**BatchResult:**
- `results: List[BatchResultItem]` - Paginated results
- `statistics: BatchStatistics` - Aggregate statistics
- `page: int`, `page_size: int`, `total_pages: int` - Pagination info

### Exceptions

- `ChemAuditError` - Base exception
- `APIError` - API returned error
- `RateLimitError` - Rate limit exceeded
- `AuthenticationError` - Invalid API key
- `ValidationError` - Request validation failed
- `BatchJobNotFoundError` - Job not found
- `TimeoutError` - Request or batch job timed out

## Rate Limits

Limits are set per endpoint (for example, 10 requests/minute for the batch status endpoint and 3 requests/minute for batch upload). Requests without an API key are counted per IP address. With an API key, requests are counted in a separate bucket for that key; the limits are the same.

On a 429 response the client waits for the number of seconds in the `Retry-After` header and retries, up to `max_retries` times. After that it raises `RateLimitError`.

## Error Handling

```python
from chemaudit import (
    ChemAuditClient,
    ValidationError,
    RateLimitError,
    AuthenticationError,
    APIError,
)

client = ChemAuditClient()

try:
    result = client.validate("invalid-smiles")
except ValidationError as e:
    print(f"Validation error: {e}")
except RateLimitError as e:
    print(f"Rate limit exceeded. Retry after {e.retry_after} seconds")
except AuthenticationError:
    print("Invalid API key")
except APIError as e:
    print(f"API error {e.status_code}: {e}")
```

## Development

Install development dependencies:

```bash
cd client
pip install -e ".[dev]"
```

Run tests:

```bash
pytest tests/ -v
```

## Examples

See the `examples/` directory for complete examples:

- `basic_usage.py` - Single molecule validation, alerts, scoring
- `batch_processing.py` - Batch file processing and export

## License

MIT License

## Links

- **Documentation:** https://github.com/Kohulan/ChemAudit#readme
- **Issue Tracker:** https://github.com/Kohulan/ChemAudit/issues
- **API Documentation:** http://localhost:8000/api/v1/docs (when server is running; http://localhost:8001/api/v1/docs with the development Docker setup)
