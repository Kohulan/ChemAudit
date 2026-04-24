"""
Prometheus metrics for ChemAudit monitoring.

Provides custom metrics for validation performance, cache efficiency,
batch processing, and alert matching. These metrics enable:

- Performance monitoring (latency, throughput)
- Resource utilization (batch sizes, active jobs)
- Cache efficiency (hit/miss rates)
- Alert pattern analysis (match counts by type)

Usage:
    from app.core.metrics import (
        VALIDATION_DURATION,
        MOLECULES_PROCESSED,
    )

    # Record validation timing
    with VALIDATION_DURATION.time():
        result = validate(mol)

    # Increment molecules processed counter
    MOLECULES_PROCESSED.labels(status="success").inc()
"""

from prometheus_client import Counter, Gauge, Histogram, Info

from app.core.config import settings

# Application info metric
APP_INFO = Info(
    "chemaudit",
    "ChemAudit application information",
)
APP_INFO.info(
    {
        "version": settings.APP_VERSION,
        "app_name": settings.APP_NAME,
    }
)

# Validation timing histogram
# Buckets chosen for typical validation times: 10ms to 30s
VALIDATION_DURATION = Histogram(
    "chemaudit_validation_duration_seconds",
    "Time spent validating molecules",
    ["validation_type"],  # single, batch
    buckets=(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0),
)

# Molecules processed counter with status label
MOLECULES_PROCESSED = Counter(
    "chemaudit_molecules_processed_total",
    "Total number of molecules processed",
    ["status"],  # success, error, invalid
)

# Batch size histogram
# Buckets for typical batch sizes: 1 to 10000
BATCH_SIZE = Histogram(
    "chemaudit_batch_size",
    "Distribution of batch job sizes",
    buckets=(1, 10, 50, 100, 250, 500, 1000, 2500, 5000, 10000),
)

# Cache metrics
CACHE_HITS = Counter(
    "chemaudit_cache_hits_total",
    "Total number of cache hits",
)

CACHE_MISSES = Counter(
    "chemaudit_cache_misses_total",
    "Total number of cache misses",
)

# Active batch jobs gauge
ACTIVE_BATCH_JOBS = Gauge(
    "chemaudit_active_batch_jobs",
    "Number of currently active batch processing jobs",
)

# Alert matches counter by alert type
ALERT_MATCHES = Counter(
    "chemaudit_alert_matches_total",
    "Total number of structural alert matches",
    ["alert_type"],  # PAINS, BRENK, NIH, ZINC
)

# Standardization counter
STANDARDIZATIONS_PERFORMED = Counter(
    "chemaudit_standardizations_total",
    "Total number of molecule standardizations performed",
    ["status"],  # success, error
)

# External API calls
EXTERNAL_API_CALLS = Counter(
    "chemaudit_external_api_calls_total",
    "Total number of external API calls",
    ["api", "status"],  # api: pubchem, chembl, coconut; status: success, error
)

EXTERNAL_API_DURATION = Histogram(
    "chemaudit_external_api_duration_seconds",
    "Time spent calling external APIs",
    ["api"],
    buckets=(0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0),
)
