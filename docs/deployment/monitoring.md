# Monitoring Setup Guide

Monitor ChemStructVal performance and health with Prometheus and Grafana.

## Overview

The monitoring stack provides:

- **Prometheus**: Metrics collection, storage, and alerting
- **Grafana**: Visualization, dashboards, and notifications
- **Pre-built Dashboard**: ChemStructVal-specific metrics
- **Alerting**: Proactive issue detection

## Architecture

```
┌─────────────────────────────────────────────────────┐
│                   Monitoring Flow                   │
├─────────────────────────────────────────────────────┤
│                                                     │
│  [Backend]                                          │
│      │                                              │
│      ├── Exposes /metrics endpoint                 │
│      │   (prometheus-fastapi-instrumentator)       │
│      │                                              │
│      ▼                                              │
│  [Prometheus :9090]                                 │
│      │                                              │
│      ├── Scrapes metrics every 10s                 │
│      ├── Stores time-series data (30 days)         │
│      ├── Evaluates alert rules                     │
│      │                                              │
│      ▼                                              │
│  [Grafana :3001]                                    │
│      │                                              │
│      ├── Queries Prometheus                        │
│      ├── Renders dashboards                        │
│      └── Sends alerts (email, Slack, etc.)         │
│                                                     │
└─────────────────────────────────────────────────────┘
```

## Quick Start

### Accessing Dashboards

After deployment:

**Grafana UI:**
```
URL: http://localhost:3001
Username: admin
Password: (from GRAFANA_PASSWORD in .env)
```

**Prometheus UI:**
```
URL: http://localhost:9090
(Development only - not exposed in production)
```

### Default Dashboard

Navigate to: Dashboards → ChemStructVal Overview

## Prometheus Configuration

### prometheus.yml

Located at project root:

```yaml
global:
  scrape_interval: 10s      # How often to scrape targets
  evaluation_interval: 10s  # How often to evaluate rules

# Scrape targets
scrape_configs:
  - job_name: 'chemstructval-backend'
    metrics_path: '/metrics'
    static_configs:
      - targets: ['backend:8000']
        labels:
          service: 'chemstructval'
          environment: 'production'

  - job_name: 'prometheus'
    static_configs:
      - targets: ['localhost:9090']
```

### Storage and Retention

In `docker-compose.prod.yml`:

```yaml
prometheus:
  command:
    - '--config.file=/etc/prometheus/prometheus.yml'
    - '--storage.tsdb.path=/prometheus'
    - '--storage.tsdb.retention.time=30d'  # Keep 30 days of data
    - '--storage.tsdb.retention.size=10GB' # Max 10GB storage
```

Adjust based on your needs:
- High-traffic: 15-30 days
- Low-traffic: 60-90 days
- Compliance: 1-2 years (requires more storage)

## Available Metrics

### HTTP Metrics (Automatic)

Provided by `prometheus-fastapi-instrumentator`:

#### http_requests_total
Total HTTP requests by method, endpoint, and status.

```promql
# Request rate (requests per second)
rate(http_requests_total[5m])

# Requests by status code
sum by (status) (rate(http_requests_total[5m]))

# 5xx error rate
rate(http_requests_total{status=~"5.."}[5m])
```

#### http_request_duration_seconds
HTTP request latency histogram.

```promql
# P95 latency
histogram_quantile(0.95, rate(http_request_duration_seconds_bucket[5m]))

# P99 latency
histogram_quantile(0.99, rate(http_request_duration_seconds_bucket[5m]))

# Average latency
rate(http_request_duration_seconds_sum[5m]) / rate(http_request_duration_seconds_count[5m])
```

#### http_requests_inprogress
Currently processing requests.

```promql
# Active requests
http_requests_inprogress

# Average concurrent requests
avg_over_time(http_requests_inprogress[5m])
```

### Custom ChemStructVal Metrics

Located in `backend/app/monitoring/metrics.py`:

#### chemstructval_molecules_processed_total
Counter of molecules validated.

**Labels:**
- `status`: success | failed
- `source`: single | batch

```promql
# Molecules validated per second
rate(chemstructval_molecules_processed_total[5m])

# Success rate
rate(chemstructval_molecules_processed_total{status="success"}[5m])
  /
rate(chemstructval_molecules_processed_total[5m])
```

#### chemstructval_validation_duration_seconds
Histogram of validation latency.

```promql
# P95 validation time
histogram_quantile(0.95, rate(chemstructval_validation_duration_seconds_bucket[5m]))

# Average validation time
rate(chemstructval_validation_duration_seconds_sum[5m])
  /
rate(chemstructval_validation_duration_seconds_count[5m])
```

#### chemstructval_cache_hits_total
Counter of cache hits vs misses.

**Labels:**
- `result`: hit | miss

```promql
# Cache hit rate
rate(chemstructval_cache_hits_total{result="hit"}[5m])
  /
rate(chemstructval_cache_hits_total[5m])

# Cache effectiveness (should be >70%)
100 * (
  rate(chemstructval_cache_hits_total{result="hit"}[5m])
  /
  rate(chemstructval_cache_hits_total[5m])
)
```

#### chemstructval_batch_size
Histogram of batch job sizes.

```promql
# Average batch size
rate(chemstructval_batch_size_sum[5m]) / rate(chemstructval_batch_size_count[5m])

# P95 batch size
histogram_quantile(0.95, rate(chemstructval_batch_size_bucket[5m]))
```

#### chemstructval_active_batch_jobs
Gauge of currently processing batch jobs.

```promql
# Current batch jobs
chemstructval_active_batch_jobs

# Max concurrent batches
max_over_time(chemstructval_active_batch_jobs[1h])
```

#### chemstructval_alert_matches_total
Counter of structural alerts found.

**Labels:**
- `alert_type`: PAINS_A | PAINS_B | PAINS_C | BRENK

```promql
# PAINS alerts per second
rate(chemstructval_alert_matches_total{alert_type=~"PAINS.*"}[5m])

# Most common alert type
topk(3, rate(chemstructval_alert_matches_total[5m]))
```

#### chemstructval_external_api_requests_total
Counter of external API calls.

**Labels:**
- `service`: pubchem | chembl | decimer | coconut
- `status`: success | failure

```promql
# External API success rate
rate(chemstructval_external_api_requests_total{status="success"}[5m])
  /
rate(chemstructval_external_api_requests_total[5m])
```

## Grafana Dashboard

### Dashboard Layout

The pre-configured dashboard has 8 panels:

#### Row 1: Overview
1. **Request Rate**: Total requests per second
2. **P95 Latency**: 95th percentile response time
3. **Error Rate**: Percentage of 5xx responses

#### Row 2: Validation Metrics
4. **Molecules Processed**: Validation throughput
5. **Validation Duration**: P50/P95/P99 latency
6. **Cache Hit Rate**: Cache effectiveness

#### Row 3: Batch and Details
7. **Active Batch Jobs**: Current batch processing
8. **Structural Alerts**: PAINS/BRENK detection rate

### Dashboard Configuration

Located at `grafana/dashboards/chemstructval.json`:

```json
{
  "dashboard": {
    "title": "ChemStructVal Overview",
    "panels": [
      {
        "title": "Request Rate",
        "targets": [
          {
            "expr": "rate(http_requests_total[5m])"
          }
        ]
      }
      // ... more panels
    ]
  }
}
```

### Dashboard Provisioning

Grafana auto-loads dashboards from `grafana/provisioning/`:

```yaml
# grafana/provisioning/dashboards/dashboard.yml
apiVersion: 1

providers:
  - name: 'ChemStructVal'
    orgId: 1
    folder: ''
    type: file
    disableDeletion: false
    updateIntervalSeconds: 10
    allowUiUpdates: true
    options:
      path: /var/lib/grafana/dashboards
```

### Data Source Configuration

```yaml
# grafana/provisioning/datasources/datasource.yml
apiVersion: 1

datasources:
  - name: Prometheus
    type: prometheus
    access: proxy
    url: http://prometheus:9090
    isDefault: true
    editable: true
```

## Setting Up Alerts

### Alert Rules

Create `prometheus/alerts.yml`:

```yaml
groups:
  - name: chemstructval_alerts
    interval: 30s
    rules:
      # High error rate
      - alert: HighErrorRate
        expr: |
          (
            rate(http_requests_total{status=~"5.."}[5m])
            /
            rate(http_requests_total[5m])
          ) > 0.05
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High error rate detected"
          description: "Error rate is {{ $value | humanizePercentage }} (threshold: 5%)"

      # High latency
      - alert: HighLatency
        expr: |
          histogram_quantile(0.95, rate(http_request_duration_seconds_bucket[5m])) > 3
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High API latency"
          description: "P95 latency is {{ $value }}s (threshold: 3s)"

      # Cache hit rate low
      - alert: LowCacheHitRate
        expr: |
          (
            rate(chemstructval_cache_hits_total{result="hit"}[5m])
            /
            rate(chemstructval_cache_hits_total[5m])
          ) < 0.5
        for: 10m
        labels:
          severity: info
        annotations:
          summary: "Cache hit rate below 50%"
          description: "Cache effectiveness is {{ $value | humanizePercentage }}"

      # Too many batch jobs
      - alert: TooManyBatchJobs
        expr: chemstructval_active_batch_jobs > 10
        for: 10m
        labels:
          severity: warning
        annotations:
          summary: "Too many concurrent batch jobs"
          description: "{{ $value }} active batch jobs (threshold: 10)"

      # External API failures
      - alert: ExternalAPIFailures
        expr: |
          (
            rate(chemstructval_external_api_requests_total{status="failure"}[5m])
            /
            rate(chemstructval_external_api_requests_total[5m])
          ) > 0.2
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High external API failure rate"
          description: "{{ $labels.service }} failure rate: {{ $value | humanizePercentage }}"

      # Service down
      - alert: ServiceDown
        expr: up{job="chemstructval-backend"} == 0
        for: 1m
        labels:
          severity: critical
        annotations:
          summary: "ChemStructVal backend is down"
          description: "Backend service has been down for 1 minute"
```

Update `prometheus.yml`:

```yaml
global:
  scrape_interval: 10s
  evaluation_interval: 10s

rule_files:
  - '/etc/prometheus/alerts.yml'

scrape_configs:
  # ... existing config
```

Update `docker-compose.prod.yml`:

```yaml
prometheus:
  volumes:
    - ./prometheus.yml:/etc/prometheus/prometheus.yml:ro
    - ./prometheus/alerts.yml:/etc/prometheus/alerts.yml:ro
```

### Grafana Alert Channels

#### Email Notifications

In Grafana UI:

1. Navigate to: Alerting → Contact points
2. Click **New contact point**
3. Configure:
   ```
   Name: Email
   Type: Email
   Addresses: ops@yourcompany.com
   ```
4. Test and save

#### Slack Notifications

1. Create Slack incoming webhook: https://api.slack.com/messaging/webhooks
2. In Grafana: Alerting → Contact points
3. Configure:
   ```
   Name: Slack
   Type: Slack
   URL: https://hooks.slack.com/services/YOUR/WEBHOOK/URL
   Channel: #chemstructval-alerts
   ```

#### Alert Rules in Grafana

1. Navigate to panel → Edit
2. Click **Alert** tab
3. Create rule:
   ```
   Alert name: High Error Rate
   Evaluate every: 1m
   For: 5m

   Condition:
   WHEN avg() OF query(A, 5m, now) IS ABOVE 0.05

   Notification: Email, Slack
   Message: Error rate is {{ $value }}%
   ```

## Monitoring Best Practices

### Dashboard Design

✅ **Do:**
- Use consistent time ranges
- Include P50, P95, P99 for latency
- Show rate over time (not raw counters)
- Use appropriate visualization (graph, gauge, stat)
- Group related metrics

❌ **Don't:**
- Overcrowd dashboards (max 12 panels)
- Use instant queries without rate()
- Mix different time ranges
- Use pie charts for time-series

### Query Performance

✅ **Do:**
- Use recording rules for complex queries
- Aggregate before querying
- Use appropriate time ranges
- Limit cardinality of labels

```yaml
# Recording rules for expensive queries
groups:
  - name: chemstructval_recordings
    interval: 10s
    rules:
      - record: job:validation_success_rate:5m
        expr: |
          rate(chemstructval_molecules_processed_total{status="success"}[5m])
          /
          rate(chemstructval_molecules_processed_total[5m])
```

❌ **Don't:**
- Query long time ranges (>24h) at high resolution
- Use regex unnecessarily
- Create high-cardinality labels
- Run expensive queries in alerts

### Alert Fatigue Prevention

✅ **Do:**
- Set appropriate thresholds
- Use `for` clause to avoid flapping
- Group related alerts
- Include actionable information
- Test alerts before deploying

❌ **Don't:**
- Alert on everything
- Set unrealistic thresholds
- Ignore repeated alerts
- Send alerts without context

## Performance Tuning

### Prometheus Resource Usage

**Memory:**
```
Memory = (samples/sec) × (retention days) × (bytes/sample)
Typical: 2-8 GB for moderate load
```

Reduce memory:
```yaml
prometheus:
  command:
    - '--storage.tsdb.retention.time=15d'  # Shorter retention
    - '--storage.tsdb.min-block-duration=2h'  # Larger blocks
```

**Disk:**
```
Disk usage = memory usage × 1.5 to 2
Typical: 5-20 GB for 30 days retention
```

### Grafana Optimization

**Dashboard queries:**
- Use $__interval for dynamic step
- Enable min interval: 10s
- Use caching where appropriate

**Example:**
```promql
# Good: Uses dashboard interval
rate(http_requests_total[$__interval])

# Bad: Fixed interval doesn't scale
rate(http_requests_total[1m])
```

## Troubleshooting

### Grafana shows "No Data"

**Check:**
1. Prometheus is scraping:
   ```
   Visit: http://localhost:9090/targets
   Status should be: UP
   ```

2. Metrics exist in Prometheus:
   ```
   Visit: http://localhost:9090/graph
   Query: http_requests_total
   Should show data
   ```

3. Grafana data source configured:
   ```
   Settings → Data Sources → Prometheus
   Test: Should succeed
   ```

4. Backend exposes metrics:
   ```bash
   curl http://localhost:8000/metrics
   # Should return Prometheus metrics
   ```

### Prometheus Target Down

**Check backend health:**
```bash
# Container running?
docker-compose ps backend

# Health endpoint accessible?
curl http://localhost:8000/api/v1/health

# Metrics endpoint accessible?
curl http://localhost:8000/metrics
```

**Check network:**
```bash
# From Prometheus container
docker exec chemstructval-prometheus wget -O- http://backend:8000/metrics
```

**Check configuration:**
```bash
# Verify prometheus.yml syntax
docker exec chemstructval-prometheus promtool check config /etc/prometheus/prometheus.yml
```

### High Memory Usage

**Prometheus:**
```bash
# Check TSDB size
docker exec chemstructval-prometheus du -sh /prometheus

# Check cardinality
curl http://localhost:9090/api/v1/status/tsdb

# Reduce retention
# Edit docker-compose.prod.yml:
# - '--storage.tsdb.retention.time=15d'
```

**Grafana:**
```bash
# Clear query cache
docker exec chemstructval-grafana rm -rf /var/lib/grafana/cache/*
docker-compose restart grafana
```

### Missing Metrics

**Check metric registration:**
```python
# In backend code
from app.monitoring.metrics import (
    molecules_processed,
    validation_duration,
)

# Ensure metrics are incremented
molecules_processed.labels(status="success", source="single").inc()
```

**Check Prometheus scrape:**
```bash
# View raw metrics
curl http://localhost:8000/metrics | grep chemstructval
```

## Advanced Configuration

### Multi-Environment Monitoring

```yaml
# prometheus.yml
scrape_configs:
  - job_name: 'chemstructval-prod'
    static_configs:
      - targets: ['prod-backend:8000']
        labels:
          environment: 'production'

  - job_name: 'chemstructval-staging'
    static_configs:
      - targets: ['staging-backend:8000']
        labels:
          environment: 'staging'
```

### Federation

For multiple Prometheus instances:

```yaml
# Central Prometheus
scrape_configs:
  - job_name: 'federate'
    honor_labels: true
    metrics_path: '/federate'
    params:
      'match[]':
        - '{job="chemstructval-backend"}'
    static_configs:
      - targets:
        - 'prom-region1:9090'
        - 'prom-region2:9090'
```

### Custom Exporters

Create custom metrics exporter:

```python
# backend/app/monitoring/custom_exporter.py
from prometheus_client import Gauge

# Custom gauge
molecules_in_db = Gauge(
    'chemstructval_molecules_in_database',
    'Total molecules stored in database'
)

# Update periodically
async def update_database_metrics():
    count = await get_molecule_count()
    molecules_in_db.set(count)
```

## Resources

- [Prometheus Documentation](https://prometheus.io/docs/)
- [Grafana Documentation](https://grafana.com/docs/)
- [PromQL Cheat Sheet](https://promlabs.com/promql-cheat-sheet/)
- [Grafana Dashboard Examples](https://grafana.com/grafana/dashboards/)
- [Prometheus Best Practices](https://prometheus.io/docs/practices/)

## Getting Help

- Check [FAQ](../troubleshooting/faq.md)
- Review Prometheus logs: `docker-compose logs prometheus`
- Review Grafana logs: `docker-compose logs grafana`
- Test queries in Prometheus UI
- Verify metrics endpoint
- Open GitHub issue with details
