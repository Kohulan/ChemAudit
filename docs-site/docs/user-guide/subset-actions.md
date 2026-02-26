---
sidebar_position: 4
title: Subset Actions & Sharing
description: Work with molecule subsets, share results via permalinks, and receive notifications
---

# Subset Actions & Sharing

Work with subsets of your batch results without re-uploading files. Select molecules, perform actions, share results via permalinks, and receive notifications when processing completes.

## Selecting Molecules

Select molecules from the batch results table using row checkboxes:

- **Individual selection**: Click the checkbox on any row
- **Page select-all**: Click the header checkbox to select all molecules on the current page
- **Cross-page selection**: Selections persist across pages — navigate freely without losing selections

The selection count appears in a floating **Actions** button at the bottom of the results table.

:::tip Selection from Charts
You can also select molecules by clicking on chart elements in the [Batch Analytics](/docs/user-guide/batch-analytics) visualizations. Selections sync between the table and all charts.
:::

## Subset Action Panel

Click the **Actions** button to open the slide-over panel with two tabs:

### Validation Tab

Displays the existing validation results for your selected molecules without making any API calls. Includes:

- Summary strip showing pass/warn/fail counts for the selection
- Per-molecule validation cards with score, key metrics, and issues
- **Download CSV** button to export the subset's validation data

### Scoring Tab

Re-score selected molecules against a different scoring profile — computed inline (no background job required):

1. Select a profile from the dropdown (shows all presets + custom profiles)
2. Results appear immediately with per-molecule profile scores
3. Each scored card shows the score bar and per-property desirability pills
4. A mini distribution bar shows the pass/moderate/fail breakdown
5. **Download CSV** button exports scores and per-property breakdowns

### Re-validate Subset

Create a new batch job containing only the selected molecules:

```bash
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/subset/revalidate \
  -H "Content-Type: application/json" \
  -d '{"indices": [0, 5, 12, 23, 42]}'
```

**Response:**

```json
{
  "new_job_id": "a1b2c3d4-...",
  "source_job_id": "550e8400-...",
  "molecule_count": 5,
  "action": "revalidate"
}
```

### Re-score Subset

Re-score selected molecules with a different profile as a new batch job:

```bash
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/subset/rescore \
  -H "Content-Type: application/json" \
  -d '{"indices": [0, 5, 12, 23, 42], "profile_id": 3}'
```

### Inline Score (No Background Job)

For quick re-scoring without creating a new batch job, use the inline endpoint:

```bash
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/subset/score-inline \
  -H "Content-Type: application/json" \
  -d '{"indices": [0, 5, 12], "profile_id": 4}'
```

**Response:**

```json
{
  "profile_name": "CNS-penetrant",
  "profile_id": 4,
  "molecules": [
    {
      "index": 0,
      "name": "Aspirin",
      "smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "profile": {
        "profile_id": 4,
        "profile_name": "CNS-penetrant",
        "score": 72.5,
        "properties": {
          "mw": {"value": 180.16, "min": 0, "max": 400, "in_range": true, "desirability": 1.0},
          "logp": {"value": 1.2, "min": 1, "max": 5, "in_range": true, "desirability": 1.0},
          "tpsa": {"value": 63.6, "min": 0, "max": 90, "in_range": true, "desirability": 1.0},
          "hbd": {"value": 1, "min": 0, "max": 2, "in_range": true, "desirability": 1.0}
        }
      }
    }
  ]
}
```

### Export Subset

Export only the selected molecules in any supported format:

```bash
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/subset/export \
  -H "Content-Type: application/json" \
  -d '{"indices": [0, 5, 12, 23, 42], "format": "csv"}' \
  -o subset.csv
```

Supports all export formats: `csv`, `excel`, `sdf`, `json`, `pdf`, `fingerprint`, `dedup`, `scaffold`, `property_matrix`.

## Batch to Single Navigation

Each row in the batch results table has an **Open in Single Validation** button. Clicking it navigates to the single validation page with the molecule pre-loaded, showing:

- A back-to-batch navigation bar with the molecule name and index
- Full single validation view with all tabs (Validation, Alerts, Scoring, Standardization, Database Lookup, Scoring Profiles)

This lets you explore a batch molecule in full detail, then return to the batch results.

## Sharing & Permalinks

### Batch Report Permalinks

Share your batch results with others using short-lived permalinks:

```bash
# Create a permalink
curl -X POST http://localhost:8001/api/v1/permalinks \
  -H "Content-Type: application/json" \
  -d '{"job_id": "550e8400-..."}'
```

**Response:**

```json
{
  "short_id": "a1b2c3d4e5f6g7h8",
  "job_id": "550e8400-...",
  "url": "https://chemaudit.example.com/report/a1b2c3d4e5f6g7h8",
  "created_at": "2026-02-26T10:30:00Z",
  "expires_at": "2026-03-28T10:30:00Z"
}
```

The URL is automatically copied to your clipboard in the web interface. Permalinks expire after 30 days.

### Resolving Permalinks

```bash
# Resolve a permalink
curl http://localhost:8001/api/v1/report/a1b2c3d4e5f6g7h8
```

Returns the job data and optional snapshot. Returns **410 Gone** if the permalink has expired.

### Single Molecule Permalinks

Share a single molecule validation via URL parameter — no database storage, no expiry:

```
https://chemaudit.example.com/?smiles=CC(=O)Oc1ccccc1C(=O)O
```

This triggers a fresh validation when the link is opened.

## Notifications

ChemAudit can notify you when batch processing completes via email or webhook.

### Email Notifications

Provide an email address when uploading a batch to receive a completion notification:

```bash
curl -X POST http://localhost:8001/api/v1/batch/upload \
  -F "file=@molecules.sdf" \
  -F "notification_email=researcher@example.com"
```

The email includes molecule counts, pass/fail statistics, and a direct link to the batch results.

:::info SMTP Configuration
Email notifications require SMTP server configuration. See [Configuration](/docs/getting-started/configuration) for the required environment variables.
:::

### Webhook Notifications

Configure webhook notifications to integrate with your pipeline:

**Environment variables** (set in `.env` or `docker-compose.yml`):

| Variable | Description |
|----------|-------------|
| `WEBHOOK_URL` | HTTP endpoint to receive POST notifications |
| `WEBHOOK_SECRET` | Secret key for HMAC-SHA256 signature verification |

**Webhook payload:**

```json
{
  "event": "batch_complete",
  "job_id": "550e8400-...",
  "status": "complete",
  "molecule_count": 1000,
  "pass_count": 850,
  "fail_count": 150,
  "avg_score": 78.5,
  "report_url": "https://chemaudit.example.com/batch/550e8400-...",
  "timestamp": "2026-02-26T10:35:00Z"
}
```

**Signature verification:**

The webhook request includes an `X-ChemAudit-Signature` header containing an HMAC-SHA256 signature of the request body:

```
X-ChemAudit-Signature: sha256=<hex_digest>
```

Verify the signature by computing `HMAC-SHA256(WEBHOOK_SECRET, request_body)` and comparing.

**Retry policy:** Failed webhook deliveries are retried up to 3 times with exponential backoff.

## Next Steps

- **[Batch Analytics](/docs/user-guide/batch-analytics)** — Explore analytics that drive subset selection
- **[Exporting Results](/docs/user-guide/exporting-results)** — Export full or subset results
- **[Scoring Profiles](/docs/user-guide/scoring/profiles)** — Create profiles for inline scoring
- **[Batch Processing](/docs/user-guide/batch-processing)** — Upload and process datasets
