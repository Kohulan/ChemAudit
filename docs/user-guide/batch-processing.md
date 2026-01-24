# Batch Processing Guide

Process thousands of molecules at once with batch validation.

## Overview

Batch processing allows you to:
- Validate entire chemical libraries
- Process up to 10,000 molecules in under 5 minutes
- Track progress in real-time
- Export results in multiple formats
- Handle individual failures gracefully

## Supported File Formats

| Format | Extension | Max Size | Notes |
|--------|-----------|----------|-------|
| SDF | .sdf, .sd | 100MB | Multi-molecule, includes coordinates |
| CSV | .csv | 50MB | Requires SMILES column mapping |

### SDF Format

Standard Structure Data File with multiple MOL blocks:

```
Example.sdf
-----------
Molecule 1
...
$$$$
Molecule 2
...
$$$$
```

**Advantages:**
- Preserves 3D coordinates
- Includes metadata fields
- Standard in chemistry

**Considerations:**
- Larger file size than CSV
- Slower parsing than SMILES

### CSV Format

Comma-separated values with a SMILES column:

```csv
id,name,smiles,activity
1,Aspirin,CC(=O)Oc1ccccc1C(=O)O,active
2,Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C,active
3,Ethanol,CCO,inactive
```

**Requirements:**
- Must have a column with SMILES strings
- Header row recommended
- UTF-8 encoding

**Advantages:**
- Compact file size
- Fast processing
- Easy to create from databases

## Uploading a Batch

### Step 1: Navigate to Batch Processing

1. Click the **Batch Processing** tab in the navigation
2. You'll see the batch upload interface

### Step 2: Select File

Two ways to upload:

**Drag and Drop:**
1. Drag your SDF or CSV file to the upload area
2. Drop when the area highlights

**File Browser:**
1. Click the upload area
2. Select file from your computer

### Step 3: Configure (CSV only)

For CSV files, you'll see a preview:

1. **Select SMILES column** from dropdown
2. **Preview** shows first 5 molecules
3. **Verify** structures are parsed correctly

Example preview:
```
Row 1: CCO          → Ethanol (3 atoms)
Row 2: c1ccccc1     → Benzene (6 atoms)
Row 3: CC(=O)O      → Acetic acid (4 atoms)
```

### Step 4: Start Processing

1. Click **Start Processing**
2. Job submits to queue
3. Processing begins immediately

## Progress Tracking

During processing, you'll see:

### Progress Bar
```
Processing: [████████████--------] 65% (650/1000)
```

Shows:
- Completion percentage
- Molecules processed / total
- Visual progress indicator

### Status Updates
```
Status: Processing...
ETA: 2 minutes remaining
Rate: ~100 molecules/second
```

Shows:
- Current status
- Estimated time remaining
- Processing rate

### Real-time Updates

Progress updates every 0.5 seconds via WebSocket connection.

**If WebSocket disconnects:**
- Progress bar may freeze
- Processing continues on server
- Refresh page to reconnect
- Results are not lost

## Handling Failures

Individual molecule failures don't stop the batch:

### Failure Types

**Parse Failures:**
```
Molecule 42: Invalid SMILES syntax
Error: "Unmatched parenthesis in SMILES"
```

**Validation Failures:**
```
Molecule 108: Valence error
Error: "Atom 3 (C) has invalid valence 5"
```

**Processing Errors:**
```
Molecule 256: Internal error
Error: "Descriptor calculation failed"
```

### Failure Handling

1. **Logged**: Error message captured
2. **Marked**: Status set to "failed"
3. **Continued**: Next molecule processed
4. **Reported**: Included in final results

### Reviewing Failures

In exported results:
- Failed molecules have status="failed"
- Error message in "error" column
- Original input preserved
- Easy to filter and retry

## Canceling a Batch

You can cancel at any time:

1. Click **Cancel** button
2. Processing stops
3. Completed results are preserved
4. Partial results available for export

**What happens:**
- Molecules already processed: Saved
- Currently processing: May complete
- Not yet started: Skipped

**Use case:** Realized input needs correction, want partial results.

## Exporting Results

After processing completes, export in multiple formats:

### CSV Export
```
id,smiles,quality_score,ml_readiness,status,issues
1,CCO,100,99.5,success,0
2,CC(O)N,90,98.2,success,1
3,invalid,0,0,failed,Parse error
```

**Best for:**
- Spreadsheet analysis
- Database import
- Filtering and sorting
- Scripting/automation

### Excel Export
```
Similar to CSV but:
- Formatted headers
- Conditional formatting (score colors)
- Multiple sheets (summary + details)
- Freeze panes
```

**Color coding:**
- Green: Score ≥ 80
- Yellow: Score 50-79
- Red: Score < 50

**Best for:**
- Reports
- Manual review
- Stakeholder presentations

### SDF Export
```
Original SDF with added fields:
> <quality_score>
95

> <ml_readiness>
98.5

> <validation_issues>
0
```

**Best for:**
- Loading into other chemistry tools
- Maintaining 3D coordinates
- Preserving metadata

### JSON Export
```json
[
  {
    "id": 1,
    "smiles": "CCO",
    "quality_score": 100,
    "ml_readiness": 99.5,
    "issues": [],
    "metadata": {...}
  }
]
```

**Best for:**
- Programmatic access
- Web services integration
- Custom analysis scripts

### PDF Export
```
Multi-page report with:
- Summary statistics
- Top issues
- Structure images
- Detailed results table
```

**Best for:**
- Documentation
- Archival
- Regulatory submissions
- Non-technical stakeholders

## Performance

### Expected Throughput

| Metric | Target | Typical |
|--------|--------|---------|
| Processing rate | >100 mol/s | ~127 mol/s |
| 1,000 molecules | <10 seconds | ~8 seconds |
| 10,000 molecules | <5 minutes | ~79 seconds |

### Performance Factors

**Fast processing:**
- Simple molecules (MW < 500)
- Valid structures
- No structural alerts
- Cached validation results

**Slow processing:**
- Large molecules (MW > 1000)
- Complex ring systems
- Many stereocenters
- First-time validation (no cache)

### Optimization Tips

1. **Remove duplicates** before uploading
2. **Standardize input** format (canonical SMILES)
3. **Use CSV** for SMILES-only data (faster than SDF)
4. **Split very large batches** (>20,000 molecules)
5. **Reuse cache** (same molecules validate faster)

## Architecture

Understanding the backend helps troubleshoot:

```
Upload → Redis Queue → Celery Workers (4) → PostgreSQL
                              ↓
                    Progress WebSocket
```

**Components:**
- **Redis**: Job queue
- **Celery**: Parallel processing with 4 workers
- **WebSocket**: Real-time progress updates
- **PostgreSQL**: Results storage

**Chunk processing:**
- Batch split into chunks of 100 molecules
- Each chunk processed in parallel
- Progress aggregated across chunks

## Common Workflows

### Workflow 1: Library Cleanup

```
1. Export compound library as CSV
2. Upload to ChemStructVal batch
3. Export results as Excel
4. Filter for quality_score >= 80
5. Identify and fix low-quality entries
6. Re-export cleaned library
```

### Workflow 2: Database Migration

```
1. Export old database as SDF
2. Batch validate in ChemStructVal
3. Review failures (parse errors, valence)
4. Standardize with ChEMBL pipeline
5. Export standardized SDF
6. Import to new database
```

### Workflow 3: ML Dataset Preparation

```
1. Upload candidate molecules
2. Batch validate and score
3. Export as CSV
4. Filter ml_readiness >= 95
5. Check for structural alerts
6. Remove PAINS if needed
7. Final dataset ready for training
```

### Workflow 4: Hit Triage

```
1. Upload screening hits (CSV from plate reader)
2. Batch validate
3. Filter quality_score >= 70
4. Flag PAINS alerts
5. Export PDF report for chemists
6. Prioritize for synthesis
```

## API Access

For programmatic batch processing, use the Python client:

```python
from chemstructval import ChemStructValClient

client = ChemStructValClient(
    base_url="http://localhost:8000",
    api_key="your-api-key"  # Optional, higher rate limits
)

# Submit batch job
job = client.submit_batch_file("compounds.csv", smiles_column="smiles")

# Monitor progress
for progress in client.monitor_batch(job.job_id):
    print(f"Progress: {progress.percent}%")

# Get results when complete
results = client.get_batch_results(job.job_id)

# Export
client.export_batch(job.job_id, format="csv", output="results.csv")
```

See [Python Client Documentation](../api/python-client.md) for details.

## Troubleshooting

### Upload fails
- Check file size limits (CSV: 50MB, SDF: 100MB)
- Verify file format (valid CSV/SDF)
- Check file encoding (use UTF-8)

### Progress stuck at 0%
- Check Celery worker is running: `docker-compose logs celery-worker`
- Verify Redis connection: `redis-cli ping`
- Check backend logs: `docker-compose logs backend`

### WebSocket disconnects
- Normal for long-running jobs (>10 min)
- Refresh page to reconnect
- Results are not lost
- Can still export after completion

### High failure rate
- Review first 10 failures for patterns
- Check input format consistency
- Verify SMILES are valid
- Consider standardizing input first

### Slow processing
- Check system resources (CPU, memory)
- Reduce batch size if OOM errors
- Check Prometheus metrics for bottlenecks
- Consider scaling Celery workers

See [FAQ](../troubleshooting/faq.md) for more troubleshooting.

## Best Practices

1. **Start small**: Test with 10-100 molecules first
2. **Check preview**: Verify column mapping for CSV
3. **Monitor progress**: Watch for high failure rates
4. **Export early**: Download results before closing browser
5. **Document settings**: Note any custom configurations
6. **Keep originals**: Always maintain original input files
7. **Version results**: Track validation runs with timestamps

## Limits and Quotas

### Anonymous Users
- 1 batch job at a time
- Maximum 1,000 molecules per batch
- Rate limit: 10 requests/minute

### Authenticated Users (API Key)
- 5 concurrent batch jobs
- Maximum 10,000 molecules per batch
- Rate limit: 300 requests/minute

### Enterprise
- Unlimited concurrent jobs
- No batch size limit
- Custom rate limits
- Contact for details

To create an API key:
```bash
curl -X POST http://localhost:8000/api/v1/api-keys \
  -H "Content-Type: application/json" \
  -d '{"name": "My API Key"}'
```

Returns:
```json
{
  "key": "csv_abc123...",
  "name": "My API Key",
  "created_at": "2026-01-24T10:00:00Z"
}
```

**Important**: Save the key immediately - it cannot be retrieved later.
