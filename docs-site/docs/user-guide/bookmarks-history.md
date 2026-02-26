---
sidebar_position: 8
title: Bookmarks & History
description: Save result snapshots and track your validation audit trail
---

# Bookmarks & History

ChemAudit lets you bookmark validation results for later reference and automatically tracks every validation in an audit trail. Both features are session-scoped for privacy and support GDPR-compliant data erasure.

## Bookmarks

Bookmarks save a snapshot of a molecule's validation results so you can return to them later without re-validating.

### Creating a Bookmark

After validating a molecule on the Single Validation page, click the **Bookmark** button (star icon) in the results header. A confirmation toast appears and the star fills to indicate the molecule is saved.

The bookmark automatically captures:

- SMILES string
- Molecule name (if available)
- InChIKey (computed from SMILES)
- Tags and notes (editable)
- Source context (e.g., `single_validation` or `batch`)
- Full result snapshot stored locally in the browser

### Bookmarks Page

Navigate to `/bookmarks` to manage your saved molecules.

**Search and Filter:**

- Search by SMILES substring using the search bar
- Filter by tag using the tag chips displayed above the table
- Filter by source (single validation, batch)

**Table Columns:**

| Column | Description |
|--------|-------------|
| **Molecule** | Name + truncated SMILES (click chevron to expand full SMILES + notes) |
| **InChIKey** | First 14 characters for quick identification |
| **Tags** | Up to 3 visible tags with "+N more" indicator |
| **Actions** | View (opens in validation) and Delete |

A "Results saved" badge appears on rows that have a locally cached result snapshot.

### Managing Bookmarks

**Edit tags and notes:** Click a bookmark to expand it, then modify tags or notes.

**Delete:** Click the trash icon on individual bookmarks, or select multiple and use **Delete Selected** in the bulk action bar.

**Bulk actions:** Select bookmarks using checkboxes, then use the floating action bar:

- **Submit as Batch** — creates a new batch job from selected bookmarks
- **Delete Selected** — removes all selected bookmarks
- **Clear** — deselects all

### Submit Bookmarks as Batch

Select one or more bookmarks and click **Submit as Batch**. ChemAudit creates a new batch job containing those molecules and redirects you to the batch results page.

### API Reference

```bash
# List bookmarks (paginated, with optional filters)
curl "http://localhost:8001/api/v1/bookmarks?page=1&page_size=50&tag=kinase&search=CCO"

# Create a bookmark
curl -X POST http://localhost:8001/api/v1/bookmarks \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "name": "Aspirin",
    "tags": ["analgesic", "nsaid"],
    "notes": "Good validation score",
    "source": "single_validation"
  }'

# Update bookmark tags/notes
curl -X PUT http://localhost:8001/api/v1/bookmarks/42 \
  -H "Content-Type: application/json" \
  -d '{
    "tags": ["analgesic", "nsaid", "approved"],
    "notes": "Updated notes"
  }'

# Delete a single bookmark
curl -X DELETE http://localhost:8001/api/v1/bookmarks/42

# Bulk delete bookmarks
curl -X DELETE "http://localhost:8001/api/v1/bookmarks/bulk?ids=1&ids=2&ids=3"

# Submit bookmarks as a batch job
curl -X POST http://localhost:8001/api/v1/bookmarks/batch-submit \
  -H "Content-Type: application/json" \
  -d '{"bookmark_ids": [1, 2, 3, 4, 5]}'
```

**Bookmark creation response:**

```json
{
  "id": 42,
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "name": "Aspirin",
  "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "tags": ["analgesic", "nsaid"],
  "notes": "Good validation score",
  "source": "single_validation",
  "job_id": null,
  "created_at": "2026-02-26T10:30:00Z"
}
```

**Batch submit response:**

```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "molecule_count": 5,
  "message": "Job submitted. Processing 5 molecules."
}
```

## Validation History

Every validation — single or batch — is recorded in an immutable audit trail. This provides a complete record of what was validated, when, and the outcome.

### What Gets Logged

Each audit entry captures:

| Field | Description |
|-------|-------------|
| **SMILES** | The validated molecule |
| **InChIKey** | Computed identifier |
| **Outcome** | `pass`, `warn`, or `fail` |
| **Score** | Validation score (0–100) |
| **Source** | `single` or `batch` |
| **Job ID** | Batch job ID (for batch validations) |
| **Timestamp** | When the validation occurred |

:::info Immutable Audit Trail
History entries are append-only. Individual entries cannot be edited or deleted — only a full session purge removes them.
:::

### History Page

Navigate to `/history` to browse your validation history.

**Statistics Grid:**

The top of the page shows summary cards:

- **Total Validations** — total count across all sessions
- **Passed** — validations with `pass` outcome (green)
- **Warnings** — validations with `warn` outcome (amber)
- **Failed** — validations with `fail` outcome (red)

**Filtering:**

| Filter | Description |
|--------|-------------|
| **Date Range** | From/To date pickers |
| **Outcome** | `pass`, `warn`, or `fail` |
| **Source** | `single` or `batch` |
| **SMILES Search** | Substring match on SMILES |

Click **Apply Filters** to update results or **Reset** to clear all filters.

**Results Table:**

Results are paginated (25 per page, newest first) with columns:

| Column | Description |
|--------|-------------|
| **Date** | Validation timestamp |
| **SMILES** | Truncated SMILES string |
| **Outcome** | Color-coded badge (green/amber/red) |
| **Score** | Color-coded score (green ≥70, amber ≥40, red below 40) |
| **Source** | `single` or `batch` badge |
| **Job ID** | Truncated batch job ID (if applicable) |

### API Reference

```bash
# Get paginated history with filters
curl "http://localhost:8001/api/v1/history?page=1&page_size=50&outcome=pass&source=single&date_from=2026-02-01T00:00:00Z&date_to=2026-02-28T23:59:59Z"

# Get summary statistics
curl http://localhost:8001/api/v1/history/stats
```

**History response:**

```json
{
  "entries": [
    {
      "id": 1,
      "smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
      "outcome": "pass",
      "score": 95,
      "job_id": null,
      "source": "single",
      "created_at": "2026-02-26T10:30:00Z"
    }
  ],
  "total": 150,
  "page": 1,
  "page_size": 50
}
```

**Statistics response:**

```json
{
  "total_validations": 150,
  "outcome_distribution": {"pass": 120, "warn": 20, "fail": 10},
  "source_distribution": {"single": 80, "batch": 70}
}
```

## Session Scope

Bookmarks and history are isolated per user:

- **Anonymous users**: Data is scoped to a session cookie (`chemaudit_sid`), which lasts 30 days
- **API key users**: Data is scoped to the API key hash (SHA-256)
- **Database-level isolation**: PostgreSQL Row-Level Security ensures users can only access their own data

API key-scoped data takes precedence when both a session cookie and API key are present.

## Data Retention & Privacy

### Auto-Purge

Anonymous session data is automatically purged after 30 days of inactivity via a scheduled Celery Beat task. API key-scoped data is not auto-purged.

### GDPR Erasure

The **Privacy** page includes a "Purge My Data" button that immediately deletes all your data:

```bash
# Delete all bookmarks and history for current session
curl -X DELETE http://localhost:8001/api/v1/me/data
```

**Response:**

```json
{
  "status": "purged",
  "deleted": {
    "bookmarks": 42,
    "history": 150
  }
}
```

This removes:

- All bookmarks (server-side database records)
- All validation history entries
- Locally cached result snapshots (IndexedDB cleared in the browser)

:::warning Irreversible
Data purge is permanent and cannot be undone. Make sure to export any data you need before purging.
:::

## Next Steps

- **[Single Validation](/docs/user-guide/single-validation)** — Validate molecules and bookmark results
- **[Batch Processing](/docs/user-guide/batch-processing)** — Process datasets and track history
- **[Exporting Results](/docs/user-guide/exporting-results)** — Export before purging data
- **[Scoring Profiles](/docs/user-guide/scoring/profiles)** — Apply custom scoring to bookmarked molecules
