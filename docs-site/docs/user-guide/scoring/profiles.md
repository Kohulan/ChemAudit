---
sidebar_position: 9
title: Scoring Profiles
description: Customizable property thresholds and weights for compound filtering
---

# Scoring Profiles

Scoring profiles let you define custom property thresholds and weights for evaluating molecules against your specific criteria. ChemAudit includes 8 built-in presets and a visual profile builder for creating your own.

## Built-in Presets

Eight immutable presets cover common drug discovery use cases. Presets cannot be edited or deleted, but can be duplicated as a starting point for custom profiles.

| Preset | Key Thresholds | Description |
|--------|----------------|-------------|
| **Drug-like (Lipinski)** | MW 0–500, LogP −5 to 5, HBD 0–5, HBA 0–10 | Lipinski Rule of Five |
| **Lead-like** | MW 200–450, LogP −1 to 4, RotB 0–7 | Lead-like chemical space |
| **Fragment-like (Ro3)** | MW 0–300, LogP −3 to 3, HBD 0–3, HBA 0–3, RotB 0–3 | Fragment screening |
| **CNS-penetrant** | MW 0–400, LogP 1–5, TPSA 0–90, HBD 0–2 | Blood-brain barrier permeable |
| **Ghose (Amgen)** | MW 160–480, LogP −0.4 to 5.6, Atoms 20–70, Refractivity 40–130 | Ghose filter |
| **Veber (GSK)** | RotB 0–10, TPSA 0–140 | Oral bioavailability |
| **PPI-like** | MW 400–800, LogP 2–6, RotB 5–15 | Protein-protein interaction space |
| **NP-like** | MW 200–800, RotB 0–10, HBD 0–5, Rings 2–20 | Natural product space |

### Using Presets

Navigate to the **Profiles** page (`/profiles`) and click the **Presets** tab. Each preset shows its name, description, and up to 4 key thresholds. Click **Apply** to set it as your active profile, or **Customize** to duplicate it as a custom profile.

## Custom Profile Builder

The **Custom Builder** tab provides a visual editor for creating profiles with precise control over thresholds and weights.

### Threshold Properties

Eight molecular properties can be configured with minimum and maximum thresholds:

| Property | Slider Range | Step | Unit |
|----------|-------------|------|------|
| **Molecular Weight** | 0–1000 | 10 | Da |
| **LogP** | −5 to 10 | 0.5 | — |
| **TPSA** | 0–300 | 5 | A² |
| **H-Bond Donors** | 0–10 | 1 | — |
| **H-Bond Acceptors** | 0–20 | 1 | — |
| **Rotatable Bonds** | 0–20 | 1 | — |
| **Aromatic Rings** | 0–10 | 1 | — |
| **Fsp3** | 0–1 | 0.05 | — |

Each property has a dual-range slider for setting min/max bounds and a weight input (0–2, step 0.1).

### Weights

Weights control how much each property contributes to the composite score. The default weight is 1.0. Higher weights increase the property's influence on the final score. Set a weight to 0 to exclude a property from scoring entirely.

### Saving and Sharing

- **Save** — creates or updates a custom profile
- **Reset** — restores default thresholds (Lipinski-like)
- **Export** — downloads the profile as JSON for sharing
- **Import** — loads a profile from a JSON file

## How Profile Scoring Works

### Desirability Function

Each property is evaluated using a desirability function that returns a value between 0 and 1:

- **In range** (between min and max): desirability = 1.0
- **Outside range**: linear falloff from 1.0 to 0.0 over a distance equal to the range width

For example, with MW thresholds of 200–500 (range width = 300):

| MW Value | Desirability | Explanation |
|----------|-------------|-------------|
| 350 | 1.0 | Within range |
| 200 | 1.0 | At boundary |
| 650 | 0.5 | 50% beyond upper bound |
| 800 | 0.0 | 100% beyond upper bound |

### Composite Score

The composite profile score uses a QED-style **weighted geometric mean** of per-property desirabilities, scaled to 0–100:

```
Score = (∏ desirability_i ^ weight_i) ^ (1 / Σ weight_i) × 100
```

This means:

- A single property with desirability 0 drives the entire score to 0
- All properties must be reasonably in range for a high score
- Properties with higher weights have more influence

### Per-Property Breakdown

The scoring result includes a per-property breakdown showing:

```json
{
  "profile_id": 1,
  "profile_name": "Drug-like (Lipinski)",
  "score": 85.3,
  "properties": {
    "mw": {
      "value": 342.1,
      "min": 0,
      "max": 500,
      "in_range": true,
      "desirability": 1.0
    },
    "logp": {
      "value": 2.3,
      "min": -5,
      "max": 5,
      "in_range": true,
      "desirability": 1.0
    }
  }
}
```

## Using Profiles in Batch Processing

### Apply During Upload

In the batch upload flow, expand the **Scoring Profile** sidebar to select a profile before uploading. The selected profile is applied to every molecule during processing.

```bash
# Upload batch with a scoring profile
curl -X POST http://localhost:8001/api/v1/batch/upload \
  -F "file=@molecules.sdf" \
  -F "profile_id=1"
```

### Profile Score in Results

When a profile is applied, batch results include a **Profile Score** column. The batch statistics also report the average profile score across all molecules.

### Re-Score with Different Profile

After processing, you can re-score a subset of molecules with a different profile using the **Subset Actions** panel — no need to re-upload the file. See [Subset Actions](/docs/user-guide/subset-actions) for details.

## API Reference

```bash
# List all profiles (presets + custom)
curl http://localhost:8001/api/v1/profiles

# Get a single profile
curl http://localhost:8001/api/v1/profiles/1

# Create a custom profile
curl -X POST http://localhost:8001/api/v1/profiles \
  -H "Content-Type: application/json" \
  -d '{
    "name": "My Custom Profile",
    "description": "Optimized for CNS targets",
    "thresholds": {
      "mw": {"min": 150, "max": 400},
      "logp": {"min": 1, "max": 4},
      "tpsa": {"min": 20, "max": 90},
      "hbd": {"min": 0, "max": 2}
    },
    "weights": {
      "mw": 1.0,
      "logp": 1.5,
      "tpsa": 1.2,
      "hbd": 1.0
    }
  }'

# Update a custom profile (presets return 400)
curl -X PUT http://localhost:8001/api/v1/profiles/9 \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Updated Profile",
    "thresholds": {"mw": {"min": 100, "max": 450}}
  }'

# Delete a custom profile (presets return 400)
curl -X DELETE http://localhost:8001/api/v1/profiles/9

# Duplicate any profile (including presets)
curl -X POST http://localhost:8001/api/v1/profiles/1/duplicate \
  -H "Content-Type: application/json" \
  -d '{"name": "My Lipinski Variant"}'

# Export profile as JSON
curl http://localhost:8001/api/v1/profiles/1/export -o profile.json

# Import profile from JSON
curl -X POST http://localhost:8001/api/v1/profiles/import \
  -H "Content-Type: application/json" \
  -d @profile.json
```

**Profile response:**

```json
{
  "id": 1,
  "name": "Drug-like (Lipinski)",
  "description": "Lipinski Rule of Five: MW<500, LogP<5, HBD<=5, HBA<=10",
  "thresholds": {
    "mw": {"min": 0, "max": 500},
    "logp": {"min": -5, "max": 5},
    "hbd": {"min": 0, "max": 5},
    "hba": {"min": 0, "max": 10}
  },
  "weights": {
    "mw": 1.0,
    "logp": 1.0,
    "hbd": 1.0,
    "hba": 1.0
  },
  "is_preset": true,
  "is_active": true,
  "created_at": "2026-02-01T00:00:00Z",
  "updated_at": "2026-02-01T00:00:00Z"
}
```

:::info Preset Protection
Built-in presets cannot be modified or deleted. Attempting to update or delete a preset returns HTTP 400. Use the duplicate endpoint to create an editable copy.
:::

## Next Steps

- **[Batch Processing](/docs/user-guide/batch-processing)** — Apply profiles to batch jobs
- **[Subset Actions](/docs/user-guide/subset-actions)** — Re-score subsets with different profiles
- **[Exporting Results](/docs/user-guide/exporting-results)** — Export profile scores with results
- **[Scoring Overview](/docs/user-guide/scoring/overview)** — All scoring systems
