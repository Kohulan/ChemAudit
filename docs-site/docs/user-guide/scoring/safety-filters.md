---
sidebar_position: 4
title: Safety Filters
description: Screen for problematic substructures using PAINS, BRENK, NIH, ZINC, and ChEMBL filters
---

# Safety Filter Scoring

Safety filter scoring evaluates molecules against structural alert databases to identify potentially problematic compounds. This is a scored summary of the detailed [Structural Alerts](/docs/user-guide/structural-alerts) screening.

## What It Measures

Safety filter scoring tests molecules against multiple catalogs:

| Catalog | Patterns | Focus |
|---------|----------|-------|
| **PAINS A/B/C** | ~480 | Pan-assay interference compounds — frequent hitters in HTS |
| **Brenk** | ~105 | Unfavorable chemical moieties for drug development |
| **NIH** | ~180 | NIH MLSMR screening exclusion filters |
| **ZINC** | ~95 | Drug-likeness and reactivity filters |
| **ChEMBL** | ~700+ | Combined from 7 pharma sub-catalogs (see below) |

### ChEMBL Sub-Catalogs

The ChEMBL structural alerts combine filters from 7 pharmaceutical industry sources:

| Sub-Catalog | Source | Focus |
|-------------|--------|-------|
| **BMS** | Bristol-Myers Squibb | Reactive functional groups and chemical liabilities |
| **Dundee** | University of Dundee | Promiscuous compound filters for screening |
| **Glaxo** | GlaxoSmithKline | Undesirable moieties and toxicophores |
| **Inpharmatica** | Inpharmatica Ltd. | Chemical liabilities and ADMET flags |
| **LINT** | Lead Identification Noise | Noise-causing patterns in lead identification |
| **MLSMR** | NIH/MLSMR | Molecular Libraries Screening Center Network filters |
| **SureChEMBL** | EMBL-EBI | Patent literature structural alerts |

Each sub-catalog targets different aspects of compound quality. When ChEMBL alerts are enabled, all 7 sub-catalogs are screened simultaneously.

All pattern matching uses RDKit's `FilterCatalog` module with SMARTS substructure matching.

## Scoring Output

The safety filter score provides a pass/fail summary:

```json
{
  "safety_filters": {
    "pains": {
      "passed": true,
      "alerts": [],
      "alert_count": 0
    },
    "brenk": {
      "passed": true,
      "alerts": [],
      "alert_count": 0
    },
    "nih": {
      "passed": true,
      "alerts": [],
      "alert_count": 0
    },
    "zinc": {
      "passed": true,
      "alerts": [],
      "alert_count": 0
    },
    "chembl": {
      "passed": true,
      "total_alerts": 0
    },
    "all_passed": true,
    "total_alerts": 0,
    "interpretation": "No safety alerts detected"
  }
}
```

## Interpretation

| Result | Interpretation | Recommendation |
|--------|---------------|----------------|
| **All passed** | No alerts found | Good safety profile for screening |
| **PAINS only** | Assay interference risk | Review assay compatibility before screening |
| **Multiple catalogs** | Multiple issues | Investigate further, consider excluding |

## Safety Filter vs. Structural Alerts

**Safety Filter Scoring (this page):**
- Part of comprehensive scoring API
- Pass/fail summary only
- Quick screening
- Included in batch processing

**[Structural Alerts](/docs/user-guide/structural-alerts):**
- Detailed alert information
- Matched atoms identified
- Pattern names and descriptions
- Dedicated API endpoint

Use safety filter scoring for quick filtering, then use detailed structural alerts for investigation.

## Use in Batch Processing

In batch mode, safety filter pass rate is reported:

```json
{
  "statistics": {
    "safety_pass_rate": 0.91
  }
}
```

This indicates 91% of molecules pass all safety filters.

## Best Practices

1. **Use for prioritization**: Not absolute rejection criteria
2. **Review alerts**: Always investigate specific alerts for hits
3. **Consider context**: 87 FDA-approved drugs match PAINS patterns
4. **Check assay type**: PAINS alerts are assay-specific
5. **Document decisions**: Record why you accepted/rejected flagged molecules

## Next Steps

- **[Structural Alerts](/docs/user-guide/structural-alerts)** - Detailed alert information
- **[Scoring Overview](/docs/user-guide/scoring/overview)** - All scoring systems
- **[API Reference](/docs/api/endpoints)** - Full API documentation
