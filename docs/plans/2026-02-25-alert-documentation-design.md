# Design: ChEMBL Structural Alert Documentation & Citations

**Date**: 2026-02-25
**Branch**: features_v1

## Problem

When a structural alert fires (e.g., morphine → "Dundee: isolated alkene"), the UI shows:
- A cryptic catalog code ("CHEMBL_DUNDEE")
- The raw RDKit pattern name as the "description" (just echoes the name)
- No explanation of what the filter set is or why the pattern matters
- No citation to the original research
- No functional group category to help chemists assess concern level

A medicinal chemist seeing "Dundee — 1 alert" has no actionable context.

## Solution

### 1. Extract RDKit Built-in Metadata

RDKit `FilterCatalogEntry` objects already store per-entry properties that the current code ignores:
- `Reference` — paper citation or URL (e.g., "Brenk R et al. ChemMedChem 3 (2008) 435-444")
- `Scope` — what the filter set screens for (e.g., "unwanted functionality due to potential tox reasons")
- `FilterSet` — which rule set it belongs to (e.g., "ChEMBL23_Dundee")

Extract these in `alert_manager.py` and pass them through to the API response.

### 2. Enrich Catalog Metadata with Proper Citations

Update `AVAILABLE_CATALOGS` in `filter_catalog.py` with:
- Full human-readable names
- What each filter set screens for
- Original paper citations with DOIs/PMIDs
- Pattern counts (from actual RDKit catalog, not hardcoded)

| Catalog | Full Name | Citation |
|---------|-----------|----------|
| PAINS | Pan-Assay Interference Compounds | Baell & Holloway, J Med Chem 53 (2010) 2719-2740 |
| BRENK | Brenk Structural Alerts | Brenk et al., ChemMedChem 3 (2008) 435-444 |
| NIH | NIH/NCGC Alerts | Jadhav et al., J Med Chem 53 (2009) 37-51 |
| ZINC | ZINC Drug-likeness Filters | Irwin & Shoichet, J Chem Inf Model 45 (2005) 177-182 |
| CHEMBL_BMS | Bristol-Myers Squibb HTS Deck Filters | Pearce et al., J Chem Inf Model 46 (2006) 1060-1068 |
| CHEMBL_DUNDEE | U. of Dundee NTD Screening Filters | Brenk et al., ChemMedChem 3 (2008) 435-444 |
| CHEMBL_GLAXO | Glaxo Wellcome Hard Filters | Hann et al., J Chem Inf Comput Sci 39 (1999) 897-902 |
| CHEMBL_INPHARMATICA | Inpharmatica Unwanted Fragments | Inpharmatica Ltd (internal, no peer-reviewed pub) |
| CHEMBL_LINT | Lilly MedChem Rules | Bruns & Watson, J Med Chem 55 (2012) 9763-9772 |
| CHEMBL_MLSMR | NIH MLSMR Excluded Functionality | NIH Molecular Libraries Program |
| CHEMBL_SURECHEMBL | SureChEMBL Patent Alerts | Papadatos et al., Nucleic Acids Res 44 (2016) D1220-D1228 |

### 3. Pattern Category Classification

Map pattern names to broad concern categories using keyword rules on the pattern name string:

**Categories:**
- **Reactive Group** — electrophilic/nucleophilic reactivity (Michael acceptors, acyl halides, epoxides, etc.)
- **Metabolic Liability** — metabolically labile groups (isolated alkenes, aldehydes, esters, etc.)
- **Toxicophore** — known toxic substructures (nitro, aniline, hydrazine, heavy metals, etc.)
- **Assay Interference** — HTS assay artifacts (PAINS patterns, fluorescent compounds, aggregators)
- **Physicochemical** — property-based filters (molecular weight, atom counts, rotatable bonds)
- **Unwanted Functionality** — general undesirable groups that don't fit above categories

Implementation: ~40 keyword-to-category rules covering the pattern namespace. PAINS patterns default to "Assay Interference". Unmatched patterns default to "Unwanted Functionality".

### 4. Schema Changes

**Backend `AlertResultSchema`** — add fields:
- `reference: Optional[str]` — paper citation from RDKit entry
- `scope: Optional[str]` — what the filter screens for from RDKit entry
- `filter_set: Optional[str]` — filter set identifier
- `catalog_description: Optional[str]` — human-readable catalog name + description
- `category: Optional[str]` — functional group concern category

**Backend `CatalogInfoSchema`** — add fields:
- `reference: Optional[str]` — paper citation
- `scope: Optional[str]` — what it screens for
- `doi: Optional[str]` — DOI link
- `pmid: Optional[str]` — PubMed ID

**Backend `AlertScreenRequest`** — extend catalog validator to accept ChEMBL catalogs.

### 5. Frontend Changes

**TypeScript types** — mirror new schema fields.

**AlertCard** — redesign to show:
- Pattern name (formatted)
- Human-readable catalog name (not "CHEMBL_DUNDEE")
- Category badge (color-coded by concern type)
- Filter set description (scope)
- Matched atoms with hover
- Reference citation (collapsible or subtle)

**AlertResults** — catalog group headers show full context.

**SingleValidation** — ChEMBL catalog selector: collapsible group with master toggle + individual checkboxes.

### 6. Batch Integration

Batch alert results in `tasks.py` will also include the enriched fields so batch export/display benefits too.

## Non-Goals

- Individual per-pattern prose descriptions for all ~1265 patterns
- Changing the underlying RDKit screening logic
- Adding new alert catalogs beyond what RDKit provides
