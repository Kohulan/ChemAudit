# Alert Documentation & Citations — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Enrich structural alerts with catalog-level context, citations, functional group categories, and RDKit metadata so chemists understand what each alert means and why it was flagged.

**Architecture:** Extract Reference/Scope/FilterSet from RDKit FilterCatalogEntry (already embedded, currently ignored). Build keyword-based pattern category classifier. Propagate enriched data through schemas to frontend AlertCard. Add collapsible ChEMBL catalog selector.

**Tech Stack:** Python/FastAPI backend, RDKit FilterCatalog API, Pydantic v2 schemas, React/TypeScript frontend, Tailwind CSS

---

### Task 1: Enrich AVAILABLE_CATALOGS with citations and metadata

**Files:**
- Modify: `backend/app/services/alerts/filter_catalog.py`

**Step 1: Update AVAILABLE_CATALOGS dictionary**

Replace the current `AVAILABLE_CATALOGS` dict with enriched metadata. Each entry gets: `name`, `description`, `scope`, `reference`, `doi`, `pmid`, `pattern_count`, `severity`.

```python
AVAILABLE_CATALOGS: Dict[str, Dict[str, str]] = {
    "PAINS": {
        "name": "PAINS (Pan-Assay Interference Compounds)",
        "description": "Substructure filters for removal of frequent hitters in high-throughput screens",
        "scope": "Identifies compounds that appear active in multiple assay types due to non-specific mechanisms such as aggregation, redox cycling, fluorescence interference, or chemical reactivity rather than genuine target binding",
        "reference": "Baell JB, Holloway GA. New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. J Med Chem 53 (2010) 2719-2740.",
        "doi": "10.1021/jm901137j",
        "pmid": "20131845",
        "pattern_count": "480",
        "severity": "warning",
        "note": "87 FDA-approved drugs contain PAINS patterns — these are warnings, not rejections",
    },
    "PAINS_A": {
        "name": "PAINS Class A (Most Frequent Hitters)",
        "description": "16 most promiscuous PAINS patterns — highest hit rates across multiple assay technologies",
        "scope": "PAINS filters (family A)",
        "reference": "Baell JB, Holloway GA. J Med Chem 53 (2010) 2719-2740.",
        "doi": "10.1021/jm901137j",
        "pmid": "20131845",
        "pattern_count": "16",
        "severity": "warning",
    },
    "PAINS_B": {
        "name": "PAINS Class B (Moderate Hitters)",
        "description": "55 moderately promiscuous PAINS patterns",
        "scope": "PAINS filters (family B)",
        "reference": "Baell JB, Holloway GA. J Med Chem 53 (2010) 2719-2740.",
        "doi": "10.1021/jm901137j",
        "pmid": "20131845",
        "pattern_count": "55",
        "severity": "warning",
    },
    "PAINS_C": {
        "name": "PAINS Class C (Lower Frequency Hitters)",
        "description": "409 less frequent but still notable PAINS patterns",
        "scope": "PAINS filters (family C)",
        "reference": "Baell JB, Holloway GA. J Med Chem 53 (2010) 2719-2740.",
        "doi": "10.1021/jm901137j",
        "pmid": "20131845",
        "pattern_count": "409",
        "severity": "warning",
    },
    "BRENK": {
        "name": "Brenk Structural Alerts",
        "description": "Filters for potentially problematic functional groups identified during assembly of screening libraries for neglected tropical disease drug discovery",
        "scope": "Unwanted functionality due to potential toxicity reasons or unfavourable pharmacokinetic properties",
        "reference": "Brenk R, Schipani A, James D, et al. Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases. ChemMedChem 3 (2008) 435-444.",
        "doi": "10.1002/cmdc.200700139",
        "pmid": "18064617",
        "pattern_count": "105",
        "severity": "warning",
    },
    "NIH": {
        "name": "NIH/NCGC Alerts",
        "description": "Structural alert patterns from the NIH Chemical Genomics Center to annotate compounds with problematic functional groups",
        "scope": "Annotate compounds with problematic functional groups causing aggregation, autofluorescence, and reactivity artifacts",
        "reference": "Jadhav A, Ferreira RS, Klumpp C, et al. Quantitative Analyses of Aggregation, Autofluorescence, and Reactivity Artifacts in a Screen for Inhibitors of a Thiol Protease. J Med Chem 53 (2010) 37-51.",
        "doi": "10.1021/jm901070c",
        "pmid": "19908840",
        "pattern_count": "180",
        "severity": "warning",
    },
    "ZINC": {
        "name": "ZINC Drug-likeness Filters",
        "description": "Property and functional group filters from the ZINC database used to ensure drug-likeness of screening compounds",
        "scope": "Drug-likeness and unwanted functional group filters",
        "reference": "Irwin JJ, Shoichet BK. ZINC — A Free Database of Commercially Available Compounds for Virtual Screening. J Chem Inf Model 45 (2005) 177-182.",
        "doi": "10.1021/ci049714+",
        "pmid": "15667143",
        "pattern_count": "50",
        "severity": "warning",
    },
    "CHEMBL_BMS": {
        "name": "Bristol-Myers Squibb HTS Deck Filters",
        "description": "Empirically derived filters to remove promiscuous compounds from BMS high-throughput screening decks",
        "scope": "Identifies reactive and promiscuous functional groups that cause high false-positive rates in HTS",
        "reference": "Pearce BC, Sofia MJ, Good AC, et al. An Empirical Process for the Design of High-Throughput Screening Deck Filters. J Chem Inf Model 46 (2006) 1060-1068.",
        "doi": "10.1021/ci050504m",
        "pmid": "16711725",
        "pattern_count": "180",
        "severity": "warning",
    },
    "CHEMBL_DUNDEE": {
        "name": "University of Dundee NTD Screening Filters",
        "description": "Compound filters developed at the University of Dundee for neglected tropical disease screening library curation",
        "scope": "Screens for unwanted functionality due to potential toxicity or unfavourable pharmacokinetic properties",
        "reference": "Brenk R, Schipani A, James D, et al. Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases. ChemMedChem 3 (2008) 435-444.",
        "doi": "10.1002/cmdc.200700139",
        "pmid": "18064617",
        "pattern_count": "105",
        "severity": "warning",
    },
    "CHEMBL_GLAXO": {
        "name": "Glaxo Wellcome Hard Filters",
        "description": "Compound quality filters developed at GlaxoSmithKline to remove reactive and undesirable compounds before HTS pooling",
        "scope": "Identifies reactive, unstable, or otherwise undesirable compounds to exclude from screening pools",
        "reference": "Hann M, Hudson B, Lewell X, et al. Strategic Pooling of Compounds for High-Throughput Screening. J Chem Inf Comput Sci 39 (1999) 897-902.",
        "doi": "10.1021/ci990423o",
        "pmid": "10529988",
        "pattern_count": "55",
        "severity": "warning",
    },
    "CHEMBL_INPHARMATICA": {
        "name": "Inpharmatica Unwanted Fragment Filters",
        "description": "Filters developed at Inpharmatica Ltd to identify unwanted fragments including mutagenic groups, groups with unfavourable pharmacokinetic properties, and reactive functionalities",
        "scope": "Identifies unwanted fragments: mutagenic groups (e.g. nitro), poor PK groups (e.g. sulfates), and reactive groups (e.g. thiols, 2-halopyridines)",
        "reference": "Inpharmatica Ltd. (no peer-reviewed publication; included in ChEMBL structural alerts since ChEMBL 23)",
        "pattern_count": "91",
        "severity": "warning",
    },
    "CHEMBL_LINT": {
        "name": "Lilly MedChem Rules (LINT)",
        "description": "275 rules developed over 18 years at Eli Lilly to identify reactive or promiscuous compounds that interfere with biological assays",
        "scope": "Identifies reactive compounds, assay interference (fluorescence/absorbance), protein-damaging agents, unstable compounds, and non-drug-like molecules",
        "reference": "Bruns RF, Watson IA. Rules for Identifying Potentially Reactive or Promiscuous Compounds. J Med Chem 55 (2012) 9763-9772.",
        "doi": "10.1021/jm301008n",
        "pmid": "23061697",
        "pattern_count": "57",
        "severity": "warning",
    },
    "CHEMBL_MLSMR": {
        "name": "NIH MLSMR Excluded Functionality Filters",
        "description": "Structural filters from the NIH Molecular Libraries Small Molecule Repository to exclude compounds with problematic functionalities",
        "scope": "Excludes compounds with reactive, toxic, or otherwise undesirable functional groups from the NIH screening collection",
        "reference": "NIH Molecular Libraries Program. MLSMR Excluded Functionality Filters (included in ChEMBL structural alerts since ChEMBL 23).",
        "pattern_count": "116",
        "severity": "warning",
    },
    "CHEMBL_SURECHEMBL": {
        "name": "SureChEMBL Patent-derived Alerts",
        "description": "Structural alerts used in the SureChEMBL patent chemistry database for compound quality assessment",
        "scope": "Identifies known problematic substructures in patent-extracted chemical entities",
        "reference": "Papadatos G, Davies M, Dedman N, et al. SureChEMBL: a large-scale, chemically annotated patent document database. Nucleic Acids Res 44 (2016) D1220-D1228.",
        "doi": "10.1093/nar/gkv1253",
        "pmid": "26582922",
        "pattern_count": "166",
        "severity": "warning",
    },
    "ALL": {
        "name": "All Catalogs Combined",
        "description": "Screen against all available structural alert patterns",
        "scope": "Combined screening across all filter sets",
        "pattern_count": "variable",
        "severity": "warning",
    },
}
```

**Step 2: Verify catalog loads correctly**

Run: `cd backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -c "from app.services.alerts.filter_catalog import AVAILABLE_CATALOGS; print(len(AVAILABLE_CATALOGS)); print(AVAILABLE_CATALOGS['CHEMBL_DUNDEE']['reference'])"`

Expected: `15` and the Brenk citation string.

**Step 3: Commit**

```bash
git add backend/app/services/alerts/filter_catalog.py
git commit -m "feat: enrich AVAILABLE_CATALOGS with citations, DOIs, and scope descriptions"
```

---

### Task 2: Build pattern category classifier

**Files:**
- Create: `backend/app/services/alerts/pattern_categories.py`

**Step 1: Create the category classifier module**

This module maps pattern names to broad concern categories using keyword rules. It does NOT try to describe individual patterns — it classifies them by concern type.

```python
"""
Pattern Category Classifier for Structural Alerts

Maps structural alert pattern names to broad concern categories
using keyword-based rules. Categories help chemists quickly assess
what type of concern an alert represents.

Categories:
- Reactive Group: Electrophilic/nucleophilic reactivity
- Metabolic Liability: Metabolically labile groups
- Toxicophore: Known toxic substructures
- Assay Interference: HTS assay artifacts (PAINS default)
- Physicochemical: Property-based filters
- Unwanted Functionality: General undesirable groups
"""

from typing import Optional


# Ordered list of (keywords, category) — first match wins.
# Keywords are matched case-insensitively against the pattern name.
_CATEGORY_RULES: list[tuple[list[str], str]] = [
    # --- Reactive Groups ---
    (
        [
            "michael acceptor", "acyl halide", "acyl_halide", "acid halide", "acid_halide",
            "acyl chloride", "acyl fluoride", "acyl cyanide", "acyl_cyanide",
            "epoxide", "aziridine", "isocyanate", "isothiocyanate",
            "ketene", "acyl_phosph", "sulfonyl_halide", "sulfonyl halide",
            "allyl_halide", "alkyl halide", "alkyl_halide", "vinyl_halide",
            "activated_acetylene", "activated_vinyl", "activated_4mem",
            "activated_diazo", "activated_S#O",
            "beta_lactam", "beta-keto", "anhydride", "peroxide",
            "acetal", "halo_ether", "halo ether", "N-halo", "N_halo",
            "triflate", "pentafluorophenyl ester",
            "phosphorane", "phosphor", "silicon halogen",
            "reactive", "Reactive",
        ],
        "Reactive Group",
    ),
    # --- Toxicophores ---
    (
        [
            "nitro", "nitroso", "N-nitroso", "N_nitroso",
            "azide", "azido", "diazo",
            "hydrazine", "hydrazide", "acyl hydrazine", "acyl_hydrazine",
            "nitrogen_mustard", "nitrogen mustard",
            "aniline", "diaminobenzene",
            "heavy metal",
            "benzidine", "polycyclic aromatic",
            "cyanamide", "cyanohydrin",
            "hydroxamic acid", "hydroxamic_acid",
        ],
        "Toxicophore",
    ),
    # --- Metabolic Liability ---
    (
        [
            "isolated alkene", "isolated_alkene", "acyclic C=C",
            "acyclic_C=C", "stilbene", "polyene",
            "aldehyde", "thiol", "disulphide", "disulfide",
            "enamine", "imine", "oxime",
            "thioester", "thiocarbonyl", "thio_ester",
            "ester_of_HOBT", "ester of HOBT", "phenol ester", "phenol_ester",
            "phenyl carbonate", "phenyl_carbonate",
            "four member lactone", "four_member_lactone",
            "catechol", "hydroquinone", "quinone", "chinone",
            "N oxide", "N_oxide",
            "crown ether", "crown_ether",
            "Oxygen-nitrogen single bond", "oxygen_nitrogen",
            "sulfinic acid", "sulfonic acid", "Sulfonic_acid",
            "sulphate", "sulfur oxygen", "sulfur_oxygen",
            "Sulphur-nitrogen", "sulphur_nitrogen",
            "thiol", "triple bond", "triple_bond",
        ],
        "Metabolic Liability",
    ),
    # --- Physicochemical ---
    (
        [
            "Non-Hydrogen_atoms", "Non-Hydrogen atoms",
            "carbons", "N,O,S", "halogens",
            "rings", "rotatable", "chiral",
            "molecular_weight", "Molecular weight",
            "Aliphatic long chain", "aliphatic_long_chain",
            "perfluorinated", "Polycyclic aromatic hydrocarbon",
        ],
        "Physicochemical",
    ),
    # --- Assay Interference (PAINS-like) ---
    (
        [
            "azo", "Azo", "rhodanine",
            "thiazole_amine", "thiaz_ene",
            "pyrrole_A", "catechol_A",
            "ene_six_het", "hzone", "anil_di_alk",
            "mannich", "het_5_pyrazole",
            "quinone_A", "imine_one", "keto_phenone",
            "azo_A", "sulfonamide_A",
            "indol_3yl_alk", "phenol_sulfite",
            "biotin",
        ],
        "Assay Interference",
    ),
]

# Catalog-level defaults: if catalog is PAINS*, default to Assay Interference
_CATALOG_DEFAULTS: dict[str, str] = {
    "PAINS": "Assay Interference",
    "PAINS_A": "Assay Interference",
    "PAINS_B": "Assay Interference",
    "PAINS_C": "Assay Interference",
}


def classify_pattern(pattern_name: str, catalog_source: str) -> str:
    """
    Classify a structural alert pattern into a concern category.

    Args:
        pattern_name: The pattern name from RDKit FilterCatalog
        catalog_source: The catalog source (PAINS, BRENK, CHEMBL_DUNDEE, etc.)

    Returns:
        Category string: one of "Reactive Group", "Metabolic Liability",
        "Toxicophore", "Assay Interference", "Physicochemical",
        or "Unwanted Functionality" (default)
    """
    name_lower = pattern_name.lower()

    for keywords, category in _CATEGORY_RULES:
        for kw in keywords:
            if kw.lower() in name_lower:
                return category

    # Catalog-level default
    catalog_upper = catalog_source.upper()
    for prefix, default_cat in _CATALOG_DEFAULTS.items():
        if catalog_upper.startswith(prefix):
            return default_cat

    return "Unwanted Functionality"
```

**Step 2: Verify classification works**

Run:
```bash
cd backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -c "
from app.services.alerts.pattern_categories import classify_pattern
print(classify_pattern('isolated alkene', 'CHEMBL_DUNDEE'))
print(classify_pattern('Michael acceptor', 'BRENK'))
print(classify_pattern('ene_six_het_A(483)', 'PAINS_A'))
print(classify_pattern('nitro group', 'CHEMBL_DUNDEE'))
print(classify_pattern('Non-Hydrogen_atoms', 'ZINC'))
print(classify_pattern('2-halo pyridine', 'CHEMBL_DUNDEE'))
"
```

Expected:
```
Metabolic Liability
Reactive Group
Assay Interference
Toxicophore
Physicochemical
Unwanted Functionality
```

**Step 3: Commit**

```bash
git add backend/app/services/alerts/pattern_categories.py
git commit -m "feat: add pattern category classifier for structural alerts"
```

---

### Task 3: Enrich AlertResult and AlertManager with RDKit metadata

**Files:**
- Modify: `backend/app/services/alerts/alert_manager.py`

**Step 1: Add new fields to AlertResult dataclass**

Add `reference`, `scope`, `filter_set`, `catalog_description`, and `category` to the `AlertResult` dataclass. Update the `screen()` method to extract `Reference`, `Scope`, `FilterSet` from each `FilterCatalogEntry` and classify the pattern.

In `alert_manager.py`, add the import at the top:

```python
from .pattern_categories import classify_pattern
from .filter_catalog import AVAILABLE_CATALOGS
```

Update the `AlertResult` dataclass to add:

```python
@dataclass
class AlertResult:
    """Result from a structural alert pattern match."""

    pattern_name: str
    description: str
    severity: AlertSeverity
    matched_atoms: List[int]
    catalog_source: str
    smarts: Optional[str] = None
    reference: Optional[str] = None
    scope: Optional[str] = None
    filter_set: Optional[str] = None
    catalog_description: Optional[str] = None
    category: Optional[str] = None
```

**Step 2: Update `screen()` to extract RDKit entry properties and classify**

In the `screen()` method, after `pattern_name = entry.GetDescription() or "unknown_pattern"`, extract the entry properties:

```python
# Extract RDKit entry metadata
reference = None
scope = None
filter_set_id = None
try:
    if "Reference" in entry.GetPropList():
        reference = entry.GetProp("Reference")
    if "Scope" in entry.GetPropList():
        scope = entry.GetProp("Scope")
    if "FilterSet" in entry.GetPropList():
        filter_set_id = entry.GetProp("FilterSet")
except Exception:
    pass

# Get catalog description from AVAILABLE_CATALOGS
cat_key = catalog_type.upper()
cat_meta = AVAILABLE_CATALOGS.get(cat_key, {})
catalog_desc = cat_meta.get("name", cat_key)

# If entry-level reference is just a GitHub URL, prefer catalog-level reference
if reference and "github.com" in reference:
    catalog_ref = cat_meta.get("reference")
    if catalog_ref:
        reference = catalog_ref

# Classify pattern into concern category
category = classify_pattern(pattern_name, catalog_type)
```

Then update the `AlertResult` construction to pass these new fields:

```python
alert = AlertResult(
    pattern_name=pattern_name,
    description=description,
    severity=severity,
    matched_atoms=matched_atoms,
    catalog_source=catalog_type.upper(),
    smarts=smarts,
    reference=reference,
    scope=scope,
    filter_set=filter_set_id,
    catalog_description=catalog_desc,
    category=category,
)
```

**Step 3: Verify with morphine**

Run:
```bash
cd backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -c "
from rdkit import Chem
from app.services.alerts.alert_manager import AlertManager
mgr = AlertManager()
mol = Chem.MolFromSmiles('CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@@H]3[C@@H]1C5')
result = mgr.screen(mol, catalogs=['PAINS', 'BRENK', 'CHEMBL_DUNDEE'])
for a in result.alerts:
    print(f'Pattern: {a.pattern_name}')
    print(f'  Category: {a.category}')
    print(f'  Catalog: {a.catalog_description}')
    print(f'  Scope: {a.scope}')
    print(f'  Reference: {a.reference}')
    print()
"
```

Expected: Two alerts (BRENK + CHEMBL_DUNDEE) both showing "isolated alkene" / "isolated_alkene" with category "Metabolic Liability", proper catalog names, scope text, and Brenk citation.

**Step 4: Commit**

```bash
git add backend/app/services/alerts/alert_manager.py
git commit -m "feat: extract RDKit entry metadata and classify alert patterns"
```

---

### Task 4: Update schemas to include enriched alert fields

**Files:**
- Modify: `backend/app/schemas/alerts.py`

**Step 1: Add new fields to AlertResultSchema**

Add after the `smarts` field:

```python
reference: Optional[str] = Field(
    None, description="Citation for the original paper defining this alert"
)
scope: Optional[str] = Field(
    None, description="What aspect of compound quality this filter screens for"
)
filter_set: Optional[str] = Field(
    None, description="Filter set identifier from RDKit (e.g. ChEMBL23_Dundee)"
)
catalog_description: Optional[str] = Field(
    None, description="Human-readable name of the source catalog"
)
category: Optional[str] = Field(
    None,
    description="Concern category: Reactive Group, Metabolic Liability, Toxicophore, "
    "Assay Interference, Physicochemical, or Unwanted Functionality",
)
```

**Step 2: Add new fields to CatalogInfoSchema**

Add after the `note` field:

```python
reference: Optional[str] = Field(None, description="Citation for the original paper")
scope: Optional[str] = Field(None, description="What the filter set screens for")
doi: Optional[str] = Field(None, description="DOI link for the reference paper")
pmid: Optional[str] = Field(None, description="PubMed ID for the reference paper")
```

**Step 3: Extend catalog validator in AlertScreenRequest**

Update the `valid_catalogs` set in `validate_catalogs` to include ChEMBL catalogs:

```python
valid_catalogs = {
    "PAINS",
    "PAINS_A",
    "PAINS_B",
    "PAINS_C",
    "BRENK",
    "NIH",
    "ZINC",
    "CHEMBL_BMS",
    "CHEMBL_DUNDEE",
    "CHEMBL_GLAXO",
    "CHEMBL_INPHARMATICA",
    "CHEMBL_LINT",
    "CHEMBL_MLSMR",
    "CHEMBL_SURECHEMBL",
    "ALL",
}
```

**Step 4: Update the alerts route to pass new fields**

In `backend/app/api/routes/alerts.py`, update the `AlertResultSchema` construction in `screen_alerts()` to include the new fields:

```python
alerts = [
    AlertResultSchema(
        pattern_name=alert.pattern_name,
        description=alert.description,
        severity=AlertSeverity(alert.severity.value),
        matched_atoms=alert.matched_atoms,
        catalog_source=alert.catalog_source,
        smarts=alert.smarts,
        reference=alert.reference,
        scope=alert.scope,
        filter_set=alert.filter_set,
        catalog_description=alert.catalog_description,
        category=alert.category,
    )
    for alert in screening_result.alerts
]
```

Also update the `list_catalogs()` endpoint to pass the new catalog fields:

```python
catalogs[cat_type] = CatalogInfoSchema(
    name=cat_info["name"],
    description=cat_info["description"],
    pattern_count=cat_info["pattern_count"],
    severity=cat_info["severity"],
    note=cat_info.get("note"),
    reference=cat_info.get("reference"),
    scope=cat_info.get("scope"),
    doi=cat_info.get("doi"),
    pmid=cat_info.get("pmid"),
)
```

**Step 5: Verify API response**

Run the backend server and test:
```bash
cd backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -c "
from app.schemas.alerts import AlertResultSchema, CatalogInfoSchema
# Verify schemas accept new fields without error
a = AlertResultSchema(
    pattern_name='test', description='test', matched_atoms=[], catalog_source='PAINS',
    reference='Some ref', scope='Some scope', category='Reactive Group',
    catalog_description='PAINS Full Name'
)
print(a.model_dump())
"
```

**Step 6: Commit**

```bash
git add backend/app/schemas/alerts.py backend/app/api/routes/alerts.py
git commit -m "feat: add reference, scope, category fields to alert schemas and route"
```

---

### Task 5: Update batch tasks to include enriched alert fields

**Files:**
- Modify: `backend/app/services/batch/tasks.py`

**Step 1: Update alert dict construction in both batch functions**

In both `_process_single_molecule_simple()` (line ~132) and `_process_single_molecule()` (line ~335), update the alert dict to include the new fields. Change:

```python
{
    "catalog": a.catalog_source,
    "rule_name": a.pattern_name,
    "severity": (
        a.severity.value
        if hasattr(a.severity, "value")
        else str(a.severity)
    ),
    "matched_atoms": a.matched_atoms,
}
```

To:

```python
{
    "catalog": a.catalog_source,
    "rule_name": a.pattern_name,
    "severity": (
        a.severity.value
        if hasattr(a.severity, "value")
        else str(a.severity)
    ),
    "matched_atoms": a.matched_atoms,
    "reference": a.reference,
    "scope": a.scope,
    "catalog_description": a.catalog_description,
    "category": a.category,
}
```

**Step 2: Commit**

```bash
git add backend/app/services/batch/tasks.py
git commit -m "feat: include enriched alert metadata in batch processing results"
```

---

### Task 6: Update frontend TypeScript types

**Files:**
- Modify: `frontend/src/types/alerts.ts`

**Step 1: Add new fields to AlertResult interface**

Add after `smarts`:

```typescript
reference?: string | null;
scope?: string | null;
filter_set?: string | null;
catalog_description?: string | null;
category?: string | null;
```

**Step 2: Add new fields to CatalogInfo interface**

Add after `note`:

```typescript
reference?: string | null;
scope?: string | null;
doi?: string | null;
pmid?: string | null;
```

**Step 3: Commit**

```bash
git add frontend/src/types/alerts.ts
git commit -m "feat: add enriched alert fields to TypeScript types"
```

---

### Task 7: Redesign AlertCard to show enriched context

**Files:**
- Modify: `frontend/src/components/alerts/AlertCard.tsx`

**Step 1: Replace AlertCard implementation**

The new AlertCard shows:
1. Pattern name (formatted) + severity badge + category badge (color-coded)
2. Human-readable catalog name (not "CHEMBL_DUNDEE")
3. Scope description (what the filter set screens for)
4. Matched atoms with hover highlight
5. Reference citation (subtle, at bottom)
6. Approved drug examples (kept from existing code)

Category color mapping:
- Reactive Group → red
- Toxicophore → rose/pink
- Metabolic Liability → amber
- Assay Interference → purple
- Physicochemical → slate
- Unwanted Functionality → gray

Replace the full file with:

```tsx
import { useState } from 'react';
import { ChevronDown, ChevronUp } from 'lucide-react';
import type { AlertResult, AlertSeverity } from '../../types/alerts';

interface AlertCardProps {
  alert: AlertResult;
  onAtomHover?: (atoms: number[]) => void;
  className?: string;
}

const CATEGORY_STYLES: Record<string, { bg: string; text: string; label: string }> = {
  'Reactive Group': { bg: 'bg-red-100', text: 'text-red-700', label: 'Reactive Group' },
  'Toxicophore': { bg: 'bg-rose-100', text: 'text-rose-700', label: 'Toxicophore' },
  'Metabolic Liability': { bg: 'bg-amber-100', text: 'text-amber-700', label: 'Metabolic Liability' },
  'Assay Interference': { bg: 'bg-purple-100', text: 'text-purple-700', label: 'Assay Interference' },
  'Physicochemical': { bg: 'bg-slate-100', text: 'text-slate-700', label: 'Physicochemical' },
  'Unwanted Functionality': { bg: 'bg-gray-100', text: 'text-gray-600', label: 'Unwanted Functionality' },
};

/**
 * Known PAINS patterns that appear in FDA-approved drugs.
 * Provides educational context to help users understand alerts are warnings.
 */
const APPROVED_DRUG_EXAMPLES: Record<string, string[]> = {
  rhodanine: ['Methotrexate', 'Epalrestat'],
  catechol: ['Dopamine', 'Epinephrine', 'Levodopa'],
  quinone: ['Doxorubicin', 'Vitamin K'],
  michael_acceptor: ['Ibrutinib', 'Afatinib'],
  azo: ['Sulfasalazine', 'Phenazopyridine'],
  thiourea: ['Methimazole', 'Thiouracil'],
};

function getApprovedDrugNote(patternName: string): string | null {
  const patternLower = patternName.toLowerCase();
  for (const [key, drugs] of Object.entries(APPROVED_DRUG_EXAMPLES)) {
    if (patternLower.includes(key)) {
      return `Found in approved drugs: ${drugs.slice(0, 2).join(', ')}`;
    }
  }
  return null;
}

export function AlertCard({ alert, onAtomHover, className = '' }: AlertCardProps) {
  const [expanded, setExpanded] = useState(false);

  const getSeverityStyles = (severity: AlertSeverity) => {
    switch (severity) {
      case 'critical':
        return { bg: 'bg-red-50', border: 'border-red-200', badge: 'bg-red-100 text-red-800', icon: '!!!' };
      case 'warning':
        return { bg: 'bg-amber-50', border: 'border-amber-200', badge: 'bg-amber-100 text-amber-800', icon: '!' };
      case 'info':
        return { bg: 'bg-blue-50', border: 'border-blue-200', badge: 'bg-blue-100 text-blue-800', icon: 'i' };
      default:
        return { bg: 'bg-gray-50', border: 'border-gray-200', badge: 'bg-gray-100 text-gray-800', icon: '?' };
    }
  };

  const styles = getSeverityStyles(alert.severity);
  const categoryStyle = alert.category ? CATEGORY_STYLES[alert.category] || CATEGORY_STYLES['Unwanted Functionality'] : null;

  const formatPatternName = (name: string) => {
    const cleaned = name.replace(/\(\d+\)$/, '').trim();
    return cleaned
      .split('_')
      .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
      .join(' ');
  };

  const approvedDrugNote = getApprovedDrugNote(alert.pattern_name);

  return (
    <div
      className={`${styles.bg} border ${styles.border} rounded-lg p-4 transition-all hover:shadow-md ${className}`}
      onMouseEnter={() => onAtomHover?.(alert.matched_atoms)}
      onMouseLeave={() => onAtomHover?.([])}
    >
      <div className="flex items-start gap-3">
        <div
          className={`flex items-center justify-center w-8 h-8 rounded-full ${styles.badge} font-bold text-sm flex-shrink-0`}
        >
          {styles.icon}
        </div>
        <div className="flex-1 min-w-0">
          {/* Row 1: Pattern name + badges */}
          <div className="flex items-center gap-2 mb-1.5 flex-wrap">
            <h4 className="font-medium text-gray-900">
              {formatPatternName(alert.pattern_name)}
            </h4>
            <span className={`px-2 py-0.5 text-xs font-medium rounded-full ${styles.badge}`}>
              {alert.severity.toUpperCase()}
            </span>
            {categoryStyle && (
              <span className={`px-2 py-0.5 text-xs font-medium rounded-full ${categoryStyle.bg} ${categoryStyle.text}`}>
                {categoryStyle.label}
              </span>
            )}
          </div>

          {/* Row 2: Human-readable catalog name */}
          <p className="text-sm text-gray-600 mb-1.5">
            {alert.catalog_description || alert.catalog_source}
          </p>

          {/* Row 3: Scope — what this filter screens for */}
          {alert.scope && (
            <p className="text-sm text-gray-500 italic mb-2">{alert.scope}</p>
          )}

          {/* Matched atoms */}
          {alert.matched_atoms.length > 0 && (
            <div className="text-xs text-gray-500 mb-2">
              <span className="font-medium">Matched atoms:</span>{' '}
              {alert.matched_atoms.join(', ')}
              <span className="ml-2 text-amber-600">(hover to highlight)</span>
            </div>
          )}

          {/* Approved drug note */}
          {approvedDrugNote && (
            <div className="text-xs text-amber-700 bg-yellow-50 rounded px-2 py-1 inline-block mb-2">
              {approvedDrugNote}
            </div>
          )}

          {/* Expandable reference section */}
          {alert.reference && (
            <button
              onClick={(e) => { e.stopPropagation(); setExpanded(!expanded); }}
              className="flex items-center gap-1 text-xs text-gray-400 hover:text-gray-600 transition-colors"
            >
              {expanded ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
              {expanded ? 'Hide reference' : 'Show reference'}
            </button>
          )}
          {expanded && alert.reference && (
            <div className="mt-2 text-xs text-gray-500 bg-white/50 rounded p-2 border border-gray-100">
              <span className="font-medium">Reference:</span> {alert.reference}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
```

**Step 2: Verify TypeScript compiles**

Run: `cd frontend && npx tsc --noEmit`

Expected: No errors related to AlertCard.

**Step 3: Commit**

```bash
git add frontend/src/components/alerts/AlertCard.tsx
git commit -m "feat: redesign AlertCard with catalog context, category badges, and citations"
```

---

### Task 8: Add collapsible ChEMBL catalog selector to SingleValidation

**Files:**
- Modify: `frontend/src/pages/SingleValidation.tsx`

**Step 1: Add ChEMBL catalog state and toggle logic**

Near the existing `selectedCatalogs` state (line ~226), add state for ChEMBL group expansion:

```typescript
const [chemblExpanded, setChemblExpanded] = useState(false);
```

Define the ChEMBL catalogs:

```typescript
const CHEMBL_CATALOGS = [
  { id: 'CHEMBL_BMS', label: 'BMS HTS Filters' },
  { id: 'CHEMBL_DUNDEE', label: 'Dundee NTD Filters' },
  { id: 'CHEMBL_GLAXO', label: 'Glaxo Hard Filters' },
  { id: 'CHEMBL_INPHARMATICA', label: 'Inpharmatica' },
  { id: 'CHEMBL_LINT', label: 'Lilly MedChem (LINT)' },
  { id: 'CHEMBL_MLSMR', label: 'NIH MLSMR' },
  { id: 'CHEMBL_SURECHEMBL', label: 'SureChEMBL' },
];
```

Add a helper to toggle all ChEMBL catalogs:

```typescript
const toggleAllChembl = (enabled: boolean) => {
  const chemblIds = CHEMBL_CATALOGS.map((c) => c.id);
  if (enabled) {
    setSelectedCatalogs((prev) => [...new Set([...prev, ...chemblIds])]);
  } else {
    setSelectedCatalogs((prev) => prev.filter((c) => !chemblIds.includes(c)));
  }
};

const allChemblSelected = CHEMBL_CATALOGS.every((c) => selectedCatalogs.includes(c.id));
const someChemblSelected = CHEMBL_CATALOGS.some((c) => selectedCatalogs.includes(c.id));
```

**Step 2: Replace the catalog selector JSX**

Replace the existing catalog selector section (lines ~1194-1211) with:

```tsx
{/* Catalog selector */}
<div>
  <p className="text-xs text-[var(--color-text-muted)] mb-2">Select catalogs to screen:</p>
  {/* Core catalogs */}
  <div className="flex flex-wrap gap-2 mb-2">
    {['PAINS', 'BRENK', 'NIH', 'ZINC'].map((catalog) => (
      <button
        key={catalog}
        onClick={() => toggleCatalog(catalog)}
        className={cn(
          'px-3 py-1.5 text-sm rounded-lg transition-all',
          selectedCatalogs.includes(catalog)
            ? 'bg-[var(--color-primary)]/15 text-[var(--color-primary)] border border-[var(--color-primary)]/30'
            : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] border border-transparent hover:border-[var(--color-border)]'
        )}
      >
        {catalog}
      </button>
    ))}
  </div>
  {/* ChEMBL group toggle */}
  <div className="border border-[var(--color-border)] rounded-lg overflow-hidden">
    <button
      onClick={() => setChemblExpanded(!chemblExpanded)}
      className="w-full flex items-center justify-between px-3 py-2 text-sm hover:bg-[var(--color-surface-sunken)] transition-colors"
    >
      <div className="flex items-center gap-2">
        <input
          type="checkbox"
          checked={allChemblSelected}
          ref={(el) => { if (el) el.indeterminate = someChemblSelected && !allChemblSelected; }}
          onChange={(e) => toggleAllChembl(e.target.checked)}
          onClick={(e) => e.stopPropagation()}
          className="rounded border-gray-300 text-[var(--color-primary)] focus:ring-[var(--color-primary)]"
        />
        <span className="text-[var(--color-text-secondary)]">ChEMBL Pharma Filters</span>
        <span className="text-xs text-[var(--color-text-muted)]">(7 filter sets)</span>
      </div>
      <svg
        className={cn('w-4 h-4 text-[var(--color-text-muted)] transition-transform', chemblExpanded && 'rotate-180')}
        fill="none" stroke="currentColor" viewBox="0 0 24 24"
      >
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
      </svg>
    </button>
    {chemblExpanded && (
      <div className="px-3 pb-2 pt-1 border-t border-[var(--color-border)] flex flex-wrap gap-2">
        {CHEMBL_CATALOGS.map((cat) => (
          <button
            key={cat.id}
            onClick={() => toggleCatalog(cat.id)}
            className={cn(
              'px-2.5 py-1 text-xs rounded-md transition-all',
              selectedCatalogs.includes(cat.id)
                ? 'bg-[var(--color-primary)]/15 text-[var(--color-primary)] border border-[var(--color-primary)]/30'
                : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] border border-transparent hover:border-[var(--color-border)]'
            )}
          >
            {cat.label}
          </button>
        ))}
      </div>
    )}
  </div>
</div>
```

**Step 3: Verify TypeScript compiles**

Run: `cd frontend && npx tsc --noEmit`

**Step 4: Commit**

```bash
git add frontend/src/pages/SingleValidation.tsx
git commit -m "feat: add collapsible ChEMBL pharma filter selector to alerts tab"
```

---

### Task 9: Update AlertResults group headers with catalog context

**Files:**
- Modify: `frontend/src/components/alerts/AlertResults.tsx`

**Step 1: Add catalog context to group headers**

When alerts are grouped by catalog, the group header should show the human-readable catalog name and scope from the first alert in that group. Update the grouped view section to use `catalog_description` and `scope` from the alert data:

In the grouped view map (around line 153), replace the group header:

```tsx
{Object.entries(alertsByCatalog).map(([catalog, catalogAlerts]) => {
  const firstAlert = catalogAlerts[0];
  return (
    <div key={catalog}>
      <div className="mb-3">
        <h4 className="font-medium text-gray-700 flex items-center gap-2">
          <span>{firstAlert?.catalog_description || catalog}</span>
          <span className="text-sm font-normal text-gray-500">
            ({catalogAlerts.length} alert{catalogAlerts.length !== 1 ? 's' : ''})
          </span>
        </h4>
        {firstAlert?.scope && (
          <p className="text-xs text-gray-500 mt-0.5 italic">{firstAlert.scope}</p>
        )}
      </div>
      <div className="space-y-3">
        {catalogAlerts.map((alert, index) => (
          <AlertCard
            key={`${alert.pattern_name}-${index}`}
            alert={alert}
            onAtomHover={onHighlightAtoms}
          />
        ))}
      </div>
    </div>
  );
})}
```

**Step 2: Verify TypeScript compiles**

Run: `cd frontend && npx tsc --noEmit`

**Step 3: Commit**

```bash
git add frontend/src/components/alerts/AlertResults.tsx
git commit -m "feat: show catalog context in alert group headers"
```

---

### Task 10: Run full lint and type checks

**Files:** None (verification only)

**Step 1: Backend lint**

Run: `cd backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m ruff check .`

Fix any issues found.

**Step 2: Frontend type check**

Run: `cd frontend && npx tsc --noEmit`

Fix any issues found.

**Step 3: Manual verification — test with morphine**

Start the backend and hit the alerts endpoint with morphine and multiple catalogs to verify the full enriched response:

```bash
cd backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -c "
from rdkit import Chem
from app.services.alerts.alert_manager import AlertManager
mgr = AlertManager()
mol = Chem.MolFromSmiles('CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@@H]3[C@@H]1C5')
result = mgr.screen(mol, catalogs=['PAINS', 'BRENK', 'CHEMBL_DUNDEE', 'CHEMBL_GLAXO', 'CHEMBL_LINT'])
print(f'Total alerts: {result.total_alerts}')
for a in result.alerts:
    print(f'  {a.pattern_name} | {a.category} | {a.catalog_description}')
    print(f'    Scope: {a.scope}')
    print(f'    Ref: {a.reference}')
"
```

**Step 4: Final commit if any fixes were needed**

```bash
git add -A && git commit -m "fix: lint and type check fixes for alert enrichment"
```
