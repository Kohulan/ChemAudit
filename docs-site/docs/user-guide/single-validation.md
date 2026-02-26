---
sidebar_position: 1
title: Single Validation
description: Validate individual molecules with comprehensive structural checks
---

# Single Molecule Validation

ChemAudit's single molecule validation provides comprehensive structural analysis for individual molecules. This feature is ideal for quick checks, exploring new compounds, or validating structures before batch processing.

## Supported Input Formats

ChemAudit automatically detects and supports multiple input formats:

| Format | Example | Auto-Detected |
|--------|---------|---------------|
| **SMILES** | `CCO` | Yes |
| **InChI** | `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3` | Yes |
| **MOL Block** | V2000/V3000 format | Yes |
| **IUPAC Name** | `aspirin`, `2-acetoxybenzoic acid` | Yes |

:::tip Auto-Detection
Simply paste your molecule in any format — including IUPAC names like "aspirin" or "2-acetoxybenzoic acid". ChemAudit automatically detects the input type and converts names to SMILES. See [IUPAC Name Conversion](/docs/user-guide/iupac-conversion) for details.
:::

## How to Validate

### Using the Web Interface

1. Navigate to the **Single Validation** page (home)
2. Enter or paste your molecule in the input field
3. Click **Validate**
4. Review results across all tabs:
   - **Validation**: Structural checks and overall score
   - **Alerts**: PAINS, BRENK, NIH, ZINC, ChEMBL screening
   - **Scoring**: ML-readiness, drug-likeness, ADMET, NP-likeness
   - **Scoring Profiles**: Consensus score, lead/fragment-likeness, property breakdowns, bioavailability radar
   - **Standardization**: ChEMBL-compatible cleanup with provenance timeline
   - **Database Lookup**: PubChem, ChEMBL, COCONUT cross-references

### Using the API

```bash
curl -X POST http://localhost:8001/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CCO",
    "format": "auto"
  }'
```

## Validation Checks Explained

ChemAudit runs 5 basic checks on every molecule plus 17 deep validation checks organized into three domains. All check severities can be customized through the severity configuration panel.

### Basic Checks

These fundamental checks assess structural validity:

| Check | Severity | Description | Common Causes of Failure |
|-------|----------|-------------|-------------------------|
| **Parsability** | Critical | Can the input be parsed into a valid molecule? | Malformed SMILES, invalid characters |
| **Sanitization** | Error | Does the molecule pass RDKit sanitization? | Structural inconsistencies, invalid atom types |
| **Valence** | Critical | Do all atoms have chemically valid bond counts? | Typos in SMILES, incorrect charges |
| **Aromaticity** | Error | Can aromatic systems be kekulized? | Invalid aromatic systems, wrong electron count |
| **Connectivity** | Warning | Is the molecule a single connected component? | Mixtures, salts, disconnected fragments |

:::danger Critical Failures
If a critical check fails, the molecule structure is invalid and cannot be used for further analysis.
:::

### Deep Validation — Chemical Composition

Six checks examining what the molecule is made of:

| Check | Severity | Description |
|-------|----------|-------------|
| **Mixture Detection** | Warning | Identifies disconnected fragments and classifies each as drug, salt, solvent, or unknown |
| **Solvent Contamination** | Warning | Matches fragments against 15+ known solvents (water, DMSO, DMF, acetonitrile, methanol, ethanol, etc.) |
| **Inorganic Filter** | Warning/Error | Detects inorganic (no carbon → Error) or organometallic (carbon + metal → Warning) compounds |
| **Radical Detection** | Warning | Flags atoms with unpaired electrons |
| **Isotope Labels** | Info | Detects isotope-labeled atoms (deuterium, ¹³C, tritium, etc.) |
| **Trivial Molecule** | Error | Flags molecules with ≤ 3 heavy atoms as too small for meaningful analysis |

### Deep Validation — Structural Complexity

Six checks assessing structural features that may complicate analysis:

| Check | Severity | Description |
|-------|----------|-------------|
| **Hypervalent Atoms** | Warning | Atoms exceeding maximum allowed valence for their element |
| **Polymer Detection** | Info | Detected via SGroup markers, MW > 1500 Da, or dummy atoms |
| **Ring Strain** | Warning | 3- or 4-membered rings with significant angle strain |
| **Macrocycle Detection** | Info | Rings with > 12 atoms |
| **Charged Species** | Info | Net charge calculation, positive/negative atoms, zwitterion detection |
| **Explicit Hydrogen Audit** | Info | Mixed explicit/implicit hydrogen representation |

### Deep Validation — Stereo & Tautomers

Five checks related to stereochemistry and tautomeric forms:

| Check | Severity | Description |
|-------|----------|-------------|
| **Stereoisomer Enumeration** | Warning | Enumerates possible stereoisomers from undefined centers (cap: 128) |
| **Undefined Stereocenters** | Warning | Counts chiral centers without R/S specification |
| **Tautomer Detection** | Info | Enumerates tautomers and checks if input matches canonical form |
| **Aromatic System Validation** | Warning | Flags unusual aromatic ring sizes (not 5 or 6) and charged aromatic atoms |
| **Coordinate Dimension** | Info | Detects 2D, 3D, or no coordinate data |

:::tip Severity Customization
All deep validation check severities can be overridden through the severity configuration panel. This lets you adjust which checks are treated as errors vs. warnings vs. informational based on your specific use case.
:::

## Validation Options

### Preserve Aromatic Notation

By default, ChemAudit outputs canonical SMILES in kekulized form (explicit single/double bonds). Enable this option to preserve aromatic notation:

**Input:**
```
c1ccccc1
```

**Output (default - kekulized):**
```
C1=CC=CC=C1
```

**Output (preserve_aromatic=true):**
```
c1ccccc1
```

## Understanding Results

### Overall Score

The overall validation score ranges from 0-100 and indicates structure quality:

| Score Range | Quality | Interpretation |
|-------------|---------|---------------|
| **90-100** | Excellent | Structure is valid and ready for use |
| **70-89** | Good | Minor issues detected, review recommended |
| **50-69** | Fair | Significant issues need attention |
| **0-49** | Poor | Critical problems, structure likely invalid |

### Molecule Information

Every validated molecule returns comprehensive information:

- **Input SMILES**: Original input
- **Canonical SMILES**: Standardized SMILES representation
- **InChI**: International Chemical Identifier
- **InChIKey**: Hashed InChI for database lookups
- **Molecular Formula**: Element composition
- **Molecular Weight**: Exact mass
- **Atom Count**: Number of heavy atoms

### Issue Details

Failed checks appear in the issues list with:

- **Check name**: Which validation failed
- **Severity**: Critical, Warning, or Info
- **Message**: Human-readable description
- **Affected atoms**: Atom indices involved (if applicable)
- **Details**: Additional technical information

## Scoring Profiles Tab

The Scoring Profiles tab provides advanced drug-likeness analysis beyond the standard Scoring tab:

- **Consensus Score** — Drug-likeness consensus across 5 rule sets (Lipinski, Ghose, Veber, Egan, Muegge), scored 0–5
- **Lead & Fragment Likeness** — Lead-likeness assessment, Rule of 3 compliance, salt inventory, ligand efficiency
- **Property Breakdown** — Per-atom TPSA and LogP breakdowns, Bertz complexity index, Fsp3 detail
- **Bioavailability Radar** — 6-axis radar chart for oral bioavailability assessment plus BOILED-Egg scatter for GI absorption and BBB permeation predictions
- **Atom Contribution Viewer** — Per-atom property contribution heatmap

See [Scoring Profiles](/docs/user-guide/scoring/profiles) for full details on profile scoring and the custom profile builder.

## Standardization Provenance

When viewing the Standardization tab, a provenance timeline shows exactly what changed at each pipeline stage:

- **Vertical timeline** with per-stage cards that auto-expand when changes occurred
- **Change types tracked**: charge normalization, bond normalization, ring aromaticity, radical changes, fragment removal
- **DVAL cross-references**: Links to deep validation findings (e.g., "DVAL-01: 2 undefined stereocenters detected")
- **Tautomer summary**: Number of tautomers enumerated, canonical form, complexity flag
- **Stereo detail**: Per-center before/after configuration and reason for change

Enable provenance in the API with `include_provenance=true` on the `/standardize` endpoint. See [Standardization](/docs/user-guide/standardization) for details.

## Bookmarking Results

After validation, click the **Bookmark** button (star icon) in the results header to save the molecule and its full result snapshot. Bookmarked molecules appear on the [Bookmarks & History](/docs/user-guide/bookmarks-history) page where you can search, tag, and revisit them.

## Severity Configuration

Click the gear icon to open the severity configuration panel. This lets you override the default severity of any validation check:

- Choose **Error**, **Warning**, or **Info** for each check
- Reset individual checks back to their defaults
- Overrides persist in your browser's localStorage
- The overall verdict dynamically recomputes based on your effective severities

## RDKit Version Source

The canonical SMILES shown in molecule info includes a tooltip indicating which RDKit version was used for canonicalization. This is useful for reproducibility when sharing results.

## Common Validation Errors

### Valence Errors

**Problem**: Atom has incorrect number of bonds

**Examples:**

| Invalid | Valid | Explanation |
|---------|-------|-------------|
| `CC(C)(C)(C)C` | `CC(C)(C)C` | Carbon can have max 4 bonds |
| `CN` (with 4 bonds to N) | `C[N+]` | Quaternary nitrogen needs charge |

### Kekulization Failures

**Problem**: Aromatic ring doesn't follow Hückel's rule (4n+2 electrons)

**Examples:**

| Invalid | Valid | Explanation |
|---------|-------|-------------|
| `c1cccc1` | `c1ccccc1` | Benzene needs 6 carbons (6 pi electrons) |
| `c1cccccccc1` | `C1=CC=CC=CC=CC=C1` | 8 pi electrons - not aromatic |

### Unclosed Rings

**Problem**: Ring opening doesn't have matching closing

**Examples:**

| Invalid | Valid | Explanation |
|---------|-------|-------------|
| `C1CCC` | `C1CCCC1` | Ring 1 not closed |
| `C1CCC2` | `C1CCC2CC2C1` | Rings 1 and 2 must both close |

## Best Practices

1. **Validate early**: Check structures before starting experiments
2. **Review warnings**: Don't ignore stereo or representation warnings
3. **Use canonical forms**: Work with canonical SMILES for consistency
4. **Cross-check databases**: Use database lookup to verify compound identity
5. **Document issues**: Record any validation warnings in your data

## Next Steps

- **[Batch Processing](/docs/user-guide/batch-processing)** — Validate thousands of molecules at once
- **[Bookmarks & History](/docs/user-guide/bookmarks-history)** — Save and revisit validation results
- **[IUPAC Name Conversion](/docs/user-guide/iupac-conversion)** — Enter chemical names directly
- **[Scoring Profiles](/docs/user-guide/scoring/profiles)** — Custom property scoring
- **[Structural Alerts](/docs/user-guide/structural-alerts)** — Screen for problematic patterns
- **[Scoring](/docs/user-guide/scoring/overview)** — Comprehensive molecular scoring
- **[API Reference](/docs/api/endpoints)** — Detailed API documentation
