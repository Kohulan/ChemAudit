---
sidebar_position: 16
title: Diagnostics
description: Low-level chemical structure analysis tools for debugging and pre-screening
---

# Diagnostics

The Diagnostics page provides low-level chemical structure analysis tools for debugging parsing errors, investigating structural discrepancies, and pre-screening files before batch processing.

Access Diagnostics from the **Data Preparation** dropdown in the header navigation.

## Available Tools

### SMILES Diagnostics

Parse and analyze a SMILES string to understand its structure at the atom level:

- **Atom-by-atom breakdown**: Each atom with its element, charge, isotope, and bonding
- **Ring closures**: Identification of ring opening/closing tokens
- **Branch points**: Parenthetical branch structure
- **Parse errors**: Detailed error messages for invalid SMILES with the position of failure

:::tip Debugging Invalid SMILES
When a molecule fails validation with "Cannot parse input," paste it into SMILES Diagnostics to see exactly where and why parsing fails.
:::

### InChI Layer Diff

Compare InChI layers between two molecules to identify exact structural differences:

- **Formula layer** — Molecular formula comparison
- **Connectivity layer** — Atom connectivity differences
- **Hydrogen layer** — Hydrogen attachment differences
- **Charge layer** — Formal charge differences
- **Stereo layers** — Stereochemistry differences (tetrahedral, double bond)

This is useful for understanding why two molecules that look similar have different InChIKeys or why standardization changed a structure.

### Round-Trip Validation

Validate SMILES → MOL → SMILES round-trip fidelity:

1. Parse the input SMILES into an RDKit molecule object
2. Convert to MOL block (3D coordinates if possible)
3. Convert back to SMILES
4. Compare the original and round-tripped SMILES

Differences indicate potential representation loss — the molecule may not survive serialization/deserialization cycles in some workflows.

### File Pre-Validation

Check a file for parseable SMILES before uploading to batch processing:

- **File parsing**: Validates that the file format is correct (CSV, SDF)
- **Column detection**: Identifies SMILES and name columns
- **Parse rate**: Reports what percentage of SMILES can be parsed by RDKit
- **Error summary**: Lists unparseable entries with row numbers and error messages

:::info Pre-Screening Saves Time
Running file pre-validation before a large batch upload helps you fix format issues and unparseable entries upfront, rather than discovering them during processing.
:::

### Coordinate Dimension Analysis

Detect the coordinate dimensionality of molecules in MOL blocks or SDF files:

| Dimension | Meaning |
|-----------|---------|
| **2D** | X/Y coordinates only (Z = 0) — suitable for 2D depiction |
| **3D** | X/Y/Z coordinates — includes conformational information |
| **No coordinates** | No atom coordinate data present |
| **Mixed** | File contains molecules with different dimensionalities |

This helps verify that SDF files have the expected coordinate data before using them for 3D analysis or docking workflows.

## Use Cases

### Debugging Validation Failures

1. A molecule fails validation → copy the SMILES
2. Open Diagnostics → SMILES Diagnostics
3. Paste the SMILES → see the exact parsing error or structural issue
4. Fix the SMILES and re-validate

### Comparing Standardized vs. Original

1. Standardize a molecule → note the original and standardized SMILES
2. Open Diagnostics → InChI Layer Diff
3. Enter both SMILES → see exactly which structural layers changed
4. Understand whether stereochemistry, connectivity, or charge was affected

### Pre-Screening Batch Files

1. Open Diagnostics → File Pre-Validation
2. Upload your CSV or SDF
3. Review the parse rate and error summary
4. Fix problematic entries
5. Upload the cleaned file to Batch Processing

## Next Steps

- **[Single Validation](/docs/user-guide/single-validation)** — Validate individual molecules
- **[Batch Processing](/docs/user-guide/batch-processing)** — Process large datasets
- **[Standardization](/docs/user-guide/standardization)** — Understand the standardization pipeline
- **[QSAR-Ready Pipeline](/docs/user-guide/qsar-ready)** — Prepare datasets for ML
