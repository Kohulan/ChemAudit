---
sidebar_position: 7
title: IUPAC Name Conversion
description: Convert IUPAC chemical names to SMILES automatically
---

# IUPAC Name Conversion

ChemAudit can accept IUPAC names, common drug names, and trade names as input — converting them to SMILES automatically before validation. This makes it easy to validate molecules without knowing their SMILES notation.

## How It Works

### Auto-Detection

When you enter text in the input field, ChemAudit classifies it as one of three types:

| Detection | Indicators | Example |
|-----------|-----------|---------|
| **SMILES** | Contains `()`, `[]`, `=`, `#`, `@`, `/`, `\` | `CC(=O)Oc1ccccc1C(=O)O` |
| **IUPAC** | Contains spaces, hyphens, or chemical suffixes (-ol, -ane, -ene, -acid, -one, -al, -amine, -amide, -yl, -ide, -oxy, -oic) | `2-acetoxybenzoic acid` |
| **Ambiguous** | Purely alphabetic, ≥4 characters, no clear indicators | `aspirin` |

An input type badge appears below the input field showing the detected type (e.g., "Detected: IUPAC name"). When detection is ambiguous, you can use force-override buttons to specify the format manually.

### Conversion Pipeline

IUPAC names are converted to SMILES using a two-tier approach:

1. **OPSIN** (primary) — Open Parser for Systematic IUPAC Nomenclature. Runs locally via Java, fast and works offline. Handles systematic IUPAC names with high accuracy.

2. **PubChem PUG REST** (fallback) — If OPSIN cannot parse the name, ChemAudit queries the PubChem API to resolve common names and trade names.

The conversion source (OPSIN or PubChem) is reported in the results.

### What Happens After Conversion

Once a name is converted to SMILES, validation proceeds normally — all tabs work as expected (Validation, Alerts, Scoring, Standardization, Database Lookup, Scoring Profiles). A conversion badge shows the original name and the converted SMILES for reference.

## Using IUPAC Input

### Web Interface

1. Enter a chemical name in the input field (e.g., `aspirin` or `2-acetoxybenzoic acid`)
2. The input type badge shows "Detected: IUPAC name"
3. Click **Validate**
4. The conversion badge shows: `aspirin → CC(=O)Oc1ccccc1C(=O)O (via PubChem)`
5. All validation results are based on the converted SMILES

:::tip Placeholder Help
The input field shows: "Enter SMILES, InChI, or IUPAC name (e.g., aspirin)"
:::

### API

Send the name directly to the validation endpoint:

```bash
curl -X POST http://localhost:8001/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "aspirin",
    "format": "auto"
  }'
```

The response includes an `input_interpretation` field with conversion details:

```json
{
  "status": "completed",
  "input_interpretation": {
    "detected_as": "iupac",
    "original_input": "aspirin",
    "converted_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "conversion_source": "pubchem"
  },
  "molecule_info": {
    "input_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "canonical_smiles": "CC(=O)OC1=CC=CC=C1C(O)=O",
    "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
  },
  "overall_score": 95,
  "issues": []
}
```

You can also force the input type:

```bash
curl -X POST http://localhost:8001/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "morphine",
    "input_type": "iupac"
  }'
```

## Supported Names

| Name Type | Examples | Conversion Source |
|-----------|---------|-------------------|
| **Systematic IUPAC** | 2-acetoxybenzoic acid, ethanol, propan-2-one | OPSIN |
| **Common names** | aspirin, caffeine, morphine | PubChem |
| **Trade names** | Tylenol, Advil (if in PubChem) | PubChem |
| **CAS-style** | Names resolvable by PubChem | PubChem |

## Limitations

- **OPSIN requires Java**: The Docker image includes a Java JRE. For local development, ensure Java is installed and `OPSIN_JAR_PATH` is configured.
- **PubChem fallback requires internet**: If OPSIN cannot parse the name and no internet connection is available, conversion will fail.
- **Obscure names may fail**: Very specialized or proprietary names not in PubChem's database cannot be resolved.
- **Ambiguous names**: Some short names may be ambiguous between SMILES and IUPAC. Use the force-override buttons or set `input_type` in the API to resolve.
- **PubChem rate limits**: The PubChem API has usage limits. For batch conversion of many names, consider pre-converting using a dedicated tool.

:::warning Verify Conversions
Always verify that the converted SMILES matches your intended molecule, especially for ambiguous or common names that may map to unexpected structures.
:::

## Next Steps

- **[Single Validation](/docs/user-guide/single-validation)** — Full validation features for converted molecules
- **[Database Integrations](/docs/user-guide/database-integrations)** — Cross-reference converted molecules against databases
- **[Batch Processing](/docs/user-guide/batch-processing)** — Process datasets (batch input uses SMILES columns)
