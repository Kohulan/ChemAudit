<div align="center">

<img src="assets/logo.png" alt="ChemVault" width="80" />

# Getting Started

### Your First Steps with ChemVault

</div>

---

## 📋 Table of Contents

- [Installation](#-installation)
- [Your First Validation](#-your-first-validation)
- [Understanding Results](#-understanding-results)
- [Next Steps](#-next-steps)

---

## 💻 Installation

### Option 1: Docker (Recommended)

The fastest way to get started:

```bash
# Clone the repository
git clone https://github.com/yourusername/chemvault.git
cd chemvault

# Start all services (development mode)
docker-compose up -d

# Wait for services to be ready (~30 seconds)
docker-compose logs -f
```

🌐 **Access Points:**

| Service | URL |
|---------|-----|
| **Web Interface** | http://localhost:3000 |
| **API Documentation** | http://localhost:8000/docs |

### Option 1b: Production Deployment

For production use with configurable batch limits:

```bash
# Interactive deployment - select a profile
./deploy.sh

# Or specify profile directly
./deploy.sh medium   # 10K molecules, 4 workers
./deploy.sh large    # 50K molecules, 8 workers
./deploy.sh coconut  # 1M molecules, 16 workers
```

See [Deployment Guide](DEPLOYMENT.md) for full production setup.

### Option 2: Manual Installation

<details>
<summary><b>Backend Setup</b></summary>

**Requirements:** Python 3.11+, Poetry

```bash
cd backend

# Install dependencies
poetry install

# Start development server
poetry run uvicorn app.main:app --reload --port 8000
```

</details>

<details>
<summary><b>Frontend Setup</b></summary>

**Requirements:** Node.js 18+, npm

```bash
cd frontend

# Install dependencies
npm install

# Start development server
npm run dev
```

</details>

---

## 🧪 Your First Validation

### Using the Web Interface

1. **Open** http://localhost:3000
2. **Enter** a SMILES string in the input field:
   ```
   CC(=O)Oc1ccccc1C(=O)O
   ```
   *(This is Aspirin)*
3. **Click** "Validate"
4. **View** your results!

### Using the API

```bash
curl -X POST http://localhost:8000/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)O",
    "format": "smiles"
  }'
```

<details>
<summary><b>📄 Example Response</b></summary>

```json
{
  "valid": true,
  "validation_score": 95,
  "checks": [
    {
      "name": "valence",
      "passed": true,
      "severity": "critical",
      "message": "All atoms have valid valence"
    },
    {
      "name": "aromaticity",
      "passed": true,
      "severity": "warning",
      "message": "Aromatic system is valid"
    }
  ],
  "standardized_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
  "inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
}
```

</details>

### Using Python

```python
import requests

response = requests.post(
    "http://localhost:8000/api/v1/validate",
    json={
        "molecule": "CC(=O)Oc1ccccc1C(=O)O",
        "format": "smiles"
    }
)

result = response.json()
print(f"Valid: {result['valid']}")
print(f"Score: {result['validation_score']}")
```

---

## 📊 Understanding Results

### Validation Score

| Score Range | Quality | Recommendation |
|-------------|---------|----------------|
| **90-100** | 🟢 Excellent | Ready for use |
| **70-89** | 🟡 Good | Minor issues, review recommended |
| **50-69** | 🟠 Fair | Needs attention |
| **0-49** | 🔴 Poor | Significant issues |

### Check Severities

| Severity | Icon | Meaning |
|----------|------|---------|
| **Critical** | 🔴 | Must be fixed - structure is invalid |
| **Warning** | 🟡 | Should be reviewed - may affect results |
| **Info** | 🔵 | Informational - no action required |

### Common Validation Checks

| Check | What It Does |
|-------|--------------|
| **Valence** | Verifies all atoms have correct number of bonds |
| **Aromaticity** | Validates aromatic ring systems |
| **Stereochemistry** | Checks stereo centers and E/Z bonds |
| **Kekulization** | Tests if aromatic rings can be kekulized |
| **Connectivity** | Ensures molecule is properly connected |

---

## 🎓 Next Steps

Now that you've completed your first validation, explore more features:

<table>
<tr>
<td width="50%" valign="top">

### 📚 Learn More

- [User Guide](USER_GUIDE.md) - Complete feature walkthrough
- [API Reference](API_REFERENCE.md) - Full API documentation
- [Troubleshooting](TROUBLESHOOTING.md) - Common issues & solutions

</td>
<td width="50%" valign="top">

### 🚀 Try Features

- **Batch Processing** - Validate thousands of molecules
- **Structural Alerts** - Screen for PAINS/BRENK
- **ML-Readiness** - Assess ML suitability
- **Database Lookup** - Search PubChem, ChEMBL

</td>
</tr>
</table>

---

## 🧬 Sample Molecules to Try

| Name | SMILES | Description |
|------|--------|-------------|
| Aspirin | `CC(=O)Oc1ccccc1C(=O)O` | Common pain reliever |
| Caffeine | `Cn1cnc2c1c(=O)n(c(=O)n2C)C` | Stimulant |
| Ibuprofen | `CC(C)Cc1ccc(cc1)C(C)C(=O)O` | Anti-inflammatory |
| Penicillin G | `CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O` | Antibiotic |
| Morphine | `CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O` | Analgesic |

---

<div align="center">

**Questions?** Check the [FAQ](TROUBLESHOOTING.md) or [Open an Issue](https://github.com/yourusername/chemvault/issues)

</div>
