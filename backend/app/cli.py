"""
ChemAudit CLI - Chemical structure validation, scoring, and standardization.

Usage:
    chemaudit validate --smiles "CCO"
    chemaudit score --smiles "CCO"
    chemaudit standardize --smiles "CCO"
    chemaudit profile --smiles "CCO"
    echo "CCO" | chemaudit validate
    chemaudit validate --smiles "CCO" --local    # offline, no server
    chemaudit validate --smiles "CCO" --format json
"""

import json
import sys
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(
    name="chemaudit",
    help="ChemAudit CLI - Chemical structure validation, scoring, and standardization",
    no_args_is_help=True,
)
console = Console()


def _get_smiles(smiles: Optional[str], file: Optional[str]) -> str:
    """Resolve SMILES input from --smiles flag, --file, or stdin pipe.

    Args:
        smiles: Direct SMILES string from --smiles flag.
        file: Path to an input file (CSV or SDF).

    Returns:
        Resolved SMILES or molecule string.

    Raises:
        typer.BadParameter: If no input source provides a valid molecule string.
    """
    if smiles:
        return smiles

    if file:
        with open(file) as f:
            lines = f.readlines()
        if not lines:
            raise typer.BadParameter(f"File {file} is empty")
        first_line = lines[0].strip()
        # If first line looks like a CSV header containing "smiles", use second line
        if first_line.lower().find("smiles") != -1 and len(lines) > 1:
            return lines[1].strip().split(",")[0].strip().strip('"').strip("'")
        # For SDF files (contain V2000/V3000 markers), read entire content
        content = "".join(lines)
        if "V2000" in content or "V3000" in content:
            return content
        return first_line

    if not sys.stdin.isatty():
        data = sys.stdin.read().strip()
        if data:
            return data

    raise typer.BadParameter("Provide --smiles, --file, or pipe input via stdin")


def _is_json_output(format_override: Optional[str]) -> bool:
    """Determine whether output should be JSON.

    Returns True when --format json is given, or when stdout is piped
    (non-interactive).
    """
    if format_override is not None:
        return format_override.lower() == "json"
    return not sys.stdout.isatty()


def _output_json(data: dict) -> None:
    """Print data as pretty-printed JSON."""
    typer.echo(json.dumps(data, indent=2, default=str))


def _output_validation_table(data: dict) -> None:
    """Print validation result as a Rich table."""
    table = Table(title="Validation Result")
    table.add_column("Field", style="cyan")
    table.add_column("Value", style="white")

    table.add_row("Status", str(data.get("status", "completed")))
    table.add_row("Overall Score", str(data.get("overall_score", "N/A")))
    mol_info = data.get("molecule_info", {})
    table.add_row("SMILES", str(mol_info.get("canonical_smiles", "N/A")))
    table.add_row("InChIKey", str(mol_info.get("inchikey", "N/A")))
    table.add_row("Issues", str(len(data.get("issues", []))))
    console.print(table)


def _output_scoring_table(data: dict) -> None:
    """Print scoring result as a Rich table."""
    table = Table(title="Scoring Result")
    table.add_column("Field", style="cyan")
    table.add_column("Value", style="white")

    # Extract druglikeness info if present
    dl = data.get("druglikeness", {})
    if dl:
        lipinski = dl.get("lipinski", {})
        table.add_row("Lipinski Passed", str(lipinski.get("passed", "N/A")))
        qed = dl.get("qed", {})
        table.add_row("QED Score", str(qed.get("score", "N/A")))

    # ML readiness
    ml = data.get("ml_readiness", {})
    if ml:
        table.add_row("ML Readiness", str(ml.get("score", "N/A")))

    # NP-likeness
    np_res = data.get("np_likeness", {})
    if np_res:
        table.add_row("NP-likeness", str(np_res.get("score", "N/A")))

    # ADMET
    admet = data.get("admet", {})
    if admet:
        sa = admet.get("synthetic_accessibility", {})
        table.add_row("SA Score", str(sa.get("score", "N/A")))

    console.print(table)


def _output_standardization_table(data: dict) -> None:
    """Print standardization result as a Rich table."""
    table = Table(title="Standardization Result")
    table.add_column("Field", style="cyan")
    table.add_column("Value", style="white")

    result = data.get("result", data)
    table.add_row("Original SMILES", str(result.get("original_smiles", "N/A")))
    table.add_row("Standardized SMILES", str(result.get("standardized_smiles", "N/A")))
    changes = result.get("steps_applied", result.get("changes", []))
    table.add_row("Changes", str(len(changes) if isinstance(changes, list) else changes))
    console.print(table)


def _output_profile_table(data: dict) -> None:
    """Print compound profile as a Rich table."""
    table = Table(title="Compound Profile")
    table.add_column("Property", style="cyan")
    table.add_column("Value", style="white")

    # If there's a properties dict, iterate over it
    properties = data.get("properties", data)
    if isinstance(properties, dict):
        for key, value in properties.items():
            if key in ("molecule_info", "execution_time_ms"):
                continue
            if isinstance(value, dict):
                # Flatten nested dicts one level
                for sub_key, sub_val in value.items():
                    table.add_row(f"{key}.{sub_key}", str(sub_val))
            else:
                table.add_row(key, str(value))
    console.print(table)


def _http_post(server: str, path: str, payload: dict) -> dict:
    """Send a POST request to the ChemAudit API server.

    Args:
        server: Base URL of the server (e.g., http://localhost:8000).
        path: API path (e.g., /api/v1/validate).
        payload: JSON-serializable request body.

    Returns:
        Parsed JSON response as a dict.
    """
    import httpx

    response = httpx.post(f"{server}{path}", json=payload, timeout=60)
    response.raise_for_status()
    return response.json()


def _local_validate(smiles: str) -> dict:
    """Validate a molecule using direct service imports (offline mode)."""
    from app.services.parser.molecule_parser import parse_molecule
    from app.services.validation.engine import validation_engine

    parsed = parse_molecule(smiles)
    if parsed is None or not parsed.success or parsed.mol is None:
        return {"status": "error", "message": "Invalid molecule"}

    results, score = validation_engine.validate(parsed.mol)
    issues = [
        {
            "check_name": r.check_name,
            "passed": r.passed,
            "severity": r.severity.value if hasattr(r.severity, "value") else str(r.severity),
            "message": r.message,
        }
        for r in results
        if not r.passed
    ]
    return {
        "status": "completed",
        "overall_score": score,
        "molecule_info": {
            "canonical_smiles": smiles,
        },
        "issues": issues,
    }


def _local_score(smiles: str) -> dict:
    """Score a molecule using direct service imports (offline mode)."""
    from rdkit import Chem

    from app.services.scoring import (
        calculate_admet,
        calculate_druglikeness,
        calculate_ml_readiness,
        calculate_np_likeness,
    )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "error", "message": "Invalid molecule"}

    ml = calculate_ml_readiness(mol)
    dl = calculate_druglikeness(mol, include_extended=True)
    np_res = calculate_np_likeness(mol)
    admet = calculate_admet(mol)

    return {
        "ml_readiness": {"score": ml.score, "label": ml.label},
        "druglikeness": {
            "lipinski": {"passed": dl.lipinski.passed, "violations": dl.lipinski.violations},
            "qed": {"score": dl.qed.score},
        },
        "np_likeness": {"score": np_res.score, "interpretation": np_res.interpretation},
        "admet": {
            "synthetic_accessibility": {
                "score": admet.synthetic_accessibility.score,
                "classification": admet.synthetic_accessibility.classification,
            },
        },
    }


def _local_standardize(smiles: str) -> dict:
    """Standardize a molecule using direct service imports (offline mode)."""
    from rdkit import Chem

    from app.services.standardization.chembl_pipeline import standardize_molecule

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "error", "message": "Invalid molecule"}

    result = standardize_molecule(mol)
    return {
        "result": {
            "original_smiles": result.original_smiles,
            "standardized_smiles": result.standardized_smiles,
            "success": result.success,
            "steps_applied": [
                {"step_name": s.step_name, "applied": s.applied, "description": s.description}
                for s in result.steps_applied
            ],
        }
    }


def _local_profile(smiles: str) -> dict:
    """Profile a molecule using direct service imports (offline mode)."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "error", "message": "Invalid molecule"}

    return {
        "properties": {
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "num_aromatic_rings": Descriptors.NumAromaticRings(mol),
            "fsp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
            "num_atoms": mol.GetNumAtoms(),
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        }
    }


@app.command()
def validate(
    smiles: Optional[str] = typer.Option(None, "--smiles", "-s", help="SMILES string to validate"),
    file: Optional[str] = typer.Option(None, "--file", "-f", help="Input file (CSV/SDF)"),
    server: str = typer.Option(
        "http://localhost:8000", "--server", help="API server URL"
    ),
    local: bool = typer.Option(False, "--local", help="Use direct service import (no server)"),
    format: Optional[str] = typer.Option(None, "--format", help="Output format: json|table"),
) -> None:
    """Validate a chemical structure."""
    molecule = _get_smiles(smiles, file)
    try:
        if local:
            data = _local_validate(molecule)
        else:
            data = _http_post(server, "/api/v1/validate", {"molecule": molecule})
        if _is_json_output(format):
            _output_json(data)
        else:
            _output_validation_table(data)
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(code=1)


@app.command()
def score(
    smiles: Optional[str] = typer.Option(None, "--smiles", "-s", help="SMILES string to score"),
    file: Optional[str] = typer.Option(None, "--file", "-f", help="Input file (CSV/SDF)"),
    server: str = typer.Option(
        "http://localhost:8000", "--server", help="API server URL"
    ),
    local: bool = typer.Option(False, "--local", help="Use direct service import (no server)"),
    format: Optional[str] = typer.Option(None, "--format", help="Output format: json|table"),
) -> None:
    """Score a chemical structure for drug-likeness."""
    molecule = _get_smiles(smiles, file)
    try:
        if local:
            data = _local_score(molecule)
        else:
            data = _http_post(server, "/api/v1/score", {"molecule": molecule})
        if _is_json_output(format):
            _output_json(data)
        else:
            _output_scoring_table(data)
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(code=1)


@app.command()
def standardize(
    smiles: Optional[str] = typer.Option(
        None, "--smiles", "-s", help="SMILES string to standardize"
    ),
    file: Optional[str] = typer.Option(None, "--file", "-f", help="Input file (CSV/SDF)"),
    server: str = typer.Option(
        "http://localhost:8000", "--server", help="API server URL"
    ),
    local: bool = typer.Option(False, "--local", help="Use direct service import (no server)"),
    format: Optional[str] = typer.Option(None, "--format", help="Output format: json|table"),
) -> None:
    """Standardize a chemical structure."""
    molecule = _get_smiles(smiles, file)
    try:
        if local:
            data = _local_standardize(molecule)
        else:
            data = _http_post(server, "/api/v1/standardize", {"molecule": molecule})
        if _is_json_output(format):
            _output_json(data)
        else:
            _output_standardization_table(data)
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(code=1)


@app.command()
def profile(
    smiles: Optional[str] = typer.Option(
        None, "--smiles", "-s", help="SMILES string to profile"
    ),
    file: Optional[str] = typer.Option(None, "--file", "-f", help="Input file (CSV/SDF)"),
    server: str = typer.Option(
        "http://localhost:8000", "--server", help="API server URL"
    ),
    local: bool = typer.Option(False, "--local", help="Use direct service import (no server)"),
    format: Optional[str] = typer.Option(None, "--format", help="Output format: json|table"),
) -> None:
    """Profile a chemical structure (properties, drug-likeness, etc.)."""
    molecule = _get_smiles(smiles, file)
    try:
        if local:
            data = _local_profile(molecule)
        else:
            data = _http_post(
                server,
                "/api/v1/score",
                {"molecule": molecule, "include": ["druglikeness", "admet", "np_likeness"]},
            )
        if _is_json_output(format):
            _output_json(data)
        else:
            _output_profile_table(data)
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(code=1)


def main():
    """CLI entry point."""
    app()
