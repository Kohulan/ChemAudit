"""
File Parser Module

Parses SDF and CSV files into molecule data for batch processing.
Handles errors per-molecule without crashing the entire batch.
"""
import io
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any

import pandas as pd
from rdkit import Chem


@dataclass
class MoleculeData:
    """Data extracted from a molecule in a batch file."""

    smiles: str
    name: Optional[str] = None
    index: int = 0
    properties: Dict[str, Any] = field(default_factory=dict)
    parse_error: Optional[str] = None


def parse_sdf(file_content: bytes) -> List[MoleculeData]:
    """
    Parse an SDF file into a list of molecule data.

    Args:
        file_content: Raw bytes of the SDF file

    Returns:
        List of MoleculeData objects (includes entries with parse_error for invalid molecules)

    Note:
        Individual parse errors are captured per-molecule rather than failing the entire batch.
    """
    molecules = []
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(file_content))

    for idx, mol in enumerate(supplier):
        if mol is None:
            # Failed to parse this molecule
            molecules.append(
                MoleculeData(
                    smiles="",
                    name=f"mol_{idx}",
                    index=idx,
                    parse_error=f"Failed to parse molecule at index {idx}",
                )
            )
            continue

        try:
            # Get SMILES
            smiles = Chem.MolToSmiles(mol)

            # Extract name from _Name property or generate index-based
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{idx}"
            if not name.strip():
                name = f"mol_{idx}"

            # Extract all other properties
            properties = {}
            for prop_name in mol.GetPropsAsDict():
                if prop_name != "_Name":
                    try:
                        properties[prop_name] = mol.GetProp(prop_name)
                    except Exception:
                        pass

            molecules.append(
                MoleculeData(
                    smiles=smiles,
                    name=name,
                    index=idx,
                    properties=properties,
                )
            )
        except Exception as e:
            molecules.append(
                MoleculeData(
                    smiles="",
                    name=f"mol_{idx}",
                    index=idx,
                    parse_error=f"Error extracting data: {str(e)}",
                )
            )

    return molecules


def parse_csv(
    file_content: bytes,
    smiles_column: str = "SMILES",
    name_column: Optional[str] = None,
) -> List[MoleculeData]:
    """
    Parse a CSV file into a list of molecule data.

    Args:
        file_content: Raw bytes of the CSV file
        smiles_column: Name of the column containing SMILES strings
        name_column: Optional name of the column containing molecule names

    Returns:
        List of MoleculeData objects

    Raises:
        ValueError: If the SMILES column is not found in the CSV

    Note:
        Individual parse errors are captured per-molecule.
    """
    try:
        df = pd.read_csv(io.BytesIO(file_content))
    except Exception as e:
        raise ValueError(f"Failed to parse CSV file: {str(e)}")

    # Check for SMILES column (case-insensitive search)
    smiles_col_actual = None
    for col in df.columns:
        if col.lower() == smiles_column.lower():
            smiles_col_actual = col
            break

    if smiles_col_actual is None:
        available = ", ".join(df.columns.tolist())
        raise ValueError(
            f"SMILES column '{smiles_column}' not found. Available columns: {available}"
        )

    # Determine name column
    name_col_actual = None
    if name_column:
        for col in df.columns:
            if col.lower() == name_column.lower():
                name_col_actual = col
                break
    else:
        # Try to auto-detect name column
        for candidate in ["Name", "name", "NAME", "ID", "id", "Compound", "compound"]:
            if candidate in df.columns:
                name_col_actual = candidate
                break

    molecules = []
    for idx, row in df.iterrows():
        smiles = str(row[smiles_col_actual]).strip()

        # Validate SMILES is not empty or NaN
        if not smiles or smiles.lower() == "nan":
            molecules.append(
                MoleculeData(
                    smiles="",
                    name=f"row_{idx}",
                    index=idx,
                    parse_error="Empty or missing SMILES",
                )
            )
            continue

        # Get name
        if name_col_actual and pd.notna(row[name_col_actual]):
            name = str(row[name_col_actual]).strip()
        else:
            name = f"row_{idx}"

        # Collect other columns as properties
        properties = {}
        for col in df.columns:
            if col not in [smiles_col_actual, name_col_actual]:
                val = row[col]
                if pd.notna(val):
                    properties[col] = val

        molecules.append(
            MoleculeData(
                smiles=smiles,
                name=name,
                index=idx,
                properties=properties,
            )
        )

    return molecules


def detect_csv_columns(file_content: bytes) -> Dict[str, List[str]]:
    """
    Detect column names in a CSV file for user selection.

    Args:
        file_content: Raw bytes of the CSV file

    Returns:
        Dictionary with 'columns' list and 'suggested_smiles' (best guess for SMILES column)
    """
    try:
        df = pd.read_csv(io.BytesIO(file_content), nrows=5)
        columns = df.columns.tolist()

        # Try to detect SMILES column
        suggested_smiles = None
        for col in columns:
            col_lower = col.lower()
            if "smiles" in col_lower or col_lower == "smi":
                suggested_smiles = col
                break

        return {
            "columns": columns,
            "suggested_smiles": suggested_smiles,
            "row_count_estimate": len(pd.read_csv(io.BytesIO(file_content))),
        }
    except Exception as e:
        raise ValueError(f"Failed to read CSV: {str(e)}")
