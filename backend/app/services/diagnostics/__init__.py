"""
Diagnostics package for structure quality diagnostics.

Provides five diagnostic tools:
- SMILES error diagnostics with position extraction and fix suggestions (DIAG-01)
- InChI layer-by-layer comparison (DIAG-02)
- Format round-trip lossiness detection (DIAG-03)
- Cross-pipeline standardization comparison (DIAG-04)
- SDF/CSV file pre-validation (DIAG-05)
"""

from .file_prevalidator import prevalidate_csv, prevalidate_sdf
from .format_roundtrip import check_roundtrip
from .inchi_diff import diff_inchi_layers, parse_inchi_layers
from .smiles_diagnostics import diagnose_smiles
from .std_comparison import compare_pipelines

__all__ = [
    "diagnose_smiles",
    "diff_inchi_layers",
    "parse_inchi_layers",
    "check_roundtrip",
    "compare_pipelines",
    "prevalidate_sdf",
    "prevalidate_csv",
]
