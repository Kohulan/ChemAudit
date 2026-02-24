"""
IUPAC Name to SMILES Converter

Uses JPype + OPSIN JAR for converting IUPAC chemical names to SMILES notation.
Includes a heuristic input type detector to distinguish SMILES from IUPAC names.
"""

import logging
import re
from pathlib import Path

from app.core.config import settings

logger = logging.getLogger(__name__)

# Module-level singleton for OPSIN NameToStructure
_nts = None


def init_opsin() -> None:
    """Initialize the OPSIN converter via JPype.

    Starts the JVM (if not already running) and loads the OPSIN NameToStructure
    singleton. This should be called once during application startup.

    Raises:
        FileNotFoundError: If OPSIN JAR not found at configured path
        ImportError: If jpype1 is not installed
        RuntimeError: If JVM fails to start
    """
    global _nts

    if _nts is not None:
        return  # Already initialized

    jar_path = Path(settings.OPSIN_JAR_PATH)
    if not jar_path.exists():
        raise FileNotFoundError(f"OPSIN JAR not found at {jar_path}")

    import jpype

    if not jpype.isJVMStarted():
        jpype.startJVM(
            jpype.getDefaultJVMPath(),
            "-Djava.awt.headless=true",
            f"-Djava.class.path={jar_path}",
        )

    nts_class = jpype.JClass("uk.ac.cam.ch.wwmm.opsin.NameToStructure")
    _nts = nts_class.getInstance()
    logger.info("OPSIN initialized successfully from %s", jar_path)


def iupac_to_smiles(name: str) -> str | None:
    """Convert an IUPAC chemical name to SMILES.

    Args:
        name: IUPAC chemical name (e.g., 'aspirin', '2-methylpropan-1-ol')

    Returns:
        SMILES string if conversion succeeds, None otherwise

    Raises:
        RuntimeError: If OPSIN has not been initialized
    """
    if _nts is None:
        raise RuntimeError("OPSIN not initialized. Call init_opsin() first.")

    try:
        result = _nts.parseChemicalName(name)
        smiles = str(result.getSmiles())
        return smiles if smiles and smiles != "null" else None
    except Exception:
        logger.debug("OPSIN conversion failed for '%s'", name, exc_info=True)
        return None


def is_opsin_available() -> bool:
    """Check if OPSIN is initialized and available.

    Returns:
        True if OPSIN is ready for use
    """
    return _nts is not None


# SMILES-specific characters that distinguish SMILES from text
_SMILES_CHARS = re.compile(r"[()[\]=#@/\\]")

# IUPAC suffixes and patterns
_IUPAC_SUFFIXES = (
    "-yl", "-ane", "-ol", "-ene", "-yne", "acid", "-one", "-al",
    "-amine", "-amide", "-ester", "-ether", "-oxy", "-oic",
    "-ic", "-ous", "-ate", "-ide", "-ite",
)


def detect_input_type(value: str) -> str:
    """Detect whether input is likely SMILES or an IUPAC name.

    Uses heuristics based on character patterns:
    - SMILES typically contains: (), [], =, #, @, /, \\
    - IUPAC names contain: spaces, hyphens, chemical suffixes, mostly alphabetic

    Args:
        value: Input string to classify

    Returns:
        One of 'smiles', 'iupac', or 'ambiguous'
    """
    stripped = value.strip()
    if not stripped:
        return "ambiguous"

    # Check for SMILES-specific characters
    if _SMILES_CHARS.search(stripped):
        return "smiles"

    # Check for IUPAC indicators
    lower = stripped.lower()

    # Spaces strongly indicate IUPAC (SMILES never has spaces)
    if " " in stripped:
        return "iupac"

    # IUPAC suffixes
    for suffix in _IUPAC_SUFFIXES:
        if lower.endswith(suffix):
            return "iupac"

    # Purely alphabetic with hyphens (common IUPAC pattern)
    if re.match(r"^[a-zA-Z][a-zA-Z0-9,\-]*$", stripped) and "-" in stripped:
        return "iupac"

    # Short purely alphabetic strings like "aspirin", "caffeine" are likely IUPAC
    if re.match(r"^[a-zA-Z]{4,}$", stripped):
        return "iupac"

    # Could be simple SMILES like "CCO" or abbreviation
    return "ambiguous"
