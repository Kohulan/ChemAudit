"""SA Comparison service: SA Score + SCScore + SYBA side-by-side.

Computes three independent synthesizability scores and returns them together
for comparison. Handles all graceful-fallback scenarios:
  - SA Score: Always available via RDKit Contrib sascorer (bundled with RDKit)
  - SCScore: Optional pip install; gracefully falls back on ImportError or
    weight file NumPy 2.x incompatibility
  - SYBA: Requires subprocess isolation; gracefully falls back on timeout or
    subprocess failure

LICENSING NOTE: SYBA is GPL-3.0. It MUST ONLY be called via subprocess
isolation. Never import syba at module level in this Apache-2.0 codebase.
"""

import logging
import os
import subprocess
import sys
from typing import Any, Dict, Optional

from rdkit import Chem, RDConfig

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore  # noqa: E402

logger = logging.getLogger(__name__)

# Module-level cache for SCScore scorer (loaded once, reused across calls)
_SCSCORE_SCORER: Optional[Any] = None
_SCSCORE_LOAD_ATTEMPTED: bool = False


# ---------------------------------------------------------------------------
# SA Score (RDKit Contrib — always available)
# ---------------------------------------------------------------------------


def _compute_sa_score(mol: Chem.Mol) -> Dict[str, Any]:
    """Compute SA Score via RDKit Contrib sascorer.

    SA Score ranges from 1 (easy) to 10 (difficult to synthesize).
    This always succeeds for valid RDKit molecules — sascorer is bundled
    with every RDKit installation.

    Args:
        mol: RDKit molecule object (must be valid)

    Returns:
        dict with keys: score (float), scale (str), classification (str), available (bool)
    """
    score = sascorer.calculateScore(mol)
    if score < 3:
        classification = "easy"
    elif score < 5:
        classification = "moderate"
    else:
        classification = "difficult"

    return {
        "score": round(score, 2),
        "scale": "1-10",
        "classification": classification,
        "available": True,
    }


# ---------------------------------------------------------------------------
# SCScore (optional pip install — graceful fallback)
# ---------------------------------------------------------------------------


def _load_scscore() -> Optional[Any]:
    """Attempt to load and restore the SCScore model.

    Tries two paths:
    1. Standard SCScorer().restore() — works when scscore is pip-installed
       and NumPy compat is OK
    2. Vendor path fallback — looks for scscore_weights.npz in the vendor/
       subdirectory alongside this module (pre-converted for NumPy 2.x)

    Returns:
        Loaded SCScorer instance, or None if unavailable
    """
    global _SCSCORE_SCORER, _SCSCORE_LOAD_ATTEMPTED

    if _SCSCORE_LOAD_ATTEMPTED:
        return _SCSCORE_SCORER

    _SCSCORE_LOAD_ATTEMPTED = True

    try:
        from scscore.standalone_model_numpy import SCScorer  # type: ignore

        scorer = SCScorer()
        try:
            scorer.restore()
            _SCSCORE_SCORER = scorer
            return _SCSCORE_SCORER
        except Exception as restore_err:
            logger.debug("SCScore restore() failed (%s), trying vendor fallback", restore_err)
            # Try vendor-supplied .npz weight file (pre-converted for NumPy 2.x)
            vendor_path = os.path.join(os.path.dirname(__file__), "vendor", "scscore_weights.npz")
            if os.path.isfile(vendor_path):
                try:
                    scorer2 = SCScorer()
                    scorer2.restore(weight_path=vendor_path)
                    _SCSCORE_SCORER = scorer2
                    return _SCSCORE_SCORER
                except Exception as vendor_err:
                    logger.debug("SCScore vendor fallback failed: %s", vendor_err)

    except ImportError:
        logger.debug("scscore not installed — SCScore will be unavailable")
    except Exception as err:
        logger.debug("Unexpected error loading SCScore: %s", err)

    _SCSCORE_SCORER = None
    return None


def _compute_scscore(smiles: str) -> Dict[str, Any]:
    """Compute SCScore for a SMILES string.

    SCScore ranges from 1 (purchasable starting material) to 5
    (requires many complex synthesis steps).

    Args:
        smiles: Canonical SMILES string

    Returns:
        dict with keys: score/error, scale, classification, available
    """
    scorer = _load_scscore()
    if scorer is None:
        return {
            "error": "scscore not available",
            "available": False,
        }

    try:
        _, score = scorer.apply(smiles)
        if score < 2:
            classification = "easy"
        elif score < 3:
            classification = "moderate"
        else:
            classification = "difficult"

        return {
            "score": round(score, 2),
            "scale": "1-5",
            "classification": classification,
            "available": True,
        }
    except Exception as err:
        logger.debug("SCScore apply() failed for SMILES %r: %s", smiles, err)
        return {
            "error": f"SCScore computation failed: {err}",
            "available": False,
        }


# ---------------------------------------------------------------------------
# SYBA (GPL-3.0 — subprocess isolation only, never imported at module level)
# ---------------------------------------------------------------------------


def _syba_via_subprocess(smiles: str) -> Optional[float]:
    """Run SYBA prediction via subprocess to maintain GPL-3.0 isolation.

    SYBA (SYnthesizaBility Assessment) is GPL-3.0 licensed. This Apache-2.0
    codebase must never import it directly. Instead, we spawn a fresh Python
    subprocess that imports and runs SYBA, then captures the numeric output.

    Args:
        smiles: SMILES string to score

    Returns:
        SYBA score as float, or None on any failure (timeout, import error,
        subprocess error, parse error)
    """
    script = (
        "from syba.syba import SybaClassifier; "
        "s = SybaClassifier(); "
        "s.fitDefaultScore(); "
        f"score = s.predict({smiles!r}); "
        "print(score)"
    )
    try:
        proc = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if proc.returncode != 0:
            logger.debug("SYBA subprocess failed (rc=%d): %s", proc.returncode, proc.stderr[:200])
            return None

        return float(proc.stdout.strip())

    except subprocess.TimeoutExpired:
        logger.debug("SYBA subprocess timed out after 30s")
        return None
    except (ValueError, OSError) as err:
        logger.debug("SYBA subprocess error: %s", err)
        return None
    except Exception as err:
        logger.debug("Unexpected SYBA subprocess error: %s", err)
        return None


def _compute_syba(smiles: str) -> Dict[str, Any]:
    """Compute SYBA score for a SMILES string via subprocess isolation.

    SYBA score ranges from approximately -166 (very hard) to +275 (very easy).
    Positive scores indicate synthesizable, negative indicate difficult.

    Args:
        smiles: Canonical SMILES string

    Returns:
        dict with keys: score/error, scale, classification, available
    """
    score = _syba_via_subprocess(smiles)
    if score is None:
        return {
            "error": "SYBA unavailable",
            "available": False,
        }

    if score > 50:
        classification = "easy"
    elif score > 0:
        classification = "moderate"
    else:
        classification = "difficult"

    return {
        "score": round(score, 1),
        "scale": "-166 to +275",
        "classification": classification,
        "available": True,
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_sa_comparison(mol: Chem.Mol, smiles: str) -> Dict[str, Any]:
    """Compute all three synthesizability scores side-by-side.

    Returns SA Score (always available), SCScore (optional dep), SYBA
    (subprocess isolation), and a RAscore placeholder slot (out of scope
    for v3.0 — requires Python 3.7 + TensorFlow 2.5).

    Args:
        mol: RDKit molecule object (must be valid; used for SA Score)
        smiles: Canonical SMILES string (used for SCScore and SYBA)

    Returns:
        dict with keys: sa_score, scscore, syba, rascore
        Each value is a dict with at minimum "available" (bool) and either
        "score" (float) or "error" (str).
    """
    return {
        "sa_score": _compute_sa_score(mol),
        "scscore": _compute_scscore(smiles),
        "syba": _compute_syba(smiles),
        "rascore": {
            "available": False,
            "note": "Coming in a future release",
        },
    }
