"""SA Comparison service: SA Score + SCScore + SYBA side-by-side.

Computes three independent synthesizability scores and returns them together
for comparison. Handles all graceful-fallback scenarios:
  - SA Score: Always available via RDKit Contrib sascorer (bundled with RDKit)
  - SCScore: Self-contained reimplementation over vendored weights (no external
    package); gracefully falls back if the weight file is absent
  - SYBA: Requires subprocess isolation; gracefully falls back on timeout or
    subprocess failure

LICENSING NOTE: SYBA is GPL-3.0. It MUST ONLY be called via subprocess
isolation. Never import syba at module level in this Apache-2.0 codebase.
"""

import json
import logging
import os
import queue
import subprocess
import sys
import threading
from functools import lru_cache
from typing import Any, Dict, Optional

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import rdFingerprintGenerator

from app.core.error_sanitizer import safe_error_detail

logger = logging.getLogger(__name__)


@lru_cache(maxsize=1)
def _get_sascorer():
    """Lazily import RDKit Contrib sascorer; return None if unavailable."""
    try:
        sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_path not in sys.path:
            sys.path.insert(0, sa_path)
        import sascorer  # type: ignore[import-untyped]

        return sascorer
    except Exception:
        return None


# Module-level cache for the vendored SCScore weights (loaded once, reused).
_SCSCORE_WEIGHTS: Optional[list] = None
_SCSCORE_LOAD_ATTEMPTED: bool = False
_SCSCORE_WEIGHTS_PATH = os.path.join(os.path.dirname(__file__), "vendor", "scscore_weights.npz")

# SCScore model architecture (Coley et al., 2018) — must match the trained
# weights. The model is a small MLP over a Morgan-2 1024-bit fingerprint;
# scores range 1 (trivial starting material) to 5 (many complex steps).
_SCSCORE_FP_LEN = 1024
_SCSCORE_FP_RAD = 2
_SCSCORE_SCALE = 5.0


@lru_cache(maxsize=1)
def _get_scscore_fp_gen():
    """Reusable Morgan generator matching SCScore's training featurisation."""
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=_SCSCORE_FP_RAD, fpSize=_SCSCORE_FP_LEN, includeChirality=True
    )


# Upper bound on SMILES length accepted by the SYBA subprocess path. Scoring a
# SMILES longer than this is pointless and only widens the input surface.
_MAX_SYBA_SMILES_LEN = 2000


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
    _sa = _get_sascorer()
    if _sa is None:
        return {
            "score": None,
            "scale": "1-10",
            "classification": None,
            "available": False,
            "error": "sascorer not available",
        }

    score = _sa.calculateScore(mol)
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
# SCScore (self-contained reimplementation over vendored weights)
# ---------------------------------------------------------------------------


def _load_scscore_weights() -> Optional[list]:
    """Load the vendored SCScore MLP weights once (list of 12 numpy arrays).

    The weights are the upstream ``model.ckpt-10654`` (full_reaxys_model_1024bool,
    Coley et al.) re-saved from the project-agnostic ``json.gz`` form into a
    portable, pickle-free ``.npz`` (plain ``.npy`` arrays — readable across NumPy
    versions). Returns None if the file is absent, so SCScore degrades
    gracefully without ever raising.

    Returns:
        Ordered list of weight/bias arrays [W0, b0, W1, b1, ...], or None.
    """
    global _SCSCORE_WEIGHTS, _SCSCORE_LOAD_ATTEMPTED

    if _SCSCORE_LOAD_ATTEMPTED:
        return _SCSCORE_WEIGHTS

    _SCSCORE_LOAD_ATTEMPTED = True

    if not os.path.isfile(_SCSCORE_WEIGHTS_PATH):
        logger.debug(
            "SCScore weights not present at %s — SCScore unavailable (see vendor/README.md)",
            _SCSCORE_WEIGHTS_PATH,
        )
        _SCSCORE_WEIGHTS = None
        return None

    try:
        with np.load(_SCSCORE_WEIGHTS_PATH) as data:
            keys = sorted(data.files, key=lambda s: int(s.split("_")[1]))
            _SCSCORE_WEIGHTS = [data[k] for k in keys]
    except Exception as err:
        logger.warning("Failed to load vendored SCScore weights (%s); SCScore unavailable", err)
        _SCSCORE_WEIGHTS = None

    return _SCSCORE_WEIGHTS


def _scscore_fingerprint(mol: Chem.Mol) -> np.ndarray:
    """Morgan fingerprint (radius 2, 1024 bits, chirality) as a float64 vector.

    Reproduces the exact featurisation the upstream SCScore model was trained
    with. Bit values are 0/1, exactly representable, so float64 here matches the
    model's original float32 featurisation bit-for-bit (verified against the
    upstream standalone model to < 1e-6).
    """
    return _get_scscore_fp_gen().GetFingerprintAsNumPy(mol).astype(np.float64)


def _scscore_forward(weights: list, fp: np.ndarray) -> float:
    """Run the SCScore MLP forward pass: dense layers + ReLU, sigmoid-scaled to 1-5.

    Mirrors upstream ``SCScorer.apply``: weights are (W, b) pairs applied in
    sequence with ReLU between hidden layers, then ``1 + (scale-1)*sigmoid(x)``.
    """
    x = fp
    n = len(weights)
    for i in range(0, n, 2):
        x = np.matmul(x, weights[i]) + weights[i + 1]
        if i != n - 2:  # ReLU on every layer except the last
            x = x * (x > 0)
    x = 1.0 + (_SCSCORE_SCALE - 1.0) * (1.0 / (1.0 + np.exp(-x)))
    return float(np.asarray(x).ravel()[0])


def _compute_scscore(smiles: str) -> Dict[str, Any]:
    """Compute SCScore for a SMILES string.

    SCScore ranges from 1 (purchasable starting material) to 5
    (requires many complex synthesis steps).

    Args:
        smiles: Canonical SMILES string

    Returns:
        dict with keys: score/error, scale, classification, available
    """
    weights = _load_scscore_weights()
    if weights is None:
        return {
            "error": "scscore not available",
            "available": False,
        }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                "error": "Invalid SMILES for SCScore",
                "available": False,
            }

        fp = _scscore_fingerprint(mol)
        # Upstream returns 0 when the fingerprint is empty (no bits set).
        score = 0.0 if float(fp.sum()) == 0.0 else _scscore_forward(weights, fp)

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
        return {
            "error": safe_error_detail(err, "SCScore computation failed"),
            "available": False,
        }


# ---------------------------------------------------------------------------
# SYBA (GPL-3.0 — subprocess isolation only, never imported at module level)
# ---------------------------------------------------------------------------


_SYBA_WORKER_PATH = os.path.join(os.path.dirname(__file__), "syba_worker.py")
_SYBA_READY_TIMEOUT = 60.0  # seconds to allow for one-time SYBA model load
_SYBA_PREDICT_TIMEOUT = 30.0  # seconds per prediction


class _SybaWorker:
    """Manages a long-lived SYBA subprocess to amortise the model-load cost.

    SYBA (GPL-3.0) is confined to a separate process (``syba_worker.py``); this
    Apache-2.0 module only exchanges JSON lines with it over stdin/stdout. The
    worker is started lazily on first use and reused for all later predictions.
    Access is serialised by a lock (a single worker, one request at a time), a
    background reader thread feeds responses to a queue for timeout handling, and
    a crashed/hung worker is killed and transparently restarted on the next call.
    A permanent failure (SYBA not importable) is sticky so we don't respawn.
    """

    def __init__(self, command: Optional[list] = None):
        self._command = command or [sys.executable, "-u", _SYBA_WORKER_PATH]
        self._proc: Optional[subprocess.Popen] = None
        self._queue: "queue.Queue[Optional[str]]" = queue.Queue()
        self._lock = threading.Lock()
        self._fatal = False  # SYBA permanently unavailable (e.g. not installed)

    def _read_loop(self, proc: subprocess.Popen) -> None:
        try:
            for line in proc.stdout:  # type: ignore[union-attr]
                self._queue.put(line)
        finally:
            self._queue.put(None)  # EOF sentinel

    def _kill(self) -> None:
        if self._proc is not None:
            try:
                self._proc.kill()
            except Exception:
                pass
        self._proc = None

    def _ensure_started(self) -> bool:
        if self._fatal:
            return False
        if self._proc is not None and self._proc.poll() is None:
            return True

        try:
            self._proc = subprocess.Popen(
                self._command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                bufsize=1,
                shell=False,
            )
        except Exception as exc:
            logger.warning("SYBA worker failed to start: %s", exc)
            self._fatal = True
            return False

        self._queue = queue.Queue()
        threading.Thread(target=self._read_loop, args=(self._proc,), daemon=True).start()

        # Handshake: wait for the worker to load the model (or report fatal).
        try:
            line = self._queue.get(timeout=_SYBA_READY_TIMEOUT)
        except queue.Empty:
            logger.warning("SYBA worker did not become ready within %ss", _SYBA_READY_TIMEOUT)
            self._kill()
            return False
        try:
            msg = json.loads(line) if line else {}
        except Exception:
            msg = {}
        if msg.get("ready"):
            return True
        # Not ready. A fatal handshake (SYBA missing) is sticky; otherwise retry later.
        if msg.get("fatal"):
            logger.warning("SYBA unavailable: %s", msg.get("error", "model load failed"))
            self._fatal = True
        self._kill()
        return False

    def predict(self, smiles: str) -> Optional[float]:
        """Return the SYBA score for ``smiles``, or None if unavailable."""
        # Defense-in-depth: reject pathological input before touching the worker.
        if not smiles or len(smiles) > _MAX_SYBA_SMILES_LEN or not smiles.isprintable():
            logger.debug("SYBA input rejected (length/non-printable): %r", smiles[:60])
            return None

        with self._lock:
            if not self._ensure_started():
                return None
            try:
                self._proc.stdin.write(json.dumps({"smiles": smiles}) + "\n")  # type: ignore[union-attr]
                self._proc.stdin.flush()  # type: ignore[union-attr]
            except Exception as exc:
                logger.debug("SYBA worker write failed: %s", exc)
                self._kill()
                return None

            try:
                line = self._queue.get(timeout=_SYBA_PREDICT_TIMEOUT)
            except queue.Empty:
                logger.debug("SYBA prediction timed out; restarting worker")
                self._kill()
                return None

            if line is None:  # worker died mid-request
                self._kill()
                return None

            try:
                msg = json.loads(line)
            except Exception:
                return None
            score = msg.get("score")
            return float(score) if score is not None else None

    def shutdown(self) -> None:
        with self._lock:
            if self._proc is not None:
                try:
                    self._proc.stdin.close()  # type: ignore[union-attr]
                except Exception:
                    pass
                self._kill()


# Module-level singleton reused across all SYBA scoring calls.
_syba_worker = _SybaWorker()


def _syba_via_subprocess(smiles: str) -> Optional[float]:
    """Compute a SYBA score via the persistent worker (GPL-3.0 isolated)."""
    return _syba_worker.predict(smiles)


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
