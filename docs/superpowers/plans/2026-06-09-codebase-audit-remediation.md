# Codebase Audit Remediation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remediate the concerns surfaced by the codebase audit (`.planning/codebase/CONCERNS.md`) in risk-ordered waves, fixing the high-confidence security/reliability/performance issues first with full test coverage, and staging the larger refactors behind explicit checkpoints.

**Architecture:** Work proceeds on branch `fix/codebase-audit-concerns` off `development`. Each fix is an atomic commit with a failing-test-first (TDD) loop where behavior changes. Wave 1 (this session) is surgical and reviewable; Waves 2–4 are larger and gated behind user checkpoints. Two audit findings are reclassified as **design decisions** and deliberately NOT auto-fixed (see "Flagged for Decision").

**Tech Stack:** Python 3.11 / FastAPI / Celery / Redis / RDKit (backend); React 18 / TypeScript / Vite / Vitest (frontend). Backend tests run with the `cheminformatics` conda env: `/Users/kohulan/anaconda3/envs/cheminformatics/bin/python -m pytest` from `backend/`.

**Baseline (2026-06-09):** `2262 passed, 15 failed, 34 skipped`. All 15 failures are pre-existing **environment** failures (no local Redis → `ConnectionRefusedError` for batch/export endpoints; tSNE/openTSNE module issues). Regression bar for every wave: **2262 unit tests stay green; total failures stay ≤ 15**.

---

## Flagged for Decision (NOT auto-fixed)

These two audit findings are **deliberate, tested design decisions**, not defects. Changing them would regress availability and/or contradict sibling code. They are surfaced for product/security sign-off rather than silently "fixed."

1. **WebSocket ownership "fail-open" when Redis unavailable** (`backend/app/main.py:65-66`, `_check_ws_ownership`).
   - The audit (CONCERNS.md "CSRF Validation Bypass on WebSocket Ownership Check") recommends returning `False` (deny) when `manager._redis` is `None`.
   - **However**, the sibling HTTP batch-ownership path `app/core/ownership.verify_job_access` makes the *same* fail-open choice **by design**, and `tests/test_ownership.py::TestVerifyJobAccess::test_graceful_redis_failure` (lines 107-114) explicitly asserts "Should allow access when Redis is unavailable." The mitigation is the v4-UUID job-id (high entropy, regex-gated at `main.py:446`).
   - Flipping only the WebSocket side would (a) deny *all* live-progress connections during any Redis blip — an availability regression for every user — to defend against an attacker who already knows a UUID, and (b) make WS and HTTP ownership inconsistent.
   - **DECISION (2026-06-09): keep fail-open** — confirmed by maintainer. Consistent with the tested HTTP path; availability preserved; UUID entropy is the mitigation. No code change.

2. **`is_ip_banned` fail-open on Redis error** (`backend/app/core/rate_limit.py:47-49`).
   - Returns `False` (not banned) on Redis error, with an explicit `# Fail open for availability` comment. This is intentional: a Redis outage should not lock every user out. The genuine sub-issue (connections never closed) IS fixed in Wave 1 Task 4 via pooling. The fail-open policy itself is left as-is.

---

## Wave 1 — Surgical security, reliability & performance (THIS SESSION)

High-confidence, low-blast-radius fixes. Each is a failing-test-first TDD loop and an atomic commit. Run command prefix (from `backend/`):

```
PY=/Users/kohulan/anaconda3/envs/cheminformatics/bin/python
$PY -m pytest <targets> -q -p no:cacheprovider
```

---

### Task 1: Celery task time limits (worker-starvation backstop)

Pathological molecules (deeply-fused rings, huge macrocycles) can make RDKit run effectively forever, permanently blocking a Celery worker. Add hard + soft global time limits.

**Files:**
- Modify: `backend/app/celery_app.py:34-83` (inside `celery_app.conf.update(...)`)
- Test: `backend/tests/test_celery_config.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# backend/tests/test_celery_config.py
"""Tests for Celery application configuration hardening."""

from app.celery_app import celery_app


def test_task_time_limits_are_configured():
    """Hard and soft time limits must be set to prevent worker starvation."""
    conf = celery_app.conf
    assert conf.task_time_limit == 3600, "hard time limit should be 1 hour"
    assert conf.task_soft_time_limit == 3300, "soft limit should be 55 minutes"
    assert conf.task_soft_time_limit < conf.task_time_limit
```

- [ ] **Step 2: Run test to verify it fails**

Run: `$PY -m pytest tests/test_celery_config.py -q -p no:cacheprovider`
Expected: FAIL (task_time_limit is None by default).

- [ ] **Step 3: Implement** — add to the `celery_app.conf.update(...)` call in `backend/app/celery_app.py`, after `task_reject_on_worker_lost=True,` (line 42):

```python
    # Global task time limits — backstop against pathological molecules that
    # make RDKit ring perception run unbounded and permanently block a worker.
    task_soft_time_limit=3300,  # 55 min: raises SoftTimeLimitExceeded (graceful)
    task_time_limit=3600,       # 60 min: hard SIGKILL of the worker process
```

- [ ] **Step 4: Run test to verify it passes**

Run: `$PY -m pytest tests/test_celery_config.py -q -p no:cacheprovider` → PASS

- [ ] **Step 5: Commit**

```bash
git add backend/app/celery_app.py backend/tests/test_celery_config.py
git commit -m "fix(celery): add global task time limits to prevent worker starvation"
```

---

### Task 2: Narrow Celery `autoretry_for` to transient errors

`autoretry_for=(Exception,)` retries *every* failure 2–3× — including permanent errors. (Molecule-level errors are already caught and returned as result dicts at `tasks.py:120,163`, so they never raise; only infra errors should trigger retries.) Narrow to transient connection/timeout errors.

**Files:**
- Modify: `backend/app/services/batch/tasks.py` (imports near top; decorators at lines 59-65, 169-175, 661-668)
- Test: `backend/tests/test_batch/test_retry_policy.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# backend/tests/test_batch/test_retry_policy.py
"""Batch tasks must only auto-retry transient infrastructure errors."""

from app.services.batch import tasks


def test_transient_errors_constant_excludes_generic_exception():
    assert Exception not in tasks.TRANSIENT_ERRORS
    assert ValueError not in tasks.TRANSIENT_ERRORS


def test_transient_errors_constant_includes_connection_errors():
    import redis.exceptions

    assert ConnectionError in tasks.TRANSIENT_ERRORS
    assert TimeoutError in tasks.TRANSIENT_ERRORS
    assert redis.exceptions.ConnectionError in tasks.TRANSIENT_ERRORS


def test_tasks_use_transient_errors_tuple():
    # Celery stores autoretry_for-derived retry behaviour; assert the decorators
    # reference the shared constant rather than a blanket (Exception,).
    import inspect

    src = inspect.getsource(tasks)
    assert "autoretry_for=(Exception,)" not in src
    assert src.count("autoretry_for=TRANSIENT_ERRORS") == 3
```

- [ ] **Step 2: Run test to verify it fails**

Run: `$PY -m pytest tests/test_batch/test_retry_policy.py -q -p no:cacheprovider`
Expected: FAIL (`TRANSIENT_ERRORS` undefined).

- [ ] **Step 3: Implement**

Add near the top of `backend/app/services/batch/tasks.py` (after existing imports):

```python
import redis.exceptions

# Only transient infrastructure failures are worth retrying. Permanent errors
# (bad SMILES, schema mismatch, logic bugs) are already caught inside the task
# bodies and returned as result dicts, so a blanket (Exception,) retry only
# triples the cost of permanently-broken molecules.
TRANSIENT_ERRORS = (
    ConnectionError,
    TimeoutError,
    redis.exceptions.ConnectionError,
    redis.exceptions.TimeoutError,
    OSError,
)
```

Then replace all three occurrences of `autoretry_for=(Exception,),` with `autoretry_for=TRANSIENT_ERRORS,` (decorators at lines 59-65, 169-175, 661-668).

- [ ] **Step 4: Run test to verify it passes**

Run: `$PY -m pytest tests/test_batch/test_retry_policy.py -q -p no:cacheprovider` → PASS

- [ ] **Step 5: Commit**

```bash
git add backend/app/services/batch/tasks.py backend/tests/test_batch/test_retry_policy.py
git commit -m "fix(batch): narrow Celery autoretry_for to transient errors only"
```

---

### Task 3: Stricter insecure-default-secret gate

Defaults (`SECRET_KEY` etc. = `"CHANGE_ME_IN_PRODUCTION"`) currently only hard-fail when `DEBUG=False`; in debug they emit a `warning` and proceed. Escalate to `ERROR`-level logging and add an `ALLOW_INSECURE_DEFAULTS` switch so an operator can force-fail even in debug. Default behavior (debug accepts defaults) is preserved so existing tests/dev flows are unaffected.

**Files:**
- Modify: `backend/app/core/config.py:153-174` (the `_check_insecure_defaults` validator) + add setting near line 44
- Test: `backend/tests/test_security_startup.py` (extend — read existing first to match style)

- [ ] **Step 1: Write the failing test** (append to `tests/test_security_startup.py`)

```python
def test_insecure_defaults_blocked_when_explicitly_disallowed():
    """Even in DEBUG, ALLOW_INSECURE_DEFAULTS=False must reject default secrets."""
    import pytest
    from app.core.config import Settings

    with pytest.raises(ValueError, match="insecure default"):
        Settings(DEBUG=True, ALLOW_INSECURE_DEFAULTS=False)


def test_insecure_defaults_allowed_in_debug_by_default():
    """Default dev behaviour is preserved: DEBUG=True still accepts defaults."""
    from app.core.config import Settings

    s = Settings(DEBUG=True)  # ALLOW_INSECURE_DEFAULTS defaults True
    assert s.SECRET_KEY == "CHANGE_ME_IN_PRODUCTION"
```

- [ ] **Step 2: Run to verify it fails**

Run: `$PY -m pytest tests/test_security_startup.py -q -p no:cacheprovider`
Expected: FAIL (`ALLOW_INSECURE_DEFAULTS` unknown field → first test does not raise the expected ValueError).

- [ ] **Step 3: Implement**

Add setting after `API_KEY_MAX_EXPIRY_DAYS` (line 47) in `config.py`:

```python
    # When False, insecure default secrets are rejected even in DEBUG mode.
    # Leave True for local dev; set False in any shared/staging environment.
    ALLOW_INSECURE_DEFAULTS: bool = True
```

Replace the body of `_check_insecure_defaults` (lines 156-174) with:

```python
        secret_fields = ("SECRET_KEY", "API_KEY_ADMIN_SECRET", "CSRF_SECRET_KEY")
        insecure = [f for f in secret_fields if getattr(self, f) in _INSECURE_DEFAULTS]
        if insecure:
            must_block = (not self.DEBUG) or (not self.ALLOW_INSECURE_DEFAULTS)
            if must_block:
                raise ValueError(
                    f"Insecure default value(s) for {', '.join(insecure)}. "
                    f"Set strong secrets via environment variables before running "
                    f"outside local development."
                )
            _config_logger.error(
                "%d security setting(s) use INSECURE DEFAULTS (%s). "
                "Acceptable ONLY for local development; set ALLOW_INSECURE_DEFAULTS=False "
                "in any shared environment.",
                len(insecure),
                ", ".join(insecure),
            )
        return self
```

- [ ] **Step 4: Run to verify it passes** (and the whole security-startup file still passes)

Run: `$PY -m pytest tests/test_security_startup.py -q -p no:cacheprovider` → PASS

- [ ] **Step 5: Commit**

```bash
git add backend/app/core/config.py backend/tests/test_security_startup.py
git commit -m "fix(security): escalate insecure-default warnings and add ALLOW_INSECURE_DEFAULTS gate"
```

---

### Task 4: Pool Redis connections (fix per-call creation + leak)

`get_sync_redis_client()` (`rate_limit.py:28`) and `get_redis_client()` (`security.py:33`) create a brand-new client on every call. The sync clients in `rate_limit.py` are never closed → fd leak under unique-IP traffic. Replace both with lazily-initialized module-level singletons backed by a connection pool. Keep the function names/signatures so existing `@patch(...)` tests keep working.

**Files:**
- Modify: `backend/app/core/rate_limit.py:28-30`
- Modify: `backend/app/core/security.py:33-35`, and remove the now-wrong `await client.aclose()` calls at `:94` and `:201` (they would close the shared client)
- Test: `backend/tests/test_redis_pooling.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# backend/tests/test_redis_pooling.py
"""Redis clients must be pooled singletons, not created per call."""

import app.core.rate_limit as rl


def test_sync_redis_client_is_singleton():
    rl._SYNC_REDIS = None  # reset module cache
    c1 = rl.get_sync_redis_client()
    c2 = rl.get_sync_redis_client()
    assert c1 is c2


async def test_async_redis_client_is_singleton():
    import app.core.security as sec

    sec._ASYNC_REDIS = None
    c1 = await sec.get_redis_client()
    c2 = await sec.get_redis_client()
    assert c1 is c2
```

- [ ] **Step 2: Run to verify it fails**

Run: `$PY -m pytest tests/test_redis_pooling.py -q -p no:cacheprovider`
Expected: FAIL (`c1 is c2` false — new client each call; `_SYNC_REDIS` attr missing).

- [ ] **Step 3: Implement**

`rate_limit.py` — replace lines 28-30:

```python
_SYNC_REDIS: "redis.Redis | None" = None


def get_sync_redis_client():
    """Return a process-wide pooled synchronous Redis client (lazy singleton)."""
    global _SYNC_REDIS
    if _SYNC_REDIS is None:
        _SYNC_REDIS = redis.from_url(
            settings.REDIS_URL, decode_responses=True, max_connections=20
        )
    return _SYNC_REDIS
```

`security.py` — replace lines 33-35:

```python
_ASYNC_REDIS: "redis.Redis | None" = None


async def get_redis_client():
    """Return a process-wide pooled async Redis client (lazy singleton)."""
    global _ASYNC_REDIS
    if _ASYNC_REDIS is None:
        _ASYNC_REDIS = redis.from_url(
            settings.REDIS_URL, decode_responses=True, max_connections=20
        )
    return _ASYNC_REDIS
```

Then in `security.py`, delete the two `finally:` `await client.aclose()` lines (currently at `:93-94` in `validate_api_key` and `:200-201` in `_update_usage_stats`) — the shared client must NOT be closed per request. Keep the `try:` bodies; remove only the `finally`/`aclose`. (Closing on shutdown is handled by the existing app lifespan / acceptable to leave to process exit.)

- [ ] **Step 4: Run to verify it passes** + the touched modules' existing tests

Run: `$PY -m pytest tests/test_redis_pooling.py tests/test_rate_limit.py tests/test_api_keys.py -q -p no:cacheprovider` → PASS (existing patches target the function names, unaffected)

- [ ] **Step 5: Commit**

```bash
git add backend/app/core/rate_limit.py backend/app/core/security.py backend/tests/test_redis_pooling.py
git commit -m "perf(redis): use pooled singleton clients instead of per-call connections"
```

---

### Task 5: Vectorize Butina clustering distance matrix

`clustering.py:117-120` computes pairwise Tanimoto distances with a nested Python loop (≈500k comparisons at the 1000-mol cap). Replace with RDKit's C-level `BulkTanimotoSimilarity`. Output must be identical (`1 - similarity`, same lower-triangle order).

**Files:**
- Modify: `backend/app/services/analytics/clustering.py:116-120`
- Test: `backend/tests/test_analytics_clustering.py` (extend — existing tests are the regression net)

- [ ] **Step 1: Write the failing/guard test** (append)

```python
def test_bulk_and_loop_distances_match():
    """Vectorized BulkTanimotoSimilarity must equal the naive loop, same order."""
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    smis = ["CCO", "CCN", "c1ccccc1", "CC(=O)O", "CCCl"]
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smis]
    n = len(fps)

    naive = []
    for i in range(1, n):
        for j in range(i):
            naive.append(1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))

    bulk = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        bulk.extend(1.0 - s for s in sims)

    assert len(naive) == len(bulk)
    for a, b in zip(naive, bulk):
        assert abs(a - b) < 1e-9
```

- [ ] **Step 2: Run to verify it passes already** (this is a property guard, proves equivalence before refactor)

Run: `$PY -m pytest tests/test_analytics_clustering.py::test_bulk_and_loop_distances_match -q -p no:cacheprovider` → PASS

- [ ] **Step 3: Implement** — replace lines 116-120 in `clustering.py`:

```python
    # Step 3: Compute lower-triangle Tanimoto distance list (vectorized via
    # RDKit C-level BulkTanimotoSimilarity — ~10-50x faster than a Python loop).
    dists: list[float] = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1.0 - s for s in sims)
```

- [ ] **Step 4: Run the full clustering suite** (regression)

Run: `$PY -m pytest tests/test_analytics_clustering.py -q -p no:cacheprovider` → PASS (all)

- [ ] **Step 5: Commit**

```bash
git add backend/app/services/analytics/clustering.py backend/tests/test_analytics_clustering.py
git commit -m "perf(clustering): vectorize Tanimoto distance matrix with BulkTanimotoSimilarity"
```

---

### Task 6: Offload blocking `task.get()` off the event loop

`validation.py:307` calls the synchronous, blocking `task.get(timeout=...)` inside an `async def` handler, stalling the FastAPI event loop for up to 60s. Push it to the default thread-pool executor.

**Files:**
- Modify: `backend/app/api/routes/validation.py` (around line 305-307; confirm `import asyncio` present at top, add if missing)
- Test: `backend/tests/test_validation/` — add a focused test if the endpoint is exercisable with eager Celery; otherwise rely on existing route tests + manual verification (document which).

- [ ] **Step 1: Confirm `asyncio` import** at top of `validation.py`; add `import asyncio` if absent.

- [ ] **Step 2: Write the test** (only if a Celery-eager path exists; check `tests/test_validation/`). If not feasible as a unit test, add a regression note in the commit body and verify via the full route suite. Pattern when feasible:

```python
async def test_async_validate_does_not_block_event_loop(client, monkeypatch):
    # With CELERY eager + a fake task whose .get sleeps, the handler must await
    # via run_in_executor (event loop stays responsive). See conftest celery setup.
    ...
```

- [ ] **Step 3: Implement** — replace line 307 (`result = task.get(timeout=timeout)`):

```python
        # Run the blocking Celery result fetch in a thread so the event loop
        # stays free to serve other requests.
        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(None, lambda: task.get(timeout=timeout))
```

- [ ] **Step 4: Run the validation route suite** (regression)

Run: `$PY -m pytest tests/test_validation tests/test_api -q -p no:cacheprovider` → no new failures vs baseline.

- [ ] **Step 5: Commit**

```bash
git add backend/app/api/routes/validation.py
git commit -m "perf(validation): run blocking Celery task.get in executor to free event loop"
```

---

### Task 7: Guard module-level `sascorer` imports against missing RDKit Contrib

Three modules do `sys.path.append(RDConfig.RDContribDir/SA_Score); import sascorer` at module load. If Contrib is absent, the *entire module* fails to import and crashes the worker. Convert to a lazy `_get_sascorer()` loader mirroring the existing good pattern in `scoring/np_likeness.py:63-75`.

**Files:**
- Modify: `backend/app/services/structure_filter/scorer.py:31` (import) + `:95` (usage)
- Modify: `backend/app/services/structure_filter/filter_pipeline.py:36-37` (import) + `:271` (usage)
- Modify: `backend/app/services/profiler/sa_comparison.py:26` (import) + `:53` (usage)
- Test: `backend/tests/test_sascorer_guard.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# backend/tests/test_sascorer_guard.py
"""Module-level sascorer imports must not crash when Contrib is unavailable."""

from rdkit import Chem


def test_filter_pipeline_exposes_lazy_sascorer_loader():
    from app.services.structure_filter import filter_pipeline

    assert hasattr(filter_pipeline, "_get_sascorer")
    scorer = filter_pipeline._get_sascorer()
    # When available (normal RDKit install) it computes; the contract is "never
    # raise at import time".
    if scorer is not None:
        assert isinstance(scorer.calculateScore(Chem.MolFromSmiles("CCO")), float)


def test_no_module_level_bare_sascorer_import():
    import inspect

    from app.services.structure_filter import filter_pipeline, scorer
    from app.services.profiler import sa_comparison

    for mod in (filter_pipeline, scorer, sa_comparison):
        src = inspect.getsource(mod)
        # The bare top-level "import sascorer" (col 0) must be gone.
        assert "\nimport sascorer" not in src, mod.__name__
```

- [ ] **Step 2: Run to verify it fails**

Run: `$PY -m pytest tests/test_sascorer_guard.py -q -p no:cacheprovider` → FAIL.

- [ ] **Step 3: Implement** — in each of the three files, replace the two lines

```python
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore  # noqa: E402
```

with a cached lazy loader:

```python
from functools import lru_cache


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
```

Update each usage site to resolve the loader and degrade gracefully:
- `scorer.py:95` and `filter_pipeline.py:271`:
  ```python
  _sa = _get_sascorer()
  if _sa is None:
      raise RuntimeError("SA Score unavailable: RDKit Contrib SA_Score not installed")
  sa_raw = _sa.calculateScore(mol)
  ```
  (Adapt the variable name `sa_raw`/`sa_score` to match each call site.)
- `sa_comparison.py:53`: same pattern; if the surrounding function returns a result dict, return `{"error": "sascorer not available"}` instead of raising, matching the SCScore convention in that file.

- [ ] **Step 4: Run to verify it passes** + the structure-filter/profiler suites

Run: `$PY -m pytest tests/test_sascorer_guard.py tests/test_structure_filter -q -p no:cacheprovider` → PASS / no new failures.

- [ ] **Step 5: Commit**

```bash
git add backend/app/services/structure_filter/scorer.py backend/app/services/structure_filter/filter_pipeline.py backend/app/services/profiler/sa_comparison.py backend/tests/test_sascorer_guard.py
git commit -m "fix(scoring): guard module-level sascorer imports against missing RDKit Contrib"
```

---

### Task 8: Add coverage thresholds (ratchet, do not break CI)

No coverage gate exists for backend or frontend. Add thresholds set at/just below *measured* current coverage so they ratchet without failing the current suite.

**Files:**
- Modify: `backend/pyproject.toml:72-81` (`[tool.pytest.ini_options]`)
- Modify: `frontend/vitest.config.ts:11-20` (`coverage` block)

- [ ] **Step 1: Measure backend coverage**

Run: `$PY -m pytest --cov=app --cov-report=term-missing -q -p no:cacheprovider 2>&1 | tail -5`
Record the TOTAL %.

- [ ] **Step 2: Measure frontend coverage**

Run (from `frontend/`): `npm run test -- --coverage --run 2>&1 | tail -20`
Record the lines %.

- [ ] **Step 3: Implement backend gate** — add to `[tool.pytest.ini_options]` `addopts` (set `<N>` = floor(measured − 2), so the gate ratchets but the current suite passes; note the 15 env-failures must be excluded from the cov run or the gate measured on the green subset):

```toml
addopts = "--cov=app --cov-report=term-missing --cov-fail-under=<N>"
```

- [ ] **Step 4: Implement frontend gate** — add to the `coverage` block in `vitest.config.ts`:

```ts
      thresholds: {
        lines: <M>,
        functions: <M>,
        statements: <M>,
      },
```

- [ ] **Step 5: Verify gates pass at the chosen thresholds**, then commit:

```bash
git add backend/pyproject.toml frontend/vitest.config.ts
git commit -m "test: add coverage thresholds to ratchet backend and frontend coverage"
```

> **Note:** if enabling `--cov` by default conflicts with the env-failing integration tests, scope the gate to the unit subset or document the required `-m "not integration"` selector in the commit body. Surface the actual measured numbers — do not silently pick a low threshold.

---

### Wave 1 exit criteria

- [ ] All Wave 1 tasks committed atomically.
- [ ] `$PY -m pytest -q -p no:cacheprovider` shows ≥ 2262 passing, ≤ 15 failing (no NEW failures vs baseline).
- [ ] `cd frontend && npm test -- --run` green (unaffected by backend work, but confirm).
- [ ] Summary posted to user; **checkpoint before Wave 2**.

---

## Wave 2 — Defensive hardening (medium) — ✅ ALL DONE (2.1–2.6)

- **2.1 `progress_tracker` thread-safety** — ✅ DONE. Double-checked locking around the lazy Redis init (`threading.Lock`). Commit `fix(batch): make ProgressTracker lazy Redis init thread-safe`.
- **2.2 OPSIN JVM init-failed flag** — ✅ DONE. `_init_failed` sticky flag; first failure re-raises (startup logs), subsequent calls short-circuit. Commit `fix(iupac): skip OPSIN re-init after failure...`.
- **2.3 SYBA subprocess input hardening** — ✅ DONE. Reject oversized/non-printable SMILES before spawn (`_MAX_SYBA_SMILES_LEN=2000`), explicit `shell=False`. Commit `fix(profiler): validate SYBA SMILES input...`.
- **2.4 Admin-secret replay resistance** — ✅ DONE (backward-compatible). Added `generate_admin_token()` / `verify_admin_secret()`: HMAC-signed, time-bound tokens (`ADMIN_AUTH_MAX_SKEW_SECONDS=300`) always honoured; opt-in strict mode `ADMIN_AUTH_REQUIRE_SIGNED` (default False) rejects the static secret. Existing clients unaffected. Commit `feat(security): add opt-in HMAC-signed admin tokens...`.
- **2.5 SCScore vendor weights** — ✅ DONE (2026-06-10, superseded the block). Investigation showed the upstream `scscore` package is Python-2 / NumPy-1 era (`import cPickle`, `dtype=np.bool`) and cannot run on this stack, and its loader has no `.npz` path — so the original "supply an `.npz`" plan was unworkable. Resolved like `molvs` (2.6): vendored the `full_reaxys_model_1024bool` weights as a pickle-free, NumPy-portable `scscore_weights.npz` (from upstream's `json.gz`; MIT, attributed) and reimplemented the ~25-line inference (Morgan-2/1024 FP → MLP → sigmoid 1-5) in `sa_comparison.py`, dropping the broken import. Reproduces the upstream model bit-for-bit (validated < 1e-6); still degrades gracefully if the weight file is absent. Commit `009d221`.
- **2.6 `molvs` migration audit** — ✅ DONE (2026-06-10, superseded the deferral). Rather than swap to `rdkit.Chem.MolStandardize`, the unmaintained `molvs` dependency was dropped entirely: its 61 `REMOVE_FRAGMENTS` patterns were vendored verbatim into `services/validation/_fragment_patterns.py` with an upstream-parity test. Commit `c526c22 refactor(validation): vendor MolVS fragment patterns, drop unmaintained molvs dependency`.

---

## Wave 3 — Larger backend refactors — ✅ ALL DONE (3.1–3.7)

- **3.1 Chunked Redis batch-result storage** — ✅ DONE. Results stored as a Redis LIST (one JSON per element); default unfiltered/unsorted view paginates via `LRANGE` (no full deserialize); filtered/sorted views load-all + view-cache; backward-compatible with the legacy single-blob format. Commit `perf(batch): store results as Redis list...`.
- **3.2 WebSocket duplicate-subscriber fix + test suite** — ✅ DONE. Ownership guard: a superseded subscriber stops before broadcasting; cancel-and-await prior task in `connect`; full WS test suite added (was zero coverage). Commit `fix(websockets): prevent duplicate broadcasts...`.
- **3.3 Persistent SYBA worker** — ✅ DONE (2026-06-10, superseded the deferral). Subprocess-per-call replaced with a long-lived `services/profiler/syba_worker.py` (one request per stdin line / one response per stdout line) with restart-on-death and per-call timeout. Commit `24057f5 perf(profiler): replace per-call SYBA subprocess with persistent worker (JSON pipe, restart, timeout)`.
- **3.4 Dedicated analytics queue + time limits** — ✅ DONE. `analytics` queue + routing for `run_expensive_analytics`, per-task `soft_time_limit=600`/`time_limit=660`, dedicated `celery-worker-analytics` added to both compose files. Commit `perf(analytics): isolate expensive analytics...`.
- **3.5 `asyncio.run()` in Celery tasks** — ✅ DONE (2026-06-10, superseded the deferral). Added a dedicated sync (`psycopg`) DB engine for Celery tasks and removed `asyncio.run` from the audit/purge task code, making it robust across all Celery pool types. Commit `417555f refactor(db): add sync engine for Celery tasks, remove asyncio.run from audit/purge`.
- **3.6 Celery chord/aggregation integration test** — ✅ DONE. Eager-mode end-to-end test of `group(chunks)` → aggregate → store. Commit `test(batch): add eager Celery chord aggregation integration test`.
- **3.7 Redis memory policy** — ✅ DONE (improved over audit). Bounded `maxmemory` in dev compose; switched BOTH compose files from the risky `allkeys-lru` to `volatile-lru` so no-TTL Celery broker messages are never evicted. Commit `fix(redis): bound memory and use volatile-lru...`.

## Wave 4 — Frontend refactors & E2E — ✅ 4.1/4.2/4.3/4.4/4.5 ALL DONE (4.1 finished 2026-06-10)

- **4.1 Split `SingleValidation.tsx`** (3,115 lines) — ✅ DONE (2026-06-10). Reduced to **1,802 lines** by extracting presentational components one per commit, each verified by tsc + eslint + vitest (202) + Playwright E2E (1). State stayed in the page; components are presentational and receive props. Extractions: `ExampleMolecules` (`acff4df`), `StatusPanels` (`2facba2`), `ScoreTiles` (`e06d0a9`), `ValidationOutcomePanels` (`bd3b090`), `ValidationIssuesPanel` (`0774df7`), `MoleculeViewerPanel` (`1ac9335`), `AllChecksCard` (`1551f55`), `integrations/DatabaseLookupTab` = DatabaseLookupControls + DatabaseResultsBox (`93b5faf`), `validation/TabInputPanels` = six per-tab input panels (`b5031ed`).
- **4.2 Split monolithic backend route files** — ✅ DONE. `scoring.py` 893→392 (17 builders → `services/scoring/score_builders.py`); `qsar_ready.py` 707→520 (export/summary builders → `services/qsar_ready/builders.py`); `batch.py` 1,090→824 (subset endpoints → `routes/batch_subset.py`, shared `services/batch/parsing.py`). Verified against the full 2,341-test real-service suite.
- **4.3 Replace `window.prompt`/`window.confirm`** — ✅ DONE. Accessible, testable `ConfirmModal` + `PromptModal` wired into all three components.
- **4.4 Page-level frontend tests** — ✅ DONE. Render/characterization tests for `SingleValidation` and `BatchValidation` in the real provider tree.
- **4.5 E2E suite (Playwright)** — ✅ DONE. `frontend/playwright.config.ts` orchestrates backend (:8001) + frontend (:3002); `e2e/single-validation.spec.ts` drives a real browser through the full validate round-trip. Passes against the live stack (Redis + Postgres containers).

## Environment note (this session)

The full backend suite was validated against a **real stack**: Redis + Postgres started as Docker containers (`chemaudit-redis-test`, `chemaudit-pg-test`), and the missing-in-this-env Python deps were installed into the `cheminformatics` conda env (`asyncpg`, `psycopg[binary]`, `fakeredis`, `aiosqlite`, `pytest-cov`, `coverage`, `xlsxwriter`, `openpyxl`, `openTSNE`, plus `@playwright/test` + chromium for the frontend). Result: **2,341 backend tests pass, 0 failures** (was 15 env-only failures before the stack was up), **202 frontend tests**, **1 E2E** — all green.

---

## Wave 4 — Frontend refactors & E2E (gated, likely separate PRs)

- **4.1 Split `SingleValidation.tsx`** (3,115 lines, `frontend/src/pages/`): extract tab panels (Scoring, Alerts, Standardization, …) into `frontend/src/components/validation/` sub-components; keep the page as an orchestrator. Large refactor — own PR, incremental, snapshot/behavior tests per extraction.
- **4.2 Split monolithic backend route files**: `batch.py` (1,090), `scoring.py` (893), `qsar_ready.py` (707) → extract inline business logic into `services/`, leaving thin handlers. Own PRs.
- **4.3 Replace `window.prompt`/`window.confirm`** (`components/dataset-audit/WeightSliders.tsx:73,83`, `components/qsar-ready/PipelineConfigPanel.tsx:61,69`, `components/structure-filter/StructureFilterPresetSelector.tsx:98,107`) with accessible, testable React modal dialogs.
- **4.4 Page-level frontend tests**: `SingleValidation.tsx`, `BatchValidation.tsx` (currently untested).
- **4.5 E2E suite (Playwright)**: upload → WebSocket progress → results → export happy path. `@vitest/browser-playwright` is already present in the lockfile.

---

## Self-Review (against CONCERNS.md)

- **Tech Debt:** Redis per-call → T4; `asyncio.run()` in tasks → 3.5; monolithic routes → 4.2; `SingleValidation.tsx` → 4.1; no coverage enforcement → T8. ✅ covered.
- **Security:** insecure defaults → T3; WS ownership fail-open → deliberate (kept for parity with the tested HTTP path; rationale now documented in `_check_ws_ownership`, commit `b644e5a`); SYBA subprocess → 2.3; admin header → 2.4. ✅ covered.
- **Performance:** batch JSON blob → 3.1; blocking `task.get` → T6; O(n²) clustering → T5; SYBA subprocess/call → 3.3; t-SNE blocking → 3.4. ✅ covered.
- **Fragile areas:** WS duplicate subscriber → 3.2; no WS tests → 3.2; progress_tracker thread-safety → 2.1; `autoretry_for=(Exception,)` → T2; OPSIN JVM → 2.2. ✅ covered.
- **Scaling:** Redis SPOF/memory → 3.7; clustering O(n²) memory → vectorized via BulkTanimotoSimilarity + cap now deployment-configurable (`CLUSTERING_MAX_MOLECULES`, commit `d84c3fd`). ✅ covered.
- **Dependencies:** `sascorer`/`npscorer` import → T7; `scscore` NumPy 2.x/vendor → 2.5 (DONE — vendored weights + reimplemented inference); `molvs` → 2.6 (DONE). ✅ covered.
- **Missing features:** no E2E → 4.5; no Celery time limits → T1; Celery integration tests → 3.6. ✅ covered.
- **Test gaps:** WS manager → 3.2; frontend pages → 4.4; Celery aggregation → 3.6; `window.prompt` → 4.3; rate-limit conn close → T4. ✅ covered.

Every CONCERNS.md item is now resolved and committed (Waves 1–4, all sub-items
including 2.5/2.6/3.3/3.5 and all of 4.1). No items remain blocked or deferred.
Branch `fix/codebase-audit-concerns` is green: backend 2346 passed / 0 failed /
31 skipped; frontend tsc + eslint clean, vitest 202, Playwright E2E 1.
