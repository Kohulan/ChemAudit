# SCScore weights (optional)

SCScore (Coley et al.) is an **optional** third synthesizability score shown
alongside SA Score and SYBA in the profiler's SA-comparison view. It is *not* a
PyPI package and ships trained weights that are pickled under NumPy 1.x, which
fails to load under NumPy 2.x. ChemAudit therefore treats SCScore as optional:
when it cannot be loaded, the API returns `{"available": false, "error":
"scscore not available"}` for the SCScore field and the other scores still work.
When SCScore is installed but its weights fail to load, a one-time `WARNING` is
logged (see `_load_scscore` in `../sa_comparison.py`).

## Enabling SCScore under NumPy 2.x

1. Install SCScore from source (it is not on PyPI):

   ```bash
   pip install "git+https://github.com/connorcoley/scscore.git"
   ```

2. If `SCScorer().restore()` works in your environment (NumPy 1.x), nothing more
   is needed — the standard path is used automatically.

3. Under **NumPy 2.x**, the bundled weights fail to unpickle. Convert them once
   to a NumPy-2-compatible `.npz` and drop it next to this file as
   `scscore_weights.npz`. Using a NumPy 1.x environment (or `allow_pickle`):

   ```python
   # Run in an env where scscore's weights load (e.g. numpy<2):
   import numpy as np
   from scscore.standalone_model_numpy import SCScorer
   s = SCScorer()
   s.restore()                      # loads the bundled weights
   # s.vars is a list of weight/bias numpy arrays
   np.savez("scscore_weights.npz", *s.vars)
   ```

   Then copy `scscore_weights.npz` into this directory. `_load_scscore` will pick
   it up via `scorer.restore(weight_path=...)`.

> The weights themselves are **not vendored** here: they are a third-party
> trained-model artifact with its own licensing, and obtaining/redistributing
> them is a deployment decision for the operator.
