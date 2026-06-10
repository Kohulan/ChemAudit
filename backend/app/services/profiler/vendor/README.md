# SCScore weights (vendored)

SCScore (Coley et al., *J. Chem. Inf. Model.* 2018) is the third synthesizability
score shown alongside SA Score and SYBA in the profiler's SA-comparison view.

ChemAudit does **not** depend on the upstream `scscore` package. That package is
Python-2 / NumPy-1 era code (`import cPickle`, `dtype=np.bool`) and does not run
under this project's Python 3 + NumPy 2.x stack. Instead — mirroring how `molvs`
was handled — the trained weights are vendored here and the (tiny) inference is
reimplemented self-contained in `../sa_comparison.py`.

## What's in this directory

- **`scscore_weights.npz`** — the trained weights of the upstream
  `full_reaxys_model_1024bool` model (`model.ckpt-10654`). A 6-layer MLP over a
  1024-bit radius-2 Morgan fingerprint (12 arrays: 6 weight/bias pairs). Stored
  as a compressed, **pickle-free** `.npz` of plain `.npy` arrays, so it loads
  identically across NumPy versions.

The reimplementation (`_scscore_fingerprint` + `_scscore_forward` in
`sa_comparison.py`) reproduces the upstream `standalone_model_numpy` output
**bit-for-bit** — validated to `< 1e-6` against the reference model across a
range of molecules (see `tests/test_scscore_degradation.py`).

## License / attribution

SCScore is MIT-licensed (Connor W. Coley et al.,
<https://github.com/connorcoley/scscore>). The vendored weights are that
project's trained artifact, redistributed here under its MIT terms. Cite:

> Coley, C. W.; Rogers, L.; Green, W. H.; Jensen, K. F. *SCScore: Synthetic
> Complexity Learned from a Reaction Corpus.* J. Chem. Inf. Model. 2018, 58 (2),
> 252–261. <https://doi.org/10.1021/acs.jcim.7b00622>

## Regenerating the weights

If the model ever needs re-vendoring, convert the upstream NumPy-agnostic
`json.gz` weights (no NumPy 1.x environment required — `.json.gz` is plain text):

```python
import gzip, json, numpy as np
src = "scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz"
with gzip.GzipFile(src) as f:
    arrs = [np.asarray(x, dtype=np.float64) for x in json.loads(f.read().decode("utf-8"))]
np.savez_compressed("scscore_weights.npz", *arrs)   # keys arr_0 ... arr_11, in order
```

If `scscore_weights.npz` is absent, SCScore reports `{"available": false}` and the
other synthesizability scores continue to work — it never raises.
