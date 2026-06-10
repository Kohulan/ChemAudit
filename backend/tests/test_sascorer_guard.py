"""Module-level sascorer imports must not crash when Contrib is unavailable."""

import inspect

from rdkit import Chem


def test_filter_pipeline_exposes_lazy_sascorer_loader():
    from app.services.structure_filter import filter_pipeline

    assert hasattr(filter_pipeline, "_get_sascorer")
    scorer = filter_pipeline._get_sascorer()
    # When available (normal RDKit install) it computes; the hard contract is
    # "never raise at import time".
    if scorer is not None:
        assert isinstance(scorer.calculateScore(Chem.MolFromSmiles("CCO")), float)


def test_all_three_modules_expose_lazy_loader():
    from app.services.profiler import sa_comparison
    from app.services.structure_filter import filter_pipeline, scorer

    for mod in (filter_pipeline, scorer, sa_comparison):
        assert hasattr(mod, "_get_sascorer"), mod.__name__


def test_no_module_level_bare_sascorer_import():
    from app.services.profiler import sa_comparison
    from app.services.structure_filter import filter_pipeline, scorer

    for mod in (filter_pipeline, scorer, sa_comparison):
        src = inspect.getsource(mod)
        # The bare top-level "import sascorer" (col 0) must be gone; the import
        # now lives inside the guarded _get_sascorer() loader.
        assert "\nimport sascorer" not in src, mod.__name__
