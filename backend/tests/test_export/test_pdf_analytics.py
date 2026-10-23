"""Tests for the optional analytics sections of the batch PDF report."""

import base64

from app.services.export.pdf_analytics import ANALYTICS_SECTIONS, build_analytics_sections

SMILES = [
    "CCc1ccccc1",
    "NCCc1ccccc1",
    "OCCc1ccccc1",
    "Clc1ccccc1",
    "Oc1ccccc1",
    "Cc1ccncc1",
    "Clc1ccncc1",
    "C1CCCCC1",
    "OC1CCCCC1",
    "c1ccc2ccccc2c1",
    "Cc1ccc2ccccc2c1",
    "CC(=O)Nc1ccc(O)cc1",
]


def _results():
    """Batch-shaped results with the numeric properties the sections use."""
    return [
        {
            "index": i,
            "status": "success",
            "smiles": smi,
            "validation": {"overall_score": 60 + 3 * i},
            "scoring": {
                "druglikeness": {"qed_score": 0.3 + 0.04 * i},
                "admet": {"sa_score": 1.0 + 0.1 * (i % 5), "fsp3": 0.05 * (i % 7)},
                "ml_readiness": {"score": 70 + i},
            },
        }
        for i, smi in enumerate(SMILES)
    ]


def test_builds_every_section_with_content():
    """Each selected section yields data, and charts are valid SVG."""
    out = build_analytics_sections(_results(), set(ANALYTICS_SECTIONS))

    assert set(out) == set(ANALYTICS_SECTIONS)
    for name in ("chemical_space", "scaffold_treemap", "correlation_matrix"):
        assert "error" not in out[name], out[name]
        assert base64.b64decode(out[name]["chart"]).startswith(b"<svg")
    assert out["chemical_space"]["molecule_count"] == len(SMILES)
    assert out["scaffold_treemap"]["top"][0]["count"] >= 1
    assert {r["label"] for r in out["statistics"]["rows"]} >= {"Validation score", "QED"}
    assert out["mmp_pairs"]["rows"][0]["smiles_a"] in SMILES


def test_only_selected_sections_are_computed():
    """Unselected sections are skipped."""
    out = build_analytics_sections(_results(), {"statistics"})

    assert set(out) == {"statistics"}


def test_small_batch_reports_errors_instead_of_failing():
    """Too little data gives a readable note, not an exception."""
    out = build_analytics_sections(_results()[:2], set(ANALYTICS_SECTIONS))

    assert "error" in out["chemical_space"]
    assert "error" in out["correlation_matrix"]


def test_scaffold_labels_are_escaped_in_svg():
    """SMILES with SVG-special characters do not break the chart."""
    results = _results()
    results[0]["smiles"] = "C[C@H](N)c1ccc(cc1)C<"  # unparseable: skipped
    out = build_analytics_sections(results, {"scaffold_treemap"})

    svg = base64.b64decode(out["scaffold_treemap"]["chart"]).decode()
    assert "<svg" in svg and "C<" not in svg
