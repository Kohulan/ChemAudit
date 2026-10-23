"""
Analytics sections for the batch PDF report.

Builds template-ready data for the optional PDF sections (chemical space,
scaffolds, property statistics, correlation matrix, MMP pairs) from the same
analytics services the web UI uses. Charts are plain SVG, like the score
distribution chart, so no plotting library is needed.
"""

import base64
import logging
from html import escape
from typing import Any, Callable, Dict, List, Optional

logger = logging.getLogger(__name__)

# Validation score bands: (lower bound, color, label). Shared with the score
# distribution chart in pdf_report.py.
SCORE_BANDS = [
    (90, "#16a34a", "90-100"),
    (70, "#2563eb", "70-89"),
    (50, "#d97706", "50-69"),
    (0, "#dc2626", "0-49"),
]
_NO_SCORE_COLOR = "#9ca3af"
_BAR_COLOR = "#2563eb"

# ponytail: SVG size grows with points; sample above this to keep the PDF small
MAX_SCATTER_POINTS = 5000
TOP_SCAFFOLDS = 15
TOP_MMP_PAIRS = 20

PROPERTY_LABELS = {
    "validation_score": "Validation score",
    "qed_score": "QED",
    "sa_score": "SA score",
    "ml_readiness_score": "ML-readiness",
    "fsp3": "Fsp3",
}


def build_analytics_sections(
    results: List[Dict[str, Any]], sections: set
) -> Dict[str, Dict[str, Any]]:
    """
    Compute the selected analytics sections for the PDF report.

    Each section is computed on its own; a failure in one section is reported
    in that section (``{"error": ...}``) and does not stop the report.

    Args:
        results: Batch result dicts being exported.
        sections: Selected section ids.

    Returns:
        Dict keyed by section id with template-ready data.
    """
    out: Dict[str, Dict[str, Any]] = {}
    for name, builder in SECTION_BUILDERS.items():
        if name not in sections:
            continue
        try:
            out[name] = builder(results)
        except Exception as e:
            logger.warning("PDF section %s failed: %s", name, e)
            out[name] = {"error": "This section could not be computed for this batch."}
    return out


def _b64(svg: str) -> str:
    """Base64-encode an SVG string for an <img> data URI."""
    return base64.b64encode(svg.encode("utf-8")).decode("utf-8")


def _score_color(score: Optional[float]) -> str:
    """Color for a validation score band."""
    if score is None:
        return _NO_SCORE_COLOR
    return next(color for threshold, color, _ in SCORE_BANDS if score >= threshold)


def _label(prop: str) -> str:
    """Display name for an analytics property key."""
    return PROPERTY_LABELS.get(prop, prop)


def _truncate(text: str, limit: int) -> str:
    """Shorten long SMILES for chart labels."""
    return text if len(text) <= limit else text[: limit - 1] + "…"


# ---------------------------------------------------------------------------
# Chemical space
# ---------------------------------------------------------------------------


def _chemical_space(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """PCA of Morgan fingerprints, drawn as a scatter colored by validation score."""
    from app.services.analytics.chemical_space import compute_pca

    # Sample before PCA: the fingerprint matrix is n x 2048, too big for 1M molecules
    candidates = [r for r in results if r.get("status") == "success"]
    step = max(1, -(-len(candidates) // MAX_SCATTER_POINTS))  # ceil division
    sample = candidates[::step]

    pca = compute_pca(sample)
    if "error" in pca:
        return {"error": pca["error"]}

    scores = {r.get("index"): (r.get("validation") or {}).get("overall_score") for r in sample}
    points = list(zip(pca["molecule_indices"], pca["coordinates"]))
    variance = pca["variance_explained"]
    return {
        "chart": _b64(_scatter_svg(points, scores, variance)),
        "molecule_count": len(candidates),
        "plotted_count": len(points),
        "variance_explained": [round(v * 100, 1) for v in variance],
    }


def _scatter_svg(points, scores, variance) -> str:
    """Scatter plot of PCA coordinates."""
    width, height, pad = 600, 380, 50
    xs = [c[0] for _, c in points]
    ys = [c[1] for _, c in points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_span = (x_max - x_min) or 1.0
    y_span = (y_max - y_min) or 1.0

    def sx(x):
        return pad + (x - x_min) / x_span * (width - 2 * pad)

    def sy(y):
        return height - pad - (y - y_min) / y_span * (height - 2 * pad)

    parts = [
        f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
        f'<rect x="{pad}" y="{pad}" width="{width - 2 * pad}" height="{height - 2 * pad}" '
        'fill="none" stroke="#d1d5db" />',
    ]
    for idx, (x, y) in points:
        parts.append(
            f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="2.5" '
            f'fill="{_score_color(scores.get(idx))}" fill-opacity="0.75" />'
        )
    parts.append(
        f'<text x="{width / 2}" y="{height - 12}" text-anchor="middle" font-size="12" '
        f'fill="#374151">PC1 ({variance[0] * 100:.1f}% variance)</text>'
    )
    parts.append(
        f'<text x="16" y="{height / 2}" text-anchor="middle" font-size="12" fill="#374151" '
        f'transform="rotate(-90 16 {height / 2})">PC2 ({variance[1] * 100:.1f}% variance)</text>'
    )
    for i, (_, color, label) in enumerate(SCORE_BANDS):
        x = pad + i * 90
        parts.append(f'<circle cx="{x + 5}" cy="25" r="5" fill="{color}" />')
        parts.append(f'<text x="{x + 15}" y="29" font-size="11" fill="#374151">{label}</text>')
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Scaffolds
# ---------------------------------------------------------------------------


def _scaffolds(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Murcko scaffold frequency: top scaffolds as a horizontal bar chart."""
    from app.services.analytics.scaffold_analysis import compute_scaffold_analysis

    analysis = compute_scaffold_analysis([r for r in results if r.get("status") == "success"])
    groups = sorted(analysis["scaffolds"], key=lambda g: g["count"], reverse=True)
    if not groups:
        return {"error": "No scaffolds found (no valid molecules)."}

    total = sum(g["count"] for g in groups)
    top = [
        {
            "smiles": g["scaffold_smiles"] or "(acyclic)",
            "count": g["count"],
            "percent": round(g["count"] / total * 100, 1),
        }
        for g in groups[:TOP_SCAFFOLDS]
    ]
    return {
        "chart": _b64(_bar_svg(top)),
        "unique_count": analysis["unique_scaffold_count"],
        "shannon_entropy": round(analysis["shannon_entropy"], 2),
        "top": top,
        "molecule_count": total,
    }


def _bar_svg(rows: List[Dict[str, Any]]) -> str:
    """Horizontal bar chart: one bar per scaffold."""
    width, label_w, bar_h, gap = 600, 260, 16, 6
    height = len(rows) * (bar_h + gap) + 10
    max_count = max(r["count"] for r in rows) or 1
    bar_space = width - label_w - 70

    parts = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    for i, row in enumerate(rows):
        y = 5 + i * (bar_h + gap)
        w = max(1, row["count"] / max_count * bar_space)
        parts.append(
            f'<text x="{label_w - 8}" y="{y + 12}" text-anchor="end" font-size="9" '
            f'font-family="monospace" fill="#374151">{escape(_truncate(row["smiles"], 42))}</text>'
        )
        parts.append(
            f'<rect x="{label_w}" y="{y}" width="{w:.1f}" height="{bar_h}" fill="{_BAR_COLOR}" />'
        )
        parts.append(
            f'<text x="{label_w + w + 6:.1f}" y="{y + 12}" font-size="10" fill="#374151">'
            f"{row['count']} ({row['percent']}%)</text>"
        )
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Property statistics and correlations
# ---------------------------------------------------------------------------


def _statistics(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Descriptive statistics per numeric property."""
    from app.services.analytics.statistics import compute_property_stats

    stats = compute_property_stats(results)
    if not stats:
        return {"error": "Not enough numeric property values (need at least 2 per property)."}
    return {
        "rows": [
            {
                "label": _label(s.property_name),
                **{
                    k: round(v, 2) if isinstance(v, float) else v for k, v in s.model_dump().items()
                },
            }
            for s in stats
        ]
    }


def _correlation_matrix(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Pearson correlation heatmap between numeric properties."""
    from app.services.analytics.statistics import compute_correlations

    correlations = compute_correlations(results)
    if not correlations:
        return {
            "error": "Not enough data: needs at least 2 properties with 10 or more values each."
        }

    names: List[str] = []
    for c in correlations:
        for n in (c.property_a, c.property_b):
            if n not in names:
                names.append(n)
    matrix = {(c.property_a, c.property_b): c.pearson_r for c in correlations}
    matrix.update({(b, a): r for (a, b), r in list(matrix.items())})
    return {"chart": _b64(_heatmap_svg(names, matrix))}


def _heat_color(r: float) -> str:
    """Diverging color: blue for positive, red for negative correlation."""
    t = min(abs(r), 1.0)
    base = (37, 99, 235) if r >= 0 else (220, 38, 38)
    rgb = [round(255 + (c - 255) * t) for c in base]
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def _heatmap_svg(names: List[str], matrix: Dict[tuple, float]) -> str:
    """Square heatmap with the r value printed in each cell."""
    label_w, top, cell = 120, 24, 70
    width = label_w + cell * len(names)
    height = top + cell * len(names)
    parts = [
        f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
    ]
    for j, name in enumerate(names):
        x = label_w + j * cell + cell / 2
        parts.append(
            f'<text x="{x}" y="{top - 8}" text-anchor="middle" font-size="10" '
            f'fill="#374151">{escape(_label(name))}</text>'
        )
    for i, row_name in enumerate(names):
        y = top + i * cell
        parts.append(
            f'<text x="{label_w - 8}" y="{y + cell / 2 + 4}" text-anchor="end" font-size="10" '
            f'fill="#374151">{escape(_label(row_name))}</text>'
        )
        for j, col_name in enumerate(names):
            r = 1.0 if i == j else matrix.get((row_name, col_name), 0.0)
            x = label_w + j * cell
            text_color = "#ffffff" if abs(r) > 0.6 else "#111827"
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" '
                f'fill="{_heat_color(r)}" stroke="#ffffff" />'
            )
            parts.append(
                f'<text x="{x + cell / 2}" y="{y + cell / 2 + 4}" text-anchor="middle" '
                f'font-size="11" fill="{text_color}">{r:.2f}</text>'
            )
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Matched molecular pairs
# ---------------------------------------------------------------------------


def _mmp_pairs(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Top matched molecular pairs by Tanimoto similarity."""
    from app.services.analytics.mmp import compute_mmp_analysis

    analysis = compute_mmp_analysis(results)
    if analysis.get("status") == "refused":
        return {"error": f"{analysis['reason']}; this batch has {analysis['molecule_count']}."}

    pairs = analysis.get("pairs") or []
    if not pairs:
        return {"error": "No matched molecular pairs found."}

    smiles_by_index = {r.get("index"): r.get("smiles", "") for r in results}
    return {
        "total": len(pairs),
        "rows": [
            {
                "mol_a": p["mol_a_index"],
                "mol_b": p["mol_b_index"],
                "smiles_a": smiles_by_index.get(p["mol_a_index"], ""),
                "smiles_b": smiles_by_index.get(p["mol_b_index"], ""),
                "core": p["core_smiles"],
                "change": f"{p['rgroup_a']} → {p['rgroup_b']}",
                "tanimoto": round(p["tanimoto"], 2),
            }
            for p in pairs[:TOP_MMP_PAIRS]
        ],
    }


# Section id -> builder, in report order
SECTION_BUILDERS: Dict[str, Callable[[List[Dict[str, Any]]], Dict[str, Any]]] = {
    "chemical_space": _chemical_space,
    "scaffold_treemap": _scaffolds,
    "statistics": _statistics,
    "correlation_matrix": _correlation_matrix,
    "mmp_pairs": _mmp_pairs,
}
ANALYTICS_SECTIONS = tuple(SECTION_BUILDERS)
