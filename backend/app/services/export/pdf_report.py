"""
PDF Report Generator

Generates professional PDF reports from batch validation results using WeasyPrint.
Includes summary statistics, score distribution chart, critical issues table,
and molecule structure images for flagged compounds.
"""
import base64
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional

from jinja2 import Environment, FileSystemLoader
from rdkit import Chem
from rdkit.Chem import Draw
from weasyprint import HTML

from app.services.export.base import BaseExporter


class PDFReportGenerator(BaseExporter):
    """
    Generate PDF reports for batch validation results.

    Uses WeasyPrint to convert HTML template to PDF, with embedded
    molecule images (base64 PNG) and SVG score distribution chart.
    """

    def __init__(self):
        """Initialize PDF generator with template environment."""
        # Find templates directory relative to this file
        template_dir = Path(__file__).parent.parent.parent / "templates" / "reports"
        self.env = Environment(loader=FileSystemLoader(str(template_dir)))
        self.template = self.env.get_template("batch_report.html")

    @property
    def media_type(self) -> str:
        """MIME type for PDF files."""
        return "application/pdf"

    @property
    def file_extension(self) -> str:
        """File extension for PDF files."""
        return "pdf"

    def export(self, results: List[Dict[str, Any]]) -> BytesIO:
        """
        Generate PDF report from batch results.

        Args:
            results: List of batch result dictionaries

        Returns:
            BytesIO buffer containing PDF data
        """
        # Calculate statistics
        stats = self._calculate_statistics(results)

        # Generate chart
        chart_data = self._generate_score_distribution_chart(results)

        # Extract critical issues (top 20)
        critical_issues = self._extract_critical_issues(results, limit=20)

        # Get flagged molecules (score < 70, limit 30)
        flagged_molecules = self._get_flagged_molecules(results, score_threshold=70, limit=30)

        # Render HTML template
        html_content = self.template.render(
            job_id=self._extract_job_id(results),
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            stats=stats,
            chart_data=chart_data,
            critical_issues=critical_issues,
            flagged_molecules=flagged_molecules,
        )

        # Convert HTML to PDF
        pdf_buffer = BytesIO()
        HTML(string=html_content, base_url=str(Path(__file__).parent)).write_pdf(pdf_buffer)
        pdf_buffer.seek(0)

        return pdf_buffer

    def _calculate_statistics(self, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Calculate summary statistics from results.

        Args:
            results: List of batch result dictionaries

        Returns:
            Dictionary with statistics
        """
        total = len(results)
        successful = sum(1 for r in results if r.get("status") == "success")
        errors = total - successful

        # Calculate average scores
        validation_scores = [
            r["validation"]["overall_score"]
            for r in results
            if r.get("validation") and r["validation"].get("overall_score") is not None
        ]
        avg_validation_score = (
            sum(validation_scores) / len(validation_scores) if validation_scores else None
        )

        ml_scores = [
            r["scoring"]["ml_readiness"]["score"]
            for r in results
            if r.get("scoring")
            and r["scoring"].get("ml_readiness")
            and r["scoring"]["ml_readiness"].get("score") is not None
        ]
        avg_ml_readiness_score = sum(ml_scores) / len(ml_scores) if ml_scores else None

        # Score distribution (90-100, 70-89, 50-69, 0-49)
        excellent = sum(1 for s in validation_scores if s >= 90)
        good = sum(1 for s in validation_scores if 70 <= s < 90)
        moderate = sum(1 for s in validation_scores if 50 <= s < 70)
        poor = sum(1 for s in validation_scores if s < 50)

        # Alert summary
        alert_summary = {}
        for result in results:
            if result.get("alerts") and result["alerts"].get("alerts"):
                for alert in result["alerts"]["alerts"]:
                    catalog = alert.get("catalog", "Unknown")
                    alert_summary[catalog] = alert_summary.get(catalog, 0) + 1

        return {
            "total": total,
            "successful": successful,
            "errors": errors,
            "avg_validation_score": avg_validation_score,
            "avg_ml_readiness_score": avg_ml_readiness_score,
            "score_distribution": {
                "excellent": excellent,
                "good": good,
                "moderate": moderate,
                "poor": poor,
            },
            "alert_summary": alert_summary,
            "processing_time_seconds": None,  # Not available in results
        }

    def _generate_score_distribution_chart(self, results: List[Dict[str, Any]]) -> str:
        """
        Generate SVG bar chart for score distribution.

        Args:
            results: List of batch result dictionaries

        Returns:
            Base64-encoded SVG string
        """
        # Extract validation scores
        scores = [
            r["validation"]["overall_score"]
            for r in results
            if r.get("validation") and r["validation"].get("overall_score") is not None
        ]

        if not scores:
            return ""

        # Count distribution
        excellent = sum(1 for s in scores if s >= 90)
        good = sum(1 for s in scores if 70 <= s < 90)
        moderate = sum(1 for s in scores if 50 <= s < 70)
        poor = sum(1 for s in scores if s < 50)

        total = len(scores)

        # Generate simple SVG bar chart
        svg_width = 600
        svg_height = 300
        bar_width = 120
        bar_spacing = 30
        max_height = 220

        # Calculate bar heights (proportional to count)
        max_count = max(excellent, good, moderate, poor) or 1
        heights = {
            "excellent": int((excellent / max_count) * max_height),
            "good": int((good / max_count) * max_height),
            "moderate": int((moderate / max_count) * max_height),
            "poor": int((poor / max_count) * max_height),
        }

        # Colors matching CSS
        colors = {
            "excellent": "#16a34a",
            "good": "#2563eb",
            "moderate": "#d97706",
            "poor": "#dc2626",
        }

        # Build SVG
        svg_parts = [
            f'<svg width="{svg_width}" height="{svg_height}" xmlns="http://www.w3.org/2000/svg">',
        ]

        # Draw bars
        x_positions = [40, 190, 340, 490]
        labels = ["Excellent\n(90-100)", "Good\n(70-89)", "Moderate\n(50-69)", "Poor\n(0-49)"]
        counts = [excellent, good, moderate, poor]
        categories = ["excellent", "good", "moderate", "poor"]

        for i, (x, label, count, category) in enumerate(
            zip(x_positions, labels, counts, categories)
        ):
            bar_height = heights[category]
            bar_y = max_height - bar_height + 20

            # Draw bar
            svg_parts.append(
                f'<rect x="{x}" y="{bar_y}" width="{bar_width}" height="{bar_height}" '
                f'fill="{colors[category]}" />'
            )

            # Draw count label on bar
            svg_parts.append(
                f'<text x="{x + bar_width/2}" y="{bar_y - 10}" text-anchor="middle" '
                f'font-size="16" font-weight="bold" fill="{colors[category]}">{count}</text>'
            )

            # Draw category label
            label_lines = label.split("\n")
            for j, line in enumerate(label_lines):
                svg_parts.append(
                    f'<text x="{x + bar_width/2}" y="{max_height + 55 + j*18}" '
                    f'text-anchor="middle" font-size="12" fill="#374151">{line}</text>'
                )

        svg_parts.append("</svg>")
        svg_content = "".join(svg_parts)

        # Encode to base64
        svg_bytes = svg_content.encode("utf-8")
        return base64.b64encode(svg_bytes).decode("utf-8")

    def _extract_critical_issues(
        self, results: List[Dict[str, Any]], limit: int = 20
    ) -> List[Dict[str, Any]]:
        """
        Extract critical validation issues from results.

        Args:
            results: List of batch result dictionaries
            limit: Maximum number of issues to return

        Returns:
            List of issue dictionaries with index, smiles, score, issue, severity
        """
        issues = []

        for result in results:
            if not result.get("validation") or not result["validation"].get("issues"):
                continue

            # Look for critical or error severity issues
            for issue in result["validation"]["issues"]:
                if issue.get("severity") in ["CRITICAL", "ERROR"]:
                    issues.append(
                        {
                            "index": result.get("index", 0),
                            "smiles": result.get("smiles", ""),
                            "score": result["validation"].get("overall_score", 0),
                            "issue": issue.get("message", "Unknown issue"),
                            "severity": issue.get("severity", "UNKNOWN"),
                        }
                    )

        # Sort by severity (CRITICAL first) then by score (lowest first)
        severity_order = {"CRITICAL": 0, "ERROR": 1}
        issues.sort(key=lambda x: (severity_order.get(x["severity"], 2), x["score"]))

        return issues[:limit]

    def _get_flagged_molecules(
        self, results: List[Dict[str, Any]], score_threshold: int = 70, limit: int = 30
    ) -> List[Dict[str, Any]]:
        """
        Get molecules with low validation scores for visual inspection.

        Args:
            results: List of batch result dictionaries
            score_threshold: Score below which molecules are flagged
            limit: Maximum number of molecules to return

        Returns:
            List of molecule dictionaries with index, name, score, image_data
        """
        flagged = []

        for result in results:
            if not result.get("validation"):
                continue

            score = result["validation"].get("overall_score", 0)
            if score < score_threshold:
                # Generate molecule image
                image_data = self._mol_to_base64_png(result.get("smiles", ""), width=200, height=200)

                flagged.append(
                    {
                        "index": result.get("index", 0),
                        "name": result.get("name"),
                        "smiles": result.get("smiles", ""),
                        "score": score,
                        "image_data": image_data,
                    }
                )

        # Sort by score (lowest first)
        flagged.sort(key=lambda x: x["score"])

        return flagged[:limit]

    def _mol_to_base64_png(self, smiles: str, width: int = 200, height: int = 200) -> Optional[str]:
        """
        Convert SMILES to base64-encoded PNG image.

        Args:
            smiles: SMILES string
            width: Image width in pixels
            height: Image height in pixels

        Returns:
            Base64-encoded PNG string, or None if conversion fails
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            # Generate image
            img = Draw.MolToImage(mol, size=(width, height))

            # Convert to PNG bytes
            img_buffer = BytesIO()
            img.save(img_buffer, format="PNG")
            img_bytes = img_buffer.getvalue()

            # Encode to base64 (decode to string for HTML embedding)
            return base64.b64encode(img_bytes).decode("utf-8")

        except Exception:
            return None

    def _extract_job_id(self, results: List[Dict[str, Any]]) -> str:
        """
        Extract job ID from results (if available).

        Args:
            results: List of batch result dictionaries

        Returns:
            Job ID string or "unknown"
        """
        # Job ID not stored in individual results, use placeholder
        return "batch_results"


# Register PDF exporter
from app.services.export.base import ExportFormat, ExporterFactory

ExporterFactory.register(ExportFormat.PDF, PDFReportGenerator)
