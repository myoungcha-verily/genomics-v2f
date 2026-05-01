"""HTML → PDF renderer using weasyprint.

Reuses the existing HTML report — keeps templates consistent. weasyprint
adds page numbers, fixed margins, and a print-friendly layout.

Optional dependency: weasyprint pulls in cairo/pango at the system
level. If weasyprint isn't available at runtime we log a warning and
return None — non-fatal.
"""

from __future__ import annotations

import logging
import os
from typing import Optional

logger = logging.getLogger(__name__)


def render_html_to_pdf(html_path: str, output_path: Optional[str] = None) -> Optional[str]:
    """Convert an HTML file to PDF. Returns the PDF path, or None if
    weasyprint isn't installed or rendering failed."""
    if not os.path.exists(html_path):
        logger.warning(f"PDF render: source HTML not found: {html_path}")
        return None

    try:
        import weasyprint
    except ImportError:
        logger.warning(
            "PDF render: weasyprint not installed (system cairo/pango "
            "dependencies missing?); skipping PDF export"
        )
        return None

    if output_path is None:
        output_path = os.path.splitext(html_path)[0] + ".pdf"

    try:
        weasyprint.HTML(filename=html_path).write_pdf(output_path)
        logger.info(f"PDF rendered: {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"PDF render failed: {e}")
        return None


def is_enabled(config: dict) -> bool:
    out = (config.get("output", {}) or {})
    fmt = str(out.get("report_format", "html")).lower()
    return fmt in ("pdf", "both") or bool(out.get("pdf_export", False))
