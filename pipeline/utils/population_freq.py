"""Population frequency calculations from gnomAD data."""

import json
import logging
import os
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# Load population map
_POP_MAP = None


def _load_pop_map():
    global _POP_MAP
    if _POP_MAP is not None:
        return _POP_MAP
    ref_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "reference", "gnomad_population_map.json"
    )
    try:
        with open(ref_path) as f:
            _POP_MAP = json.load(f)
    except FileNotFoundError:
        _POP_MAP = {"populations": {}, "thresholds": {}}
    return _POP_MAP


def classify_frequency(af: Optional[float], config: dict) -> Dict[str, any]:
    """Classify variant by population frequency for ACMG criteria.

    Returns dict with:
        af: allele frequency (float or None)
        category: 'common', 'low_frequency', 'rare', 'ultra_rare', 'absent'
        ba1: bool (standalone benign, AF > 5%)
        bs1: bool (strong benign, AF > 1%)
        pm2: bool (supporting pathogenic, AF < 0.01% or absent)
    """
    acmg_cfg = config.get("acmg", {})
    ba1_thresh = acmg_cfg.get("ba1_threshold", 0.05)
    bs1_thresh = acmg_cfg.get("bs1_threshold", 0.01)
    pm2_thresh = acmg_cfg.get("pm2_threshold", 0.0001)

    if af is None or af < 0:
        return {
            "af": None,
            "category": "absent",
            "ba1": False, "bs1": False, "pm2": True,
        }

    ba1 = af > ba1_thresh
    bs1 = af > bs1_thresh
    pm2 = af < pm2_thresh

    if af > ba1_thresh:
        category = "common"
    elif af > bs1_thresh:
        category = "low_frequency"
    elif af > pm2_thresh:
        category = "rare"
    elif af > 0:
        category = "ultra_rare"
    else:
        category = "absent"

    return {
        "af": af,
        "category": category,
        "ba1": ba1, "bs1": bs1, "pm2": pm2,
    }


def get_popmax_af(gnomad_row: dict) -> Tuple[Optional[float], str]:
    """Get maximum population-specific allele frequency.

    Returns (popmax_af, population_name).
    """
    pop_map = _load_pop_map()
    max_af = 0.0
    max_pop = "none"

    for pop_code, pop_info in pop_map.get("populations", {}).items():
        af_key = f"AF_{pop_code}"
        af = gnomad_row.get(af_key)
        if af is not None and af > max_af:
            max_af = af
            max_pop = pop_info["name"]

    return (max_af if max_af > 0 else None, max_pop)


def format_frequency(af: Optional[float]) -> str:
    """Format allele frequency for display."""
    if af is None:
        return "Absent"
    if af == 0:
        return "0"
    if af < 0.0001:
        return f"{af:.2e}"
    if af < 0.01:
        return f"{af:.4f}"
    return f"{af:.3f} ({af*100:.1f}%)"
