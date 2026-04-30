"""Sequence Ontology consequence term mapping and severity ranking."""

import json
import logging
import os
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Load consequence severity from reference file
_SEVERITY_DATA = None


def _load_severity():
    global _SEVERITY_DATA
    if _SEVERITY_DATA is not None:
        return _SEVERITY_DATA

    ref_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "reference", "consequence_severity.json"
    )
    try:
        with open(ref_path) as f:
            _SEVERITY_DATA = json.load(f)
    except FileNotFoundError:
        logger.warning(f"Consequence severity file not found: {ref_path}")
        _SEVERITY_DATA = {"consequences": {}, "severity_order": []}
    return _SEVERITY_DATA


def get_severity(consequence: str) -> str:
    """Get severity level for a consequence term.

    Returns: 'HIGH', 'MODERATE', 'LOW', or 'MODIFIER'
    """
    data = _load_severity()
    info = data["consequences"].get(consequence, {})
    return info.get("severity", "MODIFIER")


def get_rank(consequence: str) -> int:
    """Get numeric rank for a consequence (lower = more severe)."""
    data = _load_severity()
    info = data["consequences"].get(consequence, {})
    return info.get("rank", 999)


def is_lof(consequence: str) -> bool:
    """Check if consequence is loss-of-function."""
    data = _load_severity()
    info = data["consequences"].get(consequence, {})
    return info.get("is_lof", False)


def get_most_severe(consequences: List[str]) -> str:
    """Given a list of consequence terms, return the most severe one."""
    if not consequences:
        return "unknown"
    return min(consequences, key=lambda c: get_rank(c))


def map_vep_consequence(consequence_str: str) -> Dict[str, any]:
    """Map a VEP consequence string to structured info.

    VEP may return multiple consequences separated by '&'.

    Returns dict with:
        consequences: list of individual terms
        most_severe: the most severe term
        severity: severity level of most severe
        is_lof: whether any consequence is LoF
        rank: numeric rank of most severe
    """
    terms = [t.strip() for t in consequence_str.split("&")]
    most_severe = get_most_severe(terms)

    return {
        "consequences": terms,
        "most_severe": most_severe,
        "severity": get_severity(most_severe),
        "is_lof": any(is_lof(t) for t in terms),
        "rank": get_rank(most_severe),
    }


LOF_CONSEQUENCES = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost",
}

MISSENSE_CONSEQUENCES = {
    "missense_variant",
}

INFRAME_CONSEQUENCES = {
    "inframe_insertion",
    "inframe_deletion",
    "stop_lost",
}

SYNONYMOUS_CONSEQUENCES = {
    "synonymous_variant",
}

SPLICE_REGION_CONSEQUENCES = {
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}
