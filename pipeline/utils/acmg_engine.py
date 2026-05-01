"""ACMG/AMP variant classification engine.

Evaluates all 28 ACMG criteria and combines evidence
to assign a 5-tier classification:
  Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign
"""

import json
import logging
import os
from typing import Dict, List, Tuple

from pipeline.utils.acmg_rules import ALL_CRITERIA

logger = logging.getLogger(__name__)


def classify_variant(variant: dict, config: dict) -> Dict:
    """Classify a single variant using ACMG/AMP criteria.

    Args:
        variant: dict with variant annotations (consequence, scores, frequencies)
        config: pipeline config

    Returns dict with:
        classification: str (Pathogenic, Likely pathogenic, VUS, Likely benign, Benign)
        criteria_met: list of (criterion_name, strength, evidence)
        pathogenic_criteria: list of triggered pathogenic criteria
        benign_criteria: list of triggered benign criteria
        evidence_summary: str
    """
    # Apply gene-specific VCEP overrides (no-op if no rule file for this gene)
    from pipeline.utils.vcep_loader import merge_vcep_overrides
    config = merge_vcep_overrides(variant, config)

    acmg_cfg = config.get("acmg", {})

    # Evaluate all criteria
    criteria_met = []
    pathogenic_criteria = []
    benign_criteria = []

    for name, eval_fn in ALL_CRITERIA.items():
        # Skip disabled criteria
        if name.startswith("PVS") and not acmg_cfg.get("enable_pvs1", True):
            continue
        if name in ("BA1", "BS1", "PM2") and not acmg_cfg.get(
                "enable_population_criteria", True):
            continue
        if name in ("PP3", "BP4", "BP7") and not acmg_cfg.get(
                "enable_computational", True):
            continue
        if name in ("PS1", "PM5", "PP5", "BP6") and not acmg_cfg.get(
                "enable_clinical_data", True):
            continue
        if name in ("PS2", "PM6", "PM3") and not acmg_cfg.get(
                "enable_segregation", True):
            continue
        if name in ("PS3", "BS3") and not acmg_cfg.get(
                "enable_functional", False):
            continue

        triggered, strength, evidence = eval_fn(variant, config)
        if triggered:
            criteria_met.append({
                "criterion": name,
                "strength": strength,
                "evidence": evidence,
            })
            if name.startswith(("PVS", "PS", "PM", "PP")):
                pathogenic_criteria.append((name, strength))
            else:
                benign_criteria.append((name, strength))

    # Combine evidence to classify
    classification = _combine_evidence(pathogenic_criteria, benign_criteria)

    # Build evidence summary
    criteria_names = [c["criterion"] for c in criteria_met]
    summary = ", ".join(criteria_names) if criteria_names else "No criteria met"

    return {
        "classification": classification,
        "criteria_met": criteria_met,
        "pathogenic_criteria": pathogenic_criteria,
        "benign_criteria": benign_criteria,
        "evidence_summary": summary,
        "n_pathogenic_criteria": len(pathogenic_criteria),
        "n_benign_criteria": len(benign_criteria),
    }


def _combine_evidence(pathogenic: List[Tuple[str, str]],
                       benign: List[Tuple[str, str]]) -> str:
    """Combine pathogenic and benign evidence per ACMG rules.

    Returns classification string.
    """
    # Count by strength
    path_counts = _count_by_strength(pathogenic)
    ben_counts = _count_by_strength(benign)

    # Check benign first (BA1 is standalone)
    if ben_counts.get("standalone", 0) >= 1:
        return "Benign"
    if ben_counts.get("strong", 0) >= 2:
        return "Benign"

    # Likely benign
    if ben_counts.get("strong", 0) >= 1 and ben_counts.get("supporting", 0) >= 1:
        return "Likely benign"
    if ben_counts.get("supporting", 0) >= 2:
        return "Likely benign"

    # Pathogenic
    pvs = path_counts.get("very_strong", 0)
    ps = path_counts.get("strong", 0)
    pm = path_counts.get("moderate", 0)
    pp = path_counts.get("supporting", 0)

    if _is_pathogenic(pvs, ps, pm, pp):
        return "Pathogenic"
    if _is_likely_pathogenic(pvs, ps, pm, pp):
        return "Likely pathogenic"

    # If both pathogenic and benign evidence, still VUS
    return "VUS"


def _is_pathogenic(pvs: int, ps: int, pm: int, pp: int) -> bool:
    """Check if evidence meets Pathogenic threshold."""
    if pvs >= 1 and ps >= 1:
        return True
    if pvs >= 1 and pm >= 2:
        return True
    if pvs >= 1 and pm >= 1 and pp >= 1:
        return True
    if pvs >= 1 and pp >= 2:
        return True
    if ps >= 2:
        return True
    if ps >= 1 and pm >= 3:
        return True
    if ps >= 1 and pm >= 2 and pp >= 2:
        return True
    if ps >= 1 and pm >= 1 and pp >= 4:
        return True
    return False


def _is_likely_pathogenic(pvs: int, ps: int, pm: int, pp: int) -> bool:
    """Check if evidence meets Likely Pathogenic threshold."""
    if pvs >= 1 and pm >= 1:
        return True
    if ps >= 1 and pm >= 1:
        return True
    if ps >= 1 and pp >= 2:
        return True
    if pm >= 3:
        return True
    if pm >= 2 and pp >= 2:
        return True
    if pm >= 1 and pp >= 4:
        return True
    return False


def _count_by_strength(criteria: List[Tuple[str, str]]) -> Dict[str, int]:
    """Count criteria by evidence strength."""
    counts = {}
    for _, strength in criteria:
        counts[strength] = counts.get(strength, 0) + 1
    return counts


def classification_tier(classification: str) -> int:
    """Map classification to display tier (1=most important)."""
    tiers = {
        "Pathogenic": 1,
        "Likely pathogenic": 1,
        "VUS": 2,
        "Likely benign": 3,
        "Benign": 3,
    }
    return tiers.get(classification, 2)


def classification_color(classification: str) -> str:
    """Map classification to display color."""
    colors = {
        "Pathogenic": "#d32f2f",
        "Likely pathogenic": "#e64a19",
        "VUS": "#f9a825",
        "Likely benign": "#388e3c",
        "Benign": "#1b5e20",
    }
    return colors.get(classification, "#757575")
