"""In-silico prediction score interpretation for ACMG PP3/BP4.

Pejaver/ClinGen 2022 calibrated tiers (PMID 36413997) — REVEL/CADD/SpliceAI/
AlphaMissense each map to BP4_strong → BP4_moderate → BP4_supporting → indeterminate
→ PP3_supporting → PP3_moderate → PP3_strong based on calibrated thresholds.
The ACMG combination engine then weights the strongest tier across predictors
when deciding whether PP3/BP4 fires and at what strength.

Set `annotation.use_calibrated_tiers: false` in pipeline_config.yaml to fall
back to the legacy single-threshold (always-supporting) behavior.
"""

import json
import logging
import os
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Tier rank for "best (most extreme) wins" merging across predictors.
# Lower index = stronger. PP3_strong (0) and BP4_strong (6) bracket the scale.
TIER_RANK = {
    "PP3_strong": 0, "PP3_moderate": 1, "PP3_supporting": 2,
    "indeterminate": 3,
    "BP4_supporting": 4, "BP4_moderate": 5, "BP4_strong": 6,
}

_CALIBRATION_CACHE: Optional[dict] = None


def _load_calibration() -> dict:
    global _CALIBRATION_CACHE
    if _CALIBRATION_CACHE is not None:
        return _CALIBRATION_CACHE
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    cal_path = os.path.join(here, "reference", "insilico_calibration.json")
    if os.path.exists(cal_path):
        with open(cal_path) as f:
            _CALIBRATION_CACHE = json.load(f)
            return _CALIBRATION_CACHE
    logger.warning(f"insilico_calibration.json not found at {cal_path}; "
                   "using legacy single-threshold mode")
    _CALIBRATION_CACHE = {}
    return _CALIBRATION_CACHE


def _classify_tier(score: float, tiers: Dict[str, float], direction: str) -> str:
    """Map a raw score to one of seven tiers using the calibrated thresholds.

    direction='ascending' (REVEL, CADD, SpliceAI, AlphaMissense): higher is more
    pathogenic. We walk from PP3_strong down; first threshold the score meets
    is the one that wins. If no PP3 threshold met, walk BP4_strong upward.
    """
    if direction != "ascending":
        # Currently all four predictors are ascending. Reserved for future.
        return "indeterminate"

    # Pathogenic ladder (highest first)
    for tier in ("PP3_strong", "PP3_moderate", "PP3_supporting"):
        thr = tiers.get(tier)
        if thr is not None and score >= thr:
            return tier

    # Benign ladder (lowest first)
    for tier in ("BP4_strong", "BP4_moderate", "BP4_supporting"):
        thr = tiers.get(tier)
        if thr is not None and score <= thr:
            return tier

    return "indeterminate"


def _best_tier(tiers: List[str]) -> str:
    """Return the strongest (most extreme) tier across predictors."""
    if not tiers:
        return "indeterminate"
    return min(tiers, key=lambda t: TIER_RANK.get(t, 3))


def interpret_scores(scores: dict, config: dict) -> Dict[str, any]:
    """Interpret in-silico prediction scores for ACMG PP3/BP4.

    Returns dict with:
        pp3 (bool):           PP3 fires (any strength)
        bp4 (bool):           BP4 fires (any strength)
        pp3_strength (str):   'supporting' | 'moderate' | 'strong' if pp3
        bp4_strength (str):   'supporting' | 'moderate' | 'strong' if bp4
        best_tier (str):      raw tier label e.g. 'PP3_moderate'
        per_predictor (dict): {predictor: tier} for transparency
        evidence (list):      'REVEL=0.95 → PP3_strong'
        summary (str)
    """
    annot = config.get("annotation", {})
    use_calibrated = annot.get("use_calibrated_tiers", True)
    cal = _load_calibration() if use_calibrated else {}

    score_map = {
        "REVEL": scores.get("revel"),
        "CADD": scores.get("cadd_phred"),
        "SpliceAI": scores.get("spliceai_max"),
        "AlphaMissense": scores.get("alphamissense"),
    }

    per_predictor: Dict[str, str] = {}
    evidence: List[str] = []

    if use_calibrated and cal:
        # Calibrated tiered classification
        for name, val in score_map.items():
            if val is None:
                continue
            tiers = cal.get(name)
            if not tiers:
                continue
            tier = _classify_tier(float(val), tiers, direction="ascending")
            per_predictor[name] = tier
            if tier != "indeterminate":
                fmt = f"{val:.1f}" if name == "CADD" else f"{val:.3f}"
                evidence.append(f"{name}={fmt} → {tier}")
    else:
        # Legacy single-threshold fallback (PP3/BP4 always 'supporting')
        thresholds = annot.get("thresholds", {})
        defaults = {
            "REVEL": (thresholds.get("revel_pathogenic", 0.644),
                       thresholds.get("revel_benign", 0.290)),
            "CADD": (thresholds.get("cadd_pathogenic", 25.3),
                      thresholds.get("cadd_benign", 15.0)),
            "SpliceAI": (thresholds.get("spliceai_pathogenic", 0.2),
                          thresholds.get("spliceai_benign", 0.1)),
            "AlphaMissense": (thresholds.get("alphamissense_pathogenic", 0.564),
                                thresholds.get("alphamissense_benign", 0.340)),
        }
        for name, val in score_map.items():
            if val is None:
                continue
            path_thr, benign_thr = defaults[name]
            if val >= path_thr:
                per_predictor[name] = "PP3_supporting"
                evidence.append(f"{name}={val:.3f} (≥{path_thr})")
            elif val <= benign_thr:
                per_predictor[name] = "BP4_supporting"
                evidence.append(f"{name}={val:.3f} (≤{benign_thr})")

    # Pick the strongest tier across all predictors
    best = _best_tier(list(per_predictor.values()))

    pp3 = best.startswith("PP3_")
    bp4 = best.startswith("BP4_")
    pp3_strength = best.split("_", 1)[1] if pp3 else None
    bp4_strength = best.split("_", 1)[1] if bp4 else None

    if pp3:
        summary = f"Damaging at {pp3_strength} ({', '.join(evidence)})"
    elif bp4:
        summary = f"Benign at {bp4_strength} ({', '.join(evidence)})"
    else:
        summary = "Indeterminate" if not evidence else f"Mixed ({', '.join(evidence)})"

    return {
        "pp3": pp3,
        "bp4": bp4,
        "pp3_strength": pp3_strength,
        "bp4_strength": bp4_strength,
        "best_tier": best,
        "per_predictor": per_predictor,
        "evidence": evidence,
        # Backward-compat fields used by existing acmg_rules code:
        "pp3_evidence": [e for e in evidence if "PP3_" in e or "≥" in e],
        "bp4_evidence": [e for e in evidence if "BP4_" in e or "≤" in e],
        "summary": summary,
        "scores": {
            "cadd_phred": score_map["CADD"],
            "revel": score_map["REVEL"],
            "spliceai_max": score_map["SpliceAI"],
            "alphamissense": score_map["AlphaMissense"],
        },
    }


def get_spliceai_max(ds_ag: float = 0, ds_al: float = 0,
                      ds_dg: float = 0, ds_dl: float = 0) -> float:
    return max(ds_ag, ds_al, ds_dg, ds_dl)


def format_score(name: str, value: Optional[float]) -> str:
    if value is None:
        return f"{name}: N/A"
    if name.lower() == "cadd":
        return f"CADD: {value:.1f}"
    return f"{name}: {value:.3f}"
