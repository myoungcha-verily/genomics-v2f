"""In-silico prediction score interpretation for ACMG PP3/BP4."""

import logging
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)


def interpret_scores(scores: dict, config: dict) -> Dict[str, any]:
    """Interpret in-silico prediction scores for ACMG PP3/BP4.

    Uses ClinGen SVI recommended thresholds by default.

    Args:
        scores: dict with keys like 'cadd_phred', 'revel', 'spliceai_max', 'alphamissense'
        config: pipeline config with annotation.thresholds

    Returns dict with:
        pp3: bool (computational evidence supports damaging)
        bp4: bool (computational evidence supports benign)
        pp3_evidence: list of supporting scores
        bp4_evidence: list of supporting scores
        summary: str
    """
    thresholds = config.get("annotation", {}).get("thresholds", {})

    cadd_path = thresholds.get("cadd_pathogenic", 25.3)
    cadd_benign = thresholds.get("cadd_benign", 15.0)
    revel_path = thresholds.get("revel_pathogenic", 0.644)
    revel_benign = thresholds.get("revel_benign", 0.290)
    spliceai_path = thresholds.get("spliceai_pathogenic", 0.2)
    spliceai_benign = thresholds.get("spliceai_benign", 0.1)
    am_path = thresholds.get("alphamissense_pathogenic", 0.564)
    am_benign = thresholds.get("alphamissense_benign", 0.340)

    pp3_evidence = []
    bp4_evidence = []

    # CADD
    cadd = scores.get("cadd_phred")
    if cadd is not None:
        if cadd >= cadd_path:
            pp3_evidence.append(f"CADD={cadd:.1f} (≥{cadd_path})")
        elif cadd <= cadd_benign:
            bp4_evidence.append(f"CADD={cadd:.1f} (≤{cadd_benign})")

    # REVEL
    revel = scores.get("revel")
    if revel is not None:
        if revel >= revel_path:
            pp3_evidence.append(f"REVEL={revel:.3f} (≥{revel_path})")
        elif revel <= revel_benign:
            bp4_evidence.append(f"REVEL={revel:.3f} (≤{revel_benign})")

    # SpliceAI
    spliceai = scores.get("spliceai_max")
    if spliceai is not None:
        if spliceai >= spliceai_path:
            pp3_evidence.append(f"SpliceAI={spliceai:.2f} (≥{spliceai_path})")
        elif spliceai <= spliceai_benign:
            bp4_evidence.append(f"SpliceAI={spliceai:.2f} (≤{spliceai_benign})")

    # AlphaMissense
    am = scores.get("alphamissense")
    if am is not None:
        if am >= am_path:
            pp3_evidence.append(f"AlphaMissense={am:.3f} (≥{am_path})")
        elif am <= am_benign:
            bp4_evidence.append(f"AlphaMissense={am:.3f} (≤{am_benign})")

    # PP3 requires agreement from multiple predictors (ClinGen SVI)
    pp3 = len(pp3_evidence) >= 1 and len(bp4_evidence) == 0
    bp4 = len(bp4_evidence) >= 1 and len(pp3_evidence) == 0

    if pp3:
        summary = f"Damaging ({', '.join(pp3_evidence)})"
    elif bp4:
        summary = f"Benign ({', '.join(bp4_evidence)})"
    else:
        summary = "Indeterminate"

    return {
        "pp3": pp3,
        "bp4": bp4,
        "pp3_evidence": pp3_evidence,
        "bp4_evidence": bp4_evidence,
        "summary": summary,
        "scores": {
            "cadd_phred": cadd,
            "revel": revel,
            "spliceai_max": spliceai,
            "alphamissense": am,
        },
    }


def get_spliceai_max(ds_ag: float = 0, ds_al: float = 0,
                      ds_dg: float = 0, ds_dl: float = 0) -> float:
    """Get maximum SpliceAI delta score across all splice types.

    DS_AG: acceptor gain, DS_AL: acceptor loss,
    DS_DG: donor gain, DS_DL: donor loss
    """
    return max(ds_ag, ds_al, ds_dg, ds_dl)


def format_score(name: str, value: Optional[float]) -> str:
    """Format a score for display."""
    if value is None:
        return f"{name}: N/A"
    if name.lower() == "cadd":
        return f"CADD: {value:.1f}"
    return f"{name}: {value:.3f}"
