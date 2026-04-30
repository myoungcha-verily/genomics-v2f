"""Stage 7: Validation & Benchmarking

Input: data/classified/acmg_results.parquet
Output: eval/validation_summary.json

Benchmarks ACMG classification accuracy against ClinVar
ground truth (variants with ≥2 review stars).
"""

import json
import logging
import os
import sys
import time
from collections import Counter

import pandas as pd
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

logger = logging.getLogger(__name__)


def run(config: dict) -> dict:
    """Execute Stage 7: Validation."""
    t0 = time.time()
    output_dir = config.get("output", {}).get("output_dir", "data")
    class_dir = os.path.join(output_dir, "classified")
    eval_dir = "eval"
    os.makedirs(eval_dir, exist_ok=True)

    logger.info("Stage 7: Validation & Benchmarking")

    # Load classified variants
    class_path = os.path.join(class_dir, "acmg_results.parquet")
    if not os.path.exists(class_path):
        raise FileNotFoundError(f"Stage 4 output not found: {class_path}")

    df = pd.read_parquet(class_path)
    logger.info(f"Loaded {len(df)} classified variants")

    # ClinVar concordance analysis
    concordance = _clinvar_concordance(df)

    # Classification distribution
    if "acmg_classification" in df.columns:
        class_dist = df["acmg_classification"].value_counts().to_dict()
    else:
        class_dist = {}

    # Criteria usage analysis
    criteria_usage = _criteria_usage_analysis(df)

    # Quality metrics summary
    quality_summary = _quality_summary(df)

    result = {
        "stage": "07_validation",
        "total_variants": len(df),
        "classification_distribution": class_dist,
        "clinvar_concordance": concordance,
        "criteria_usage": criteria_usage,
        "quality_summary": quality_summary,
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    output_path = os.path.join(eval_dir, "validation_summary.json")
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"\n{'='*60}")
    print(f"Stage 7 Complete: Validation & Benchmarking")
    print(f"{'='*60}")
    print(f"  Total variants: {len(df):,}")
    print(f"  Classification distribution:")
    for cls, count in sorted(class_dist.items()):
        print(f"    {cls:20s}: {count:,}")
    if concordance.get("n_comparable"):
        print(f"\n  ClinVar Concordance ({concordance['n_comparable']} variants):")
        print(f"    Overall agreement: {concordance.get('overall_agreement', 0):.1%}")
        print(f"    P/LP sensitivity:  {concordance.get('pathogenic_sensitivity', 0):.1%}")
        print(f"    B/LB specificity:  {concordance.get('benign_specificity', 0):.1%}")
    else:
        print(f"\n  ClinVar concordance: No comparable variants")
    print(f"  Time: {result['elapsed_seconds']}s")
    print(f"{'='*60}\n")

    return result


def _clinvar_concordance(df: pd.DataFrame) -> dict:
    """Compute concordance between our ACMG classification and ClinVar.

    Only considers ClinVar entries with ≥2 review stars.
    """
    if "clinvar_classification" not in df.columns or \
       "acmg_classification" not in df.columns:
        return {"n_comparable": 0, "note": "Missing classification columns"}

    # Filter to ClinVar variants with high confidence
    min_stars = 2
    has_clinvar = (
        df["clinvar_classification"].notna() &
        (df["clinvar_review_stars"] >= min_stars)
    )
    comparable = df[has_clinvar].copy()

    if comparable.empty:
        return {"n_comparable": 0, "note": f"No ClinVar variants with ≥{min_stars} stars"}

    n = len(comparable)

    # Map to simplified categories
    def simplify(cls):
        if cls in ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"):
            return "P/LP"
        elif cls in ("Benign", "Likely benign", "Benign/Likely benign"):
            return "B/LB"
        return "VUS"

    comparable["clinvar_simple"] = comparable["clinvar_classification"].apply(simplify)
    comparable["acmg_simple"] = comparable["acmg_classification"].apply(simplify)

    # Overall agreement
    agree = (comparable["clinvar_simple"] == comparable["acmg_simple"]).sum()

    # P/LP sensitivity: Of ClinVar P/LP, how many did we call P/LP?
    clinvar_plp = comparable[comparable["clinvar_simple"] == "P/LP"]
    n_plp = len(clinvar_plp)
    plp_correct = (clinvar_plp["acmg_simple"] == "P/LP").sum() if n_plp > 0 else 0

    # B/LB specificity: Of ClinVar B/LB, how many did we call B/LB?
    clinvar_blb = comparable[comparable["clinvar_simple"] == "B/LB"]
    n_blb = len(clinvar_blb)
    blb_correct = (clinvar_blb["acmg_simple"] == "B/LB").sum() if n_blb > 0 else 0

    # Confusion matrix
    confusion = {}
    for clinvar_cls in ["P/LP", "VUS", "B/LB"]:
        row = {}
        for acmg_cls in ["P/LP", "VUS", "B/LB"]:
            row[acmg_cls] = int(
                ((comparable["clinvar_simple"] == clinvar_cls) &
                 (comparable["acmg_simple"] == acmg_cls)).sum()
            )
        confusion[clinvar_cls] = row

    return {
        "n_comparable": n,
        "min_review_stars": min_stars,
        "overall_agreement": round(agree / n, 4) if n > 0 else 0,
        "pathogenic_sensitivity": round(plp_correct / n_plp, 4) if n_plp > 0 else None,
        "benign_specificity": round(blb_correct / n_blb, 4) if n_blb > 0 else None,
        "n_clinvar_pathogenic": n_plp,
        "n_clinvar_benign": n_blb,
        "n_clinvar_vus": len(comparable[comparable["clinvar_simple"] == "VUS"]),
        "confusion_matrix": confusion,
    }


def _criteria_usage_analysis(df: pd.DataFrame) -> dict:
    """Analyze which ACMG criteria are most frequently triggered."""
    if "acmg_criteria" not in df.columns:
        return {}

    all_criteria = []
    for criteria_str in df["acmg_criteria"].dropna():
        for c in str(criteria_str).split(", "):
            c = c.strip()
            if c and c != "No criteria met":
                all_criteria.append(c)

    counts = Counter(all_criteria)
    return {
        "total_criteria_triggered": len(all_criteria),
        "unique_criteria": len(counts),
        "most_common": dict(counts.most_common(15)),
    }


def _quality_summary(df: pd.DataFrame) -> dict:
    """Summarize variant quality metrics."""
    summary = {}
    if "read_depth" in df.columns:
        summary["mean_depth"] = round(float(df["read_depth"].mean()), 1)
        summary["median_depth"] = round(float(df["read_depth"].median()), 1)
    if "gt_quality" in df.columns:
        summary["mean_gq"] = round(float(df["gt_quality"].mean()), 1)
    if "gnomad_af" in df.columns:
        n_with_freq = df["gnomad_af"].notna().sum()
        summary["variants_with_gnomad"] = int(n_with_freq)
        summary["pct_with_gnomad"] = round(n_with_freq / len(df) * 100, 1) if len(df) > 0 else 0
    if "clinvar_classification" in df.columns:
        n_with_clinvar = df["clinvar_classification"].notna().sum()
        summary["variants_with_clinvar"] = int(n_with_clinvar)
    return summary


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    run(config)
