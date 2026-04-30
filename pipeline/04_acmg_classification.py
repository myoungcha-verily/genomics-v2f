"""Stage 4: ACMG/AMP Variant Classification

Input: data/enriched/variants_enriched.parquet
Output: data/classified/acmg_results.parquet

Evaluates all 28 ACMG/AMP criteria per variant and assigns
5-tier classification: Pathogenic, Likely pathogenic, VUS,
Likely benign, Benign.
"""

import json
import logging
import os
import sys
import time
from collections import Counter

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipeline.utils.acmg_engine import classify_variant, classification_tier

logger = logging.getLogger(__name__)


def run(config: dict) -> dict:
    """Execute Stage 4: ACMG Classification."""
    t0 = time.time()
    output_dir = config.get("output", {}).get("output_dir", "data")
    enrich_dir = os.path.join(output_dir, "enriched")
    class_dir = os.path.join(output_dir, "classified")
    os.makedirs(class_dir, exist_ok=True)

    logger.info("Stage 4: ACMG/AMP Classification")

    # Load enriched variants
    enrich_path = os.path.join(enrich_dir, "variants_enriched.parquet")
    if not os.path.exists(enrich_path):
        raise FileNotFoundError(f"Stage 3 output not found: {enrich_path}")

    df = pd.read_parquet(enrich_path)
    logger.info(f"Loaded {len(df)} enriched variants")

    # Classify each variant
    classifications = []
    criteria_details = []

    for idx, row in df.iterrows():
        variant_dict = row.to_dict()
        result = classify_variant(variant_dict, config)

        classifications.append({
            "variant_id": row["variant_id"],
            "acmg_classification": result["classification"],
            "acmg_criteria": result["evidence_summary"],
            "acmg_tier": classification_tier(result["classification"]),
            "n_pathogenic_criteria": result["n_pathogenic_criteria"],
            "n_benign_criteria": result["n_benign_criteria"],
        })

        # Store detailed criteria for reporting
        criteria_details.append({
            "variant_id": row["variant_id"],
            "criteria_met": result["criteria_met"],
            "pathogenic_criteria": [
                {"name": n, "strength": s}
                for n, s in result["pathogenic_criteria"]
            ],
            "benign_criteria": [
                {"name": n, "strength": s}
                for n, s in result["benign_criteria"]
            ],
        })

    # Merge classifications with variant data
    class_df = pd.DataFrame(classifications)
    df = df.merge(class_df, on="variant_id", how="left")

    # Sort by tier (most important first), then by classification
    tier_order = {"Pathogenic": 0, "Likely pathogenic": 1, "VUS": 2,
                  "Likely benign": 3, "Benign": 4}
    df["_sort_key"] = df["acmg_classification"].map(tier_order).fillna(5)
    df = df.sort_values("_sort_key").drop(columns=["_sort_key"])

    # Save classified variants
    output_path = os.path.join(class_dir, "acmg_results.parquet")
    df.to_parquet(output_path, index=False)
    logger.info(f"Saved classified variants: {output_path}")

    # Save detailed criteria
    criteria_path = os.path.join(class_dir, "criteria_details.json")
    with open(criteria_path, "w") as f:
        json.dump(criteria_details, f, indent=2, default=str)

    # Classification distribution
    class_counts = Counter(df["acmg_classification"])

    result = {
        "stage": "04_acmg_classification",
        "total_variants": len(df),
        "classification_distribution": dict(class_counts),
        "pathogenic": class_counts.get("Pathogenic", 0),
        "likely_pathogenic": class_counts.get("Likely pathogenic", 0),
        "vus": class_counts.get("VUS", 0),
        "likely_benign": class_counts.get("Likely benign", 0),
        "benign": class_counts.get("Benign", 0),
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    with open(os.path.join(class_dir, "classification_summary.json"), "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"\n{'='*60}")
    print(f"Stage 4 Complete: ACMG/AMP Classification")
    print(f"{'='*60}")
    print(f"  Total variants:    {len(df):,}")
    print(f"  Pathogenic:        {result['pathogenic']:,}")
    print(f"  Likely pathogenic: {result['likely_pathogenic']:,}")
    print(f"  VUS:               {result['vus']:,}")
    print(f"  Likely benign:     {result['likely_benign']:,}")
    print(f"  Benign:            {result['benign']:,}")
    print(f"  Time:              {result['elapsed_seconds']}s")
    print(f"{'='*60}\n")

    return result


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    run(config)
