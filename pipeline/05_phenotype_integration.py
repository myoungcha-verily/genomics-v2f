"""Stage 5: Phenotype Integration (Differentiator)

Input: data/classified/acmg_results.parquet + FHIR data
Output: data/phenotype/patient_phenotype.json, updated variant scores

This is the key differentiator vs competitors — using live FHIR
clinical data (conditions, observations) to prioritize variants
that match the patient's phenotype.
"""

import json
import logging
import os
import sys
import time

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipeline.utils.fhir_phenotype import (
    extract_patient_phenotype,
    score_variant_phenotype_match,
)

logger = logging.getLogger(__name__)


def run(config: dict) -> dict:
    """Execute Stage 5: Phenotype Integration."""
    t0 = time.time()
    output_dir = config.get("output", {}).get("output_dir", "data")
    class_dir = os.path.join(output_dir, "classified")
    pheno_dir = os.path.join(output_dir, "phenotype")
    os.makedirs(pheno_dir, exist_ok=True)

    pheno_config = config.get("phenotype", {})
    enabled = pheno_config.get("enabled", False)

    logger.info(f"Stage 5: Phenotype Integration (enabled={enabled})")

    # Load classified variants
    class_path = os.path.join(class_dir, "acmg_results.parquet")
    if not os.path.exists(class_path):
        raise FileNotFoundError(f"Stage 4 output not found: {class_path}")

    df = pd.read_parquet(class_path)
    logger.info(f"Loaded {len(df)} classified variants")

    # Extract patient phenotype from FHIR
    proband_id = df["sample_id"].iloc[0] if "sample_id" in df.columns else "unknown"

    if enabled:
        phenotype = extract_patient_phenotype(proband_id, config)
        logger.info(f"Phenotype: {phenotype['n_conditions']} conditions, "
                     f"{phenotype['n_candidate_genes']} candidate genes")

        # Score each variant against phenotype
        boost = pheno_config.get("boost_phenotype_match", 2.0)
        scores = []
        for _, row in df.iterrows():
            score = score_variant_phenotype_match(row.to_dict(), phenotype)
            scores.append(score)

        df["phenotype_match_score"] = scores

        # Boost phenotype-matching variants
        n_boosted = (df["phenotype_match_score"] > 0).sum()
        logger.info(f"Phenotype-matching variants: {n_boosted}")

        # Re-sort: phenotype-matching VUS get promoted
        df["_priority"] = df.apply(
            lambda r: _compute_priority(r, boost), axis=1
        )
        df = df.sort_values("_priority").drop(columns=["_priority"])

    else:
        phenotype = {
            "patient_id": proband_id,
            "conditions": [],
            "candidate_genes": [],
            "n_conditions": 0,
            "n_candidate_genes": 0,
            "note": "Phenotype integration disabled in config",
        }
        df["phenotype_match_score"] = 0.0

    # Save phenotype profile
    pheno_path = os.path.join(pheno_dir, "patient_phenotype.json")
    with open(pheno_path, "w") as f:
        json.dump(phenotype, f, indent=2, default=str)

    # Save updated variants (with phenotype scores)
    updated_path = os.path.join(class_dir, "acmg_results.parquet")
    df.to_parquet(updated_path, index=False)

    n_matching = int((df["phenotype_match_score"] > 0).sum()) if enabled else 0

    result = {
        "stage": "05_phenotype_integration",
        "enabled": enabled,
        "n_conditions": phenotype.get("n_conditions", 0),
        "n_candidate_genes": phenotype.get("n_candidate_genes", 0),
        "n_phenotype_matching_variants": n_matching,
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    with open(os.path.join(pheno_dir, "phenotype_summary.json"), "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"\n{'='*60}")
    print(f"Stage 5 Complete: Phenotype Integration")
    print(f"{'='*60}")
    print(f"  Enabled:             {enabled}")
    print(f"  Patient conditions:  {phenotype.get('n_conditions', 0)}")
    print(f"  Candidate genes:     {phenotype.get('n_candidate_genes', 0)}")
    print(f"  Matching variants:   {n_matching}")
    print(f"  Time:                {result['elapsed_seconds']}s")
    print(f"{'='*60}\n")

    return result


def _compute_priority(row, boost: float) -> float:
    """Compute sort priority (lower = higher priority)."""
    tier_map = {
        "Pathogenic": 0,
        "Likely pathogenic": 1,
        "VUS": 2,
        "Likely benign": 3,
        "Benign": 4,
    }
    tier = tier_map.get(row.get("acmg_classification", "VUS"), 2)

    # Phenotype-matching VUS get boosted
    pheno_score = row.get("phenotype_match_score", 0)
    if tier == 2 and pheno_score > 0:
        tier = 1.5  # Between LP and VUS

    return tier


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    run(config)
