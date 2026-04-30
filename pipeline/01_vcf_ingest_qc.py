"""Stage 1: VCF Ingest & Quality Control

Input: VCF file (local or GCS)
Output: data/vcf/normalized.vcf.gz, data/vcf/variants.parquet, data/vcf/qc_report.json

Steps:
1. Resolve VCF path (download from GCS if needed)
2. Normalize with bcftools (multi-allelic decomposition, left-alignment)
3. Parse variants using adapter
4. Compute QC metrics
5. Apply basic quality filters
6. Save normalized variants and QC report
"""

import json
import logging
import os
import sys
import time

import pandas as pd
import yaml

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from adapters import get_adapter
from pipeline.utils.vcf_parser import normalize_vcf, compute_qc_metrics, filter_variants

logger = logging.getLogger(__name__)


def run(config: dict) -> dict:
    """Execute Stage 1: VCF Ingest & QC.

    Returns dict with stage results and metrics.
    """
    t0 = time.time()
    input_cfg = config.get("input", {})
    output_dir = config.get("output", {}).get("output_dir", "data")
    vcf_dir = os.path.join(output_dir, "vcf")
    os.makedirs(vcf_dir, exist_ok=True)

    adapter_name = input_cfg.get("adapter", "single_sample")
    logger.info(f"Stage 1: VCF Ingest & QC (adapter={adapter_name})")

    # Step 1: Initialize adapter
    adapter = get_adapter(adapter_name, config)
    proband_id = adapter.get_proband_id()
    sample_ids = adapter.get_sample_ids()
    logger.info(f"Proband: {proband_id}, Samples: {sample_ids}")

    # Step 2: Validate VCF
    validation = adapter.validate_vcf()
    if not validation["valid"]:
        logger.error(f"VCF validation failed: {validation['issues']}")
        raise ValueError(f"VCF validation failed: {validation['issues']}")

    # Step 3: Normalize VCF (if bcftools available)
    vcf_path = adapter._resolve_vcf_path()
    normalized_vcf = os.path.join(vcf_dir, "normalized.vcf.gz")
    try:
        normalized_vcf = normalize_vcf(vcf_path, normalized_vcf)
        logger.info(f"Normalized VCF: {normalized_vcf}")
    except FileNotFoundError:
        logger.warning("bcftools not found — skipping normalization")
        normalized_vcf = vcf_path

    # Step 4: Parse variants
    variants_df = adapter.load_variants()
    if variants_df.empty:
        raise ValueError("No variants found in VCF")
    logger.info(f"Parsed {len(variants_df)} variants")

    # Step 5: QC metrics
    qc_metrics = compute_qc_metrics(variants_df)
    logger.info(f"QC: Ti/Tv={qc_metrics.get('ti_tv_ratio', 'N/A')}, "
                f"Het/Hom={qc_metrics.get('het_hom_ratio', 'N/A')}, "
                f"Mean DP={qc_metrics.get('mean_depth', 'N/A')}x")

    if qc_metrics.get("issues"):
        for issue in qc_metrics["issues"]:
            logger.warning(f"QC Warning: {issue}")

    # Step 6: Basic quality filter
    filtered_df, filter_counts = filter_variants(
        variants_df,
        min_gq=20.0,
        min_dp=10,
        min_af=0.15,
        pass_only=True,
    )
    logger.info(f"After quality filter: {len(filtered_df)} variants")

    # Step 7: Save outputs
    variants_path = os.path.join(vcf_dir, "variants.parquet")
    filtered_df.to_parquet(variants_path, index=False)
    logger.info(f"Saved variants: {variants_path}")

    # Also save unfiltered for reference
    all_variants_path = os.path.join(vcf_dir, "variants_all.parquet")
    variants_df.to_parquet(all_variants_path, index=False)

    # Save QC report
    qc_report = {
        "stage": "01_vcf_ingest_qc",
        "adapter": adapter_name,
        "proband_id": proband_id,
        "sample_ids": sample_ids,
        "reference_genome": input_cfg.get("reference_genome", "GRCh38"),
        "vcf_path": input_cfg.get("vcf_path", ""),
        "normalized_vcf": normalized_vcf,
        "total_variants_raw": len(variants_df),
        "total_variants_filtered": len(filtered_df),
        "filter_counts": filter_counts,
        "qc_metrics": qc_metrics,
        "elapsed_seconds": round(time.time() - t0, 1),
    }
    qc_path = os.path.join(vcf_dir, "qc_report.json")
    with open(qc_path, "w") as f:
        json.dump(qc_report, f, indent=2, default=str)
    logger.info(f"QC report: {qc_path}")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Stage 1 Complete: VCF Ingest & QC")
    print(f"{'='*60}")
    print(f"  Adapter:      {adapter_name}")
    print(f"  Proband:      {proband_id}")
    print(f"  Raw variants: {len(variants_df):,}")
    print(f"  After filter: {len(filtered_df):,}")
    print(f"  Ti/Tv ratio:  {qc_metrics.get('ti_tv_ratio', 'N/A')}")
    print(f"  Het/Hom:      {qc_metrics.get('het_hom_ratio', 'N/A')}")
    print(f"  Mean depth:   {qc_metrics.get('mean_depth', 'N/A')}x")
    if qc_metrics.get("issues"):
        print(f"  Warnings:     {len(qc_metrics['issues'])}")
    print(f"  Time:         {qc_report['elapsed_seconds']}s")
    print(f"{'='*60}\n")

    return qc_report


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    run(config)
