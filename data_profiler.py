#!/usr/bin/env python3
"""VCF Data Profiler — Run after Stage 1 to analyze variant quality.

Generates a quality report with:
- Variant count and type distribution
- Quality metric distributions (DP, GQ, AF)
- Consequence severity breakdown
- Per-chromosome distribution
- Recommendations for pipeline parameters

Usage:
  python3 data_profiler.py
  python3 data_profiler.py --config config/pipeline_config.yaml
"""

import argparse
import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def profile_variants(df: pd.DataFrame) -> dict:
    """Generate comprehensive variant quality profile."""
    n = len(df)
    if n == 0:
        return {"total_variants": 0, "valid": False}

    profile = {
        "total_variants": n,
        "valid": True,
    }

    # Variant types
    is_snv = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    profile["snv_count"] = int(is_snv.sum())
    profile["indel_count"] = int((~is_snv).sum())
    profile["snv_pct"] = round(is_snv.mean() * 100, 1)

    # Genotype distribution
    if "genotype" in df.columns:
        gt_counts = df["genotype"].value_counts().to_dict()
        profile["genotype_distribution"] = {str(k): int(v) for k, v in gt_counts.items()}

    # Quality distributions
    for col, name in [("read_depth", "depth"), ("gt_quality", "gq"),
                       ("allele_fraction", "af")]:
        if col in df.columns:
            vals = df[col].dropna()
            if len(vals) > 0:
                profile[f"{name}_stats"] = {
                    "mean": round(float(vals.mean()), 2),
                    "median": round(float(vals.median()), 2),
                    "std": round(float(vals.std()), 2),
                    "p5": round(float(vals.quantile(0.05)), 2),
                    "p25": round(float(vals.quantile(0.25)), 2),
                    "p75": round(float(vals.quantile(0.75)), 2),
                    "p95": round(float(vals.quantile(0.95)), 2),
                }

    # Consequence distribution (if annotated)
    if "severity" in df.columns:
        sev_counts = df["severity"].value_counts().to_dict()
        profile["severity_distribution"] = {str(k): int(v) for k, v in sev_counts.items()}

    if "consequence" in df.columns:
        cons_counts = df["consequence"].value_counts().head(20).to_dict()
        profile["top_consequences"] = {str(k): int(v) for k, v in cons_counts.items()}

    # Per-chromosome
    if "chrom" in df.columns:
        chrom_counts = df["chrom"].value_counts().to_dict()
        profile["per_chromosome"] = {str(k): int(v) for k, v in chrom_counts.items()}

    # Gene distribution (if annotated)
    if "gene" in df.columns:
        gene_counts = df["gene"].dropna().value_counts().head(20).to_dict()
        profile["top_genes"] = {str(k): int(v) for k, v in gene_counts.items()}
        profile["unique_genes"] = int(df["gene"].dropna().nunique())

    # ACMG classification distribution (if classified)
    if "acmg_classification" in df.columns:
        class_counts = df["acmg_classification"].value_counts().to_dict()
        profile["acmg_distribution"] = {str(k): int(v) for k, v in class_counts.items()}

    return profile


def generate_warnings(profile: dict) -> list:
    """Generate warnings based on profile data."""
    warnings = []

    n = profile.get("total_variants", 0)
    if n == 0:
        return [{"level": "CRITICAL", "message": "No variants found"}]

    if n < 100:
        warnings.append({"level": "WARNING", "message": f"Very few variants ({n}). Check VCF file."})
    elif n > 1000000:
        warnings.append({"level": "INFO", "message": f"Large variant set ({n:,}). Pipeline may take longer."})

    depth_stats = profile.get("depth_stats", {})
    if depth_stats.get("median", 0) < 20:
        warnings.append({"level": "WARNING", "message": f"Low median depth ({depth_stats['median']}x). Quality may be compromised."})

    gq_stats = profile.get("gq_stats", {})
    if gq_stats.get("median", 0) < 30:
        warnings.append({"level": "WARNING", "message": f"Low median GQ ({gq_stats['median']}). Consider stricter filtering."})

    return warnings


def generate_recommendations(profile: dict) -> list:
    """Generate pipeline recommendations."""
    recs = []

    n = profile.get("total_variants", 0)
    if n > 500000:
        recs.append("Large variant set — consider using a gene panel to focus analysis")

    if profile.get("snv_pct", 0) < 80:
        recs.append("High indel proportion — ensure VCF normalization (bcftools norm)")

    if profile.get("depth_stats", {}).get("median", 0) < 30:
        recs.append("Consider WES enrichment QC — median depth below 30x")

    return recs


def main():
    parser = argparse.ArgumentParser(description="VCF Data Profiler")
    parser.add_argument("--config", default="config/pipeline_config.yaml")
    parser.add_argument("--input", default=None, help="Direct path to variants parquet")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s")

    # Find variants file
    if args.input:
        path = args.input
    else:
        with open(args.config) as f:
            config = yaml.safe_load(f)
        output_dir = config.get("output", {}).get("output_dir", "data")

        # Try classified first, then enriched, then raw
        for subdir in ["classified/acmg_results", "enriched/variants_enriched",
                       "vcf/variants"]:
            path = os.path.join(output_dir, f"{subdir}.parquet")
            if os.path.exists(path):
                break
        else:
            print("No variant data found. Run the pipeline first.")
            sys.exit(1)

    print(f"\nProfiling: {path}")
    df = pd.read_parquet(path)
    profile = profile_variants(df)
    warnings = generate_warnings(profile)
    recommendations = generate_recommendations(profile)

    profile["warnings"] = warnings
    profile["recommendations"] = recommendations

    # Save profile
    profile_path = path.replace(".parquet", "_profile.json")
    with open(profile_path, "w") as f:
        json.dump(profile, f, indent=2, default=str)

    # Print summary
    print(f"\n{'='*50}")
    print(f"  Variant Data Profile")
    print(f"{'='*50}")
    print(f"  Total variants:   {profile['total_variants']:,}")
    print(f"  SNVs:             {profile.get('snv_count', 0):,} ({profile.get('snv_pct', 0)}%)")
    print(f"  Indels:           {profile.get('indel_count', 0):,}")

    if "depth_stats" in profile:
        d = profile["depth_stats"]
        print(f"  Depth (median):   {d['median']}x (mean {d['mean']}x)")

    if "unique_genes" in profile:
        print(f"  Unique genes:     {profile['unique_genes']}")

    if "acmg_distribution" in profile:
        print(f"\n  ACMG Classification:")
        for cls, count in profile["acmg_distribution"].items():
            print(f"    {cls:20s}: {count:,}")

    if warnings:
        print(f"\n  Warnings:")
        for w in warnings:
            print(f"    [{w['level']}] {w['message']}")

    if recommendations:
        print(f"\n  Recommendations:")
        for r in recommendations:
            print(f"    → {r}")

    print(f"\n  Profile saved: {profile_path}")
    print(f"{'='*50}\n")


if __name__ == "__main__":
    main()
