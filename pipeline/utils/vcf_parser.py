"""VCF parsing utilities — normalization, QC metrics, and helpers."""

import json
import logging
import os
import subprocess
from typing import Dict, Any, Optional, Tuple

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def normalize_vcf(input_vcf: str, output_vcf: str,
                  reference_fasta: Optional[str] = None) -> str:
    """Normalize VCF using bcftools norm.

    Steps:
    1. Decompose multi-allelic sites (-m -both)
    2. Left-align indels (if reference provided)
    3. Remove exact duplicates (-d exact)
    4. Index output

    Returns path to normalized VCF.
    """
    os.makedirs(os.path.dirname(output_vcf) or ".", exist_ok=True)

    # Step 1: Decompose multi-allelic
    cmd = ["bcftools", "norm", "-m", "-both"]
    if reference_fasta:
        cmd.extend(["-f", reference_fasta])
    cmd.extend(["-O", "z", "-o", output_vcf, input_vcf])

    logger.info(f"Normalizing: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.warning(f"bcftools norm failed: {result.stderr}")
        return input_vcf

    # Step 2: Remove duplicates
    dedup_vcf = output_vcf.replace(".vcf.gz", ".dedup.vcf.gz")
    cmd2 = ["bcftools", "norm", "-d", "exact",
            "-O", "z", "-o", dedup_vcf, output_vcf]
    result2 = subprocess.run(cmd2, capture_output=True, text=True)
    if result2.returncode == 0:
        os.replace(dedup_vcf, output_vcf)
    else:
        try:
            os.remove(dedup_vcf)
        except OSError:
            pass

    # Step 3: Index
    subprocess.run(["bcftools", "index", "-t", output_vcf],
                    capture_output=True, text=True)

    # Count variants
    count_result = subprocess.run(
        ["bcftools", "view", "-H", output_vcf],
        capture_output=True, text=True
    )
    n_variants = count_result.stdout.count("\n") if count_result.returncode == 0 else 0
    logger.info(f"Normalized VCF: {n_variants} variants")

    return output_vcf


def compute_qc_metrics(variants_df: pd.DataFrame) -> Dict[str, Any]:
    """Compute quality control metrics from parsed variants.

    Returns dict with:
    - total_variants: int
    - snv_count, indel_count: int
    - ti_tv_ratio: float (transition/transversion ratio; WES ~2.8-3.2, WGS ~2.0-2.1)
    - het_hom_ratio: float (typically ~1.5-2.0 for outbred populations)
    - mean_depth: float
    - mean_gq: float
    - pass_rate: float (fraction with FILTER=PASS)
    - singleton_rate: float
    """
    if variants_df.empty:
        return {"total_variants": 0, "valid": False, "issues": ["No variants"]}

    n = len(variants_df)
    metrics = {"total_variants": n}

    # SNV vs indel
    is_snv = (variants_df["ref"].str.len() == 1) & (variants_df["alt"].str.len() == 1)
    metrics["snv_count"] = int(is_snv.sum())
    metrics["indel_count"] = int((~is_snv).sum())

    # Ti/Tv ratio
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    snvs = variants_df[is_snv]
    if len(snvs) > 0:
        is_ti = snvs.apply(
            lambda r: (r["ref"], r["alt"]) in transitions, axis=1
        )
        n_ti = is_ti.sum()
        n_tv = len(snvs) - n_ti
        metrics["ti_tv_ratio"] = round(n_ti / n_tv, 3) if n_tv > 0 else float("inf")
        metrics["n_transitions"] = int(n_ti)
        metrics["n_transversions"] = int(n_tv)
    else:
        metrics["ti_tv_ratio"] = None

    # Het/Hom ratio
    is_het = variants_df["genotype"].isin(["0/1", "0|1", "1/0", "1|0"])
    is_hom_alt = variants_df["genotype"].isin(["1/1", "1|1"])
    n_het = is_het.sum()
    n_hom = is_hom_alt.sum()
    metrics["het_count"] = int(n_het)
    metrics["hom_alt_count"] = int(n_hom)
    metrics["het_hom_ratio"] = round(n_het / n_hom, 3) if n_hom > 0 else float("inf")

    # Quality metrics
    metrics["mean_depth"] = round(float(variants_df["read_depth"].mean()), 1)
    metrics["median_depth"] = round(float(variants_df["read_depth"].median()), 1)
    metrics["mean_gq"] = round(float(variants_df["gt_quality"].mean()), 1)
    metrics["mean_allele_fraction"] = round(
        float(variants_df["allele_fraction"].mean()), 3
    )

    # PASS rate
    pass_mask = variants_df["filter"] == "PASS"
    metrics["pass_rate"] = round(float(pass_mask.mean()), 4)
    metrics["pass_count"] = int(pass_mask.sum())
    metrics["filtered_count"] = int((~pass_mask).sum())

    # Per-chromosome distribution
    chrom_counts = variants_df["chrom"].value_counts().to_dict()
    metrics["per_chromosome"] = {str(k): int(v) for k, v in chrom_counts.items()}

    # Quality warnings
    issues = []
    if metrics.get("ti_tv_ratio") and metrics["ti_tv_ratio"] < 1.5:
        issues.append(f"Low Ti/Tv ratio ({metrics['ti_tv_ratio']}), expected >2.0 for WES")
    if metrics["mean_depth"] < 20:
        issues.append(f"Low mean depth ({metrics['mean_depth']}x), recommend >30x")
    if metrics["pass_rate"] < 0.8:
        issues.append(f"Low PASS rate ({metrics['pass_rate']:.1%})")
    if metrics["het_hom_ratio"] > 5.0:
        issues.append(f"High Het/Hom ratio ({metrics['het_hom_ratio']}), possible contamination")

    metrics["issues"] = issues
    metrics["valid"] = len([i for i in issues if "Low" in i]) < 2

    return metrics


def filter_variants(variants_df: pd.DataFrame,
                    min_gq: float = 20.0,
                    min_dp: int = 10,
                    min_af: float = 0.15,
                    pass_only: bool = True) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Apply quality filters to variants.

    Returns (filtered_df, filter_counts).
    """
    n_start = len(variants_df)
    filter_counts = {"input": n_start}

    df = variants_df.copy()

    if pass_only:
        mask = df["filter"] == "PASS"
        filter_counts["failed_filter"] = int((~mask).sum())
        df = df[mask]

    mask_gq = df["gt_quality"] >= min_gq
    filter_counts["low_gq"] = int((~mask_gq).sum())
    df = df[mask_gq]

    mask_dp = df["read_depth"] >= min_dp
    filter_counts["low_dp"] = int((~mask_dp).sum())
    df = df[mask_dp]

    mask_af = df["allele_fraction"] >= min_af
    filter_counts["low_af"] = int((~mask_af).sum())
    df = df[mask_af]

    filter_counts["output"] = len(df)
    filter_counts["removed"] = n_start - len(df)

    logger.info(f"Filtered: {n_start} → {len(df)} variants "
                 f"({filter_counts['removed']} removed)")
    return df, filter_counts
