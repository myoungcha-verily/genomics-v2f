#!/usr/bin/env python3
"""ClinVar Concordance Benchmark

Compares pipeline ACMG classifications against ClinVar ground truth
for variants with high-confidence ClinVar classifications (≥2 stars).

Usage:
  python3 scripts/benchmark_clinvar.py
  python3 scripts/benchmark_clinvar.py --input data/classified/acmg_results.parquet
"""

import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))


def benchmark(input_path: str):
    """Run ClinVar concordance benchmark."""
    df = pd.read_parquet(input_path)
    print(f"Loaded {len(df)} classified variants")

    if "clinvar_classification" not in df.columns or \
       "acmg_classification" not in df.columns:
        print("ERROR: Missing classification columns. Run stages 3-4 first.")
        sys.exit(1)

    # Filter to high-confidence ClinVar entries
    high_conf = df[
        df["clinvar_classification"].notna() &
        (df["clinvar_review_stars"] >= 2)
    ].copy()

    if high_conf.empty:
        print("No high-confidence ClinVar variants found.")
        return

    print(f"High-confidence ClinVar entries (≥2 stars): {len(high_conf)}")

    # Simplify classifications
    def simplify(cls):
        if pd.isna(cls):
            return "Unknown"
        cls = str(cls)
        if cls in ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"):
            return "P/LP"
        elif cls in ("Benign", "Likely benign", "Benign/Likely benign"):
            return "B/LB"
        return "VUS"

    high_conf["cv_simple"] = high_conf["clinvar_classification"].apply(simplify)
    high_conf["acmg_simple"] = high_conf["acmg_classification"].apply(simplify)

    # Concordance
    agree = (high_conf["cv_simple"] == high_conf["acmg_simple"]).sum()
    total = len(high_conf)

    # P/LP sensitivity
    cv_plp = high_conf[high_conf["cv_simple"] == "P/LP"]
    plp_sens = (cv_plp["acmg_simple"] == "P/LP").mean() if len(cv_plp) > 0 else 0

    # B/LB specificity
    cv_blb = high_conf[high_conf["cv_simple"] == "B/LB"]
    blb_spec = (cv_blb["acmg_simple"] == "B/LB").mean() if len(cv_blb) > 0 else 0

    print(f"\n{'='*50}")
    print(f"  ClinVar Concordance Benchmark")
    print(f"{'='*50}")
    print(f"  Comparable variants: {total}")
    print(f"  Overall agreement:   {agree/total:.1%}")
    print(f"  P/LP sensitivity:    {plp_sens:.1%} (n={len(cv_plp)})")
    print(f"  B/LB specificity:    {blb_spec:.1%} (n={len(cv_blb)})")

    # Confusion matrix
    print(f"\n  Confusion Matrix (ClinVar vs Pipeline):")
    print(f"  {'':15s} {'P/LP':>8s} {'VUS':>8s} {'B/LB':>8s}")
    for cv_cls in ["P/LP", "VUS", "B/LB"]:
        row = high_conf[high_conf["cv_simple"] == cv_cls]
        counts = [int((row["acmg_simple"] == c).sum()) for c in ["P/LP", "VUS", "B/LB"]]
        print(f"  ClinVar {cv_cls:5s}  {counts[0]:>8d} {counts[1]:>8d} {counts[2]:>8d}")

    print(f"{'='*50}")

    # Save results
    results = {
        "total_comparable": total,
        "overall_agreement": round(agree / total, 4),
        "pathogenic_sensitivity": round(plp_sens, 4),
        "benign_specificity": round(blb_spec, 4),
        "n_clinvar_plp": len(cv_plp),
        "n_clinvar_blb": len(cv_blb),
    }

    output_path = os.path.join("eval", "clinvar_benchmark.json")
    os.makedirs("eval", exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved: {output_path}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/classified/acmg_results.parquet")
    args = parser.parse_args()
    benchmark(args.input)
