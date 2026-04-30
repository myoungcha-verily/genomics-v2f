"""Stage 6: Report Generation

Input: data/classified/acmg_results.parquet, data/phenotype/patient_phenotype.json
Output: reports/<sample>_report.html

Generates per-proband clinical variant reports in HTML format
with tiered variant presentation (P/LP, VUS, LB/B).
"""

import json
import logging
import os
import sys
import time

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipeline.utils.report_renderer import render_proband_report

logger = logging.getLogger(__name__)


def run(config: dict) -> dict:
    """Execute Stage 6: Report Generation."""
    t0 = time.time()
    output_dir = config.get("output", {}).get("output_dir", "data")
    reports_dir = config.get("output", {}).get("reports_dir", "reports")
    class_dir = os.path.join(output_dir, "classified")
    pheno_dir = os.path.join(output_dir, "phenotype")
    os.makedirs(reports_dir, exist_ok=True)

    logger.info("Stage 6: Report Generation")

    # Load classified variants
    class_path = os.path.join(class_dir, "acmg_results.parquet")
    if not os.path.exists(class_path):
        raise FileNotFoundError(f"Stage 4 output not found: {class_path}")

    df = pd.read_parquet(class_path)
    logger.info(f"Loaded {len(df)} classified variants")

    # Load phenotype if available
    pheno_path = os.path.join(pheno_dir, "patient_phenotype.json")
    phenotype = None
    if os.path.exists(pheno_path):
        with open(pheno_path) as f:
            phenotype = json.load(f)

    # Load QC metrics if available
    qc_path = os.path.join(output_dir, "vcf", "qc_report.json")
    qc_metrics = None
    if os.path.exists(qc_path):
        with open(qc_path) as f:
            qc_report = json.load(f)
            qc_metrics = qc_report.get("qc_metrics", {})

    # Generate report per proband
    reports_generated = []

    if "sample_id" in df.columns:
        probands = df["sample_id"].unique()
    else:
        probands = ["unknown"]

    for proband_id in probands:
        proband_df = df[df["sample_id"] == proband_id] if "sample_id" in df.columns else df

        output_path = os.path.join(reports_dir, f"{proband_id}_report.html")
        try:
            report_path = render_proband_report(
                variants_df=proband_df,
                proband_id=proband_id,
                config=config,
                phenotype=phenotype,
                qc_metrics=qc_metrics,
                output_path=output_path,
            )
            reports_generated.append(report_path)
            logger.info(f"Report generated: {report_path}")
        except Exception as e:
            logger.error(f"Report generation failed for {proband_id}: {e}")

    result = {
        "stage": "06_report_generation",
        "reports_generated": len(reports_generated),
        "report_paths": reports_generated,
        "report_format": config.get("output", {}).get("report_format", "html"),
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    with open(os.path.join(reports_dir, "report_summary.json"), "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"\n{'='*60}")
    print(f"Stage 6 Complete: Report Generation")
    print(f"{'='*60}")
    print(f"  Reports generated: {len(reports_generated)}")
    for rp in reports_generated:
        print(f"    → {rp}")
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
