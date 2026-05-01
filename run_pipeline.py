#!/usr/bin/env python3
"""Variant-to-Function Reporter — Pipeline Orchestrator

Single entry point to run all 7 stages:
  1. VCF Ingest & QC
  2. Variant Annotation (VEP)
  3. Database Enrichment (ClinVar, gnomAD)
  4. ACMG/AMP Classification
  5. Phenotype Integration (FHIR)
  6. Report Generation
  7. Validation & Benchmarking

Usage:
  python3 run_pipeline.py                                    # Full pipeline
  python3 run_pipeline.py --config config/pipeline_config.yaml
  python3 run_pipeline.py --dry-run                          # Validate only
  python3 run_pipeline.py --start-stage 3                    # Resume from stage
  python3 run_pipeline.py --stages 1,4,6                     # Run specific stages
  python3 run_pipeline.py --no-gcs-sync                      # Skip GCS upload
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import time
from datetime import datetime

import yaml

# Stage definitions
STAGES = {
    1: {
        "name": "VCF Ingest & QC",
        "module": "pipeline.01_vcf_ingest_qc",
        "outputs": ["data/vcf/variants.parquet", "data/vcf/qc_report.json"],
        "prerequisites": [],
        "estimated_time": "1-5 min",
    },
    2: {
        "name": "Variant Annotation",
        "module": "pipeline.02_annotate_variants",
        "outputs": ["data/annotated/annotated_variants.parquet"],
        "prerequisites": ["data/vcf/variants.parquet"],
        "estimated_time": "5-30 min",
    },
    3: {
        "name": "Database Enrichment",
        "module": "pipeline.03_database_enrichment",
        "outputs": ["data/enriched/variants_enriched.parquet"],
        "prerequisites": ["data/annotated/annotated_variants.parquet"],
        "estimated_time": "2-10 min",
    },
    4: {
        "name": "ACMG Classification",
        "module": "pipeline.04_acmg_classification",
        "outputs": ["data/classified/acmg_results.parquet"],
        "prerequisites": ["data/enriched/variants_enriched.parquet"],
        "estimated_time": "1-5 min",
    },
    5: {
        "name": "Phenotype Integration",
        "module": "pipeline.05_phenotype_integration",
        "outputs": ["data/phenotype/patient_phenotype.json"],
        "prerequisites": ["data/classified/acmg_results.parquet"],
        "estimated_time": "1-5 min",
    },
    6: {
        "name": "Report Generation",
        "module": "pipeline.06_report_generation",
        "outputs": ["reports/"],
        "prerequisites": ["data/classified/acmg_results.parquet"],
        "estimated_time": "< 1 min",
    },
    7: {
        "name": "Validation",
        "module": "pipeline.07_validation",
        "outputs": ["eval/validation_summary.json"],
        "prerequisites": ["data/classified/acmg_results.parquet"],
        "estimated_time": "< 1 min",
    },
}


def load_config(config_path: str) -> dict:
    """Load and validate pipeline configuration."""
    if not os.path.exists(config_path):
        print(f"ERROR: Config file not found: {config_path}")
        print("Run setup_wizard.py first to generate configuration.")
        sys.exit(1)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    return config


def validate_config(config: dict) -> list:
    """Validate configuration, return list of issues."""
    issues = []

    # Check input
    input_cfg = config.get("input", {})
    if not input_cfg.get("vcf_path"):
        issues.append("input.vcf_path is empty — specify a VCF file path")

    adapter = input_cfg.get("adapter", "single_sample")
    if adapter not in ("single_sample", "trio", "panel", "gvcf"):
        issues.append(f"Unknown adapter: {adapter}")

    if adapter == "trio":
        trio = input_cfg.get("trio", {})
        if not trio.get("proband_id") and not trio.get("ped_file"):
            issues.append("Trio adapter requires proband_id or ped_file")

    # Check annotation
    ann = config.get("annotation", {})
    vep_mode = ann.get("vep_mode", "docker")
    if vep_mode == "docker":
        # Check if Docker is available
        result = subprocess.run(["docker", "--version"],
                                capture_output=True, text=True)
        if result.returncode != 0:
            issues.append("Docker not found — set annotation.vep_mode to 'skip' or 'local'")

    return issues


def check_prerequisites(stage_num: int) -> list:
    """Check if prerequisite files exist for a stage."""
    stage = STAGES[stage_num]
    missing = []
    for prereq in stage["prerequisites"]:
        if not os.path.exists(prereq):
            missing.append(prereq)
    return missing


def run_stage(stage_num: int, config: dict) -> dict:
    """Run a single pipeline stage."""
    stage = STAGES[stage_num]
    module_name = stage["module"]

    print(f"\n{'━'*60}")
    print(f"  Stage {stage_num}/7: {stage['name']}")
    print(f"  Estimated time: {stage['estimated_time']}")
    print(f"{'━'*60}")

    # Import and run the stage module
    try:
        module = __import__(module_name, fromlist=["run"])
        result = module.run(config)
        return result
    except Exception as e:
        logger.error(f"Stage {stage_num} failed: {e}")
        suggestion = _get_error_suggestion(stage_num, str(e))
        if suggestion:
            print(f"\n  SUGGESTION: {suggestion}")
        raise


def sync_to_gcs(config: dict, stage_num: int):
    """Sync stage outputs to GCS bucket."""
    gcs_bucket = config.get("output", {}).get("gcs_bucket", "")
    gcs_prefix = config.get("output", {}).get("gcs_prefix", "genomics-v2f/")

    if not gcs_bucket:
        return

    stage = STAGES[stage_num]
    for output in stage["outputs"]:
        if not os.path.exists(output):
            continue

        gcs_path = f"{gcs_bucket.rstrip('/')}/{gcs_prefix.rstrip('/')}/{output}"
        try:
            if os.path.isdir(output):
                cmd = ["gsutil", "-m", "rsync", "-r", output, gcs_path]
            else:
                cmd = ["gsutil", "cp", output, gcs_path]

            subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            logger.info(f"Synced to GCS: {output} → {gcs_path}")
        except Exception as e:
            logger.warning(f"GCS sync failed: {e}")


def dry_run(config: dict, stages_to_run: list):
    """Validate config and prerequisites without running."""
    print("\n" + "="*60)
    print("  DRY RUN — Validation Only")
    print("="*60)

    # Config validation
    issues = validate_config(config)
    if issues:
        print("\nConfiguration issues:")
        for issue in issues:
            print(f"  ✗ {issue}")
    else:
        print("\n  ✓ Configuration valid")

    # Check prerequisites
    print(f"\nStages to run: {stages_to_run}")
    for stage_num in stages_to_run:
        stage = STAGES[stage_num]
        missing = check_prerequisites(stage_num)
        status = "✓" if not missing else "✗"
        print(f"  {status} Stage {stage_num}: {stage['name']}")
        if missing:
            for m in missing:
                print(f"      Missing: {m}")

    # Check VCF file
    vcf_path = config.get("input", {}).get("vcf_path", "")
    if vcf_path:
        if vcf_path.startswith("gs://"):
            print(f"\n  VCF: {vcf_path} (GCS — will download at runtime)")
        elif os.path.exists(vcf_path):
            size_mb = os.path.getsize(vcf_path) / 1024 / 1024
            print(f"\n  ✓ VCF: {vcf_path} ({size_mb:.1f} MB)")
        else:
            print(f"\n  ✗ VCF not found: {vcf_path}")

    # Check tools
    for tool in ["bcftools", "docker"]:
        result = subprocess.run([tool, "--version"],
                                capture_output=True, text=True)
        status = "✓" if result.returncode == 0 else "✗ (optional)"
        version = result.stdout.split("\n")[0] if result.returncode == 0 else "not found"
        print(f"  {status} {tool}: {version}")

    print("\n" + "="*60)


def _get_error_suggestion(stage_num: int, error: str) -> str:
    """Get context-aware error suggestion."""
    error_lower = error.lower()

    if "not found" in error_lower and "vcf" in error_lower:
        return "Check input.vcf_path in your config"
    if "docker" in error_lower:
        return "Docker not available. Set annotation.vep_mode: skip"
    if "bigquery" in error_lower or "403" in error_lower:
        return "BigQuery access denied. Check permissions or run: gcloud auth application-default login"
    if "cyvcf2" in error_lower or "pysam" in error_lower:
        return "VCF parser not installed. Run: pip install cyvcf2 pysam"
    if stage_num == 3 and "timeout" in error_lower:
        return "BQ query timed out. Try reducing batch size or check network"

    return ""


def main():
    parser = argparse.ArgumentParser(description="Variant-to-Function Reporter Pipeline")
    parser.add_argument("--config", default="config/pipeline_config.yaml",
                        help="Path to pipeline config YAML")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate config and prerequisites without running")
    parser.add_argument("--start-stage", type=int, default=1,
                        help="Start from a specific stage (1-7)")
    parser.add_argument("--stages", type=str, default=None,
                        help="Comma-separated list of stages to run (e.g., '1,4,6')")
    parser.add_argument("--no-gcs-sync", action="store_true",
                        help="Skip GCS sync after each stage")
    parser.add_argument("--log-level", default=None,
                        help="Logging level (DEBUG, INFO, WARNING)")
    args = parser.parse_args()

    # Change to project directory
    project_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(project_dir)
    sys.path.insert(0, project_dir)

    # Load config
    config = load_config(args.config)

    # Setup logging
    log_level = args.log_level or config.get("pipeline", {}).get("log_level", "INFO")
    log_file = config.get("pipeline", {}).get("log_file", "logs/pipeline.log")
    os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)

    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file),
        ],
    )
    global logger
    logger = logging.getLogger("pipeline")

    # Determine stages to run
    if args.stages:
        stages_to_run = [int(s.strip()) for s in args.stages.split(",")]
    else:
        all_stages = list(range(args.start_stage, 8))
        stages_cfg = config.get("pipeline", {}).get("stages", "all")
        if stages_cfg != "all" and isinstance(stages_cfg, str):
            stages_to_run = [int(s.strip()) for s in stages_cfg.split(",")]
        else:
            stages_to_run = all_stages

    # Dry run?
    if args.dry_run:
        dry_run(config, stages_to_run)
        return

    # Print banner
    print("\n" + "╔" + "═"*58 + "╗")
    print("║" + " Variant-to-Function Reporter".center(58) + "║")
    print("║" + " Clinical Genomics Pipeline v1.0".center(58) + "║")
    print("╚" + "═"*58 + "╝")

    t_start = time.time()
    results = {}
    failed_stage = None

    # Record run start in run history (for the dashboard's Runs tab)
    try:
        from pipeline.utils.run_manifest import build_manifest
        from pipeline.utils.run_history import record_run_start
        manifest = build_manifest(config)
        run_id = manifest["run_id"]
        record_run_start(run_id, config, manifest)
    except Exception as e:
        print(f"  (run history record_start failed, continuing: {e})")
        manifest = None
        run_id = "unknown"

    for stage_num in stages_to_run:
        if stage_num not in STAGES:
            print(f"Unknown stage: {stage_num}")
            continue

        # Check prerequisites (skip for first stage or if resuming)
        missing = check_prerequisites(stage_num)
        if missing and stage_num != stages_to_run[0]:
            print(f"\nSkipping stage {stage_num}: missing prerequisites")
            for m in missing:
                print(f"  → {m}")
            continue

        try:
            result = run_stage(stage_num, config)
            results[stage_num] = result

            # GCS sync
            if not args.no_gcs_sync:
                sync_to_gcs(config, stage_num)

        except KeyboardInterrupt:
            print(f"\n\nInterrupted at stage {stage_num}")
            failed_stage = stage_num
            break
        except Exception as e:
            print(f"\n  ERROR in stage {stage_num}: {e}")
            failed_stage = stage_num
            break

    # Summary
    elapsed = time.time() - t_start
    print("\n" + "╔" + "═"*58 + "╗")
    print("║" + " Pipeline Summary".center(58) + "║")
    print("╚" + "═"*58 + "╝")

    for stage_num in stages_to_run:
        if stage_num in results:
            print(f"  ✓ Stage {stage_num}: {STAGES[stage_num]['name']}")
        elif stage_num == failed_stage:
            print(f"  ✗ Stage {stage_num}: {STAGES[stage_num]['name']} (FAILED)")
        else:
            print(f"  ○ Stage {stage_num}: {STAGES[stage_num]['name']} (skipped)")

    print(f"\n  Total time: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    # Record run end in run history
    try:
        from pipeline.utils.run_history import record_run_end
        status = "failure" if failed_stage else (
            "partial" if any(s not in results for s in stages_to_run) else "success")
        summary = {"completed_stages": list(results.keys()),
                    "failed_stage": failed_stage,
                    "elapsed_seconds": round(elapsed, 1)}
        record_run_end(run_id, status, summary, config)
    except Exception as e:
        print(f"  (run history record_end failed: {e})")

    if failed_stage:
        print(f"\n  Pipeline failed at stage {failed_stage}.")
        print(f"  Fix the issue and resume: python3 run_pipeline.py --start-stage {failed_stage}")
    else:
        print(f"\n  Pipeline completed successfully!")
        # Check for reports
        reports_dir = config.get("output", {}).get("reports_dir", "reports")
        if os.path.exists(reports_dir):
            reports = [f for f in os.listdir(reports_dir) if f.endswith(".html")]
            if reports:
                print(f"  Reports: {reports_dir}/")
                for r in reports:
                    print(f"    → {r}")

    # Save run log
    run_log = {
        "timestamp": datetime.now().isoformat(),
        "config_path": args.config,
        "stages_run": list(results.keys()),
        "failed_stage": failed_stage,
        "elapsed_seconds": round(elapsed, 1),
        "results": {str(k): v for k, v in results.items()},
    }
    log_path = os.path.join("logs", "run_log.json")
    os.makedirs("logs", exist_ok=True)
    with open(log_path, "w") as f:
        json.dump(run_log, f, indent=2, default=str)

    print()


if __name__ == "__main__":
    main()
