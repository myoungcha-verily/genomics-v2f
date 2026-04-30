#!/usr/bin/env python3
"""Variant-to-Function Reporter — Interactive Setup Wizard

Walks through 6 steps to generate config/pipeline_config.yaml:
  1. VCF Input Type (single_sample / trio / panel / gvcf)
  2. VCF Source (local path or GCS URI; trio: sample IDs / PED file)
  3. Reference Genome (GRCh37 or GRCh38)
  4. Annotation Sources (VEP mode, in-silico predictors, databases)
  5. Phenotype Integration (FHIR project/dataset)
  6. Output & Storage (GCS bucket, report format)

Usage:
  python3 setup_wizard.py
  python3 setup_wizard.py --output config/pipeline_config.yaml
"""

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime

import yaml


def ask(prompt: str, default: str = "", options: list = None) -> str:
    """Ask user a question with optional default and options."""
    if options:
        print(f"\n  Options: {', '.join(options)}")
    suffix = f" [{default}]" if default else ""
    while True:
        answer = input(f"  {prompt}{suffix}: ").strip()
        if not answer:
            return default
        if options and answer not in options:
            print(f"  Invalid choice. Options: {', '.join(options)}")
            continue
        return answer


def ask_yn(prompt: str, default: bool = True) -> bool:
    """Ask yes/no question."""
    suffix = " [Y/n]" if default else " [y/N]"
    answer = input(f"  {prompt}{suffix}: ").strip().lower()
    if not answer:
        return default
    return answer in ("y", "yes")


def validate_bq_table(table: str) -> bool:
    """Check if BigQuery table is accessible."""
    try:
        result = subprocess.run(
            ["bq", "show", "--format=json", table],
            capture_output=True, text=True, timeout=15,
        )
        return result.returncode == 0
    except Exception:
        return False


def run_wizard(output_path: str):
    """Run the interactive setup wizard."""
    print("\n" + "="*60)
    print("  Variant-to-Function Reporter — Setup Wizard")
    print("="*60)
    print("\n  This wizard will generate your pipeline configuration.")
    print("  Press Enter to accept [defaults].\n")

    config = build_default_config()

    # ========== Step 1: VCF Input Type ==========
    print("\n" + "─"*60)
    print("  STEP 1/6: VCF Input Type")
    print("─"*60)
    print("  • single_sample: Standard single-sample VCF/VCF.gz")
    print("  • trio: Multi-sample VCF with proband + parents")
    print("  • panel: Targeted gene panel (filters to panel genes)")
    print("  • gvcf: Genomic VCF (GATK HaplotypeCaller output)")

    adapter = ask("Input type", "single_sample",
                  ["single_sample", "trio", "panel", "gvcf"])
    config["input"]["adapter"] = adapter

    # ========== Step 2: VCF Source ==========
    print("\n" + "─"*60)
    print("  STEP 2/6: VCF Source")
    print("─"*60)

    vcf_path = ask("VCF file path (local or gs://...)")
    config["input"]["vcf_path"] = vcf_path

    if vcf_path and not vcf_path.startswith("gs://") and not os.path.exists(vcf_path):
        print(f"  WARNING: File not found: {vcf_path}")

    # Trio-specific
    if adapter == "trio":
        print("\n  Trio configuration:")
        proband = ask("Proband sample ID")
        father = ask("Father sample ID (or empty)")
        mother = ask("Mother sample ID (or empty)")
        ped = ask("PED file path (or empty)")
        config["input"]["trio"]["proband_id"] = proband
        config["input"]["trio"]["father_id"] = father
        config["input"]["trio"]["mother_id"] = mother
        config["input"]["trio"]["ped_file"] = ped

    # Panel-specific
    if adapter == "panel":
        print("\n  Available panels: cardiac, cancer, neuro, custom_template")
        panel = ask("Gene panel name or path", "cardiac")
        padding = ask("Padding (bp)", "20")
        config["input"]["panel"]["gene_panel"] = panel
        config["input"]["panel"]["padding_bp"] = int(padding)

    # ========== Step 3: Reference Genome ==========
    print("\n" + "─"*60)
    print("  STEP 3/6: Reference Genome")
    print("─"*60)

    genome = ask("Reference genome", "GRCh38", ["GRCh37", "GRCh38"])
    config["input"]["reference_genome"] = genome
    config["annotation"]["vep_assembly"] = genome

    # ========== Step 4: Annotation Sources ==========
    print("\n" + "─"*60)
    print("  STEP 4/6: Annotation Sources")
    print("─"*60)

    # VEP mode
    print("  VEP modes:")
    print("  • docker: Run VEP in Docker container (recommended)")
    print("  • local: Use locally installed VEP")
    print("  • skip: Skip VEP (use basic consequence prediction)")
    vep_mode = ask("VEP mode", "skip", ["docker", "local", "skip"])
    config["annotation"]["vep_mode"] = vep_mode

    if vep_mode == "docker":
        cache_dir = ask("VEP cache directory (empty for --database mode)")
        config["annotation"]["vep_cache_dir"] = cache_dir

    # In-silico predictors
    print("\n  In-silico predictors (used for ACMG PP3/BP4):")
    config["annotation"]["enable_cadd"] = ask_yn("Enable CADD?", True)
    config["annotation"]["enable_revel"] = ask_yn("Enable REVEL?", True)
    config["annotation"]["enable_spliceai"] = ask_yn("Enable SpliceAI?", True)

    # Databases
    print("\n  Reference databases (via BigQuery):")
    clinvar_table = config["databases"]["clinvar"]["bq_table"]
    use_clinvar = ask_yn(f"Query ClinVar ({clinvar_table})?", True)
    if not use_clinvar:
        config["databases"]["clinvar"]["bq_table"] = ""

    gnomad_table = config["databases"]["gnomad"]["bq_table"]
    use_gnomad = ask_yn(f"Query gnomAD ({gnomad_table})?", True)
    if not use_gnomad:
        config["databases"]["gnomad"]["bq_table"] = ""

    # ========== Step 5: Phenotype Integration ==========
    print("\n" + "─"*60)
    print("  STEP 5/6: Phenotype Integration (FHIR)")
    print("─"*60)
    print("  This is the key differentiator — using patient clinical data")
    print("  to prioritize variants that match their phenotype.")

    enable_pheno = ask_yn("Enable FHIR phenotype integration?", False)
    config["phenotype"]["enabled"] = enable_pheno

    if enable_pheno:
        fhir_project = ask("FHIR BigQuery project")
        fhir_dataset = ask("FHIR BigQuery dataset")
        config["phenotype"]["fhir_project"] = fhir_project
        config["phenotype"]["fhir_dataset"] = fhir_dataset

        if fhir_project and fhir_dataset:
            table = f"{fhir_project}.{fhir_dataset}.Condition"
            if validate_bq_table(table):
                print(f"  ✓ FHIR Condition table accessible")
            else:
                print(f"  ⚠ Cannot verify {table} — check permissions")

    # ========== Step 6: Output & Storage ==========
    print("\n" + "─"*60)
    print("  STEP 6/6: Output & Storage")
    print("─"*60)

    gcs_bucket = ask("GCS bucket for output persistence (empty to skip)")
    config["output"]["gcs_bucket"] = gcs_bucket

    report_format = ask("Report format", "html", ["html"])
    config["output"]["report_format"] = report_format

    include_vus = ask_yn("Include VUS in report?", True)
    config["output"]["include_vus"] = include_vus

    if include_vus:
        max_vus = ask("Max VUS in report", "50")
        config["output"]["max_vus_in_report"] = int(max_vus)

    # ========== Write Config ==========
    write_config(config, output_path)

    print("\n" + "="*60)
    print("  Setup Complete!")
    print("="*60)
    print(f"\n  Config saved: {output_path}")
    print(f"\n  Next steps:")
    print(f"    1. Review config:  cat {output_path}")
    print(f"    2. Validate:       python3 run_pipeline.py --dry-run")
    print(f"    3. Run pipeline:   python3 run_pipeline.py")
    print()


def build_default_config() -> dict:
    """Build default configuration dict."""
    return {
        "input": {
            "adapter": "single_sample",
            "vcf_path": "",
            "reference_genome": "GRCh38",
            "trio": {
                "proband_id": "",
                "father_id": "",
                "mother_id": "",
                "ped_file": "",
            },
            "panel": {
                "gene_panel": "",
                "padding_bp": 20,
            },
        },
        "annotation": {
            "vep_mode": "skip",
            "vep_docker_image": "ensemblorg/ensembl-vep:release_114.0",
            "vep_cache_dir": "",
            "vep_species": "homo_sapiens",
            "vep_assembly": "GRCh38",
            "enable_cadd": True,
            "enable_revel": True,
            "enable_spliceai": True,
            "enable_alphamissense": True,
            "thresholds": {
                "cadd_pathogenic": 25.3,
                "cadd_benign": 15.0,
                "revel_pathogenic": 0.644,
                "revel_benign": 0.290,
                "spliceai_pathogenic": 0.2,
                "spliceai_benign": 0.1,
                "alphamissense_pathogenic": 0.564,
                "alphamissense_benign": 0.340,
            },
        },
        "databases": {
            "clinvar": {
                "bq_table": "bigquery-public-data.human_variant_annotation.clinvar_hg38",
                "min_review_stars": 1,
            },
            "gnomad": {
                "bq_table": "bigquery-public-data.gnomAD.v4_1_0_exomes__variant_results",
                "bq_table_fallback": "",
            },
            "dbnsfp": {
                "enabled": False,
                "bq_table": "",
            },
            "clingen": {
                "gene_list_path": "reference/clingen_lof_genes.json",
            },
        },
        "acmg": {
            "ba1_threshold": 0.05,
            "bs1_threshold": 0.01,
            "pm2_threshold": 0.0001,
            "enable_pvs1": True,
            "enable_population_criteria": True,
            "enable_computational": True,
            "enable_clinical_data": True,
            "enable_segregation": True,
            "enable_functional": False,
        },
        "phenotype": {
            "enabled": False,
            "fhir_project": "",
            "fhir_dataset": "",
            "condition_table": "Condition",
            "observation_table": "Observation",
            "omim_mapping": True,
            "hpo_matching": True,
            "boost_phenotype_match": 2.0,
        },
        "output": {
            "output_dir": "data",
            "reports_dir": "reports",
            "gcs_bucket": "",
            "gcs_prefix": "genomics-v2f/",
            "report_format": "html",
            "include_vus": True,
            "max_vus_in_report": 50,
            "include_benign": False,
        },
        "pipeline": {
            "stages": "all",
            "resume_from": None,
            "max_workers": 4,
            "log_level": "INFO",
            "log_file": "logs/pipeline.log",
        },
    }


def write_config(config: dict, output_path: str):
    """Write config to YAML, backing up existing file."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    if os.path.exists(output_path):
        backup = output_path + f".bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        shutil.copy2(output_path, backup)
        print(f"\n  Backed up existing config: {backup}")

    with open(output_path, "w") as f:
        f.write(f"# Generated by setup_wizard.py on {datetime.now().isoformat()}\n")
        f.write("# Variant-to-Function Reporter Pipeline Configuration\n\n")
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Setup Wizard")
    parser.add_argument("--output", default="config/pipeline_config.yaml",
                        help="Output config path")
    args = parser.parse_args()

    run_wizard(args.output)
