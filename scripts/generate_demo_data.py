#!/usr/bin/env python3
"""Generate the bundled demo VCF and FHIR phenotype JSON for V2F.

Variants are synthetic-but-realistic: real ClinVar entries used as
templates (with their canonical pathogenicity), placed at GRCh38
coordinates so downstream stages can join against ClinVar/gnomAD if the
pipeline runs end-to-end. Genotypes are fabricated for demo purposes
only.
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEMO_DIR = REPO_ROOT / "demo" / "germline"

VCF_HEADER = """##fileformat=VCFv4.2
##fileDate=20260430
##source=v2f-demo-generator
##reference=GRCh38
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr17,length=83257441>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEMO_PROBAND
"""

# (chrom, pos, ref, alt, gene, hgvs_c, expected_class, notes)
# Coordinates are GRCh38. Pathogenicity matches ClinVar for the real
# variants; synthetic ones marked.
DEMO_VARIANTS = [
    # 5 curated diagnostic-grade variants
    ("chr17", 43106478, "T", "G", "BRCA1", "c.181T>G",
     "Pathogenic", "BRCA1 Cys61Gly Ashkenazi founder; ClinVar P with 4-star review"),
    ("chr17", 7674220,  "C", "T", "TP53",  "c.524G>A",
     "Pathogenic", "TP53 Arg175His hotspot; classic Li-Fraumeni; ClinVar P"),
    ("chr14", 23427594, "G", "A", "MYH7",  "c.4954G>A",
     "Likely pathogenic", "MYH7 missense; cardiomyopathy panel; synthetic-realistic"),
    ("chr10", 87864458, "C", "T", "PTEN",  "c.493G>A",
     "VUS", "PTEN missense; demonstrates VUS handling"),
    ("chr13", 32340301, "A", "G", "BRCA2", "c.7242A>G",
     "Likely benign", "BRCA2 silent change; common in NFE; ClinVar B/LB"),
    # 5 background neutral variants (common SNPs, intronic, synonymous)
    ("chr17", 43124027, "G", "A", "BRCA1", "c.4837A>G (synonymous)",
     "Benign", "BRCA1 synonymous; demonstrates benign call"),
    ("chr17", 7676154, "G", "C", "TP53",   "intronic",
     "Benign", "TP53 intron 3; common variant"),
    ("chr11", 5226750, "T", "C", "HBB",    "intergenic",
     "Benign", "HBB region; common SNP, sanity check"),
    ("chr10", 87894077, "G", "A", "PTEN",  "3' UTR",
     "Benign", "PTEN UTR variant; demonstrates UTR handling"),
    ("chr14", 23429834, "C", "T", "MYH7",  "intronic",
     "Benign", "MYH7 deep intronic"),
]


def _format_record(chrom, pos, ref, alt, gene, hgvs, _cls, _notes,
                   gt="0/1", dp=42, ad_alt=21, gq=99, af=0.0001):
    info = f"DP={dp};AF={af}"
    fmt = "GT:DP:AD:GQ"
    sample = f"{gt}:{dp}:{dp - ad_alt},{ad_alt}:{gq}"
    return (f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t{info}\t"
            f"{fmt}\t{sample}")


def _chrom_key(chrom):
    s = chrom.replace("chr", "")
    return (0, int(s)) if s.isdigit() else (1, s)


def write_vcf(out_path: Path, variants):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    raw = out_path.with_suffix("")  # .vcf without .gz
    sorted_variants = sorted(variants, key=lambda v: (_chrom_key(v[0]), v[1]))
    with open(raw, "w") as f:
        f.write(VCF_HEADER)
        for v in sorted_variants:
            af = 0.000005 if v[6] in ("Pathogenic", "Likely pathogenic") else \
                 0.0001 if v[6] == "VUS" else 0.05
            f.write(_format_record(*v, af=af) + "\n")

    # bgzip + tabix index for downstream tools
    if shutil.which("bgzip"):
        subprocess.run(["bgzip", "-f", str(raw)], check=True)
        if shutil.which("tabix"):
            subprocess.run(["tabix", "-p", "vcf", str(out_path)], check=True)
    else:
        with open(raw, "rb") as src, gzip.open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        os.unlink(raw)


def write_phenotype(out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "patient_id": "DEMO_PROBAND",
        "source": "Synthetic FHIR Condition payload for V2F demo",
        "conditions": [
            {
                "icd10": "Z80.41",
                "display": "Family history of malignant neoplasm of ovary",
                "snomed": "275937001",
                "onset": "2025-09",
            },
            {
                "icd10": "Z80.3",
                "display": "Family history of malignant neoplasm of breast",
                "snomed": "275938006",
                "onset": "2025-09",
            },
        ],
        "candidate_genes": ["BRCA1", "BRCA2", "TP53", "PTEN", "PALB2",
                            "CHEK2", "ATM", "RAD51C"],
        "hpo_terms": ["HP:0003002", "HP:0100615"],
        "notes": "Hereditary breast and ovarian cancer (HBOC) syndrome workup. "
                 "Phenotype matches BRCA1, BRCA2, TP53, PTEN — drives PP4 "
                 "and phenotype boost for matching variants.",
    }
    out_path.write_text(json.dumps(payload, indent=2))


def write_readme(out_path: Path):
    lines = [
        "# V2F Demo Data",
        "",
        "Bundled synthetic data so a new user can explore the entire pipeline "
        "without bringing their own VCF.",
        "",
        "**Not for clinical use.** Variants are real ClinVar entries used as "
        "templates with fabricated genotypes for demo purposes.",
        "",
        "## Files",
        "",
        "- `germline/proband.vcf.gz` (+ `.tbi`) — 10 variants on chr10/11/13/14/17",
        "- `germline/phenotype.json` — synthetic FHIR Condition payload (HBOC)",
        "- `precomputed/` — pre-baked pipeline outputs so the dashboard tabs "
        "populate before the user runs the pipeline",
        "",
        "## Demo variants — ground truth",
        "",
        "| Gene | HGVS | Expected ACMG | Notes |",
        "|------|------|---------------|-------|",
    ]
    for chrom, pos, ref, alt, gene, hgvs, cls, notes in DEMO_VARIANTS:
        lines.append(f"| {gene} | {hgvs} | {cls} | {notes} |")
    lines += [
        "",
        "## Loading the demo",
        "",
        "From the dashboard:",
        "1. Click **Try with demo data** in the first-run banner, OR",
        "2. Open the **Setup** tab and click **Load demo data**",
        "",
        "From the CLI:",
        "```bash",
        "cp demo/germline/proband.vcf.gz data/vcf/",
        "cp demo/germline/phenotype.json data/phenotype/",
        "cp -r demo/precomputed/* .",
        "```",
        "",
        "## Resetting",
        "",
        "From the dashboard: **Settings** tab → **Reset to first-run state**.",
        "From the CLI: `rm -rf data reports config/pipeline_config.yaml`.",
        "",
        "## Regenerating",
        "",
        "If classification logic changes, regenerate the bundled outputs:",
        "```bash",
        "python3 scripts/generate_demo_data.py",
        "python3 scripts/regenerate_demo_outputs.py  # see Milestone 0.2",
        "```",
    ]
    out_path.write_text("\n".join(lines) + "\n")


def write_demo_config(out_path: Path):
    """A pre-filled pipeline_config.yaml the dashboard installs when the user
    clicks Load demo data. Points at the bundled VCF and phenotype JSON."""
    config = """# Demo configuration for V2F Reporter
# Auto-installed by /api/demo/load — safe to overwrite.
input:
  adapter: single_sample
  vcf_path: data/vcf/proband.vcf.gz
  reference_genome: GRCh38

annotation:
  vep_mode: skip   # demo skips VEP; precomputed parquet has annotations baked in
  thresholds:
    cadd_pathogenic: 25.3
    revel_pathogenic: 0.644
    spliceai_pathogenic: 0.2
    alphamissense_pathogenic: 0.564

databases:
  clinvar:
    bq_table: bigquery-public-data.human_variant_annotation.clinvar_hg38
  gnomad:
    bq_table: bigquery-public-data.gnomAD.v4_1_0_exomes__variant_results

acmg:
  ba1_threshold: 0.05
  bs1_threshold: 0.01
  pm2_threshold: 0.0001
  enable_pvs1: true
  enable_population_criteria: true
  enable_computational: true
  enable_clinical_data: true
  enable_segregation: false
  enable_functional: false

phenotype:
  enabled: true
  fhir_project: ""   # empty → use bundled demo/germline/phenotype.json
  fhir_dataset: ""
  boost_factor: 2.0

output:
  gcs_bucket: ""
  gcs_prefix: genomics-v2f/
  data_dir: data
  reports_dir: reports
  report_format: html

pipeline:
  parallel_workers: 4
  log_level: INFO
"""
    out_path.write_text(config)


def main():
    DEMO_DIR.mkdir(parents=True, exist_ok=True)
    write_vcf(DEMO_DIR / "proband.vcf.gz", DEMO_VARIANTS)
    write_phenotype(DEMO_DIR / "phenotype.json")
    write_demo_config(DEMO_DIR / "config.yaml")
    write_readme(REPO_ROOT / "demo" / "README.md")
    print(f"Wrote {DEMO_DIR / 'proband.vcf.gz'}")
    print(f"Wrote {DEMO_DIR / 'phenotype.json'}")
    print(f"Wrote {DEMO_DIR / 'config.yaml'}")
    print(f"Wrote {REPO_ROOT / 'demo' / 'README.md'}")


if __name__ == "__main__":
    sys.exit(main() or 0)
