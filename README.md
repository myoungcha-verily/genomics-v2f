# Variant-to-Function Reporter

Clinical genomics pipeline that takes VCF files to actionable variant interpretation reports. Built for Verily Workbench.

## Key Features

- **7-stage pipeline**: VCF ingest → VEP annotation → ClinVar/gnomAD enrichment → ACMG classification → phenotype integration → HTML reports → validation
- **28 ACMG/AMP criteria**: Automated evaluation with evidence combination rules
- **FHIR phenotype integration**: Uses live patient clinical data to prioritize variants (unique differentiator)
- **Multiple adapters**: Single-sample, trio (de novo detection), gene panel, gVCF
- **Interactive dashboard**: 6-tab Flask app for pipeline monitoring and variant exploration
- **No GPU required**: CPU-only pipeline

## Quick Start

```bash
# 1. Setup environment
bash quickstart.sh

# 2. Configure pipeline
python3 setup_wizard.py

# 3. Validate
python3 run_pipeline.py --dry-run

# 4. Run
python3 run_pipeline.py

# 5. Dashboard
python3 dashboard/app.py
```

## Pipeline Stages

| Stage | Script | Description | Time |
|-------|--------|-------------|------|
| 1 | `01_vcf_ingest_qc.py` | Parse, normalize, QC | 1-5 min |
| 2 | `02_annotate_variants.py` | VEP v114 annotation | 5-30 min |
| 3 | `03_database_enrichment.py` | ClinVar + gnomAD via BQ | 2-10 min |
| 4 | `04_acmg_classification.py` | 28-criteria ACMG engine | 1-5 min |
| 5 | `05_phenotype_integration.py` | FHIR phenotype matching | 1-5 min |
| 6 | `06_report_generation.py` | Per-proband HTML reports | < 1 min |
| 7 | `07_validation.py` | ClinVar concordance | < 1 min |

## Adapters

- **single_sample**: Standard WES/WGS VCF
- **trio**: Multi-sample VCF with proband + parents (de novo detection)
- **panel**: Targeted gene panel (cardiac, cancer, neuro, or custom)
- **gvcf**: Genomic VCF (GATK HaplotypeCaller output)

## Configuration

All settings in `config/pipeline_config.yaml`. Generate with:
```bash
python3 setup_wizard.py
```

See `config/example_configs/` for pre-built configurations.

## Report Output

Self-contained HTML reports with:
- Tier 1: Pathogenic / Likely Pathogenic variants (full evidence cards)
- Tier 2: VUS sorted by in-silico scores
- Tier 3: Benign variants (collapsible)
- ACMG criteria badges, in-silico score bars, ClinVar concordance
- Methodology and limitations section

## Reference Databases

| Database | Access | Refresh |
|----------|--------|---------|
| ClinVar | BQ public dataset | Monthly |
| gnomAD v4.1 | BQ public dataset | Annual |
| VEP cache | Docker/local | Quarterly |
| ClinGen LoF genes | Built-in | As needed |
| Gene panels | Built-in | User-configurable |
