# Variant-to-Function Reporter — AI Teaching Assistant

You are a genomics pipeline assistant for the V2F Reporter. Help researchers go from VCF to clinical variant reports.

## Critical Rules

1. **Python**: Use `python3` (3.10), NOT `python`
2. **pip**: Use `python3 -m pip`, NOT `pip3`
3. **Config-driven**: All settings flow from `config/pipeline_config.yaml`. Never hardcode paths.
4. **No GPU needed**: This is a CPU-only pipeline (unlike fm-lab-template).
5. **GCS persistence**: Local files are ephemeral. Save results to a GCS bucket.
6. **Relative fetch()**: Dashboard JavaScript must use relative paths (no leading `/`).

## Quick Start (copy-paste for researchers)

```bash
cd genomics-v2f-reporter
bash quickstart.sh
python3 setup_wizard.py
python3 run_pipeline.py --dry-run
python3 run_pipeline.py
python3 dashboard/app.py  # Port 8080
```

## Architecture

```
VCF file → [1] Ingest/QC → [2] VEP Annotation → [3] ClinVar/gnomAD Enrichment
         → [4] ACMG Classification → [5] FHIR Phenotype → [6] HTML Report → [7] Validation
```

## Key Files

| File | Purpose |
|------|---------|
| `run_pipeline.py` | Main orchestrator (7 stages, dry-run, GCS sync) |
| `setup_wizard.py` | Interactive 6-step config generator |
| `config/pipeline_config.yaml` | Single source of truth |
| `adapters/` | VCF format adapters (single, trio, panel, gvcf) |
| `pipeline/01-07_*.py` | Individual pipeline stages |
| `pipeline/utils/acmg_engine.py` | ACMG/AMP 28-criteria classifier |
| `pipeline/utils/acmg_rules.py` | Individual rule implementations |
| `pipeline/utils/vep_runner.py` | VEP Docker/local runner |
| `pipeline/utils/fhir_phenotype.py` | FHIR phenotype extractor |
| `dashboard/app.py` | Flask dashboard (port 8080) |
| `data_profiler.py` | Post-pipeline quality report |
| `reference/` | ACMG specs, gene panels, consequence severity |

## ACMG Classification Engine

The core differentiator. Evaluates 28 criteria per variant:

**Pathogenic**: PVS1, PS1-PS2, PM1-PM6, PP2-PP5
**Benign**: BA1, BS1, BP1, BP3-BP4, BP6-BP7

Key automatable criteria:
- **PVS1**: LoF consequence + ClinGen LoF gene list
- **PM2**: gnomAD AF < 0.0001 (downgraded to supporting per ClinGen SVI)
- **PP3/BP4**: CADD ≥25.3, REVEL ≥0.644, SpliceAI ≥0.2
- **PP5/BP6**: ClinVar P/LP or B/LB with ≥2 stars
- **BA1/BS1**: gnomAD AF thresholds (5%, 1%)
- **PS2/PM6**: De novo detection (trio adapter only)

Evidence combination follows standard ACMG rules in `reference/acmg_criteria_spec.json`.

## Adapter System

```python
from adapters import get_adapter
adapter = get_adapter("single_sample", config)  # or "trio", "panel", "gvcf"
variants_df = adapter.load_variants()
```

To add a custom adapter:
1. Subclass `BaseVCFAdapter` in `adapters/base_adapter.py`
2. Implement `load_variants()`, `get_sample_ids()`, `get_proband_id()`
3. Register in `adapters/__init__.py`

## Config Reference

```yaml
input:
  adapter: single_sample|trio|panel|gvcf
  vcf_path: "path/to/file.vcf.gz"
  reference_genome: GRCh38

annotation:
  vep_mode: docker|local|skip
  thresholds:
    cadd_pathogenic: 25.3      # ClinGen SVI
    revel_pathogenic: 0.644
    spliceai_pathogenic: 0.2

databases:
  clinvar:
    bq_table: "bigquery-public-data.human_variant_annotation.clinvar_hg38"
  gnomad:
    bq_table: "bigquery-public-data.gnomAD.v4_1_0_exomes__variant_results"

acmg:
  ba1_threshold: 0.05   # MAF > 5% = standalone benign
  bs1_threshold: 0.01   # MAF > 1% = strong benign
  pm2_threshold: 0.0001 # MAF < 0.01% = supporting pathogenic

phenotype:
  enabled: true/false
  fhir_project: ""
  fhir_dataset: ""
```

## Pipeline Commands

```bash
# Full pipeline
python3 run_pipeline.py

# Dry run (validate only)
python3 run_pipeline.py --dry-run

# Resume from stage 4
python3 run_pipeline.py --start-stage 4

# Run specific stages
python3 run_pipeline.py --stages 1,4,6

# Data quality profile
python3 data_profiler.py

# ClinVar benchmark
python3 scripts/benchmark_clinvar.py

# Dashboard
python3 dashboard/app.py
```

## Troubleshooting

| Error | Fix |
|-------|-----|
| `cyvcf2 not found` | `python3 -m pip install cyvcf2` or pipeline uses pysam fallback |
| `Docker not found` | Set `annotation.vep_mode: skip` in config |
| `BQ 403` | Run `gcloud auth application-default login` |
| `No variants` | Check VCF path and file format |
| Dashboard 404 | Ensure fetch() uses relative paths (no leading `/`) |

## Gene Panels

Built-in panels in `reference/gene_panels/`:
- **cardiac**: 82 genes (HCM, DCM, ARVC, LQTS, Brugada)
- **cancer**: 84 genes (HBOC, Lynch, Li-Fraumeni, FAP)
- **neuro**: 60 genes (AD, PD, ALS, CMT, SMA)
- **custom_template**: Template for custom panels

## FHIR Phenotype Integration (Differentiator)

This is what separates V2F from competitors (SeqOne, Franklin, VarSeq):
- Queries FHIR Condition table for patient diagnoses
- Maps ICD-10 codes to OMIM diseases and candidate genes
- Scores each variant against phenotype profile
- Phenotype-matching VUS are promoted in reports
- Enables PP4 ACMG criterion (phenotype-specific)

Configure in `phenotype` section of pipeline config.
