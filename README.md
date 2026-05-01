# Variant-to-Function Reporter

Clinical genomics interpretation pipeline that takes raw VCF files all the way to actionable variant reports — across germline, somatic, and structural variants. Built for Verily Workbench.

[![Tests](https://github.com/myoungcha-verily/genomics-v2f/actions/workflows/test.yml/badge.svg)](https://github.com/myoungcha-verily/genomics-v2f/actions/workflows/test.yml)
[![Lint](https://github.com/myoungcha-verily/genomics-v2f/actions/workflows/lint.yml/badge.svg)](https://github.com/myoungcha-verily/genomics-v2f/actions/workflows/lint.yml)

## Highlights

- **Three classification frameworks side-by-side**, picked per variant:
  - Germline → **ACMG/AMP 2015** + Pejaver/ClinGen 2022 calibrated tiers + MaveDB PS3/BS3 + Bayesian posterior + 8-gene VCEP rules
  - Somatic → **AMP/ASCO/CAP 2017** four-tier with CIViC + OncoKB + COSMIC
  - SVs / CNVs → **ACMG/ClinGen 2019** 5-section scoring (Riggs et al.)
- **FHIR phenotype matching** — scores each variant against the patient's diagnosed conditions (key differentiator vs SeqOne / Franklin / VarSeq)
- **Manual literature curation queue** — LitVar2 + ClinGen Allele Registry; researchers add evidence per variant; engine respects overrides
- **Per-stage validation gates** — catches the failure modes (all-VUS, 0-ClinVar-match) that hit M0 in the wild
- **Three bundled demos** — germline, somatic, CNV — populate every dashboard tab without bringing your own VCF
- **Production infrastructure** — CI/CD, run-manifest provenance, structured JSON logging, run history audit trail, reference-data refresh suite
- **CPU-only** — no GPU required

## Quick start

```bash
# 1. Setup
bash quickstart.sh

# 2. Configure (or skip and load a demo from the dashboard)
python3 setup_wizard.py

# 3. Validate
python3 run_pipeline.py --dry-run

# 4. Run
python3 run_pipeline.py

# 5. Dashboard (port 8080)
python3 dashboard/app.py
```

Or, simpler: deploy the dashboard as a Workbench custom app and click **Setup → Load germline demo / Load somatic demo / Load CNV demo**.

## Pipeline stages

| Stage | Script | Description | Time |
|-------|--------|-------------|------|
| 1 | `01_vcf_ingest_qc.py` | Parse, normalize, QC + acceptance gates | 1-5 min |
| 2 | `02_annotate_variants.py` | VEP annotation (skippable) | 5-30 min |
| 3 | `03_database_enrichment.py` | ClinVar + gnomAD + LitVar + ClinGen Allele + gnomAD-SV | 2-10 min |
| 4 | `04_acmg_classification.py` | Germline ACMG / Somatic AMP / CNV — routed by `analysis_mode` | 1-5 min |
| 5 | `05_phenotype_integration.py` | FHIR phenotype matching (germline only) | 1-5 min |
| 6 | `06_report_generation.py` | HTML + optional FHIR / PDF / CSV outputs | < 1 min |
| 7 | `07_validation.py` | ClinVar concordance + run-manifest snapshot | < 1 min |

## Adapters

| Adapter | Use case |
|---|---|
| `single_sample` | Standard WES/WGS VCF |
| `trio` | Multi-sample VCF with proband + parents (de novo detection) |
| `panel` | Targeted gene panel (cardiac, cancer, neuro, custom) |
| `gvcf` | Genomic VCF (GATK HaplotypeCaller output) |
| `somatic` / `somatic_pair` | Tumor-only or tumor-normal pair |
| `sv` / `cnv` | Symbolic-allele VCFs (`<DEL>/<DUP>/<INV>/<INS>/<BND>`) |

## Analysis modes

Set `input.analysis_mode` in `config/pipeline_config.yaml` (or pick in the Setup wizard step 1.5):

- `germline` — ACMG/AMP 2015 + ClinGen SVI extensions
- `somatic` — AMP/ASCO/CAP 2017 four-tier
- `both` — run both engines, write both column families
- `cnv` — handled automatically when symbolic alleles are present, regardless of mode

## Knowledge bases

### Germline
- **ClinVar** — public BQ table or weekly refresh via `scripts/refresh_clinvar.sh`
- **gnomAD v3** (GRCh38) — per-chromosome BQ tables
- **MaveDB** — bundled fallback (BRCA1, TP53, PTEN); full mirror via `scripts/refresh_mavedb.sh`
- **ClinGen Allele Registry** — current ClinVar assertions (catches re-classifications)
- **LitVar2** (NCBI) — PubMed publications per variant

### Somatic
- **CIViC** — open API, default-on
- **OncoKB** — paid API, default-off (graceful no-op without token)
- **COSMIC** — optional BQ lookup

### Structural variants
- **ClinGen Dosage Sensitivity** — bundled subset; full mirror via `scripts/refresh_clingen_dosage.sh`
- **gnomAD-SV** — optional BQ lookup
- Recurrent loci: 22q11.2, 1q21.1, 16p11.2, 15q11-q13, 17p11.2 (PMP22)

## Dashboard tabs

| Tab | Purpose |
|---|---|
| Overview | Pipeline stages, reference DBs, reference-data freshness tile |
| Setup | 7-step wizard: adapter, **analysis_mode**, VCF, reference, VEP, DB tables, phenotype |
| Upload | Drag-drop VCF |
| Pipeline | Run / cancel + live log + demo-loaded notice |
| Variants | Filter by ACMG / AMP / CNV class; gene filter; literature column |
| Reports | HTML preview + PDF / CSV / FHIR exports + run-manifest footer |
| **Curate** | Per-variant LitVar lookup + manual curation form (PS3/BS3, PS4, PP1, BS4, CNV Section 3) |
| **Validation** | Per-stage acceptance gates (warning / hard-fail), color-coded |
| **Runs** | Pipeline run history with git SHA + config hash + duration |
| Settings | Full YAML editor + FHIR connectivity test |

## Reproducibility

Every report has a footer:

```
pipeline=abc12345 · config=sha256:xyz · mode=germline · run=01a2b3c4
```

This identifies the exact code (git SHA), the exact config (sha256 of canonical JSON), the analysis mode, and the run UUID. The full manifest (with reference-data versions) is written to `reports/run_manifest.json`.

Two pipeline runs of the same VCF + same config + same code produce byte-identical parquet outputs.

## Testing

```bash
pytest tests/                # 118 tests, all green
pytest tests/test_amp_engine.py -v
ruff check . && ruff format --check .
```

CI runs the full pytest matrix on Python 3.10 + 3.11 plus ruff / mypy / bandit on every PR.

## Configuration reference

```yaml
input:
  adapter: single_sample | trio | panel | gvcf | somatic | sv
  analysis_mode: germline | somatic | both | cnv
  vcf_path: "path/to/file.vcf.gz"
  reference_genome: GRCh38
  somatic:
    tumor_sample_id: ""
    normal_sample_id: ""

annotation:
  vep_mode: docker | local | skip
  use_calibrated_tiers: true     # Pejaver/ClinGen 2022
  thresholds:
    cadd_pathogenic: 25.3
    revel_pathogenic: 0.644
    spliceai_pathogenic: 0.2
    alphamissense_pathogenic: 0.564

databases:
  clinvar:
    bq_table: "bigquery-public-data.human_variant_annotation.ncbi_clinvar_hg38_20180701"
  gnomad:
    bq_table: "bigquery-public-data.gnomAD.v3_genomes__chr{chrom}"
  gnomad_sv:
    bq_table: ""
  civic: {enabled: true, api_url: "https://civicdb.org/api/graphql"}
  oncokb: {enabled: false, api_token: "", api_url: "https://www.oncokb.org/api/v1"}
  cosmic: {enabled: false, bq_table: ""}
  mavedb: {bq_table: ""}
  clingen_dosage: {enabled: true}
  curations: {bq_table: ""}      # local-JSONL fallback if unset
  runs: {bq_table: ""}           # local-JSONL fallback if unset

acmg:
  ba1_threshold: 0.05
  bs1_threshold: 0.01
  pm2_threshold: 0.0001
  enable_functional: false       # gates PS3/BS3
  vcep:
    enabled: true                # 8 genes: BRCA1/2, PTEN, TP53, MYH7, MYBPC3, KCNQ1, SCN5A
  bayesian:
    enabled: true
    prior_pathogenicity: 0.1
    gene_priors: {}              # per-gene overrides

acmg_amp:
  vaf_min_somatic: 0.05
  knowledge_base_priority: ["oncokb", "civic", "cosmic"]

phenotype:
  enabled: false
  fhir_project: ""
  fhir_dataset: ""

output:
  report_format: html            # html | pdf | both
  fhir_export: false
  csv_export: true

pipeline:
  parallel_workers: 4
  log_level: INFO
  gates:
    mode: warning                # warning | hard_fail
    stage_3:
      clinvar_match_rate_min: 0.0
    stage_4:
      all_vus_max_fraction: 0.95
```

## Reference data refresh

Production deployments shouldn't depend on the bundled demo fallbacks or 2018 ClinVar:

```bash
# Weekly
V2F_GCP_PROJECT=your-project bash scripts/refresh_clinvar.sh

# Monthly
V2F_GCP_PROJECT=your-project bash scripts/refresh_mavedb.sh
bash scripts/refresh_clingen_dosage.sh
```

Each script writes `reference/refresh_status/<source>.json`; the dashboard's Overview tab shows ages and warns when any source is > 30 days stale.

## What's deliberately out of scope

This is a **research tool**. Adding any of the following would require a separate, cross-departmental project:

- HIPAA controls + per-user audit log for PHI access
- CAP/CLIA accreditation + ISO 15189 QMS documentation
- FDA LDT framework / EU IVDR regulatory submission
- Two-person sign-out workflow for clinical reporting
- Per-lab variant interpretation memory database
- Multi-tenant SaaS deployment

## Repository layout

```
genomics-v2f-reporter/
├── adapters/                  # VCF format parsers
├── pipeline/                  # 7 stages + utils
│   ├── 01_vcf_ingest_qc.py
│   ├── ...
│   ├── 07_validation.py
│   └── utils/
│       ├── acmg_engine.py        # germline 28-criteria
│       ├── amp_engine.py         # somatic 4-tier (M1)
│       ├── cnv_engine.py         # SV 5-section (M3)
│       ├── bayesian_acmg.py      # Tavtigian 2018 (M2)
│       ├── in_silico_scores.py   # Pejaver tiered calibration (M2)
│       ├── functional_scores.py  # MaveDB PS3/BS3 (M2)
│       ├── vcep_loader.py        # 8-gene overrides (M2)
│       ├── litvar_client.py      # NCBI literature (P4)
│       ├── clingen_allele_client.py  # canonical CA-IDs (P4)
│       ├── curation_store.py     # manual curation queue (P4)
│       ├── quality_gates.py      # per-stage acceptance (P6)
│       ├── run_manifest.py       # provenance (P5)
│       ├── run_history.py        # audit trail (P9)
│       ├── structured_logger.py  # JSON logging (P9)
│       ├── fhir_exporter.py      # FHIR R4 (P7)
│       ├── pdf_renderer.py       # weasyprint (P7)
│       └── csv_exporter.py       # bulk CSV (P7)
├── reference/
│   ├── insilico_calibration.json     # Pejaver thresholds (M2)
│   ├── clingen_dosage_sensitivity.json  # HI=3 / TS=3 + recurrent loci (M3)
│   ├── vcep_rules/                # 8 per-gene override files (M2)
│   ├── gene_panels/               # cardiac, cancer, neuro
│   └── refresh_status/            # last-refresh timestamps (P10)
├── dashboard/                 # Flask app + index.html (8 tabs)
├── demo/
│   ├── germline/              # 10-variant VCF + precomputed
│   ├── somatic/               # 7-variant tumor-normal VCF + precomputed
│   └── sv/                    # 6-variant SV VCF + precomputed
├── workflows/                 # WDL workflows (Phase 8 — not yet shipped)
├── tests/                     # 118 pytest tests
├── scripts/                   # generators + refresh scripts
└── .github/workflows/         # CI: test, lint
```
