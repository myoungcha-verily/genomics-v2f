# Variant-to-Function Reporter — AI Assistant Context

You are a genomics pipeline assistant for the V2F Reporter. Help researchers go from VCF (or symbolic-allele VCF for SVs) to actionable variant interpretation reports across **germline**, **somatic**, and **CNV** classification frameworks.

## Critical rules

1. **Python**: Use `python3` (3.10), NOT `python`.
2. **pip**: Use `python3 -m pip`, NOT `pip3`.
3. **Config-driven**: All settings flow from `config/pipeline_config.yaml`. Never hardcode paths.
4. **No GPU needed**: This is a CPU-only pipeline.
5. **GCS persistence**: Local files are ephemeral. Save results to a GCS bucket.
6. **Relative `fetch()`**: Dashboard JavaScript must use relative paths (no leading `/`).
7. **Three classification frameworks**: pick the right one per `analysis_mode` — don't apply ACMG/AMP rules to somatic variants or to CNVs; they have their own engines.

## Quick start

```bash
cd genomics-v2f-reporter
bash quickstart.sh
python3 setup_wizard.py
python3 run_pipeline.py --dry-run
python3 run_pipeline.py
python3 dashboard/app.py     # port 8080
```

Or, faster path: deploy as a Workbench custom app and click **Setup → Load germline demo / somatic demo / CNV demo**.

## Architecture (germline, somatic, CNV — same 7-stage spine)

```
VCF → [1] Ingest/QC → [2] VEP Annotation → [3] DB Enrichment
    → [4] Classification* → [5] Phenotype → [6] Reports → [7] Validation
```

`*` Stage 4 routes per `input.analysis_mode`:
- `germline` → `acmg_engine.classify_variant()` (ACMG/AMP 2015 + Pejaver tiers + MaveDB + Bayesian + VCEP)
- `somatic` → `amp_engine.classify_somatic_variant()` (AMP/ASCO/CAP 2017 four-tier)
- `both` → run both engines, write both column families
- Rows with `svtype` set always go to `cnv_engine.classify_cnv()` (ACMG/ClinGen 2019)

## Key files

| File | Purpose |
|---|---|
| `run_pipeline.py` | Orchestrator. Records each run in run history (P9). |
| `config/pipeline_config.yaml` | Single source of truth. Use the setup wizard to write it. |
| `adapters/` | `single_sample`, `trio`, `panel`, `gvcf`, `somatic`, `sv` |
| `pipeline/utils/acmg_engine.py` | Germline 28-criteria with VCEP overrides + curation overrides |
| `pipeline/utils/amp_engine.py` | Somatic 4-tier; consults CIViC + OncoKB + COSMIC |
| `pipeline/utils/cnv_engine.py` | CNV 5-section scoring; Section 3 driven by curator entries + LitVar |
| `pipeline/utils/bayesian_acmg.py` | Tavtigian 2018 posterior probability |
| `pipeline/utils/in_silico_scores.py` | Pejaver/ClinGen 2022 calibrated PP3/BP4 tiers |
| `pipeline/utils/functional_scores.py` | MaveDB PS3/BS3 lookups |
| `pipeline/utils/vcep_loader.py` | 8 per-gene rule files in `reference/vcep_rules/` |
| `pipeline/utils/litvar_client.py` | NCBI LitVar2 API |
| `pipeline/utils/clingen_allele_client.py` | Canonical CA-IDs + current ClinVar assertions |
| `pipeline/utils/curation_store.py` | Manual curation queue (BQ + local-JSONL fallback) |
| `pipeline/utils/quality_gates.py` | Per-stage acceptance gates (warning|hard_fail) |
| `pipeline/utils/run_manifest.py` | Per-run provenance (git SHA + config hash) |
| `pipeline/utils/run_history.py` | Run audit trail |
| `pipeline/utils/structured_logger.py` | JSON formatter for cloud monitoring |
| `pipeline/utils/fhir_exporter.py` | FHIR R4 DiagnosticReport |
| `dashboard/app.py` | Flask backend (port 8080) — 14 endpoints |
| `reference/insilico_calibration.json` | Pejaver thresholds |
| `reference/clingen_dosage_sensitivity.json` | HI=3 / TS=3 + recurrent loci |
| `reference/vcep_rules/<GENE>.json` | Per-gene VCEP overrides |

## Three classification frameworks

### Germline (default) — ACMG/AMP 2015 + extensions

**28 criteria**, with calibrated tiers via Pejaver/ClinGen 2022:
- Pathogenic: PVS1, PS1-3, PM1-6, PP2-5
- Benign: BA1, BS1, BS3, BP1, BP3-4, BP6-7

PP3/BP4 strength tiers from `reference/insilico_calibration.json` (REVEL/CADD/SpliceAI/AlphaMissense).

PS3/BS3 from `pipeline/utils/functional_scores.py` (MaveDB) — gated on `acmg.enable_functional`.

VCEP overrides per `reference/vcep_rules/{BRCA1,BRCA2,PTEN,TP53,MYH7,MYBPC3,KCNQ1,SCN5A}.json`. Notable: **MYH7 PVS1 explicitly disabled** (LoF is not the disease mechanism); MYBPC3 PVS1 enabled (haploinsufficiency).

Bayesian posterior alongside the categorical call (`bayesian_posterior_prob` column). Tavtigian 2018 weights: PVS=8, S=4, M=2, P=1; BA1=hard override.

### Somatic — AMP/ASCO/CAP 2017

Four tiers: I (strong clinical) / II (potential) / III (unknown / VUS) / IV (benign).

KB priority configurable via `acmg_amp.knowledge_base_priority`. Defaults to `[oncokb, civic, cosmic]` — best-tier-wins merge.

OncoKB requires API token; without one the client logs `error="no_token"` and `amp_engine` falls back to CIViC alone.

VAF guardrail: tumor-only with VAF > 0.45 → flagged "consider germline" in evidence summary.

Stage 5 (phenotype) is automatically skipped for `analysis_mode: somatic`.

### CNV / SV — ACMG/ClinGen 2019 (Riggs)

Five sections:
1. Genomic content (n_genes affected)
2. HI/TS gene overlap + recurrent loci
3. Literature evidence — driven by curator entries (`criterion: CNV_section_3`) + LitVar fallback
4. Inheritance (de novo +0.5)
5. Population frequency (gnomAD-SV-aware)

Total score → 5-tier mapping mirroring ACMG/AMP. Recurrent pathogenic loci (22q11.2, 16p11.2, etc.) hit Section 2 hard.

## Adapter system

```python
from adapters import get_adapter
adapter = get_adapter("single_sample", config)
# Other options: "trio", "panel", "gvcf", "somatic", "somatic_pair", "sv", "cnv"
variants_df = adapter.load_variants()
```

Custom adapters subclass `BaseVCFAdapter` and implement `load_variants()`, `get_sample_ids()`, `get_proband_id()`. Register in `adapters/__init__.py`.

## Config sketch

```yaml
input:
  adapter: single_sample           # or trio, panel, gvcf, somatic, sv
  analysis_mode: germline          # or somatic, both, cnv
  vcf_path: "data/vcf/sample.vcf.gz"
  reference_genome: GRCh38

annotation:
  vep_mode: docker                 # or local, skip
  use_calibrated_tiers: true       # Pejaver/ClinGen 2022

databases:
  clinvar: {bq_table: "..."}
  gnomad:  {bq_table: "...{chrom}"}   # use {chrom} placeholder for v3 per-chrom
  civic:   {enabled: true}
  oncokb:  {enabled: false, api_token: ""}
  mavedb:  {bq_table: ""}             # bundled fallback if unset
  curations: {bq_table: ""}           # local-JSONL fallback

acmg:
  enable_functional: false           # gates PS3/BS3
  vcep: {enabled: true}              # 8 genes
  bayesian: {enabled: true, prior_pathogenicity: 0.1}

acmg_amp:                            # somatic mode
  vaf_min_somatic: 0.05

phenotype:
  enabled: false
  fhir_project: ""

pipeline:
  gates:
    mode: warning                    # or hard_fail
```

## Pipeline commands

```bash
# Full pipeline (all 7 stages)
python3 run_pipeline.py

# Validate config + paths only
python3 run_pipeline.py --dry-run

# Resume from stage 4
python3 run_pipeline.py --start-stage 4

# Just specific stages
python3 run_pipeline.py --stages 1,4,6

# Tests
pytest tests/                                    # 118 tests
pytest tests/test_amp_engine.py -v               # one suite

# Reference data refresh (production deployments)
V2F_GCP_PROJECT=your-project bash scripts/refresh_clinvar.sh
bash scripts/refresh_mavedb.sh
bash scripts/refresh_clingen_dosage.sh

# Dashboard
python3 dashboard/app.py                         # port 8080
```

## Dashboard tabs

| Tab | What's there |
|---|---|
| Overview | Pipeline stages, ref-DB freshness tile (warns >30 d stale) |
| Setup | 7-step wizard. Step 1.5 = `analysis_mode` |
| Upload | Drag-drop VCF → `data/vcf/` |
| Pipeline | Run / Cancel + log + demo-loaded notice |
| Variants | Filter by ACMG / AMP / CNV class |
| Reports | HTML + optional FHIR / PDF / CSV exports |
| **Curate** | LitVar lookup + per-variant evidence form for PS3/BS3, PS4, PP1, BS4, CNV_section_3 |
| **Validation** | Per-stage gate results (PASS / WARN / HARD FAIL header badge) |
| **Runs** | Last 50 pipeline runs with git SHA + config hash + duration |
| Settings | Full YAML editor + FHIR connectivity test button |

## Demos

Three bundled demos load curated outputs into every dashboard tab:

| Demo | Variants | Framework |
|---|---|---|
| **Germline** | 10 (BRCA1, TP53, MYH7, PTEN, BRCA2 + 5 controls) | ACMG/AMP |
| **Somatic** | 7 (BRAF V600E, EGFR L858R, TP53 R175H, KRAS G12D, PIK3CA E545K + controls) | AMP/ASCO/CAP |
| **CNV** | 6 (22q11.2 del, 16p11.2 del, PMP22 dup, BRCA1 exon del + controls) | ACMG/ClinGen 2019 |

Load via `POST /api/demo/load {"mode": "germline"|"somatic"|"cnv"}` or click the buttons in **Setup → Load <X> demo**.

## Reproducibility

Every report has a footer like:
```
pipeline=abc12345 · config=sha256:xyz12345 · mode=germline · run=01a2b3c4
```

The full manifest is `reports/run_manifest.json` and is also embedded in any FHIR DiagnosticReport extension. Same code + same config + same data → identical hashed parquet outputs.

## Quality gates (catches the M0-style failure modes)

`pipeline.gates.mode: warning|hard_fail` controls behavior.

- **stage_1.ti_tv_ratio**: in [1.5, 3.5] for WES — out-of-range often signals contamination
- **stage_3.clinvar_match_rate**: catches schema-drift failures (queries silently returning 0)
- **stage_4.vus_fraction**: > 95% all-VUS is suspicious — usually means upstream enrichment failed

Results aggregated in the **Validation** tab + per-stage summary JSON.

## Curation override flow

Researchers add literature-backed evidence per variant via the Curate tab:
1. Pick a variant from the list
2. LitVar2 publications appear automatically
3. Fill the form: criterion (PS3/BS3/PS4/PP1/BS4/CNV_section_3), action (support_pathogenic/benign/no_call), strength, evidence text, PMIDs
4. Submit → engine reads the latest entry per (variant, criterion) pair on next classification run
5. Reports tag overridden calls with `curated: true` so provenance is preserved

Storage: `databases.curations.bq_table` if set, else local `data/curations.jsonl`.

## Troubleshooting

| Error | Fix |
|---|---|
| `cyvcf2 not found` | `python3 -m pip install cyvcf2` or pipeline uses pysam fallback automatically |
| `Docker not found` | Set `annotation.vep_mode: skip` in config |
| `BQ 403` | Run `gcloud auth application-default login` |
| `No variants` after Run Pipeline | Check VCF path, then run `Validation` tab to inspect stage gates |
| Dashboard 404 on `fetch()` | Ensure paths use no leading `/` (Workbench proxy compatibility) |
| `format_frequency import error` | Already fixed; was M0 polish |
| All variants land in VUS | Stage 3 enrichment gate likely failed — check the configured BQ tables exist |
| OncoKB always returns "no_token" | Expected; set `databases.oncokb.api_token` to enable |
| Phenotype "0 conditions, 0 candidate genes" | Click "Test connection" in Setup step 6 to diagnose FHIR access |

## Gene panels

`reference/gene_panels/`:
- **cardiac**: 82 genes (HCM, DCM, ARVC, LQTS, Brugada)
- **cancer**: 84 genes (HBOC, Lynch, Li-Fraumeni, FAP)
- **neuro**: 60 genes (AD, PD, ALS, CMT, SMA)
- **custom_template**: starting point for new panels

## VCEP rules (ClinGen Variant Curation Expert Panels)

`reference/vcep_rules/`:
- **Cancer**: BRCA1, BRCA2 (ENIGMA), PTEN (PTEN-VCEP), TP53 (TP53-VCEP)
- **Cardiac**: MYH7, MYBPC3 (Cardiomyopathy VCEP), KCNQ1 (LongQT), SCN5A (Channelopathy)

Each file overrides specific criteria (PVS1 applicable/not, PM2 strength, BS1 threshold, etc.). See `pipeline/utils/vcep_loader.py` for merge logic.

## Production-readiness ceiling

This is a **research tool**. The following are deliberately out of scope:
- HIPAA / per-user audit log for PHI
- CAP/CLIA / ISO 15189 QMS
- FDA LDT / EU IVDR regulatory work
- Two-person sign-out workflow
- Per-lab interpretation memory database

These would need a separate cross-departmental project.

## What's running

| Phase | Status |
|---|---|
| M0 — Demo + onboarding | ✅ shipped |
| M0 polish — BQ schemas + banner cohesion | ✅ shipped |
| Pre-M1 polish — FHIR connectivity + report errors | ✅ shipped |
| M1 — Germline / somatic dual mode | ✅ shipped |
| M2 — Pejaver tiers + MaveDB + Bayesian + VCEP | ✅ shipped |
| M3 — Structural variants | ✅ shipped |
| Phase 4 — Literature & curation queue | ✅ shipped |
| Phase 5 — CI/CD + reproducibility | ✅ shipped |
| Phase 6 — Validation gates | ✅ shipped |
| Phase 7 — FHIR / PDF / CSV exporters | ✅ shipped |
| Phase 9 — Observability + run history | ✅ shipped |
| Phase 10 — Reference-data refresh suite | ✅ shipped |
| Phase 8 — WDL workflows + WGS scale | ⏸ deferred (real data needed to design well) |
