# INSTALL — Genomics V2F Reporter

**Snapshot date:** 2026-04-30
**Source tarball:** `gs://myoung-bucket-wb-fleeting-lemon-7624/genomics-v2f-reporter/genomics-v2f-reporter.tar.gz`
**Audience:** Researchers deploying the V2F Reporter into a fresh Workbench workspace.

This document is a **deterministic, copy-paste deployment runbook**.
Every step is a shell command. Run them in order from inside a
JupyterLab terminal in the target workspace. If any step fails,
stop and resolve before moving on — later steps assume earlier ones
succeeded.

> **Read first:**
> - `README.md` for pipeline architecture and usage.
> - `CLAUDE.md` for AI assistant integration.

---

## Section 1 — Prerequisites

### 1.1 Workspace check

```bash
wb workspace describe --format=json | jq '{name, id, googleProjectId}'
```

Confirm you are in the **target** workspace (not the source). Note
the `googleProjectId` — you will need it in Section 4.

If the workspace doesn't exist yet, create it first:

```bash
# Replace <workspace-id> and <pod-id> with your values
wb workspace create \
  --id="<workspace-id>" \
  --name="Genomics V2F Reporter" \
  --description="Variant-to-Function Reporter — clinical genomics analysis workspace" \
  --pod-id="<pod-id>"

wb workspace set --id="<workspace-id>"
```

### 1.2 Auth check

```bash
wb auth status
```

You must be authenticated with WRITER or OWNER role on the
workspace. If not, run `wb auth login`.

### 1.3 Tooling check

```bash
python3 --version            # 3.9+ required
bq version                   # BigQuery CLI
gsutil version               # GCS CLI
python3 -c 'import yaml, pandas; print("OK")'  # Core deps
```

All four must succeed. If `pyyaml` or `pandas` are missing, they
will be installed in Section 5.

---

## Section 2 — Download & Unpack

### 2.1 Download the tarball

```bash
cd /home/jupyter
gsutil cp gs://myoung-bucket-wb-fleeting-lemon-7624/genomics-v2f-reporter/genomics-v2f-reporter.tar.gz /tmp/
```

### 2.2 Unpack

```bash
cd /home/jupyter
tar xzf /tmp/genomics-v2f-reporter.tar.gz
```

This creates `/home/jupyter/genomics-v2f-reporter/` with all source
files.

### 2.3 Confirm extraction

```bash
ls genomics-v2f-reporter/run_pipeline.py \
   genomics-v2f-reporter/setup_wizard.py \
   genomics-v2f-reporter/quickstart.sh \
   genomics-v2f-reporter/manifest.txt
```

All four files must exist. If any is missing, the tarball is corrupt
— re-download and retry.

---

## Section 3 — Validate Contents

Compare extracted file count against the manifest.

```bash
cd /home/jupyter/genomics-v2f-reporter

# File count from manifest header
EXPECTED=$(grep '^# File count:' manifest.txt | awk '{print $NF}')

# Actual file count on disk
ACTUAL=$(find . -type f | wc -l)

echo "Expected: ${EXPECTED}"
echo "Actual:   ${ACTUAL}"

if [ "${ACTUAL}" -ge "${EXPECTED}" ]; then
    echo "✓ File count OK"
else
    echo "✗ File count mismatch — stop and investigate"
fi
```

For per-file integrity, spot-check a few checksums:

```bash
# Pick 3 random files from manifest and verify MD5
grep -v '^#' manifest.txt | shuf | head -3 | while read md5 size path; do
    actual_md5=$(md5sum "${path#genomics-v2f-reporter/}" 2>/dev/null | awk '{print $1}')
    if [ "${actual_md5}" = "${md5}" ]; then
        echo "✓ ${path}"
    else
        echo "✗ ${path} (expected ${md5}, got ${actual_md5})"
    fi
done
```

---

## Section 4 — Provision Resources

Get the workspace's GCP project ID:

```bash
PROJECT=$(wb workspace describe --format=json | jq -r '.googleProjectId')
echo "Project: ${PROJECT}"
```

### 4.1 Create GCS bucket

```bash
wb resource create gcs-bucket \
  --name=v2f-data \
  --description="VCF files, pipeline outputs, and reports"
```

### 4.2 Create BigQuery dataset

```bash
wb resource create bq-dataset \
  --name=v2f-results \
  --dataset-id=v2f_results \
  --description="Pipeline results and variant classifications"
```

### 4.3 Add ClinVar reference

```bash
wb resource add-ref bq-dataset \
  --name=clinvar-reference \
  --project-id=bigquery-public-data \
  --dataset-id=human_variant_annotation \
  --description="ClinVar variant annotations (public dataset)" \
  2>/dev/null || echo "⚠ ClinVar reference may require manual setup"
```

### 4.4 Add gnomAD reference

```bash
wb resource add-ref bq-dataset \
  --name=gnomad-reference \
  --project-id=bigquery-public-data \
  --dataset-id=gnomAD \
  --description="gnomAD population frequency data (public dataset)" \
  2>/dev/null || echo "⚠ gnomAD reference may require manual setup (VPC Service Controls)"
```

> **Note:** ClinVar and gnomAD references may fail if VPC Service
> Controls block access to `bigquery-public-data`. The pipeline
> still works — it will query the public tables directly via the
> BigQuery API.

### 4.5 Confirm resources

```bash
wb resource list --format=json | jq '.[] | {id: .id, type: .resourceType}'
```

You should see at least:
- `v2f-data` (GCS_BUCKET)
- `v2f-results` (BQ_DATASET)
- `clinvar-reference` (BQ_DATASET) — if it succeeded
- `gnomad-reference` (BQ_DATASET) — if it succeeded

---

## Section 5 — Install Dependencies

```bash
cd /home/jupyter/genomics-v2f-reporter
bash quickstart.sh
```

This runs `python3 -m pip install -r requirements.txt` and creates
the output directories (`data/`, `reports/`, `eval/`, `logs/`).

Confirm the install:

```bash
python3 -c '
import cyvcf2, pandas, yaml, jinja2, flask
from google.cloud import bigquery, storage
print("All dependencies OK")
'
```

If `cyvcf2` fails (requires C compiler), the pipeline falls back to
`pysam` automatically — this is fine for most use cases.

---

## Section 6 — Configure Pipeline

### Option A: Interactive wizard (recommended)

```bash
cd /home/jupyter/genomics-v2f-reporter
python3 setup_wizard.py
```

The wizard walks through 6 configuration steps:
1. Input adapter type (single_sample, trio, panel, gvcf)
2. VCF file path
3. Reference genome (GRCh37/GRCh38)
4. VEP annotation mode (docker/local/skip)
5. Database connections (ClinVar, gnomAD BQ tables)
6. Phenotype integration (FHIR project + dataset)

### Option B: Copy an example config

```bash
cp config/example_configs/single_sample_wes.yaml config/pipeline_config.yaml
```

Then edit `config/pipeline_config.yaml` to update:
- `input.vcf_path` — path to your VCF file
- `databases.clinvar.bq_table` — ClinVar BQ table path
- `databases.gnomad.bq_table` — gnomAD BQ table path

Available example configs:
- `single_sample_wes.yaml` — Single-sample whole-exome
- `trio_wes.yaml` — Trio whole-exome (proband + parents)
- `panel_targeted.yaml` — Targeted gene panel
- `research_wgs.yaml` — Research whole-genome

---

## Section 7 — Validate Setup

```bash
cd /home/jupyter/genomics-v2f-reporter
python3 run_pipeline.py --dry-run
```

The dry run validates:
- Config file is parseable
- Input VCF path exists (if configured)
- BigQuery connections are reachable
- Output directories are writable
- All pipeline stages load without import errors

Expected output: `Dry run complete — all checks passed.`

If the dry run reports errors:
- **BQ 403**: Run `gcloud auth application-default login`
- **VCF not found**: Check `input.vcf_path` in config
- **Missing module**: Re-run `bash quickstart.sh`

---

## Section 8 — Start Dashboard

### 8.1 Get the app UUID

```bash
APP_UUID=$(wb app list --format=json | jq -r '.[] | select(.status == "RUNNING") | .id' | head -1)
echo "App UUID: ${APP_UUID}"
```

### 8.2 Launch the dashboard

```bash
cd /home/jupyter/genomics-v2f-reporter
nohup python3 dashboard/app.py > /tmp/v2f_dashboard.log 2>&1 &
sleep 2
```

### 8.3 Verify it's running

```bash
curl -sf http://localhost:8080/ | head -20
```

You should see HTML content. If empty or error, check logs:

```bash
cat /tmp/v2f_dashboard.log
```

### 8.4 Access through Workbench proxy

Open in your browser:

```
https://workbench.verily.com/app/${APP_UUID}/proxy/8080/
```

Replace `${APP_UUID}` with the value from step 8.1.

### 8.5 Confirm port is listening

```bash
ss -ltnp 2>/dev/null | grep ':8080' || netstat -tlnp 2>/dev/null | grep ':8080'
```

Port 8080 must be in `LISTEN` state.

---

## Section 9 — Post-Deploy Checklist

Run through every item below to confirm the deployment is complete.

### 9.1 Resource check

```bash
wb resource list
```

- [ ] `v2f-data` bucket exists
- [ ] `v2f-results` BQ dataset exists
- [ ] ClinVar/gnomAD references exist (or noted as manual setup)

### 9.2 Config check

```bash
cat config/pipeline_config.yaml | head -20
```

- [ ] Config file exists and has correct adapter type
- [ ] VCF path points to a valid file or placeholder

### 9.3 Pipeline check

```bash
python3 run_pipeline.py --dry-run
```

- [ ] Dry run passes without errors

### 9.4 Dashboard check

```bash
curl -sf http://localhost:8080/ > /dev/null && echo "Dashboard OK" || echo "Dashboard NOT running"
```

- [ ] Dashboard responds on port 8080

### 9.5 GCS persistence reminder

> **Local storage is ephemeral.** After running the pipeline, save
> results to the workspace bucket:
>
> ```bash
> BUCKET=$(wb resource describe v2f-data --format=json | jq -r '.metadata.bucketName')
> gsutil -m cp -r reports/ "gs://${BUCKET}/reports/"
> gsutil -m cp -r data/classified/ "gs://${BUCKET}/classified/"
> ```

### 9.6 Environment variables

```bash
env | grep WORKBENCH_v2f
```

- [ ] `WORKBENCH_v2f_data` resolves to the bucket path
- [ ] `WORKBENCH_v2f_results` resolves to the BQ dataset path

---

## Done

If every section above passed without errors, the V2F Reporter is
deployed and ready to process VCF files.

Next steps:
1. Upload a VCF: `gsutil cp sample.vcf.gz gs://<v2f-data-bucket>/vcf/`
2. Configure: `python3 setup_wizard.py`
3. Run: `python3 run_pipeline.py`
4. View: Open dashboard at `https://workbench.verily.com/app/<UUID>/proxy/8080/`

If you encounter a step that fails repeatedly, capture the error
and check `CLAUDE.md` for troubleshooting guidance.

---

*Genomics V2F Reporter — snapshot 2026-04-30*
