#!/bin/bash
# ClinVar refresh — pull the latest weekly release into a BQ table.
#
# Public bigquery-public-data tables are typically frozen old snapshots
# (the V2F default points at ncbi_clinvar_hg38_20180701, from 2018).
# Production deployments mirror the weekly release into a private BQ
# table; this script does that. Run on a Sunday cron or via Workbench
# scheduled jobs.
#
# Required env vars:
#   V2F_GCP_PROJECT       — destination GCP project
#   V2F_CLINVAR_DATASET   — destination dataset (e.g. v2f_results)
#   V2F_CLINVAR_TABLE     — destination table (e.g. clinvar_weekly)
# Optional:
#   V2F_CLINVAR_BUILD     — GRCh38 (default) or GRCh37
#
# Output: prints the loaded table identifier so the V2F config can be
# updated to point at it.

set -euo pipefail

PROJECT="${V2F_GCP_PROJECT:-${GCP_PROJECT:-$(gcloud config get-value project 2>/dev/null)}}"
DATASET="${V2F_CLINVAR_DATASET:-v2f_results}"
TABLE="${V2F_CLINVAR_TABLE:-clinvar_weekly}"
BUILD="${V2F_CLINVAR_BUILD:-GRCh38}"

if [ -z "$PROJECT" ]; then
  echo "ERROR: V2F_GCP_PROJECT not set and gcloud has no default project" >&2
  exit 1
fi

case "$BUILD" in
  GRCh38) URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" ;;
  GRCh37) URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" ;;
  *) echo "Unknown build: $BUILD" >&2; exit 1 ;;
esac

WORK="${TMPDIR:-/tmp}/v2f_clinvar_refresh"
mkdir -p "$WORK"

echo "→ Downloading ClinVar VCF from NCBI…"
curl -sSf -o "$WORK/clinvar.vcf.gz" "$URL"

echo "→ Extracting release date from header…"
RELEASE=$(zcat "$WORK/clinvar.vcf.gz" | head -50 | \
  grep -oE 'fileDate=[0-9]+' | head -1 | cut -d= -f2)
RELEASE="${RELEASE:-$(date +%Y%m%d)}"
DEST_FQ="$PROJECT.$DATASET.${TABLE}_$RELEASE"
echo "  ClinVar release: $RELEASE → $DEST_FQ"

echo "→ Loading into BigQuery…"
# Using the genomics-vcf-loader pattern (vcf-to-bq) that ClinGen and
# bigquery-public-data themselves use. Caller's responsibility to have
# the loader installed; if not, fall back to a flat schema parse.
if command -v gcloud-genomics-vcf-loader &>/dev/null; then
  gcloud-genomics-vcf-loader \
    --input_pattern="$WORK/clinvar.vcf.gz" \
    --output_table="$DEST_FQ" \
    --representative_header_file="$WORK/clinvar.vcf.gz" \
    --runner=DataflowRunner \
    --project="$PROJECT"
else
  echo "  (gcloud-genomics-vcf-loader not found; using bq load with simple schema)"
  bq mk --force --dataset "$PROJECT:$DATASET" 2>/dev/null || true
  zcat "$WORK/clinvar.vcf.gz" | grep -v '^#' | awk -F'\t' '
    BEGIN { OFS="\t" }
    { print $1, $2, $4, $5, $7, $8 }
  ' > "$WORK/clinvar.tsv"
  bq load --replace --source_format=CSV --field_delimiter=tab \
    "$DEST_FQ" "$WORK/clinvar.tsv" \
    "reference_name:STRING,start_position:INT64,reference_bases:STRING,alternate_bases:STRING,filter:STRING,info:STRING"
fi

echo "→ Updating refresh status JSON…"
STATUS_DIR="${V2F_REFRESH_STATUS_DIR:-./reference/refresh_status}"
mkdir -p "$STATUS_DIR"
cat > "$STATUS_DIR/clinvar.json" <<EOF
{
  "source": "ncbi_clinvar",
  "build": "$BUILD",
  "release": "$RELEASE",
  "table": "$DEST_FQ",
  "refreshed_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
echo
echo "✓ ClinVar refreshed: $DEST_FQ"
echo "  Update databases.clinvar.bq_table in pipeline_config.yaml to use it."
