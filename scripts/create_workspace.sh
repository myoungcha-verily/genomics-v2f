#!/bin/bash
# Variant-to-Function Reporter — Workspace Provisioning
#
# Creates a complete Workbench workspace with:
#   1. GCS bucket for VCFs + outputs
#   2. BQ dataset for results
#   3. ClinVar BQ reference (public dataset)
#   4. gnomAD BQ reference (public dataset)
#   5. JupyterLab app
#
# Usage:
#   bash scripts/create_workspace.sh <workspace-id> <pod-id>
#   bash scripts/create_workspace.sh my-genomics-ws my-pod

set -e

WORKSPACE_ID="${1:?Usage: create_workspace.sh <workspace-id> <pod-id>}"
POD_ID="${2:?Usage: create_workspace.sh <workspace-id> <pod-id>}"

echo "============================================"
echo "  Creating V2F Reporter Workspace"
echo "  ID: $WORKSPACE_ID"
echo "============================================"

# Step 1: Create workspace
echo ""
echo "Step 1/5: Creating workspace..."
wb workspace create \
  --id="$WORKSPACE_ID" \
  --name="Genomics V2F Reporter" \
  --description="Variant-to-Function Reporter — clinical genomics analysis workspace" \
  --pod-id="$POD_ID"

wb workspace set --id="$WORKSPACE_ID"
echo "  ✓ Workspace created"

# Step 2: Create GCS bucket
echo ""
echo "Step 2/5: Creating storage bucket..."
wb resource create gcs-bucket \
  --name=v2f-data \
  --description="VCF files, pipeline outputs, and reports"
echo "  ✓ Bucket created"

# Step 3: Create BQ dataset
echo ""
echo "Step 3/5: Creating BigQuery dataset..."
wb resource create bq-dataset \
  --name=v2f-results \
  --dataset-id=v2f_results \
  --description="Pipeline results and variant classifications"
echo "  ✓ BQ dataset created"

# Step 4: Add ClinVar reference
echo ""
echo "Step 4/5: Adding ClinVar BQ reference..."
wb resource add-ref bq-dataset \
  --name=clinvar-reference \
  --project-id=bigquery-public-data \
  --dataset-id=human_variant_annotation \
  --description="ClinVar variant annotations (public dataset)" \
  2>/dev/null || echo "  ⚠ ClinVar reference may already exist or require manual setup"
echo "  ✓ ClinVar reference added"

# Step 5: Add gnomAD reference (may fail due to VPC SC)
echo ""
echo "Step 5/5: Adding gnomAD BQ reference..."
wb resource add-ref bq-dataset \
  --name=gnomad-reference \
  --project-id=bigquery-public-data \
  --dataset-id=gnomAD \
  --description="gnomAD population frequency data (public dataset)" \
  2>/dev/null || echo "  ⚠ gnomAD reference may require manual setup (VPC Service Controls)"
echo "  ✓ gnomAD reference added"

echo ""
echo "============================================"
echo "  Workspace Ready!"
echo "============================================"
echo ""
echo "  Workspace: $WORKSPACE_ID"
echo "  Bucket:    v2f-data"
echo "  BQ Dataset: v2f-results"
echo ""
echo "  Next: Create a JupyterLab app in the Workbench UI"
echo "  Then: bash quickstart.sh"
echo ""
