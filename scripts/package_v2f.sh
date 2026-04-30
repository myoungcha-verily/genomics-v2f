#!/bin/bash
# Variant-to-Function Reporter — Packaging Script
#
# Builds a portable tarball of the V2F Reporter for deployment
# into a new Workbench workspace.
#
# What it does:
#   1. Collects all source files (excludes data/output dirs, caches, secrets)
#   2. Generates manifest.txt with file count, sizes, and MD5 checksums
#   3. Creates /tmp/genomics-v2f-reporter.tar.gz
#   4. Uploads to GCS for transfer to the new workspace
#
# Usage:
#   cd ~/genomics-v2f-reporter
#   bash scripts/package_v2f.sh
#
# The tarball unpacks into genomics-v2f-reporter/ (one top-level dir).

set -euo pipefail

# ── Configuration ──────────────────────────────────────────────
PROJECT_DIR="/home/jupyter/genomics-v2f-reporter"
TARBALL="/tmp/genomics-v2f-reporter.tar.gz"
MANIFEST="${PROJECT_DIR}/manifest.txt"
GCS_DEST="gs://myoung-bucket-wb-fleeting-lemon-7624/genomics-v2f-reporter/genomics-v2f-reporter.tar.gz"

echo "============================================"
echo "  V2F Reporter — Package Builder"
echo "============================================"

# ── Step 1: Verify we're in the right directory ────────────────
if [ ! -f "${PROJECT_DIR}/run_pipeline.py" ]; then
    echo "ERROR: ${PROJECT_DIR}/run_pipeline.py not found."
    echo "       Run this script from the genomics-v2f-reporter directory."
    exit 1
fi

cd /home/jupyter

# ── Step 2: Build manifest with checksums ──────────────────────
echo ""
echo "Step 1/4: Building manifest..."

# Find all files to include (respecting exclusions)
find genomics-v2f-reporter \
    -path 'genomics-v2f-reporter/data' -prune -o \
    -path 'genomics-v2f-reporter/reports' -prune -o \
    -path 'genomics-v2f-reporter/eval' -prune -o \
    -path 'genomics-v2f-reporter/logs' -prune -o \
    -path 'genomics-v2f-reporter/__pycache__' -prune -o \
    -path 'genomics-v2f-reporter/venv' -prune -o \
    -path 'genomics-v2f-reporter/.venv' -prune -o \
    -path 'genomics-v2f-reporter/vep_cache' -prune -o \
    -path 'genomics-v2f-reporter/dashboard-app' -prune -o \
    -path 'genomics-v2f-reporter/.ipynb_checkpoints' -prune -o \
    -name '*.pyc' -prune -o \
    -name '*.env' -prune -o \
    -name '.DS_Store' -prune -o \
    -type f -print \
  | sort > /tmp/v2f_file_list.txt

FILE_COUNT=$(wc -l < /tmp/v2f_file_list.txt)
TOTAL_SIZE=$(xargs -a /tmp/v2f_file_list.txt du -cb 2>/dev/null | tail -1 | awk '{print $1}')
TOTAL_SIZE_KB=$((TOTAL_SIZE / 1024))

# Generate manifest header + per-file checksums
{
    echo "# V2F Reporter — Package Manifest"
    echo "# Generated: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    echo "# File count: ${FILE_COUNT}"
    echo "# Total size: ${TOTAL_SIZE_KB} KB (${TOTAL_SIZE} bytes)"
    echo "#"
    echo "# MD5  Size(bytes)  Path"
    echo "# ─────────────────────────────────────────────"
    while IFS= read -r f; do
        md5=$(md5sum "$f" | awk '{print $1}')
        size=$(stat --printf='%s' "$f" 2>/dev/null || stat -f '%z' "$f" 2>/dev/null)
        printf "%s  %8s  %s\n" "$md5" "$size" "$f"
    done < /tmp/v2f_file_list.txt
} > "${MANIFEST}"

echo "  ${FILE_COUNT} files, ${TOTAL_SIZE_KB} KB total"
echo "  Manifest written: manifest.txt"

# ── Step 3: Create tarball ─────────────────────────────────────
echo ""
echo "Step 2/4: Creating tarball..."

# Re-read file list (manifest is now included)
find genomics-v2f-reporter \
    -path 'genomics-v2f-reporter/data' -prune -o \
    -path 'genomics-v2f-reporter/reports' -prune -o \
    -path 'genomics-v2f-reporter/eval' -prune -o \
    -path 'genomics-v2f-reporter/logs' -prune -o \
    -path 'genomics-v2f-reporter/__pycache__' -prune -o \
    -path 'genomics-v2f-reporter/venv' -prune -o \
    -path 'genomics-v2f-reporter/.venv' -prune -o \
    -path 'genomics-v2f-reporter/vep_cache' -prune -o \
    -path 'genomics-v2f-reporter/dashboard-app' -prune -o \
    -path 'genomics-v2f-reporter/.ipynb_checkpoints' -prune -o \
    -name '*.pyc' -prune -o \
    -name '*.env' -prune -o \
    -name '.DS_Store' -prune -o \
    -type f -print \
  | sort > /tmp/v2f_file_list.txt

tar czf "${TARBALL}" -T /tmp/v2f_file_list.txt

TARBALL_SIZE=$(du -h "${TARBALL}" | awk '{print $1}')
TARBALL_COUNT=$(tar tzf "${TARBALL}" | wc -l)

echo "  Tarball: ${TARBALL} (${TARBALL_SIZE})"
echo "  Files in tarball: ${TARBALL_COUNT}"

# ── Step 4: Validate tarball ──────────────────────────────────
echo ""
echo "Step 3/4: Validating tarball..."

# Quick unpack test
TMPDIR=$(mktemp -d)
tar xzf "${TARBALL}" -C "${TMPDIR}"
UNPACKED_COUNT=$(find "${TMPDIR}/genomics-v2f-reporter" -type f | wc -l)

if [ "${UNPACKED_COUNT}" -eq "${TARBALL_COUNT}" ]; then
    echo "  Unpack OK: ${UNPACKED_COUNT} files extracted"
else
    echo "  WARNING: Expected ${TARBALL_COUNT} files, got ${UNPACKED_COUNT}"
fi

# Check key files exist
for f in run_pipeline.py setup_wizard.py quickstart.sh INSTALL.md manifest.txt requirements.txt; do
    if [ -f "${TMPDIR}/genomics-v2f-reporter/${f}" ]; then
        echo "  ✓ ${f}"
    else
        echo "  ✗ MISSING: ${f}"
    fi
done

rm -rf "${TMPDIR}"

# ── Step 5: Upload to GCS ─────────────────────────────────────
echo ""
echo "Step 4/4: Uploading to GCS..."

gsutil cp "${TARBALL}" "${GCS_DEST}"

echo "  Uploaded: ${GCS_DEST}"

# ── Summary ────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Package Complete!"
echo "============================================"
echo ""
echo "  Tarball:  ${TARBALL} (${TARBALL_SIZE})"
echo "  Files:    ${TARBALL_COUNT}"
echo "  GCS:      ${GCS_DEST}"
echo ""
echo "  To deploy in a new workspace:"
echo "    1. Create workspace: bash scripts/create_workspace.sh <ws-id> <pod-id>"
echo "    2. Follow INSTALL.md in the new workspace"
echo ""

# Clean up temp files
rm -f /tmp/v2f_file_list.txt
