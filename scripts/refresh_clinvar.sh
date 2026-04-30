#!/bin/bash
# Monthly ClinVar refresh check
# ClinVar is updated monthly. This script verifies
# that the BQ public dataset has the latest release.

echo "Checking ClinVar BQ table..."

bq show --format=prettyjson \
  bigquery-public-data:human_variant_annotation.clinvar_hg38 \
  2>/dev/null | python3 -c "
import sys, json
data = json.load(sys.stdin)
modified = data.get('lastModifiedTime', 'unknown')
rows = data.get('numRows', '0')
size_mb = int(data.get('numBytes', 0)) / 1024 / 1024
print(f'  ClinVar hg38:')
print(f'    Last modified: {modified}')
print(f'    Rows: {rows}')
print(f'    Size: {size_mb:.1f} MB')
" || echo "  Cannot access ClinVar table. Check permissions."
