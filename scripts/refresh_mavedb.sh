#!/bin/bash
# MaveDB refresh — pull the full MaveDB score-set catalog into BigQuery.
#
# MaveDB hosts deep mutational scanning (DMS) data used for PS3/BS3
# functional evidence. The bundled fallback in V2F covers only 3 genes
# (BRCA1, TP53, PTEN) for demo purposes. Production deployments mirror
# the full catalog so any gene with a published DMS gets coverage.
#
# Reference: https://www.mavedb.org/api/
#
# Required env vars:
#   V2F_GCP_PROJECT
#   V2F_MAVEDB_DATASET    (default v2f_results)
#   V2F_MAVEDB_TABLE      (default mavedb_full)

set -euo pipefail

PROJECT="${V2F_GCP_PROJECT:-${GCP_PROJECT:-$(gcloud config get-value project 2>/dev/null)}}"
DATASET="${V2F_MAVEDB_DATASET:-v2f_results}"
TABLE="${V2F_MAVEDB_TABLE:-mavedb_full}"
DEST_FQ="$PROJECT.$DATASET.$TABLE"

if [ -z "$PROJECT" ]; then
  echo "ERROR: V2F_GCP_PROJECT not set" >&2
  exit 1
fi

WORK="${TMPDIR:-/tmp}/v2f_mavedb_refresh"
mkdir -p "$WORK"

echo "→ Fetching MaveDB score-set list…"
curl -sSf "https://www.mavedb.org/api/v1/score-sets" -o "$WORK/score_sets.json"

echo "→ Materializing per-variant rows…"
python3 - <<'PY'
import json, os, csv
work = os.environ["WORK"]
with open(f"{work}/score_sets.json") as f:
    score_sets = json.load(f)

print(f"  Found {len(score_sets)} MaveDB score sets")
out = open(f"{work}/mavedb_variants.tsv", "w", newline="")
w = csv.writer(out, delimiter="\t")
w.writerow(["score_set_urn", "gene", "protein_change", "score",
            "function_class", "publication_doi"])

import urllib.request
n_rows = 0
for ss in score_sets:
    urn = ss.get("urn", "")
    target = (ss.get("targetGene") or {}).get("name", "")
    pubs = ss.get("publicationIdentifiers", [])
    doi = pubs[0].get("identifier") if pubs else ""
    try:
        with urllib.request.urlopen(
            f"https://www.mavedb.org/api/v1/score-sets/{urn}/scores/", timeout=30
        ) as r:
            scores = json.loads(r.read())
    except Exception as e:
        print(f"  ! skipping {urn}: {e}")
        continue
    for v in scores or []:
        score = v.get("score")
        if score is None: continue
        # MaveDB convention: damaging if score <= -1, tolerated if >= 0
        if score <= -1.0: cls = "damaging"
        elif score >= 0.0: cls = "tolerated"
        else: cls = "intermediate"
        w.writerow([urn, target, v.get("hgvs_pro", ""), f"{score:.4f}", cls, doi])
        n_rows += 1
out.close()
print(f"  Wrote {n_rows} variant rows")
PY

WORK="$WORK" python3 -c "import os; print(os.environ['WORK'])"
echo "→ Loading into BigQuery: $DEST_FQ"
bq mk --force --dataset "$PROJECT:$DATASET" 2>/dev/null || true
bq load --replace --source_format=CSV --field_delimiter=tab \
  "$DEST_FQ" "$WORK/mavedb_variants.tsv" \
  "score_set:STRING,gene:STRING,protein_change:STRING,score:FLOAT64,function_class:STRING,publication_doi:STRING"

STATUS_DIR="${V2F_REFRESH_STATUS_DIR:-./reference/refresh_status}"
mkdir -p "$STATUS_DIR"
cat > "$STATUS_DIR/mavedb.json" <<EOF
{
  "source": "mavedb",
  "table": "$DEST_FQ",
  "refreshed_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF

echo
echo "✓ MaveDB refreshed: $DEST_FQ"
echo "  Update databases.mavedb.bq_table in pipeline_config.yaml to use it."
