#!/bin/bash
# ClinGen Dosage Sensitivity refresh.
#
# Pulls the latest ClinGen DS Map (HI=3 + TS=3 genes + recurrent loci)
# and rewrites reference/clingen_dosage_sensitivity.json. The bundled
# default in V2F is a curated subset of ~15 genes for demo purposes;
# production deployments need the full ~600-gene + recurrent-locus list.
#
# Source: https://search.clinicalgenome.org/kb/gene-dosage
# Public download: https://search.clinicalgenome.org/kb/gene-dosage/download

set -euo pipefail

WORK="${TMPDIR:-/tmp}/v2f_clingen_refresh"
mkdir -p "$WORK"

echo "→ Downloading ClinGen Dosage Sensitivity Map…"
curl -sSf -o "$WORK/clingen_ds.tsv" \
  "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv" \
  || { echo "  ! ClinGen FTP download failed; please configure mirror"; exit 1; }

echo "→ Building V2F-format JSON…"
WORK="$WORK" python3 - <<'PY'
import csv, json, os, datetime

work = os.environ["WORK"]
out_path = "reference/clingen_dosage_sensitivity.json"

hi3_genes = []
ts3_genes = []
with open(f"{work}/clingen_ds.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        gene = row.get("Gene Symbol") or row.get("gene_symbol") or ""
        chrom = row.get("Chromosome") or row.get("chromosome") or ""
        if not chrom.startswith("chr") and chrom:
            chrom = f"chr{chrom}"
        try:
            start = int(row.get("Genomic Location Start") or 0)
            end = int(row.get("Genomic Location End") or 0)
        except (TypeError, ValueError):
            continue
        hi = row.get("Haploinsufficiency Score") or row.get("HI Score")
        ts = row.get("Triplosensitivity Score") or row.get("TS Score")
        if hi == "3":
            hi3_genes.append({"gene": gene, "chrom": chrom,
                                "start": start, "end": end})
        if ts == "3":
            ts3_genes.append({"gene": gene, "chrom": chrom,
                                "start": start, "end": end})

existing = {}
if os.path.exists(out_path):
    with open(out_path) as f:
        existing = json.load(f)

doc = {
    "_source": "ClinGen Dosage Sensitivity Map (refreshed)",
    "_refreshed_at": datetime.datetime.utcnow().isoformat() + "Z",
    "_haploinsufficient_3": hi3_genes,
    "_triplosensitive_3": ts3_genes,
    "_recurrent_pathogenic_loci": existing.get("_recurrent_pathogenic_loci", []),
    "_recurrent_benign_loci": existing.get("_recurrent_benign_loci", []),
}

with open(out_path, "w") as f:
    json.dump(doc, f, indent=2)

print(f"  ClinGen DS refreshed: {len(hi3_genes)} HI=3 genes, "
      f"{len(ts3_genes)} TS=3 genes")
PY

STATUS_DIR="${V2F_REFRESH_STATUS_DIR:-./reference/refresh_status}"
mkdir -p "$STATUS_DIR"
cat > "$STATUS_DIR/clingen_dosage.json" <<EOF
{
  "source": "clingen_dosage_sensitivity",
  "refreshed_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "file": "reference/clingen_dosage_sensitivity.json"
}
EOF

echo
echo "✓ ClinGen DS refreshed → reference/clingen_dosage_sensitivity.json"
