#!/usr/bin/env python3
"""Generate the bundled SV demo VCF and precomputed CNV outputs.

Curated to demonstrate the full ACMG/ClinGen 2019 tier range:
  - 22q11.2 deletion (~3 Mb DEL on chr22) — Pathogenic (DiGeorge syndrome)
  - BRCA1 exon-level deletion (~6 kb on chr17) — Likely pathogenic
  - 16p11.2 deletion (~600 kb on chr16) — Pathogenic (autism spectrum)
  - PMP22 duplication (~1.5 Mb DUP on chr17) — Pathogenic (CMT1A)
  - 200 kb VUS deletion in non-disease region (chr2) — VUS
  - Common 1.5 kb intronic deletion (chr11) — Benign
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DEMO_DIR = REPO_ROOT / "demo" / "sv"
PRECOMPUTED = REPO_ROOT / "demo" / "precomputed_sv"

VCF_HEADER = """##fileformat=VCFv4.2
##fileDate=20260501
##source=v2f-sv-demo-generator
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr22,length=50818468>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEMO_SV
"""

# (chrom, pos, end, svtype, gene_label, expected_class, n_genes, gnomad_af, notes)
SV_VARIANTS = [
    ("chr22", 18900000, 21500000, "DEL", "22q11.2",
     "Pathogenic", 50, 0.0, "DiGeorge / 22q11.2 deletion syndrome — recurrent pathogenic locus"),
    ("chr16", 29500000, 30200000, "DEL", "16p11.2",
     "Pathogenic", 28, 0.0, "16p11.2 deletion — autism spectrum"),
    ("chr17", 15159001, 16679000, "DUP", "PMP22",
     "Pathogenic", 8, 0.0, "PMP22 duplication — Charcot-Marie-Tooth 1A (triplosensitive)"),
    ("chr17", 43106000, 43112000, "DEL", "BRCA1",
     "Likely pathogenic", 1, 0.0, "BRCA1 exon-level deletion (HI=3)"),
    ("chr2",  100000000, 100200000, "DEL", "intergenic",
     "VUS", 0, 0.00001, "200 kb deletion in gene-poor region — VUS"),
    ("chr11", 5226000, 5227500, "DEL", "intronic",
     "Benign", 0, 0.07, "Common intronic deletion — Section 5 negative"),
]


def _format_record(chrom, pos, end, svtype, gene_label, expected_class,
                    n_genes, gnomad_af, notes):
    svlen = end - pos
    info = f"SVTYPE={svtype};SVLEN={svlen};END={end};AF={gnomad_af:.5g}"
    fmt = "GT:DP"
    sample = "0/1:40"
    return f"{chrom}\t{pos}\t.\tN\t<{svtype}>\t99\tPASS\t{info}\t{fmt}\t{sample}"


def _chrom_key(chrom):
    s = chrom.replace("chr", "")
    return (0, int(s)) if s.isdigit() else (1, s)


def write_vcf(out_path: Path, variants):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    raw = out_path.with_suffix("")
    sorted_v = sorted(variants, key=lambda v: (_chrom_key(v[0]), v[1]))
    with open(raw, "w") as f:
        f.write(VCF_HEADER)
        for v in sorted_v:
            f.write(_format_record(*v) + "\n")
    if shutil.which("bgzip"):
        subprocess.run(["bgzip", "-f", str(raw)], check=True)
        if shutil.which("tabix"):
            subprocess.run(["tabix", "-p", "vcf", str(out_path)], check=True)
    else:
        with open(raw, "rb") as src, gzip.open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        os.unlink(raw)


def write_demo_config(out_path: Path):
    out_path.write_text("""# Structural-variant demo configuration for V2F Reporter
input:
  adapter: sv
  analysis_mode: cnv
  vcf_path: data/vcf/sv_demo.vcf.gz
  reference_genome: GRCh38

annotation:
  vep_mode: skip   # demo skips VEP

databases:
  gnomad_sv:
    enabled: false
    bq_table: ""   # optional; bundled fallback used when empty
  clingen_dosage:
    enabled: true   # uses bundled reference/clingen_dosage_sensitivity.json
  decipher:
    enabled: false
    bq_table: ""

acmg:
  enable_pvs1: false   # ACMG/AMP not used; ACMG/ClinGen 2019 applies
  enable_population_criteria: false

phenotype:
  enabled: false   # not used for SV interpretation

output:
  data_dir: data
  reports_dir: reports
  report_format: html

pipeline:
  parallel_workers: 4
  log_level: INFO
""")


def write_precomputed_outputs():
    rows = []
    for v in SV_VARIANTS:
        chrom, pos, end, svtype, gene, expected, n_genes, af, notes = v
        rows.append({
            "variant_id": f"{chrom}_{pos}_{svtype}",
            "chrom": chrom, "pos": pos, "ref": "N", "alt": f"<{svtype}>",
            "sample_id": "DEMO_SV",
            "gene": gene,
            "consequence": "structural_variant",
            "severity": "HIGH" if expected.startswith(("Pathogenic", "Likely")) else "MODIFIER",
            "genotype": "0/1", "read_depth": 40, "allele_fraction": 0.5,
            "svtype": svtype,
            "svlen": end - pos,
            "end_pos": end,
            "n_genes": n_genes,
            "gnomad_sv_af": af,
            "gnomad_sv_ac": int(af * 100000) if af else 0,
            "cnv_classification": expected,
            "cnv_score": {
                "Pathogenic": 1.5, "Likely pathogenic": 0.95,
                "VUS": 0.0, "Likely benign": -0.95, "Benign": -1.30,
            }.get(expected, 0.0),
            "cnv_evidence_summary": notes,
            "cnv_section_scores": json.dumps({
                "section_1_genomic_content": 1.0 if n_genes >= 35 else (0.5 if n_genes >= 25 else (0.0 if n_genes else -0.5)),
                "section_2_dosage_genes": 1.0 if expected == "Pathogenic" else (0.5 if "Likely" in expected else 0.0),
                "section_3_literature": 0.0,
                "section_4_inheritance": 0.0,
                "section_5_population_frequency": -1.0 if af and af >= 0.05 else (0.1 if af and af < 1e-5 else 0.0),
            }),
            # Mirror the small-variant columns so the dashboard's existing UI still works
            "acmg_classification": expected,
            "acmg_criteria": notes,
            "acmg_tier": 1 if expected.startswith(("Pathogenic", "Likely path")) else (3 if "benign" in expected.lower() else 2),
            "gnomad_af": af,  # legacy column, used by the existing Variants table
            "clinvar_classification": "",
            "clinvar_review_stars": 0,
            "cadd_phred": None,
            "revel": None,
            "spliceai_max": None,
            "phenotype_match_score": 0.0,
        })
    df = pd.DataFrame(rows)
    base = PRECOMPUTED / "data"
    for sub in ("vcf", "annotated", "enriched", "classified"):
        (base / sub).mkdir(parents=True, exist_ok=True)

    s1_cols = ["variant_id", "chrom", "pos", "ref", "alt", "genotype",
               "read_depth", "allele_fraction", "svtype", "svlen", "end_pos"]
    df[s1_cols].to_parquet(base / "vcf" / "variants.parquet", index=False)
    df.to_parquet(base / "annotated" / "annotated_variants.parquet", index=False)
    df.to_parquet(base / "enriched" / "variants_enriched.parquet", index=False)
    df.to_parquet(base / "classified" / "acmg_results.parquet", index=False)

    (base / "phenotype").mkdir(parents=True, exist_ok=True)
    (base / "phenotype" / "patient_phenotype.json").write_text(
        json.dumps({"skipped": True, "reason": "CNV mode"}, indent=2))

    eval_dir = PRECOMPUTED / "eval"
    eval_dir.mkdir(parents=True, exist_ok=True)
    (eval_dir / "validation_summary.json").write_text(json.dumps({
        "n_variants": len(df), "demo": True,
        "framework": "ACMG/ClinGen 2019 CNV (Riggs 2019)",
    }, indent=2))

    reports_dir = PRECOMPUTED / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    rows_html = []
    for _, r in df.sort_values("cnv_score", ascending=False).iterrows():
        rows_html.append(
            f"<tr><td><strong>{r.gene}</strong></td>"
            f"<td>{r.chrom}:{r.pos:,}-{r.end_pos:,}</td>"
            f"<td>{r.svtype}</td>"
            f"<td>{r.svlen:,} bp</td>"
            f"<td><span class='cls cls-{r.cnv_classification.replace(' ', '-')}'>{r.cnv_classification}</span></td>"
            f"<td>{r.cnv_score:+.2f}</td>"
            f"<td>{r.n_genes}</td>"
            f"<td style='font-size:11px;'>{r.cnv_evidence_summary[:200]}</td></tr>"
        )
    (reports_dir / "DEMO_SV_report.html").write_text(f"""<!DOCTYPE html><html><head><meta charset='UTF-8'>
<title>V2F CNV Report — DEMO_SV</title>
<style>
body{{font-family:'Segoe UI',sans-serif;max-width:1300px;margin:24px auto;padding:0 20px;color:#333}}
h1{{color:#00838f;border-bottom:2px solid #00838f;padding-bottom:8px}}
.banner{{background:#e0f7fa;border-left:4px solid #00838f;padding:12px 18px;margin:20px 0;border-radius:4px;color:#006064}}
table{{width:100%;border-collapse:collapse;margin-top:12px;font-size:13px}}
th,td{{padding:8px 10px;border-bottom:1px solid #eee;text-align:left;vertical-align:top}}
th{{background:#e0f7fa}}
.cls{{display:inline-block;padding:2px 8px;border-radius:3px;color:#fff;font-weight:600;font-size:11px}}
.cls-Pathogenic{{background:#c62828}}.cls-Likely-pathogenic{{background:#ef6c00}}
.cls-VUS{{background:#fbc02d;color:#333}}
.cls-Likely-benign{{background:#66bb6a}}.cls-Benign{{background:#2e7d32}}
</style></head><body>
<h1>V2F CNV Report — DEMO_SV</h1>
<div class="banner"><strong>Demo report.</strong> Generated from synthetic data — not for clinical use.</div>
<p><strong>Framework:</strong> ACMG/ClinGen 2019 5-section CNV scoring (Riggs et al. PMID 31690835).<br>
<strong>Generated:</strong> {datetime.utcnow().strftime('%Y-%m-%d')}</p>
<h2>Classified Structural Variants</h2>
<table>
<thead><tr><th>Region</th><th>Coordinates</th><th>Type</th><th>Length</th><th>Classification</th><th>Score</th><th>Genes</th><th>Evidence</th></tr></thead>
<tbody>{''.join(rows_html)}</tbody>
</table>
</body></html>""")


def write_readme():
    md = ["# V2F SV / CNV Demo Data", "",
          "Bundled synthetic SV VCF for ACMG/ClinGen 2019 CNV scoring demo.",
          "**Not for clinical use.**", "",
          "## Variants (ground truth)", "",
          "| Region | Coords | Type | Size | Class | Notes |",
          "|---|---|---|---|---|---|"]
    for chrom, pos, end, svtype, gene, expected, n_genes, af, notes in SV_VARIANTS:
        md.append(f"| {gene} | {chrom}:{pos:,}-{end:,} | {svtype} | "
                  f"{end-pos:,} bp | {expected} | {notes} |")
    (DEMO_DIR / "README.md").write_text("\n".join(md) + "\n")


def main():
    DEMO_DIR.mkdir(parents=True, exist_ok=True)
    write_vcf(DEMO_DIR / "sv_demo.vcf.gz", SV_VARIANTS)
    write_demo_config(DEMO_DIR / "config.yaml")
    write_readme()
    write_precomputed_outputs()
    print(f"Wrote SV demo under {DEMO_DIR.relative_to(REPO_ROOT)}/")
    print(f"Wrote precomputed outputs under {PRECOMPUTED.relative_to(REPO_ROOT)}/")


if __name__ == "__main__":
    sys.exit(main() or 0)
