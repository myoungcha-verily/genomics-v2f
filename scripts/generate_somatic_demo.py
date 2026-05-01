#!/usr/bin/env python3
"""Generate the bundled somatic demo VCF (tumor + matched normal) and
the precomputed AMP-tier outputs used by the dashboard's somatic demo.

Variants chosen to demonstrate the full AMP/ASCO/CAP 2017 tier range:
  - BRAF V600E       — Tier I (FDA-approved targeted therapy in melanoma)
  - EGFR L858R       — Tier I (FDA-approved EGFR TKIs in NSCLC)
  - TP53 R175H       — Tier II (preclinical evidence, common driver)
  - KRAS G12D        — Tier II (off-label therapies under investigation)
  - PIK3CA E545K     — Tier II (combination therapy contexts)
  - 2 noise variants for context (Tier III/IV)
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
DEMO_DIR = REPO_ROOT / "demo" / "somatic"
PRECOMPUTED = REPO_ROOT / "demo" / "precomputed_somatic"

VCF_HEADER = """##fileformat=VCFv4.2
##fileDate=20260501
##source=v2f-somatic-demo-generator
##reference=GRCh38
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr17,length=83257441>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEMO_TUMOR\tDEMO_NORMAL
"""

# (chrom, pos, ref, alt, gene, hgvs_p, expected_tier, drug_targets, notes)
SOMATIC_VARIANTS = [
    ("chr7", 140753336, "A", "T",  "BRAF",   "p.V600E",
     "I", ["Vemurafenib", "Dabrafenib", "Encorafenib"],
     "Canonical melanoma driver; FDA-approved targeted therapy"),
    ("chr7", 55191822,  "T", "G",  "EGFR",   "p.L858R",
     "I", ["Erlotinib", "Gefitinib", "Osimertinib"],
     "NSCLC EGFR-TKI sensitizing mutation"),
    ("chr17", 7674220,  "C", "T",  "TP53",   "p.R175H",
     "II", ["APR-246 (preclinical)"],
     "Hotspot TP53 driver across many cancers"),
    ("chr12", 25245350, "C", "T",  "KRAS",   "p.G12D",
     "II", ["MRTX1133 (clinical trials)"],
     "Common pancreatic / colorectal driver; emerging therapies"),
    ("chr3", 179234297, "G", "A",  "PIK3CA", "p.E545K",
     "II", ["Alpelisib"],
     "PI3K/AKT pathway; HER2- breast cancer combination therapy"),
    ("chr17", 7676154,  "G", "C",  "TP53",   "intronic",
     "IV", [], "Background intronic variant — Tier IV"),
    ("chr12", 25257246, "G", "C",  "KRAS",   "synonymous",
     "IV", [], "Synonymous KRAS variant — Tier IV"),
]


def _chrom_key(chrom):
    s = chrom.replace("chr", "")
    return (0, int(s)) if s.isdigit() else (1, s)


def _format_record(chrom, pos, ref, alt, gene, hgvs_p, tier, drugs, notes,
                   t_vaf=0.42, n_vaf=0.0):
    """Tumor-normal genotype line. Tumor is heterozygous (subclonal-like
    VAF), matched normal is homozygous reference."""
    t_dp = 80
    t_alt = int(t_dp * t_vaf)
    t_ref = t_dp - t_alt
    n_dp = 60
    n_alt = int(n_dp * n_vaf)
    n_ref = n_dp - n_alt
    info = f"DP={t_dp + n_dp};SOMATIC"
    fmt = "GT:DP:AD:AF"
    tumor_field = f"0/1:{t_dp}:{t_ref},{t_alt}:{t_vaf:.3f}"
    normal_field = f"0/0:{n_dp}:{n_ref},{n_alt}:{n_vaf:.3f}"
    return (f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t{info}\t"
            f"{fmt}\t{tumor_field}\t{normal_field}")


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
    out_path.write_text("""# Somatic demo configuration for V2F Reporter
input:
  adapter: somatic
  analysis_mode: somatic
  vcf_path: data/vcf/tumor_normal.vcf.gz
  reference_genome: GRCh38
  somatic:
    tumor_sample_id: DEMO_TUMOR
    normal_sample_id: DEMO_NORMAL

annotation:
  vep_mode: skip   # demo skips VEP

databases:
  civic: {enabled: true, api_url: "https://civicdb.org/api/graphql"}
  oncokb: {enabled: false, api_token: "", api_url: "https://www.oncokb.org/api/v1"}
  cosmic: {enabled: false, bq_table: ""}

acmg_amp:
  vaf_min_somatic: 0.05
  require_paired_normal: false
  knowledge_base_priority: ["oncokb", "civic", "cosmic"]

phenotype:
  enabled: false   # ignored in somatic mode

output:
  data_dir: data
  reports_dir: reports
  report_format: html

pipeline:
  parallel_workers: 4
  log_level: INFO
""")


def write_precomputed_outputs():
    """Produce stage 1-7 parquets + a sample report so the dashboard tabs
    populate immediately when the user clicks 'Load somatic demo'."""
    rows = []
    for v in SOMATIC_VARIANTS:
        chrom, pos, ref, alt, gene, hgvs_p, tier, drugs, notes = v
        rows.append({
            "variant_id": f"{chrom}_{pos}_{ref}_{alt}",
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "sample_id": "DEMO_TUMOR",
            "gene": gene,
            "hgvs_p": hgvs_p,
            "consequence": "missense_variant" if "p." in hgvs_p else "intron_variant",
            "severity": "MODERATE" if "p." in hgvs_p else "MODIFIER",
            "genotype": "0/1",
            "read_depth": 80,
            "allele_fraction": 0.42,
            "tumor_vaf": 0.42,
            "tumor_alt_count": 33,
            "normal_vaf": 0.0,
            "normal_alt_count": 0,
            "is_paired": True,
            "amp_tier": tier,
            "amp_evidence": notes,
            "amp_drug_targets": "; ".join(drugs),
            "amp_kbs_consulted": "civic" if drugs else "civic",
            "amp_evidence_level": "CIViC:A" if tier == "I" else (
                "CIViC:B" if tier == "II" else None),
            "amp_evidence_sources": f"https://civicdb.org/links/genes/{gene}",
            # acmg_classification populated for somatic-only mode
            "acmg_classification": f"AMP Tier {tier}",
            "acmg_criteria": notes,
            "acmg_tier": tier,
            # gnomad / clinvar typically null for somatic drivers
            "gnomad_af": 0.0,
            "clinvar_classification": "",
            "clinvar_review_stars": 0,
            "cadd_phred": 28.0 if "p." in hgvs_p else 5.0,
            "revel": 0.85 if "p." in hgvs_p else 0.0,
            "spliceai_max": 0.0,
            "phenotype_match_score": 0.0,
        })
    df = pd.DataFrame(rows)

    # Stage parquets (mirroring germline demo layout)
    s_cols = {
        1: ["variant_id", "chrom", "pos", "ref", "alt", "genotype",
            "read_depth", "allele_fraction", "tumor_vaf", "normal_vaf"],
        2: None,  # all annotation cols
        3: None,
        4: None,
    }
    base = PRECOMPUTED / "data"
    (base / "vcf").mkdir(parents=True, exist_ok=True)
    df[s_cols[1]].to_parquet(base / "vcf" / "variants.parquet", index=False)
    (base / "annotated").mkdir(parents=True, exist_ok=True)
    df.to_parquet(base / "annotated" / "annotated_variants.parquet", index=False)
    (base / "enriched").mkdir(parents=True, exist_ok=True)
    df.to_parquet(base / "enriched" / "variants_enriched.parquet", index=False)
    (base / "classified").mkdir(parents=True, exist_ok=True)
    df.to_parquet(base / "classified" / "acmg_results.parquet", index=False)

    # Phenotype stub for somatic mode
    (base / "phenotype").mkdir(parents=True, exist_ok=True)
    (base / "phenotype" / "patient_phenotype.json").write_text(
        json.dumps({"skipped": True, "reason": "somatic mode"}, indent=2))

    # Validation summary
    eval_dir = PRECOMPUTED / "eval"
    eval_dir.mkdir(parents=True, exist_ok=True)
    (eval_dir / "validation_summary.json").write_text(json.dumps({
        "n_variants": len(df), "demo": True, "framework": "AMP/ASCO/CAP 2017",
    }, indent=2))

    # Sample report
    reports_dir = PRECOMPUTED / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    rows_html = []
    for _, r in df.sort_values("amp_tier").iterrows():
        rows_html.append(
            f"<tr><td><strong>{r.gene}</strong></td>"
            f"<td>{r.chrom}:{r.pos}</td>"
            f"<td>{r.hgvs_p}</td>"
            f"<td><span class='tier tier-{r.amp_tier}'>Tier {r.amp_tier}</span></td>"
            f"<td>{r.tumor_vaf:.2f}</td>"
            f"<td style='font-size:11px;'>{r.amp_drug_targets}</td>"
            f"<td style='font-size:11px;'>{r.amp_evidence}</td></tr>"
        )
    (reports_dir / "DEMO_TUMOR_report.html").write_text(f"""<!DOCTYPE html><html><head><meta charset='UTF-8'>
<title>V2F Somatic Report — DEMO_TUMOR</title>
<style>
body{{font-family:'Segoe UI',sans-serif;max-width:1200px;margin:24px auto;padding:0 20px;color:#333}}
h1{{color:#7b1fa2;border-bottom:2px solid #7b1fa2;padding-bottom:8px}}
.banner{{background:#fff8e1;border-left:4px solid #f9a825;padding:12px 18px;margin:20px 0;border-radius:4px;color:#5d4037}}
table{{width:100%;border-collapse:collapse;margin-top:12px;font-size:13px}}
th,td{{padding:8px 10px;border-bottom:1px solid #eee;text-align:left;vertical-align:top}}
th{{background:#f3e5f5}}
.tier{{display:inline-block;padding:2px 8px;border-radius:3px;color:#fff;font-weight:600;font-size:11px}}
.tier-I{{background:#c62828}}.tier-II{{background:#ef6c00}}.tier-III{{background:#fbc02d;color:#333}}.tier-IV{{background:#2e7d32}}
</style></head><body>
<h1>V2F Somatic Variant Report — DEMO_TUMOR</h1>
<div class="banner"><strong>Demo report.</strong> Generated from synthetic data — not for clinical use.</div>
<p><strong>Framework:</strong> AMP/ASCO/CAP 2017 four-tier somatic interpretation.<br>
<strong>Generated:</strong> {datetime.utcnow().strftime('%Y-%m-%d')}</p>
<h2>Classified Somatic Variants</h2>
<table>
<thead><tr><th>Gene</th><th>Position</th><th>Protein change</th><th>AMP tier</th><th>VAF</th><th>Drug targets</th><th>Evidence</th></tr></thead>
<tbody>{''.join(rows_html)}</tbody>
</table>
</body></html>""")


def write_readme():
    md = ["# V2F Somatic Demo Data", "",
          "Bundled synthetic tumor-normal data for AMP/ASCO/CAP 2017 demo.",
          "**Not for clinical use.**", "",
          "## Variants (ground truth)", "",
          "| Gene | Protein change | AMP tier | Drugs | Notes |",
          "|---|---|---|---|---|"]
    for chrom, pos, ref, alt, gene, hgvs_p, tier, drugs, notes in SOMATIC_VARIANTS:
        md.append(f"| {gene} | {hgvs_p} | {tier} | {', '.join(drugs) or '—'} | {notes} |")
    (DEMO_DIR / "README.md").write_text("\n".join(md) + "\n")


def main():
    DEMO_DIR.mkdir(parents=True, exist_ok=True)
    write_vcf(DEMO_DIR / "tumor_normal.vcf.gz", SOMATIC_VARIANTS)
    write_demo_config(DEMO_DIR / "config.yaml")
    write_readme()
    write_precomputed_outputs()
    print(f"Wrote somatic demo VCF, config, README under {DEMO_DIR.relative_to(REPO_ROOT)}/")
    print(f"Wrote precomputed outputs under {PRECOMPUTED.relative_to(REPO_ROOT)}/")


if __name__ == "__main__":
    sys.exit(main() or 0)
