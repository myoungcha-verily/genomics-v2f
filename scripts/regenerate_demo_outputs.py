#!/usr/bin/env python3
"""Generate pre-computed pipeline outputs for the demo VCF.

So a new user clicking "Try with demo data" sees populated tabs (Variants,
Reports, etc.) immediately, without having to wait for the pipeline to run.

Mirrors the schema produced by stages 1-6, but writes plausible values
directly instead of running the actual annotation/enrichment/classification
code. Re-run after any classification logic change so the bundled outputs
match the live pipeline.
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
PRECOMPUTED = REPO_ROOT / "demo" / "precomputed"

# Mirror of DEMO_VARIANTS in generate_demo_data.py + per-variant
# annotation values that would come out of stages 2/3/4.
ROWS = [
    {
        "variant_id": "chr17_43106478_T_G",
        "chrom": "chr17", "pos": 43106478, "ref": "T", "alt": "G",
        "gene": "BRCA1", "consequence": "missense_variant", "severity": "MODERATE",
        "hgvs_c": "c.181T>G", "hgvs_p": "p.Cys61Gly",
        "genotype": "0/1", "read_depth": 42, "allele_fraction": 0.50,
        "gnomad_af": 0.000005, "gnomad_af_popmax": 0.00012,
        "clinvar_classification": "Pathogenic", "clinvar_review_stars": 4,
        "cadd_phred": 28.4, "revel": 0.842, "spliceai_max": 0.04,
        "alphamissense": 0.91,
        "acmg_classification": "Pathogenic",
        "acmg_criteria": "PVS1, PM2, PP3, PP5",
        "phenotype_match_score": 0.95,
    },
    {
        "variant_id": "chr17_7674220_C_T",
        "chrom": "chr17", "pos": 7674220, "ref": "C", "alt": "T",
        "gene": "TP53", "consequence": "missense_variant", "severity": "MODERATE",
        "hgvs_c": "c.524G>A", "hgvs_p": "p.Arg175His",
        "genotype": "0/1", "read_depth": 38, "allele_fraction": 0.47,
        "gnomad_af": 0.000003, "gnomad_af_popmax": 0.00008,
        "clinvar_classification": "Pathogenic", "clinvar_review_stars": 3,
        "cadd_phred": 32.0, "revel": 0.927, "spliceai_max": 0.02,
        "alphamissense": 0.97,
        "acmg_classification": "Pathogenic",
        "acmg_criteria": "PS1, PM1, PM2, PP3, PP5",
        "phenotype_match_score": 0.78,
    },
    {
        "variant_id": "chr14_23427594_G_A",
        "chrom": "chr14", "pos": 23427594, "ref": "G", "alt": "A",
        "gene": "MYH7", "consequence": "missense_variant", "severity": "MODERATE",
        "hgvs_c": "c.4954G>A", "hgvs_p": "p.Glu1652Lys",
        "genotype": "0/1", "read_depth": 51, "allele_fraction": 0.49,
        "gnomad_af": 0.0,
        "clinvar_classification": "Likely pathogenic", "clinvar_review_stars": 2,
        "cadd_phred": 26.7, "revel": 0.683, "spliceai_max": 0.13,
        "alphamissense": 0.71,
        "acmg_classification": "Likely pathogenic",
        "acmg_criteria": "PM1, PM2, PP3, PP5",
        "phenotype_match_score": 0.05,
    },
    {
        "variant_id": "chr10_87864458_C_T",
        "chrom": "chr10", "pos": 87864458, "ref": "C", "alt": "T",
        "gene": "PTEN", "consequence": "missense_variant", "severity": "MODERATE",
        "hgvs_c": "c.493G>A", "hgvs_p": "p.Glu165Lys",
        "genotype": "0/1", "read_depth": 35, "allele_fraction": 0.45,
        "gnomad_af": 0.0001,
        "clinvar_classification": "Uncertain significance", "clinvar_review_stars": 1,
        "cadd_phred": 22.1, "revel": 0.512, "spliceai_max": 0.08,
        "alphamissense": 0.45,
        "acmg_classification": "VUS",
        "acmg_criteria": "PM2, PP3",
        "phenotype_match_score": 0.62,
    },
    {
        "variant_id": "chr13_32340301_A_G",
        "chrom": "chr13", "pos": 32340301, "ref": "A", "alt": "G",
        "gene": "BRCA2", "consequence": "synonymous_variant", "severity": "LOW",
        "hgvs_c": "c.7242A>G", "hgvs_p": "p.Lys2414=",
        "genotype": "0/1", "read_depth": 47, "allele_fraction": 0.48,
        "gnomad_af": 0.058,
        "clinvar_classification": "Benign", "clinvar_review_stars": 3,
        "cadd_phred": 8.2, "revel": 0.045, "spliceai_max": 0.01,
        "alphamissense": 0.08,
        "acmg_classification": "Likely benign",
        "acmg_criteria": "BS1, BP4, BP6",
        "phenotype_match_score": 0.00,
    },
    {
        "variant_id": "chr17_43124027_G_A",
        "chrom": "chr17", "pos": 43124027, "ref": "G", "alt": "A",
        "gene": "BRCA1", "consequence": "synonymous_variant", "severity": "LOW",
        "hgvs_c": "c.4837A>G (synonymous)", "hgvs_p": "p.Lys1613=",
        "genotype": "0/1", "read_depth": 44, "allele_fraction": 0.48,
        "gnomad_af": 0.072,
        "clinvar_classification": "Benign", "clinvar_review_stars": 4,
        "cadd_phred": 5.1, "revel": 0.012, "spliceai_max": 0.00,
        "alphamissense": 0.04,
        "acmg_classification": "Benign",
        "acmg_criteria": "BA1, BP4, BP6, BP7",
        "phenotype_match_score": 0.00,
    },
    {
        "variant_id": "chr17_7676154_G_C",
        "chrom": "chr17", "pos": 7676154, "ref": "G", "alt": "C",
        "gene": "TP53", "consequence": "intron_variant", "severity": "MODIFIER",
        "hgvs_c": "intronic", "hgvs_p": "",
        "genotype": "0/1", "read_depth": 41, "allele_fraction": 0.49,
        "gnomad_af": 0.13,
        "clinvar_classification": "Benign", "clinvar_review_stars": 2,
        "cadd_phred": 2.3, "revel": 0.0, "spliceai_max": 0.00,
        "alphamissense": 0.0,
        "acmg_classification": "Benign",
        "acmg_criteria": "BA1, BP4",
        "phenotype_match_score": 0.00,
    },
    {
        "variant_id": "chr11_5226750_T_C",
        "chrom": "chr11", "pos": 5226750, "ref": "T", "alt": "C",
        "gene": "HBB", "consequence": "intergenic_variant", "severity": "MODIFIER",
        "hgvs_c": "intergenic", "hgvs_p": "",
        "genotype": "0/1", "read_depth": 39, "allele_fraction": 0.50,
        "gnomad_af": 0.21,
        "clinvar_classification": "", "clinvar_review_stars": 0,
        "cadd_phred": 1.1, "revel": 0.0, "spliceai_max": 0.00,
        "alphamissense": 0.0,
        "acmg_classification": "Benign",
        "acmg_criteria": "BA1",
        "phenotype_match_score": 0.00,
    },
    {
        "variant_id": "chr10_87894077_G_A",
        "chrom": "chr10", "pos": 87894077, "ref": "G", "alt": "A",
        "gene": "PTEN", "consequence": "3_prime_UTR_variant", "severity": "MODIFIER",
        "hgvs_c": "3' UTR", "hgvs_p": "",
        "genotype": "0/1", "read_depth": 40, "allele_fraction": 0.48,
        "gnomad_af": 0.085,
        "clinvar_classification": "Likely benign", "clinvar_review_stars": 1,
        "cadd_phred": 4.4, "revel": 0.0, "spliceai_max": 0.00,
        "alphamissense": 0.0,
        "acmg_classification": "Likely benign",
        "acmg_criteria": "BS1, BP6",
        "phenotype_match_score": 0.05,
    },
    {
        "variant_id": "chr14_23429834_C_T",
        "chrom": "chr14", "pos": 23429834, "ref": "C", "alt": "T",
        "gene": "MYH7", "consequence": "intron_variant", "severity": "MODIFIER",
        "hgvs_c": "intronic", "hgvs_p": "",
        "genotype": "0/1", "read_depth": 46, "allele_fraction": 0.50,
        "gnomad_af": 0.092,
        "clinvar_classification": "", "clinvar_review_stars": 0,
        "cadd_phred": 3.0, "revel": 0.0, "spliceai_max": 0.00,
        "alphamissense": 0.0,
        "acmg_classification": "Benign",
        "acmg_criteria": "BA1",
        "phenotype_match_score": 0.00,
    },
]


def write_parquets():
    df = pd.DataFrame(ROWS)
    # Stage 1: minimal variant table
    s1_cols = ["variant_id", "chrom", "pos", "ref", "alt", "genotype",
               "read_depth", "allele_fraction"]
    s1 = (PRECOMPUTED / "data" / "vcf" / "variants.parquet")
    s1.parent.mkdir(parents=True, exist_ok=True)
    df[s1_cols].to_parquet(s1, index=False)
    # Stage 2: + annotation
    s2_cols = s1_cols + ["gene", "consequence", "severity", "hgvs_c", "hgvs_p",
                         "cadd_phred", "revel", "spliceai_max", "alphamissense"]
    s2 = (PRECOMPUTED / "data" / "annotated" / "annotated_variants.parquet")
    s2.parent.mkdir(parents=True, exist_ok=True)
    df[s2_cols].to_parquet(s2, index=False)
    # Stage 3: + enrichment
    s3_cols = s2_cols + ["gnomad_af", "clinvar_classification",
                         "clinvar_review_stars"]
    if "gnomad_af_popmax" in df.columns:
        s3_cols.append("gnomad_af_popmax")
    s3 = (PRECOMPUTED / "data" / "enriched" / "variants_enriched.parquet")
    s3.parent.mkdir(parents=True, exist_ok=True)
    df[s3_cols].to_parquet(s3, index=False)
    # Stage 4: + ACMG (full dashboard schema)
    s4 = (PRECOMPUTED / "data" / "classified" / "acmg_results.parquet")
    s4.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(s4, index=False)
    return df


def write_phenotype_summary():
    out = PRECOMPUTED / "data" / "phenotype" / "patient_phenotype.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "patient_id": "DEMO_PROBAND",
        "conditions": [
            {"icd10": "Z80.41",
             "display": "Family history of malignant neoplasm of ovary"},
            {"icd10": "Z80.3",
             "display": "Family history of malignant neoplasm of breast"},
        ],
        "candidate_genes": ["BRCA1", "BRCA2", "TP53", "PTEN", "PALB2",
                            "CHEK2", "ATM", "RAD51C"],
        "n_variants_phenotype_matching": 4,
        "phenotype_match_threshold": 0.5,
        "generated_at": datetime.utcnow().isoformat() + "Z",
    }
    out.write_text(json.dumps(payload, indent=2))
    return out


def write_report(df):
    out = PRECOMPUTED / "reports" / "DEMO_PROBAND_report.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = df.sort_values("phenotype_match_score", ascending=False)
    rank_order = {"Pathogenic": 0, "Likely pathogenic": 1, "VUS": 2,
                  "Likely benign": 3, "Benign": 4}
    rows = rows.assign(_rank=rows["acmg_classification"].map(rank_order))
    rows = rows.sort_values(["_rank", "phenotype_match_score"],
                            ascending=[True, False])

    table_rows = []
    for _, r in rows.iterrows():
        table_rows.append(
            f"<tr><td><strong>{r.gene}</strong></td>"
            f"<td>{r.chrom}:{r.pos}</td>"
            f"<td>{r.ref}&gt;{r.alt}</td>"
            f"<td>{r.hgvs_c}</td>"
            f"<td>{r.acmg_classification}</td>"
            f"<td>{r.acmg_criteria}</td>"
            f"<td>{r.gnomad_af:.2e}</td>"
            f"<td>{r.cadd_phred:.1f}</td>"
            f"<td>{r.phenotype_match_score:.2f}</td></tr>"
        )

    html = f"""<!DOCTYPE html><html><head><meta charset="UTF-8">
<title>V2F Report — DEMO_PROBAND</title>
<style>
body{{font-family:'Segoe UI',sans-serif;max-width:1100px;margin:24px auto;padding:0 20px;color:#333}}
h1{{color:#1565c0;border-bottom:2px solid #1565c0;padding-bottom:8px}}
h2{{color:#0d47a1;margin-top:28px}}
table{{width:100%;border-collapse:collapse;margin-top:12px;font-size:13px}}
th,td{{padding:8px 10px;border-bottom:1px solid #eee;text-align:left}}
th{{background:#f5f5f5}}
.banner{{background:#fff8e1;border-left:4px solid #f9a825;padding:12px 18px;margin:20px 0;border-radius:4px;color:#5d4037}}
.summary{{display:flex;gap:14px;margin-top:14px}}
.card{{flex:1;background:#f8f9fa;border-radius:6px;padding:14px;text-align:center}}
.card .v{{font-size:24px;font-weight:700;color:#1565c0}}
.card .l{{font-size:11px;text-transform:uppercase;color:#888;letter-spacing:1px}}
</style></head><body>
<h1>V2F Variant Interpretation Report — DEMO_PROBAND</h1>
<div class="banner"><strong>Demo report.</strong> Generated from synthetic data — not for clinical use.</div>

<h2>Patient Summary</h2>
<p><strong>Patient ID:</strong> DEMO_PROBAND<br>
<strong>Reference:</strong> GRCh38 &middot; <strong>Adapter:</strong> single_sample &middot;
<strong>Phenotype:</strong> HBOC (Z80.3, Z80.41)</p>

<div class="summary">
  <div class="card"><div class="v">{(rows.acmg_classification == 'Pathogenic').sum()}</div><div class="l">Pathogenic</div></div>
  <div class="card"><div class="v">{(rows.acmg_classification == 'Likely pathogenic').sum()}</div><div class="l">Likely path.</div></div>
  <div class="card"><div class="v">{(rows.acmg_classification == 'VUS').sum()}</div><div class="l">VUS</div></div>
  <div class="card"><div class="v">{(rows.acmg_classification == 'Likely benign').sum()}</div><div class="l">Likely benign</div></div>
  <div class="card"><div class="v">{(rows.acmg_classification == 'Benign').sum()}</div><div class="l">Benign</div></div>
</div>

<h2>Classified Variants</h2>
<table>
<thead><tr><th>Gene</th><th>Position</th><th>Change</th><th>HGVS</th><th>Classification</th><th>ACMG criteria</th><th>gnomAD AF</th><th>CADD</th><th>Phenotype score</th></tr></thead>
<tbody>
{''.join(table_rows)}
</tbody>
</table>

<h2>Notes</h2>
<p>This report was generated from V2F demo data on {datetime.utcnow().strftime('%Y-%m-%d')}.
Phenotype score reflects match between the patient's diagnosed conditions
(BRCA1/2-related cancer family history) and each variant's gene.</p>
</body></html>"""
    out.write_text(html)
    return out


def write_validation_summary():
    """Stage 7 indicator file so the dashboard shows 7/7 stages complete."""
    out = PRECOMPUTED / "eval" / "validation_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "n_variants": 10,
        "n_pathogenic_concordant": 2,
        "n_benign_concordant": 5,
        "concordance_rate": 0.95,
        "demo": True,
    }, indent=2))


def main():
    df = write_parquets()
    write_phenotype_summary()
    write_report(df)
    write_validation_summary()
    n = sum(1 for _ in PRECOMPUTED.rglob("*") if _.is_file())
    print(f"Wrote {n} precomputed files under {PRECOMPUTED.relative_to(REPO_ROOT)}/")


if __name__ == "__main__":
    sys.exit(main() or 0)
