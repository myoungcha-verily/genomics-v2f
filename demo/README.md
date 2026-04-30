# V2F Demo Data

Bundled synthetic data so a new user can explore the entire pipeline without bringing their own VCF.

**Not for clinical use.** Variants are real ClinVar entries used as templates with fabricated genotypes for demo purposes.

## Files

- `germline/proband.vcf.gz` (+ `.tbi`) — 10 variants on chr10/11/13/14/17
- `germline/phenotype.json` — synthetic FHIR Condition payload (HBOC)
- `precomputed/` — pre-baked pipeline outputs so the dashboard tabs populate before the user runs the pipeline

## Demo variants — ground truth

| Gene | HGVS | Expected ACMG | Notes |
|------|------|---------------|-------|
| BRCA1 | c.181T>G | Pathogenic | BRCA1 Cys61Gly Ashkenazi founder; ClinVar P with 4-star review |
| TP53 | c.524G>A | Pathogenic | TP53 Arg175His hotspot; classic Li-Fraumeni; ClinVar P |
| MYH7 | c.4954G>A | Likely pathogenic | MYH7 missense; cardiomyopathy panel; synthetic-realistic |
| PTEN | c.493G>A | VUS | PTEN missense; demonstrates VUS handling |
| BRCA2 | c.7242A>G | Likely benign | BRCA2 silent change; common in NFE; ClinVar B/LB |
| BRCA1 | c.4837A>G (synonymous) | Benign | BRCA1 synonymous; demonstrates benign call |
| TP53 | intronic | Benign | TP53 intron 3; common variant |
| HBB | intergenic | Benign | HBB region; common SNP, sanity check |
| PTEN | 3' UTR | Benign | PTEN UTR variant; demonstrates UTR handling |
| MYH7 | intronic | Benign | MYH7 deep intronic |

## Loading the demo

From the dashboard:
1. Click **Try with demo data** in the first-run banner, OR
2. Open the **Setup** tab and click **Load demo data**

From the CLI:
```bash
cp demo/germline/proband.vcf.gz data/vcf/
cp demo/germline/phenotype.json data/phenotype/
cp -r demo/precomputed/* .
```

## Resetting

From the dashboard: **Settings** tab → **Reset to first-run state**.
From the CLI: `rm -rf data reports config/pipeline_config.yaml`.

## Regenerating

If classification logic changes, regenerate the bundled outputs:
```bash
python3 scripts/generate_demo_data.py
python3 scripts/regenerate_demo_outputs.py  # see Milestone 0.2
```
