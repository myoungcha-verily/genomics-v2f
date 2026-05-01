"""Stage 6: Report Generation

Input: data/classified/acmg_results.parquet, data/phenotype/patient_phenotype.json
Output: reports/<sample>_report.html

Generates per-proband clinical variant reports in HTML format
with tiered variant presentation (P/LP, VUS, LB/B).
"""

import json
import logging
import os
import sys
import time

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipeline.utils.report_renderer import render_proband_report

logger = logging.getLogger(__name__)


def _render_somatic_inline(proband_df, proband_id, config, output_path):
    """Minimal inline somatic report. Used when analysis_mode is 'somatic'.

    Avoids the germline Jinja template (which assumes phenotype + ACMG
    columns) and emphasizes the somatic-specific data: AMP tiers, drug
    targets, VAF.
    """
    from datetime import datetime
    rows_html = []
    df_sorted = proband_df.copy()
    tier_rank = {"I": 0, "II": 1, "III": 2, "IV": 3}
    if "amp_tier" in df_sorted.columns:
        df_sorted["_rank"] = df_sorted["amp_tier"].map(tier_rank).fillna(4)
        df_sorted = df_sorted.sort_values("_rank")

    for _, r in df_sorted.iterrows():
        amp_tier = r.get("amp_tier", "?")
        drugs = r.get("amp_drug_targets", "") or ""
        vaf = float(r.get("tumor_vaf") or r.get("allele_fraction") or 0.0)
        gene = r.get("gene", "")
        hgvs_p = r.get("hgvs_p", "")
        rows_html.append(
            f"<tr><td><strong>{gene}</strong></td>"
            f"<td>{r.get('chrom','')}:{r.get('pos','')}</td>"
            f"<td>{hgvs_p}</td>"
            f"<td><span class='tier tier-{amp_tier}'>Tier {amp_tier}</span></td>"
            f"<td>{vaf:.2f}</td>"
            f"<td style='font-size:11px;'>{drugs[:200]}</td>"
            f"<td style='font-size:11px;'>{(r.get('amp_evidence','') or '')[:160]}</td></tr>"
        )

    n_total = len(df_sorted)
    n_t1 = int((df_sorted.get("amp_tier") == "I").sum()) if "amp_tier" in df_sorted.columns else 0
    n_t2 = int((df_sorted.get("amp_tier") == "II").sum()) if "amp_tier" in df_sorted.columns else 0

    # Manifest footer is appended to the report body just before </body>
    html = f"""<!DOCTYPE html><html><head><meta charset='UTF-8'>
<title>V2F Somatic Report — {proband_id}</title>
<style>
body{{font-family:'Segoe UI',sans-serif;max-width:1200px;margin:24px auto;padding:0 20px;color:#333}}
h1{{color:#7b1fa2;border-bottom:2px solid #7b1fa2;padding-bottom:8px}}
h2{{color:#4a148c;margin-top:28px}}
table{{width:100%;border-collapse:collapse;margin-top:12px;font-size:13px}}
th,td{{padding:8px 10px;border-bottom:1px solid #eee;text-align:left;vertical-align:top}}
th{{background:#f3e5f5}}
.tier{{display:inline-block;padding:2px 8px;border-radius:3px;color:#fff;font-weight:600;font-size:11px}}
.tier-I{{background:#c62828}}.tier-II{{background:#ef6c00}}.tier-III{{background:#fbc02d;color:#333}}.tier-IV{{background:#2e7d32}}
.summary{{display:flex;gap:14px;margin-top:14px}}
.card{{flex:1;background:#fafafa;border-radius:6px;padding:14px;text-align:center}}
.card .v{{font-size:24px;font-weight:700;color:#7b1fa2}}
.card .l{{font-size:11px;text-transform:uppercase;color:#888;letter-spacing:1px}}
</style></head><body>
<h1>V2F Somatic Variant Report — {proband_id}</h1>
<p style='color:#666'>Framework: AMP/ASCO/CAP 2017 four-tier somatic interpretation guidelines.<br>
Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}</p>

<div class="summary">
  <div class="card"><div class="v">{n_total}</div><div class="l">Total variants</div></div>
  <div class="card"><div class="v">{n_t1}</div><div class="l">Tier I</div></div>
  <div class="card"><div class="v">{n_t2}</div><div class="l">Tier II</div></div>
</div>

<h2>Classified Somatic Variants</h2>
<table>
<thead><tr><th>Gene</th><th>Position</th><th>Protein change</th><th>AMP tier</th><th>VAF</th><th>Drug targets</th><th>Evidence</th></tr></thead>
<tbody>
{''.join(rows_html) if rows_html else '<tr><td colspan="7" style="text-align:center;color:#999;">No variants</td></tr>'}
</tbody>
</table>
{config.get('_run_manifest_footer_html', '')}
</body></html>"""
    with open(output_path, "w") as f:
        f.write(html)
    return output_path


def run(config: dict) -> dict:
    """Execute Stage 6: Report Generation.

    Routes per analysis_mode:
      'germline' / 'both' -> existing Jinja-based germline report
      'somatic'           -> inline minimal somatic report
                             (emphasizes AMP tiers + drug targets + VAF)
    """
    t0 = time.time()
    output_dir = config.get("output", {}).get("output_dir", "data")
    reports_dir = config.get("output", {}).get("reports_dir", "reports")
    class_dir = os.path.join(output_dir, "classified")
    pheno_dir = os.path.join(output_dir, "phenotype")
    os.makedirs(reports_dir, exist_ok=True)

    mode = (config.get("input", {}) or {}).get("analysis_mode", "germline")
    logger.info(f"Stage 6: Report Generation (mode={mode})")

    # Snapshot the run manifest once per pipeline invocation
    from pipeline.utils.run_manifest import write_manifest, manifest_footer_html
    manifest = write_manifest(config, reports_dir)
    # Make the rendered footer available to per-renderer code paths via config
    config = dict(config)
    config["_run_manifest_footer_html"] = manifest_footer_html(manifest)

    # Load classified variants
    class_path = os.path.join(class_dir, "acmg_results.parquet")
    if not os.path.exists(class_path):
        raise FileNotFoundError(f"Stage 4 output not found: {class_path}")

    df = pd.read_parquet(class_path)
    logger.info(f"Loaded {len(df)} classified variants")

    # Load phenotype if available
    pheno_path = os.path.join(pheno_dir, "patient_phenotype.json")
    phenotype = None
    if os.path.exists(pheno_path):
        with open(pheno_path) as f:
            phenotype = json.load(f)

    # Load QC metrics if available
    qc_path = os.path.join(output_dir, "vcf", "qc_report.json")
    qc_metrics = None
    if os.path.exists(qc_path):
        with open(qc_path) as f:
            qc_report = json.load(f)
            qc_metrics = qc_report.get("qc_metrics", {})

    # Generate report per proband
    reports_generated = []
    report_errors = []  # list of {sample_id, error_type, message}

    if "sample_id" in df.columns:
        probands = df["sample_id"].unique()
    else:
        probands = ["unknown"]

    for proband_id in probands:
        proband_df = df[df["sample_id"] == proband_id] if "sample_id" in df.columns else df

        output_path = os.path.join(reports_dir, f"{proband_id}_report.html")
        try:
            if mode == "somatic":
                report_path = _render_somatic_inline(
                    proband_df, proband_id, config, output_path)
            else:
                # germline or both — use the existing Jinja germline template
                report_path = render_proband_report(
                    variants_df=proband_df,
                    proband_id=proband_id,
                    config=config,
                    phenotype=phenotype,
                    qc_metrics=qc_metrics,
                    output_path=output_path,
                )
            reports_generated.append(report_path)
            logger.info(f"Report generated: {report_path}")
        except (FileNotFoundError, OSError) as e:
            err = {"sample_id": proband_id, "error_type": "io_error", "message": str(e)}
            report_errors.append(err)
            logger.error(f"Report I/O failed for {proband_id}: {e}")
        except (KeyError, ValueError, TypeError) as e:
            err = {"sample_id": proband_id, "error_type": "data_error", "message": str(e)}
            report_errors.append(err)
            logger.error(f"Report data error for {proband_id}: {e}")
        except ImportError as e:
            err = {"sample_id": proband_id, "error_type": "import_error", "message": str(e)}
            report_errors.append(err)
            logger.error(f"Report import error for {proband_id}: {e}")
        except Exception as e:
            # Unknown error — still capture so dashboard can surface it
            err = {"sample_id": proband_id, "error_type": "unknown",
                   "message": f"{type(e).__name__}: {e}"}
            report_errors.append(err)
            logger.exception(f"Report rendering failed for {proband_id}")

    result = {
        "stage": "06_report_generation",
        "reports_generated": len(reports_generated),
        "report_paths": reports_generated,
        "report_errors": report_errors,
        "n_errors": len(report_errors),
        "report_format": config.get("output", {}).get("report_format", "html"),
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    with open(os.path.join(reports_dir, "report_summary.json"), "w") as f:
        json.dump(result, f, indent=2, default=str)

    # Separate errors file the dashboard reads on /api/reports
    if report_errors:
        with open(os.path.join(reports_dir, "report_errors.json"), "w") as f:
            json.dump({"errors": report_errors}, f, indent=2)
    else:
        # Remove a stale errors file from a prior failed run so the
        # Reports tab doesn't keep showing a red banner after a fix.
        stale = os.path.join(reports_dir, "report_errors.json")
        if os.path.exists(stale):
            os.unlink(stale)

    print(f"\n{'='*60}")
    print(f"Stage 6 Complete: Report Generation")
    print(f"{'='*60}")
    print(f"  Reports generated: {len(reports_generated)}")
    for rp in reports_generated:
        print(f"    → {rp}")
    print(f"  Time:              {result['elapsed_seconds']}s")
    print(f"{'='*60}\n")

    return result


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    run(config)
