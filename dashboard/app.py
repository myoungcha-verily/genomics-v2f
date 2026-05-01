#!/usr/bin/env python3
"""Variant-to-Function Reporter — Dashboard

Flask app serving on port 8080 with 6 tabs:
  Overview, Upload, Pipeline, Variants, Reports, Settings

All fetch() calls use RELATIVE paths for Workbench proxy compatibility.
"""

import json
import logging
import os
import shutil
import subprocess
import sys
import threading
import time

import pandas as pd
import yaml
from flask import Flask, jsonify, request, render_template, send_from_directory

app = Flask(__name__)
logger = logging.getLogger(__name__)

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_PATH = os.path.join(PROJECT_DIR, "config", "pipeline_config.yaml")
DEMO_DIR = os.path.join(PROJECT_DIR, "demo")

# Make pipeline.* importable regardless of cwd
if PROJECT_DIR not in sys.path:
    sys.path.insert(0, PROJECT_DIR)

# Pipeline process tracking
pipeline_process = None
pipeline_log = []


def _load_config():
    """Load pipeline config."""
    try:
        with open(CONFIG_PATH) as f:
            return yaml.safe_load(f)
    except Exception:
        return {}


def _check_stage(stage_num):
    """Check if a pipeline stage has completed."""
    indicators = {
        1: "data/vcf/variants.parquet",
        2: "data/annotated/annotated_variants.parquet",
        3: "data/enriched/variants_enriched.parquet",
        4: "data/classified/acmg_results.parquet",
        5: "data/phenotype/patient_phenotype.json",
        6: "reports/",
        7: "eval/validation_summary.json",
    }
    path = os.path.join(PROJECT_DIR, indicators.get(stage_num, ""))
    if os.path.isdir(path):
        return len(os.listdir(path)) > 0
    return os.path.exists(path)


def _get_artifact_info(path):
    """Get file/directory info."""
    full = os.path.join(PROJECT_DIR, path)
    if os.path.isdir(full):
        files = os.listdir(full)
        return {"type": "directory", "files": len(files), "path": path}
    elif os.path.exists(full):
        size = os.path.getsize(full)
        return {"type": "file", "size_kb": round(size / 1024, 1), "path": path}
    return None


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/status")
def api_status():
    """Pipeline status — which stages are complete."""
    stages = []
    stage_names = {
        1: "VCF Ingest & QC",
        2: "Variant Annotation",
        3: "Database Enrichment",
        4: "ACMG Classification",
        5: "Phenotype Integration",
        6: "Report Generation",
        7: "Validation",
    }
    for i in range(1, 8):
        complete = _check_stage(i)
        stages.append({
            "stage": i,
            "name": stage_names[i],
            "complete": complete,
        })

    n_complete = sum(1 for s in stages if s["complete"])
    running = pipeline_process is not None and pipeline_process.poll() is None

    return jsonify({
        "stages": stages,
        "completed": n_complete,
        "total": 7,
        "running": running,
    })


@app.route("/api/databases")
def api_databases():
    """Check database connectivity."""
    config = _load_config()
    databases = []

    clinvar = config.get("databases", {}).get("clinvar", {})
    if clinvar.get("bq_table"):
        databases.append({
            "name": "ClinVar",
            "table": clinvar["bq_table"],
            "status": "configured",
        })

    gnomad = config.get("databases", {}).get("gnomad", {})
    if gnomad.get("bq_table"):
        databases.append({
            "name": "gnomAD",
            "table": gnomad["bq_table"],
            "status": "configured",
        })

    return jsonify({"databases": databases})


@app.route("/api/config", methods=["GET"])
def api_get_config():
    """Return current config."""
    return jsonify(_load_config())


@app.route("/api/config", methods=["POST"])
def api_save_config():
    """Save updated config."""
    try:
        new_config = request.json
        os.makedirs(os.path.dirname(CONFIG_PATH), exist_ok=True)
        with open(CONFIG_PATH, "w") as f:
            yaml.dump(new_config, f, default_flow_style=False, sort_keys=False)
        return jsonify({"status": "ok"})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/upload", methods=["POST"])
def api_upload():
    """Handle VCF file upload."""
    if "file" not in request.files:
        return jsonify({"error": "No file provided"}), 400

    file = request.files["file"]
    if not file.filename:
        return jsonify({"error": "No filename"}), 400

    # Save to data/vcf/
    vcf_dir = os.path.join(PROJECT_DIR, "data", "vcf")
    os.makedirs(vcf_dir, exist_ok=True)
    save_path = os.path.join(vcf_dir, file.filename)
    file.save(save_path)

    # Update config
    config = _load_config()
    config.setdefault("input", {})["vcf_path"] = save_path
    with open(CONFIG_PATH, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    return jsonify({
        "status": "ok",
        "path": save_path,
        "size_mb": round(os.path.getsize(save_path) / 1024 / 1024, 2),
    })


@app.route("/api/pipeline/run", methods=["POST"])
def api_pipeline_run():
    """Start pipeline execution."""
    global pipeline_process, pipeline_log

    if pipeline_process is not None and pipeline_process.poll() is None:
        return jsonify({"error": "Pipeline already running"}), 409

    pipeline_log = []
    cmd = ["python3", os.path.join(PROJECT_DIR, "run_pipeline.py"),
           "--config", CONFIG_PATH]

    # Optional: specific stages
    stages = request.json.get("stages") if request.json else None
    if stages:
        cmd.extend(["--stages", stages])

    pipeline_process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, cwd=PROJECT_DIR,
    )

    # Stream output in background
    def _read_output():
        for line in pipeline_process.stdout:
            pipeline_log.append(line.rstrip())
            if len(pipeline_log) > 1000:
                pipeline_log.pop(0)

    thread = threading.Thread(target=_read_output, daemon=True)
    thread.start()

    return jsonify({"status": "started", "pid": pipeline_process.pid})


@app.route("/api/pipeline/status")
def api_pipeline_status():
    """Get pipeline execution status and logs."""
    running = pipeline_process is not None and pipeline_process.poll() is None
    return jsonify({
        "running": running,
        "returncode": pipeline_process.returncode if pipeline_process and not running else None,
        "log": pipeline_log[-100:],  # Last 100 lines
    })


@app.route("/api/pipeline/cancel", methods=["POST"])
def api_pipeline_cancel():
    """Cancel running pipeline."""
    global pipeline_process
    if pipeline_process and pipeline_process.poll() is None:
        pipeline_process.terminate()
        return jsonify({"status": "cancelled"})
    return jsonify({"status": "not_running"})


@app.route("/api/variants")
def api_variants():
    """Return classified variants for interactive table."""
    parquet_path = os.path.join(PROJECT_DIR, "data", "classified",
                                "acmg_results.parquet")
    if not os.path.exists(parquet_path):
        return jsonify({"variants": [], "total": 0})

    df = pd.read_parquet(parquet_path)

    # Apply filters from query params
    classification = request.args.get("classification")
    gene = request.args.get("gene")
    min_cadd = request.args.get("min_cadd", type=float)

    if classification:
        df = df[df["acmg_classification"] == classification]
    if gene:
        df = df[df["gene"].str.contains(gene, case=False, na=False)]
    if min_cadd is not None and "cadd_phred" in df.columns:
        df = df[df["cadd_phred"] >= min_cadd]

    # Pagination
    page = request.args.get("page", 1, type=int)
    per_page = request.args.get("per_page", 50, type=int)
    total = len(df)
    df = df.iloc[(page-1)*per_page : page*per_page]

    # Convert to records
    cols = ["variant_id", "chrom", "pos", "ref", "alt", "gene",
            "consequence", "severity", "genotype", "acmg_classification",
            "acmg_criteria", "gnomad_af", "clinvar_classification",
            "cadd_phred", "revel", "spliceai_max", "read_depth",
            "allele_fraction"]
    cols = [c for c in cols if c in df.columns]
    records = df[cols].fillna("").to_dict("records")

    return jsonify({
        "variants": records,
        "total": total,
        "page": page,
        "per_page": per_page,
    })


@app.route("/api/reports/<sample_id>")
def api_report(sample_id):
    """Serve a generated report."""
    reports_dir = os.path.join(PROJECT_DIR, "reports")
    filename = f"{sample_id}_report.html"
    if os.path.exists(os.path.join(reports_dir, filename)):
        return send_from_directory(reports_dir, filename)
    return jsonify({"error": "Report not found"}), 404


@app.route("/api/reports")
def api_list_reports():
    """List generated reports.

    Also surfaces any per-sample report-generation errors written by stage 6
    to reports/report_errors.json so the Reports tab can warn the user when
    a run silently produced fewer reports than expected.
    """
    reports_dir = os.path.join(PROJECT_DIR, "reports")
    if not os.path.exists(reports_dir):
        return jsonify({"reports": [], "errors": []})

    reports = []
    for f in os.listdir(reports_dir):
        if f.endswith("_report.html"):
            sample = f.replace("_report.html", "")
            path = os.path.join(reports_dir, f)
            reports.append({
                "sample_id": sample,
                "filename": f,
                "size_kb": round(os.path.getsize(path) / 1024, 1),
                "modified": os.path.getmtime(path),
            })

    errors = []
    err_path = os.path.join(reports_dir, "report_errors.json")
    if os.path.exists(err_path):
        try:
            with open(err_path) as fh:
                errors = (json.load(fh) or {}).get("errors", [])
        except Exception:
            pass

    return jsonify({"reports": reports, "errors": errors})


@app.route("/api/gene_panels")
def api_gene_panels():
    """List available gene panels."""
    panel_dir = os.path.join(PROJECT_DIR, "reference", "gene_panels")
    panels = []
    if os.path.exists(panel_dir):
        for f in os.listdir(panel_dir):
            if f.endswith(".json") and f != "custom_template.json":
                with open(os.path.join(panel_dir, f)) as fp:
                    data = json.load(fp)
                panels.append({
                    "name": f.replace(".json", ""),
                    "display_name": data.get("panel_name", f),
                    "gene_count": data.get("gene_count", 0),
                })
    return jsonify({"panels": panels})


@app.route("/api/phenotype/test", methods=["POST"])
def api_phenotype_test():
    """Probe whether the configured FHIR project + dataset are reachable.

    Body (optional): {fhir_project, fhir_dataset, condition_table}.
    If body is empty, reads values from the saved config.
    """
    body = request.get_json(silent=True) or {}
    if "fhir_project" not in body or "fhir_dataset" not in body:
        cfg = _load_config()
        pheno = cfg.get("phenotype", {}) if cfg else {}
        body.setdefault("fhir_project", pheno.get("fhir_project", ""))
        body.setdefault("fhir_dataset", pheno.get("fhir_dataset", ""))
        body.setdefault("condition_table", pheno.get("condition_table", "Condition"))

    try:
        from pipeline.utils.fhir_phenotype import test_fhir_connectivity
    except Exception as e:
        return jsonify({"ok": False, "error": f"import failed: {e}"}), 500

    result = test_fhir_connectivity(
        body.get("fhir_project", ""),
        body.get("fhir_dataset", ""),
        body.get("condition_table", "Condition"),
    )
    return jsonify(result)


@app.route("/api/validation")
def api_validation():
    """Aggregate quality-gate results from each stage's summary JSON.

    Returns:
        {
          "stages": [
            {"stage": "stage_1", "gates": [{"name", "passed", ...}], ...},
            {"stage": "stage_3", "gates": [...], ...},
            {"stage": "stage_4", "gates": [...], ...}
          ],
          "overall_passed": bool,
          "n_warnings": int,
          "n_hard_fails": int,
        }
    """
    summaries = [
        ("stage_1", "data/vcf/qc_report.json"),
        ("stage_3", "data/enriched/enrichment_summary.json"),
        ("stage_4", "data/classified/classification_summary.json"),
    ]
    out = {"stages": [], "overall_passed": True,
           "n_warnings": 0, "n_hard_fails": 0}
    for stage, rel_path in summaries:
        path = os.path.join(PROJECT_DIR, rel_path)
        if not os.path.exists(path):
            continue
        try:
            with open(path) as f:
                data = json.load(f)
        except Exception:
            continue
        gates = data.get("quality_gates", []) or []
        if not isinstance(gates, list):
            gates = []
        for g in gates:
            if not g.get("passed"):
                if g.get("severity") == "hard_fail":
                    out["n_hard_fails"] += 1
                    out["overall_passed"] = False
                else:
                    out["n_warnings"] += 1
        out["stages"].append({
            "stage": stage,
            "label": data.get("stage", stage),
            "gates": gates,
            "summary": {k: data.get(k) for k in
                         ("total_variants", "total_variants_filtered",
                          "clinvar_matches", "gnomad_matches",
                          "pathogenic", "likely_pathogenic", "vus",
                          "likely_benign", "benign")
                         if k in data},
        })
    return jsonify(out)


@app.route("/api/curation/<variant_id>", methods=["GET"])
def api_curation_list(variant_id):
    """List curation entries for a variant."""
    try:
        from pipeline.utils.curation_store import list_curations
    except Exception as e:
        return jsonify({"error": f"import failed: {e}"}), 500
    cfg = _load_config()
    entries = list_curations(variant_id, cfg)
    return jsonify({"variant_id": variant_id, "entries": entries})


@app.route("/api/curation/<variant_id>", methods=["POST"])
def api_curation_add(variant_id):
    """Add a curation entry for a variant."""
    try:
        from pipeline.utils.curation_store import add_curation
    except Exception as e:
        return jsonify({"error": f"import failed: {e}"}), 500
    body = request.get_json(silent=True) or {}
    body["variant_id"] = variant_id
    cfg = _load_config()
    try:
        stored = add_curation(body, cfg)
        return jsonify({"stored": True, "entry": stored})
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": f"unexpected: {e}"}), 500


@app.route("/api/literature/<variant_id>", methods=["GET"])
def api_literature(variant_id):
    """Look up literature for a variant via LitVar2.

    Variant id is chrom-pos-ref-alt; the gene + hgvs_p are taken from the
    classified parquet so the caller doesn't have to resend them.
    """
    pq = os.path.join(PROJECT_DIR, "data", "classified", "acmg_results.parquet")
    if not os.path.exists(pq):
        return jsonify({"error": "no classified variants — run pipeline first"}), 404

    df = pd.read_parquet(pq)
    matches = df[df["variant_id"] == variant_id]
    if matches.empty:
        return jsonify({"error": f"variant {variant_id} not found"}), 404
    row = matches.iloc[0]
    gene = row.get("gene", "")
    hgvs_p = row.get("hgvs_p", "")

    try:
        from pipeline.utils.litvar_client import query_litvar
    except Exception as e:
        return jsonify({"error": f"import failed: {e}"}), 500
    return jsonify(query_litvar(gene, hgvs_p))


@app.route("/api/demo/load", methods=["POST"])
def api_demo_load():
    """Copy bundled demo VCF + precomputed outputs into the live data dirs.

    Body (optional): {"mode": "germline" | "somatic"}. Defaults to germline.
    """
    body = request.get_json(silent=True) or {}
    mode = body.get("mode", "germline")
    if mode not in ("germline", "somatic", "cnv"):
        return jsonify({"error": f"Unknown demo mode: {mode}"}), 400

    if mode == "germline":
        src_demo_dir = os.path.join(DEMO_DIR, "germline")
        precomputed_dir = os.path.join(DEMO_DIR, "precomputed")
        vcf_basename = "proband.vcf.gz"
    elif mode == "somatic":
        src_demo_dir = os.path.join(DEMO_DIR, "somatic")
        precomputed_dir = os.path.join(DEMO_DIR, "precomputed_somatic")
        vcf_basename = "tumor_normal.vcf.gz"
    else:  # cnv
        src_demo_dir = os.path.join(DEMO_DIR, "sv")
        precomputed_dir = os.path.join(DEMO_DIR, "precomputed_sv")
        vcf_basename = "sv_demo.vcf.gz"

    if not os.path.isdir(src_demo_dir) or not os.path.isdir(precomputed_dir):
        return jsonify({"error": f"{mode} demo data not bundled in this build"}), 500

    # 1. Copy demo VCF + tabix index into data/vcf/
    vcf_dst_dir = os.path.join(PROJECT_DIR, "data", "vcf")
    os.makedirs(vcf_dst_dir, exist_ok=True)
    for fname in (vcf_basename, vcf_basename + ".tbi"):
        src = os.path.join(src_demo_dir, fname)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(vcf_dst_dir, fname))

    # 2. Install demo config
    demo_config = os.path.join(src_demo_dir, "config.yaml")
    if os.path.exists(demo_config):
        os.makedirs(os.path.dirname(CONFIG_PATH), exist_ok=True)
        shutil.copy2(demo_config, CONFIG_PATH)

    # 3. Copy precomputed pipeline outputs into the live data/ + reports/ + eval/ trees
    for sub in ("data", "reports", "eval"):
        src = os.path.join(precomputed_dir, sub)
        if os.path.isdir(src):
            dst = os.path.join(PROJECT_DIR, sub)
            for root, _, files in os.walk(src):
                rel = os.path.relpath(root, src)
                tgt = os.path.join(dst, rel) if rel != "." else dst
                os.makedirs(tgt, exist_ok=True)
                for f in files:
                    shutil.copy2(os.path.join(root, f), os.path.join(tgt, f))

    n_variants = 0
    pq = os.path.join(PROJECT_DIR, "data", "classified", "acmg_results.parquet")
    if os.path.exists(pq):
        try:
            n_variants = len(pd.read_parquet(pq))
        except Exception:
            n_variants = 0

    return jsonify({
        "loaded": True,
        "mode": mode,
        "vcf": f"data/vcf/{vcf_basename}",
        "config": "config/pipeline_config.yaml",
        "variants": n_variants,
    })


@app.route("/api/demo/reset", methods=["POST"])
def api_demo_reset():
    """Wipe data/, reports/, eval/, logs/ and pipeline_config.yaml — reverts
    to first-run state for repeatable demos."""
    for sub in ("data", "reports", "eval", "logs"):
        path = os.path.join(PROJECT_DIR, sub)
        if os.path.isdir(path):
            shutil.rmtree(path)
    if os.path.exists(CONFIG_PATH):
        os.unlink(CONFIG_PATH)
    return jsonify({"reset": True})


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    app.run(host="0.0.0.0", port=8080, debug=False)
