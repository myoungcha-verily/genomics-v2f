"""Pre-M1 polish tests: stage 6 writes report_errors.json on failure and
/api/reports surfaces them."""

import json
import os
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest


def _write_classified_parquet(tmp_path):
    df = pd.DataFrame([
        {"sample_id": "S1", "gene": "BRCA1", "acmg_classification": "Pathogenic"},
        {"sample_id": "S2", "gene": "TP53",  "acmg_classification": "VUS"},
    ])
    out = tmp_path / "data" / "classified" / "acmg_results.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    return out


def test_stage6_writes_report_errors_when_render_fails(repo_root, tmp_path,
                                                       monkeypatch):
    """If render_proband_report raises, stage 6 must capture the error per-sample
    in report_errors.json instead of silently producing zero reports."""
    monkeypatch.chdir(tmp_path)
    _write_classified_parquet(tmp_path)
    (tmp_path / "data" / "phenotype").mkdir(parents=True, exist_ok=True)
    (tmp_path / "data" / "phenotype" / "patient_phenotype.json").write_text("{}")

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "pipeline.06_report_generation",
        repo_root / "pipeline" / "06_report_generation.py",
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    # Patch render_proband_report to always raise
    def _raise(*a, **kw):
        raise KeyError("missing template var 'gene_panel'")
    monkeypatch.setattr(mod, "render_proband_report", _raise)

    config = {
        "output": {"reports_dir": "reports", "report_format": "html"},
    }
    result = mod.run(config)

    # Stage 6 should still complete (return result) but with errors recorded
    assert result["reports_generated"] == 0
    assert result["n_errors"] >= 1
    assert any(e["error_type"] == "data_error"
               for e in result.get("report_errors", []))

    # report_errors.json should be on disk
    err_path = tmp_path / "reports" / "report_errors.json"
    assert err_path.exists(), "report_errors.json missing after failed render"
    payload = json.loads(err_path.read_text())
    assert "errors" in payload and len(payload["errors"]) >= 1


def test_dashboard_reports_endpoint_surfaces_errors(repo_root, tmp_path,
                                                     monkeypatch):
    """/api/reports must return errors[] when reports/report_errors.json exists."""
    reports_dir = tmp_path / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    (reports_dir / "report_errors.json").write_text(json.dumps({
        "errors": [
            {"sample_id": "S1", "error_type": "data_error",
             "message": "missing column gene_panel"},
        ]
    }))

    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()
    rv = client.get("/api/reports")
    assert rv.status_code == 200
    body = rv.get_json()
    assert "errors" in body
    assert len(body["errors"]) == 1
    assert body["errors"][0]["sample_id"] == "S1"


def test_dashboard_reports_endpoint_clean_when_no_errors_file(repo_root, tmp_path,
                                                               monkeypatch):
    """/api/reports returns empty errors[] when no report_errors.json exists."""
    reports_dir = tmp_path / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()
    rv = client.get("/api/reports")
    body = rv.get_json()
    assert body.get("errors", []) == []
