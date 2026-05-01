"""P6.2: /api/validation aggregates per-stage gate results."""

import json
from pathlib import Path


def _seed_summary(tmp_path, stage_dir, fname, content):
    """Write a stage summary JSON with embedded gate results."""
    out = tmp_path / "data" / stage_dir / fname
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(content))


def test_validation_endpoint_aggregates_stage_gates(tmp_path, monkeypatch):
    _seed_summary(tmp_path, "vcf", "qc_report.json", {
        "stage": "01_vcf_ingest_qc",
        "total_variants_filtered": 100,
        "quality_gates": [
            {"name": "stage_1.min_variants", "passed": True,
             "severity": "warning", "actual": 100, "expected": ">= 1",
             "message": "ok"},
            {"name": "stage_1.ti_tv_ratio", "passed": False,
             "severity": "warning", "actual": 0.5, "expected": "in [1.5, 3.5]",
             "message": "Ti/Tv low"},
        ],
    })
    _seed_summary(tmp_path, "enriched", "enrichment_summary.json", {
        "stage": "03_database_enrichment",
        "total_variants": 100, "clinvar_matches": 0,
        "quality_gates": [
            {"name": "stage_3.clinvar_match_rate", "passed": False,
             "severity": "hard_fail", "actual": "0.0%", "expected": ">= 30%",
             "message": "Likely BQ misconfig"},
        ],
    })

    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()
    rv = client.get("/api/validation")
    assert rv.status_code == 200
    body = rv.get_json()
    assert body["overall_passed"] is False
    assert body["n_hard_fails"] == 1
    assert body["n_warnings"] == 1
    assert len(body["stages"]) == 2
    stage_names = [s["stage"] for s in body["stages"]]
    assert "stage_1" in stage_names
    assert "stage_3" in stage_names


def test_validation_endpoint_no_runs_yet(tmp_path, monkeypatch):
    """Empty data tree → empty stages list, overall_passed True (vacuous)."""
    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()
    rv = client.get("/api/validation")
    body = rv.get_json()
    assert body["stages"] == []
    assert body["overall_passed"] is True
    assert body["n_hard_fails"] == 0


def test_validation_endpoint_all_pass(tmp_path, monkeypatch):
    _seed_summary(tmp_path, "classified", "classification_summary.json", {
        "stage": "04_acmg_classification",
        "pathogenic": 5, "vus": 20, "benign": 30,
        "quality_gates": [
            {"name": "stage_4.vus_fraction", "passed": True,
             "severity": "warning", "actual": "36.4%",
             "expected": "<= 95%", "message": "ok"},
        ],
    })

    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()
    body = client.get("/api/validation").get_json()
    assert body["overall_passed"] is True
    assert body["n_hard_fails"] == 0
    assert body["n_warnings"] == 0
