"""M3.4: bundled SV demo data parses + classifies as expected."""

import json
import pandas as pd
import pytest


def test_sv_demo_vcf_present(repo_root):
    vcf = repo_root / "demo" / "sv" / "sv_demo.vcf.gz"
    assert vcf.exists()
    assert vcf.stat().st_size < 5000


def test_sv_demo_config_yaml(repo_root):
    cfg = (repo_root / "demo" / "sv" / "config.yaml").read_text()
    assert "analysis_mode: cnv" in cfg
    assert "adapter: sv" in cfg


def test_sv_precomputed_parquet_has_cnv_classifications(repo_root):
    pq = (repo_root / "demo" / "precomputed_sv"
          / "data" / "classified" / "acmg_results.parquet")
    df = pd.read_parquet(pq)
    assert "cnv_classification" in df.columns
    assert "svtype" in df.columns
    assert "cnv_score" in df.columns
    # 22q11.2 deletion should be Pathogenic
    twenty_two = df[df.gene == "22q11.2"]
    assert len(twenty_two) == 1
    assert twenty_two.iloc[0]["cnv_classification"] == "Pathogenic"


def test_dashboard_demo_load_cnv_mode(repo_root, tmp_path, monkeypatch):
    """POST /api/demo/load with mode=cnv populates SV data tree."""
    import shutil

    shutil.copytree(repo_root / "demo", tmp_path / "demo")
    (tmp_path / "config").mkdir()

    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()
    rv = client.post("/api/demo/load", json={"mode": "cnv"})
    assert rv.status_code == 200
    body = rv.get_json()
    assert body["loaded"] is True
    assert body["mode"] == "cnv"
    assert body["variants"] == 6
    assert (tmp_path / "data" / "vcf" / "sv_demo.vcf.gz").exists()
    cfg = (tmp_path / "config" / "pipeline_config.yaml").read_text()
    assert "analysis_mode: cnv" in cfg

    rv2 = client.get("/api/reports")
    sample_ids = [r["sample_id"] for r in (rv2.get_json().get("reports") or [])]
    assert "DEMO_SV" in sample_ids
