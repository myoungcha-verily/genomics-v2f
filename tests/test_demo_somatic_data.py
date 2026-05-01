"""M1: bundled somatic demo files exist and parse correctly."""

import json
from pathlib import Path

import pandas as pd
import pytest


def test_somatic_demo_vcf_present(repo_root):
    vcf = repo_root / "demo" / "somatic" / "tumor_normal.vcf.gz"
    assert vcf.exists()
    assert vcf.stat().st_size < 5000


def test_somatic_demo_config_yaml(repo_root):
    cfg_path = repo_root / "demo" / "somatic" / "config.yaml"
    assert cfg_path.exists()
    text = cfg_path.read_text()
    assert "analysis_mode: somatic" in text
    assert "DEMO_TUMOR" in text
    assert "DEMO_NORMAL" in text


def test_somatic_precomputed_parquet_has_amp_tier(repo_root):
    pq = (repo_root / "demo" / "precomputed_somatic"
          / "data" / "classified" / "acmg_results.parquet")
    assert pq.exists()
    df = pd.read_parquet(pq)
    assert "amp_tier" in df.columns
    assert "amp_drug_targets" in df.columns
    # BRAF V600E should be Tier I in the curated demo
    braf = df[df.gene == "BRAF"]
    assert len(braf) >= 1
    assert braf.iloc[0]["amp_tier"] == "I"
    assert "Vemurafenib" in braf.iloc[0]["amp_drug_targets"]


def test_dashboard_demo_load_somatic_mode(repo_root, tmp_path, monkeypatch):
    """POST /api/demo/load with mode=somatic populates somatic data tree."""
    import shutil

    shutil.copytree(repo_root / "demo", tmp_path / "demo")
    (tmp_path / "config").mkdir()

    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()
    rv = client.post("/api/demo/load", json={"mode": "somatic"})
    assert rv.status_code == 200
    body = rv.get_json()
    assert body["loaded"] is True
    assert body["mode"] == "somatic"
    assert body["variants"] >= 5
    assert (tmp_path / "data" / "vcf" / "tumor_normal.vcf.gz").exists()
    cfg = (tmp_path / "config" / "pipeline_config.yaml").read_text()
    assert "analysis_mode: somatic" in cfg

    # Reports endpoint should now show DEMO_TUMOR
    rv2 = client.get("/api/reports")
    assert rv2.status_code == 200
    sample_ids = [r["sample_id"] for r in (rv2.get_json().get("reports") or [])]
    assert "DEMO_TUMOR" in sample_ids
