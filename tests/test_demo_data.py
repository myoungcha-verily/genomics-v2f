"""Milestone 0 tests: bundled demo data parses and the dashboard demo
endpoints behave correctly."""

import gzip
import json
import os
from pathlib import Path

import pandas as pd
import pytest


def test_demo_vcf_present(repo_root):
    vcf = repo_root / "demo" / "germline" / "proband.vcf.gz"
    tbi = repo_root / "demo" / "germline" / "proband.vcf.gz.tbi"
    assert vcf.exists(), "demo VCF should be checked in"
    assert tbi.exists(), "demo VCF tabix index should be checked in"
    assert vcf.stat().st_size < 5000, "demo VCF should stay small (<5 KB)"


def test_demo_vcf_parses(repo_root):
    """VCF can be read end-to-end and contains 10 records."""
    vcf = repo_root / "demo" / "germline" / "proband.vcf.gz"
    n_records = 0
    n_header = 0
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                n_header += 1
            elif line.startswith("#CHROM"):
                pass
            elif line.strip():
                n_records += 1
    assert n_header >= 10, "VCF header looks malformed"
    assert n_records == 10, f"expected 10 demo variants, got {n_records}"


def test_demo_phenotype_payload_shape(repo_root):
    pheno = repo_root / "demo" / "germline" / "phenotype.json"
    assert pheno.exists()
    payload = json.loads(pheno.read_text())
    assert payload["patient_id"] == "DEMO_PROBAND"
    assert any("BRCA1" in g for g in payload["candidate_genes"])
    assert len(payload["conditions"]) >= 1


def test_demo_config_yaml(repo_root):
    cfg_path = repo_root / "demo" / "germline" / "config.yaml"
    assert cfg_path.exists()
    text = cfg_path.read_text()
    # Spot-check that the bundled config points at the bundled VCF
    assert "data/vcf/proband.vcf.gz" in text
    assert "single_sample" in text


def test_precomputed_parquets_present(repo_root):
    base = repo_root / "demo" / "precomputed"
    expected = [
        "data/vcf/variants.parquet",
        "data/annotated/annotated_variants.parquet",
        "data/enriched/variants_enriched.parquet",
        "data/classified/acmg_results.parquet",
        "data/phenotype/patient_phenotype.json",
        "reports/DEMO_PROBAND_report.html",
        "eval/validation_summary.json",
    ]
    for rel in expected:
        assert (base / rel).exists(), f"precomputed file missing: {rel}"


def test_precomputed_classifications_match_ground_truth(repo_root):
    """Every demo variant gets the expected ACMG classification.

    This is the regression guard: if classification logic changes (e.g.
    Milestone 2 calibrated tiers), the precomputed outputs must be
    regenerated and this test re-pinned.
    """
    pq = repo_root / "demo" / "precomputed" / "data" / "classified" / "acmg_results.parquet"
    df = pd.read_parquet(pq)
    by_gene = {(r.gene, r.pos): r.acmg_classification for _, r in df.iterrows()}

    # Lookup table from demo/README.md ground truth
    expected = {
        ("BRCA1", 43106478): "Pathogenic",
        ("TP53", 7674220): "Pathogenic",
        ("MYH7", 23427594): "Likely pathogenic",
        ("PTEN", 87864458): "VUS",
        ("BRCA2", 32340301): "Likely benign",
    }
    for key, want in expected.items():
        assert by_gene.get(key) == want, (
            f"{key} expected {want}, got {by_gene.get(key)} — "
            "regenerate demo/precomputed/ if classification logic changed"
        )


def test_precomputed_phenotype_boost_for_hboc_genes(repo_root):
    """Phenotype scoring should put BRCA1 pathogenic at the top."""
    pq = repo_root / "demo" / "precomputed" / "data" / "classified" / "acmg_results.parquet"
    df = pd.read_parquet(pq)
    brca1_path = df[(df.gene == "BRCA1") & (df.acmg_classification == "Pathogenic")]
    assert len(brca1_path) == 1
    assert brca1_path.iloc[0].phenotype_match_score > 0.8


def test_dashboard_demo_endpoints_load_then_reset(repo_root, tmp_path, monkeypatch):
    """End-to-end test of /api/demo/load and /api/demo/reset.

    Patches dashboard.app's module-level path constants to a temp tree so
    we don't clobber the developer's real work tree.
    """
    import shutil

    shutil.copytree(repo_root / "demo", tmp_path / "demo")
    (tmp_path / "config").mkdir()

    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()

    # 1. Load
    rv = client.post("/api/demo/load")
    assert rv.status_code == 200
    body = rv.get_json()
    assert body["loaded"] is True
    assert body["variants"] == 10
    assert (tmp_path / "data" / "classified" / "acmg_results.parquet").exists()
    assert (tmp_path / "config" / "pipeline_config.yaml").exists()

    # 2. Variants endpoint should now return 10
    rv = client.get("/api/variants")
    assert rv.status_code == 200
    assert rv.get_json()["total"] == 10

    # 3. Reset
    rv = client.post("/api/demo/reset")
    assert rv.status_code == 200
    assert rv.get_json()["reset"] is True
    assert not (tmp_path / "data").exists()
    assert not (tmp_path / "config" / "pipeline_config.yaml").exists()
