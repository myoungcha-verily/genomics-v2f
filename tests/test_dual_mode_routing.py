"""M1: pipeline_config.yaml input.analysis_mode routes stage 4 + 5 + 6."""

import json
import os
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest


def _seed_enriched_parquet(tmp_path):
    df = pd.DataFrame([
        {"variant_id": "chr7_140753336_A_T", "chrom": "chr7", "pos": 140753336,
         "ref": "A", "alt": "T", "sample_id": "S1",
         "gene": "BRAF", "hgvs_p": "p.V600E",
         "consequence": "missense_variant", "severity": "MODERATE",
         "tumor_vaf": 0.42, "is_paired": True,
         "gnomad_af": 0.0, "clinvar_classification": "",
         "clinvar_review_stars": 0, "cadd_phred": 28.0, "revel": 0.85,
         "spliceai_max": 0.0, "alphamissense": 0.91,
         "genotype": "0/1", "read_depth": 80,
         "allele_fraction": 0.42, "allele_depth_alt": 33,
         "allele_depth_ref": 47, "qual": 99.0, "filter": "PASS",
         "gt_quality": 99.0},
    ])
    out = tmp_path / "data" / "enriched" / "variants_enriched.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)


def _load_stage4(repo_root):
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "pipeline.stage4",
        repo_root / "pipeline" / "04_acmg_classification.py",
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_germline_mode_populates_acmg_classification(tmp_path, repo_root,
                                                      monkeypatch):
    monkeypatch.chdir(tmp_path)
    _seed_enriched_parquet(tmp_path)
    stage4 = _load_stage4(repo_root)
    config = {
        "input": {"analysis_mode": "germline"},
        "output": {"output_dir": "data"},
        "acmg": {},
    }
    result = stage4.run(config)
    out = pd.read_parquet(tmp_path / "data" / "classified" / "acmg_results.parquet")
    assert "acmg_classification" in out.columns
    assert "amp_tier" not in out.columns  # somatic-only column
    assert result["analysis_mode"] == "germline"


def test_somatic_mode_populates_amp_tier(tmp_path, repo_root, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _seed_enriched_parquet(tmp_path)
    stage4 = _load_stage4(repo_root)

    # Mock CIViC to return Tier I for BRAF V600E
    with patch("pipeline.utils.civic_client.is_enabled", return_value=True), \
         patch("pipeline.utils.oncokb_client.is_enabled", return_value=False), \
         patch("pipeline.utils.cosmic_client.is_enabled", return_value=False), \
         patch("pipeline.utils.civic_client.query_civic") as q_civic:
        q_civic.return_value = {
            "found": True, "highest_level": "A", "amp_tier": "I",
            "drug_targets": ["Vemurafenib"], "evidence_items": [],
            "source_url": None, "error": None, "latency_ms": 50,
        }
        config = {
            "input": {"analysis_mode": "somatic"},
            "output": {"output_dir": "data"},
            "databases": {"civic": {"enabled": True}},
            "acmg_amp": {"knowledge_base_priority": ["civic"]},
        }
        result = stage4.run(config)

    out = pd.read_parquet(tmp_path / "data" / "classified" / "acmg_results.parquet")
    assert "amp_tier" in out.columns
    assert out.iloc[0]["amp_tier"] == "I"
    # acmg_classification gets a fallback "AMP Tier I" string in somatic mode
    assert "AMP" in str(out.iloc[0]["acmg_classification"])
    assert result["amp_tier_i"] == 1


def test_both_mode_populates_both_columns(tmp_path, repo_root, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _seed_enriched_parquet(tmp_path)
    stage4 = _load_stage4(repo_root)

    with patch("pipeline.utils.civic_client.is_enabled", return_value=True), \
         patch("pipeline.utils.oncokb_client.is_enabled", return_value=False), \
         patch("pipeline.utils.cosmic_client.is_enabled", return_value=False), \
         patch("pipeline.utils.civic_client.query_civic") as q_civic:
        q_civic.return_value = {
            "found": True, "highest_level": "A", "amp_tier": "I",
            "drug_targets": ["Vemurafenib"], "evidence_items": [],
            "source_url": None, "error": None, "latency_ms": 50,
        }
        config = {
            "input": {"analysis_mode": "both"},
            "output": {"output_dir": "data"},
            "databases": {"civic": {"enabled": True}},
            "acmg_amp": {"knowledge_base_priority": ["civic"]},
            "acmg": {},
        }
        result = stage4.run(config)

    out = pd.read_parquet(tmp_path / "data" / "classified" / "acmg_results.parquet")
    assert "acmg_classification" in out.columns
    assert "amp_tier" in out.columns
    # ACMG classification is a real string (not 'AMP Tier ...'), AMP tier is 'I'
    assert "AMP" not in str(out.iloc[0]["acmg_classification"])
    assert out.iloc[0]["amp_tier"] == "I"


def test_stage5_skipped_in_somatic_mode(tmp_path, repo_root, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Seed a classified parquet (stage 5 reads from there)
    df = pd.DataFrame([
        {"variant_id": "x", "sample_id": "S1", "gene": "BRAF",
         "acmg_classification": "AMP Tier I"},
    ])
    out = tmp_path / "data" / "classified" / "acmg_results.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "pipeline.stage5", repo_root / "pipeline" / "05_phenotype_integration.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    config = {
        "input": {"analysis_mode": "somatic"},
        "output": {"output_dir": "data"},
        "phenotype": {"enabled": True},
    }
    result = mod.run(config)
    assert result.get("skipped_reason") == "analysis_mode=somatic"
    # The phenotype JSON sentinel is written so dashboard knows stage was intentional
    assert (tmp_path / "data" / "phenotype" / "patient_phenotype.json").exists()
