"""P5.2: Pipeline run manifest provenance."""

import json
from pathlib import Path


def test_build_manifest_shape():
    from pipeline.utils.run_manifest import build_manifest
    cfg = {"input": {"adapter": "single_sample", "analysis_mode": "germline",
                      "reference_genome": "GRCh38"}}
    m = build_manifest(cfg)
    assert m["schema_version"] == "1.0"
    assert m["analysis_mode"] == "germline"
    assert m["adapter"] == "single_sample"
    assert m["reference_genome"] == "GRCh38"
    assert m["config_hash"].startswith("sha256:")
    assert m["run_id"]
    assert "started_at" in m
    assert "reference_data" in m


def test_config_hash_stable():
    """Same config → same hash."""
    from pipeline.utils.run_manifest import build_manifest
    cfg = {"input": {"adapter": "single_sample"}, "acmg": {"ba1_threshold": 0.05}}
    m1 = build_manifest(cfg)
    m2 = build_manifest(cfg)
    assert m1["config_hash"] == m2["config_hash"]


def test_config_hash_changes():
    """Different config → different hash."""
    from pipeline.utils.run_manifest import build_manifest
    m1 = build_manifest({"input": {"adapter": "single_sample"}})
    m2 = build_manifest({"input": {"adapter": "trio"}})
    assert m1["config_hash"] != m2["config_hash"]


def test_write_manifest_to_disk(tmp_path):
    from pipeline.utils.run_manifest import write_manifest
    cfg = {"input": {"adapter": "single_sample"}}
    m = write_manifest(cfg, str(tmp_path))
    out = tmp_path / "run_manifest.json"
    assert out.exists()
    loaded = json.loads(out.read_text())
    assert loaded["schema_version"] == m["schema_version"]
    assert loaded["config_hash"] == m["config_hash"]


def test_manifest_footer_html():
    from pipeline.utils.run_manifest import build_manifest, manifest_footer_html
    cfg = {"input": {"analysis_mode": "germline"}}
    m = build_manifest(cfg)
    html = manifest_footer_html(m)
    assert "run-manifest" in html
    assert "germline" in html
    assert "config=" in html
