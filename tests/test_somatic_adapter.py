"""M1: Somatic adapter parses tumor-normal VCFs with VAF columns."""

import pytest


def test_somatic_adapter_registered():
    from adapters import ADAPTER_REGISTRY
    assert "somatic" in ADAPTER_REGISTRY
    assert "somatic_pair" in ADAPTER_REGISTRY


def test_somatic_adapter_loads_demo_vcf(repo_root):
    """The bundled somatic demo VCF parses through the SomaticAdapter."""
    pytest.importorskip("cyvcf2", reason="cyvcf2 not available in test env")

    from adapters import get_adapter
    config = {
        "input": {
            "adapter": "somatic",
            "analysis_mode": "somatic",
            "vcf_path": str(repo_root / "demo" / "somatic" / "tumor_normal.vcf.gz"),
            "reference_genome": "GRCh38",
            "somatic": {
                "tumor_sample_id": "DEMO_TUMOR",
                "normal_sample_id": "DEMO_NORMAL",
            },
        },
        "acmg_amp": {"vaf_min_somatic": 0.05},
    }
    adapter = get_adapter("somatic", config)
    samples = adapter.get_sample_ids()
    assert "DEMO_TUMOR" in samples
    assert "DEMO_NORMAL" in samples
    assert adapter.get_proband_id() == "DEMO_TUMOR"

    df = adapter.load_variants()
    # Demo VCF has 7 variants
    assert len(df) >= 5
    # Required somatic columns present
    for col in ("tumor_vaf", "tumor_alt_count", "normal_vaf", "is_paired"):
        assert col in df.columns
    # All rows are from tumor sample
    assert (df["sample_id"] == "DEMO_TUMOR").all()
    # is_paired flag set since matched normal exists
    assert df["is_paired"].all()
    # Tumor VAFs above threshold
    assert (df["tumor_vaf"] >= 0.05).all()


def test_somatic_adapter_tumor_only_inference(repo_root, tmp_path):
    """Single-sample VCF (tumor only) should still parse, with is_paired=False."""
    pytest.importorskip("cyvcf2", reason="cyvcf2 not available")
    # We can use the same demo VCF but config without normal_sample_id
    from adapters import get_adapter
    config = {
        "input": {
            "adapter": "somatic",
            "analysis_mode": "somatic",
            "vcf_path": str(repo_root / "demo" / "somatic" / "tumor_normal.vcf.gz"),
            "reference_genome": "GRCh38",
            "somatic": {
                "tumor_sample_id": "DEMO_TUMOR",
                # normal omitted intentionally
            },
        },
        "acmg_amp": {"vaf_min_somatic": 0.05},
    }
    adapter = get_adapter("somatic", config)
    df = adapter.load_variants()
    assert len(df) >= 5
    # is_paired should be False since we didn't specify a normal
    assert (df["is_paired"] == False).all()
