"""M3.1: SV adapter parses symbolic VCFs."""

import pytest


def test_sv_adapter_registered():
    from adapters import ADAPTER_REGISTRY
    assert "sv" in ADAPTER_REGISTRY
    assert "cnv" in ADAPTER_REGISTRY  # alias


def test_sv_adapter_classify_alt_symbolic():
    from adapters.sv_adapter import SVAdapter
    a = SVAdapter({"input": {}})
    assert a._classify_alt("<DEL>") == ("DEL", None)
    assert a._classify_alt("<DUP>") == ("DUP", None)
    assert a._classify_alt("<INV>") == ("INV", None)
    assert a._classify_alt("<INS:ME>") == ("INS", None)
    assert a._classify_alt("A") == (None, None)


def test_sv_adapter_classify_alt_bnd():
    from adapters.sv_adapter import SVAdapter
    a = SVAdapter({"input": {}})
    svtype, bnd = a._classify_alt("N[chr2:321682[")
    assert svtype == "BND"
    assert bnd == "chr2:321682"


def test_sv_adapter_loads_demo_vcf(repo_root):
    pytest.importorskip("cyvcf2")
    from adapters import get_adapter
    config = {
        "input": {
            "adapter": "sv",
            "vcf_path": str(repo_root / "demo" / "sv" / "sv_demo.vcf.gz"),
            "reference_genome": "GRCh38",
        }
    }
    adapter = get_adapter("sv", config)
    df = adapter.load_variants()
    assert len(df) == 6
    for col in ("svtype", "svlen", "end_pos"):
        assert col in df.columns
    # All rows have an svtype
    assert df["svtype"].notna().all()
    # Specific demo entries present
    assert (df["svtype"] == "DEL").sum() >= 4
    assert (df["svtype"] == "DUP").sum() >= 1
