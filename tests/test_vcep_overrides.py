"""M2.4: Gene-specific VCEP rule overrides."""

import pytest


def test_supported_genes_loaded():
    from pipeline.utils.vcep_loader import list_supported_genes
    genes = list_supported_genes({})
    for expected in ("BRCA1", "BRCA2", "PTEN", "TP53", "MYH7", "MYBPC3",
                       "KCNQ1", "SCN5A"):
        assert expected in genes, f"VCEP rule file missing for {expected}"


def test_vcep_loader_returns_none_for_unknown_gene():
    from pipeline.utils.vcep_loader import load_vcep_rules
    assert load_vcep_rules("OBSCURE_GENE", {}) is None


def test_pten_pm2_threshold_overridden():
    """PTEN VCEP tightens PM2 threshold to 5e-6."""
    from pipeline.utils.vcep_loader import merge_vcep_overrides
    base = {"acmg": {"pm2_threshold": 0.0001}}
    merged = merge_vcep_overrides({"gene": "PTEN"}, base)
    assert merged["acmg"]["pm2_threshold"] == 0.000005


def test_myh7_pvs1_disabled_via_vcep():
    """MYH7 VCEP rules say PVS1 is not applicable. eval_pvs1 must respect this."""
    from pipeline.utils.vcep_loader import merge_vcep_overrides
    from pipeline.utils.acmg_rules import eval_pvs1

    # A real LoF variant in MYH7 (frameshift) — would normally trigger PVS1
    variant = {"gene": "MYH7", "consequence": "frameshift_variant"}
    base_config = {"acmg": {}}

    # Without VCEP merge — PVS1 would still depend on is_lof_gene; but
    # critically: WITH VCEP merge, PVS1 should not fire even if it would
    # have otherwise.
    merged = merge_vcep_overrides(variant, base_config)
    triggered, strength, evidence = eval_pvs1(variant, merged)
    assert triggered is False
    assert "not applicable" in evidence.lower() or evidence == ""


def test_brca1_pm2_strength_supporting():
    """ENIGMA VCEP downgrades BRCA1 PM2 to supporting."""
    from pipeline.utils.vcep_loader import merge_vcep_overrides, get_override
    merged = merge_vcep_overrides({"gene": "BRCA1"}, {"acmg": {}})
    pm2_ov = get_override("PM2", merged)
    assert pm2_ov is not None
    assert pm2_ov["strength"] == "supporting"


def test_unknown_gene_passes_config_through_unchanged():
    """No VCEP file for the gene → config returned unchanged."""
    from pipeline.utils.vcep_loader import merge_vcep_overrides
    base = {"acmg": {"pm2_threshold": 0.0001}}
    merged = merge_vcep_overrides({"gene": "OBSCURE_GENE"}, base)
    assert merged["acmg"]["pm2_threshold"] == 0.0001
    assert "vcep_overrides" not in merged["acmg"]


def test_vcep_disabled_via_config():
    from pipeline.utils.vcep_loader import merge_vcep_overrides
    base = {"acmg": {"pm2_threshold": 0.0001, "vcep": {"enabled": False}}}
    merged = merge_vcep_overrides({"gene": "PTEN"}, base)
    # Should NOT have applied PTEN's stricter PM2 threshold
    assert merged["acmg"]["pm2_threshold"] == 0.0001
