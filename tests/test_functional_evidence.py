"""M2.2: PS3/BS3 wired through MaveDB lookup."""

import pytest


def test_mavedb_lookup_brca1_pathogenic():
    from pipeline.utils.functional_scores import lookup_mavedb
    cfg = {"acmg": {"enable_functional": True}}
    res = lookup_mavedb("BRCA1", "p.Cys61Gly", cfg)
    assert res["found"] is True
    assert res["tier"] == "PS3"
    assert res["function"] == "damaging"
    assert "Findlay" in (res["source"] or "")


def test_mavedb_lookup_brca1_tolerated():
    from pipeline.utils.functional_scores import lookup_mavedb
    cfg = {"acmg": {"enable_functional": True}}
    res = lookup_mavedb("BRCA1", "p.Asp67Tyr", cfg)
    assert res["found"] is True
    assert res["tier"] == "BS3"
    assert res["function"] == "tolerated"


def test_mavedb_lookup_unknown_gene_returns_not_found():
    from pipeline.utils.functional_scores import lookup_mavedb
    cfg = {"acmg": {"enable_functional": True}}
    res = lookup_mavedb("OBSCURE_GENE", "p.X1Y", cfg)
    assert res["found"] is False
    assert res["tier"] is None


def test_eval_ps3_fires_when_enabled():
    from pipeline.utils.acmg_rules import eval_ps3
    variant = {"gene": "TP53", "hgvs_p": "p.Arg175His"}
    cfg = {"acmg": {"enable_functional": True}}
    triggered, strength, evidence = eval_ps3(variant, cfg)
    assert triggered is True
    assert strength == "strong"
    assert "PS3" in evidence


def test_eval_ps3_no_op_when_disabled():
    """enable_functional=false: PS3 must not fire even on a known damaging variant."""
    from pipeline.utils.acmg_rules import eval_ps3
    variant = {"gene": "TP53", "hgvs_p": "p.Arg175His"}
    cfg = {"acmg": {"enable_functional": False}}
    triggered, strength, evidence = eval_ps3(variant, cfg)
    assert triggered is False


def test_eval_bs3_fires_for_tolerated():
    from pipeline.utils.acmg_rules import eval_bs3
    variant = {"gene": "TP53", "hgvs_p": "p.Pro72Arg"}
    cfg = {"acmg": {"enable_functional": True}}
    triggered, strength, evidence = eval_bs3(variant, cfg)
    assert triggered is True
    assert strength == "strong"


def test_ps3_in_all_criteria_registry():
    from pipeline.utils.acmg_rules import ALL_CRITERIA
    assert "PS3" in ALL_CRITERIA
    assert "BS3" in ALL_CRITERIA
