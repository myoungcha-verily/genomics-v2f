"""M2.1: Pejaver/ClinGen 2022 calibrated PP3/BP4 strength tiers."""

import pytest


def _interpret(scores, use_calibrated=True):
    from pipeline.utils.in_silico_scores import interpret_scores
    cfg = {"annotation": {"use_calibrated_tiers": use_calibrated}}
    return interpret_scores(scores, cfg)


def test_revel_strong_pathogenic():
    """REVEL >= 0.932 should land in PP3_strong tier."""
    res = _interpret({"revel": 0.95})
    assert res["pp3"] is True
    assert res["pp3_strength"] == "strong"
    assert res["best_tier"] == "PP3_strong"


def test_revel_moderate_pathogenic():
    """REVEL between 0.773 and 0.932 should be PP3_moderate."""
    res = _interpret({"revel": 0.80})
    assert res["pp3"] is True
    assert res["pp3_strength"] == "moderate"


def test_revel_supporting_pathogenic():
    """REVEL between 0.644 and 0.773 should be PP3_supporting."""
    res = _interpret({"revel": 0.70})
    assert res["pp3"] is True
    assert res["pp3_strength"] == "supporting"


def test_revel_indeterminate():
    """REVEL between 0.290 and 0.644 should be indeterminate."""
    res = _interpret({"revel": 0.45})
    assert res["pp3"] is False
    assert res["bp4"] is False
    assert res["best_tier"] == "indeterminate"


def test_revel_strong_benign():
    """REVEL <= 0.016 should be BP4_strong (Pejaver 2022)."""
    res = _interpret({"revel": 0.005})
    assert res["bp4"] is True
    assert res["best_tier"] == "BP4_strong"


def test_best_tier_picks_strongest_across_predictors():
    """When predictors disagree, the strongest tier wins."""
    # CADD says supporting-pathogenic, REVEL says strong-pathogenic
    res = _interpret({"cadd_phred": 26.0, "revel": 0.95})
    assert res["pp3"] is True
    assert res["pp3_strength"] == "strong"  # REVEL wins


def test_legacy_mode_always_supporting():
    """When use_calibrated_tiers=false, PP3/BP4 always fire as 'supporting'."""
    res = _interpret({"revel": 0.95}, use_calibrated=False)
    assert res["pp3"] is True
    assert res["pp3_strength"] == "supporting"


def test_eval_pp3_returns_calibrated_strength():
    """eval_pp3 propagates the calibrated strength to the engine."""
    from pipeline.utils.acmg_rules import eval_pp3
    cfg = {"annotation": {"use_calibrated_tiers": True}}
    triggered, strength, evidence = eval_pp3({"revel": 0.95}, cfg)
    assert triggered is True
    assert strength == "strong"


def test_eval_bp4_clamps_moderate_to_supporting():
    """ACMG/AMP combination grid has no 'moderate' benign tier;
    eval_bp4 clamps Pejaver BP4_moderate down to 'supporting'."""
    from pipeline.utils.acmg_rules import eval_bp4
    cfg = {"annotation": {"use_calibrated_tiers": True}}
    triggered, strength, evidence = eval_bp4({"revel": 0.10}, cfg)  # in BP4_moderate range
    assert triggered is True
    assert strength == "supporting"
