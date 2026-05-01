"""P6.1: Per-stage acceptance gates."""

import pytest


def _cfg(mode="warning", **stage_overrides):
    cfg = {"pipeline": {"gates": {"mode": mode, **stage_overrides}}}
    return cfg


def test_stage_1_titv_in_range_passes():
    from pipeline.utils.quality_gates import check_stage_1
    qc = {"n_variants": 10000, "ti_tv_ratio": 2.1}
    results = check_stage_1(qc, _cfg())
    titv = next(r for r in results if r.name == "stage_1.ti_tv_ratio")
    assert titv.passed is True


def test_stage_1_titv_out_of_range_warns():
    from pipeline.utils.quality_gates import check_stage_1
    qc = {"n_variants": 10000, "ti_tv_ratio": 0.5}
    results = check_stage_1(qc, _cfg())
    titv = next(r for r in results if r.name == "stage_1.ti_tv_ratio")
    assert titv.passed is False
    assert titv.severity == "warning"


def test_stage_1_hard_fail_raises():
    from pipeline.utils.quality_gates import evaluate, GateFailure
    qc = {"n_variants": 1, "ti_tv_ratio": 2.0}
    cfg = _cfg(mode="hard_fail", stage_1={"min_variants": 100})
    with pytest.raises(GateFailure) as exc:
        evaluate("stage_1", qc, cfg)
    assert exc.value.gate.name == "stage_1.min_variants"


def test_stage_3_clinvar_match_rate_below_threshold():
    """The exact failure mode that hit M0: queries ran but returned 0 matches."""
    from pipeline.utils.quality_gates import check_stage_3
    summary = {"total_variants": 10, "clinvar_matches": 0, "gnomad_matches": 0}
    cfg = _cfg(stage_3={"clinvar_match_rate_min": 0.30})
    results = check_stage_3(summary, cfg)
    cv = next(r for r in results if r.name == "stage_3.clinvar_match_rate")
    assert cv.passed is False
    assert "0.0%" in cv.actual
    assert "BQ table misconfigured" in cv.message


def test_stage_3_no_threshold_no_gate():
    """If no threshold is configured, no gate is generated."""
    from pipeline.utils.quality_gates import check_stage_3
    summary = {"total_variants": 10, "clinvar_matches": 0, "gnomad_matches": 0}
    cfg = _cfg()  # no threshold
    results = check_stage_3(summary, cfg)
    assert all(r.name != "stage_3.clinvar_match_rate" for r in results)


def test_stage_4_all_vus_warning():
    """Catches the all-VUS regression — M0 produced 10/10 VUS due to broken
    upstream enrichment."""
    from pipeline.utils.quality_gates import check_stage_4
    counts = {"VUS": 10, "Pathogenic": 0, "Benign": 0}
    results = check_stage_4(counts, _cfg())
    vus = next(r for r in results if r.name == "stage_4.vus_fraction")
    assert vus.passed is False  # 100% VUS > default 95% max


def test_stage_4_normal_distribution_passes():
    from pipeline.utils.quality_gates import check_stage_4
    counts = {"VUS": 20, "Pathogenic": 5, "Benign": 25,
              "Likely benign": 30, "Likely pathogenic": 5}
    results = check_stage_4(counts, _cfg())
    vus = next(r for r in results if r.name == "stage_4.vus_fraction")
    assert vus.passed is True


def test_evaluate_returns_dict_form():
    from pipeline.utils.quality_gates import evaluate
    qc = {"n_variants": 10000, "ti_tv_ratio": 2.1}
    out = evaluate("stage_1", qc, _cfg())
    assert isinstance(out, list)
    assert all(isinstance(r, dict) for r in out)
    assert all("passed" in r and "name" in r for r in out)


def test_evaluate_unknown_stage_returns_empty():
    from pipeline.utils.quality_gates import evaluate
    assert evaluate("stage_99", {}, _cfg()) == []
