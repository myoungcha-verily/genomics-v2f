"""M2.3: Bayesian ACMG (Tavtigian 2018) posterior probability."""

import pytest


def test_no_evidence_returns_prior():
    from pipeline.utils.bayesian_acmg import compute_posterior
    res = compute_posterior([], [], prior=0.1)
    # No evidence → posterior ≈ prior
    assert abs(res["posterior_prob"] - 0.1) < 0.001
    assert res["log2_odds_ratio"] == 0.0


def test_pvs1_alone_gives_high_posterior():
    """PVS1 = 8 log2 OR. With prior 0.1 → very high posterior."""
    from pipeline.utils.bayesian_acmg import compute_posterior
    res = compute_posterior([("PVS1", "very_strong")], [], prior=0.1)
    # Prior odds = 0.1/0.9 ≈ 0.111
    # Posterior odds = 0.111 * 256 = 28.4
    # Posterior = 28.4 / 29.4 ≈ 0.966
    assert res["posterior_prob"] > 0.95


def test_two_strong_pathogenic_classified_pathogenic():
    """ACMG rule: 2 strong → pathogenic. Tavtigian Bayes should mirror this."""
    from pipeline.utils.bayesian_acmg import compute_posterior, classify_from_posterior
    res = compute_posterior([("PS1", "strong"), ("PS3", "strong")], [],
                              prior=0.1)
    # Posterior odds = 0.111 * 16 * 16 = 28.4 → 0.966
    assert res["posterior_prob"] > 0.95
    assert classify_from_posterior(res["posterior_prob"]) in (
        "Likely pathogenic", "Pathogenic")


def test_ba1_override_zeroes_posterior():
    """Standalone benign → posterior = 0 regardless of pathogenic evidence."""
    from pipeline.utils.bayesian_acmg import compute_posterior
    res = compute_posterior(
        [("PVS1", "very_strong")],
        [("BA1", "standalone")],
        prior=0.1,
    )
    assert res["posterior_prob"] == 0.0
    assert res["ba1_override"] is True


def test_classify_from_posterior_buckets():
    from pipeline.utils.bayesian_acmg import classify_from_posterior
    assert classify_from_posterior(0.995) == "Pathogenic"
    assert classify_from_posterior(0.95) == "Likely pathogenic"
    assert classify_from_posterior(0.50) == "VUS"
    assert classify_from_posterior(0.05) == "Likely benign"
    assert classify_from_posterior(0.001) == "Benign"


def test_two_supporting_benign_gives_likely_benign_range():
    from pipeline.utils.bayesian_acmg import compute_posterior, classify_from_posterior
    res = compute_posterior([], [("BP4", "supporting"), ("BP6", "supporting")],
                              prior=0.1)
    # Prior odds 0.111 * 0.5 * 0.5 = 0.0277 → posterior 0.027
    assert classify_from_posterior(res["posterior_prob"]) == "Likely benign"


def test_per_gene_prior_lookup():
    from pipeline.utils.bayesian_acmg import get_prior
    cfg = {"acmg": {"bayesian": {
        "prior_pathogenicity": 0.1,
        "gene_priors": {"BRCA1": 0.30},  # higher prior for BRCA1
    }}}
    assert get_prior("BRCA1", cfg) == 0.30
    assert get_prior("UNKNOWN", cfg) == 0.10
