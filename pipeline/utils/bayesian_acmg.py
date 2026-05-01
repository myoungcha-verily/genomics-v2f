"""Bayesian ACMG/AMP per Tavtigian et al. 2018 (PMID 29300386).

Converts the categorical ACMG/AMP combination rules into a quantitative
posterior probability of pathogenicity. Useful for borderline VUSes and
for downstream tools that need a continuous score rather than a bucket.

Formula
-------
Each evidence type contributes a multiplicative odds-ratio:

    Strength            Pathogenic OR    Benign OR
    -----------------   --------------   ----------
    very_strong (PVS)   2^8 = 256        2^-8 = 1/256
    strong (PS)         2^4 = 16         2^-4 = 1/16
    moderate (PM)       2^2 = 4          2^-2 = 1/4
    supporting (PP/BP)  2^1 = 2          2^-1 = 1/2

The standalone benign criterion (BA1) acts as a hard override → posterior 0.

Posterior probability:

    posterior = (prior_odds × product_of_evidence_ORs)
              / (1 + prior_odds × product_of_evidence_ORs)

Default prior is 0.1 (10% prior probability of pathogenicity for an
unselected variant). Per-gene priors can be configured.

Reference: Tavtigian SV et al. Modeling the ACMG/AMP variant classification
guidelines as a Bayesian classification framework. Genet Med. 2018.
"""

import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

# Evidence weights as exponents of 2 (Tavtigian 2018 Table 1)
PATHOGENIC_WEIGHT = {
    "very_strong": 8,   # PVS1
    "strong": 4,        # PS1-PS4
    "moderate": 2,      # PM1-PM6
    "supporting": 1,    # PP1-PP5
}
BENIGN_WEIGHT = {
    "standalone": 1000,  # BA1 — hard override (handled separately)
    "strong": -4,        # BS1-BS4
    "supporting": -1,    # BP1-BP7
    "moderate": -2,      # Pejaver-style (only used if BP4 emits moderate)
}

DEFAULT_PRIOR = 0.1


def _odds_from_prob(p: float) -> float:
    if p >= 1.0:
        return float("inf")
    if p <= 0:
        return 0.0
    return p / (1.0 - p)


def _prob_from_odds(o: float) -> float:
    if o == float("inf"):
        return 1.0
    return o / (1.0 + o)


def compute_posterior(pathogenic_criteria: List[Tuple[str, str]],
                       benign_criteria: List[Tuple[str, str]],
                       prior: float = DEFAULT_PRIOR) -> Dict:
    """Compute posterior probability of pathogenicity given a list of
    triggered ACMG criteria.

    Inputs are lists of (criterion_name, strength) tuples — same shape as
    `acmg_engine.classify_variant()` produces in `pathogenic_criteria` and
    `benign_criteria`.

    Returns:
        {
          "posterior_prob": float (0..1),
          "log2_odds_ratio": float,
          "ba1_override": bool,
          "components": list[{name, strength, log2_or}],
        }
    """
    components = []

    # BA1 override (any standalone benign → automatically benign)
    if any(s == "standalone" for _, s in benign_criteria):
        return {
            "posterior_prob": 0.0,
            "log2_odds_ratio": float("-inf"),
            "ba1_override": True,
            "components": [{"name": n, "strength": "standalone", "log2_or": "-inf"}
                            for n, s in benign_criteria if s == "standalone"],
        }

    log2_or_total = 0.0
    for name, strength in pathogenic_criteria:
        w = PATHOGENIC_WEIGHT.get(strength, 0)
        log2_or_total += w
        components.append({"name": name, "strength": strength, "log2_or": w})

    for name, strength in benign_criteria:
        w = BENIGN_WEIGHT.get(strength, 0)
        log2_or_total += w  # already negative
        components.append({"name": name, "strength": strength, "log2_or": w})

    prior_odds = _odds_from_prob(max(min(prior, 1.0), 0.0))
    posterior_odds = prior_odds * (2.0 ** log2_or_total)
    posterior = _prob_from_odds(posterior_odds)

    return {
        "posterior_prob": posterior,
        "log2_odds_ratio": log2_or_total,
        "ba1_override": False,
        "components": components,
    }


def classify_from_posterior(posterior: float) -> str:
    """Tavtigian 2018 Table 2 classification thresholds (mirror ACMG buckets).

    Pathogenic >=  0.99
    Likely pathogenic  0.90 - 0.99
    VUS                0.10 - 0.90
    Likely benign      0.01 - 0.10
    Benign            <  0.01
    """
    if posterior >= 0.99:
        return "Pathogenic"
    if posterior >= 0.90:
        return "Likely pathogenic"
    if posterior < 0.01:
        return "Benign"
    if posterior < 0.10:
        return "Likely benign"
    return "VUS"


def get_prior(gene: str, config: dict) -> float:
    """Return the prior probability of pathogenicity for a gene.

    Looks up acmg.bayesian.gene_priors[gene] if present, else falls back to
    acmg.bayesian.prior_pathogenicity (default 0.1).
    """
    bayes_cfg = (config.get("acmg", {}) or {}).get("bayesian", {}) or {}
    gene_priors = bayes_cfg.get("gene_priors", {}) or {}
    if gene and gene in gene_priors:
        return float(gene_priors[gene])
    return float(bayes_cfg.get("prior_pathogenicity", DEFAULT_PRIOR))


def is_enabled(config: dict) -> bool:
    bayes_cfg = (config.get("acmg", {}) or {}).get("bayesian", {}) or {}
    return bool(bayes_cfg.get("enabled", True))  # default-on
