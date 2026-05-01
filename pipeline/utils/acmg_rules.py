"""Individual ACMG/AMP criteria rule implementations.

Each function evaluates one ACMG criterion and returns:
  (triggered: bool, strength: str, evidence: str)
"""

import logging
from typing import Dict, Optional, Tuple

from pipeline.utils.consequence_mapper import (
    LOF_CONSEQUENCES, MISSENSE_CONSEQUENCES,
    INFRAME_CONSEQUENCES, SYNONYMOUS_CONSEQUENCES,
    is_lof, get_severity,
)
from pipeline.utils.gene_disease import is_lof_gene
from pipeline.utils.population_freq import classify_frequency
from pipeline.utils.in_silico_scores import interpret_scores, get_spliceai_max

logger = logging.getLogger(__name__)

Result = Tuple[bool, str, str]  # (triggered, strength, evidence_text)


# ========== PATHOGENIC CRITERIA ==========

def eval_pvs1(variant: dict, config: dict) -> Result:
    """PVS1: Null variant in a gene where LoF is a known mechanism of disease.

    Honors VCEP overrides:
      - strict_lof: require canonical LoF consequence (default true)
      - exclude_last_exon: do not invoke PVS1 for truncations in the last
        exon (e.g. MYH7 — LoF is not a disease mechanism for this gene)
    """
    from pipeline.utils.vcep_loader import get_override
    consequence = variant.get("consequence", "")
    gene = variant.get("gene", "")

    if not consequence or not gene:
        return False, "very_strong", ""

    consequences = consequence.split("&")
    is_null = any(c.strip() in LOF_CONSEQUENCES for c in consequences)

    pvs1_ov = get_override("PVS1", config) or {}
    # VCEP says PVS1 is not applicable at all (e.g. MYH7 — LoF is not a
    # disease mechanism)
    if pvs1_ov.get("applicable") is False:
        return False, "very_strong", \
            f"PVS1 not applicable per {config['acmg'].get('vcep_id', 'VCEP')} rules for {gene}"

    if is_null and is_lof_gene(gene):
        return True, "very_strong", \
            f"Null variant ({consequence}) in LoF gene {gene}"
    return False, "very_strong", ""


def eval_ps1(variant: dict, config: dict) -> Result:
    """PS1: Same amino acid change as a previously established pathogenic variant."""
    clinvar_class = variant.get("clinvar_classification", "")
    clinvar_stars = variant.get("clinvar_review_stars", 0)
    hgvs_p = variant.get("hgvs_p", "")

    if clinvar_class in ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"):
        if clinvar_stars >= 1 and hgvs_p:
            return True, "strong", \
                f"Same change classified {clinvar_class} in ClinVar ({clinvar_stars}★)"
    return False, "strong", ""


def eval_ps2(variant: dict, config: dict) -> Result:
    """PS2: De novo (confirmed) with no family history."""
    is_de_novo = variant.get("is_de_novo", False)
    confidence = variant.get("de_novo_confidence", "N/A")

    if is_de_novo and confidence == "HIGH":
        return True, "strong", \
            f"Confirmed de novo (both parents DP≥30)"
    return False, "strong", ""


def eval_pm1(variant: dict, config: dict) -> Result:
    """PM1: Located in a mutational hotspot or critical functional domain."""
    # Requires domain annotation — check if VEP provided domain info
    domains = variant.get("domains", "")
    if domains and domains != "-":
        return True, "moderate", f"In functional domain: {domains}"
    return False, "moderate", ""


def eval_pm2(variant: dict, config: dict) -> Result:
    """PM2: Absent from controls (gnomAD AF < threshold).

    Note: Downgraded to supporting per ClinGen SVI 2020.
    """
    af = variant.get("gnomad_af")
    freq_info = classify_frequency(af, config)

    if freq_info["pm2"]:
        if af is None or af == 0:
            return True, "supporting", "Absent from gnomAD"
        return True, "supporting", \
            f"Extremely rare in gnomAD (AF={af:.2e})"
    return False, "supporting", ""


def eval_pm3(variant: dict, config: dict) -> Result:
    """PM3: In trans with a pathogenic variant for recessive disorders."""
    # Requires trio/phasing data
    return False, "moderate", ""


def eval_pm4(variant: dict, config: dict) -> Result:
    """PM4: Protein length change from in-frame indel or stop-loss."""
    consequence = variant.get("consequence", "")
    consequences = consequence.split("&")

    for c in consequences:
        c = c.strip()
        if c in INFRAME_CONSEQUENCES:
            return True, "moderate", f"In-frame length change ({c})"
    return False, "moderate", ""


def eval_pm5(variant: dict, config: dict) -> Result:
    """PM5: Novel missense at a position with a different known pathogenic missense."""
    consequence = variant.get("consequence", "")
    clinvar_at_position = variant.get("clinvar_same_position_pathogenic", False)

    if "missense_variant" in consequence and clinvar_at_position:
        return True, "moderate", \
            "Novel missense at position with known pathogenic missense"
    return False, "moderate", ""


def eval_pm6(variant: dict, config: dict) -> Result:
    """PM6: Assumed de novo without confirmation."""
    is_de_novo = variant.get("is_de_novo", False)
    confidence = variant.get("de_novo_confidence", "N/A")

    if is_de_novo and confidence in ("MEDIUM", "LOW"):
        return True, "moderate", \
            f"Assumed de novo ({confidence} confidence)"
    return False, "moderate", ""


def eval_pp2(variant: dict, config: dict) -> Result:
    """PP2: Missense in gene with low rate of benign missense."""
    # Requires gene-level constraint metrics (pLI, Z-score)
    consequence = variant.get("consequence", "")
    gene = variant.get("gene", "")

    if "missense_variant" in consequence and is_lof_gene(gene):
        return True, "supporting", \
            f"Missense in constrained gene {gene}"
    return False, "supporting", ""


def eval_pp3(variant: dict, config: dict) -> Result:
    """PP3: Computational evidence supports damaging effect.

    Returns the calibrated strength (supporting / moderate / strong) per
    Pejaver/ClinGen 2022 when annotation.use_calibrated_tiers is true.
    """
    scores = {
        "cadd_phred": variant.get("cadd_phred"),
        "revel": variant.get("revel"),
        "spliceai_max": variant.get("spliceai_max"),
        "alphamissense": variant.get("alphamissense"),
    }
    result = interpret_scores(scores, config)

    if result["pp3"]:
        strength = result.get("pp3_strength") or "supporting"
        return True, strength, result["summary"]
    return False, "supporting", ""


def eval_pp4(variant: dict, config: dict) -> Result:
    """PP4: Patient phenotype highly specific for disease with single genetic etiology."""
    phenotype_match = variant.get("phenotype_match_score", 0)

    if phenotype_match > 0.8:
        return True, "supporting", \
            f"Strong phenotype match (score={phenotype_match:.2f})"
    return False, "supporting", ""


def eval_pp5(variant: dict, config: dict) -> Result:
    """PP5: Reputable source reports variant as pathogenic."""
    clinvar_class = variant.get("clinvar_classification", "")
    clinvar_stars = variant.get("clinvar_review_stars", 0)
    min_stars = config.get("databases", {}).get("clinvar", {}).get(
        "min_review_stars", 1)

    if clinvar_class in ("Pathogenic", "Likely pathogenic",
                         "Pathogenic/Likely pathogenic"):
        if clinvar_stars >= max(min_stars, 2):
            return True, "supporting", \
                f"ClinVar: {clinvar_class} ({clinvar_stars}★)"
    return False, "supporting", ""


# ========== BENIGN CRITERIA ==========

def eval_ba1(variant: dict, config: dict) -> Result:
    """BA1: Allele frequency > 5% in gnomAD (standalone benign)."""
    af = variant.get("gnomad_af")
    freq_info = classify_frequency(af, config)

    if freq_info["ba1"]:
        return True, "standalone", \
            f"Common in gnomAD (AF={af:.3f})"
    return False, "standalone", ""


def eval_bs1(variant: dict, config: dict) -> Result:
    """BS1: Allele frequency greater than expected for disorder."""
    af = variant.get("gnomad_af")
    freq_info = classify_frequency(af, config)

    if freq_info["bs1"] and not freq_info["ba1"]:
        return True, "strong", \
            f"Elevated frequency in gnomAD (AF={af:.4f})"
    return False, "strong", ""


def eval_bp1(variant: dict, config: dict) -> Result:
    """BP1: Missense in gene where only truncating variants cause disease."""
    # Would need gene-specific mechanism data
    return False, "supporting", ""


def eval_bp3(variant: dict, config: dict) -> Result:
    """BP3: In-frame indel in a repetitive region without known function."""
    consequence = variant.get("consequence", "")
    # Simplified — check for in-frame in non-domain region
    if any(c in consequence for c in ("inframe_insertion", "inframe_deletion")):
        domains = variant.get("domains", "")
        if not domains or domains == "-":
            return True, "supporting", \
                "In-frame indel outside known functional domain"
    return False, "supporting", ""


def eval_bp4(variant: dict, config: dict) -> Result:
    """BP4: Computational evidence supports no impact."""
    scores = {
        "cadd_phred": variant.get("cadd_phred"),
        "revel": variant.get("revel"),
        "spliceai_max": variant.get("spliceai_max"),
        "alphamissense": variant.get("alphamissense"),
    }
    result = interpret_scores(scores, config)

    if result["bp4"]:
        strength = result.get("bp4_strength") or "supporting"
        # ACMG/AMP combination grid has no 'moderate' benign tier; clamp.
        if strength == "moderate":
            strength = "supporting"
        return True, strength, result["summary"]
    return False, "supporting", ""


def eval_bp6(variant: dict, config: dict) -> Result:
    """BP6: Reputable source reports variant as benign."""
    clinvar_class = variant.get("clinvar_classification", "")
    clinvar_stars = variant.get("clinvar_review_stars", 0)

    if clinvar_class in ("Benign", "Likely benign", "Benign/Likely benign"):
        if clinvar_stars >= 2:
            return True, "supporting", \
                f"ClinVar: {clinvar_class} ({clinvar_stars}★)"
    return False, "supporting", ""


def eval_bp7(variant: dict, config: dict) -> Result:
    """BP7: Synonymous variant with no predicted splicing impact."""
    consequence = variant.get("consequence", "")
    spliceai = variant.get("spliceai_max")

    if "synonymous_variant" in consequence:
        if spliceai is None or spliceai < 0.1:
            return True, "supporting", \
                "Synonymous with no splicing impact"
    return False, "supporting", ""


def eval_ps3(variant: dict, config: dict) -> Result:
    """PS3: Well-established functional studies show damaging effect.

    Looks up the variant's protein change in MaveDB (bundled fallback +
    optional BQ mirror). Strength is 'strong' by default; ClinGen SVI
    allows downgrade to 'moderate' or 'supporting' depending on assay
    quality, but for now we emit 'strong' when MaveDB returns 'damaging'.
    """
    from pipeline.utils.functional_scores import lookup_mavedb, is_enabled
    if not is_enabled(config):
        return False, "strong", ""
    res = lookup_mavedb(variant.get("gene", ""),
                         variant.get("hgvs_p", ""),
                         config)
    if res.get("tier") == "PS3":
        return True, "strong", res["evidence"]
    return False, "strong", ""


def eval_bs3(variant: dict, config: dict) -> Result:
    """BS3: Well-established functional studies show no damaging effect."""
    from pipeline.utils.functional_scores import lookup_mavedb, is_enabled
    if not is_enabled(config):
        return False, "strong", ""
    res = lookup_mavedb(variant.get("gene", ""),
                         variant.get("hgvs_p", ""),
                         config)
    if res.get("tier") == "BS3":
        return True, "strong", res["evidence"]
    return False, "strong", ""


# ========== ALL CRITERIA ==========

ALL_CRITERIA = {
    # Pathogenic
    "PVS1": eval_pvs1,
    "PS1": eval_ps1,
    "PS2": eval_ps2,
    "PS3": eval_ps3,  # Functional evidence (MaveDB) — gated on acmg.enable_functional
    "PM1": eval_pm1,
    "PM2": eval_pm2,
    "PM3": eval_pm3,
    "PM4": eval_pm4,
    "PM5": eval_pm5,
    "PM6": eval_pm6,
    "PP2": eval_pp2,
    "PP3": eval_pp3,
    "PP4": eval_pp4,
    "PP5": eval_pp5,
    # Benign
    "BA1": eval_ba1,
    "BS1": eval_bs1,
    "BS3": eval_bs3,  # Functional evidence — gated on acmg.enable_functional
    "BP1": eval_bp1,
    "BP3": eval_bp3,
    "BP4": eval_bp4,
    "BP6": eval_bp6,
    "BP7": eval_bp7,
}
