"""AMP/ASCO/CAP 2017 four-tier somatic variant interpretation engine.

Reference: Li et al. 2017 (PMID: 27993330). Standards and Guidelines for
the Interpretation and Reporting of Sequence Variants in Cancer.

Tiers:
  I    — Variants of strong clinical significance (FDA-approved therapy,
         professional guidelines, well-powered studies)
  II   — Variants of potential clinical significance (off-label therapy,
         multiple small studies, preclinical evidence with consensus)
  III  — Variants of unknown clinical significance (VUS for the somatic
         framework — not yet annotated in databases)
  IV   — Benign / likely benign in the somatic context (e.g. germline
         polymorphisms, recurrent neutral SNPs)

Inputs to the classifier:
  - Variant gene + protein change
  - VAF (tumor allele fraction)
  - Whether matched normal is available
  - Knowledge-base hits from CIViC / OncoKB / COSMIC

Output:
  {
    "amp_tier": "I" | "II" | "III" | "IV",
    "evidence_summary": str,
    "drug_targets": List[str],
    "knowledge_bases_consulted": List[str],
    "highest_evidence_level": str,
    "vaf": float,
    "is_paired": bool,
  }

The engine merges evidence from whichever knowledge bases are enabled in
config. OncoKB takes precedence when present (commercial-grade curation);
CIViC is the open-source baseline.
"""

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


# Rough AMP-tier ordering for "best wins"
TIER_RANK = {"I": 0, "II": 1, "III": 2, "IV": 3, None: 4}


def _best_tier(*tiers: Optional[str]) -> Optional[str]:
    """Return the highest-confidence (lowest-rank) non-None tier."""
    valid = [t for t in tiers if t]
    if not valid:
        return None
    return min(valid, key=lambda t: TIER_RANK.get(t, 4))


def classify_somatic_variant(variant: dict, config: dict) -> Dict:
    """Apply AMP/ASCO/CAP 2017 tiering to a single variant.

    `variant` should carry at minimum: gene, hgvs_p (or protein_change),
    chrom, pos, ref, alt, tumor_vaf, is_paired.
    """
    db_config = config.get("databases", {}) or {}
    amp_config = config.get("acmg_amp", {}) or {}
    priority = amp_config.get("knowledge_base_priority",
                                ["oncokb", "civic", "cosmic"])

    gene = variant.get("gene") or ""
    protein_change = variant.get("hgvs_p") or variant.get("protein_change") or ""
    chrom = variant.get("chrom", "")
    pos = variant.get("pos", 0)
    ref = variant.get("ref", "")
    alt = variant.get("alt", "")
    vaf = float(variant.get("tumor_vaf") or variant.get("allele_fraction") or 0.0)
    is_paired = bool(variant.get("is_paired", False))

    consulted: List[str] = []
    drugs = set()
    levels = []
    tiers = []
    notes = []
    sources = []

    # Iterate KBs in priority order (caller-configurable)
    for kb in priority:
        if kb == "oncokb":
            from pipeline.utils import oncokb_client
            if oncokb_client.is_enabled(db_config):
                token = oncokb_client.get_token(db_config)
                api_url = (db_config.get("oncokb", {}) or {}).get("api_url",
                                                                    "https://www.oncokb.org/api/v1")
                res = oncokb_client.query_oncokb(gene, protein_change,
                                                  api_token=token, api_url=api_url)
                consulted.append("oncokb")
                if res.get("found"):
                    if res.get("amp_tier"):
                        tiers.append(res["amp_tier"])
                    if res.get("highest_level"):
                        levels.append(f"OncoKB:{res['highest_level']}")
                    drugs.update(res.get("drug_targets") or [])
                    if res.get("source_url"):
                        sources.append(res["source_url"])
                    if res.get("oncogenic"):
                        notes.append(f"OncoKB: {res['oncogenic']}")
                elif res.get("error") and res["error"] != "no_token":
                    notes.append(f"OncoKB error: {res['error'][:120]}")

        elif kb == "civic":
            from pipeline.utils import civic_client
            if civic_client.is_enabled(db_config):
                api_url = (db_config.get("civic", {}) or {}).get("api_url",
                                                                   "https://civicdb.org/api/graphql")
                res = civic_client.query_civic(gene, protein_change, api_url=api_url)
                consulted.append("civic")
                if res.get("found"):
                    if res.get("amp_tier"):
                        tiers.append(res["amp_tier"])
                    if res.get("highest_level"):
                        levels.append(f"CIViC:{res['highest_level']}")
                    drugs.update(res.get("drug_targets") or [])
                    if res.get("source_url"):
                        sources.append(res["source_url"])
                elif res.get("error"):
                    notes.append(f"CIViC error: {res['error'][:120]}")

        elif kb == "cosmic":
            from pipeline.utils import cosmic_client
            if cosmic_client.is_enabled(db_config):
                bq_table = (db_config.get("cosmic", {}) or {}).get("bq_table", "")
                res = cosmic_client.query_cosmic(chrom, pos, ref, alt,
                                                   bq_table=bq_table)
                consulted.append("cosmic")
                if res.get("found"):
                    notes.append(
                        f"COSMIC: {res['occurrences']} occurrences across "
                        f"{len(res['tissue_distribution'])} tissues"
                    )
                    if res.get("source_url"):
                        sources.append(res["source_url"])
                    # COSMIC alone doesn't determine tier, but recurrence
                    # in many tissues nudges toward Tier II
                    if res["occurrences"] >= 10 and not tiers:
                        tiers.append("II")

    # Final tier selection
    best = _best_tier(*tiers)

    # Tumor-only VAF guardrail: very high VAF (>0.45) without paired normal
    # is suspicious for germline contamination — note it in the summary
    if not is_paired and vaf > 0.45 and best in ("I", "II"):
        notes.append(
            f"High VAF ({vaf:.2f}) without matched normal — consider germline"
        )

    if best is None:
        # No KB hit. Default to Tier III (VUS for somatic) if anything was
        # consulted; if nothing was consulted at all we still report III
        # rather than failing.
        best = "III"
        if not consulted:
            notes.append("No somatic knowledge bases enabled")

    summary = "; ".join(notes) if notes else f"Tier {best} from {','.join(consulted) or 'no KBs'}"

    return {
        "amp_tier": best,
        "evidence_summary": summary,
        "drug_targets": sorted(drugs),
        "knowledge_bases_consulted": consulted,
        "highest_evidence_level": ", ".join(levels) if levels else None,
        "evidence_sources": sources,
        "vaf": vaf,
        "is_paired": is_paired,
    }


def tier_label(tier: str) -> str:
    """Human-readable label for an AMP tier."""
    return {
        "I": "Tier I — Strong clinical significance",
        "II": "Tier II — Potential clinical significance",
        "III": "Tier III — Unknown clinical significance",
        "IV": "Tier IV — Benign / Likely benign (somatic)",
    }.get(tier, tier)
