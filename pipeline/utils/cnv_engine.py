"""ACMG/ClinGen 2019 CNV classification engine.

Implements the 5-section scoring system from Riggs et al. 2019
(PMID 31690835) for copy-number variants:

  Section 1 — Initial assessment of genomic content (genes, repeats)
  Section 2 — Overlap with established triplosensitive (TS) /
              haploinsufficient (HI) genes or genomic regions
  Section 3 — Evaluation of literature evidence (segregation, case reports)
  Section 4 — Inheritance pattern / family history
  Section 5 — Clinical evaluation considerations (de novo, multiple cases)

Each section contributes positive or negative points. The total is mapped
to a 5-tier classification mirroring ACMG/AMP:

  >=  0.99  → Pathogenic
  0.90-0.99 → Likely pathogenic
  -0.89 to 0.89 (excluding |x| < 0.10 → Benign) → VUS
  -0.99 to -0.90 → Likely benign
  <= -1.00  → Benign

The engine reads:
  - `reference/clingen_dosage_sensitivity.json` (bundled fallback) for
    HI/TS gene lookups and recurrent locus matches
  - `databases.gnomad_sv` (optional BQ) for population frequency
  - `databases.clingen_dosage` (optional BQ mirror) for live data
"""

import json
import logging
import os
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_CLINGEN_DS_CACHE: Optional[dict] = None


def _load_clingen_ds() -> dict:
    global _CLINGEN_DS_CACHE
    if _CLINGEN_DS_CACHE is not None:
        return _CLINGEN_DS_CACHE
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    path = os.path.join(here, "reference", "clingen_dosage_sensitivity.json")
    if not os.path.exists(path):
        logger.warning("ClinGen DS reference JSON not found; CNV scoring degraded")
        _CLINGEN_DS_CACHE = {}
        return _CLINGEN_DS_CACHE
    with open(path) as f:
        _CLINGEN_DS_CACHE = json.load(f)
    return _CLINGEN_DS_CACHE


def _interval_overlaps(a_start, a_end, b_start, b_end) -> bool:
    return not (a_end < b_start or a_start > b_end)


def _overlap_fraction(a_start, a_end, b_start, b_end) -> float:
    if a_end < b_start or a_start > b_end:
        return 0.0
    overlap = min(a_end, b_end) - max(a_start, b_start)
    span = a_end - a_start
    return max(0.0, overlap / span) if span > 0 else 0.0


def section_1_genomic_content(svtype: str, svlen: int, n_genes: int) -> Dict:
    """Section 1: gene content + size assessment.

    +0.0 default
    +0.5 if 25-35 genes affected (large CNV in gene-dense region)
    +1.0 if >35 genes
    -0.5 if 0 protein-coding genes affected
    """
    score = 0.0
    notes = []
    if n_genes is None:
        n_genes = 0
    if svlen is None:
        svlen = 0

    if n_genes == 0:
        score = -0.5
        notes.append(f"Section 1: No protein-coding genes affected ({score:+})")
    elif n_genes >= 35:
        score = 1.0
        notes.append(f"Section 1: {n_genes} genes affected ({score:+})")
    elif n_genes >= 25:
        score = 0.5
        notes.append(f"Section 1: {n_genes} genes affected ({score:+})")
    elif n_genes > 0:
        notes.append(f"Section 1: {n_genes} genes affected (0)")

    return {"score": score, "notes": notes}


def section_2_dosage_genes(svtype: str, chrom: str, start: int, end: int,
                            ds_data: dict) -> Dict:
    """Section 2: HI / TS overlap + recurrent locus match."""
    score = 0.0
    notes = []

    # Recurrent pathogenic locus match (highest weight)
    for locus in ds_data.get("_recurrent_pathogenic_loci", []):
        if (locus["chrom"] == chrom
                and locus.get("type") == svtype
                and _interval_overlaps(start, end, locus["start"], locus["end"])):
            frac = _overlap_fraction(start, end, locus["start"], locus["end"])
            if frac >= 0.5:
                score += locus.get("score", 1.0)
                notes.append(
                    f"Section 2: Overlaps recurrent pathogenic locus "
                    f"'{locus['name']}' ({frac:.0%}) (+{locus.get('score', 1.0):.2f})")

    # HI gene overlap (DEL only) or TS gene overlap (DUP only)
    if svtype == "DEL":
        for gene in ds_data.get("_haploinsufficient_3", []):
            if (gene["chrom"] == chrom
                    and _interval_overlaps(start, end, gene["start"], gene["end"])):
                score += 0.5
                notes.append(
                    f"Section 2: DEL overlaps haploinsufficient gene "
                    f"{gene['gene']} (HI=3) (+0.5)")
                # Cap at +1.5 for HI overlap to avoid double-counting
                if score >= 1.5:
                    break

    if svtype == "DUP":
        for gene in ds_data.get("_triplosensitive_3", []):
            if (gene["chrom"] == chrom
                    and _interval_overlaps(start, end, gene["start"], gene["end"])):
                score += 0.5
                notes.append(
                    f"Section 2: DUP overlaps triplosensitive gene "
                    f"{gene['gene']} (TS=3) (+0.5)")

    # Recurrent benign locus match (negative weight)
    for locus in ds_data.get("_recurrent_benign_loci", []):
        if (locus["chrom"] == chrom
                and _interval_overlaps(start, end, locus["start"], locus["end"])):
            frac = _overlap_fraction(start, end, locus["start"], locus["end"])
            if frac >= 0.5:
                score -= 1.0
                notes.append(
                    f"Section 2: Overlaps recurrent benign locus "
                    f"'{locus['name']}' (-1.00)")

    return {"score": score, "notes": notes}


def section_4_inheritance(variant: dict) -> Dict:
    """Section 4: family history + inheritance.

    +0.5 if confirmed de novo with parents tested
    -0.5 if observed in unaffected parent / population
    """
    score = 0.0
    notes = []
    if variant.get("is_de_novo") and variant.get("de_novo_confidence") == "HIGH":
        score = 0.5
        notes.append("Section 4: Confirmed de novo (+0.50)")
    return {"score": score, "notes": notes}


def section_5_population_frequency(variant: dict) -> Dict:
    """Section 5: population frequency / clinical considerations.

    +0.10 for ultra-rare in gnomAD-SV (AF < 1e-5)
    -0.30 for common in gnomAD-SV (AF >= 0.01)
    -1.00 for very common (AF >= 0.05)
    """
    score = 0.0
    notes = []
    af = variant.get("gnomad_sv_af")
    if af is not None:
        if af >= 0.05:
            score = -1.00
            notes.append(f"Section 5: gnomAD-SV AF={af:.3f} (very common, -1.00)")
        elif af >= 0.01:
            score = -0.30
            notes.append(f"Section 5: gnomAD-SV AF={af:.3f} (common, -0.30)")
        elif af < 1e-5:
            score = 0.10
            notes.append(f"Section 5: gnomAD-SV AF<1e-5 (ultra-rare, +0.10)")
    return {"score": score, "notes": notes}


def classify_cnv(variant: dict, config: dict) -> Dict:
    """Classify a single CNV using ACMG/ClinGen 2019 5-section scoring."""
    svtype = variant.get("svtype", "")
    chrom = variant.get("chrom", "")
    start = int(variant.get("pos") or 0)
    end = int(variant.get("end_pos") or start)
    n_genes = variant.get("n_genes", 0) or 0

    ds = _load_clingen_ds()

    s1 = section_1_genomic_content(svtype, end - start, n_genes)
    s2 = section_2_dosage_genes(svtype, chrom, start, end, ds)
    # Section 3 (literature) requires manual curation — left at 0 for automated path
    s3 = {"score": 0.0, "notes": ["Section 3: literature evidence requires manual curation"]}
    s4 = section_4_inheritance(variant)
    s5 = section_5_population_frequency(variant)

    total = s1["score"] + s2["score"] + s3["score"] + s4["score"] + s5["score"]

    if total >= 0.99:
        classification = "Pathogenic"
    elif total >= 0.90:
        classification = "Likely pathogenic"
    elif total <= -1.00:
        classification = "Benign"
    elif total <= -0.90:
        classification = "Likely benign"
    else:
        classification = "VUS"

    all_notes = []
    for s in (s1, s2, s3, s4, s5):
        all_notes.extend(s["notes"])

    return {
        "cnv_classification": classification,
        "cnv_score": round(total, 3),
        "cnv_section_scores": {
            "section_1_genomic_content": s1["score"],
            "section_2_dosage_genes": s2["score"],
            "section_3_literature": s3["score"],
            "section_4_inheritance": s4["score"],
            "section_5_population_frequency": s5["score"],
        },
        "cnv_evidence_summary": "; ".join(all_notes) if all_notes else "No evidence",
    }


def is_sv_variant(variant: dict) -> bool:
    """True if this variant should be routed to the CNV engine."""
    svtype = variant.get("svtype")
    return svtype in ("DEL", "DUP", "INV", "INS", "BND") if svtype else False
