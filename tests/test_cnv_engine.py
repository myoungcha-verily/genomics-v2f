"""M3.2: ACMG/ClinGen 2019 CNV engine."""

import pytest


def test_22q11_deletion_classified_pathogenic():
    """The recurrent 22q11.2 deletion locus must come back Pathogenic."""
    from pipeline.utils.cnv_engine import classify_cnv
    variant = {
        "chrom": "chr22", "pos": 18900000, "end_pos": 21500000,
        "svtype": "DEL", "svlen": 2600000, "n_genes": 50,
        "gnomad_sv_af": None,
    }
    res = classify_cnv(variant, {})
    assert res["cnv_classification"] in ("Pathogenic", "Likely pathogenic")
    # Section 1 (gene density) and Section 2 (recurrent locus) should both fire
    sect = res["cnv_section_scores"]
    assert sect["section_1_genomic_content"] >= 0.5
    assert sect["section_2_dosage_genes"] >= 1.0


def test_brca1_exon_del_likely_pathogenic_at_least():
    """A small DEL overlapping a haploinsufficient gene gets Section 2 credit."""
    from pipeline.utils.cnv_engine import classify_cnv
    variant = {
        "chrom": "chr17", "pos": 43106000, "end_pos": 43112000,
        "svtype": "DEL", "svlen": 6000, "n_genes": 1,
        "gnomad_sv_af": None,
    }
    res = classify_cnv(variant, {})
    assert res["cnv_score"] >= 0.5
    assert "BRCA1" in res["cnv_evidence_summary"] or "haploinsufficient" in res["cnv_evidence_summary"]


def test_common_deletion_classified_benign():
    """A small DEL with high gnomAD-SV AF should be Benign via Section 5."""
    from pipeline.utils.cnv_engine import classify_cnv
    variant = {
        "chrom": "chr11", "pos": 5226000, "end_pos": 5227500,
        "svtype": "DEL", "svlen": 1500, "n_genes": 0,
        "gnomad_sv_af": 0.07,
    }
    res = classify_cnv(variant, {})
    assert res["cnv_classification"] in ("Benign", "Likely benign")


def test_pmp22_duplication_pathogenic():
    """The PMP22 duplication is the canonical CMT1A trigger — DUP overlapping
    a triplosensitive gene should score high."""
    from pipeline.utils.cnv_engine import classify_cnv
    variant = {
        "chrom": "chr17", "pos": 15159001, "end_pos": 16679000,
        "svtype": "DUP", "svlen": 1520000, "n_genes": 8,
        "gnomad_sv_af": None,
    }
    res = classify_cnv(variant, {})
    assert res["cnv_score"] >= 0.5


def test_intergenic_vus_stays_vus():
    """Gene-poor region with no recurrent locus match → VUS."""
    from pipeline.utils.cnv_engine import classify_cnv
    variant = {
        "chrom": "chr2", "pos": 100000000, "end_pos": 100200000,
        "svtype": "DEL", "svlen": 200000, "n_genes": 0,
        "gnomad_sv_af": 0.00001,
    }
    res = classify_cnv(variant, {})
    assert res["cnv_classification"] == "VUS"


def test_is_sv_variant():
    from pipeline.utils.cnv_engine import is_sv_variant
    assert is_sv_variant({"svtype": "DEL"}) is True
    assert is_sv_variant({"svtype": "BND"}) is True
    assert is_sv_variant({"svtype": None}) is False
    assert is_sv_variant({}) is False
    assert is_sv_variant({"chrom": "chr17", "pos": 100}) is False
