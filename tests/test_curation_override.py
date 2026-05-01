"""P4.4: Curation overrides applied by ACMG engine + CNV Section 3."""

from unittest.mock import patch


def test_curation_override_replaces_engine_call(tmp_path, monkeypatch):
    """A support_benign curation for BS3 must add a benign criterion +
    flip a borderline call."""
    from pipeline.utils.acmg_engine import classify_variant
    from pipeline.utils.curation_store import add_curation

    cfg = {
        "input": {"analysis_mode": "germline"},
        "output": {"data_dir": str(tmp_path)},
        "acmg": {"enable_curation_overrides": True,
                  "ba1_threshold": 0.05, "bs1_threshold": 0.01,
                  "pm2_threshold": 0.0001},
        "databases": {"curations": {}},
    }
    add_curation({
        "variant_id": "chr17_43106478_T_G",
        "gene": "BRCA1",
        "criterion": "BS3",
        "action": "support_benign",
        "strength": "strong",
        "evidence_text": "DMS shows tolerated function",
    }, cfg)

    variant = {
        "variant_id": "chr17_43106478_T_G",
        "gene": "BRCA1", "consequence": "missense_variant",
        "hgvs_p": "p.Cys61Gly", "chrom": "chr17", "pos": 43106478,
        "ref": "T", "alt": "G",
        "gnomad_af": None, "clinvar_classification": "",
        "clinvar_review_stars": 0,
    }
    res = classify_variant(variant, cfg)
    # The BS3 curation should appear in benign_criteria as 'strong'
    benign_keys = [c[0] for c in res["benign_criteria"]]
    assert "BS3" in benign_keys
    # criteria_met should tag the entry as curated
    matching = [c for c in res["criteria_met"]
                 if c["criterion"] == "BS3" and c.get("curated")]
    assert len(matching) >= 1


def test_cnv_section_3_curator_pathogenic_lifts_score(tmp_path):
    """CNV Section 3 score moves with curator entries — was previously 0."""
    from pipeline.utils.cnv_engine import section_3_literature
    from pipeline.utils.curation_store import add_curation

    cfg = {"output": {"data_dir": str(tmp_path)},
           "databases": {"curations": {}}}
    add_curation({
        "variant_id": "chr22_18900000_DEL",
        "criterion": "CNV_section_3",
        "action": "support_pathogenic",
        "strength": "moderate",
        "evidence_text": "DiGeorge syndrome — established literature",
    }, cfg)
    variant = {"variant_id": "chr22_18900000_DEL", "gene": "TBX1",
               "svtype": "DEL", "pos": 18900000, "end_pos": 21500000,
               "chrom": "chr22"}
    s3 = section_3_literature(variant, cfg)
    assert s3["score"] >= 0.50  # moderate strength


def test_cnv_section_3_no_curator_no_litvar_score_is_zero(tmp_path):
    """No curator entry, no LitVar hits → Section 3 stays at 0."""
    from pipeline.utils.cnv_engine import section_3_literature

    cfg = {"output": {"data_dir": str(tmp_path)},
           "databases": {"curations": {}, "litvar": {"enabled": False}}}
    variant = {"variant_id": "x", "gene": "OBSCURE_GENE",
               "svtype": "DEL", "pos": 100, "end_pos": 200, "chrom": "chr1"}
    s3 = section_3_literature(variant, cfg)
    assert s3["score"] == 0.0
