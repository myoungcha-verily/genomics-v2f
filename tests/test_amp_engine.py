"""M1: AMP/ASCO/CAP 2017 four-tier classifier picks correct tier per variant.

Tests use mocked CIViC/OncoKB clients to avoid external API dependency.
"""

from unittest.mock import patch


def _braf_v600e():
    return {
        "gene": "BRAF", "hgvs_p": "p.V600E",
        "chrom": "chr7", "pos": 140753336, "ref": "A", "alt": "T",
        "tumor_vaf": 0.42, "is_paired": True,
    }


def test_amp_classifier_tier_i_with_civic_a_level():
    from pipeline.utils import amp_engine

    with patch("pipeline.utils.civic_client.is_enabled", return_value=True), \
         patch("pipeline.utils.oncokb_client.is_enabled", return_value=False), \
         patch("pipeline.utils.cosmic_client.is_enabled", return_value=False), \
         patch("pipeline.utils.civic_client.query_civic") as q_civic:
        q_civic.return_value = {
            "found": True, "evidence_count": 15,
            "highest_level": "A", "amp_tier": "I",
            "drug_targets": ["Vemurafenib", "Dabrafenib"],
            "evidence_items": [], "source_url": "https://civicdb.org/x",
            "error": None, "latency_ms": 100,
        }
        result = amp_engine.classify_somatic_variant(
            _braf_v600e(), {"databases": {"civic": {"enabled": True}}})
    assert result["amp_tier"] == "I"
    assert "Vemurafenib" in result["drug_targets"]
    assert "civic" in result["knowledge_bases_consulted"]


def test_amp_classifier_oncokb_overrides_civic():
    """When OncoKB is enabled and returns a tier, it should be in the consulted list."""
    from pipeline.utils import amp_engine

    with patch("pipeline.utils.civic_client.is_enabled", return_value=True), \
         patch("pipeline.utils.oncokb_client.is_enabled", return_value=True), \
         patch("pipeline.utils.oncokb_client.get_token", return_value="fake-token"), \
         patch("pipeline.utils.cosmic_client.is_enabled", return_value=False), \
         patch("pipeline.utils.oncokb_client.query_oncokb") as q_oncokb, \
         patch("pipeline.utils.civic_client.query_civic") as q_civic:
        q_oncokb.return_value = {
            "found": True, "highest_level": "LEVEL_1", "amp_tier": "I",
            "drug_targets": ["Vemurafenib"], "evidence_items": [],
            "source_url": None, "error": None, "latency_ms": 100,
            "oncogenic": "Oncogenic",
        }
        q_civic.return_value = {
            "found": True, "highest_level": "B", "amp_tier": "II",
            "drug_targets": [], "evidence_items": [], "source_url": None,
            "error": None, "latency_ms": 50,
        }
        result = amp_engine.classify_somatic_variant(
            _braf_v600e(),
            {"databases": {"oncokb": {"enabled": True}, "civic": {"enabled": True}},
             "acmg_amp": {"knowledge_base_priority": ["oncokb", "civic", "cosmic"]}})

    # Both KBs returned a tier; best (Tier I from OncoKB) wins
    assert result["amp_tier"] == "I"
    assert "oncokb" in result["knowledge_bases_consulted"]
    assert "civic" in result["knowledge_bases_consulted"]


def test_amp_classifier_no_kb_hits_returns_tier_iii():
    """No KB hit → default to Tier III (somatic VUS)."""
    from pipeline.utils import amp_engine

    with patch("pipeline.utils.civic_client.is_enabled", return_value=True), \
         patch("pipeline.utils.oncokb_client.is_enabled", return_value=False), \
         patch("pipeline.utils.cosmic_client.is_enabled", return_value=False), \
         patch("pipeline.utils.civic_client.query_civic") as q_civic:
        q_civic.return_value = {
            "found": False, "evidence_count": 0, "highest_level": None,
            "amp_tier": None, "drug_targets": [], "evidence_items": [],
            "source_url": None, "error": None, "latency_ms": 100,
        }
        result = amp_engine.classify_somatic_variant(
            {"gene": "OBSCURE_GENE", "hgvs_p": "p.X1Y", "tumor_vaf": 0.3,
             "is_paired": True},
            {"databases": {"civic": {"enabled": True}}})
    assert result["amp_tier"] == "III"


def test_amp_classifier_high_vaf_unpaired_warns():
    """Tumor-only VAF > 0.45 warns about possible germline contamination."""
    from pipeline.utils import amp_engine

    with patch("pipeline.utils.civic_client.is_enabled", return_value=True), \
         patch("pipeline.utils.oncokb_client.is_enabled", return_value=False), \
         patch("pipeline.utils.cosmic_client.is_enabled", return_value=False), \
         patch("pipeline.utils.civic_client.query_civic") as q_civic:
        q_civic.return_value = {
            "found": True, "highest_level": "A", "amp_tier": "I",
            "drug_targets": ["Vemurafenib"], "evidence_items": [],
            "source_url": None, "error": None, "latency_ms": 50,
        }
        v = _braf_v600e()
        v["tumor_vaf"] = 0.50
        v["is_paired"] = False
        result = amp_engine.classify_somatic_variant(
            v, {"databases": {"civic": {"enabled": True}}})
    assert result["amp_tier"] == "I"
    assert "germline" in result["evidence_summary"].lower()


def test_oncokb_no_token_graceful_no_op():
    """OncoKB client without a token returns error='no_token' and the engine
    should still produce a valid tier (from CIViC alone)."""
    from pipeline.utils import oncokb_client
    res = oncokb_client.query_oncokb("BRAF", "p.V600E", api_token="")
    assert res["found"] is False
    assert res["error"] == "no_token"


def test_tier_label_human_readable():
    from pipeline.utils.amp_engine import tier_label
    assert "Strong" in tier_label("I")
    assert "Potential" in tier_label("II")
    assert "Unknown" in tier_label("III")
    assert "Benign" in tier_label("IV")
