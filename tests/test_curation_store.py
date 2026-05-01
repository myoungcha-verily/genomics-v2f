"""P4.3: Curation store with local-JSONL fallback."""

import json
from pathlib import Path


def test_add_then_list_local_backend(tmp_path):
    from pipeline.utils.curation_store import add_curation, list_curations
    cfg = {"output": {"data_dir": str(tmp_path)},
           "databases": {"curations": {}}}
    e = add_curation({
        "variant_id": "chr17_43106478_T_G",
        "gene": "BRCA1",
        "criterion": "PS3",
        "action": "support_pathogenic",
        "strength": "strong",
        "evidence_text": "Findlay 2018 saturation editing damaging",
        "pubmed_ids": [30209399],
        "curator_email": "test@v2f.io",
    }, cfg)
    assert e["entry_id"]
    assert e["created_at"]

    items = list_curations("chr17_43106478_T_G", cfg)
    assert len(items) == 1
    assert items[0]["criterion"] == "PS3"
    assert items[0]["action"] == "support_pathogenic"


def test_missing_required_fields_raises(tmp_path):
    from pipeline.utils.curation_store import add_curation
    cfg = {"output": {"data_dir": str(tmp_path)},
           "databases": {"curations": {}}}
    try:
        add_curation({"gene": "BRCA1"}, cfg)
        assert False, "should have raised"
    except ValueError as e:
        assert "missing" in str(e)


def test_overrides_for_variant_picks_latest_per_criterion(tmp_path):
    from pipeline.utils.curation_store import (
        add_curation, get_overrides_for_variant
    )
    cfg = {"output": {"data_dir": str(tmp_path)},
           "databases": {"curations": {}}}
    add_curation({"variant_id": "v1", "criterion": "PS3",
                   "action": "support_pathogenic", "strength": "strong",
                   "created_at": 1000}, cfg)
    add_curation({"variant_id": "v1", "criterion": "PS3",
                   "action": "support_benign", "strength": "supporting",
                   "created_at": 2000}, cfg)
    overrides = get_overrides_for_variant("v1", cfg)
    # Only one entry per criterion, the latest wins
    assert len(overrides) == 1
    assert overrides[0]["action"] == "support_benign"
    assert overrides[0]["created_at"] == 2000


def test_list_filtered_by_variant(tmp_path):
    from pipeline.utils.curation_store import add_curation, list_curations
    cfg = {"output": {"data_dir": str(tmp_path)},
           "databases": {"curations": {}}}
    add_curation({"variant_id": "v1", "criterion": "PS3",
                   "action": "support_pathogenic"}, cfg)
    add_curation({"variant_id": "v2", "criterion": "PS3",
                   "action": "support_pathogenic"}, cfg)
    assert len(list_curations("v1", cfg)) == 1
    assert len(list_curations("v2", cfg)) == 1
    assert len(list_curations(None, cfg)) == 2
