"""P4.1: NCBI LitVar2 client (offline mocks)."""

from unittest.mock import patch, MagicMock


def test_no_input_returns_clear_error():
    from pipeline.utils.litvar_client import query_litvar
    res = query_litvar("", "", "")
    assert res["found"] is False
    assert "must provide" in (res["error"] or "")


def test_litvar_success_path():
    from pipeline.utils import litvar_client
    fake_response = MagicMock()
    fake_response.raise_for_status = MagicMock()
    fake_response.json.return_value = {
        "variant_id": "litvar@chr7:140753336",
        "pmids": [29667901, 30209399, 31690835, 36413997, 37733941],
        "publication_count": 142,
    }
    with patch("requests.get", return_value=fake_response):
        res = litvar_client.query_litvar("BRAF", "p.V600E")

    assert res["found"] is True
    assert res["n_pubmed_ids"] == 5
    assert res["publication_count"] == 142
    assert 29667901 in res["pubmed_ids"]
    assert "litvar2" in (res["source_url"] or "")


def test_litvar_truncates_pmid_list():
    from pipeline.utils import litvar_client
    fake_response = MagicMock()
    fake_response.raise_for_status = MagicMock()
    fake_response.json.return_value = {
        "variant_id": "x",
        "pmids": list(range(1, 100)),
        "publication_count": 99,
    }
    with patch("requests.get", return_value=fake_response):
        res = litvar_client.query_litvar("BRAF", "p.V600E", max_pmids=10)
    assert len(res["pubmed_ids"]) == 10
    assert res["publication_count"] == 99   # full count preserved


def test_litvar_failure_captured():
    from pipeline.utils import litvar_client
    with patch("requests.get", side_effect=RuntimeError("connection refused")):
        res = litvar_client.query_litvar("BRAF", "p.V600E")
    assert res["found"] is False
    assert "connection refused" in (res["error"] or "")
