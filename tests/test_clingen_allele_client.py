"""P4.2: ClinGen Allele Registry client (offline mocks)."""

from unittest.mock import patch, MagicMock


def test_unknown_chrom_returns_error():
    from pipeline.utils.clingen_allele_client import query_clingen_allele
    res = query_clingen_allele("chrXX", 100, "A", "T")
    assert res["found"] is False
    assert "unsupported" in (res["error"] or "")


def test_missing_input_returns_error():
    from pipeline.utils.clingen_allele_client import query_clingen_allele
    res = query_clingen_allele("", 0, "", "")
    assert res["found"] is False


def test_clingen_success_path():
    from pipeline.utils import clingen_allele_client
    fake = MagicMock()
    fake.raise_for_status = MagicMock()
    fake.json.return_value = {
        "@id": "http://reg.clinicalgenome.org/allele/CA127401",
        "externalRecords": {
            "ClinVarVariations": [
                {"RCV": [{"clinicalSignificance": "Pathogenic"}],
                 "reviewStatus": "criteria_provided"},
            ],
        },
    }
    with patch("requests.get", return_value=fake):
        res = clingen_allele_client.query_clingen_allele(
            "chr17", 43106478, "T", "G", reference_genome="GRCh38")
    assert res["found"] is True
    assert res["ca_id"] == "CA127401"
    assert res["current_clinvar_class"] == "Pathogenic"
    assert res["n_assertions"] >= 1


def test_refseq_for_chrom():
    from pipeline.utils.clingen_allele_client import _refseq_for_chrom
    assert _refseq_for_chrom("17", "GRCh38") == "NC_000017.11"
    assert _refseq_for_chrom("chr17", "GRCh38") == "NC_000017.11"
    assert _refseq_for_chrom("17", "GRCh37") == "NC_000017.10"
    assert _refseq_for_chrom("chrZZ", "GRCh38") is None
