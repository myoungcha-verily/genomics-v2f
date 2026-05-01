"""P7: FHIR + PDF + CSV exporters."""

import json
from pathlib import Path

import pandas as pd
import pytest


def _sample_df():
    return pd.DataFrame([
        {"variant_id": "chr17_43106478_T_G", "chrom": "chr17", "pos": 43106478,
         "ref": "T", "alt": "G", "gene": "BRCA1", "hgvs_p": "p.Cys61Gly",
         "hgvs_c": "c.181T>G",
         "consequence": "missense_variant", "severity": "MODERATE",
         "acmg_classification": "Pathogenic",
         "acmg_criteria": "PVS1, PM2, PP3, PP5",
         "acmg_tier": 1,
         "bayesian_posterior_prob": 0.97,
         "bayesian_classification": "Pathogenic",
         "amp_tier": None, "amp_drug_targets": None,
         "cnv_classification": None, "cnv_score": None,
         "clinvar_classification": "Pathogenic", "clinvar_review_stars": 4,
         "gnomad_af": 0.000005,
         "cadd_phred": 28.4, "revel": 0.842,
         "phenotype_match_score": 0.95},
        {"variant_id": "chr13_32340301_A_G", "chrom": "chr13", "pos": 32340301,
         "ref": "A", "alt": "G", "gene": "BRCA2", "hgvs_p": "",
         "hgvs_c": "c.7242A>G",
         "consequence": "synonymous", "severity": "LOW",
         "acmg_classification": "Likely benign", "acmg_criteria": "BS1, BP4",
         "acmg_tier": 3,
         "bayesian_posterior_prob": 0.05,
         "bayesian_classification": "Likely benign",
         "amp_tier": None, "amp_drug_targets": None,
         "cnv_classification": None, "cnv_score": None,
         "clinvar_classification": "Benign", "clinvar_review_stars": 3,
         "gnomad_af": 0.058, "cadd_phred": 8.2, "revel": 0.045,
         "phenotype_match_score": 0.0},
    ])


def test_csv_exporter_writes_expected_columns(tmp_path):
    from pipeline.utils.csv_exporter import export_csv
    df = _sample_df()
    out = tmp_path / "x_variants.csv"
    export_csv(df, str(out))
    assert out.exists()
    head = out.read_text().splitlines()[0]
    # Spot-check key columns are in the header
    for col in ("variant_id", "gene", "acmg_classification",
                "bayesian_posterior_prob"):
        assert col in head


def test_csv_exporter_handles_missing_columns(tmp_path):
    """Different analysis modes produce different columns; missing ones
    should be skipped silently."""
    from pipeline.utils.csv_exporter import export_csv
    df = pd.DataFrame([{"variant_id": "x", "chrom": "chr1", "pos": 100,
                         "ref": "A", "alt": "T",
                         "acmg_classification": "VUS"}])
    out = tmp_path / "minimal.csv"
    export_csv(df, str(out))
    assert out.exists()
    # No 'gene', 'amp_tier', etc. — but the file should still be written
    head = out.read_text().splitlines()[0]
    assert "variant_id" in head
    assert "gene" not in head  # not present in df


def test_fhir_exporter_produces_valid_bundle(tmp_path):
    from pipeline.utils.fhir_exporter import export_fhir_diagnostic_report
    df = _sample_df()
    out = tmp_path / "x_fhir.json"
    export_fhir_diagnostic_report(df, "PATIENT_X", str(out),
                                    manifest={"git_sha": "abc"})
    payload = json.loads(out.read_text())
    assert payload["resourceType"] == "Bundle"
    assert payload["type"] == "collection"

    # Find the DiagnosticReport
    reports = [e["resource"] for e in payload["entry"]
                if e["resource"]["resourceType"] == "DiagnosticReport"]
    assert len(reports) == 1
    rep = reports[0]
    assert rep["status"] == "preliminary"
    assert "PRELIMINARY" in rep["conclusion"]
    assert "extension" in rep  # manifest stamped

    # Only Pathogenic gets reported (LB does not)
    obs = [e["resource"] for e in payload["entry"]
            if e["resource"]["resourceType"] == "Observation"]
    assert len(obs) == 1
    classification = obs[0]["valueCodeableConcept"]["text"]
    assert classification == "Pathogenic"


def test_fhir_exporter_includes_amp_and_cnv_reportables(tmp_path):
    from pipeline.utils.fhir_exporter import export_fhir_diagnostic_report
    df = pd.DataFrame([
        {"variant_id": "v1", "chrom": "chr7", "pos": 100, "ref": "A", "alt": "T",
         "gene": "BRAF", "hgvs_p": "p.V600E",
         "amp_tier": "I", "amp_drug_targets": "Vemurafenib",
         "amp_evidence": "FDA-approved", "acmg_classification": None,
         "cnv_classification": None},
        {"variant_id": "v2", "chrom": "chr22", "pos": 18900000, "ref": "N",
         "alt": "<DEL>", "gene": "22q11.2", "svtype": "DEL",
         "cnv_classification": "Pathogenic", "cnv_evidence_summary": "DiGeorge",
         "amp_tier": None, "acmg_classification": None},
    ])
    out = tmp_path / "fhir.json"
    export_fhir_diagnostic_report(df, "PATIENT_Y", str(out))
    payload = json.loads(out.read_text())
    obs = [e["resource"] for e in payload["entry"]
            if e["resource"]["resourceType"] == "Observation"]
    assert len(obs) == 2  # both reportable


def test_fhir_no_reportable_variants_yields_empty_results(tmp_path):
    from pipeline.utils.fhir_exporter import export_fhir_diagnostic_report
    df = pd.DataFrame([{"variant_id": "x", "chrom": "chr1", "pos": 1,
                         "ref": "A", "alt": "T",
                         "acmg_classification": "Benign"}])
    out = tmp_path / "fhir.json"
    export_fhir_diagnostic_report(df, "PATIENT_Z", str(out))
    payload = json.loads(out.read_text())
    reports = [e["resource"] for e in payload["entry"]
                if e["resource"]["resourceType"] == "DiagnosticReport"]
    assert len(reports[0]["result"]) == 0
    assert "0 reportable" in reports[0]["conclusion"]


def test_pdf_renderer_no_op_when_weasyprint_missing(tmp_path):
    """The PDF renderer is best-effort — without weasyprint it returns
    None and logs, but never raises."""
    from pipeline.utils.pdf_renderer import render_html_to_pdf
    html = tmp_path / "report.html"
    html.write_text("<html><body>Test</body></html>")
    # weasyprint is not in the test deps — this should return None gracefully
    result = render_html_to_pdf(str(html))
    # In environments where weasyprint IS installed, the test is still valid
    # (returns the PDF path, not None). Either way it should not raise.
    assert result is None or result.endswith(".pdf")


def test_fhir_exporter_is_enabled_flag():
    from pipeline.utils.fhir_exporter import is_enabled
    assert is_enabled({"output": {"fhir_export": True}}) is True
    assert is_enabled({"output": {}}) is False
    assert is_enabled({}) is False


def test_pdf_renderer_is_enabled_flag():
    from pipeline.utils.pdf_renderer import is_enabled
    assert is_enabled({"output": {"report_format": "pdf"}}) is True
    assert is_enabled({"output": {"report_format": "both"}}) is True
    assert is_enabled({"output": {"pdf_export": True}}) is True
    assert is_enabled({"output": {"report_format": "html"}}) is False
