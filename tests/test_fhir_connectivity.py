"""Pre-M1 polish tests: FHIR connectivity probe + dashboard endpoint."""

from unittest.mock import MagicMock, patch


def test_unconfigured_returns_clear_error(repo_root):
    from pipeline.utils.fhir_phenotype import test_fhir_connectivity
    res = test_fhir_connectivity("", "")
    assert res["ok"] is False
    assert "must both be set" in (res["error"] or "")
    assert res["latency_ms"] == 0


def test_partial_config_treated_as_unconfigured(repo_root):
    from pipeline.utils.fhir_phenotype import test_fhir_connectivity
    res = test_fhir_connectivity("only-project", "")
    assert res["ok"] is False
    assert "must both be set" in (res["error"] or "")


def test_query_failure_captured_with_error(repo_root, monkeypatch):
    """When BQ raises, the error message lands in the result without raising."""
    from pipeline.utils import fhir_phenotype as fp

    fake_client = MagicMock()
    fake_client.query.side_effect = RuntimeError("table not found")

    fake_module = MagicMock()
    fake_module.Client.return_value = fake_client
    fake_module.QueryJobConfig.return_value = MagicMock()

    monkeypatch.setattr(
        "google.cloud.bigquery", fake_module, raising=False
    )
    import sys
    sys.modules["google.cloud.bigquery"] = fake_module

    res = fp.test_fhir_connectivity("p", "d", "Condition")
    assert res["ok"] is False
    assert "table not found" in (res["error"] or "")
    assert res["latency_ms"] >= 0


def test_dashboard_endpoint_returns_payload(repo_root, tmp_path, monkeypatch):
    """POST /api/phenotype/test wires through to the connectivity probe."""
    import shutil

    shutil.copytree(repo_root / "demo", tmp_path / "demo")
    (tmp_path / "config").mkdir()

    import dashboard.app as dapp
    monkeypatch.setattr(dapp, "PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr(dapp, "CONFIG_PATH",
                        str(tmp_path / "config" / "pipeline_config.yaml"))
    monkeypatch.setattr(dapp, "DEMO_DIR", str(tmp_path / "demo"))

    client = dapp.app.test_client()
    rv = client.post("/api/phenotype/test", json={"fhir_project": "", "fhir_dataset": ""})
    assert rv.status_code == 200
    body = rv.get_json()
    assert body["ok"] is False
    assert "must both be set" in (body["error"] or "")
    assert body["project"] == ""
    assert body["dataset"] == ""


def test_extract_phenotype_status_field(repo_root):
    """extract_patient_phenotype now tags the result with a status field."""
    from pipeline.utils.fhir_phenotype import extract_patient_phenotype

    # Disabled
    p1 = extract_patient_phenotype("X", {"phenotype": {"enabled": False}})
    assert p1["status"] == "disabled"

    # Unconfigured (enabled=true but no project/dataset)
    p2 = extract_patient_phenotype("X", {"phenotype": {"enabled": True}})
    assert p2["status"] == "unconfigured"
