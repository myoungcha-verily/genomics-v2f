"""P9: Structured logger + run history."""

import json
import logging
import time
from pathlib import Path


def test_json_formatter_produces_valid_json():
    """Each log line is one JSON object with the expected fields."""
    from pipeline.utils.structured_logger import JsonFormatter
    fmt = JsonFormatter(service="v2f-reporter",
                          default_fields={"run_id": "abc123"})
    rec = logging.LogRecord(
        name="test", level=logging.INFO, pathname="x.py", lineno=1,
        msg="hello", args=(), exc_info=None,
    )
    rec.stage = "stage_4"
    rec.variant_id = "chr17_43106478_T_G"
    out = fmt.format(rec)
    parsed = json.loads(out)
    assert parsed["msg"] == "hello"
    assert parsed["level"] == "INFO"
    assert parsed["service"] == "v2f-reporter"
    assert parsed["run_id"] == "abc123"
    assert parsed["stage"] == "stage_4"
    assert parsed["variant_id"] == "chr17_43106478_T_G"


def test_setup_structured_logging_replaces_handlers():
    """Calling setup twice should not double-add handlers."""
    from pipeline.utils.structured_logger import setup_structured_logging
    setup_structured_logging(run_id="x")
    n1 = len(logging.getLogger().handlers)
    setup_structured_logging(run_id="y")
    n2 = len(logging.getLogger().handlers)
    assert n1 == n2 == 1


def test_run_history_local_round_trip(tmp_path):
    from pipeline.utils.run_history import (
        record_run_start, record_run_end, list_runs
    )
    cfg = {"output": {"data_dir": str(tmp_path)}, "databases": {}}
    manifest = {"analysis_mode": "germline", "adapter": "single_sample",
                 "git_sha": "abc1234567", "config_hash": "sha256:xyz"}
    record_run_start("run-1", cfg, manifest)
    time.sleep(0.01)
    record_run_end("run-1", "success",
                    {"completed_stages": [1, 2, 3, 4, 5, 6, 7]}, cfg)

    runs = list_runs(cfg, limit=10)
    assert len(runs) == 1
    r = runs[0]
    assert r["run_id"] == "run-1"
    assert r["status"] == "success"
    assert r["analysis_mode"] == "germline"
    assert r["duration_seconds"] is not None
    assert r["duration_seconds"] >= 0


def test_run_history_orders_descending(tmp_path):
    from pipeline.utils.run_history import record_run_start, list_runs
    cfg = {"output": {"data_dir": str(tmp_path)}, "databases": {}}
    for i in range(3):
        record_run_start(f"run-{i}", cfg,
                          {"analysis_mode": "germline"})
        time.sleep(0.01)
    runs = list_runs(cfg, limit=10)
    assert [r["run_id"] for r in runs] == ["run-2", "run-1", "run-0"]
