"""Structured JSON logging for cloud-monitoring compatibility.

Replaces the default plain-text logger with a JSON formatter so every
log line carries structured fields: timestamp, level, run_id, stage,
sample, message. Cloud Logging / Loki / Datadog all parse this directly.

Usage:
    from pipeline.utils.structured_logger import setup_structured_logging
    setup_structured_logging(run_id="abc-123", stage="04_acmg_classification")
    logger.info("classified", extra={"variant_id": "...", "duration_ms": 42})

The formatter merges any `extra` dict into the JSON payload so callers
can attach per-event context without ceremony.
"""

from __future__ import annotations

import json
import logging
import sys
import time
from typing import Optional


class JsonFormatter(logging.Formatter):
    """Emits one JSON object per log line."""

    # Standard LogRecord attributes we don't want to repeat
    _RESERVED = {
        "name", "msg", "args", "levelname", "levelno", "pathname", "filename",
        "module", "exc_info", "exc_text", "stack_info", "lineno", "funcName",
        "created", "msecs", "relativeCreated", "thread", "threadName",
        "processName", "process", "getMessage",
    }

    def __init__(self, *, service: str = "v2f-reporter",
                  default_fields: Optional[dict] = None):
        super().__init__()
        self.service = service
        self.default_fields = default_fields or {}

    def format(self, record: logging.LogRecord) -> str:
        payload = {
            "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(record.created)),
            "level": record.levelname,
            "logger": record.name,
            "service": self.service,
            "msg": record.getMessage(),
        }
        payload.update(self.default_fields)

        # Pull extra fields the caller attached via extra={}
        for k, v in record.__dict__.items():
            if k in self._RESERVED or k.startswith("_"):
                continue
            if k in payload:
                continue
            try:
                json.dumps(v, default=str)
                payload[k] = v
            except (TypeError, ValueError):
                payload[k] = str(v)

        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)

        return json.dumps(payload, default=str)


def setup_structured_logging(level: str = "INFO", service: str = "v2f-reporter",
                              **default_fields) -> None:
    """Install the JSON formatter on the root logger.

    `default_fields` are merged into every emitted record (e.g. run_id,
    stage). Idempotent — calling twice replaces the handler.
    """
    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Drop any existing handlers (avoid duplicate output)
    for h in list(root.handlers):
        root.removeHandler(h)

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(JsonFormatter(service=service,
                                          default_fields=default_fields))
    root.addHandler(handler)


def stage_logger(name: str, run_id: str, stage: str) -> logging.LoggerAdapter:
    """Get a per-stage logger adapter. Pre-tags every record with run_id +
    stage so caller doesn't have to."""
    base = logging.getLogger(name)
    return logging.LoggerAdapter(base, {"run_id": run_id, "stage": stage})
