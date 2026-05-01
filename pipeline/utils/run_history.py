"""Per-run history tracking — one row per pipeline run.

Stores run metadata so the dashboard can show "the last 100 runs" with
filterable columns (sample, mode, status, duration). Two backends:

  - Local JSONL (`data/run_history.jsonl`) — used by default, works
    offline, persists across pipeline runs in the same workspace
  - BigQuery (`databases.runs.bq_table` configured) — used in production
    deployments where multiple researchers share state

This is a thin parallel of `curation_store.py` — same shape, different
schema.
"""

from __future__ import annotations

import json
import logging
import os
import time
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def _local_path(config: dict) -> str:
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    data_dir = (config.get("output", {}) or {}).get("data_dir", "data")
    if not os.path.isabs(data_dir):
        data_dir = os.path.join(here, data_dir)
    return os.path.join(data_dir, "run_history.jsonl")


def _bq_table(config: dict) -> Optional[str]:
    return ((config.get("databases", {}) or {})
             .get("runs", {}) or {}).get("bq_table")


def record_run_start(run_id: str, config: dict, manifest: dict) -> None:
    """Append a 'started' entry. Call this at the top of run_pipeline.py."""
    entry = {
        "run_id": run_id,
        "event": "started",
        "ts": time.time(),
        "analysis_mode": manifest.get("analysis_mode"),
        "adapter": manifest.get("adapter"),
        "git_sha": manifest.get("git_sha"),
        "config_hash": manifest.get("config_hash"),
    }
    _append(entry, config)


def record_run_end(run_id: str, status: str, summary: Dict, config: dict) -> None:
    """Append a 'finished' entry once all stages complete (or one fails)."""
    entry = {
        "run_id": run_id,
        "event": "finished",
        "ts": time.time(),
        "status": status,  # success | failure | partial
        "summary": summary,
    }
    _append(entry, config)


def list_runs(config: dict, limit: int = 100) -> List[Dict]:
    """Return the most-recent runs (descending by ts). Aggregates start +
    end events into one row per run_id."""
    bq_table = _bq_table(config)
    raw = (_read_from_bq(bq_table, limit * 4) if bq_table
            else _read_from_local(_local_path(config), limit * 4))
    # Aggregate
    by_run: Dict[str, Dict] = {}
    for e in raw:
        rid = e.get("run_id")
        if not rid:
            continue
        if rid not in by_run:
            by_run[rid] = {"run_id": rid}
        if e.get("event") == "started":
            by_run[rid]["started_at"] = e.get("ts")
            for k in ("analysis_mode", "adapter", "git_sha", "config_hash"):
                if e.get(k):
                    by_run[rid][k] = e[k]
        elif e.get("event") == "finished":
            by_run[rid]["finished_at"] = e.get("ts")
            by_run[rid]["status"] = e.get("status")
            by_run[rid]["summary"] = e.get("summary")
    runs = sorted(by_run.values(),
                   key=lambda r: r.get("started_at", 0), reverse=True)
    # Compute duration
    for r in runs:
        if r.get("started_at") and r.get("finished_at"):
            r["duration_seconds"] = r["finished_at"] - r["started_at"]
        else:
            r["duration_seconds"] = None
    return runs[:limit]


def _append(entry: Dict, config: dict) -> None:
    bq_table = _bq_table(config)
    if bq_table:
        _write_to_bq(entry, bq_table)
    else:
        _append_local(entry, _local_path(config))


def _append_local(entry: Dict, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(entry, default=str) + "\n")


def _read_from_local(path: str, limit: int) -> List[Dict]:
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as f:
        for line in f:
            try:
                out.append(json.loads(line))
            except Exception:
                continue
    out.sort(key=lambda x: x.get("ts", 0), reverse=True)
    return out[:limit]


def _write_to_bq(entry: Dict, bq_table: str) -> None:
    try:
        from google.cloud import bigquery
        client = bigquery.Client()
        # summary is a dict; serialize for the JSON column
        row = dict(entry)
        if "summary" in row and not isinstance(row["summary"], str):
            row["summary"] = json.dumps(row["summary"], default=str)
        errors = client.insert_rows_json(bq_table, [row])
        if errors:
            logger.warning(f"Run history BQ insert failed: {errors}")
    except Exception as e:
        logger.warning(f"Run history BQ insert error: {e}")


def _read_from_bq(bq_table: str, limit: int) -> List[Dict]:
    try:
        from google.cloud import bigquery
        client = bigquery.Client()
        query = f"SELECT * FROM `{bq_table}` ORDER BY ts DESC LIMIT {limit}"
        rows = list(client.query(query).result())
        return [dict(r) for r in rows]
    except Exception as e:
        logger.warning(f"Run history BQ read error: {e}")
        return []
