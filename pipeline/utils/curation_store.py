"""Per-workspace curation store for manual literature evidence.

Researchers add literature-backed evidence per variant via the Curate
tab in the dashboard. Entries are stored in a BigQuery table when one
is configured; otherwise in a local JSONL file (`data/curations.jsonl`)
so the demo flow works without BQ setup.

Each entry has:
  - variant_id (chrom-pos-ref-alt)
  - gene + hgvs_p
  - criterion (PS3, PS4, PP1, BS3, BS4, BP6, "CNV_section_3", ...)
  - action (support_pathogenic | support_benign | no_call)
  - strength (very_strong | strong | moderate | supporting | standalone)
  - evidence_text
  - pubmed_ids (array)
  - curator_email
  - created_at (epoch seconds)
  - replaces (CA-ID of a prior entry this amends)

Curator entries override the automated engine's call for that variant
+ criterion. The classifier tags overridden calls with `curated: true`.
"""

from __future__ import annotations

import json
import logging
import os
import time
import uuid
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def _bq_table(config: dict) -> Optional[str]:
    return ((config.get("databases", {}) or {}).get("curations", {}) or {}).get("bq_table")


def _local_path(config: dict) -> str:
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    data_dir = (config.get("output", {}) or {}).get("data_dir", "data")
    if not os.path.isabs(data_dir):
        data_dir = os.path.join(here, data_dir)
    return os.path.join(data_dir, "curations.jsonl")


def add_curation(entry: Dict, config: dict) -> Dict:
    """Persist a curation entry. Returns the stored entry (with id + ts filled)."""
    entry = dict(entry)
    entry.setdefault("entry_id", str(uuid.uuid4()))
    entry.setdefault("created_at", int(time.time()))

    required = ("variant_id", "criterion", "action")
    missing = [k for k in required if not entry.get(k)]
    if missing:
        raise ValueError(f"missing required curation fields: {missing}")

    bq_table = _bq_table(config)
    if bq_table:
        _write_to_bq(entry, bq_table)
    else:
        _append_to_local(entry, _local_path(config))

    return entry


def list_curations(variant_id: Optional[str], config: dict) -> List[Dict]:
    """List curation entries, optionally filtered to a single variant."""
    bq_table = _bq_table(config)
    if bq_table:
        return _read_from_bq(bq_table, variant_id)
    return _read_from_local(_local_path(config), variant_id)


def get_overrides_for_variant(variant_id: str, config: dict) -> List[Dict]:
    """Get the most-recent curation override per criterion for one variant."""
    entries = list_curations(variant_id, config)
    by_criterion: Dict[str, Dict] = {}
    for e in sorted(entries, key=lambda x: x.get("created_at", 0)):
        by_criterion[e["criterion"]] = e   # latest wins
    return list(by_criterion.values())


# ---- BQ backend ---------------------------------------------------

_BQ_SCHEMA = [
    ("entry_id", "STRING"),
    ("variant_id", "STRING"),
    ("gene", "STRING"),
    ("hgvs_p", "STRING"),
    ("criterion", "STRING"),
    ("action", "STRING"),
    ("strength", "STRING"),
    ("evidence_text", "STRING"),
    ("pubmed_ids", "STRING"),  # comma-joined for simplicity
    ("curator_email", "STRING"),
    ("created_at", "INT64"),
    ("replaces", "STRING"),
]


def _write_to_bq(entry: Dict, bq_table: str) -> None:
    from google.cloud import bigquery
    client = bigquery.Client()
    row = {
        "entry_id": entry.get("entry_id"),
        "variant_id": entry.get("variant_id"),
        "gene": entry.get("gene", ""),
        "hgvs_p": entry.get("hgvs_p", ""),
        "criterion": entry.get("criterion", ""),
        "action": entry.get("action", ""),
        "strength": entry.get("strength", "supporting"),
        "evidence_text": entry.get("evidence_text", ""),
        "pubmed_ids": ",".join(str(p) for p in (entry.get("pubmed_ids") or [])),
        "curator_email": entry.get("curator_email", ""),
        "created_at": entry.get("created_at"),
        "replaces": entry.get("replaces", ""),
    }
    errors = client.insert_rows_json(bq_table, [row])
    if errors:
        raise RuntimeError(f"BQ insert failed: {errors}")


def _read_from_bq(bq_table: str, variant_id: Optional[str]) -> List[Dict]:
    from google.cloud import bigquery
    client = bigquery.Client()
    where = f"WHERE variant_id = @vid" if variant_id else ""
    params = ([bigquery.ScalarQueryParameter("vid", "STRING", variant_id)]
               if variant_id else [])
    query = f"SELECT * FROM `{bq_table}` {where} ORDER BY created_at DESC"
    job_config = bigquery.QueryJobConfig(query_parameters=params)
    rows = list(client.query(query, job_config=job_config).result())
    return [dict(r) for r in rows]


# ---- Local-JSONL backend (fallback) -------------------------------

def _append_to_local(entry: Dict, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(entry, default=str) + "\n")


def _read_from_local(path: str, variant_id: Optional[str]) -> List[Dict]:
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as f:
        for line in f:
            try:
                e = json.loads(line)
            except Exception:
                continue
            if variant_id and e.get("variant_id") != variant_id:
                continue
            out.append(e)
    out.sort(key=lambda x: x.get("created_at", 0), reverse=True)
    return out
