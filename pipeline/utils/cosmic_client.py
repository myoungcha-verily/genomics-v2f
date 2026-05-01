"""COSMIC client (BigQuery-based).

COSMIC's mutation census is mirrored on internal BQ tables for some
deployments. Default-disabled; if a `cosmic.bq_table` is configured,
this client looks up the variant by chrom/pos/ref/alt and returns
recurrence + tissue distribution metadata.
"""

import logging
import time
from typing import Dict

logger = logging.getLogger(__name__)


def query_cosmic(chrom: str, pos: int, ref: str, alt: str,
                  bq_table: str,
                  timeout_s: float = 10.0) -> Dict:
    out = {
        "found": False, "occurrences": 0, "tissue_distribution": [],
        "cosmic_id": None, "source_url": None, "error": None,
        "latency_ms": 0,
    }
    if not bq_table:
        out["error"] = "no_table"
        return out

    try:
        from google.cloud import bigquery
    except ImportError as e:
        out["error"] = f"bigquery client not installed: {e}"
        return out

    chrom_bare = str(chrom).replace("chr", "")
    query = f"""
    SELECT genomic_mutation_id AS cosmic_id,
           COUNT(*) AS occurrences,
           ARRAY_AGG(DISTINCT primary_site IGNORE NULLS LIMIT 10) AS tissues
    FROM `{bq_table}`
    WHERE chromosome = @chrom AND start_position = @pos
      AND mutation_cds_ref = @ref AND mutation_cds_alt = @alt
    GROUP BY genomic_mutation_id
    """

    t0 = time.time()
    try:
        client = bigquery.Client()
        params = [
            bigquery.ScalarQueryParameter("chrom", "STRING", chrom_bare),
            bigquery.ScalarQueryParameter("pos", "INT64", int(pos)),
            bigquery.ScalarQueryParameter("ref", "STRING", ref),
            bigquery.ScalarQueryParameter("alt", "STRING", alt),
        ]
        job_config = bigquery.QueryJobConfig(query_parameters=params)
        rows = list(client.query(query, job_config=job_config,
                                  timeout=timeout_s).result())
    except Exception as e:
        out["error"] = str(e)
        out["latency_ms"] = int((time.time() - t0) * 1000)
        return out

    out["latency_ms"] = int((time.time() - t0) * 1000)
    if not rows:
        return out
    row = rows[0]
    out["found"] = True
    out["occurrences"] = int(row["occurrences"])
    out["tissue_distribution"] = list(row["tissues"]) if row["tissues"] else []
    out["cosmic_id"] = row["cosmic_id"]
    if out["cosmic_id"]:
        out["source_url"] = (
            f"https://cancer.sanger.ac.uk/cosmic/search?q={out['cosmic_id']}"
        )
    return out


def is_enabled(databases_config: dict) -> bool:
    cosmic = (databases_config or {}).get("cosmic", {}) or {}
    return bool(cosmic.get("enabled", False)) and bool(cosmic.get("bq_table"))
