"""NCBI LitVar2 client (https://www.ncbi.nlm.nih.gov/research/litvar2-api/).

LitVar2 links variants to PubMed articles. Free, no auth required.
Used to surface a "literature" column on every variant in the V2F
output, plus to feed CNV Section 3 evidence (Riggs 2019) automatically
where prior code left it at zero.

Lookup paths in this client:
  - by HGVS expression (gene + protein change), e.g. BRAF p.V600E
  - by RSID, e.g. rs113488022
  - by variant ID via ClinGen Allele Registry (caller resolves first)
"""

from __future__ import annotations

import logging
import time
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

DEFAULT_API = "https://www.ncbi.nlm.nih.gov/research/litvar2-api"


def query_litvar(gene: str, hgvs_p: str = "", rsid: str = "",
                  api_url: str = DEFAULT_API,
                  timeout_s: float = 8.0,
                  max_pmids: int = 25) -> Dict:
    """Query LitVar2 for PubMed articles linking to a variant.

    Provide one of: (gene + hgvs_p) or rsid.

    Returns:
        {
          "found": bool,
          "litvar_id": str,
          "n_pubmed_ids": int,
          "pubmed_ids": list[int],   # truncated to max_pmids
          "publication_count": int,  # total in LitVar
          "search_term": str,        # what we asked for
          "source_url": str,
          "error": Optional[str],
          "latency_ms": int,
        }
    """
    out = {
        "found": False, "litvar_id": None, "n_pubmed_ids": 0,
        "pubmed_ids": [], "publication_count": 0,
        "search_term": "", "source_url": None, "error": None,
        "latency_ms": 0,
    }

    if not gene and not rsid:
        out["error"] = "must provide gene+hgvs_p or rsid"
        return out

    try:
        import requests
    except ImportError as e:
        out["error"] = f"requests not installed: {e}"
        return out

    # Step 1: variant lookup by HGVS or RSID
    if rsid:
        search = rsid if rsid.startswith("rs") else f"rs{rsid}"
        out["search_term"] = search
        url = f"{api_url}/variant/get/litvar/%23{search}"
    elif hgvs_p:
        # LitVar accepts HGVS protein change with gene prefix
        search = f"{gene}:{hgvs_p}" if hgvs_p.startswith("p.") else f"{gene}:p.{hgvs_p}"
        out["search_term"] = search
        url = f"{api_url}/variant/get?query={search}"
    else:
        out["search_term"] = gene
        url = f"{api_url}/variant/autocomplete?query={gene}"

    t0 = time.time()
    try:
        r = requests.get(url, timeout=timeout_s)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        out["error"] = str(e)
        out["latency_ms"] = int((time.time() - t0) * 1000)
        return out

    out["latency_ms"] = int((time.time() - t0) * 1000)

    # LitVar returns either a single hit or a list; normalize
    hit = None
    if isinstance(data, list) and data:
        hit = data[0]
    elif isinstance(data, dict):
        hit = data

    if not hit:
        return out

    out["found"] = True
    out["litvar_id"] = hit.get("variant_id") or hit.get("id")
    pmids = hit.get("pmids") or hit.get("pubmed_ids") or []
    if not isinstance(pmids, list):
        pmids = []
    out["pubmed_ids"] = [int(p) for p in pmids[:max_pmids] if str(p).isdigit()]
    out["n_pubmed_ids"] = len(out["pubmed_ids"])
    out["publication_count"] = (
        hit.get("publication_count") or hit.get("publications_count") or
        hit.get("nb_publications") or len(pmids)
    )
    if out["litvar_id"]:
        out["source_url"] = (
            f"https://www.ncbi.nlm.nih.gov/research/litvar2/docsum/?variant={out['litvar_id']}"
        )

    return out


def is_enabled(config: dict) -> bool:
    db = (config.get("databases", {}) or {})
    litvar = db.get("litvar", {}) or {}
    return bool(litvar.get("enabled", True))  # default-on (free API)
