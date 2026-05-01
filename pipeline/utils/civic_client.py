"""CIViC (https://civicdb.org/) client for somatic variant interpretation.

CIViC is open-data, no auth required. We query the GraphQL API by
gene symbol + protein change (e.g. BRAF p.V600E) and return the highest-
confidence evidence items grouped by clinical significance.

Used by `amp_engine.classify_somatic_variant()` to assign AMP/ASCO/CAP 2017
tiers based on community-curated clinical evidence.
"""

import logging
import time
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

DEFAULT_API = "https://civicdb.org/api/graphql"

# How CIViC evidence levels map to AMP tiers.
# CIViC level A (validated) → high-quality professional guidelines / FDA-approved → Tier I
# CIViC level B (clinical, multiple studies) → Tier I or II depending on context
# CIViC level C (case studies) → Tier II
# CIViC level D (preclinical) → Tier III
# CIViC level E (inferential) → Tier IV/III
EVIDENCE_LEVEL_TO_TIER = {
    "A": "I",
    "B": "I",
    "C": "II",
    "D": "III",
    "E": "III",
}


def query_civic(gene: str, protein_change: str,
                api_url: str = DEFAULT_API,
                timeout_s: float = 8.0) -> Dict:
    """Look up a gene + protein change in CIViC.

    Returns:
        {
          "found": bool,
          "evidence_count": int,
          "highest_level": str,            # 'A'..'E' or None
          "amp_tier": str,                 # 'I'..'IV' or None
          "drug_targets": List[str],
          "evidence_items": List[dict],    # truncated summary
          "source_url": str,
          "error": Optional[str],
          "latency_ms": int,
        }
    """
    out = {
        "found": False, "evidence_count": 0, "highest_level": None,
        "amp_tier": None, "drug_targets": [], "evidence_items": [],
        "source_url": None, "error": None, "latency_ms": 0,
    }
    if not gene:
        return out

    try:
        import requests
    except ImportError as e:
        out["error"] = f"requests not installed: {e}"
        return out

    query = """
    query VariantSearch($gene: String!, $name: String!) {
      variants(geneName: $gene, name: $name, first: 20) {
        nodes {
          id
          name
          gene { name }
          evidenceItems(first: 50, status: ACCEPTED) {
            totalCount
            nodes {
              evidenceLevel
              clinicalSignificance
              evidenceDirection
              evidenceType
              drugs { name }
              disease { name }
            }
          }
        }
      }
    }
    """

    t0 = time.time()
    try:
        resp = requests.post(api_url,
                              json={"query": query,
                                    "variables": {"gene": gene,
                                                  "name": protein_change}},
                              timeout=timeout_s)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        out["error"] = str(e)
        out["latency_ms"] = int((time.time() - t0) * 1000)
        return out

    out["latency_ms"] = int((time.time() - t0) * 1000)
    nodes = (((data or {}).get("data", {}) or {}).get("variants", {}) or {}).get("nodes", []) or []
    if not nodes:
        return out

    # Pick the first (best) match
    var = nodes[0]
    evidence = var.get("evidenceItems", {}).get("nodes", []) or []
    out["found"] = True
    out["evidence_count"] = var.get("evidenceItems", {}).get("totalCount") or len(evidence)
    out["source_url"] = f"https://civicdb.org/links/variants/{var.get('id')}"

    drugs = set()
    levels = []
    items_summary = []
    for ev in evidence:
        lvl = ev.get("evidenceLevel")
        if lvl:
            levels.append(lvl)
        for d in ev.get("drugs") or []:
            n = d.get("name")
            if n:
                drugs.add(n)
        items_summary.append({
            "level": lvl,
            "significance": ev.get("clinicalSignificance"),
            "type": ev.get("evidenceType"),
            "disease": (ev.get("disease") or {}).get("name"),
        })

    if levels:
        # Pick the alphabetically lowest (A < B < C < D < E => most clinically actionable)
        out["highest_level"] = sorted(levels)[0]
        out["amp_tier"] = EVIDENCE_LEVEL_TO_TIER.get(out["highest_level"])

    out["drug_targets"] = sorted(drugs)
    out["evidence_items"] = items_summary[:10]  # truncate for parquet rows
    return out


def is_enabled(databases_config: dict) -> bool:
    civic = (databases_config or {}).get("civic", {}) or {}
    return bool(civic.get("enabled", True))  # default-enabled
