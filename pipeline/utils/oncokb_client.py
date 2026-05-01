"""OncoKB (https://www.oncokb.org/) client for somatic variant interpretation.

OncoKB requires a commercial API token for full access. Without a token
this client logs once and returns no-op results — the AMP engine then
falls back to CIViC alone. With a token it returns OncoKB's curated
oncogenic / clinical-actionability classification.
"""

import logging
import time
from typing import Dict, Optional

logger = logging.getLogger(__name__)

DEFAULT_API = "https://www.oncokb.org/api/v1"

# OncoKB level → AMP tier. OncoKB has its own level system (1-4 + R1/R2)
# — these are rough mappings used when AMP-tier output is needed.
ONCOKB_LEVEL_TO_TIER = {
    "LEVEL_1": "I",
    "LEVEL_2": "I",
    "LEVEL_3A": "II",
    "LEVEL_3B": "II",
    "LEVEL_4": "III",
    "LEVEL_R1": "I",
    "LEVEL_R2": "II",
}


def query_oncokb(gene: str, protein_change: str,
                  api_token: str,
                  api_url: str = DEFAULT_API,
                  timeout_s: float = 8.0) -> Dict:
    """Look up a gene + protein change in OncoKB.

    Returns the same shape as `civic_client.query_civic()` for
    compatibility. If api_token is empty, returns a no-op result with
    error='no_token' (not raised — graceful degradation).
    """
    out = {
        "found": False, "evidence_count": 0, "highest_level": None,
        "amp_tier": None, "drug_targets": [], "evidence_items": [],
        "source_url": None, "error": None, "latency_ms": 0,
        "oncogenic": None,
    }
    if not gene:
        return out
    if not api_token:
        # Default-disabled path: log once at module level (caller logs it),
        # return cleanly so amp_engine can use CIViC alone.
        out["error"] = "no_token"
        return out

    try:
        import requests
    except ImportError as e:
        out["error"] = f"requests not installed: {e}"
        return out

    url = f"{api_url}/annotate/mutations/byProteinChange"
    params = {
        "hugoSymbol": gene,
        "alteration": protein_change.replace("p.", "") if protein_change else "",
    }
    headers = {"Authorization": f"Bearer {api_token}"}

    t0 = time.time()
    try:
        resp = requests.get(url, params=params, headers=headers,
                             timeout=timeout_s)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        out["error"] = str(e)
        out["latency_ms"] = int((time.time() - t0) * 1000)
        return out

    out["latency_ms"] = int((time.time() - t0) * 1000)
    if not data:
        return out

    out["found"] = True
    out["oncogenic"] = data.get("oncogenic")  # Oncogenic | Likely Oncogenic | Likely Neutral | Unknown
    level = data.get("highestSensitiveLevel") or data.get("highestResistanceLevel")
    if level:
        out["highest_level"] = level
        out["amp_tier"] = ONCOKB_LEVEL_TO_TIER.get(level)

    drugs = set()
    for tx in data.get("treatments", []) or []:
        for d in tx.get("drugs") or []:
            n = d.get("drugName") or d.get("name")
            if n:
                drugs.add(n)
    out["drug_targets"] = sorted(drugs)
    out["evidence_count"] = len(data.get("treatments", []) or [])
    out["source_url"] = (
        f"https://www.oncokb.org/gene/{gene}/{protein_change.replace('p.', '')}"
    )
    return out


def is_enabled(databases_config: dict) -> bool:
    oncokb = (databases_config or {}).get("oncokb", {}) or {}
    return bool(oncokb.get("enabled", False)) and bool(oncokb.get("api_token"))


def get_token(databases_config: dict) -> str:
    return (databases_config or {}).get("oncokb", {}).get("api_token", "")
