"""Gene-specific VCEP rule loader.

ClinGen Variant Curation Expert Panels (VCEPs) publish gene-specific
ACMG rule overrides — e.g. MYH7 disabling PVS1, PTEN tightening PM2,
TP53 promoting hotspot PM5 to strong. Generic ACMG underperforms VCEP
rules on covered genes.

This loader reads `reference/vcep_rules/<GENE>.json` files and merges
the overrides into the per-variant config that's passed to each
`eval_*` function.

Override schema (per file):
{
  "gene": "MYH7",
  "overrides": {
    "PVS1": {"strict_lof": true, "exclude_last_exon": true},
    "PM2": {"threshold": 0.00004, "strength": "supporting"},
    "BS1": {"threshold": 0.00004}
  }
}

The merged config keeps the existing acmg.* keys and adds
`acmg.vcep_overrides_applied: <gene>` (informational) plus per-criterion
strength/threshold tweaks under `acmg.<criterion>_override`. Existing
eval_* functions read overrides via `get_override(criterion, config)`.
"""

import json
import logging
import os
from typing import Dict, Optional

logger = logging.getLogger(__name__)

_VCEP_CACHE: Dict[str, dict] = {}


def is_enabled(config: dict) -> bool:
    return bool(((config.get("acmg", {}) or {}).get("vcep", {}) or {})
                .get("enabled", True))  # default-on


def _vcep_dir(config: dict) -> str:
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    rel = ((config.get("acmg", {}) or {}).get("vcep", {}) or {}).get(
        "panels_dir", "reference/vcep_rules")
    return os.path.join(here, rel)


def load_vcep_rules(gene: str, config: dict) -> Optional[dict]:
    """Load the VCEP rule file for a gene, or return None if no file exists."""
    if not gene:
        return None
    if gene in _VCEP_CACHE:
        return _VCEP_CACHE[gene]
    path = os.path.join(_vcep_dir(config), f"{gene}.json")
    if not os.path.exists(path):
        _VCEP_CACHE[gene] = None
        return None
    try:
        with open(path) as f:
            data = json.load(f)
        _VCEP_CACHE[gene] = data
        return data
    except Exception as e:
        logger.warning(f"VCEP rule load failed for {gene}: {e}")
        _VCEP_CACHE[gene] = None
        return None


def merge_vcep_overrides(variant: dict, config: dict) -> dict:
    """Return a per-variant config with VCEP overrides merged in.

    Caller can use this as the `config` arg to eval_* functions in place
    of the global config. Non-destructive — original config is unchanged.
    """
    if not is_enabled(config):
        return config

    gene = variant.get("gene", "")
    rules = load_vcep_rules(gene, config)
    if not rules:
        return config

    # Shallow copy then deep-update the acmg section
    merged = dict(config)
    merged["acmg"] = dict(config.get("acmg", {}) or {})
    overrides = rules.get("overrides", {}) or {}

    # PM2 / BS1 / BA1 thresholds: patch the corresponding acmg.* keys
    threshold_map = {
        "PM2": "pm2_threshold",
        "BS1": "bs1_threshold",
        "BA1": "ba1_threshold",
    }
    for criterion, ov in overrides.items():
        if criterion in threshold_map and "threshold" in ov:
            merged["acmg"][threshold_map[criterion]] = ov["threshold"]

    # Per-criterion overrides under acmg.vcep_overrides for eval_* functions
    merged["acmg"]["vcep_overrides"] = overrides
    merged["acmg"]["vcep_gene"] = gene
    merged["acmg"]["vcep_id"] = rules.get("vcep_id", "")
    return merged


def get_override(criterion: str, config: dict) -> Optional[dict]:
    """Retrieve VCEP override for a specific criterion (None if absent)."""
    return ((config.get("acmg", {}) or {}).get("vcep_overrides", {}) or {}).get(criterion)


def list_supported_genes(config: dict) -> list:
    """List genes that have a VCEP rule file."""
    d = _vcep_dir(config)
    if not os.path.isdir(d):
        return []
    return sorted(os.path.splitext(f)[0] for f in os.listdir(d) if f.endswith(".json"))
