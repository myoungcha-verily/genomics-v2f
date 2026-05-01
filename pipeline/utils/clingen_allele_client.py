"""ClinGen Allele Registry client (http://reg.clinicalgenome.org/).

The Allele Registry assigns canonical CA-IDs to every variant, links to
ClinVar / dbSNP / OMIM, and surfaces *current* pathogenicity assertions —
catching variants that were re-classified after the V2F pipeline's bundled
ClinVar snapshot.

Free, no auth required.
"""

from __future__ import annotations

import logging
import time
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

DEFAULT_API = "http://reg.clinicalgenome.org/allele"


def query_clingen_allele(chrom: str, pos: int, ref: str, alt: str,
                          reference_genome: str = "GRCh38",
                          api_url: str = DEFAULT_API,
                          timeout_s: float = 8.0) -> Dict:
    """Resolve a variant to a ClinGen canonical allele identifier and
    return any current pathogenicity assertions.

    Returns:
        {
          "found": bool,
          "ca_id": str,                   # ClinGen Canonical Allele ID
          "current_clinvar_class": str,   # most recent ClinVar classification
          "n_assertions": int,
          "assertions": list[dict],       # truncated summary
          "source_url": str,
          "error": Optional[str],
          "latency_ms": int,
        }
    """
    out = {
        "found": False, "ca_id": None,
        "current_clinvar_class": None,
        "n_assertions": 0, "assertions": [],
        "source_url": None, "error": None, "latency_ms": 0,
    }

    if not chrom or pos is None or not ref or not alt:
        out["error"] = "chrom, pos, ref, alt required"
        return out

    try:
        import requests
    except ImportError as e:
        out["error"] = f"requests not installed: {e}"
        return out

    # ClinGen Allele Registry HGVS lookup format:
    #   NC_000017.11:g.43106478T>G   (GRCh38 chr17 example)
    refseq_chrom = _refseq_for_chrom(chrom, reference_genome)
    if not refseq_chrom:
        out["error"] = f"unsupported chrom for {reference_genome}: {chrom}"
        return out

    hgvs = f"{refseq_chrom}:g.{pos}{ref}>{alt}"
    url = f"{api_url}?hgvs={hgvs}"

    t0 = time.time()
    try:
        r = requests.get(url, timeout=timeout_s,
                          headers={"Accept": "application/json"})
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        out["error"] = str(e)
        out["latency_ms"] = int((time.time() - t0) * 1000)
        return out

    out["latency_ms"] = int((time.time() - t0) * 1000)
    if not data or not isinstance(data, dict):
        return out

    # ClinGen response shape:
    #   {"@id": ".../CA000123", "communityStandardTitle": [...],
    #    "externalRecords": {"ClinVarVariations": [...]}, ...}
    full_id = data.get("@id", "")
    if "/" in full_id:
        out["ca_id"] = full_id.rsplit("/", 1)[-1]
    out["found"] = bool(out["ca_id"])

    # Pull ClinVar assertions if present
    cv_recs = (data.get("externalRecords", {}) or {}).get("ClinVarVariations") or []
    assertions = []
    for cv in cv_recs:
        clas = cv.get("RCV", [{}])
        if isinstance(clas, list) and clas:
            clas = clas[0]
        sig = (clas or {}).get("clinicalSignificance") or cv.get("clinicalSignificance")
        if sig:
            assertions.append({
                "source": "ClinVar",
                "rcv": cv.get("RCV"),
                "significance": sig,
                "review_status": cv.get("reviewStatus"),
            })
    out["n_assertions"] = len(assertions)
    out["assertions"] = assertions[:5]
    if assertions:
        out["current_clinvar_class"] = assertions[0]["significance"]

    if out["ca_id"]:
        out["source_url"] = f"http://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid={out['ca_id']}"

    return out


# Minimal RefSeq map for the common chromosomes we care about (GRCh38).
# A full deployment would source this from a checked-in reference.
_REFSEQ_GRCH38 = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
}
_REFSEQ_GRCH37 = {
    "1": "NC_000001.10", "2": "NC_000002.11", "3": "NC_000003.11",
    "4": "NC_000004.11", "5": "NC_000005.9", "6": "NC_000006.11",
    "7": "NC_000007.13", "17": "NC_000017.10", "13": "NC_000013.10",
    "10": "NC_000010.10", "11": "NC_000011.9", "14": "NC_000014.8",
    "X": "NC_000023.10", "Y": "NC_000024.9",
}


def _refseq_for_chrom(chrom: str, reference_genome: str) -> Optional[str]:
    bare = str(chrom).replace("chr", "")
    if reference_genome == "GRCh37" or reference_genome == "hg19":
        return _REFSEQ_GRCH37.get(bare)
    return _REFSEQ_GRCH38.get(bare)


def is_enabled(config: dict) -> bool:
    db = (config.get("databases", {}) or {})
    cg = db.get("clingen_allele", {}) or {}
    return bool(cg.get("enabled", True))  # default-on (free API)
