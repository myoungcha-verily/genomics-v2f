"""Pipeline run manifest — provenance for every report.

Goal: a researcher (or auditor) running this pipeline a year from now
should be able to identify exactly which code, config, and reference
data produced any given result. Two runs of the same VCF + same config
+ same code should produce byte-identical outputs (modulo the
manifest's own timestamp).

The manifest is written alongside each report and embedded as a footer
in the HTML, a structured field in FHIR output, and the run-history
table.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import subprocess
import time
import uuid
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# Bumped manually when the manifest schema changes
MANIFEST_SCHEMA_VERSION = "1.0"


def _git_sha() -> Optional[str]:
    """Best-effort git SHA. Returns None if not in a git checkout
    (e.g. running from a deployed image without .git/)."""
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=here, stderr=subprocess.DEVNULL
        ).decode().strip()
        return out
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        # Fallback: env var (set in Dockerfile via LABEL → docker build arg)
        return os.environ.get("V2F_GIT_SHA")


def _git_dirty() -> Optional[bool]:
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=here, stderr=subprocess.DEVNULL
        ).decode().strip()
        return bool(out)
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None


def _config_hash(config: dict) -> str:
    canonical = json.dumps(config, sort_keys=True, default=str)
    return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


def _ref_data_versions(config: dict) -> Dict[str, Optional[str]]:
    """Snapshot of the reference data version each external table is at.

    Best-effort — returns table names + last_modified when BQ is reachable;
    bundled file paths + mtime otherwise. Never raises.
    """
    versions: Dict[str, Optional[str]] = {}

    db = (config.get("databases", {}) or {})
    versions["clinvar_table"] = (db.get("clinvar", {}) or {}).get("bq_table")
    versions["gnomad_table_pattern"] = (db.get("gnomad", {}) or {}).get("bq_table")
    versions["gnomad_sv_table"] = (db.get("gnomad_sv", {}) or {}).get("bq_table")
    versions["civic_api"] = (db.get("civic", {}) or {}).get("api_url")
    versions["oncokb_enabled"] = bool((db.get("oncokb", {}) or {}).get("api_token"))
    versions["mavedb_table"] = (db.get("mavedb", {}) or {}).get("bq_table")

    # Bundled reference files (their mtime is a proxy for their version)
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    for rel in ("reference/insilico_calibration.json",
                "reference/clingen_dosage_sensitivity.json",
                "reference/acmg_criteria_spec.json"):
        p = os.path.join(here, rel)
        if os.path.exists(p):
            versions[os.path.basename(rel)] = (
                f"mtime={int(os.path.getmtime(p))}")

    return versions


def build_manifest(config: dict, run_id: Optional[str] = None) -> Dict:
    """Build the per-run manifest dict. Cheap; no external network calls."""
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "run_id": run_id or str(uuid.uuid4()),
        "started_at": int(time.time()),
        "git_sha": _git_sha(),
        "git_dirty": _git_dirty(),
        "config_hash": _config_hash(config),
        "analysis_mode": (config.get("input", {}) or {}).get("analysis_mode", "germline"),
        "adapter": (config.get("input", {}) or {}).get("adapter", ""),
        "reference_genome": (config.get("input", {}) or {}).get("reference_genome", ""),
        "reference_data": _ref_data_versions(config),
    }


def write_manifest(config: dict, output_dir: str,
                    run_id: Optional[str] = None) -> Dict:
    """Write the manifest as JSON to `<output_dir>/run_manifest.json` and
    return it. Idempotent; safe to call multiple times in one pipeline run."""
    os.makedirs(output_dir, exist_ok=True)
    manifest = build_manifest(config, run_id=run_id)
    out_path = os.path.join(output_dir, "run_manifest.json")
    with open(out_path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True, default=str)
    return manifest


def manifest_footer_html(manifest: Dict) -> str:
    """Render the run manifest as an HTML footer for inclusion in reports."""
    sha = (manifest.get("git_sha") or "unknown")[:12]
    dirty = "+" if manifest.get("git_dirty") else ""
    cfg = (manifest.get("config_hash") or "")[7:19]  # strip 'sha256:' prefix
    rid = manifest.get("run_id", "?")[:8]
    mode = manifest.get("analysis_mode", "?")
    return (
        '<div class="run-manifest" '
        'style="margin-top:32px;padding:10px 14px;border-top:1px solid #ccc;'
        'font-size:11px;color:#888;font-family:monospace;">'
        f'pipeline={sha}{dirty} · config={cfg} · '
        f'mode={mode} · run={rid}'
        '</div>'
    )
