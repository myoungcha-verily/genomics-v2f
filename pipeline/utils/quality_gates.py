"""Per-stage quality gates.

Acceptance criteria that fire at the end of each pipeline stage. Each gate
returns a `GateResult` describing whether it passed, the actual values, the
expected range, and a severity (warning | hard_fail).

Configured via `pipeline.gates.<stage>.<gate>` in pipeline_config.yaml. By
default gates are warnings only — researchers can opt into hard-fail mode
once they trust the gates for their specific data type.

Example config:
    pipeline:
      gates:
        mode: warning           # warning | hard_fail
        stage_1:
          min_variants: 10
          ti_tv_ratio: [1.5, 3.5]
        stage_3:
          clinvar_match_rate_min: 0.30
        stage_4:
          all_vus_max_fraction: 0.95   # > 95% all-VUS is suspicious
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class GateResult:
    name: str
    passed: bool
    severity: str  # 'warning' | 'hard_fail'
    actual: Any
    expected: Any
    message: str

    def to_dict(self) -> Dict:
        return asdict(self)


class GateFailure(Exception):
    """Raised when a hard_fail gate is breached. Stage halts; pipeline aborts."""
    def __init__(self, gate: GateResult):
        super().__init__(gate.message)
        self.gate = gate


def _gate_config(config: dict, stage: str) -> dict:
    return ((config.get("pipeline", {}) or {}).get("gates", {}) or {}).get(stage, {}) or {}


def _gate_mode(config: dict) -> str:
    return ((config.get("pipeline", {}) or {})
             .get("gates", {}) or {}).get("mode", "warning")


def _result(name: str, passed: bool, actual, expected, message: str,
             config: dict) -> GateResult:
    severity = _gate_mode(config) if not passed else "warning"
    return GateResult(name=name, passed=passed, severity=severity,
                       actual=actual, expected=expected, message=message)


def _enforce(result: GateResult) -> None:
    """Log + optionally raise based on severity."""
    if result.passed:
        logger.debug(f"gate[{result.name}] passed: {result.actual}")
        return
    log_fn = logger.error if result.severity == "hard_fail" else logger.warning
    log_fn(f"gate[{result.name}] {result.severity}: {result.message}")
    if result.severity == "hard_fail":
        raise GateFailure(result)


# ---------- Stage 1: VCF ingest QC ---------------------------------

def check_stage_1(qc_metrics: Dict, config: dict) -> List[GateResult]:
    """Gates that apply after stage 1 (VCF ingest)."""
    gates_cfg = _gate_config(config, "stage_1")
    results = []

    n = int(qc_metrics.get("n_variants", 0))
    min_var = int(gates_cfg.get("min_variants", 1))
    results.append(_result(
        "stage_1.min_variants",
        passed=n >= min_var,
        actual=n, expected=f">= {min_var}",
        message=f"Stage 1: only {n} variants; expected ≥ {min_var}",
        config=config,
    ))

    titv = qc_metrics.get("ti_tv_ratio")
    titv_range = gates_cfg.get("ti_tv_ratio", [1.5, 3.5])
    if titv is not None and titv > 0:
        passed = titv_range[0] <= titv <= titv_range[1]
        results.append(_result(
            "stage_1.ti_tv_ratio",
            passed=passed,
            actual=round(titv, 2), expected=f"in {titv_range}",
            message=(f"Stage 1: Ti/Tv = {titv:.2f}, expected in "
                     f"{titv_range}. Out-of-range Ti/Tv often signals "
                     f"contamination, mis-mapping, or wrong sequencing assay."),
            config=config,
        ))

    return results


# ---------- Stage 3: enrichment ------------------------------------

def check_stage_3(enrichment_summary: Dict, config: dict) -> List[GateResult]:
    """Gates after stage 3 (ClinVar/gnomAD enrichment).

    The most common stage-3 failure mode (which we hit live in M0!) is
    that the configured BQ tables don't exist or use a different schema —
    queries silently return zero matches. The enrichment-rate gate
    catches this fast.
    """
    gates_cfg = _gate_config(config, "stage_3")
    results = []

    n_total = int(enrichment_summary.get("total_variants", 0))
    n_clinvar = int(enrichment_summary.get("clinvar_matches", 0))
    n_gnomad = int(enrichment_summary.get("gnomad_matches", 0))
    rate_clinvar = (n_clinvar / n_total) if n_total else 0
    rate_gnomad = (n_gnomad / n_total) if n_total else 0

    cv_min = float(gates_cfg.get("clinvar_match_rate_min", 0.0))
    if cv_min > 0:
        results.append(_result(
            "stage_3.clinvar_match_rate",
            passed=rate_clinvar >= cv_min,
            actual=f"{rate_clinvar:.1%}", expected=f">= {cv_min:.0%}",
            message=(f"Stage 3: ClinVar match rate {rate_clinvar:.1%} below "
                     f"threshold {cv_min:.0%}. Likely causes: BQ table "
                     f"misconfigured, schema drift, or VPC SC blocking. "
                     f"Run `bq show {enrichment_summary.get('clinvar_table', '?')}` "
                     f"to verify."),
            config=config,
        ))

    g_min = float(gates_cfg.get("gnomad_match_rate_min", 0.0))
    if g_min > 0:
        results.append(_result(
            "stage_3.gnomad_match_rate",
            passed=rate_gnomad >= g_min,
            actual=f"{rate_gnomad:.1%}", expected=f">= {g_min:.0%}",
            message=(f"Stage 3: gnomAD match rate {rate_gnomad:.1%} below "
                     f"threshold {g_min:.0%}."),
            config=config,
        ))

    return results


# ---------- Stage 4: classification distribution -------------------

def check_stage_4(class_counts: Dict[str, int], config: dict) -> List[GateResult]:
    """Gates after stage 4 (ACMG classification).

    Catches the 'pipeline produced all VUS' failure mode that hits when
    upstream enrichment was broken — the engine has no positive evidence,
    nothing fires, and every variant lands in VUS.
    """
    gates_cfg = _gate_config(config, "stage_4")
    results = []

    total = sum(class_counts.values())
    if total == 0:
        return results

    vus = class_counts.get("VUS", 0)
    vus_frac = vus / total
    max_frac = float(gates_cfg.get("all_vus_max_fraction", 0.95))
    results.append(_result(
        "stage_4.vus_fraction",
        passed=vus_frac <= max_frac,
        actual=f"{vus_frac:.1%} ({vus}/{total})",
        expected=f"<= {max_frac:.0%}",
        message=(f"Stage 4: {vus_frac:.1%} of variants ended up VUS "
                 f"(>{max_frac:.0%} threshold). This usually signals "
                 f"upstream enrichment failure — the engine had no positive "
                 f"or negative evidence to fire on."),
        config=config,
    ))

    return results


# ---------- Aggregate runner ---------------------------------------

def evaluate(stage: str, payload: Dict, config: dict) -> List[Dict]:
    """Run all gates for `stage` and return their dict-form results.

    Hard-fail gates raise GateFailure; warning gates only log + return.
    Pipeline stages call this at the end of their work; the result list
    is included in the per-stage summary JSON.
    """
    if stage == "stage_1":
        results = check_stage_1(payload, config)
    elif stage == "stage_3":
        results = check_stage_3(payload, config)
    elif stage == "stage_4":
        results = check_stage_4(payload, config)
    else:
        return []

    for r in results:
        _enforce(r)  # may raise GateFailure for hard_fail
    return [r.to_dict() for r in results]
