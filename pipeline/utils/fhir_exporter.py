"""FHIR R4 DiagnosticReport exporter.

Produces a `DiagnosticReport` resource with one `Observation` per
Pathogenic / Likely-pathogenic variant (germline) or Tier I / II
variant (somatic). Uses the HL7 FHIR `genetics-VariantGenetics`
profile.

Targets ResearchReport / non-clinical use — `status` is set to
'preliminary' since the V2F pipeline does not produce final clinical
sign-out.

Reference:
  https://www.hl7.org/fhir/R4/diagnosticreport.html
  https://www.hl7.org/fhir/R4/genetics.html
"""

from __future__ import annotations

import json
import logging
import os
import uuid
from datetime import datetime, timezone
from typing import Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _is_reportable(row: pd.Series) -> bool:
    """Which variants make it into the FHIR report? Pathogenic/LP germline,
    AMP Tier I/II somatic, Pathogenic/LP CNV."""
    cls = str(row.get("acmg_classification") or "")
    amp = str(row.get("amp_tier") or "")
    cnv = str(row.get("cnv_classification") or "")
    return (cls in ("Pathogenic", "Likely pathogenic")
            or amp in ("I", "II")
            or cnv in ("Pathogenic", "Likely pathogenic"))


def _variant_observation(row: pd.Series, patient_ref: str,
                          report_ref: str) -> Dict:
    """Build a single FHIR Observation for one reportable variant."""
    obs_id = f"v2f-obs-{uuid.uuid4()}"
    coding_significance = str(row.get("acmg_classification")
                                or f"AMP Tier {row.get('amp_tier')}"
                                or row.get("cnv_classification") or "")

    components = []
    if row.get("gene"):
        components.append({
            "code": {"coding": [{
                "system": "http://loinc.org",
                "code": "48018-6",
                "display": "Gene studied",
            }]},
            "valueCodeableConcept": {"coding": [{
                "system": "http://www.genenames.org",
                "code": str(row["gene"]),
                "display": str(row["gene"]),
            }]},
        })
    if row.get("hgvs_p"):
        components.append({
            "code": {"coding": [{
                "system": "http://loinc.org",
                "code": "48005-3",
                "display": "Amino acid change (HGVS.p)",
            }]},
            "valueString": str(row["hgvs_p"]),
        })
    if row.get("hgvs_c"):
        components.append({
            "code": {"coding": [{
                "system": "http://loinc.org",
                "code": "48004-6",
                "display": "DNA change (HGVS.c)",
            }]},
            "valueString": str(row["hgvs_c"]),
        })
    # Bayesian posterior, if computed
    if row.get("bayesian_posterior_prob") is not None:
        try:
            components.append({
                "code": {"coding": [{
                    "system": "http://v2f-reporter.io/codes",
                    "code": "bayesian-posterior-prob",
                    "display": "Bayesian posterior probability of pathogenicity",
                }]},
                "valueQuantity": {"value": float(row["bayesian_posterior_prob"]),
                                    "unit": "probability"},
            })
        except (TypeError, ValueError):
            pass
    # AMP drug targets, if any
    if row.get("amp_drug_targets"):
        components.append({
            "code": {"coding": [{
                "system": "http://v2f-reporter.io/codes",
                "code": "amp-drug-targets",
                "display": "AMP-tier drug targets",
            }]},
            "valueString": str(row["amp_drug_targets"]),
        })

    return {
        "resourceType": "Observation",
        "id": obs_id,
        "meta": {
            "profile": [
                "http://hl7.org/fhir/StructureDefinition/Variant"
            ],
        },
        "status": "preliminary",
        "category": [{
            "coding": [{
                "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                "code": "laboratory",
            }],
        }],
        "code": {"coding": [{
            "system": "http://loinc.org",
            "code": "69548-6",
            "display": "Genetic variant assessment",
        }]},
        "subject": {"reference": patient_ref},
        "valueCodeableConcept": {
            "coding": [{
                "system": "http://v2f-reporter.io/codes/classification",
                "code": coding_significance.replace(" ", "_"),
                "display": coding_significance,
            }],
            "text": coding_significance,
        },
        "component": components,
        "note": [{
            "text": str(row.get("acmg_criteria") or row.get("amp_evidence")
                         or row.get("cnv_evidence_summary") or "")[:500],
        }] if (row.get("acmg_criteria") or row.get("amp_evidence")
                or row.get("cnv_evidence_summary")) else [],
    }


def export_fhir_diagnostic_report(df: pd.DataFrame, proband_id: str,
                                    output_path: str,
                                    config: Optional[dict] = None,
                                    manifest: Optional[dict] = None) -> str:
    """Build a FHIR R4 DiagnosticReport from a classified-variant dataframe.

    Each Pathogenic/LP germline variant + AMP I/II somatic variant + Pathogenic/LP
    CNV becomes an Observation linked from the report's `result` array.
    """
    config = config or {}
    patient_id = f"v2f-patient-{proband_id}"
    patient_ref = f"Patient/{patient_id}"
    report_id = f"v2f-report-{uuid.uuid4()}"
    report_ref = f"DiagnosticReport/{report_id}"

    reportable = df[df.apply(_is_reportable, axis=1)] if not df.empty else df
    observations = [_variant_observation(row, patient_ref, report_ref)
                     for _, row in reportable.iterrows()]

    diagnostic_report = {
        "resourceType": "DiagnosticReport",
        "id": report_id,
        "meta": {
            "profile": [
                "http://hl7.org/fhir/StructureDefinition/DiagnosticReport"
            ],
        },
        "status": "preliminary",
        "category": [{
            "coding": [{
                "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                "code": "GE",
                "display": "Genetics",
            }],
        }],
        "code": {"coding": [{
            "system": "http://loinc.org",
            "code": "81247-9",
            "display": "Master HL7 genetic variant reporting panel",
        }]},
        "subject": {"reference": patient_ref},
        "issued": _now_iso(),
        "performer": [{"display": "V2F Reporter (research-tool)"}],
        "result": [{"reference": f"Observation/{o['id']}"} for o in observations],
        "conclusion": (
            f"V2F automated interpretation of {len(df)} variants — "
            f"{len(reportable)} reportable. "
            "PRELIMINARY: this report is generated by an automated research "
            "tool and is not a clinical sign-out."
        ),
    }
    # Attach the run manifest as an extension so consumers can audit
    if manifest:
        diagnostic_report["extension"] = [{
            "url": "http://v2f-reporter.io/extensions/run-manifest",
            "valueString": json.dumps(manifest, sort_keys=True, default=str),
        }]

    bundle = {
        "resourceType": "Bundle",
        "id": f"v2f-bundle-{uuid.uuid4()}",
        "type": "collection",
        "entry": [
            {"fullUrl": f"urn:uuid:{patient_id}",
             "resource": {"resourceType": "Patient", "id": patient_id,
                           "active": True,
                           "name": [{"family": proband_id, "use": "anonymous"}]}},
            {"fullUrl": f"urn:uuid:{report_id}",
             "resource": diagnostic_report},
        ] + [
            {"fullUrl": f"urn:uuid:{o['id']}", "resource": o}
            for o in observations
        ],
    }

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(bundle, f, indent=2, sort_keys=False)
    logger.info(f"FHIR DiagnosticReport exported: {output_path} "
                f"({len(observations)} observations)")
    return output_path


def is_enabled(config: dict) -> bool:
    out = (config.get("output", {}) or {})
    return bool(out.get("fhir_export", False))
