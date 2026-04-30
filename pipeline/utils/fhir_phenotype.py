"""FHIR phenotype extraction for variant-phenotype correlation.

Queries FHIR Condition and Observation resources from BigQuery
to build a patient phenotype profile for variant prioritization.
"""

import logging
import re
from typing import Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)

# ICD-10 to OMIM disease mapping (subset of common genetic conditions)
ICD10_TO_OMIM = {
    # Cardiac
    "I42.1": {"omim": "115200", "disease": "Hypertrophic cardiomyopathy", "genes": ["MYH7", "MYBPC3", "TNNT2", "TNNI3"]},
    "I42.0": {"omim": "115200", "disease": "Dilated cardiomyopathy", "genes": ["TTN", "LMNA", "MYH7", "TNNT2"]},
    "I45.81": {"omim": "192500", "disease": "Long QT syndrome", "genes": ["KCNQ1", "KCNH2", "SCN5A"]},
    "I49.8": {"omim": "601144", "disease": "Brugada syndrome", "genes": ["SCN5A"]},
    # Cancer predisposition
    "Z15.01": {"omim": "604370", "disease": "BRCA1-related cancer", "genes": ["BRCA1"]},
    "Z15.02": {"omim": "612555", "disease": "BRCA2-related cancer", "genes": ["BRCA2"]},
    "Z80.0": {"omim": "120435", "disease": "Lynch syndrome", "genes": ["MLH1", "MSH2", "MSH6", "PMS2"]},
    # Neuro
    "G10": {"omim": "143100", "disease": "Huntington disease", "genes": ["HTT"]},
    "G12.21": {"omim": "253300", "disease": "Spinal muscular atrophy", "genes": ["SMN1"]},
    "G71.0": {"omim": "310200", "disease": "Duchenne muscular dystrophy", "genes": ["DMD"]},
    # Metabolic
    "E84": {"omim": "219700", "disease": "Cystic fibrosis", "genes": ["CFTR"]},
    "E80.0": {"omim": "176000", "disease": "Acute intermittent porphyria", "genes": ["HMBS"]},
    # Hematologic
    "D57": {"omim": "603903", "disease": "Sickle cell disease", "genes": ["HBB"]},
    "D56": {"omim": "141900", "disease": "Thalassemia", "genes": ["HBA1", "HBA2", "HBB"]},
    # Connective tissue
    "Q87.40": {"omim": "154700", "disease": "Marfan syndrome", "genes": ["FBN1"]},
    "Q79.6": {"omim": "130050", "disease": "Ehlers-Danlos syndrome", "genes": ["COL5A1", "COL5A2"]},
}

# HPO term mapping for common phenotypic features
PHENOTYPE_TO_HPO = {
    "cardiomyopathy": "HP:0001638",
    "arrhythmia": "HP:0011675",
    "seizure": "HP:0001250",
    "intellectual disability": "HP:0001249",
    "hearing loss": "HP:0000365",
    "vision loss": "HP:0000572",
    "muscle weakness": "HP:0001324",
    "short stature": "HP:0004322",
    "renal failure": "HP:0000083",
}


def extract_patient_phenotype(patient_id: str, config: dict) -> Dict:
    """Extract phenotype profile from FHIR data for a patient.

    Returns dict with:
        patient_id: str
        conditions: list of ICD-10 codes with descriptions
        observations: list of relevant clinical observations
        omim_diseases: list of matched OMIM diseases
        candidate_genes: list of genes associated with phenotype
        hpo_terms: list of HPO terms (if mappable)
        status: 'disabled' | 'unconfigured' | 'failed' | 'ok' — distinguishes
                quiet skip from real failure so the dashboard can surface it.
    """
    pheno_config = config.get("phenotype", {})
    if not pheno_config.get("enabled", False):
        logger.info("Phenotype matching disabled in config (phenotype.enabled=false)")
        return _empty_phenotype(patient_id, status="disabled")

    fhir_project = pheno_config.get("fhir_project", "")
    fhir_dataset = pheno_config.get("fhir_dataset", "")

    if not fhir_project or not fhir_dataset:
        logger.warning(
            "Phenotype enabled but FHIR project/dataset not configured — "
            "set phenotype.fhir_project and phenotype.fhir_dataset to enable"
        )
        return _empty_phenotype(patient_id, status="unconfigured")

    try:
        from google.cloud import bigquery
        client = bigquery.Client()
    except Exception as e:
        logger.error(f"BigQuery client init failed (phenotype broken): {e}")
        return _empty_phenotype(patient_id, status="failed", error=str(e))

    # Query conditions
    conditions = _query_conditions(client, patient_id, fhir_project,
                                    fhir_dataset, pheno_config)

    # Query observations
    observations = _query_observations(client, patient_id, fhir_project,
                                        fhir_dataset, pheno_config)

    # Map to OMIM diseases and candidate genes
    omim_diseases = []
    candidate_genes = set()

    for cond in conditions:
        icd_code = cond.get("code", "")
        # Try exact match first, then prefix
        omim_match = ICD10_TO_OMIM.get(icd_code)
        if not omim_match:
            # Try prefix (e.g., I42 for I42.1)
            prefix = icd_code.split(".")[0]
            omim_match = ICD10_TO_OMIM.get(prefix)

        if omim_match:
            omim_diseases.append(omim_match)
            candidate_genes.update(omim_match.get("genes", []))

    return {
        "patient_id": patient_id,
        "conditions": conditions,
        "observations": observations,
        "omim_diseases": omim_diseases,
        "candidate_genes": sorted(candidate_genes),
        "hpo_terms": [],  # Would need full HPO mapping
        "n_conditions": len(conditions),
        "n_candidate_genes": len(candidate_genes),
        "status": "ok",
    }


def _query_conditions(client, patient_id: str, project: str,
                       dataset: str, pheno_config: dict) -> List[Dict]:
    """Query FHIR Condition table for patient diagnoses."""
    table = pheno_config.get("condition_table", "Condition")
    query = f"""
    SELECT
        code.coding[SAFE_OFFSET(0)].code AS code,
        code.coding[SAFE_OFFSET(0)].display AS display,
        code.coding[SAFE_OFFSET(0)].system AS system,
        onsetDateTime,
        clinicalStatus.coding[SAFE_OFFSET(0)].code AS clinical_status
    FROM `{project}.{dataset}.{table}`
    WHERE subject.reference LIKE '%{patient_id}%'
    ORDER BY onsetDateTime DESC
    LIMIT 100
    """
    try:
        df = client.query(query).to_dataframe()
        return df.to_dict("records")
    except Exception as e:
        logger.warning(f"Condition query failed: {e}")
        return []


def _query_observations(client, patient_id: str, project: str,
                         dataset: str, pheno_config: dict) -> List[Dict]:
    """Query FHIR Observation table for relevant clinical observations."""
    table = pheno_config.get("observation_table", "Observation")
    query = f"""
    SELECT
        code.coding[SAFE_OFFSET(0)].code AS code,
        code.coding[SAFE_OFFSET(0)].display AS display,
        valueQuantity.value AS value,
        valueQuantity.unit AS unit,
        effectiveDateTime
    FROM `{project}.{dataset}.{table}`
    WHERE subject.reference LIKE '%{patient_id}%'
    ORDER BY effectiveDateTime DESC
    LIMIT 50
    """
    try:
        df = client.query(query).to_dataframe()
        return df.to_dict("records")
    except Exception as e:
        logger.warning(f"Observation query failed: {e}")
        return []


def score_variant_phenotype_match(variant: dict,
                                   phenotype: dict) -> float:
    """Score how well a variant matches the patient phenotype.

    Returns score 0.0-1.0 where:
    - 1.0: Variant gene directly matches a phenotype-associated gene
    - 0.5: Variant gene is in the same pathway
    - 0.0: No phenotype match

    This is the key differentiator — using FHIR clinical data
    to prioritize variants.
    """
    gene = variant.get("gene", "")
    if not gene:
        return 0.0

    candidate_genes = set(phenotype.get("candidate_genes", []))
    if not candidate_genes:
        return 0.0

    # Direct match
    if gene in candidate_genes:
        return 1.0

    # Could add pathway-level matching here
    return 0.0


def _empty_phenotype(patient_id: str, status: str = "disabled",
                      error: Optional[str] = None) -> Dict:
    """Return empty phenotype profile.

    `status` distinguishes:
      - 'disabled' — phenotype.enabled=false in config (silent skip)
      - 'unconfigured' — enabled but project/dataset missing
      - 'failed' — enabled and configured but query/auth broke
    """
    return {
        "patient_id": patient_id,
        "conditions": [],
        "observations": [],
        "omim_diseases": [],
        "candidate_genes": [],
        "hpo_terms": [],
        "n_conditions": 0,
        "n_candidate_genes": 0,
        "status": status,
        "error": error,
    }


def test_fhir_connectivity(fhir_project: str, fhir_dataset: str,
                            condition_table: str = "Condition") -> Dict:
    """One-shot connectivity probe: SELECT 1 FROM <project>.<dataset>.Condition.

    Returns {ok: bool, error: Optional[str], latency_ms: int, project, dataset, table}.
    Used by the dashboard to surface configuration problems early instead of
    letting stage 5 silently report 0 conditions.
    """
    import time
    out = {
        "project": fhir_project,
        "dataset": fhir_dataset,
        "table": condition_table,
        "ok": False,
        "error": None,
        "latency_ms": 0,
    }
    if not fhir_project or not fhir_dataset:
        out["error"] = "FHIR project and dataset must both be set"
        return out
    try:
        from google.cloud import bigquery
    except ImportError as e:
        out["error"] = f"google-cloud-bigquery not installed: {e}"
        return out

    t0 = time.time()
    try:
        client = bigquery.Client()
        query = f"SELECT 1 AS ok FROM `{fhir_project}.{fhir_dataset}.{condition_table}` LIMIT 1"
        # Use dry run first to avoid running the actual query if syntax/access fails fast
        job_config = bigquery.QueryJobConfig(dry_run=True, use_query_cache=False)
        client.query(query, job_config=job_config).result()
        out["ok"] = True
    except Exception as e:
        out["error"] = str(e)
    finally:
        out["latency_ms"] = int((time.time() - t0) * 1000)
    return out
