"""Bulk CSV export of classified variants.

Produces a flat CSV per sample with all classification + curation +
literature columns. Useful for downstream R / Python / Excel analysis.
"""

from __future__ import annotations

import logging
import os
from typing import Iterable, Optional

import pandas as pd

logger = logging.getLogger(__name__)

# Curated set of columns to export; order is meaningful.
_DEFAULT_COLUMNS = [
    "variant_id", "chrom", "pos", "ref", "alt", "gene", "hgvs_p", "hgvs_c",
    "consequence", "severity",
    "acmg_classification", "acmg_criteria", "acmg_tier",
    "bayesian_posterior_prob", "bayesian_classification",
    "amp_tier", "amp_drug_targets", "amp_evidence_level",
    "cnv_classification", "cnv_score", "svtype", "svlen",
    "clinvar_classification", "clinvar_review_stars",
    "gnomad_af", "gnomad_sv_af",
    "cadd_phred", "revel", "spliceai_max", "alphamissense",
    "phenotype_match_score",
    "n_pubmed_ids", "clingen_ca_id",
    "tumor_vaf", "normal_vaf", "is_paired",
    "genotype", "read_depth", "allele_fraction",
]


def export_csv(df: pd.DataFrame, output_path: str,
                columns: Optional[Iterable[str]] = None) -> str:
    """Write `df` to CSV at `output_path`. Returns the path.

    Only columns that actually exist in the dataframe are written; missing
    columns are skipped silently (different analysis modes produce different
    column families).
    """
    cols = list(columns) if columns else _DEFAULT_COLUMNS
    available = [c for c in cols if c in df.columns]
    extras = [c for c in df.columns if c not in cols and not c.startswith("_")]
    final_cols = available + sorted(extras)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df[final_cols].to_csv(output_path, index=False)
    logger.info(f"CSV exported: {output_path} ({len(df)} rows, {len(final_cols)} cols)")
    return output_path
