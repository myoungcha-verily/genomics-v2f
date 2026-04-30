"""Jinja2 HTML report rendering engine."""

import logging
import os
from datetime import datetime
from typing import Dict, List, Optional

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from pipeline.utils.acmg_engine import classification_color, classification_tier
from pipeline.utils.hgvs_formatter import variant_display_name, format_frequency
from pipeline.utils.population_freq import format_frequency as fmt_freq

logger = logging.getLogger(__name__)


def render_proband_report(
    variants_df: pd.DataFrame,
    proband_id: str,
    config: dict,
    phenotype: Optional[Dict] = None,
    qc_metrics: Optional[Dict] = None,
    output_path: Optional[str] = None,
) -> str:
    """Render a clinical variant report for a proband.

    Args:
        variants_df: Classified variants DataFrame
        proband_id: Sample ID
        config: Pipeline config
        phenotype: Patient phenotype profile (optional)
        qc_metrics: QC metrics from Stage 1 (optional)
        output_path: Where to save the report

    Returns:
        Path to rendered HTML report
    """
    template_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "templates"
    )

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(["html"]),
    )

    # Add custom filters
    env.filters["format_af"] = fmt_freq
    env.filters["classification_color"] = classification_color

    template = env.get_template("proband_report.html")

    # Prepare variant data by tier
    tier1 = _filter_tier(variants_df, ["Pathogenic", "Likely pathogenic"])
    tier2 = _filter_tier(variants_df, ["VUS"])
    tier3 = _filter_tier(variants_df, ["Likely benign", "Benign"])

    # Limit VUS in report
    max_vus = config.get("output", {}).get("max_vus_in_report", 50)
    if len(tier2) > max_vus:
        tier2 = tier2.head(max_vus)

    # Include benign?
    include_benign = config.get("output", {}).get("include_benign", False)
    if not include_benign:
        tier3 = pd.DataFrame()

    # Prepare template context
    context = {
        "proband_id": proband_id,
        "report_date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "pipeline_version": "1.0.0",
        "reference_genome": config.get("input", {}).get("reference_genome", "GRCh38"),
        "adapter_type": config.get("input", {}).get("adapter", "single_sample"),
        # Variant counts
        "total_variants": len(variants_df),
        "n_pathogenic": len(variants_df[variants_df["acmg_classification"] == "Pathogenic"])
            if "acmg_classification" in variants_df.columns else 0,
        "n_likely_pathogenic": len(variants_df[variants_df["acmg_classification"] == "Likely pathogenic"])
            if "acmg_classification" in variants_df.columns else 0,
        "n_vus": len(variants_df[variants_df["acmg_classification"] == "VUS"])
            if "acmg_classification" in variants_df.columns else 0,
        "n_likely_benign": len(variants_df[variants_df["acmg_classification"] == "Likely benign"])
            if "acmg_classification" in variants_df.columns else 0,
        "n_benign": len(variants_df[variants_df["acmg_classification"] == "Benign"])
            if "acmg_classification" in variants_df.columns else 0,
        # Variant tiers
        "tier1_variants": _variants_to_dicts(tier1),
        "tier2_variants": _variants_to_dicts(tier2),
        "tier3_variants": _variants_to_dicts(tier3),
        "max_vus_shown": max_vus,
        "include_benign": include_benign,
        # Phenotype
        "phenotype": phenotype,
        "has_phenotype": phenotype is not None and phenotype.get("n_conditions", 0) > 0,
        # QC
        "qc_metrics": qc_metrics,
    }

    html = template.render(**context)

    if output_path is None:
        reports_dir = config.get("output", {}).get("reports_dir", "reports")
        os.makedirs(reports_dir, exist_ok=True)
        output_path = os.path.join(reports_dir, f"{proband_id}_report.html")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)

    logger.info(f"Report saved: {output_path}")
    return output_path


def _filter_tier(df: pd.DataFrame, classifications: List[str]) -> pd.DataFrame:
    """Filter variants by ACMG classification."""
    if "acmg_classification" not in df.columns:
        return pd.DataFrame()
    return df[df["acmg_classification"].isin(classifications)].copy()


def _variants_to_dicts(df: pd.DataFrame) -> List[Dict]:
    """Convert variant DataFrame to list of dicts for template rendering."""
    if df.empty:
        return []

    records = []
    for _, row in df.iterrows():
        record = {
            "variant_id": row.get("variant_id", ""),
            "gene": row.get("gene", ""),
            "display_name": variant_display_name(
                row.get("gene", ""),
                row.get("hgvs_p", ""),
                row.get("hgvs_c", ""),
                row.get("chrom", ""),
                row.get("pos", 0),
                row.get("ref", ""),
                row.get("alt", ""),
            ),
            "chrom": row.get("chrom", ""),
            "pos": row.get("pos", 0),
            "ref": row.get("ref", ""),
            "alt": row.get("alt", ""),
            "genotype": row.get("genotype", ""),
            "consequence": row.get("consequence", ""),
            "severity": row.get("severity", ""),
            "hgvs_c": row.get("hgvs_c", ""),
            "hgvs_p": row.get("hgvs_p", ""),
            "classification": row.get("acmg_classification", "VUS"),
            "classification_color": classification_color(
                row.get("acmg_classification", "VUS")
            ),
            "criteria": row.get("acmg_criteria", ""),
            "gnomad_af": row.get("gnomad_af"),
            "gnomad_af_display": fmt_freq(row.get("gnomad_af")),
            "clinvar": row.get("clinvar_classification", ""),
            "clinvar_stars": row.get("clinvar_review_stars", 0),
            "cadd": row.get("cadd_phred"),
            "revel": row.get("revel"),
            "spliceai": row.get("spliceai_max"),
            "alphamissense": row.get("alphamissense"),
            "read_depth": row.get("read_depth", 0),
            "allele_fraction": row.get("allele_fraction", 0),
            "is_de_novo": row.get("is_de_novo", False),
            "phenotype_match": row.get("phenotype_match_score", 0),
        }
        records.append(record)

    return records
