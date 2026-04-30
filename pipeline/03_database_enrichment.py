"""Stage 3: Reference Database Enrichment

Input: data/annotated/annotated_variants.parquet
Output: data/enriched/variants_enriched.parquet

Steps:
1. Query ClinVar via BigQuery for clinical significance
2. Query gnomAD via BigQuery for population frequencies
3. Merge enrichment data with annotated variants
"""

import json
import logging
import os
import sys
import time

import pandas as pd
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

logger = logging.getLogger(__name__)


def run(config: dict) -> dict:
    """Execute Stage 3: Database Enrichment."""
    t0 = time.time()
    output_dir = config.get("output", {}).get("output_dir", "data")
    ann_dir = os.path.join(output_dir, "annotated")
    enrich_dir = os.path.join(output_dir, "enriched")
    os.makedirs(enrich_dir, exist_ok=True)

    logger.info("Stage 3: Database Enrichment")

    # Load annotated variants
    ann_path = os.path.join(ann_dir, "annotated_variants.parquet")
    if not os.path.exists(ann_path):
        raise FileNotFoundError(f"Stage 2 output not found: {ann_path}")

    df = pd.read_parquet(ann_path)
    logger.info(f"Loaded {len(df)} annotated variants")

    db_config = config.get("databases", {})

    # Step 1: ClinVar enrichment
    clinvar_config = db_config.get("clinvar", {})
    if clinvar_config.get("bq_table"):
        df = _enrich_clinvar(df, clinvar_config)
    else:
        logger.info("ClinVar BQ table not configured — skipping")
        df["clinvar_classification"] = None
        df["clinvar_review_stars"] = 0
        df["clinvar_id"] = None

    # Step 2: gnomAD enrichment
    gnomad_config = db_config.get("gnomad", {})
    if gnomad_config.get("bq_table"):
        df = _enrich_gnomad(df, gnomad_config)
    else:
        logger.info("gnomAD BQ table not configured — skipping")
        df["gnomad_af"] = None
        df["gnomad_af_popmax"] = None

    # Use VEP gnomAD AF as fallback
    if "gnomad_af" in df.columns and "gnomad_af_vep" in df.columns:
        df["gnomad_af"] = df["gnomad_af"].fillna(df["gnomad_af_vep"])

    # Save enriched variants
    output_path = os.path.join(enrich_dir, "variants_enriched.parquet")
    df.to_parquet(output_path, index=False)

    # Stats
    n_clinvar = df["clinvar_classification"].notna().sum() if "clinvar_classification" in df.columns else 0
    n_gnomad = df["gnomad_af"].notna().sum() if "gnomad_af" in df.columns else 0

    result = {
        "stage": "03_database_enrichment",
        "total_variants": len(df),
        "clinvar_matches": int(n_clinvar),
        "gnomad_matches": int(n_gnomad),
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    with open(os.path.join(enrich_dir, "enrichment_summary.json"), "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"\n{'='*60}")
    print(f"Stage 3 Complete: Database Enrichment")
    print(f"{'='*60}")
    print(f"  Total variants:  {len(df):,}")
    print(f"  ClinVar matches: {n_clinvar:,} ({n_clinvar/len(df)*100:.1f}%)")
    print(f"  gnomAD matches:  {n_gnomad:,} ({n_gnomad/len(df)*100:.1f}%)")
    print(f"  Time:            {result['elapsed_seconds']}s")
    print(f"{'='*60}\n")

    return result


def _enrich_clinvar(df: pd.DataFrame, clinvar_config: dict) -> pd.DataFrame:
    """Enrich variants with ClinVar clinical significance."""
    from google.cloud import bigquery

    bq_table = clinvar_config["bq_table"]
    min_stars = clinvar_config.get("min_review_stars", 1)

    logger.info(f"Querying ClinVar: {bq_table}")
    client = bigquery.Client()

    # Build variant list for query
    # ClinVar uses chrom, pos, ref, alt
    variants_for_query = df[["chrom", "pos", "ref", "alt"]].drop_duplicates()
    n_variants = len(variants_for_query)
    logger.info(f"Querying {n_variants} unique variants against ClinVar")

    # For large variant sets, query in batches
    batch_size = 1000
    clinvar_results = []

    for start in range(0, n_variants, batch_size):
        batch = variants_for_query.iloc[start:start + batch_size]

        # Build WHERE clause
        conditions = []
        for _, row in batch.iterrows():
            chrom = str(row["chrom"]).replace("chr", "")
            conditions.append(
                f"(chromosome = '{chrom}' AND start_position = {row['pos']} "
                f"AND reference_bases = '{row['ref']}' AND alternate_bases = '{row['alt']}')"
            )

        if not conditions:
            continue

        where = " OR ".join(conditions)
        query = f"""
        SELECT
            chromosome,
            start_position AS pos,
            reference_bases AS ref,
            alternate_bases AS alt,
            clinical_significance AS classification,
            review_status,
            variation_id AS clinvar_id
        FROM `{bq_table}`
        WHERE {where}
        """

        try:
            batch_df = client.query(query).to_dataframe()
            clinvar_results.append(batch_df)
        except Exception as e:
            logger.warning(f"ClinVar query batch failed: {e}")
            continue

    if not clinvar_results:
        logger.info("No ClinVar matches found")
        df["clinvar_classification"] = None
        df["clinvar_review_stars"] = 0
        df["clinvar_id"] = None
        return df

    clinvar_df = pd.concat(clinvar_results, ignore_index=True)

    # Map review status to stars
    from pipeline.utils._clinvar_stars import map_review_status
    clinvar_df["clinvar_review_stars"] = clinvar_df["review_status"].apply(
        map_review_status
    )

    # Normalize chromosome format
    clinvar_df["chrom"] = clinvar_df["chromosome"].apply(
        lambda x: f"chr{x}" if not str(x).startswith("chr") else str(x)
    )

    # Merge
    clinvar_df = clinvar_df.rename(columns={
        "classification": "clinvar_classification",
        "clinvar_id": "clinvar_id",
    })
    clinvar_df = clinvar_df[["chrom", "pos", "ref", "alt",
                             "clinvar_classification", "clinvar_review_stars",
                             "clinvar_id"]].drop_duplicates(
        subset=["chrom", "pos", "ref", "alt"], keep="first"
    )

    # Ensure pos types match
    df["pos"] = df["pos"].astype(int)
    clinvar_df["pos"] = clinvar_df["pos"].astype(int)

    df = df.merge(clinvar_df, on=["chrom", "pos", "ref", "alt"], how="left")
    n_matches = df["clinvar_classification"].notna().sum()
    logger.info(f"ClinVar: {n_matches} variants matched")

    return df


def _enrich_gnomad(df: pd.DataFrame, gnomad_config: dict) -> pd.DataFrame:
    """Enrich variants with gnomAD population frequencies."""
    from google.cloud import bigquery

    bq_table = gnomad_config["bq_table"]
    fallback = gnomad_config.get("bq_table_fallback", "")

    logger.info(f"Querying gnomAD: {bq_table}")
    client = bigquery.Client()

    variants_for_query = df[["chrom", "pos", "ref", "alt"]].drop_duplicates()
    n_variants = len(variants_for_query)

    batch_size = 1000
    gnomad_results = []

    for start in range(0, n_variants, batch_size):
        batch = variants_for_query.iloc[start:start + batch_size]

        conditions = []
        for _, row in batch.iterrows():
            chrom = str(row["chrom"]).replace("chr", "")
            conditions.append(
                f"(reference_name = '{chrom}' AND start_position = {row['pos']} "
                f"AND reference_bases = '{row['ref']}' AND alternate_bases = '{row['alt']}')"
            )

        if not conditions:
            continue

        where = " OR ".join(conditions)
        query = f"""
        SELECT
            reference_name AS chromosome,
            start_position AS pos,
            reference_bases AS ref,
            alternate_bases AS alt,
            AF AS gnomad_af,
            AF_popmax AS gnomad_af_popmax,
            AN AS gnomad_an,
            AC AS gnomad_ac
        FROM `{bq_table}`
        WHERE {where}
        """

        try:
            batch_df = client.query(query).to_dataframe()
            gnomad_results.append(batch_df)
        except Exception as e:
            logger.warning(f"gnomAD query failed: {e}")
            if fallback:
                logger.info(f"Trying fallback: {fallback}")
                try:
                    query = query.replace(bq_table, fallback)
                    batch_df = client.query(query).to_dataframe()
                    gnomad_results.append(batch_df)
                except Exception as e2:
                    logger.warning(f"gnomAD fallback also failed: {e2}")
            continue

    if not gnomad_results:
        logger.info("No gnomAD matches found")
        df["gnomad_af"] = None
        df["gnomad_af_popmax"] = None
        return df

    gnomad_df = pd.concat(gnomad_results, ignore_index=True)
    gnomad_df["chrom"] = gnomad_df["chromosome"].apply(
        lambda x: f"chr{x}" if not str(x).startswith("chr") else str(x)
    )
    gnomad_df = gnomad_df[["chrom", "pos", "ref", "alt",
                           "gnomad_af", "gnomad_af_popmax"]].drop_duplicates(
        subset=["chrom", "pos", "ref", "alt"], keep="first"
    )

    df["pos"] = df["pos"].astype(int)
    gnomad_df["pos"] = gnomad_df["pos"].astype(int)

    df = df.merge(gnomad_df, on=["chrom", "pos", "ref", "alt"], how="left")
    n_matches = df["gnomad_af"].notna().sum()
    logger.info(f"gnomAD: {n_matches} variants matched")

    return df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    run(config)
