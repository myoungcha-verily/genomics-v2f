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

    # Step 2: gnomAD enrichment (small variants)
    gnomad_config = db_config.get("gnomad", {})
    sv_mask = df["svtype"].notna() if "svtype" in df.columns else None
    if sv_mask is not None and sv_mask.any():
        small_df = df[~sv_mask].copy()
    else:
        small_df = df

    if not small_df.empty and (gnomad_config.get("bq_table_pattern")
                                 or gnomad_config.get("bq_table")):
        small_df = _enrich_gnomad(small_df, gnomad_config)
    else:
        logger.info("gnomAD BQ table not configured / no small variants — skipping")
        small_df["gnomad_af"] = None
        small_df["gnomad_af_popmax"] = None

    # SV-specific enrichment (gnomAD-SV + gene-content count)
    if sv_mask is not None and sv_mask.any():
        sv_df = df[sv_mask].copy()
        sv_df = _enrich_svs(sv_df, db_config)
        df = pd.concat([small_df, sv_df], ignore_index=True)
    else:
        df = small_df

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
        "clinvar_table": (config.get("databases", {}).get("clinvar", {}) or {}).get("bq_table"),
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    # Run stage 3 acceptance gates
    from pipeline.utils import quality_gates
    result["quality_gates"] = quality_gates.evaluate("stage_3", result, config)

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
    """Enrich variants with ClinVar clinical significance.

    Targets the ncbi_clinvar_hg38_* schema on bigquery-public-data:
      reference_name STRING (no 'chr' prefix)
      start_position INTEGER (0-based)
      reference_bases STRING
      alternate_bases ARRAY<STRUCT<alt STRING>>
      CLNSIG ARRAY<STRING>
      CLNREVSTAT ARRAY<STRING>
      ALLELEID INTEGER
    """
    from google.cloud import bigquery

    bq_table = clinvar_config["bq_table"]
    logger.info(f"Querying ClinVar: {bq_table}")
    client = bigquery.Client()

    variants_for_query = df[["chrom", "pos", "ref", "alt"]].drop_duplicates()
    n_variants = len(variants_for_query)
    logger.info(f"Querying {n_variants} unique variants against ClinVar")

    batch_size = 500
    clinvar_results = []

    for start in range(0, n_variants, batch_size):
        batch = variants_for_query.iloc[start:start + batch_size]

        conditions = []
        for _, row in batch.iterrows():
            chrom = str(row["chrom"]).replace("chr", "")
            # ClinVar table uses 0-based start_position; pipeline pos is 1-based
            conditions.append(
                f"(reference_name = '{chrom}' AND start_position = {int(row['pos']) - 1} "
                f"AND reference_bases = '{row['ref']}' AND alt_struct.alt = '{row['alt']}')"
            )

        if not conditions:
            continue

        where = " OR ".join(conditions)
        query = f"""
        SELECT
            reference_name,
            start_position + 1 AS pos,
            reference_bases AS ref,
            alt_struct.alt AS alt,
            ARRAY_TO_STRING(CLNSIG, '|') AS classification,
            ARRAY_TO_STRING(CLNREVSTAT, '|') AS review_status,
            ALLELEID AS clinvar_id
        FROM `{bq_table}`,
             UNNEST(alternate_bases) AS alt_struct
        WHERE {where}
        """

        try:
            batch_df = client.query(query).to_dataframe()
            clinvar_results.append(batch_df)
        except Exception as e:
            logger.warning(f"ClinVar query batch failed: {e}")
            continue

    if not clinvar_results or all(d.empty for d in clinvar_results):
        logger.info("No ClinVar matches found")
        df["clinvar_classification"] = None
        df["clinvar_review_stars"] = 0
        df["clinvar_id"] = None
        return df

    clinvar_df = pd.concat(clinvar_results, ignore_index=True)

    # Map CLNREVSTAT ('criteria_provided|_multiple_submitters|_no_conflicts'
    # -> 2 stars, 'reviewed_by_expert_panel' -> 3, etc.)
    from pipeline.utils._clinvar_stars import map_review_status
    clinvar_df["clinvar_review_stars"] = clinvar_df["review_status"].apply(
        map_review_status
    )

    # Normalize CLNSIG: pick the most-pathogenic value if multiple are present
    def _pick_sig(sig_str):
        if not sig_str:
            return None
        parts = [p.strip() for p in str(sig_str).split("|") if p.strip()]
        priority = ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
                    "Uncertain_significance", "Likely_benign", "Benign",
                    "Benign/Likely_benign"]
        for p in priority:
            if p in parts:
                return p.replace("_", " ")
        return parts[0].replace("_", " ") if parts else None

    clinvar_df["clinvar_classification"] = clinvar_df["classification"].apply(_pick_sig)

    clinvar_df["chrom"] = clinvar_df["reference_name"].apply(
        lambda x: f"chr{x}" if not str(x).startswith("chr") else str(x)
    )

    clinvar_df = clinvar_df[["chrom", "pos", "ref", "alt",
                             "clinvar_classification", "clinvar_review_stars",
                             "clinvar_id"]].drop_duplicates(
        subset=["chrom", "pos", "ref", "alt"], keep="first"
    )

    df["pos"] = df["pos"].astype(int)
    clinvar_df["pos"] = clinvar_df["pos"].astype(int)

    df = df.merge(clinvar_df, on=["chrom", "pos", "ref", "alt"], how="left")
    n_matches = df["clinvar_classification"].notna().sum()
    logger.info(f"ClinVar: {n_matches} variants matched")

    return df


def _enrich_gnomad(df: pd.DataFrame, gnomad_config: dict) -> pd.DataFrame:
    """Enrich variants with gnomAD population frequencies.

    Targets the per-chromosome v3 genomes layout on bigquery-public-data:
        bigquery-public-data.gnomAD.v3_genomes__chr{N}
    Schema: reference_name STRING ('chr17'), start_position INTEGER (0-based),
    reference_bases STRING, alternate_bases ARRAY<STRUCT<alt, AF, AC, ...>>.

    Config supports a {chrom} placeholder in `bq_table_pattern`, e.g.
        bq_table_pattern: bigquery-public-data.gnomAD.v3_genomes__chr{chrom}
    Falls back to the legacy single-table `bq_table` for backward compat.
    """
    from google.cloud import bigquery

    pattern = gnomad_config.get("bq_table_pattern") or gnomad_config.get("bq_table", "")
    if "{chrom}" not in pattern:
        logger.info("gnomAD table has no {chrom} placeholder — skipping enrichment")
        df["gnomad_af"] = None
        df["gnomad_af_popmax"] = None
        return df

    logger.info(f"Querying gnomAD pattern: {pattern}")
    client = bigquery.Client()

    df["pos"] = df["pos"].astype(int)
    gnomad_results = []

    # Group queries by chromosome — one BQ query per chrom
    for chrom_prefixed, group in df[["chrom", "pos", "ref", "alt"]].drop_duplicates().groupby("chrom"):
        chrom_bare = str(chrom_prefixed).replace("chr", "")
        bq_table = pattern.format(chrom=chrom_bare)

        conditions = []
        for _, row in group.iterrows():
            # gnomAD v3 start_position is 0-based
            conditions.append(
                f"(start_position = {int(row['pos']) - 1} "
                f"AND reference_bases = '{row['ref']}' AND ab.alt = '{row['alt']}')"
            )
        if not conditions:
            continue

        where = " OR ".join(conditions)
        query = f"""
        SELECT
            reference_name,
            start_position + 1 AS pos,
            reference_bases AS ref,
            ab.alt AS alt,
            ab.AF AS gnomad_af,
            ab.AC AS gnomad_ac,
            AN AS gnomad_an
        FROM `{bq_table}`,
             UNNEST(alternate_bases) AS ab
        WHERE {where}
        """

        try:
            batch_df = client.query(query).to_dataframe()
            gnomad_results.append(batch_df)
        except Exception as e:
            logger.warning(f"gnomAD query for {chrom_prefixed} ({bq_table}) failed: {e}")
            continue

    if not gnomad_results or all(d.empty for d in gnomad_results):
        logger.info("No gnomAD matches found")
        df["gnomad_af"] = None
        df["gnomad_af_popmax"] = None
        return df

    gnomad_df = pd.concat(gnomad_results, ignore_index=True)
    gnomad_df["chrom"] = gnomad_df["reference_name"].apply(
        lambda x: x if str(x).startswith("chr") else f"chr{x}"
    )
    # popmax not present in this schema; leave as None
    gnomad_df["gnomad_af_popmax"] = None
    gnomad_df = gnomad_df[["chrom", "pos", "ref", "alt",
                           "gnomad_af", "gnomad_af_popmax"]].drop_duplicates(
        subset=["chrom", "pos", "ref", "alt"], keep="first"
    )

    gnomad_df["pos"] = gnomad_df["pos"].astype(int)

    df = df.merge(gnomad_df, on=["chrom", "pos", "ref", "alt"], how="left")
    n_matches = df["gnomad_af"].notna().sum()
    logger.info(f"gnomAD: {n_matches} variants matched")

    return df


def _enrich_svs(df: pd.DataFrame, db_config: dict) -> pd.DataFrame:
    """Enrich structural variants with gnomAD-SV AF + gene-content count.

    Rolls SV rows through:
      1. gnomAD-SV BQ lookup if configured (databases.gnomad_sv.bq_table)
      2. Gene count from bundled ClinGen DS reference (counts HI/TS gene
         overlaps as a rough lower bound — production deployments should
         join against a refseq/Ensembl gene-coords table)
    """
    gnomad_sv_cfg = db_config.get("gnomad_sv", {}) or {}
    df["gnomad_sv_af"] = None
    df["gnomad_sv_ac"] = None

    if gnomad_sv_cfg.get("bq_table"):
        df = _enrich_gnomad_sv(df, gnomad_sv_cfg)

    # Count genes overlapped per SV using bundled ClinGen DS table as a
    # lower-bound approximation. (A full implementation would query refseq
    # exons or Ensembl regulatory tables.)
    from pipeline.utils.cnv_engine import _load_clingen_ds, _interval_overlaps
    ds = _load_clingen_ds()
    hi_genes = ds.get("_haploinsufficient_3", []) + ds.get("_triplosensitive_3", [])

    n_genes_list = []
    for _, row in df.iterrows():
        chrom = row.get("chrom")
        start = int(row.get("pos") or 0)
        end_pos = row.get("end_pos")
        end = int(end_pos) if end_pos is not None else start
        n = sum(1 for g in hi_genes
                if g["chrom"] == chrom
                and _interval_overlaps(start, end, g["start"], g["end"]))
        n_genes_list.append(n)
    df["n_genes"] = n_genes_list

    logger.info(f"SV enrichment: {len(df)} SVs annotated with gene counts + gnomAD-SV AF")
    return df


def _enrich_gnomad_sv(df: pd.DataFrame, gnomad_sv_cfg: dict) -> pd.DataFrame:
    """Query gnomAD-SV BigQuery table for SV allele frequencies.

    The schema varies by deployment; this implementation assumes a
    `(chrom, start, end, svtype, AF, AC)` projection. Falls back gracefully
    on schema mismatches.
    """
    from google.cloud import bigquery
    bq_table = gnomad_sv_cfg["bq_table"]
    logger.info(f"Querying gnomAD-SV: {bq_table}")
    client = bigquery.Client()

    n = len(df)
    if n == 0:
        return df

    # Build a UNION of small WHERE clauses per row (small N expected for SV)
    conditions = []
    for _, row in df.iterrows():
        chrom = str(row["chrom"]).replace("chr", "")
        start = int(row["pos"])
        end = int(row["end_pos"] or start)
        svtype = row["svtype"] or ""
        # Match by chromosome + interval overlap fraction
        conditions.append(
            f"(reference_name='{chrom}' AND svtype='{svtype}' "
            f"AND start_position <= {end} AND end_position >= {start})"
        )

    where = " OR ".join(conditions) if conditions else "FALSE"
    query = f"""
    SELECT reference_name AS chrom, start_position AS pos,
           end_position AS end_pos, svtype, AF, AC
    FROM `{bq_table}`
    WHERE {where}
    """
    try:
        sv_df = client.query(query).to_dataframe()
        if not sv_df.empty:
            # Take max AF across overlapping SVs (conservative: assume the
            # CNV in our VCF could be the same as any of these)
            for idx, row in df.iterrows():
                chrom = str(row["chrom"]).replace("chr", "")
                matches = sv_df[(sv_df["chrom"] == chrom)
                                 & (sv_df["svtype"] == row["svtype"])]
                if not matches.empty:
                    df.at[idx, "gnomad_sv_af"] = float(matches["AF"].max())
                    df.at[idx, "gnomad_sv_ac"] = int(matches["AC"].max())
    except Exception as e:
        logger.warning(f"gnomAD-SV query failed: {e}")
    return df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    run(config)
