"""Stage 2: Variant Annotation with VEP

Input: data/vcf/variants.parquet (or normalized VCF)
Output: data/annotated/vep_output.tsv, data/annotated/annotated_variants.parquet

Steps:
1. Run VEP (Docker/local/skip)
2. Parse VEP output
3. Extract key annotations (consequence, gene, HGVS, scores)
4. Merge with variant DataFrame
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
from pipeline.utils.vep_runner import run_vep, parse_vep_output, create_mock_vep_output
from pipeline.utils.consequence_mapper import map_vep_consequence
from pipeline.utils.in_silico_scores import get_spliceai_max

logger = logging.getLogger(__name__)


def run(config: dict) -> dict:
    """Execute Stage 2: Variant Annotation."""
    t0 = time.time()
    output_dir = config.get("output", {}).get("output_dir", "data")
    vcf_dir = os.path.join(output_dir, "vcf")
    ann_dir = os.path.join(output_dir, "annotated")
    os.makedirs(ann_dir, exist_ok=True)

    logger.info("Stage 2: Variant Annotation")

    # Load variants from Stage 1
    variants_path = os.path.join(vcf_dir, "variants.parquet")
    if not os.path.exists(variants_path):
        raise FileNotFoundError(f"Stage 1 output not found: {variants_path}")

    variants_df = pd.read_parquet(variants_path)
    logger.info(f"Loaded {len(variants_df)} variants from Stage 1")

    # Run VEP annotation
    ann_config = config.get("annotation", {})
    vep_mode = ann_config.get("vep_mode", "docker")
    vep_output = os.path.join(ann_dir, "vep_output.tsv")

    if vep_mode == "skip":
        logger.info("VEP mode: skip — creating mock annotations")
        vep_output = create_mock_vep_output(variants_df, vep_output)
    else:
        # Need VCF file for VEP
        normalized_vcf = os.path.join(vcf_dir, "normalized.vcf.gz")
        input_vcf = normalized_vcf if os.path.exists(normalized_vcf) else \
            config.get("input", {}).get("vcf_path", "")

        try:
            vep_output = run_vep(input_vcf, vep_output, config)
        except (RuntimeError, FileNotFoundError) as e:
            logger.warning(f"VEP failed ({e}), falling back to mock annotations")
            vep_output = create_mock_vep_output(variants_df, vep_output)

    # Parse VEP output
    vep_df = parse_vep_output(vep_output)
    logger.info(f"Parsed {len(vep_df)} VEP annotations")

    # Extract and merge annotations
    annotated_df = _merge_annotations(variants_df, vep_df, config)
    logger.info(f"Annotated {len(annotated_df)} variants")

    # Save annotated variants
    output_path = os.path.join(ann_dir, "annotated_variants.parquet")
    annotated_df.to_parquet(output_path, index=False)

    # Summary stats
    consequence_counts = annotated_df["severity"].value_counts().to_dict() \
        if "severity" in annotated_df.columns else {}

    result = {
        "stage": "02_annotate_variants",
        "vep_mode": vep_mode,
        "input_variants": len(variants_df),
        "annotated_variants": len(annotated_df),
        "consequence_distribution": consequence_counts,
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    with open(os.path.join(ann_dir, "annotation_summary.json"), "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"\n{'='*60}")
    print(f"Stage 2 Complete: Variant Annotation")
    print(f"{'='*60}")
    print(f"  VEP mode:       {vep_mode}")
    print(f"  Input variants: {len(variants_df):,}")
    print(f"  Annotated:      {len(annotated_df):,}")
    if consequence_counts:
        for sev, count in sorted(consequence_counts.items()):
            print(f"  {sev:16s}: {count:,}")
    print(f"  Time:           {result['elapsed_seconds']}s")
    print(f"{'='*60}\n")

    return result


def _merge_annotations(variants_df: pd.DataFrame,
                        vep_df: pd.DataFrame,
                        config: dict) -> pd.DataFrame:
    """Merge VEP annotations with variant data."""
    # Build annotation lookup from VEP
    annotation_rows = []

    for _, vep_row in vep_df.iterrows():
        # Parse variant ID from VEP
        location = str(vep_row.get("Location", ""))
        allele = str(vep_row.get("Allele", ""))
        uploaded = str(vep_row.get("Uploaded_variation", ""))

        # Extract key fields
        consequence = str(vep_row.get("Consequence", ""))
        cons_info = map_vep_consequence(consequence)

        # SpliceAI max score
        spliceai_max = get_spliceai_max(
            vep_row.get("SpliceAI_pred_DS_AG", 0) or 0,
            vep_row.get("SpliceAI_pred_DS_AL", 0) or 0,
            vep_row.get("SpliceAI_pred_DS_DG", 0) or 0,
            vep_row.get("SpliceAI_pred_DS_DL", 0) or 0,
        )

        annotation_rows.append({
            "variant_id": uploaded,
            "gene": str(vep_row.get("SYMBOL", "") or ""),
            "transcript": str(vep_row.get("Feature", "") or ""),
            "consequence": cons_info["most_severe"],
            "all_consequences": consequence,
            "severity": cons_info["severity"],
            "is_lof": cons_info["is_lof"],
            "impact": str(vep_row.get("IMPACT", "") or ""),
            "hgvs_c": str(vep_row.get("HGVSc", "") or ""),
            "hgvs_p": str(vep_row.get("HGVSp", "") or ""),
            "canonical": str(vep_row.get("CANONICAL", "") or ""),
            "mane_select": str(vep_row.get("MANE_SELECT", "") or ""),
            "biotype": str(vep_row.get("BIOTYPE", "") or ""),
            "sift": str(vep_row.get("SIFT", "") or ""),
            "polyphen": str(vep_row.get("PolyPhen", "") or ""),
            "cadd_phred": vep_row.get("CADD_PHRED"),
            "revel": vep_row.get("REVEL"),
            "spliceai_max": spliceai_max if spliceai_max > 0 else None,
            "alphamissense": vep_row.get("am_pathogenicity"),
            "gnomad_af_vep": vep_row.get("gnomADe_AF"),
            "max_af": vep_row.get("MAX_AF"),
            "existing_variation": str(vep_row.get("Existing_variation", "") or ""),
            "domains": str(vep_row.get("DOMAINS", "") or ""),
        })

    ann_df = pd.DataFrame(annotation_rows)

    if ann_df.empty:
        logger.warning("No VEP annotations to merge")
        # Add empty annotation columns
        for col in ["gene", "consequence", "severity", "hgvs_c", "hgvs_p",
                     "cadd_phred", "revel", "spliceai_max"]:
            variants_df[col] = None
        return variants_df

    # Prefer canonical / MANE Select transcripts
    if "canonical" in ann_df.columns:
        canonical = ann_df[ann_df["canonical"] == "YES"]
        if not canonical.empty:
            # Keep canonical, fall back to any for variants without canonical
            canonical_ids = set(canonical["variant_id"])
            non_canonical = ann_df[~ann_df["variant_id"].isin(canonical_ids)]
            non_canonical = non_canonical.drop_duplicates(subset=["variant_id"],
                                                          keep="first")
            ann_df = pd.concat([canonical, non_canonical]).drop_duplicates(
                subset=["variant_id"], keep="first")
        else:
            ann_df = ann_df.drop_duplicates(subset=["variant_id"], keep="first")
    else:
        ann_df = ann_df.drop_duplicates(subset=["variant_id"], keep="first")

    # Merge
    merged = variants_df.merge(ann_df, on="variant_id", how="left")
    return merged


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    run(config)
