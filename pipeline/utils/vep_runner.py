"""Ensembl VEP (Variant Effect Predictor) runner.

Supports:
- Docker mode: runs VEP in official Docker container
- Local mode: uses locally installed VEP
- Skip mode: assumes pre-annotated VCF
"""

import json
import logging
import os
import subprocess
from typing import Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)

# VEP output fields we need
VEP_FIELDS = [
    "Uploaded_variation", "Location", "Allele",
    "Gene", "Feature", "Feature_type", "Consequence",
    "cDNA_position", "CDS_position", "Protein_position",
    "Amino_acids", "Codons", "Existing_variation",
    "IMPACT", "DISTANCE", "STRAND", "FLAGS",
    "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID",
    "BIOTYPE", "CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL",
    "TSL", "APPRIS", "SIFT", "PolyPhen",
    "HGVSc", "HGVSp", "HGVSg",
    "CLIN_SIG", "SOMATIC", "PHENO",
    "CADD_PHRED", "CADD_RAW",
    "REVEL",
    "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
    "am_class", "am_pathogenicity",
    "gnomADe_AF", "gnomADg_AF",
    "MAX_AF", "MAX_AF_POPS",
]


def run_vep(input_vcf: str, output_tsv: str, config: dict) -> str:
    """Run VEP annotation on a VCF file.

    Args:
        input_vcf: Path to input VCF
        output_tsv: Path for output TSV
        config: Pipeline config dict

    Returns:
        Path to VEP output file
    """
    ann_config = config.get("annotation", {})
    mode = ann_config.get("vep_mode", "docker")

    if mode == "skip":
        logger.info("VEP mode is 'skip' — using input VCF as-is")
        return input_vcf
    elif mode == "docker":
        return _run_vep_docker(input_vcf, output_tsv, ann_config)
    elif mode == "local":
        return _run_vep_local(input_vcf, output_tsv, ann_config)
    else:
        raise ValueError(f"Unknown VEP mode: {mode}")


def _run_vep_docker(input_vcf: str, output_tsv: str,
                     ann_config: dict) -> str:
    """Run VEP using Docker container."""
    image = ann_config.get("vep_docker_image",
                           "ensemblorg/ensembl-vep:release_114.0")
    species = ann_config.get("vep_species", "homo_sapiens")
    assembly = ann_config.get("vep_assembly", "GRCh38")
    cache_dir = ann_config.get("vep_cache_dir", "")

    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)

    # Build VEP command
    vep_cmd = [
        "vep",
        "--input_file", f"/data/{os.path.basename(input_vcf)}",
        "--output_file", f"/output/{os.path.basename(output_tsv)}",
        "--format", "vcf",
        "--tab", "--no_stats",
        "--species", species,
        "--assembly", assembly,
        "--offline" if cache_dir else "--database",
        "--everything",
        "--force_overwrite",
        "--canonical",
        "--mane_select",
        "--hgvs",
        "--symbol",
        "--numbers",
        "--domains",
        "--regulatory",
        "--var_synonyms",
        "--show_ref_allele",
        "--pick_allele_gene",  # Pick one consequence per gene
    ]

    # Add cache if available
    if cache_dir:
        vep_cmd.extend(["--cache", "--dir_cache", "/cache"])

    # Add prediction plugins
    plugins = []
    if ann_config.get("enable_cadd", True):
        plugins.append("CADD")
    if ann_config.get("enable_revel", True):
        plugins.append("REVEL")
    if ann_config.get("enable_spliceai", True):
        plugins.append("SpliceAI")
    if ann_config.get("enable_alphamissense", True):
        plugins.append("AlphaMissense")

    # Docker command
    input_dir = os.path.dirname(os.path.abspath(input_vcf))
    output_dir = os.path.dirname(os.path.abspath(output_tsv))

    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{input_dir}:/data:ro",
        "-v", f"{output_dir}:/output",
    ]

    if cache_dir:
        docker_cmd.extend(["-v", f"{cache_dir}:/cache:ro"])

    docker_cmd.append(image)
    docker_cmd.extend(vep_cmd)

    logger.info(f"Running VEP Docker: {image}")
    logger.debug(f"Command: {' '.join(docker_cmd)}")

    try:
        result = subprocess.run(
            docker_cmd,
            capture_output=True, text=True,
            timeout=3600  # 1 hour timeout
        )
        if result.returncode != 0:
            logger.error(f"VEP Docker failed: {result.stderr}")
            raise RuntimeError(f"VEP failed: {result.stderr[:500]}")
        logger.info("VEP annotation completed successfully")
        return output_tsv
    except FileNotFoundError:
        logger.error("Docker not found. Install Docker or use vep_mode: local")
        raise
    except subprocess.TimeoutExpired:
        logger.error("VEP timed out after 1 hour")
        raise


def _run_vep_local(input_vcf: str, output_tsv: str,
                    ann_config: dict) -> str:
    """Run locally installed VEP."""
    species = ann_config.get("vep_species", "homo_sapiens")
    assembly = ann_config.get("vep_assembly", "GRCh38")
    cache_dir = ann_config.get("vep_cache_dir", "")

    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)

    cmd = [
        "vep",
        "--input_file", input_vcf,
        "--output_file", output_tsv,
        "--format", "vcf",
        "--tab", "--no_stats",
        "--species", species,
        "--assembly", assembly,
        "--everything",
        "--force_overwrite",
        "--canonical",
        "--mane_select",
        "--hgvs",
        "--symbol",
        "--pick_allele_gene",
    ]

    if cache_dir:
        cmd.extend(["--cache", "--dir_cache", cache_dir, "--offline"])
    else:
        cmd.append("--database")

    logger.info(f"Running VEP locally: {' '.join(cmd[:10])}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

    if result.returncode != 0:
        raise RuntimeError(f"VEP failed: {result.stderr[:500]}")

    logger.info("VEP annotation completed")
    return output_tsv


def parse_vep_output(vep_tsv: str) -> pd.DataFrame:
    """Parse VEP tab-delimited output into DataFrame.

    Handles VEP header lines (starting with ##) and column names (starting with #).
    """
    logger.info(f"Parsing VEP output: {vep_tsv}")

    # Find header line
    header_line = None
    skip_lines = 0
    with open(vep_tsv) as f:
        for i, line in enumerate(f):
            if line.startswith("#Uploaded"):
                header_line = line.strip().lstrip("#").split("\t")
                skip_lines = i
                break
            elif not line.startswith("#"):
                skip_lines = i
                break

    if header_line is None:
        # Try reading as TSV directly
        df = pd.read_csv(vep_tsv, sep="\t", comment="#", low_memory=False)
    else:
        df = pd.read_csv(vep_tsv, sep="\t", skiprows=skip_lines,
                          header=0, low_memory=False)

    # Clean column names
    df.columns = [c.lstrip("#").strip() for c in df.columns]

    # Parse numeric scores
    for col in ["CADD_PHRED", "CADD_RAW", "REVEL",
                "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
                "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
                "am_pathogenicity", "gnomADe_AF", "gnomADg_AF", "MAX_AF"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    logger.info(f"Parsed {len(df)} VEP annotations")
    return df


def create_mock_vep_output(variants_df: pd.DataFrame,
                            output_path: str) -> str:
    """Create mock VEP-like output for testing without VEP installed.

    Uses basic consequence prediction based on variant type.
    """
    logger.warning("Creating MOCK VEP output — install VEP for real annotations")

    records = []
    for _, v in variants_df.iterrows():
        ref_len = len(v["ref"])
        alt_len = len(v["alt"])

        # Basic consequence guess
        if ref_len == 1 and alt_len == 1:
            consequence = "missense_variant"
            impact = "MODERATE"
        elif ref_len > alt_len:
            if (ref_len - alt_len) % 3 == 0:
                consequence = "inframe_deletion"
                impact = "MODERATE"
            else:
                consequence = "frameshift_variant"
                impact = "HIGH"
        elif alt_len > ref_len:
            if (alt_len - ref_len) % 3 == 0:
                consequence = "inframe_insertion"
                impact = "MODERATE"
            else:
                consequence = "frameshift_variant"
                impact = "HIGH"
        else:
            consequence = "coding_sequence_variant"
            impact = "MODIFIER"

        records.append({
            "Uploaded_variation": v["variant_id"],
            "Location": f"{v['chrom']}:{v['pos']}",
            "Allele": v["alt"],
            "Gene": "",
            "Feature": "",
            "Feature_type": "Transcript",
            "Consequence": consequence,
            "IMPACT": impact,
            "SYMBOL": "",
            "CANONICAL": "YES",
            "HGVSc": "",
            "HGVSp": "",
            "CADD_PHRED": None,
            "REVEL": None,
            "SpliceAI_pred_DS_AG": None,
            "SpliceAI_pred_DS_AL": None,
            "SpliceAI_pred_DS_DG": None,
            "SpliceAI_pred_DS_DL": None,
            "am_pathogenicity": None,
            "gnomADe_AF": None,
            "MAX_AF": None,
        })

    df = pd.DataFrame(records)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    logger.info(f"Mock VEP output: {len(df)} annotations → {output_path}")
    return output_path
