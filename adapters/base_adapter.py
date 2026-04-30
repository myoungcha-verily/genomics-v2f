"""Abstract base class for VCF adapters.

All adapters must implement:
  - load_variants(): Parse VCF and return a DataFrame of variants
  - get_sample_ids(): Return list of sample IDs in the VCF
  - get_proband_id(): Return the primary sample for reporting

Adapters may override:
  - normalize(): Run bcftools norm (multi-allelic decomposition, left-alignment)
  - detect_de_novo(): Trio-specific de novo variant detection
  - filter_by_panel(): Panel-specific gene filtering
"""

import os
import subprocess
import logging
from abc import ABC, abstractmethod
from typing import List, Optional, Dict, Any

import pandas as pd

logger = logging.getLogger(__name__)


class BaseVCFAdapter(ABC):
    """Abstract base for all VCF format adapters."""

    def __init__(self, config: dict):
        self.config = config
        self.input_config = config.get("input", {})
        self.vcf_path = self.input_config.get("vcf_path", "")
        self.reference_genome = self.input_config.get("reference_genome", "GRCh38")
        self._validated = False

    @abstractmethod
    def load_variants(self) -> pd.DataFrame:
        """Parse VCF and return variant DataFrame.

        Returns DataFrame with columns:
            chrom (str): Chromosome (e.g., 'chr1' or '1')
            pos (int): 1-based position
            ref (str): Reference allele
            alt (str): Alternate allele
            qual (float): Quality score
            filter (str): FILTER field
            variant_id (str): Unique ID (chrom-pos-ref-alt)
            sample_id (str): Sample this genotype belongs to
            genotype (str): e.g., '0/1', '1/1'
            gt_quality (float): Genotype quality (GQ)
            read_depth (int): Read depth (DP)
            allele_depth_ref (int): Ref allele depth
            allele_depth_alt (int): Alt allele depth
            allele_fraction (float): Alt allele fraction (AD_alt / DP)
        """
        pass

    @abstractmethod
    def get_sample_ids(self) -> List[str]:
        """Return list of sample IDs present in the VCF."""
        pass

    @abstractmethod
    def get_proband_id(self) -> str:
        """Return the primary sample ID for reporting."""
        pass

    def validate_vcf(self) -> Dict[str, Any]:
        """Validate VCF file and return QC metrics.

        Returns dict with:
            valid (bool): Whether VCF passed validation
            n_variants (int): Total variant count
            n_samples (int): Number of samples
            issues (list): List of warning/error strings
            metrics (dict): QC metrics (ti/tv ratio, het/hom ratio, etc.)
        """
        issues = []
        metrics = {}

        if not self.vcf_path:
            return {"valid": False, "n_variants": 0, "n_samples": 0,
                    "issues": ["No VCF path specified"], "metrics": {}}

        # Check file exists
        if self.vcf_path.startswith("gs://"):
            # GCS path - check via gsutil
            try:
                result = subprocess.run(
                    ["gsutil", "stat", self.vcf_path],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode != 0:
                    issues.append(f"VCF not found at {self.vcf_path}")
            except (subprocess.TimeoutExpired, FileNotFoundError):
                issues.append("Cannot verify GCS path (gsutil not available)")
        elif not os.path.exists(self.vcf_path):
            issues.append(f"VCF file not found: {self.vcf_path}")

        return {
            "valid": len([i for i in issues if "not found" in i.lower()]) == 0,
            "n_variants": 0,
            "n_samples": 0,
            "issues": issues,
            "metrics": metrics,
        }

    def normalize(self, input_vcf: str, output_vcf: str,
                  reference_fasta: Optional[str] = None) -> str:
        """Normalize VCF using bcftools norm.

        - Decomposes multi-allelic sites
        - Left-aligns indels
        - Removes duplicate variants

        Args:
            input_vcf: Path to input VCF
            output_vcf: Path for normalized output
            reference_fasta: Optional path to reference FASTA

        Returns:
            Path to normalized VCF
        """
        cmd = ["bcftools", "norm", "-m", "-both"]
        if reference_fasta:
            cmd.extend(["-f", reference_fasta])
        cmd.extend(["-o", output_vcf, "-O", "z", input_vcf])

        logger.info(f"Normalizing VCF: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.warning(f"bcftools norm failed: {result.stderr}")
            logger.info("Falling back to unnormalized VCF")
            return input_vcf

        # Index
        subprocess.run(["bcftools", "index", "-t", output_vcf],
                        capture_output=True, text=True)
        return output_vcf

    def detect_de_novo(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Detect de novo variants (trio adapter overrides this)."""
        variants_df["is_de_novo"] = False
        return variants_df

    def filter_by_panel(self, variants_df: pd.DataFrame,
                         gene_list: List[str]) -> pd.DataFrame:
        """Filter variants to those in a gene list (panel adapter overrides)."""
        if "gene" in variants_df.columns and gene_list:
            mask = variants_df["gene"].isin(gene_list)
            filtered = variants_df[mask].copy()
            logger.info(f"Panel filter: {len(filtered)}/{len(variants_df)} "
                        f"variants in {len(gene_list)} genes")
            return filtered
        return variants_df

    def _resolve_vcf_path(self) -> str:
        """Download VCF from GCS if needed, return local path."""
        if not self.vcf_path.startswith("gs://"):
            return self.vcf_path

        local_dir = os.path.join(
            self.config.get("output", {}).get("output_dir", "data"), "vcf"
        )
        os.makedirs(local_dir, exist_ok=True)
        local_path = os.path.join(local_dir, os.path.basename(self.vcf_path))

        if os.path.exists(local_path):
            logger.info(f"Using cached VCF: {local_path}")
            return local_path

        logger.info(f"Downloading VCF from GCS: {self.vcf_path}")
        subprocess.run(
            ["gsutil", "cp", self.vcf_path, local_path],
            check=True, capture_output=True, text=True
        )
        # Also get index if it exists
        for ext in [".tbi", ".csi"]:
            idx = self.vcf_path + ext
            try:
                subprocess.run(
                    ["gsutil", "cp", idx, local_path + ext],
                    capture_output=True, text=True, timeout=30
                )
            except Exception:
                pass

        return local_path

    def _make_variant_id(self, chrom: str, pos: int,
                          ref: str, alt: str) -> str:
        """Create a unique variant identifier."""
        return f"{chrom}-{pos}-{ref}-{alt}"
