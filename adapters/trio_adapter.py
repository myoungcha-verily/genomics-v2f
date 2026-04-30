"""Trio VCF adapter — proband + parents.

Handles multi-sample VCF with family members for:
- De novo variant detection (PS2/PM6 ACMG criteria)
- Compound heterozygote detection (PM3)
- Inheritance pattern analysis
"""

import logging
from typing import List, Dict, Optional, Tuple

import pandas as pd

from adapters.base_adapter import BaseVCFAdapter

logger = logging.getLogger(__name__)


class TrioAdapter(BaseVCFAdapter):
    """Adapter for trio (proband + parents) VCF files."""

    def __init__(self, config: dict):
        super().__init__(config)
        trio_cfg = self.input_config.get("trio", {})
        self.proband_id = trio_cfg.get("proband_id", "")
        self.father_id = trio_cfg.get("father_id", "")
        self.mother_id = trio_cfg.get("mother_id", "")
        self.ped_file = trio_cfg.get("ped_file", "")

        if self.ped_file:
            self._parse_ped_file()

    def _parse_ped_file(self):
        """Parse PED file to get family relationships."""
        try:
            with open(self.ped_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        # PED format: family_id, individual_id, father_id, mother_id, sex, phenotype
                        ind_id = parts[1]
                        father = parts[2]
                        mother = parts[3]
                        phenotype = parts[5]
                        # Affected individual is the proband
                        if phenotype == "2":  # Affected
                            self.proband_id = ind_id
                            if father != "0":
                                self.father_id = father
                            if mother != "0":
                                self.mother_id = mother
            logger.info(f"PED: proband={self.proband_id}, "
                        f"father={self.father_id}, mother={self.mother_id}")
        except Exception as e:
            logger.warning(f"Failed to parse PED file: {e}")

    def get_sample_ids(self) -> List[str]:
        """Return all family member sample IDs."""
        local_path = self._resolve_vcf_path()
        try:
            from cyvcf2 import VCF
            vcf = VCF(local_path)
            samples = list(vcf.samples)
            vcf.close()
            return samples
        except ImportError:
            import pysam
            vcf = pysam.VariantFile(local_path)
            samples = list(vcf.header.samples)
            vcf.close()
            return samples

    def get_proband_id(self) -> str:
        """Return the proband sample ID."""
        if not self.proband_id:
            raise ValueError("Proband ID not specified. Set input.trio.proband_id "
                             "in config or provide a PED file.")
        return self.proband_id

    def load_variants(self) -> pd.DataFrame:
        """Parse trio VCF into variant DataFrame with all family members."""
        local_path = self._resolve_vcf_path()
        logger.info(f"Loading trio variants from: {local_path}")

        try:
            from cyvcf2 import VCF
            return self._load_trio_cyvcf2(local_path)
        except ImportError:
            return self._load_trio_pysam(local_path)

    def _load_trio_cyvcf2(self, vcf_path: str) -> pd.DataFrame:
        """Load trio variants using cyvcf2."""
        from cyvcf2 import VCF

        vcf = VCF(vcf_path)
        samples = list(vcf.samples)
        self._validate_trio_samples(samples)

        # Map sample names to indices
        sample_idx = {s: i for i, s in enumerate(samples)}
        proband_idx = sample_idx[self.proband_id]

        records = []
        for variant in vcf:
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            qual = variant.QUAL if variant.QUAL is not None else 0.0
            filt = variant.FILTER if variant.FILTER else "PASS"

            for alt in variant.ALT:
                # Get proband genotype
                gt = variant.genotypes[proband_idx]
                gt_str = f"{gt[0]}/{gt[1]}"

                if gt[0] == 0 and gt[1] == 0:
                    continue

                gq = variant.gt_quals[proband_idx] if variant.gt_quals is not None else 0.0
                dp = variant.gt_depths[proband_idx] if variant.gt_depths is not None else 0
                ad_alt = variant.gt_alt_depths[proband_idx] if variant.gt_alt_depths is not None else 0
                ad_ref = variant.gt_ref_depths[proband_idx] if variant.gt_ref_depths is not None else 0
                af = ad_alt / dp if dp > 0 else 0.0

                # Get parent genotypes
                father_gt, mother_gt = ".", "."
                father_dp, mother_dp = 0, 0

                if self.father_id and self.father_id in sample_idx:
                    fi = sample_idx[self.father_id]
                    fgt = variant.genotypes[fi]
                    father_gt = f"{fgt[0]}/{fgt[1]}"
                    father_dp = variant.gt_depths[fi] if variant.gt_depths is not None else 0

                if self.mother_id and self.mother_id in sample_idx:
                    mi = sample_idx[self.mother_id]
                    mgt = variant.genotypes[mi]
                    mother_gt = f"{mgt[0]}/{mgt[1]}"
                    mother_dp = variant.gt_depths[mi] if variant.gt_depths is not None else 0

                records.append({
                    "chrom": str(chrom),
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "qual": float(qual),
                    "filter": filt if isinstance(filt, str) else "PASS",
                    "variant_id": self._make_variant_id(chrom, pos, ref, alt),
                    "sample_id": self.proband_id,
                    "genotype": gt_str,
                    "gt_quality": float(gq),
                    "read_depth": int(dp),
                    "allele_depth_ref": int(ad_ref),
                    "allele_depth_alt": int(ad_alt),
                    "allele_fraction": float(af),
                    "father_genotype": father_gt,
                    "mother_genotype": mother_gt,
                    "father_depth": int(father_dp),
                    "mother_depth": int(mother_dp),
                })

        vcf.close()
        df = pd.DataFrame(records)
        logger.info(f"Loaded {len(df)} proband variants from trio")

        # Detect de novo variants
        df = self.detect_de_novo(df)
        return df

    def _load_trio_pysam(self, vcf_path: str) -> pd.DataFrame:
        """Fallback trio loader using pysam."""
        import pysam

        vcf = pysam.VariantFile(vcf_path)
        samples = list(vcf.header.samples)
        self._validate_trio_samples(samples)

        records = []
        for rec in vcf:
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            qual = rec.qual if rec.qual is not None else 0.0
            filt_set = rec.filter
            filt = "PASS" if not filt_set or "PASS" in filt_set else ";".join(filt_set)

            proband_data = rec.samples[self.proband_id]
            gt_alleles = proband_data.get("GT", (None, None))
            if gt_alleles == (None, None) or gt_alleles == (0, 0):
                continue

            for alt in rec.alts or []:
                gt_str = "/".join(str(a) if a is not None else "." for a in gt_alleles)
                gq = proband_data.get("GQ", 0) or 0
                dp = proband_data.get("DP", 0) or 0
                ad = proband_data.get("AD", (0, 0))
                ad_ref = ad[0] if len(ad) > 0 else 0
                ad_alt = ad[1] if len(ad) > 1 else 0

                father_gt, mother_gt = ".", "."
                father_dp, mother_dp = 0, 0

                if self.father_id in samples:
                    fd = rec.samples[self.father_id]
                    fgt = fd.get("GT", (None, None))
                    father_gt = "/".join(str(a) if a is not None else "." for a in fgt)
                    father_dp = fd.get("DP", 0) or 0

                if self.mother_id in samples:
                    md = rec.samples[self.mother_id]
                    mgt = md.get("GT", (None, None))
                    mother_gt = "/".join(str(a) if a is not None else "." for a in mgt)
                    mother_dp = md.get("DP", 0) or 0

                records.append({
                    "chrom": str(chrom), "pos": int(pos),
                    "ref": ref, "alt": str(alt),
                    "qual": float(qual), "filter": filt,
                    "variant_id": self._make_variant_id(chrom, pos, ref, str(alt)),
                    "sample_id": self.proband_id,
                    "genotype": gt_str,
                    "gt_quality": float(gq), "read_depth": int(dp),
                    "allele_depth_ref": int(ad_ref),
                    "allele_depth_alt": int(ad_alt),
                    "allele_fraction": ad_alt / dp if dp > 0 else 0.0,
                    "father_genotype": father_gt,
                    "mother_genotype": mother_gt,
                    "father_depth": int(father_dp),
                    "mother_depth": int(mother_dp),
                })

        vcf.close()
        df = pd.DataFrame(records)
        df = self.detect_de_novo(df)
        return df

    def detect_de_novo(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Detect de novo variants in proband.

        A variant is de novo if:
        - Present in proband (0/1 or 1/1)
        - Absent in both parents (0/0)
        - Both parents have adequate read depth (DP >= 10)
        """
        if "father_genotype" not in variants_df.columns:
            variants_df["is_de_novo"] = False
            variants_df["de_novo_confidence"] = "N/A"
            return variants_df

        min_parent_dp = 10

        is_de_novo = (
            (variants_df["father_genotype"] == "0/0") &
            (variants_df["mother_genotype"] == "0/0") &
            (variants_df["father_depth"] >= min_parent_dp) &
            (variants_df["mother_depth"] >= min_parent_dp)
        )

        # Confidence based on parent depth
        confidence = []
        for _, row in variants_df.iterrows():
            if not is_de_novo.loc[row.name]:
                confidence.append("N/A")
            elif row.get("father_depth", 0) >= 30 and row.get("mother_depth", 0) >= 30:
                confidence.append("HIGH")
            elif row.get("father_depth", 0) >= 20 and row.get("mother_depth", 0) >= 20:
                confidence.append("MEDIUM")
            else:
                confidence.append("LOW")

        variants_df["is_de_novo"] = is_de_novo
        variants_df["de_novo_confidence"] = confidence

        n_de_novo = is_de_novo.sum()
        logger.info(f"Detected {n_de_novo} de novo variants")
        return variants_df

    def _validate_trio_samples(self, vcf_samples: List[str]):
        """Validate that trio sample IDs are present in VCF."""
        missing = []
        if self.proband_id not in vcf_samples:
            missing.append(f"proband '{self.proband_id}'")
        if self.father_id and self.father_id not in vcf_samples:
            missing.append(f"father '{self.father_id}'")
        if self.mother_id and self.mother_id not in vcf_samples:
            missing.append(f"mother '{self.mother_id}'")

        if missing:
            available = ", ".join(vcf_samples)
            raise ValueError(
                f"Trio samples not found in VCF: {', '.join(missing)}. "
                f"Available samples: {available}"
            )
