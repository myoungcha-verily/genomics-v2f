"""Single-sample VCF adapter.

Handles standard single-sample VCF/VCF.gz files from WES or WGS.
"""

import logging
from typing import List

import pandas as pd

from adapters.base_adapter import BaseVCFAdapter

logger = logging.getLogger(__name__)


class SingleSampleAdapter(BaseVCFAdapter):
    """Adapter for single-sample VCF files."""

    def get_sample_ids(self) -> List[str]:
        """Return sample IDs from VCF header."""
        local_path = self._resolve_vcf_path()
        try:
            from cyvcf2 import VCF
            vcf = VCF(local_path)
            samples = list(vcf.samples)
            vcf.close()
            return samples
        except ImportError:
            return self._get_samples_pysam(local_path)

    def get_proband_id(self) -> str:
        """Return the single sample ID."""
        samples = self.get_sample_ids()
        if not samples:
            raise ValueError("No samples found in VCF")
        if len(samples) > 1:
            logger.warning(f"Multiple samples in VCF ({len(samples)}), "
                           f"using first: {samples[0]}")
        return samples[0]

    def load_variants(self) -> pd.DataFrame:
        """Parse single-sample VCF into variant DataFrame."""
        local_path = self._resolve_vcf_path()
        logger.info(f"Loading variants from: {local_path}")

        try:
            from cyvcf2 import VCF
            return self._load_cyvcf2(local_path)
        except ImportError:
            logger.info("cyvcf2 not available, falling back to pysam")
            return self._load_pysam(local_path)

    def _load_cyvcf2(self, vcf_path: str) -> pd.DataFrame:
        """Load variants using cyvcf2 (preferred, faster)."""
        from cyvcf2 import VCF

        vcf = VCF(vcf_path)
        samples = list(vcf.samples)
        if not samples:
            raise ValueError("No samples in VCF")

        records = []
        n_skipped = 0

        for variant in vcf:
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            qual = variant.QUAL if variant.QUAL is not None else 0.0
            filt = variant.FILTER if variant.FILTER else "PASS"

            for i, alt in enumerate(variant.ALT):
                # Get genotype info for first sample
                gt = variant.genotypes[0]  # [allele1, allele2, phased]
                gt_str = f"{gt[0]}/{gt[1]}" if not gt[2] else f"{gt[0]}|{gt[1]}"

                # Skip homozygous reference
                if gt[0] == 0 and gt[1] == 0:
                    n_skipped += 1
                    continue

                # Genotype quality
                gq = variant.gt_quals[0] if variant.gt_quals is not None else 0.0

                # Read depth
                dp = variant.gt_depths[0] if variant.gt_depths is not None else 0

                # Allele depths
                ad_ref, ad_alt = 0, 0
                if variant.gt_alt_depths is not None:
                    ad_alt = variant.gt_alt_depths[0]
                if variant.gt_ref_depths is not None:
                    ad_ref = variant.gt_ref_depths[0]

                af = ad_alt / dp if dp > 0 else 0.0

                records.append({
                    "chrom": str(chrom),
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "qual": float(qual),
                    "filter": filt if isinstance(filt, str) else "PASS",
                    "variant_id": self._make_variant_id(chrom, pos, ref, alt),
                    "sample_id": samples[0],
                    "genotype": gt_str,
                    "gt_quality": float(gq),
                    "read_depth": int(dp),
                    "allele_depth_ref": int(ad_ref),
                    "allele_depth_alt": int(ad_alt),
                    "allele_fraction": float(af),
                })

        vcf.close()
        df = pd.DataFrame(records)
        logger.info(f"Loaded {len(df)} variants ({n_skipped} hom-ref skipped) "
                     f"from {samples[0]}")
        return df

    def _load_pysam(self, vcf_path: str) -> pd.DataFrame:
        """Fallback loader using pysam."""
        import pysam

        vcf = pysam.VariantFile(vcf_path)
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError("No samples in VCF")

        records = []
        sample_name = samples[0]

        for rec in vcf:
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            qual = rec.qual if rec.qual is not None else 0.0

            filt_set = rec.filter
            filt = "PASS" if not filt_set or "PASS" in filt_set else ";".join(filt_set)

            sample_data = rec.samples[sample_name]
            gt_alleles = sample_data.get("GT", (None, None))

            if gt_alleles == (None, None) or gt_alleles == (0, 0):
                continue

            for alt in rec.alts or []:
                gt_str = "/".join(str(a) if a is not None else "." for a in gt_alleles)
                gq = sample_data.get("GQ", 0) or 0
                dp = sample_data.get("DP", 0) or 0
                ad = sample_data.get("AD", (0, 0))
                ad_ref = ad[0] if len(ad) > 0 else 0
                ad_alt = ad[1] if len(ad) > 1 else 0
                af = ad_alt / dp if dp > 0 else 0.0

                records.append({
                    "chrom": str(chrom),
                    "pos": int(pos),
                    "ref": ref,
                    "alt": str(alt),
                    "qual": float(qual),
                    "filter": filt,
                    "variant_id": self._make_variant_id(chrom, pos, ref, str(alt)),
                    "sample_id": sample_name,
                    "genotype": gt_str,
                    "gt_quality": float(gq),
                    "read_depth": int(dp),
                    "allele_depth_ref": int(ad_ref),
                    "allele_depth_alt": int(ad_alt),
                    "allele_fraction": float(af),
                })

        vcf.close()
        df = pd.DataFrame(records)
        logger.info(f"Loaded {len(df)} variants from {sample_name}")
        return df

    def _get_samples_pysam(self, vcf_path: str) -> List[str]:
        """Get sample names using pysam."""
        import pysam
        vcf = pysam.VariantFile(vcf_path)
        samples = list(vcf.header.samples)
        vcf.close()
        return samples
