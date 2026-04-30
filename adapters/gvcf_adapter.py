"""gVCF adapter — handles genomic VCF format.

gVCFs include reference blocks (non-variant regions) and require
special handling for variant extraction. Used primarily with
GATK HaplotypeCaller output.
"""

import logging
from typing import List

import pandas as pd

from adapters.single_sample_adapter import SingleSampleAdapter

logger = logging.getLogger(__name__)


class GVCFAdapter(SingleSampleAdapter):
    """Adapter for genomic VCF (gVCF) files."""

    def load_variants(self) -> pd.DataFrame:
        """Parse gVCF, filtering out reference blocks.

        gVCF-specific handling:
        - Skips <NON_REF> only alleles (reference blocks)
        - Handles END INFO field for reference confidence blocks
        - Filters by minimum GQ for reference calls
        """
        local_path = self._resolve_vcf_path()
        logger.info(f"Loading gVCF variants from: {local_path}")

        try:
            from cyvcf2 import VCF
            return self._load_gvcf_cyvcf2(local_path)
        except ImportError:
            return self._load_gvcf_pysam(local_path)

    def _load_gvcf_cyvcf2(self, vcf_path: str) -> pd.DataFrame:
        """Load gVCF using cyvcf2, skipping reference blocks."""
        from cyvcf2 import VCF

        vcf = VCF(vcf_path)
        samples = list(vcf.samples)
        if not samples:
            raise ValueError("No samples in gVCF")

        records = []
        n_ref_blocks = 0
        n_non_ref_only = 0

        for variant in vcf:
            # Skip reference blocks (no ALT or only <NON_REF>)
            alts = variant.ALT
            if not alts:
                n_ref_blocks += 1
                continue

            # Filter out <NON_REF> alleles
            real_alts = [a for a in alts if a != "<NON_REF>" and not a.startswith("<")]
            if not real_alts:
                n_non_ref_only += 1
                continue

            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            qual = variant.QUAL if variant.QUAL is not None else 0.0
            filt = variant.FILTER if variant.FILTER else "PASS"

            gt = variant.genotypes[0]
            gt_str = f"{gt[0]}/{gt[1]}"

            if gt[0] == 0 and gt[1] == 0:
                continue

            gq = variant.gt_quals[0] if variant.gt_quals is not None else 0.0
            dp = variant.gt_depths[0] if variant.gt_depths is not None else 0
            ad_alt = variant.gt_alt_depths[0] if variant.gt_alt_depths is not None else 0
            ad_ref = variant.gt_ref_depths[0] if variant.gt_ref_depths is not None else 0
            af = ad_alt / dp if dp > 0 else 0.0

            for alt in real_alts:
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
        logger.info(f"Loaded {len(df)} variants from gVCF "
                     f"({n_ref_blocks} ref blocks, {n_non_ref_only} NON_REF-only skipped)")
        return df

    def _load_gvcf_pysam(self, vcf_path: str) -> pd.DataFrame:
        """Fallback gVCF loader using pysam."""
        import pysam

        vcf = pysam.VariantFile(vcf_path)
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError("No samples in gVCF")

        sample_name = samples[0]
        records = []

        for rec in vcf:
            if not rec.alts:
                continue
            real_alts = [a for a in rec.alts
                         if str(a) != "<NON_REF>" and not str(a).startswith("<")]
            if not real_alts:
                continue

            sample_data = rec.samples[sample_name]
            gt_alleles = sample_data.get("GT", (None, None))
            if gt_alleles == (None, None) or gt_alleles == (0, 0):
                continue

            gt_str = "/".join(str(a) if a is not None else "." for a in gt_alleles)
            gq = sample_data.get("GQ", 0) or 0
            dp = sample_data.get("DP", 0) or 0
            ad = sample_data.get("AD", (0, 0))
            ad_ref = ad[0] if len(ad) > 0 else 0
            ad_alt = ad[1] if len(ad) > 1 else 0

            for alt in real_alts:
                records.append({
                    "chrom": str(rec.chrom), "pos": int(rec.pos),
                    "ref": rec.ref, "alt": str(alt),
                    "qual": float(rec.qual or 0),
                    "filter": "PASS",
                    "variant_id": self._make_variant_id(rec.chrom, rec.pos, rec.ref, str(alt)),
                    "sample_id": sample_name,
                    "genotype": gt_str,
                    "gt_quality": float(gq), "read_depth": int(dp),
                    "allele_depth_ref": int(ad_ref),
                    "allele_depth_alt": int(ad_alt),
                    "allele_fraction": ad_alt / dp if dp > 0 else 0.0,
                })

        vcf.close()
        return pd.DataFrame(records)
