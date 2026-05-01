"""Somatic VCF adapter for tumor-only or tumor-normal-pair calls.

Differences from the germline adapters:
- The "proband" concept is replaced by a tumor sample (and optionally a
  matched normal). Sample assignment can be explicit (config) or inferred
  from VCF sample names (heuristic: 'tumor' / 'normal' substrings).
- Each row carries somatic-specific fields: tumor_vaf, normal_vaf,
  t_alt_count, n_alt_count. These let the AMP engine apply VAF cutoffs
  and identify likely tumor-only vs germline contamination calls.
- Skips the trio de-novo path; somatic mutations are by definition
  not transmitted.
"""

import logging
from typing import List, Optional, Tuple

import pandas as pd

from adapters.base_adapter import BaseVCFAdapter

logger = logging.getLogger(__name__)


class SomaticAdapter(BaseVCFAdapter):
    """Adapter for somatic VCFs (tumor-only or tumor-normal pair)."""

    def __init__(self, config: dict):
        super().__init__(config)
        somatic_cfg = (config.get("input", {}) or {}).get("somatic", {}) or {}
        self.explicit_tumor = somatic_cfg.get("tumor_sample_id", "")
        self.explicit_normal = somatic_cfg.get("normal_sample_id", "")
        amp_cfg = config.get("acmg_amp", {}) or {}
        self.vaf_min = float(amp_cfg.get("vaf_min_somatic", 0.05))

    def get_sample_ids(self) -> List[str]:
        local_path = self._resolve_vcf_path()
        try:
            from cyvcf2 import VCF
            vcf = VCF(local_path)
            samples = list(vcf.samples)
            vcf.close()
            return samples
        except ImportError:
            try:
                import pysam
                with pysam.VariantFile(local_path) as vf:
                    return list(vf.header.samples)
            except Exception:
                return []

    def _resolve_tumor_normal(self) -> Tuple[Optional[str], Optional[str]]:
        """Decide which sample is tumor vs normal.

        Precedence:
          1. config.input.somatic.tumor_sample_id / normal_sample_id
          2. Heuristic: any sample whose name matches 'tumor'/'normal'
             (case-insensitive) substrings
          3. Fallback (single sample): treat the only sample as tumor
        """
        samples = self.get_sample_ids()
        if not samples:
            return None, None

        if self.explicit_tumor and self.explicit_tumor in samples:
            normal = self.explicit_normal if self.explicit_normal in samples else None
            return self.explicit_tumor, normal

        # Heuristic match
        tumor = next((s for s in samples if "tumor" in s.lower()
                      or s.lower().startswith("t_")
                      or s.lower().endswith("_t")), None)
        normal = next((s for s in samples if "normal" in s.lower()
                       or s.lower().startswith("n_")
                       or s.lower().endswith("_n")), None)

        if not tumor and len(samples) == 1:
            tumor = samples[0]

        if not tumor and len(samples) == 2:
            # Two samples, no obvious labels — pick first as tumor, second as normal
            tumor, normal = samples[0], samples[1]
            logger.warning(
                "Could not infer tumor/normal from sample names "
                f"({samples}); picking '{tumor}' as tumor and '{normal}' as normal. "
                "Set input.somatic.tumor_sample_id explicitly to override."
            )
        return tumor, normal

    def get_proband_id(self) -> str:
        tumor, _ = self._resolve_tumor_normal()
        if not tumor:
            raise ValueError(
                "Could not identify tumor sample in VCF. Set "
                "input.somatic.tumor_sample_id in config."
            )
        return tumor

    def load_variants(self) -> pd.DataFrame:
        local_path = self._resolve_vcf_path()
        logger.info(f"Loading somatic variants from: {local_path}")
        try:
            from cyvcf2 import VCF  # noqa: F401
            return self._load_cyvcf2(local_path)
        except ImportError:
            logger.info("cyvcf2 not available, falling back to pysam")
            return self._load_pysam(local_path)

    def _load_cyvcf2(self, vcf_path: str) -> pd.DataFrame:
        from cyvcf2 import VCF

        vcf = VCF(vcf_path)
        samples = list(vcf.samples)
        if not samples:
            raise ValueError("No samples in VCF")

        tumor, normal = self._resolve_tumor_normal()
        t_idx = samples.index(tumor) if tumor in samples else 0
        n_idx = samples.index(normal) if normal and normal in samples else None

        records = []
        n_skipped = 0
        n_below_vaf = 0
        for variant in vcf:
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            qual = variant.QUAL if variant.QUAL is not None else 0.0
            filt = variant.FILTER if variant.FILTER else "PASS"

            ad = variant.format("AD")
            dp = variant.format("DP")
            gts = variant.genotypes  # list of [a1, a2, phased]

            for i, alt in enumerate(variant.ALT):
                t_gt = gts[t_idx] if t_idx < len(gts) else [0, 0, False]
                if t_gt[0] == 0 and t_gt[1] == 0:
                    n_skipped += 1
                    continue
                t_dp = int(dp[t_idx][0]) if dp is not None and dp[t_idx][0] is not None else 0
                t_alt = int(ad[t_idx][i + 1]) if ad is not None and len(ad[t_idx]) > i + 1 else 0
                t_ref = int(ad[t_idx][0]) if ad is not None and len(ad[t_idx]) > 0 else 0
                t_vaf = (t_alt / t_dp) if t_dp else 0.0

                if t_vaf < self.vaf_min:
                    n_below_vaf += 1
                    continue

                n_dp = n_alt = 0
                n_vaf = None
                if n_idx is not None:
                    n_dp = int(dp[n_idx][0]) if dp is not None and dp[n_idx][0] is not None else 0
                    n_alt = int(ad[n_idx][i + 1]) if ad is not None and len(ad[n_idx]) > i + 1 else 0
                    n_vaf = (n_alt / n_dp) if n_dp else 0.0

                t_gt_str = f"{t_gt[0]}/{t_gt[1]}" if not t_gt[2] else f"{t_gt[0]}|{t_gt[1]}"

                records.append({
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                    "qual": qual, "filter": filt,
                    "variant_id": self._make_variant_id(chrom, pos, ref, alt),
                    "sample_id": tumor,
                    "genotype": t_gt_str,
                    "gt_quality": 99.0,  # not always available in somatic callers
                    "read_depth": t_dp,
                    "allele_depth_ref": t_ref,
                    "allele_depth_alt": t_alt,
                    "allele_fraction": t_vaf,
                    "tumor_vaf": t_vaf,
                    "tumor_alt_count": t_alt,
                    "normal_vaf": n_vaf,
                    "normal_alt_count": n_alt,
                    "is_paired": n_idx is not None,
                })

        vcf.close()
        df = pd.DataFrame(records)
        logger.info(f"Loaded {len(df)} somatic variants "
                    f"(skipped {n_skipped} hom-ref, {n_below_vaf} below VAF threshold)")
        return df

    def _load_pysam(self, vcf_path: str) -> pd.DataFrame:
        import pysam
        records = []
        n_skipped = 0
        n_below_vaf = 0
        with pysam.VariantFile(vcf_path) as vf:
            samples = list(vf.header.samples)
            tumor, normal = self._resolve_tumor_normal()
            for rec in vf:
                for i, alt in enumerate(rec.alts or []):
                    t_call = rec.samples.get(tumor)
                    if t_call is None:
                        continue
                    gt = t_call.get("GT")
                    if not gt or all(a == 0 or a is None for a in gt):
                        n_skipped += 1
                        continue
                    ad = t_call.get("AD") or (0, 0)
                    dp = t_call.get("DP") or 0
                    t_alt = ad[i + 1] if len(ad) > i + 1 else 0
                    t_ref = ad[0] if ad else 0
                    t_vaf = (t_alt / dp) if dp else 0.0
                    if t_vaf < self.vaf_min:
                        n_below_vaf += 1
                        continue

                    n_dp = n_alt = 0
                    n_vaf = None
                    if normal:
                        n_call = rec.samples.get(normal)
                        if n_call:
                            n_ad = n_call.get("AD") or (0, 0)
                            n_dp = n_call.get("DP") or 0
                            n_alt = n_ad[i + 1] if len(n_ad) > i + 1 else 0
                            n_vaf = (n_alt / n_dp) if n_dp else 0.0

                    records.append({
                        "chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref,
                        "alt": alt, "qual": rec.qual or 0.0,
                        "filter": ";".join(rec.filter.keys()) or "PASS",
                        "variant_id": self._make_variant_id(rec.chrom, rec.pos, rec.ref, alt),
                        "sample_id": tumor,
                        "genotype": "/".join(str(a) for a in gt if a is not None),
                        "gt_quality": float(t_call.get("GQ") or 0),
                        "read_depth": dp,
                        "allele_depth_ref": t_ref,
                        "allele_depth_alt": t_alt,
                        "allele_fraction": t_vaf,
                        "tumor_vaf": t_vaf,
                        "tumor_alt_count": t_alt,
                        "normal_vaf": n_vaf,
                        "normal_alt_count": n_alt,
                        "is_paired": normal is not None,
                    })
        df = pd.DataFrame(records)
        logger.info(f"Loaded {len(df)} somatic variants "
                    f"(skipped {n_skipped} hom-ref, {n_below_vaf} below VAF threshold)")
        return df
