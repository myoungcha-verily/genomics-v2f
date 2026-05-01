"""Structural variant (SV) VCF adapter.

Parses VCFs with symbolic alleles (<DEL>, <DUP>, <INV>, <INS>, <BND>),
which the small-variant adapters drop. Surfaces SV-specific columns:

    svtype       (str)   : DEL | DUP | INV | INS | BND
    svlen        (int)   : Length in bp (positive for DEL/INV, negative for DUP
                           in some callers; we normalize to abs value here)
    end_pos      (int)   : End coordinate (from INFO.END or POS+SVLEN)
    bnd_partner  (str)   : For BND, the partner breakend (chrom:pos)
    n_genes      (int)   : Number of protein-coding genes affected (filled at
                           stage 3 enrichment via VEP / refseq overlap)

The adapter focuses on parsing — the ACMG/ClinGen 2019 (Riggs 2019) CNV
scoring lives in `pipeline/utils/cnv_engine.py` and runs at stage 4.
"""

import logging
import re
from typing import List

import pandas as pd

from adapters.base_adapter import BaseVCFAdapter

logger = logging.getLogger(__name__)

SV_ALT_PATTERN = re.compile(r"<([A-Z]+)(?::[A-Z]+)?>")
BND_PATTERN = re.compile(r"[\[\]]([\w\d.]+):(\d+)[\[\]]")


class SVAdapter(BaseVCFAdapter):
    """Adapter for VCFs containing symbolic structural-variant alleles."""

    def get_sample_ids(self) -> List[str]:
        local = self._resolve_vcf_path()
        try:
            from cyvcf2 import VCF
            v = VCF(local)
            samples = list(v.samples)
            v.close()
            return samples
        except ImportError:
            try:
                import pysam
                with pysam.VariantFile(local) as vf:
                    return list(vf.header.samples)
            except Exception:
                return []

    def get_proband_id(self) -> str:
        samples = self.get_sample_ids()
        if not samples:
            raise ValueError("No samples in SV VCF")
        return samples[0]

    def load_variants(self) -> pd.DataFrame:
        local = self._resolve_vcf_path()
        logger.info(f"Loading SVs from: {local}")
        try:
            from cyvcf2 import VCF  # noqa: F401
            return self._load_cyvcf2(local)
        except ImportError:
            return self._load_pysam(local)

    def _classify_alt(self, alt: str):
        """Return (svtype, bnd_partner) for an ALT allele.

        Returns (None, None) for non-symbolic alleles so caller can skip.
        """
        if alt and alt.startswith("<") and alt.endswith(">"):
            m = SV_ALT_PATTERN.match(alt)
            if m:
                return m.group(1), None
        # Breakend notation: e.g. N[chr2:321682[
        m = BND_PATTERN.search(alt or "")
        if m:
            return "BND", f"{m.group(1)}:{m.group(2)}"
        return None, None

    def _load_cyvcf2(self, vcf_path: str) -> pd.DataFrame:
        from cyvcf2 import VCF
        v = VCF(vcf_path)
        samples = list(v.samples)
        proband = samples[0] if samples else "unknown"

        records = []
        n_skipped_nonsv = 0
        for variant in v:
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            qual = variant.QUAL if variant.QUAL is not None else 0.0
            filt = variant.FILTER if variant.FILTER else "PASS"

            for alt in variant.ALT:
                svtype, bnd_partner = self._classify_alt(alt)
                if svtype is None:
                    n_skipped_nonsv += 1
                    continue

                # END position / SVLEN
                end_pos = (variant.INFO.get("END") if hasattr(variant, "INFO")
                            else None)
                svlen = variant.INFO.get("SVLEN") if hasattr(variant, "INFO") else None
                if isinstance(svlen, (list, tuple)):
                    svlen = svlen[0] if svlen else None
                if svlen is not None:
                    svlen = abs(int(svlen))
                if end_pos is None and svlen is not None and svtype != "BND":
                    end_pos = pos + svlen

                gt_field = "0/1"
                read_depth = 0
                if variant.genotypes:
                    g = variant.genotypes[0]
                    gt_field = f"{g[0]}/{g[1]}" if not g[2] else f"{g[0]}|{g[1]}"
                if hasattr(variant, "format") and variant.format("DP") is not None:
                    try:
                        read_depth = int(variant.format("DP")[0][0])
                    except (TypeError, IndexError):
                        read_depth = 0

                records.append({
                    "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                    "qual": qual, "filter": filt,
                    "variant_id": self._make_variant_id(chrom, pos, ref, alt),
                    "sample_id": proband,
                    "genotype": gt_field,
                    "gt_quality": 99.0,
                    "read_depth": read_depth,
                    "allele_depth_ref": 0,
                    "allele_depth_alt": 0,
                    "allele_fraction": 0.5,
                    "svtype": svtype,
                    "svlen": svlen,
                    "end_pos": end_pos,
                    "bnd_partner": bnd_partner,
                })
        v.close()
        df = pd.DataFrame(records)
        logger.info(f"Loaded {len(df)} structural variants "
                    f"(skipped {n_skipped_nonsv} non-symbolic alleles)")
        return df

    def _load_pysam(self, vcf_path: str) -> pd.DataFrame:
        import pysam
        records = []
        with pysam.VariantFile(vcf_path) as vf:
            samples = list(vf.header.samples)
            proband = samples[0] if samples else "unknown"
            for rec in vf:
                for alt in rec.alts or []:
                    svtype, bnd_partner = self._classify_alt(alt)
                    if svtype is None:
                        continue
                    end_pos = rec.info.get("END")
                    svlen = rec.info.get("SVLEN")
                    if isinstance(svlen, (list, tuple)):
                        svlen = svlen[0] if svlen else None
                    if svlen is not None:
                        svlen = abs(int(svlen))
                    if end_pos is None and svlen and svtype != "BND":
                        end_pos = rec.pos + svlen
                    records.append({
                        "chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref,
                        "alt": alt, "qual": rec.qual or 0.0,
                        "filter": ";".join(rec.filter.keys()) or "PASS",
                        "variant_id": self._make_variant_id(rec.chrom, rec.pos, rec.ref, alt),
                        "sample_id": proband,
                        "genotype": "0/1",
                        "gt_quality": 99.0, "read_depth": 0,
                        "allele_depth_ref": 0, "allele_depth_alt": 0,
                        "allele_fraction": 0.5,
                        "svtype": svtype, "svlen": svlen,
                        "end_pos": end_pos, "bnd_partner": bnd_partner,
                    })
        df = pd.DataFrame(records)
        logger.info(f"Loaded {len(df)} structural variants (pysam)")
        return df
