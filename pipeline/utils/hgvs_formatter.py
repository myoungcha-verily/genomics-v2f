"""HGVS nomenclature formatting utilities."""

import re
import logging
from typing import Optional

logger = logging.getLogger(__name__)


def format_hgvs_c(transcript: str, hgvs_c: str) -> str:
    """Format coding HGVS notation.

    Example: NM_000059.4:c.5946del
    """
    if not hgvs_c:
        return ""
    if ":" in hgvs_c:
        return hgvs_c
    if transcript:
        return f"{transcript}:{hgvs_c}"
    return hgvs_c


def format_hgvs_p(transcript: str, hgvs_p: str) -> str:
    """Format protein HGVS notation.

    Example: NP_000050.3:p.Leu1982fs
    """
    if not hgvs_p:
        return ""
    if ":" in hgvs_p:
        return hgvs_p
    if transcript:
        return f"{transcript}:{hgvs_p}"
    return hgvs_p


def format_genomic(chrom: str, pos: int, ref: str, alt: str,
                    assembly: str = "GRCh38") -> str:
    """Format genomic HGVS notation.

    Example: NC_000013.11:g.32316461del
    """
    # Simple representation for display
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom}:g.{pos}{ref}>{alt}"
    elif len(ref) > len(alt):
        deleted = ref[len(alt):]
        if len(deleted) == 1:
            return f"{chrom}:g.{pos + len(alt)}del"
        return f"{chrom}:g.{pos + len(alt)}_{pos + len(ref) - 1}del"
    elif len(alt) > len(ref):
        inserted = alt[len(ref):]
        return f"{chrom}:g.{pos + len(ref) - 1}_{pos + len(ref)}ins{inserted}"
    return f"{chrom}:g.{pos}{ref}>{alt}"


def shorten_hgvs_p(hgvs_p: str) -> str:
    """Convert 3-letter amino acid codes to 1-letter for compact display.

    Example: p.Leu1982fs → p.L1982fs
    """
    aa_map = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*",
    }

    result = hgvs_p
    for three, one in aa_map.items():
        result = result.replace(three, one)
    return result


def extract_gene_from_hgvs(hgvs: str) -> Optional[str]:
    """Extract gene symbol from HGVS if present in parentheses.

    Example: NM_000059.4(BRCA2):c.5946del → BRCA2
    """
    match = re.search(r'\(([^)]+)\)', hgvs)
    return match.group(1) if match else None


def variant_display_name(gene: str, hgvs_p: Optional[str] = None,
                          hgvs_c: Optional[str] = None,
                          chrom: str = "", pos: int = 0,
                          ref: str = "", alt: str = "") -> str:
    """Create a human-readable variant name for display.

    Priority: Gene + protein change > Gene + coding change > genomic
    """
    if gene and hgvs_p and hgvs_p != "":
        short_p = shorten_hgvs_p(hgvs_p.split(":")[-1])
        return f"{gene} {short_p}"
    if gene and hgvs_c and hgvs_c != "":
        short_c = hgvs_c.split(":")[-1]
        return f"{gene} {short_c}"
    if gene:
        return f"{gene} {chrom}:{pos} {ref}>{alt}"
    return f"{chrom}:{pos} {ref}>{alt}"
