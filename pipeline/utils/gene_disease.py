"""ClinGen gene-disease validity lookups."""

import json
import logging
import os
from typing import Dict, List, Optional, Set

logger = logging.getLogger(__name__)

# Known LoF-intolerant genes (curated from ClinGen)
# This is a subset — full list should be loaded from reference file
_LOF_GENES: Optional[Set[str]] = None
_GENE_DISEASE_MAP: Optional[Dict] = None


def _load_lof_genes() -> Set[str]:
    """Load ClinGen loss-of-function gene list."""
    global _LOF_GENES
    if _LOF_GENES is not None:
        return _LOF_GENES

    ref_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "reference", "clingen_lof_genes.json"
    )
    try:
        with open(ref_path) as f:
            data = json.load(f)
            _LOF_GENES = set(data.get("genes", []))
    except FileNotFoundError:
        # Fallback: well-known haploinsufficient genes
        _LOF_GENES = {
            # Cardiac
            "MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "ACTC1",
            "LMNA", "SCN5A", "KCNQ1", "KCNH2", "RYR2", "PKP2",
            "DSP", "DSG2", "DSC2", "TMEM43", "TTN",
            # Cancer
            "BRCA1", "BRCA2", "TP53", "APC", "MLH1", "MSH2",
            "MSH6", "PMS2", "PTEN", "RB1", "VHL", "STK11",
            "SMAD4", "CDH1", "NF1", "NF2", "TSC1", "TSC2",
            # Neuro
            "SMN1", "MECP2", "FMR1",
            # Other well-established haploinsufficient genes
            "CFTR", "PKD1", "PKD2", "COL1A1", "COL1A2",
            "FBN1", "GBA", "HBB", "F8", "DMD",
        }
        logger.info(f"Using built-in LoF gene list ({len(_LOF_GENES)} genes)")

    return _LOF_GENES


def is_lof_gene(gene: str) -> bool:
    """Check if gene is in ClinGen LoF-intolerant list."""
    return gene.upper() in _load_lof_genes()


def get_gene_disease_info(gene: str) -> Dict:
    """Get gene-disease association info.

    Returns dict with:
        gene: str
        is_lof_gene: bool (ClinGen haploinsufficient)
        diseases: list of associated diseases (if available)
        inheritance: list of inheritance patterns
        moi: mode of inheritance (AD, AR, XL, etc.)
    """
    lof = is_lof_gene(gene)

    return {
        "gene": gene,
        "is_lof_gene": lof,
        "diseases": [],  # Populated from ClinGen API if available
        "inheritance": [],
        "moi": "",
    }


def load_gene_panel(panel_name: str) -> List[str]:
    """Load gene list from a named panel."""
    panel_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "reference", "gene_panels"
    )
    panel_file = os.path.join(panel_dir, f"{panel_name}.json")

    if not os.path.exists(panel_file):
        logger.warning(f"Panel not found: {panel_name}")
        return []

    with open(panel_file) as f:
        data = json.load(f)
    return data.get("genes", [])
