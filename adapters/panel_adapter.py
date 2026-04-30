"""Panel adapter — targeted gene panel analysis.

Extends single-sample adapter with gene panel filtering.
Loads gene lists from reference/gene_panels/ and filters variants
to only those within panel genes.
"""

import json
import logging
import os
from typing import List

import pandas as pd

from adapters.single_sample_adapter import SingleSampleAdapter

logger = logging.getLogger(__name__)


class PanelAdapter(SingleSampleAdapter):
    """Adapter for targeted gene panel VCF files."""

    def __init__(self, config: dict):
        super().__init__(config)
        panel_cfg = self.input_config.get("panel", {})
        self.panel_name = panel_cfg.get("gene_panel", "")
        self.padding_bp = panel_cfg.get("padding_bp", 20)
        self.gene_list = self._load_gene_list()

    def _load_gene_list(self) -> List[str]:
        """Load gene list from panel definition or BED file."""
        if not self.panel_name:
            return []

        # Check if it's a built-in panel name
        panel_dir = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "reference", "gene_panels"
        )
        panel_file = os.path.join(panel_dir, f"{self.panel_name}.json")

        if os.path.exists(panel_file):
            with open(panel_file) as f:
                panel_data = json.load(f)
            genes = panel_data.get("genes", [])
            logger.info(f"Loaded panel '{self.panel_name}': "
                        f"{len(genes)} genes")
            return genes

        # Check if it's a path to a custom file
        if os.path.exists(self.panel_name):
            if self.panel_name.endswith(".json"):
                with open(self.panel_name) as f:
                    panel_data = json.load(f)
                return panel_data.get("genes", [])
            elif self.panel_name.endswith(".bed"):
                return self._load_bed_genes(self.panel_name)
            elif self.panel_name.endswith(".txt"):
                with open(self.panel_name) as f:
                    return [line.strip() for line in f
                            if line.strip() and not line.startswith("#")]

        logger.warning(f"Panel '{self.panel_name}' not found. "
                        f"Available: {self._list_available_panels()}")
        return []

    def _load_bed_genes(self, bed_path: str) -> List[str]:
        """Extract gene names from BED file (assumes 4th column is gene name)."""
        genes = set()
        with open(bed_path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    genes.add(parts[3])
        return sorted(genes)

    def _list_available_panels(self) -> str:
        """List available built-in panels."""
        panel_dir = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "reference", "gene_panels"
        )
        if not os.path.exists(panel_dir):
            return "none"
        panels = [f.replace(".json", "") for f in os.listdir(panel_dir)
                  if f.endswith(".json") and f != "custom_template.json"]
        return ", ".join(panels) if panels else "none"

    def load_variants(self) -> pd.DataFrame:
        """Load variants and filter to panel genes."""
        df = super().load_variants()

        if self.gene_list:
            # Gene filtering happens after annotation (stage 2)
            # For now, store the gene list for later use
            logger.info(f"Panel loaded with {len(self.gene_list)} genes. "
                        f"Gene filtering will be applied after annotation.")

        return df

    def get_panel_info(self) -> dict:
        """Return panel metadata for reporting."""
        return {
            "panel_name": self.panel_name,
            "gene_count": len(self.gene_list),
            "genes": self.gene_list,
            "padding_bp": self.padding_bp,
        }
