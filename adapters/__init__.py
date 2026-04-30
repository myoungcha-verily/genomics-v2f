"""VCF adapter registry — maps adapter names to classes."""

from adapters.single_sample_adapter import SingleSampleAdapter
from adapters.trio_adapter import TrioAdapter
from adapters.panel_adapter import PanelAdapter
from adapters.gvcf_adapter import GVCFAdapter

ADAPTER_REGISTRY = {
    "single_sample": SingleSampleAdapter,
    "trio": TrioAdapter,
    "panel": PanelAdapter,
    "gvcf": GVCFAdapter,
}


def get_adapter(name: str, config: dict):
    """Get adapter instance by name."""
    if name not in ADAPTER_REGISTRY:
        available = ", ".join(ADAPTER_REGISTRY.keys())
        raise ValueError(f"Unknown adapter '{name}'. Available: {available}")
    return ADAPTER_REGISTRY[name](config)
