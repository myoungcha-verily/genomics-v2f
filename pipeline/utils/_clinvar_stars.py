"""ClinVar review status to star mapping helper."""

import json
import os

_STARS_MAP = None


def _load_map():
    global _STARS_MAP
    if _STARS_MAP is not None:
        return _STARS_MAP

    ref_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "reference", "clinvar_stars_mapping.json"
    )
    try:
        with open(ref_path) as f:
            data = json.load(f)
        _STARS_MAP = data.get("review_status_to_stars", {})
    except FileNotFoundError:
        _STARS_MAP = {}
    return _STARS_MAP


def map_review_status(status: str) -> int:
    """Map ClinVar review status string to number of review stars (0-4)."""
    if not status:
        return 0
    stars_map = _load_map()
    return stars_map.get(str(status).strip().lower(),
                          stars_map.get(str(status).strip(), 0))
