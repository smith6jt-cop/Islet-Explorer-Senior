"""Native 6551 sheet regression. Run in the original SCIMAP 2.3.5 environment."""
import os
from pathlib import Path

from scimap_all_cells import check_scoring


def test_native_scimap_scoring():
    default = Path(__file__).resolve().parents[3] / "core/panels/data/tissues/pancreatic_lymph_node/manifest.yml"
    manifest = Path(os.environ.get("PLN_PANEL_MANIFEST", str(default)))
    check_scoring(manifest)
