"""Code, data, and output roots for the absorbed spleen-follicle lane.

The tracked implementation lives below ``analysis/follicle`` in the monorepo,
while the multi-gigabyte QuPath project and measurement exports remain outside
Git.  Keeping those axes explicit prevents a checkout location from silently
changing which biological files are read.
"""

from __future__ import annotations

import os
from pathlib import Path


CODE_ROOT = Path(__file__).resolve().parent
DATA_ROOT = (
    Path(os.environ.get("FOLLICLE_DATA_ROOT", "/data/SpleenFollicleCounterQP"))
    .expanduser()
    .resolve()
)
OUTPUT_ROOT = (
    Path(os.environ.get("FOLLICLE_OUTPUT_ROOT", str(DATA_ROOT))).expanduser().resolve()
)

REQUIRED_ANALYSIS_DATA = (
    Path("Measurements/AnnotationsFinal.csv"),
    Path("Measurements/ForSH2B3.csv"),
    Path("Measurements/Cells.csv"),
    Path("Measurements/Spleen_Upload2_2024Sept17_v2_distances.csv"),
)
REQUIRED_QUPATH_DATA = (Path("Images"), Path("data"))


def missing_data_paths(*relative_paths: str | Path) -> tuple[Path, ...]:
    """Return every absent path below :data:`DATA_ROOT`."""

    paths = tuple(DATA_ROOT / Path(relative) for relative in relative_paths)
    return tuple(path for path in paths if not path.exists())


def require_data_paths(*relative_paths: str | Path) -> tuple[Path, ...]:
    """Return required paths below :data:`DATA_ROOT`, or fail as one preflight."""

    paths = tuple(DATA_ROOT / Path(relative) for relative in relative_paths)
    missing = missing_data_paths(*relative_paths)
    if missing:
        rendered = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(
            "spleen-follicle data preflight failed; set FOLLICLE_DATA_ROOT if "
            f"the external project moved:\n{rendered}"
        )
    return paths
