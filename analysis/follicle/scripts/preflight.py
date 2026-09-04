#!/usr/bin/env python3
"""Check the external data contract without reading large measurements.

Exit 0 means all paths for the selected profile exist. Exit 2 is a complete,
machine-readable missing-path report; nothing is restored or generated.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


CODE_ROOT = Path(__file__).resolve().parent.parent
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from paths import (  # noqa: E402
    DATA_ROOT,
    OUTPUT_ROOT,
    REQUIRED_ANALYSIS_DATA,
    REQUIRED_QUPATH_DATA,
    missing_data_paths,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=("analysis", "qupath", "all"),
        default="all",
        help="path group to validate (default: all)",
    )
    args = parser.parse_args()

    required = []
    if args.profile in ("analysis", "all"):
        required.extend(REQUIRED_ANALYSIS_DATA)
    if args.profile in ("qupath", "all"):
        required.extend(REQUIRED_QUPATH_DATA)
    missing = missing_data_paths(*required)
    payload = {
        "status": "PASS" if not missing else "FAIL",
        "profile": args.profile,
        "code_root": str(CODE_ROOT),
        "data_root": str(DATA_ROOT),
        "output_root": str(OUTPUT_ROOT),
        "required": [str(DATA_ROOT / path) for path in required],
        "missing": [str(path) for path in missing],
    }
    print(json.dumps(payload, indent=2))
    return 0 if not missing else 2


if __name__ == "__main__":
    raise SystemExit(main())
