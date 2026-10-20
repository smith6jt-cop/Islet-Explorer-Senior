#!/usr/bin/env python3
"""Recompute one native typing scale from one declared alignment arm."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
import alignment_block as AB  # noqa: E402
import panel  # noqa: E402
import run_step3 as STEP3  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--aligned-root", required=True, type=Path)
    ap.add_argument("--shifts", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--block-sidecars", type=Path,
                    default=Path("data/development/celltype/alignment_block"))
    ap.add_argument("--n-per-donor", type=int, default=200_000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    shift = json.loads(args.shifts.read_text())
    donors = list(map(str, shift["donors"]))
    seen = {p.name.split("=", 1)[1] for p in args.aligned_root.glob("donor_id=*") if p.is_dir()}
    if seen != set(donors):
        raise ValueError(f"alignment donor mismatch; missing={sorted(set(donors)-seen)}, "
                         f"extra={sorted(seen-set(donors))}")
    markers = panel.markers(include_variant=True)
    STEP3.ALIGNED = args.aligned_root
    AB.OUT = args.block_sidecars
    scale = STEP3.build_scale(donors, markers, n_per_donor=args.n_per_donor,
                              seed=args.seed, ratio=False)
    if scale.index.duplicated().any() or not {"floor", "top", "span"}.issubset(scale.columns):
        raise AssertionError("scale producer returned an invalid scale table")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    scale.sort_index().reset_index().to_parquet(args.out, index=False)
    meta = {
        "producer": "run_step3.build_scale",
        "aligned_root": str(args.aligned_root.resolve()),
        "shifts": str(args.shifts.resolve()), "shifts_sha256": sha256(args.shifts),
        "donors": donors, "markers_requested": markers,
        "markers_written": list(scale.index), "n_per_donor": args.n_per_donor,
        "seed": args.seed, "nuclear_ratio": False,
        "block_sidecars": str(args.block_sidecars.resolve()),
        "scale_file": str(args.out.resolve()), "scale_sha256": sha256(args.out),
        "values": "floors and p99.9 tops from the arm's aligned values; per donor at most 200,000 eligible pancreas-block cells",
        "code_revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=HERE.parents[3],
                                         check=True, capture_output=True, text=True).stdout.strip(),
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(json.dumps(meta, indent=2))
    print(f"{len(scale)} markers -> {args.out}")
    print(f"provenance -> {args.out.with_suffix(args.out.suffix + '.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
