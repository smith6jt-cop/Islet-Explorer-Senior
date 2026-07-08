#!/usr/bin/env python3
"""
Apply a hormone-strength floor to the RESTORE positivity calls — the false-endocrine fix.

RESTORE's per-image threshold lands IN THE NOISE for β-loss (T1D) donors (no separated bright INS
population), so it calls thousands of acinar cells INS⁺ at `_norm ~1.1` (barely over threshold), 100%
in tissue, no islet coherence — while real β-cells sit at `_norm` 8–25. Requiring `_norm >= min_norm`
(threshold-relative, so per-image-adaptive) rejects the noise-floor false positives while keeping the
real, separated positives. Validated per-donor: 6380 (visually β-null) INS⁺ 11,905 → ~5 at K=5, while
6539 (ND) and 6593 (β-retained T1D) keep their real β. See scripts/senior/redsea_reassess_diag.py for
the yardstick and data/redsea_reassess/ for the figures.

This rewrites ONLY the `{m}_pos` boolean of the listed markers (default the three endocrine hormones);
`_norm` / `_log2r` and every other column are copied through unchanged. Point it at the current gated
dir and write to a new dir, then swap.

  python scripts/senior/apply_hormone_floor.py \
      --gated-dir data/restore_gated_redsea --out-dir data/restore_gated_redsea_k5 --min-norm 5
"""
from __future__ import annotations
import argparse
import glob
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[2] / "data"
DEFAULT_MARKERS = ["INS", "GCG", "SST"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gated-dir", type=Path, default=BASE / "restore_gated_redsea")
    ap.add_argument("--out-dir", type=Path, default=BASE / "restore_gated_redsea_k5")
    ap.add_argument("--min-norm", type=float, default=5.0,
                    help="require {m}_norm >= this (threshold-relative) for {m}_pos to stay True")
    ap.add_argument("--markers", default=",".join(DEFAULT_MARKERS),
                    help="comma-separated markers to floor (default the endocrine hormones INS,GCG,SST)")
    a = ap.parse_args(argv)
    markers = [m.strip() for m in a.markers.split(",") if m.strip()]

    files = sorted(glob.glob(str(a.gated_dir / "donor_id=*" / "*.parquet")))
    if not files:
        raise SystemExit(f"no parquet under {a.gated_dir}")
    for f in files:
        donor = f.split("donor_id=")[1].split("/")[0]
        df = pd.read_parquet(f)
        report = []
        for m in markers:
            pc, nc = f"{m}_pos", f"{m}_norm"
            if pc not in df.columns or nc not in df.columns:
                continue
            before = int(df[pc].sum())
            df[pc] = (df[nc].to_numpy() >= a.min_norm)
            report.append(f"{m} {before:,}->{int(df[pc].sum()):,}")
        out_dir = a.out_dir / f"donor_id={donor}"
        out_dir.mkdir(parents=True, exist_ok=True)
        df.to_parquet(out_dir / "data_0.parquet", index=False)
        print(f"[{donor}] {len(df):,} cells | " + " | ".join(report), flush=True)
    print(f"[done] floored {markers} at _norm>={a.min_norm} -> {a.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
