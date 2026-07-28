#!/usr/bin/env python3
"""
Assert that a screen reproduces the production apply exactly on the accepted pairs.

The per-image review bundle is rebuilt from a SCREEN, not from the production apply — the bundle re-runs
``evaluate_locked_pair`` using the config it reconstructs from that screen's ``method_spec.json``. If the
screen and the apply disagree, the maintainer reviews images of a threshold production never used, and
the review is worthless without anyone noticing.

Same frozen math + same inputs + same policy must give the same divisor. This asserts it rather than
assuming it, and reports what differs when it does not. Exit 1 on any mismatch.

    python scripts/gate5/check_screen_matches_production.py \
        --screen data/restore_pair_validation/screen_v12_production/ctrl-reference_only_div-p99
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
from phenocycler.config import load_config                                     # noqa: E402
from phenocycler.restore_validation import ACCEPTED_PAIRS, METHOD_VERSION      # noqa: E402

# Columns that must agree to the bit. `canonical_divisor` is the one the calls depend on; the others
# would each independently explain a divisor difference, so a mismatch localises the cause.
EXACT = ["canonical_divisor", "candidate_full_maximum", "reference_separator", "target_input_floor",
         "reference_input_floor"]
CATEGORICAL = ["state", "feature_compartments", "control_definition", "divisor_statistic"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=REPO / "config.ini")
    ap.add_argument("--screen", type=Path, required=True)
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    screen_csv = Path(a.screen) / "pair_evaluations.csv"
    spec_path = Path(a.screen) / "method_spec.json"
    for p in (screen_csv, spec_path, cfg.restore_apply_manifest_csv):
        if not p.is_file():
            raise SystemExit(f"[err] not found: {p}")

    spec = json.loads(spec_path.read_text())
    if spec.get("method_version") != METHOD_VERSION:
        raise SystemExit(f"[FAIL] screen is {spec.get('method_version')}, package is {METHOD_VERSION}")

    screen = pd.read_csv(screen_csv, dtype={"donor": str}, low_memory=False)
    prod = pd.read_csv(cfg.restore_apply_manifest_csv, dtype={"donor": str}, low_memory=False)
    accepted = set(ACCEPTED_PAIRS)
    key = ["donor", "target", "reference"]
    s = screen[[tuple(r) in {(t, r_) for t, r_ in accepted}
                for r in screen[["target", "reference"]].itertuples(index=False)]]
    m = s.merge(prod, on=key, suffixes=("_screen", "_prod"), validate="1:1")

    missing = len(prod) - len(m)
    print(f"[check] method {METHOD_VERSION} | accepted pairs {len(accepted)} | "
          f"production rows {len(prod)} | matched {len(m)}")
    failures = []
    if missing:
        failures.append(f"{missing} production rows have no screen counterpart")

    for col in EXACT:
        a_, b_ = m[f"{col}_screen"].to_numpy(float), m[f"{col}_prod"].to_numpy(float)
        both_nan = np.isnan(a_) & np.isnan(b_)
        differ = ~(both_nan | (a_ == b_))
        if differ.any():
            worst = m.loc[differ].assign(
                delta=np.abs(a_[differ] - b_[differ]),
                rel=np.abs(a_[differ] - b_[differ]) / np.maximum(np.abs(b_[differ]), 1e-12))
            failures.append(f"{col}: {int(differ.sum())}/{len(m)} rows differ "
                            f"(max |Δ| {worst['delta'].max():.6g}, max rel {worst['rel'].max():.3%})")
            print(worst.nlargest(5, "rel")[key + [f"{col}_screen", f"{col}_prod", "delta", "rel"]]
                  .to_string(index=False))
        else:
            print(f"  [ok] {col}: identical on all {len(m)} rows")

    for col in CATEGORICAL:
        if f"{col}_screen" not in m or f"{col}_prod" not in m:
            print(f"  [skip] {col}: not in both tables")
            continue
        differ = m[f"{col}_screen"].astype(str) != m[f"{col}_prod"].astype(str)
        if differ.any():
            failures.append(f"{col}: {int(differ.sum())}/{len(m)} rows differ")
            print(m.loc[differ, key + [f"{col}_screen", f"{col}_prod"]].head(10).to_string(index=False))
        else:
            print(f"  [ok] {col}: identical on all {len(m)} rows")

    if failures:
        print("\n[FAIL] the screen does NOT reproduce production:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"\n[PASS] the screen reproduces the production apply exactly on all {len(m)} accepted "
          "donor-pairs — the review bundle built from it shows the thresholds production used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
