"""Measure every candidate divisor against each donor-marker's own control distribution.

For each accepted pair this runs the production `evaluate_locked_pair`, then reports
what each fixed divisor rule would have chosen on that pair's own reference-negative
control, and what fraction of cells each rule would call positive. It also keeps a
subsample of the control and of the full population so the figures can show the
distribution the rules are being applied to.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from phenocycler.config import load_config  # noqa: E402
from phenocycler.marker_registry import load_registry  # noqa: E402
from phenocycler.restore_stage import (  # noqa: E402
    RestoreStageConfig,
    _read_columns,
    build_restore_input_frame,
    production_pairs,
)
from phenocycler.restore_validation import (  # noqa: E402
    evaluate_locked_pair,
    load_expression_reviews,
    negative_control_statistics,
)

STATS = ("max", "p999", "p99", "p95")


def donor_frame(run_root: Path, donor: str, registry, cfg):
    pairs = production_pairs(registry)
    markers = tuple(dict.fromkeys(m for p in pairs for m in p))
    comps = [c for c in cfg.feature_compartments if c != "Cell"]
    cells = _read_columns(
        run_root / "cells",
        donor,
        ["object_id", "image", "X_centroid", "Y_centroid", *markers,
         *(f"{c}__{m}" for c in comps for m in markers)],
        allow_missing=True,
    )
    redsea = _read_columns(
        run_root / "redsea",
        donor,
        ["object_id", *(f"Cell__{m}" for m in markers),
         *(f"{c}__{m}" for c in comps for m in markers)],
        allow_missing=True,
    )
    geometry = _read_columns(
        run_root / "geometry_qc", donor, ["object_id", "analysis_eligible"]
    )
    return build_restore_input_frame(
        cells, redsea, geometry,
        donor_id=donor, pairs=pairs,
        feature_compartments=cfg.feature_compartments,
    ), pairs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="3bf4fc8e391d007766d1")
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pipeline_cfg = load_config(REPO / "config.ini")
    registry = load_registry(pipeline_cfg.marker_registry)
    stage_cfg = RestoreStageConfig.from_pipeline(pipeline_cfg)
    pair_cfg = stage_cfg.pair_config()
    reviews, _ = load_expression_reviews(pipeline_cfg.restore_expression_reviews)

    out = Path(args.out)
    (out / "dist").mkdir(parents=True, exist_ok=True)
    rows = []
    rng = np.random.default_rng(0)

    for donor in args.donors:
        frame, pairs = donor_frame(
            REPO / "data" / "runs" / args.run, donor, registry, stage_cfg
        )
        print(f"[{donor}] {len(frame):,} cells", flush=True)
        for target, reference in pairs:
            review = reviews.get((donor, target))
            kw = {"expression_review": review} if review is not None else {}
            try:
                ev = evaluate_locked_pair(
                    frame, donor, target, reference, config=pair_cfg, **kw
                )
            except Exception as exc:  # noqa: BLE001
                rows.append({"donor": donor, "target": target,
                             "reference": reference, "state": f"ERROR:{exc}"})
                continue

            # `target_redsea` is the QC-retained subset; the `full_*` masks and
            # `call_idx` index the whole donor frame. Keep everything full-length.
            tgt_cor = frame[f"redsea__{target}"].to_numpy(dtype=float)
            control = np.asarray(ev["full_reference_control"], dtype=bool)
            call_idx = ev.get("call_idx")
            callable_mask = np.zeros(len(tgt_cor), dtype=bool)
            if call_idx is not None:
                callable_mask[np.asarray(call_idx, dtype=int)] = True
            else:
                callable_mask[np.asarray(ev["valid_idx"], dtype=int)] = True
            ctrl_vals = tgt_cor[control & np.isfinite(tgt_cor)]

            row = {
                "donor": donor, "target": target, "reference": reference,
                "state": ev["state"].value if hasattr(ev["state"], "value") else str(ev["state"]),
                "target_input_floor": float(ev["target_input_floor"]),
                "control_n": int(ctrl_vals.size),
                "callable_n": int(callable_mask.sum()),
                "target_n": int(ev["target_n"]),
                "reference_n": int(ev["reference_n"]),
                "target_anchor_n": int(ev["target_anchor_n"]),
                "target_arm_saturated": bool(ev["target_arm_saturated"]),
                "reference_arm_saturated": bool(ev["reference_arm_saturated"]),
                "double_high_fraction": float(ev["double_high_fraction"]),
                "production_divisor": float(ev["canonical_divisor"])
                if ev["canonical_divisor"] is not None and np.isfinite(ev["canonical_divisor"] or np.nan)
                else np.nan,
            }
            if ctrl_vals.size:
                cand = {}
                for stat in STATS:
                    try:
                        cand[stat] = float(
                            negative_control_statistics(ctrl_vals, stat)["divisor"]
                        )
                    except Exception:  # noqa: BLE001
                        cand[stat] = np.nan
                cand["mean3sd"] = float(ctrl_vals.mean() + 3 * ctrl_vals.std())
                callv = tgt_cor[callable_mask]
                callv = callv[np.isfinite(callv)]
                for name, dv in cand.items():
                    row[f"div_{name}"] = dv
                    row[f"pos_{name}"] = (
                        float((callv > dv).mean()) if callv.size and np.isfinite(dv) else np.nan
                    )
                    row[f"ctrl_le1_{name}"] = (
                        float((ctrl_vals <= dv).mean()) if np.isfinite(dv) else np.nan
                    )
                    row[f"anchor_ratio_{name}"] = (
                        float((callv > dv).sum() / row["target_anchor_n"])
                        if row["target_anchor_n"] and np.isfinite(dv) else np.nan
                    )
                keep_c = ctrl_vals if ctrl_vals.size <= 60000 else rng.choice(ctrl_vals, 60000, replace=False)
                allv = tgt_cor[callable_mask & np.isfinite(tgt_cor)]
                keep_a = allv if allv.size <= 120000 else rng.choice(allv, 120000, replace=False)
                np.savez_compressed(
                    out / "dist" / f"{donor}__{target}.npz",
                    control=keep_c.astype(np.float32),
                    population=keep_a.astype(np.float32),
                    floor=np.float32(row["target_input_floor"]),
                    **{f"div_{k}": np.float32(v) for k, v in cand.items()},
                )
            rows.append(row)
            print(f"  {target:<16} <- {reference:<16} {row['state']:<12} "
                  f"ctrl={row['control_n']:>7,}", flush=True)

    table = pd.DataFrame(rows)
    table.to_parquet(out / "divisor_survey.parquet", index=False)
    print(f"\nwrote {out/'divisor_survey.parquet'}  ({len(table)} pairs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
