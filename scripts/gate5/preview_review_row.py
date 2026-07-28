#!/usr/bin/env python3
"""
Render ONE review row through the exact production path, so a display change can be judged before
committing to a ~50-minute 34-pair rebuild.

This exists because two display defects — the Step-2 line drawn at the manuscript maximum, and marker
panels windowed on a percentile that lit up background — both shipped in a full bundle and were caught
only by a reviewer opening the finished PDFs. Rendering one row takes about a minute.

It calls the same functions ``build_review_bundle`` calls, in the same order, with the config read from
the screen's own ``method_spec.json``. Nothing here is a reimplementation.

    python scripts/gate5/preview_review_row.py --donor 6591 --target CD11b --reference Pan_Cytokeratin
    -> data/gate5/preview/<target>_from_<reference>__<donor>.png
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))

import matplotlib                                                              # noqa: E402
matplotlib.use("Agg")
from phenocycler import restore_pair_review as rpr                             # noqa: E402
from phenocycler.config import load_config                                     # noqa: E402
from phenocycler.restore_validation import (                                   # noqa: E402
    ExpressionReview, PairReview, evaluate_locked_pair, load_expression_reviews,
    load_pair_reviews, load_validation_sample,
)

DEFAULT_SCREEN = ("data/restore_pair_validation/screen_v12_production/"
                  "ctrl-reference_only_div-p99")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=REPO / "config.ini")
    ap.add_argument("--screen", type=Path, default=REPO / DEFAULT_SCREEN)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--sample", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=REPO / "data" / "gate5" / "preview")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    sample = a.sample or (cfg.restore_pair_validation_dir / "qupath_compartment_full_20donor")
    spec = json.loads((a.screen / "method_spec.json").read_text())
    config = rpr._config_from_specification(spec)
    reviews, _ = load_expression_reviews(a.screen / "expression_reviews.csv")
    pairs, _ = load_pair_reviews(a.screen / "pair_reviews.csv")

    print(f"[preview] loading donor {a.donor} …", flush=True)
    donor_df, _ = load_validation_sample(sample, cfg.cells_dir, cfg.cells_redsea_dir, a.donor)
    evaluation = evaluate_locked_pair(
        donor_df, a.donor, a.target, a.reference, config=config,
        expression_review=reviews.get((a.donor, a.target), ExpressionReview.UNREVIEWED),
        pair_review=pairs.get((a.target, a.reference), PairReview.UNREVIEWED),
    )
    stats = evaluation["target_full_stats"] or {}
    lo, hi = evaluation["display_target_bounds"]
    print(f"[preview] state={evaluation['state']}  divisor={stats.get('divisor')}  "
          f"manuscript_max={stats.get('maximum')}")
    print(f"[preview] y-axis [{lo:.1f}, {hi:.1f}] -> Step 2 sits at "
          f"{(float(stats['divisor']) - lo) / (hi - lo):.3f} of axis height" if stats else "")
    print(f"[preview] display floors: target={evaluation['target_input_floor']} "
          f"reference={evaluation['reference_input_floor']}")

    categories = rpr.review_categories(evaluation, len(donor_df))
    zoom_x, zoom_y = rpr.representative_zoom_bounds(donor_df, categories)
    marker_zoom = rpr.load_marker_zoom(cfg.images_dir, donor_df, a.target, a.reference, zoom_x, zoom_y)
    projected = rpr.project_review_crop(donor_df, evaluation,
                                        marker_zoom.x_bounds_um, marker_zoom.y_bounds_um)

    fig, axes, legend_axis = rpr._create_review_figure(1)
    rpr._draw_review_row(axes[0], donor_df, evaluation, categories, "preview",
                         marker_zoom, projected)
    rpr._add_review_legend(legend_axis)
    rpr._add_review_suptitle(fig, a.target, a.reference)

    a.outdir.mkdir(parents=True, exist_ok=True)
    out = a.outdir / f"{a.target}_from_{a.reference}__{a.donor}.png"
    # NOT bbox_inches="tight": _create_review_figure already reserves a title band, and tightening
    # expands the canvas around the unwrapped suptitle into a uselessly wide image.
    fig.savefig(out, dpi=110)
    print(f"[saved] {out}")
    meta = rpr._marker_zoom_metadata(marker_zoom, evaluation)
    print(f"[preview] target window {meta['target_display_window']} "
          f"| reference window {meta['reference_display_window']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
