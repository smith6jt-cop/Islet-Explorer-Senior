"""Step 2: derive and apply the UniFORM shift for every marker, then check it.

Two passes.  The first walks markers, because a shift is a per-marker object estimated from
all 20 donors at once; the second walks donors, because that is how the output is written.

    python scripts/celltype/run_step2.py             # from core/phenocycler/
    python scripts/celltype/run_step2.py --dry-run   # derive and check, write no parquet
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import pathlib
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
import pyarrow.parquet as pq             # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                             # noqa: E402
import usability as U                    # noqa: E402
import uniform_adapter as UA             # noqa: E402

OUT = panel.OUT / "uniform"
FIGS = panel.OUT / "figures_step2"

#: Acceptance, declared before the run.  Stated on the quantity the shift acts on -- the
#: background mode -- and never on the positive fraction, which must not be tuned toward
#: cross-donor uniformity.
BAR_MODE_SD = 0.10


def marker_values(marker: str, donors: list) -> dict:
    """One marker's analysis-eligible values, per donor."""
    out = {}
    for d in donors:
        t = pq.read_table(panel.EXPRESSION / f"donor_id={d}" / "data.parquet",
                          columns=["qc_analysis_eligible", marker]).to_pandas()
        out[d] = t.loc[t["qc_analysis_eligible"], marker].to_numpy("float64")
    return out


def derive(markers: list, donors: list, usab: pd.DataFrame, save_figs: bool = True) -> dict:
    reg, norm, _ = UA._import_uniform()
    FIGS.mkdir(parents=True, exist_ok=True)
    flagged = set(usab.loc[usab["acq_ratio"] >= U.BARS["acq_ratio_flag"], "marker"])
    result = {}

    for m in markers:
        t = time.time()
        vals = marker_values(m, donors)
        counts, ranges = UA.marker_histograms(vals)
        before = UA.background_mode_spread(counts, ranges)

        hist = {m: {"counts": [counts[d] for d in donors],
                    "bin_edges": [None] * len(donors)}}
        with contextlib.redirect_stdout(io.StringIO()):     # UniFORM narrates per marker
            shift_map, refs, _ = reg.automatic_registration(
                hist, [m], [m], donors, None, None, UA.NUM_BINS)
        if save_figs:
            plt.gcf().savefig(FIGS / f"{m}.png", bbox_inches="tight", dpi=130)
        plt.close("all")

        bins = dict(zip(donors, shift_map[m]))
        after = UA.background_mode_spread(counts, ranges, bins)
        log_shift = dict(zip(donors, norm.calculate_shift_in_log_pixels(
            {m: ranges}, [m], UA.NUM_BINS, shift_map)[m]))
        scale = {d: float(np.exp(-log_shift[d])) for d in donors}

        # A rigid shift models a scale difference.  It does not model a different antibody
        # clone: that changes specificity and background too, so aligning the modes would
        # glue two different measurements together behind a good-looking statistic.
        refused = m in flagged
        if refused:
            scale = {d: 1.0 for d in donors}
        over_bar = not refused and after["sd_log"] > BAR_MODE_SD

        result[m] = {
            "reference_donor": refs[0],
            "shift_bins": {d: int(bins[d]) for d in donors},
            "shift_log": {d: float(log_shift[d]) for d in donors},
            "scale": scale,
            "provenance": ("refused_reagent_change" if refused else
                           "automatic_registration_over_bar" if over_bar else
                           "automatic_registration"),
            "over_bar": bool(over_bar),
            "acq_ratio": float(usab.loc[usab.marker == m, "acq_ratio"].iloc[0]),
            "mode_sd_before": before["sd_log"], "mode_sd_after": after["sd_log"],
            "mode_fold_before": before["fold"], "mode_fold_after": after["fold"],
            "global_min": ranges["global_min"], "global_max": ranges["global_max"],
            "gmm_cutoff": UA.gmm_cutoff(np.concatenate(list(vals.values())), ranges),
        }
        tag = "REFUSED" if refused else ("ok" if after["sd_log"] <= BAR_MODE_SD else "OVER BAR")
        print(f"  {m:<18} ref={refs[0]:<6} mode fold {before['fold']:8.2f} -> "
              f"{after['fold']:6.2f}   sd {before['sd_log']:.3f} -> {after['sd_log']:.3f}"
              f"   {tag}   {time.time() - t:.0f}s", flush=True)
        del vals
    return result


def apply_shift(shifts: dict, markers: list, donors: list) -> None:
    """Write the normalised value per cell, keyed on object_id."""
    OUT.mkdir(parents=True, exist_ok=True)
    keep = ["donor_id", "object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
            "cell_region", "islet_num"]
    for d in donors:
        df = pq.read_table(panel.EXPRESSION / f"donor_id={d}" / "data.parquet",
                           columns=keep + markers).to_pandas()
        for m in markers:
            df[m] = df[m].to_numpy("float64") * shifts[m]["scale"][d]
        p = OUT / f"donor_id={d}"
        p.mkdir(parents=True, exist_ok=True)
        df.to_parquet(p / "data.parquet", index=False)
        print(f"  {d}: {len(df):,} rows -> {p/'data.parquet'}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--markers", nargs="*", default=None)
    args = ap.parse_args()

    usab = pd.read_parquet(panel.OUT / "marker_donor_usability.parquet")
    per_marker = usab.groupby("marker").agg(
        acq_ratio=("acq_ratio", "first"),
        n_unusable=("verdict", lambda s: (s == "unusable").sum())).reset_index()
    donors = panel.cohort()["donor_id"].tolist()
    markers = args.markers or sorted(per_marker["marker"])

    print(f"Step 2 — UniFORM @ {UA.UNIFORM_COMMIT[:12]}, {len(markers)} markers x "
          f"{len(donors)} donors, {UA.NUM_BINS} bins")
    flagged = per_marker.loc[per_marker.acq_ratio >= U.BARS["acq_ratio_flag"], "marker"].tolist()
    if flagged:
        print(f"  refusing to shift (acquisition-group structure): {', '.join(flagged)}")

    shifts = derive(markers, donors, per_marker)

    over = {m: v["mode_sd_after"] for m, v in shifts.items() if v.get("over_bar")}
    sd = [v["mode_sd_after"] for v in shifts.values()
          if v["provenance"].startswith("automatic_registration")]
    print(f"\nbackground-mode SD after: median {np.median(sd):.4f}, "
          f"{len(sd) - len(over)}/{len(sd)} at or under the {BAR_MODE_SD} bar")
    if over:
        print("  over the bar, carried with the residual recorded (three attempted fixes were")
        print("  measured and refused -- see uniform_adapter.py):")
        for m, v in sorted(over.items(), key=lambda kv: -kv[1]):
            b = shifts[m]["mode_sd_before"]
            print(f"    {m:<18} {b:.3f} -> {v:.3f}   ({b / v:.1f}x better, still over "
                  f"{BAR_MODE_SD})")

    panel.OUT.mkdir(parents=True, exist_ok=True)
    meta = {"uniform_commit": UA.UNIFORM_COMMIT, "num_bins": UA.NUM_BINS,
            "bar_mode_sd": BAR_MODE_SD, "donors": donors, "markers": markers,
            "shifts": shifts}
    (panel.OUT / "shifts.json").write_text(json.dumps(meta, indent=2, sort_keys=True))
    print(f"\nwrote {panel.OUT / 'shifts.json'}")

    if args.dry_run:
        print("dry run: no normalised parquet written")
        return 0
    apply_shift(shifts, markers, donors)
    return 0


if __name__ == "__main__":
    sys.exit(main())
