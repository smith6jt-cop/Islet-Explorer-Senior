"""Step 2: derive and apply the UniFORM shift for every marker, then check it.

Two passes.  The first walks markers, because a shift is a per-marker object estimated from
all selected donors at once; the second walks donors, because that is how output is written.

    python scripts/celltype/run_step2.py             # from core/phenocycler/
    python scripts/celltype/run_step2.py --dry-run   # derive and check, write no parquet
    python scripts/celltype/run_step2.py --out ../../data/development/celltype/uniform22 \
        --shifts-out ../../data/development/celltype/shifts_22donor.json
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import pathlib
import sys
import time
import warnings

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


def marker_values(marker: str, donors: list, unavailable: dict | None = None) -> dict:
    """One marker's analysis-eligible values, for the donors that carry it.

    A panel-variant marker is aligned across its own donors -- CD38 over 13, b_Catenin1 over
    4 -- and left alone on the rest.  Aligning it against donors that never acquired it would
    be aligning against nothing. Acquired columns without any finite QC-eligible value
    at or above UniFORM's log floor are unavailable, not a synthetic zero histogram.
    """
    out = {}
    for d in [x for x in donors if marker in set(panel.donor_markers(x))]:
        t = pq.ParquetFile(panel.expression_path(d)).read(
            columns=["qc_analysis_eligible", marker]).to_pandas()
        values = t.loc[t["qc_analysis_eligible"], marker].to_numpy("float64")
        finite = np.isfinite(values)
        if not np.any(finite & (values >= 1.0)):
            reason = (f"No finite values among {len(values)} QC-eligible cells." if not finite.any()
                      else "No finite QC-eligible values at or above UniFORM's log floor (1.0).")
            if unavailable is not None:
                unavailable[d] = reason
            else:
                warnings.warn(f"{marker}, donor {d}: unavailable: {reason}",
                              RuntimeWarning, stacklevel=2)
            continue
        out[d] = values
    return out


def derive(markers: list, donors: list, usab: pd.DataFrame, save_figs: bool = True,
           figs: pathlib.Path = FIGS) -> dict:
    reg, norm, _ = UA._import_uniform()
    if save_figs:
        figs.mkdir(parents=True, exist_ok=True)
    flagged = set(usab.loc[usab["acq_ratio"] >= U.BARS["acq_ratio_flag"], "marker"])
    result = {}

    for m in markers:
        t = time.time()
        acq_ratio = float(usab.loc[usab.marker == m, "acq_ratio"].iloc[0])
        unavailable = {}
        vals = marker_values(m, donors, unavailable)
        for d, reason in unavailable.items():
            print(f"  {m}, donor {d}: unavailable: {reason}", flush=True)
        have = sorted(vals)
        if len(have) < 2:
            print(f"  {m:<18} measurable on {len(have)} donor(s) — no alignment possible")
            result[m] = {"reference_donor": have[0] if have else "unavailable",
                         "shift_bins": {}, "shift_log": {},
                         "scale": {d: "unavailable" if d in unavailable else 1.0 for d in donors},
                         "provenance": "too_few_donors_to_align", "over_bar": False,
                         "acq_ratio": acq_ratio, "n_donors": len(have),
                         "n_acquired_donors": len(have) + len(unavailable),
                         "unavailable_donors": unavailable,
                         "unavailable_reasons": ({} if have else {
                             "reference_donor": "No donor has a measurable histogram."}),
                         "mode_sd_before": float("nan"), "mode_sd_after": float("nan"),
                         "mode_fold_before": float("nan"), "mode_fold_after": float("nan"),
                         "global_min": float("nan"), "global_max": float("nan"),
                         "gmm_cutoff": float("nan")}
            continue
        counts, ranges = UA.marker_histograms(vals)
        before = UA.background_mode_spread(counts, ranges)

        hist = {m: {"counts": [counts[d] for d in have],
                    "bin_edges": [None] * len(have)}}
        with contextlib.redirect_stdout(io.StringIO()):     # UniFORM narrates per marker
            shift_map, refs, _ = reg.automatic_registration(
                hist, [m], [m], have, None, None, UA.NUM_BINS)
        if save_figs:
            plt.gcf().savefig(figs / f"{m}.png", bbox_inches="tight", dpi=130)
        plt.close("all")

        bins = dict(zip(have, shift_map[m]))
        after = UA.background_mode_spread(counts, ranges, bins)
        log_shift = dict(zip(have, norm.calculate_shift_in_log_pixels(
            {m: ranges}, [m], UA.NUM_BINS, shift_map)[m]))
        scale = {d: float(np.exp(-log_shift[d])) if d in log_shift else 1.0 for d in donors}

        # A rigid shift models a scale difference.  It does not model a different antibody
        # clone: that changes specificity and background too, so aligning the modes would
        # glue two different measurements together behind a good-looking statistic.
        refused = m in flagged
        if refused:
            scale = {d: 1.0 for d in donors}
        for d in unavailable:
            scale[d] = "unavailable"
        over_bar = not refused and after["sd_log"] > BAR_MODE_SD

        result[m] = {
            "reference_donor": refs[0],
            # Only the donors that carry the marker have a shift; the rest carry scale 1.0
            # and are listed in `scale` so every donor is accounted for either way.
            "shift_bins": {d: int(bins[d]) for d in have},
            "shift_log": {d: float(log_shift[d]) for d in have},
            "scale": scale,
            "provenance": ("refused_reagent_change" if refused else
                           "automatic_registration_over_bar" if over_bar else
                           "automatic_registration"),
            "over_bar": bool(over_bar),
            "acq_ratio": acq_ratio,
            "n_donors": len(have),
            "n_acquired_donors": len(have) + len(unavailable),
            "unavailable_donors": unavailable,
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


def apply_shift(shifts: dict, markers: list, donors: list, out: pathlib.Path = OUT) -> None:
    """Write the normalised value per cell, keyed on object_id."""
    out.mkdir(parents=True, exist_ok=True)
    keep = ["donor_id", "object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
            "cell_region", "islet_num"]
    for d in donors:
        got = [m for m in markers if m in set(panel.donor_markers(d))]
        df = pq.ParquetFile(panel.expression_path(d)).read(columns=keep + got).to_pandas()
        df["donor_id"] = d               # match the cohort and output partition (0112 == 112)
        for m in markers:
            if m not in df.columns:
                continue                       # this donor never acquired it
            if d in shifts[m].get("unavailable_donors", {}):
                continue                       # retain missing/raw input; no estimated scale
            df[m] = df[m].to_numpy("float64") * shifts[m]["scale"][d]
        p = out / f"donor_id={d}"
        p.mkdir(parents=True, exist_ok=True)
        df.to_parquet(p / "data.parquet", index=False)
        print(f"  {d}: {len(df):,} rows -> {p/'data.parquet'}", flush=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--markers", nargs="*", default=None)
    ap.add_argument("--donors", nargs="+", default=None,
                    help="Eligible donor IDs; defaults to the full eligible cohort. "
                         "Leading-zero aliases are accepted and deduplicated.")
    ap.add_argument("--out", type=pathlib.Path, default=OUT,
                    help="Normalised parquet directory (default: %(default)s). "
                         "Custom output also isolates figures and the default shift record.")
    ap.add_argument("--shifts-out", type=pathlib.Path, default=None,
                    help="Shift JSON path (default: shifts.json beside the historical uniform/ "
                         "output, or inside a custom --out directory).")
    args = ap.parse_args(argv)

    cohort = panel.cohort()["donor_id"].tolist()
    donors = cohort
    if args.donors is not None:
        donors = list(dict.fromkeys(str(int(d)) if d.isdecimal() else d for d in args.donors))
        unknown = sorted(set(donors) - set(cohort))
        if unknown:
            ap.error(f"Donors are not in the eligible cohort: {', '.join(unknown)}")
    historical = args.out.resolve() == OUT.resolve()
    shifts_out = args.shifts_out or (panel.OUT / "shifts.json" if historical
                                    else args.out / "shifts.json")
    figs = FIGS if historical else args.out / "figures_step2"

    usab = pd.read_parquet(panel.OUT / "marker_donor_usability.parquet")
    per_marker = usab.groupby("marker").agg(acq_ratio=("acq_ratio", "first")).reset_index()
    markers = args.markers or sorted(per_marker["marker"])
    markers = [m for m in markers if any(m in set(panel.donor_markers(d)) for d in donors)]
    # Resolve before importing UniFORM, whose module suppresses Python warnings globally.
    expression_inputs = {d: str(panel.expression_path(d)) for d in donors}

    print(f"Step 2 — UniFORM @ {UA.UNIFORM_COMMIT[:12]}, {len(markers)} markers x "
          f"{len(donors)} donors, {UA.NUM_BINS} bins")
    flagged = per_marker.loc[per_marker.acq_ratio >= U.BARS["acq_ratio_flag"], "marker"].tolist()
    if flagged:
        print(f"  refusing to shift (acquisition-group structure): {', '.join(flagged)}")

    shifts = derive(markers, donors, per_marker, figs=figs)

    over = {m: v["mode_sd_after"] for m, v in shifts.items() if v.get("over_bar")}
    sd = [v["mode_sd_after"] for v in shifts.values()
          if v["provenance"].startswith("automatic_registration")]
    median_sd = f"{np.median(sd):.4f}" if sd else "unavailable (no automatically aligned markers)"
    print(f"\nbackground-mode SD after: median {median_sd}, "
          f"{len(sd) - len(over)}/{len(sd)} at or under the {BAR_MODE_SD} bar")
    if over:
        print("  over the bar, carried with the residual recorded (three attempted fixes were")
        print("  measured and refused -- see uniform_adapter.py):")
        for m, v in sorted(over.items(), key=lambda kv: -kv[1]):
            b = shifts[m]["mode_sd_before"]
            print(f"    {m:<18} {b:.3f} -> {v:.3f}   ({b / v:.1f}x better, still over "
                  f"{BAR_MODE_SD})")

    shifts_out.parent.mkdir(parents=True, exist_ok=True)
    recorded_shifts = {}
    for m, shift in shifts.items():
        record = dict(shift)
        unavailable = dict(shift.get("unavailable_reasons", {}))
        for field, value in shift.items():
            if isinstance(value, float) and not np.isfinite(value):
                record[field] = "unavailable"
                unavailable[field] = (
                    "No finite acquisition-group ratio in the Step-1 usability input; "
                    "Step 2 does not recompute that measurement." if field == "acq_ratio" else
                    "Fewer than two donors have measurable histograms." if shift["n_donors"] < 2 else
                    "The estimator returned no finite measurement."
                )
        if unavailable:
            record["unavailable_reasons"] = unavailable
        recorded_shifts[m] = record
    meta = {"uniform_commit": UA.UNIFORM_COMMIT, "num_bins": UA.NUM_BINS,
            "bar_mode_sd": BAR_MODE_SD, "donors": donors, "markers": markers,
            "expression_inputs": expression_inputs, "shifts": recorded_shifts}
    shifts_out.write_text(json.dumps(meta, indent=2, sort_keys=True, allow_nan=False))
    print(f"\nwrote {shifts_out}")

    if args.dry_run:
        print("dry run: no normalised parquet written")
        return 0
    apply_shift(shifts, markers, donors, out=args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
