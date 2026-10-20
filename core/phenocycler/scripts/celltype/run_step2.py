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
import functools
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
import alignment_block as AB             # noqa: E402

OUT = panel.OUT / "uniform"
FIGS = panel.OUT / "figures_step2"

#: Acceptance, declared before the run.  Stated on the quantity the shift acts on -- the
#: background mode -- and never on the positive fraction, which must not be tuned toward
#: cross-donor uniformity.
BAR_MODE_SD = 0.10

#: WHEN A MARKER IS LEFT UNALIGNED (adjusted 2026-10-02).
#:
#: Until then a marker whose median steps with the acquisition group (Step 1's ratio at or over
#: ``usability.BARS['acq_ratio_flag']``) was refused outright, "however good its shift statistic
#: looks", on the reading that such a step is a reagent change.  With all 25 donors in one cohort
#: that rule refused CD8 for everyone -- and UniFORM's own registration puts CD8's background
#: modes within 0.065 of each other, inside the bar above.  The maintainer, 2026-10-02: "UniFORM
#: should be able to handle the different pixel intensity distributions, check if there are
#: levers that need to be adjusted."  This rule was the lever: it is the lane's, not UniFORM's,
#: and what it read as a reagent change is not something the files can establish (the scan
#: files' clone entries are unreliable -- maintainer, same day).
#:
#: So the refusal now needs BOTH: the acquisition-group ratio at or over its flag AND UniFORM's
#: shift leaving the background modes farther apart than ``BAR_MODE_SD``.  A flagged marker that
#: UniFORM aligns within the bar is aligned like any other (and recorded ``acq_ratio_flagged``).
#: CD57 is still refused under it: 0.345 after the shift, and under the floor of 1.0 on up to
#: 91 % of a donor's cells.  The label written for a refusal is unchanged because
#: ``typing/refine.py`` and ``kintsugi/features.py`` read it; it names the rule, not a finding.
REFUSED = "refused_reagent_change"


@functools.lru_cache(maxsize=None)
def _sample_rows(donor: str, path: str) -> np.ndarray:
    """Per row of a donor's expression file: analysis-eligible AND in the pancreas block.

    The sample UniFORM aligns is the BLOCK (maintainer, 2026-10-02: "UniFORM should be done on a
    same-block basis"): tissue plus any same-block lymph node, never a different-block node
    that shares the slide.  Keyed on the file as well as the donor, and read once: the mask is
    the same for every marker.
    """
    t = pq.ParquetFile(path).read(
        columns=["object_id", "qc_analysis_eligible", "cell_region"]).to_pandas()
    return (t["qc_analysis_eligible"].to_numpy(bool)
            & AB.population(donor, t["object_id"], t["cell_region"]))


def marker_values(marker: str, donors: list, unavailable: dict | None = None) -> dict:
    """One marker's analysis-eligible values in each donor's pancreas block, for the donors that carry it.

    A panel-variant marker is aligned across its own donors -- CD38 over 15, b_Catenin1 over
    6 (2 of those donors carry both) -- and left alone on the rest.  Aligning it against
    donors that never acquired it would be aligning against nothing. Acquired columns without
    any finite QC-eligible value at or above UniFORM's log floor are unavailable, not a
    synthetic zero histogram.
    """
    out = {}
    for d in [x for x in donors if marker in set(panel.donor_markers(x))]:
        path = panel.expression_path(d)
        values = pq.ParquetFile(path).read(columns=[marker]).to_pandas()[marker] \
            .to_numpy("float64")[_sample_rows(str(d), str(path))]
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
           figs: pathlib.Path = FIGS, references: dict | None = None) -> dict:
    """The shift per marker and donor.

    ``references`` pins a marker's reference donor (UniFORM's own ``reference_samples``): the
    donor every other donor is shifted onto, and so the unit the aligned values are in.  Left
    to itself UniFORM picks the histogram closest to the cohort's mean one, which changes when
    the cohort does -- and a floor or a top measured in one reference donor's unit says nothing
    about values aligned onto another's.  A marker without a pinned donor, or whose pinned donor
    carries no measurable histogram, falls back to UniFORM's choice and says so.
    """
    reg, norm, _ = UA._import_uniform()
    references = references or {}
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
        pinned = references.get(m)
        if pinned is not None and pinned not in have:
            print(f"  {m}: pinned reference donor {pinned} has no measurable histogram here -- "
                  "UniFORM chooses", flush=True)
            pinned = None
        with contextlib.redirect_stdout(io.StringIO()):     # UniFORM narrates per marker
            shift_map, refs, _ = reg.automatic_registration(
                hist, [m], [m], have, [pinned] if pinned is not None else None, None, UA.NUM_BINS)
        if save_figs:
            plt.gcf().savefig(figs / f"{m}.png", bbox_inches="tight", dpi=130)
        plt.close("all")

        bins = dict(zip(have, shift_map[m]))
        after = UA.background_mode_spread(counts, ranges, bins)
        log_shift = dict(zip(have, norm.calculate_shift_in_log_pixels(
            {m: ranges}, [m], UA.NUM_BINS, shift_map)[m]))
        scale = {d: float(np.exp(-log_shift[d])) if d in log_shift else 1.0 for d in donors}

        # A marker is left unaligned only when its level steps with the acquisition group AND
        # UniFORM's shift does not bring the background modes within the bar (see REFUSED).
        refused = m in flagged and after["sd_log"] > BAR_MODE_SD
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
            "provenance": (REFUSED if refused else
                           "automatic_registration_over_bar" if over_bar else
                           "automatic_registration"),
            "reference_source": "pinned" if pinned is not None else "uniform_choice",
            "over_bar": bool(over_bar),
            "acq_ratio": acq_ratio,
            "acq_ratio_flagged": bool(m in flagged),
            "n_donors": len(have),
            "n_acquired_donors": len(have) + len(unavailable),
            "unavailable_donors": unavailable,
            "mode_sd_before": before["sd_log"], "mode_sd_after": after["sd_log"],
            "mode_fold_before": before["fold"], "mode_fold_after": after["fold"],
            "global_min": ranges["global_min"], "global_max": ranges["global_max"],
            "gmm_cutoff": UA.gmm_cutoff(np.concatenate(list(vals.values())), ranges),
        }
        tag = "REFUSED" if refused else ("ok" if after["sd_log"] <= BAR_MODE_SD else "OVER BAR")
        if m in flagged and not refused:
            tag += " (group ratio over its flag; UniFORM within the bar -- aligned)"
        print(f"  {m:<18} ref={refs[0]:<6} mode fold {before['fold']:8.2f} -> "
              f"{after['fold']:6.2f}   sd {before['sd_log']:.3f} -> {after['sd_log']:.3f}"
              f"   {tag}   {time.time() - t:.0f}s", flush=True)
        del vals
    return result


def apply_shift(shifts: dict, markers: list, donors: list, out: pathlib.Path = OUT) -> None:
    """Write the normalised value per cell, keyed on object_id.

    A cell of another block on the same slide (a different-block lymph node) was not part of the
    sample the shift was estimated from, so it has NO aligned value here: every marker is written
    unavailable (NaN) for it and ``in_alignment_block`` is False.  The row stays, so the table
    still has one row per cell of the expression stage; what it must never carry is that cell's
    raw value standing in a column of aligned ones.
    """
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
        inside = AB.population(d, df["object_id"], df["cell_region"])
        if not inside.all():
            df.loc[~inside, got] = np.nan
        df["in_alignment_block"] = inside
        p = out / f"donor_id={d}"
        p.mkdir(parents=True, exist_ok=True)
        df.to_parquet(p / "data.parquet", index=False)
        print(f"  {d}: {len(df):,} rows ({int(inside.sum()):,} in the pancreas block) -> "
              f"{p/'data.parquet'}", flush=True)


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
    ap.add_argument("--references-from", type=pathlib.Path, default=None,
                    help="A shifts.json of an earlier run: each marker keeps the reference donor "
                         "recorded there, so the aligned values stay in that run's unit (and a "
                         "scale measured on it still applies). Markers it does not carry, or whose "
                         "donor is not in this cohort, are left to UniFORM's choice.")
    args = ap.parse_args(argv)

    cohort = panel.cohort()["donor_id"].tolist()
    donors = cohort
    if args.donors is not None:
        def identity(value: str) -> str:
            return str(int(value)) if value.isdecimal() else value
        # Keep the cohort's authoritative spelling (e.g. `0112`) after resolving numeric
        # aliases. Normalizing both sides to `112` used to reject the documented leading-zero
        # alias because the roster deliberately preserves `0112`.
        cohort_spelling = {identity(str(d)): str(d) for d in cohort}
        donors = list(dict.fromkeys(cohort_spelling.get(identity(str(d)), str(d))
                                    for d in args.donors))
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
        print(f"  level steps with the acquisition group: {', '.join(flagged)} -- left unaligned "
              f"only if UniFORM's shift also leaves the modes over the {BAR_MODE_SD} bar")

    references = {}
    if args.references_from is not None:
        canon = {str(int(d)) if d.isdecimal() else d: d for d in donors}       # 112 == 0112
        earlier = json.loads(args.references_from.read_text())["shifts"]
        for m in markers:
            ref = str(earlier.get(m, {}).get("reference_donor", ""))
            ref = canon.get(str(int(ref)) if ref.isdecimal() else ref)
            if ref is not None:
                references[m] = ref
        print(f"  reference donor pinned from {args.references_from} for {len(references)} of "
              f"{len(markers)} markers; UniFORM chooses for: "
              f"{', '.join(m for m in markers if m not in references) or 'none'}")

    shifts = (derive(markers, donors, per_marker, figs=figs, references=references) if references
              else derive(markers, donors, per_marker, figs=figs))

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
            "references_from": str(args.references_from) if args.references_from else None,
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
