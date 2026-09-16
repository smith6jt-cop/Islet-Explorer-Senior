#!/usr/bin/env python3
"""Per-islet MEAN MARKER VALUES -- the thing Phase 5 promised and never built.

The Trajectory tab's Multi-Feature Heatmap is meant to show how MARKERS move along pseudotime.
It was showing the eight cell classes instead, and faithfully so: `adata_ins_root.h5ad` is
(7302 islets x 8 cell classes), its `X` IS the composition, and there is no marker anywhere in
the file -- no layer, and `obsm` holds only the ALR composition, the diffmap and the UMAP. The
heatmap read `var_names` and got what was there. `aggregate_islets.py` never computed marker
means, so nothing downstream could have shown them.

This builds them from the per-cell UniFORM-aligned values -- the REDSEA-corrected, aligned
intensities that are the authoritative one-value-per-marker -- averaged within each islet.

TWO CHOICES THAT ARE NOT ARBITRARY:

  * ELIGIBLE CELLS ONLY. `qc_analysis_eligible` excludes detections with no nucleus, which are
    not known to be cells at all, and averaging their intensities into an islet's mean would
    import measurements from objects the rest of the analysis refuses to count.

  * CORE AND PERI ARE SEPARATE COLUMNS, never pooled. An islet's core mean and the mean of
    the 20 um ring around it are different quantities -- that is the entire point of the region
    split -- and a pooled mean would move with how much ring a given islet happens to have.

A marker absent from a donor's panel yields NaN, not zero: unavailable is not absence.

    python core/phenocycler/scripts/typing/build_islet_marker_means.py [--donors 115]
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time

import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

OUT = S.OUT / "app"
META = {"donor_id", "object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
        "cell_region", "islet_num"}


def donor_means(donor: str) -> pd.DataFrame | None:
    path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    if not path.exists():
        return None
    table = pq.read_table(path)
    df = table.to_pandas()
    markers = [c for c in df.columns if c not in META]

    df = df[df["qc_analysis_eligible"].fillna(False).astype(bool)]
    # THE ALIGNED TABLE SAYS "peri", NOT "expansion". Checked, not assumed: the first version
    # filtered on ["core", "expansion"], matched nothing for the ring, and silently produced a
    # core-only table -- 61 columns where 116 were expected, with no error anywhere.
    df = df[df["islet_num"].notna() & df["cell_region"].isin(["core", "peri"])]
    if df.empty:
        return None

    df["islet_key"] = "Islet_" + df["islet_num"].astype("Int64").astype(str)
    out = []
    for region, suffix in (("core", ""), ("peri", "_peri")):
        part = df[df["cell_region"] == region]
        if part.empty:
            continue
        grouped = part.groupby("islet_key")[markers].mean()
        grouped.columns = [f"marker_mean{suffix}_{m}" for m in markers]
        out.append(grouped)
    if not out:
        return None
    if len(out) != 2:
        raise RuntimeError(
            f"{donor}: expected core AND peri means, built {len(out)}. A region label that "
            f"matches nothing yields a silently half-built table, which is how the first run "
            f"produced core-only columns without erroring.")
    merged = pd.concat(out, axis=1)
    merged.insert(0, "combined_islet_id", [f"{donor}_{k}" for k in merged.index])
    merged.insert(1, "case_id", donor)
    return merged.reset_index()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--include-marginal", action="store_true")
    args = ap.parse_args()

    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    frames = []
    for donor in donors:
        started = time.time()
        frame = donor_means(donor)
        if frame is None:
            print(f"  {donor}: no aligned per-cell values"); continue
        frames.append(frame)
        n_core = sum(c.startswith("marker_mean_") for c in frame.columns)
        print(f"  {donor}: {len(frame):>6,} islets | {n_core} markers x core+expansion | "
              f"{time.time()-started:.0f}s", flush=True)
    if not frames:
        return 1
    everything = pd.concat(frames, ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "islet_marker_means.parquet"
    everything.to_parquet(target, index=False)
    print(f"\n{len(everything):,} islets, {everything.shape[1]} columns -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
