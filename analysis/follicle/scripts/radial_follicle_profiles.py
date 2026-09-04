#!/usr/bin/env python
"""Radial marker composition around follicles, in 30 µm bands.

Walks outward from the follicles in fixed-width bands and records the
composition of every marker common to all 14 donors. The perifollicular zone is
read off where that composition stops looking like follicle and settles to red
pulp, rather than being asserted as a fixed annulus.

Distance is measured to the follicle *boundary*, signed negative inside. That is
the radial coordinate generalised to non-circular follicles: bands from a
centroid would sit at different radii relative to the edge for every follicle,
smearing the transition when pooled. At distance 0 every follicle is aligned.

Marker values are expressed as fold-change over each donor's own red-pulp
baseline (cells beyond RED_PULP_FROM µm of any follicle), which makes donors
comparable despite differing dynamic range and the clipped-zero floor.

The perifollicular width is then read off the SMA arm. SMA is suppressed inside
the follicle, peaks just outside it, and decays back to red-pulp baseline; the
width is where that peak has fallen halfway back. Because the profiles are
fold-change over red pulp, baseline is 1.0 by construction, so the half-max
level is simply 1 + (peak - 1) / 2.

Outputs
-------
analysis/tables/radial_profiles.csv        donor x band x marker
analysis/tables/radial_pfz_width.csv       derived perifollicular width
                                           (consumed by build_region_layer_qupath.py)

The P2_F11 figure is drawn by plot_region_layer_results.py, which reads both
tables back.

    python scripts/radial_follicle_profiles.py
"""
import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT / "analysis"))
import data_utils as du  # noqa: E402

BAND = 30                 # µm
D_MIN, D_MAX = -300, 600  # µm, relative to the follicle boundary
RED_PULP_FROM = 400       # µm — baseline is taken beyond this
BOUNDARY_STEP = 8         # µm between sampled boundary points
CD20_BASELINE = 1.2       # fold over red pulp; where the follicle signal ends


def signed_distance(coords, polys):
    """Distance to the nearest follicle boundary, negative inside."""
    return du.signed_distance_to_regions(coords, polys,
                                         boundary_step=BOUNDARY_STEP)


def pfz_width(prof):
    """Per-donor perifollicular width, read off the SMA arm of the profiles.

    Only bands outside the follicle count. SMA peaks there; the width is the
    first band at or below half its height above the red-pulp baseline. A donor
    with no SMA peak — HDL070 and HDL073, whose maxima land in the far tail —
    gets a width the caller's `<= 200` guard rejects, and falls back to the
    cohort median. CD20_baseline_um is reported alongside as the independent
    check: the follicle signal should be gone by roughly the same distance.
    """
    rows = []
    for sample, g in prof.groupby("Sample", observed=True):
        out = g[g["distance_um"] > 0].sort_values("distance_um")
        if out.empty or "SMA" not in out or out["SMA"].isna().all():
            continue
        i = int(np.nanargmax(out["SMA"].values))
        peak_um = float(out["distance_um"].values[i])
        peak = float(out["SMA"].values[i])

        half = 1.0 + (peak - 1.0) / 2
        beyond = out[(out["distance_um"] > peak_um) & (out["SMA"] <= half)]
        width = float(beyond["distance_um"].iloc[0]) if len(beyond) else np.nan

        cd20 = out[out["CD20"] <= CD20_BASELINE] if "CD20" in out else out.iloc[:0]
        rows.append({"Sample": sample,
                     "Platform": out["Platform"].iloc[0],
                     "SMA_peak_um": peak_um,
                     "SMA_peak": peak,
                     "PFZ_width_um": width,
                     "CD20_baseline_um": (float(cd20["distance_um"].iloc[0])
                                          if len(cd20) else np.nan)})
    return pd.DataFrame(rows).sort_values(["PFZ_width_um", "Sample"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", nargs="*", default=list(du.MARKERS_USABLE_ALL_DONORS))
    args = ap.parse_args()
    markers = [m for m in args.markers if m != "DAPI"]

    layer = du.load_region_layer()
    by_sample = {du.extract_sample_id(k): v for k, v in layer.items()}

    print(f"loading cells for {len(markers)} markers …", flush=True)
    cells = du.load_single_cells(markers=markers)
    print(f"  {len(cells):,} cells, {cells['Sample'].nunique()} donors", flush=True)

    edges = np.arange(D_MIN, D_MAX + BAND, BAND)
    centres = edges[:-1] + BAND / 2
    rows = []

    for sample, g in cells.groupby("Sample", observed=True):
        fol = by_sample.get(sample, {}).get("Follicle", [])
        if not fol:
            print(f"  {sample:<14} no follicles — skipped", flush=True)
            continue
        coords = g[["x", "y"]].values.astype(float)
        d = signed_distance(coords, fol)

        base = d > RED_PULP_FROM
        if base.sum() < 500:
            base = d > np.nanpercentile(d, 90)
        idx = np.digitize(d, edges) - 1

        for b, centre in enumerate(centres):
            m = idx == b
            if m.sum() < 30:
                continue
            row = {"Sample": sample, "Platform": g["Platform"].iloc[0],
                   "Genotype": str(g["Genotype"].iloc[0]),
                   "distance_um": centre, "n_cells": int(m.sum())}
            for mk in markers:
                v = g[mk].values.astype(float)
                ref = np.nanmean(v[base])
                row[mk] = (np.nanmean(v[m]) / ref) if ref and ref > 0 else np.nan
            rows.append(row)
        print(f"  {sample:<14} {len(fol):4d} follicles, "
              f"{int(base.sum()):>9,} baseline cells", flush=True)

    prof = pd.DataFrame(rows)
    du.save_table(prof, "radial_profiles")

    width = pfz_width(prof)
    du.save_table(width, "radial_pfz_width")
    usable = width.loc[width["PFZ_width_um"] <= 200, "PFZ_width_um"]
    print(f"\nperifollicular width: {len(usable)}/{len(width)} donors measured, "
          f"median {usable.median():.0f} µm "
          f"(the rest fall back to it)", flush=True)
    return prof, width


if __name__ == "__main__":
    main()
