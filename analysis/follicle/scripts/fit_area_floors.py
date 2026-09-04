#!/usr/bin/env python
"""Derive the QuPath Follicle and PALS area floors instead of asserting them.

data_utils.MIN_QUPATH_FOLLICLE_AREA and MIN_QUPATH_PALS_AREA drop fragments
before the region layer is assembled. They arrived in the codebase as bare
constants with the fitting code absent, so nothing could check them or re-derive
them when the geometry changed. This script is that missing derivation.

Method
------
Pool every QuPath Follicle (or PALS) object across the cohort, take log10 of its
area in µm², and fit a 2-component Gaussian mixture. The floor is the crossover
— the area at which an object is equally likely to belong to the fragment
component and to the real-region component. Areas come from the annotation CSVs
rather than the GeoJSON polygons: the CSV is what `validate_region_layer.py`
checks against, so both ends of the pipeline agree on what an object's area is.

The fit is deterministic (fixed seed, and in practice seed-independent — the two
modes are far enough apart that EM lands in the same place from any start).

This script REPORTS; it does not edit data_utils. A floor is a judgement about
what counts as a real anatomical structure, so a moved number is something to
look at, not something to apply automatically.

    python scripts/fit_area_floors.py

Writes analysis/tables/area_floor_fit.csv.
"""
import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT / "analysis"))
import data_utils as du  # noqa: E402

# What data_utils currently ships, so the report is a comparison rather than a
# bare number.
COMMITTED = {"Follicle": du.MIN_QUPATH_FOLLICLE_AREA,
             "PALS": du.MIN_QUPATH_PALS_AREA}


def load_areas(include_codex=True):
    """Per-class QuPath object areas in µm², pooled over the cohort."""
    frames = []
    af = pd.read_csv(du.DATA_ROOT / "Measurements" / "AnnotationsFinal.csv",
                     usecols=["Image", "Classification", "Area µm^2"],
                     low_memory=False)
    frames.append(af.rename(columns={"Area µm^2": "area"}))

    sh_path = du.DATA_ROOT / "Measurements" / "ForSH2B3.csv"
    if include_codex and sh_path.exists():
        sh = pd.read_csv(sh_path, low_memory=False)
        acol = next(c for c in sh.columns if c.startswith("Area"))
        frames.append(sh[["Image", "Classification", acol]]
                      .rename(columns={acol: "area"}))

    df = pd.concat(frames, ignore_index=True)
    df["Sample"] = df["Image"].apply(du.extract_sample_id)
    df = df[~df["Sample"].isin(du.EXCLUDE_SAMPLES)]
    return df[np.isfinite(df["area"]) & (df["area"] > 0)]


def fit_floor(areas, seed=0):
    """Crossover of a 2-component GMM on log10 area, in µm²."""
    x = np.log10(np.asarray(areas, float)).reshape(-1, 1)
    gmm = GaussianMixture(2, random_state=seed).fit(x)
    lo, hi = np.sort(gmm.means_.ravel())
    grid = np.linspace(lo, hi, 200_001).reshape(-1, 1)
    p = gmm.predict_proba(grid)
    floor = 10 ** grid[int(np.argmin(np.abs(p[:, 0] - p[:, 1]))), 0]
    return floor, 10 ** lo, 10 ** hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", nargs="*", default=["Follicle", "PALS"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--phenocycler-only", action="store_true",
                    help="exclude the ForSH2B3 CODEX sections from the pool")
    args = ap.parse_args()

    df = load_areas(include_codex=not args.phenocycler_only)
    rows = []
    for name in args.classes:
        a = df.loc[df["Classification"] == name, "area"].values
        if len(a) < 50:
            print(f"{name}: only {len(a)} objects — skipped", flush=True)
            continue
        floor, mode_lo, mode_hi = fit_floor(a, seed=args.seed)
        below = a < floor
        committed = COMMITTED.get(name)
        rows.append({
            "Class": name,
            "n_objects": len(a),
            "fragment_mode_um2": mode_lo,
            "region_mode_um2": mode_hi,
            "fitted_floor_um2": floor,
            "committed_floor_um2": committed,
            "ratio_fitted_to_committed": (floor / committed) if committed else np.nan,
            "pct_objects_below_fitted": below.mean() * 100,
            "pct_area_below_fitted": a[below].sum() / a.sum() * 100,
        })
        c_below = a < committed if committed else None
        print(f"{name}: n={len(a):,}  modes {mode_lo:,.0f} / {mode_hi:,.0f} µm²")
        print(f"  fitted floor    {floor:>9,.0f} µm²  "
              f"→ drops {below.mean()*100:5.1f}% of objects, "
              f"{a[below].sum()/a.sum()*100:5.2f}% of area")
        if committed:
            print(f"  committed floor {committed:>9,.0f} µm²  "
                  f"→ drops {c_below.mean()*100:5.1f}% of objects, "
                  f"{a[c_below].sum()/a.sum()*100:5.2f}% of area"
                  f"   ({floor/committed:.2f}× apart)")
        print(flush=True)

    out = pd.DataFrame(rows)
    du.save_table(out, "area_floor_fit")
    return out


if __name__ == "__main__":
    main()
