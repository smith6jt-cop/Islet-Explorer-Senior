#!/usr/bin/env python3
"""The galectin layer -- correction J.

Three classes per islet and per expansion, all of them QuPath's own object-classifier verdicts read
from the annotation export:

    galectin-positive        `Num M2Galectin3`
    galectin+ macrophage     `Num Iba1_CD68: M2Galectin3`
    metabolically-active     `Num M2Galectin3` - `Num Iba1_CD68: M2Galectin3`
                             -- galectin without macrophage markers, the maintainer's definition
    macrophage               `Num Iba1_CD68`          (the denominator for the split above)

**Why the counts are read and not recomputed, against what the plan says.**  The plan calls for
recomputing per cell from the per-donor-normalised object classifications, because correction A
found the batch-2 annotation counts inflated ~2.5x by a duplicated detection pass.  That reason
does not apply here: **all ten clean-cohort donors are batch 1**, where `Num Detections` matches the
post-QC eligible count within 0.1% (115: 59,577 vs 59,541).  And recomputation is not available
anyway -- the classifiers are OpenCV `RTrees` over `Cell: M2Gal3: Mean`, `Max` and `Std.Dev.`, and
**the export carries only the mean**, so two of the three input features do not exist outside
QuPath.  Reproducing them would mean substituting a different method and calling it the same one.
Track B, which adds the batch-2 donors, needs either a corrected annotation export or those two
measurement columns; it cannot be solved on this side.

**Two denominators, both reported, because they are not the same.**  `n_cells_qupath` is
`Num Detections`, which QuPath counts without the nucleus rule; `n_cells_eligible` is this lane's
denominator, which excludes a detection with no nucleus. `pos_frac` uses QuPath's own, so numerator
and denominator come from one source; the eligible count sits beside it so the gap is visible.

    python core/phenocycler/scripts/typing/build_galectin_layer.py
"""
from __future__ import annotations

import argparse
import re
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

APP = S.OUT / "app"
ANNOTATIONS = S.REPO.parent / "results/AnnotationsMeasurements.csv"

#: annotation column -> the class name this build uses.
COUNT_COLUMNS = {
    "Num M2Galectin3": "Galectin_positive",
    "Num Iba1_CD68: M2Galectin3": "Galectin_macrophage",
    "Num Iba1_CD68": "Macrophage",
}
#: The app's word for the 20 um zone, from the annotation's own name.
REGION_OF = {"": "core", "_exp20um": "peri"}


def parse_annotations(donors: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(ANNOTATIONS)
    frame["case_id"] = frame["Image"].str.extract(r"(\d{3,4})")[0]
    frame = frame[frame["case_id"].isin(donors)].copy()

    parsed = frame["Name"].str.extract(r"^(Islet_\d+)(_exp20um)?$")
    frame["islet_key"] = parsed[0]
    frame["region"] = parsed[1].fillna("").map(REGION_OF)
    frame = frame[frame["islet_key"].notna() & frame["region"].notna()]
    return frame


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    missing = [c for c in COUNT_COLUMNS if c not in pd.read_csv(ANNOTATIONS, nrows=1).columns]
    if missing:
        raise KeyError(f"the annotation export lacks {missing}")

    frame = parse_annotations(donors)
    print(f"{len(frame):,} annotation rows over {frame['case_id'].nunique()} donors "
          f"({(frame['region'] == 'core').sum():,} core, {(frame['region'] == 'peri').sum():,} peri)")

    out = frame[["case_id", "islet_key", "region"]].copy()
    out["n_cells_qupath"] = frame["Num Detections"]
    for column, label in COUNT_COLUMNS.items():
        out[f"count_{label}"] = frame[column].fillna(0).astype(int)
    # The maintainer's definition: galectin without macrophage markers.
    out["count_Metabolically_active"] = (out["count_Galectin_positive"]
                                         - out["count_Galectin_macrophage"]).clip(lower=0)

    # This lane's own denominator beside QuPath's, so the difference is visible rather than assumed.
    eligible = []
    for donor in donors:
        cells = pq.read_table(S.OUT / f"cells/donor_id={donor}/data.parquet",
                              columns=["cell_class", "cell_region_20um", "islet_key"]).to_pandas()
        cells = cells[(cells["cell_class"] != "Excluded") & cells["islet_key"].notna()]
        cells["region"] = cells["cell_region_20um"].map({"core": "core", "expansion": "peri"})
        counted = (cells.dropna(subset=["region"]).groupby(["islet_key", "region"])
                   .size().rename("n_cells_eligible").reset_index())
        counted["case_id"] = donor
        eligible.append(counted)
    out = out.merge(pd.concat(eligible, ignore_index=True),
                    on=["case_id", "islet_key", "region"], how="left")

    for label in list(COUNT_COLUMNS.values()) + ["Metabolically_active"]:
        out[f"prop_{label}"] = np.where(out["n_cells_qupath"] > 0,
                                        out[f"count_{label}"] / out["n_cells_qupath"], np.nan)

    target = APP / "islet_galectin.parquet"
    out.to_parquet(target, index=False)

    gap = (out["n_cells_qupath"] - out["n_cells_eligible"]).dropna()
    print(f"\nQuPath vs eligible cell count: median difference {gap.median():+.0f} "
          f"({100 * (gap / out['n_cells_qupath'].replace(0, np.nan)).median():+.1f}%)")
    print("\nmedian proportion per class, by region:")
    for region in ("core", "peri"):
        here = out[out["region"] == region]
        stats = {label: here[f"prop_{label}"].median()
                 for label in list(COUNT_COLUMNS.values()) + ["Metabolically_active"]}
        print(f"  {region:>5}: " + "  ".join(f"{k} {v:.4f}" for k, v in stats.items()))
    print(f"\n{len(out):,} islet-regions -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
