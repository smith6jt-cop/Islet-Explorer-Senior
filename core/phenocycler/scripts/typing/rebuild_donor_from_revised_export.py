#!/usr/bin/env python3
"""Rebuild one donor's per-cell table from a revised QuPath export, REUSING REDSEA.

WHY THIS EXISTS INSTEAD OF A RE-INGEST. The pipeline is content-addressed: `run_id` is a hash
over the whole cohort manifest, so repointing three donors at their revised exports re-keys the
run and forces all twenty donors through ingest, geometry, REDSEA and expression again -- about
2.3 hours to recompute byte-identical results for seventeen donors that did not change. That
contract earns its keep for a release. It does not earn it while the lane is still being built.

THE KEY FACT THAT MAKES REUSE CORRECT: the revised export did not re-segment anything. QuPath
scripts 0b and 1-4 re-ran against a NEW TISSUE ANNOTATION -- 0b removed the duplicate detections,
1_panc_objects re-derived the Islet annotations inside the new tissue boundary, and 3_signed_distance
recomputed the distances. Cell DETECTIONS keep their UUIDs from the original detection run (measured:
99.5-99.8% of revised object ids are already in the run, and the cell outline geojson is a SUPERSET
of every revised export). A cell that survives de-duplication is the same object, with the same
polygon, so its REDSEA-corrected and UniFORM-aligned marker values are unchanged. They are joined
here by `object_id` -- never centroid, never row order.

WHAT IS RECOMPUTED, because the new tissue annotation changed it:
  parent / islet_num     from the export's `Parent`
  cell_region            from the export's re-run signed distance, at the 20 um rule
  qc_analysis_eligible   from the geometry thresholds in config.ini
  qc_estimation_eligible the same, plus the no-cytoplasm rule, which gates ESTIMATION only

WHAT CANNOT BE RECOMPUTED HERE, and is therefore dropped from the eligibility rule: the three
raster checks (`raster_area_mismatch`, `raster_empty`, `nucleus_raster_empty`) need the whole-slide
label raster. On 6534 they excluded 51,517 cells in the old run, and `raster_area_mismatch` was
overwhelmingly the DUPLICATE artefact -- the loser of each pair lost every pixel to the survivor.
With the duplicates gone at source that reason should be largely empty, but "should be" is not
"measured", so the column records which rules were actually applied.

    python core/phenocycler/scripts/typing/rebuild_donor_from_revised_export.py 6505 6534 6623
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import subpanel_io as S

REVISED = S.REPO / "data/RevisedQuPathExports"
OUT = S.OUT / "revised_cells"

#: from core/phenocycler/config.ini [geometry_qc]. Read, not guessed -- see load_thresholds().
CONFIG = S.REPO / "core/phenocycler/config.ini"

#: Every islet expansion is 20 um. Derived from the signed distance rather than read from an
#: annotation, because 6482's annotated expansion measures 23.1 um and that is an error.
EXPANSION_UM = 20.0

IDENT = {
    "Object ID": "object_id",
    "Image": "image",
    "Parent": "parent_raw",
    "Classification": "qupath_classification",
    "Centroid X µm": "X_centroid",
    "Centroid Y µm": "Y_centroid",
}
SHAPE = {
    "Nucleus: Area µm^2": "nucleus_area_um2",
    "Cell: Area µm^2": "cell_area",
    "Cell: Solidity": "cell_solidity",
    "Nucleus/Cell area ratio": "nucleus_cell_area_ratio",
}
#: 3_signed_distance.groovy writes "with Islet"; an older run of it wrote the form without "with".
#: Prefer the new one and fall back, rather than assuming which export is in hand.
DIST = ["Signed distance to annotation with Islet µm", "Signed distance to annotation Islet µm"]


def load_thresholds() -> dict[str, float]:
    import configparser
    cp = configparser.ConfigParser()
    cp.read(CONFIG)
    g = cp["geometry_qc"]
    return {
        "min_cell_area": g.getfloat("min_cell_area"),
        "max_area_median_multiple": g.getfloat("max_area_median_multiple"),
        "min_cell_solidity": g.getfloat("min_cell_solidity"),
        "no_cytoplasm_nc_ratio": g.getfloat("no_cytoplasm_nc_ratio"),
    }


def read_export(path: pathlib.Path) -> pd.DataFrame:
    """Only the columns needed. The exports are 8-17 GB and 994 columns wide."""
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    index = {n: i for i, n in enumerate(header)}
    dist_col = next((c for c in DIST if c in index), None)
    if dist_col is None:
        raise KeyError(f"{path.name}: no signed-distance-to-Islet column; 3_signed_distance did not run")

    wanted = {**IDENT, **SHAPE, dist_col: "signed_dist_islet_um"}
    missing = [c for c in wanted if c not in index]
    if missing:
        raise KeyError(f"{path.name}: missing {missing}")

    picks = [index[c] for c in wanted]
    names = [wanted[c] for c in wanted]
    rows = []
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            rows.append([row[i] for i in picks])
    df = pd.DataFrame(rows, columns=names)
    for c in list(SHAPE.values()) + ["X_centroid", "Y_centroid", "signed_dist_islet_um"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def islet_number(parent: pd.Series) -> pd.Series:
    """`Islet_269` -> 269. Anything that is not an islet annotation is NA, not 0."""
    return pd.to_numeric(parent.str.extract(r"Islet[_ ]?(\d+)", expand=False), errors="coerce")


def region(dist: pd.Series) -> pd.Series:
    """core / expansion / tissue, from the signed distance alone.

    Negative or zero is inside the islet; (0, 20] is the peri-islet expansion; beyond that is
    ordinary tissue. Built this way so every expansion is exactly 20 um by construction.
    """
    out = pd.Series(pd.NA, index=dist.index, dtype="object")
    finite = np.isfinite(dist)
    out[finite & (dist <= 0)] = "core"
    out[finite & (dist > 0) & (dist <= EXPANSION_UM)] = "expansion"
    out[finite & (dist > EXPANSION_UM)] = "tissue"
    return out


def eligibility(df: pd.DataFrame, thr: dict[str, float]):
    """(analysis_eligible, analysis_reason, estimation_eligible, estimation_reason).

    THE PIPELINE KEEPS THREE SEPARATE ROLES -- analysis, estimation, spillover context -- and a
    rule that disqualifies a cell for one does NOT disqualify it for the others. Collapsing them
    into one flag is not a conservative simplification, it is a different answer: the first
    version of this function applied `no_cytoplasm` to ANALYSIS and excluded 100,059 cells for
    6534 that the run counts as analysis-eligible. Measured in the run's own geometry QC:

        analysis_reason            no_cytoplasm = 0
        estimation_reason          no_cytoplasm = 99,024
        spillover_context_reason   no_cytoplasm = 0

    A nucleus that fills its whole cell gives you nothing to ESTIMATE a background level from;
    it is still a cell, and still counts.

    A cell with no nucleus is excluded from every role -- we cannot be sure it is a cell at all.
    Order matters: the reason recorded is the FIRST rule the cell failed, matching how the
    geometry QC reports one reason per cell.
    """
    area = df["cell_area"]
    med = np.nanmedian(area[np.isfinite(area)])
    shared = [
        ("missing_nucleus_geometry", ~np.isfinite(df["nucleus_area_um2"])),
        ("cell_area_below_min", area < thr["min_cell_area"]),
        ("cell_area_above_max", area > thr["max_area_median_multiple"] * med),
        ("low_solidity", df["cell_solidity"] < thr["min_cell_solidity"]),
    ]

    def apply(checks):
        reason = pd.Series("", index=df.index, dtype="object")
        for name, failed in checks:
            reason[(reason == "") & failed.fillna(False)] = name
        return reason == "", reason

    a_ok, a_why = apply(shared)
    e_ok, e_why = apply(shared + [
        ("no_cytoplasm", df["nucleus_cell_area_ratio"] > thr["no_cytoplasm_nc_ratio"]),
    ])
    return a_ok, a_why, e_ok, e_why


def rebuild(donor: str, thr: dict[str, float]) -> None:
    started = time.time()
    export = REVISED / f"{donor}.tsv"
    if not export.exists():
        print(f"  {donor}: no revised export at {export}"); return

    df = read_export(export)
    df.insert(0, "donor_id", donor)
    df["islet_num"] = islet_number(df["parent_raw"])
    df["cell_region"] = region(df["signed_dist_islet_um"])
    (df["qc_analysis_eligible"], df["analysis_reason"],
     df["qc_estimation_eligible"], df["estimation_reason"]) = eligibility(df, thr)

    # ---- REDSEA-corrected, UniFORM-aligned marker values, BY object_id -------------------
    aligned_path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    aligned = pq.read_table(aligned_path).to_pandas()
    markers = [c for c in aligned.columns
               if c not in {"donor_id", "object_id", "qc_analysis_eligible",
                            "X_centroid", "Y_centroid", "cell_region", "islet_num"}]
    merged = df.merge(aligned[["object_id"] + markers], on="object_id",
                      how="left", validate="one_to_one", indicator="_merge")

    # A revised cell with no aligned row has NO marker measurement. That is unavailable, not
    # zero -- the markers stay NaN and the count is reported rather than quietly dropped.
    #
    # COUNTED FROM THE MERGE INDICATOR, not from whether a marker column is NaN. The first
    # attempt probed markers[0] and reported 191,322 missing rows for 6534 when the true
    # number is 6,182: a marker is legitimately NaN for a cell that HAS a row, because not
    # every marker was acquired on every donor and REDSEA can leave one unset. "This column
    # is empty" and "this cell has no measurements" are different questions.
    unmeasured = int((merged["_merge"] == "left_only").sum())
    merged = merged.drop(columns="_merge")

    target = OUT / f"donor_id={donor}"
    target.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(target / "data.parquet", index=False)

    n = len(merged)
    elig = int(merged["qc_analysis_eligible"].sum())
    by_region = merged.loc[merged["qc_analysis_eligible"], "cell_region"].value_counts()
    print(f"  {donor}: {n:>9,} cells | eligible {elig:>9,} ({100*elig/n:4.1f}%) | "
          f"core {by_region.get('core', 0):>7,} expansion {by_region.get('expansion', 0):>7,} "
          f"tissue {by_region.get('tissue', 0):>9,} | "
          f"no marker row {unmeasured:>6,} | islets {int(merged['islet_num'].nunique()):>5,} | "
          f"{time.time()-started:.0f}s", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("donors", nargs="*", default=["6505", "6534", "6623"])
    args = ap.parse_args()
    thr = load_thresholds()
    print(f"geometry thresholds from config.ini: {thr}")
    print(f"expansion fixed at {EXPANSION_UM:.0f} um, derived from the signed distance\n")
    OUT.mkdir(parents=True, exist_ok=True)
    for donor in args.donors:
        rebuild(donor, thr)
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
