#!/usr/bin/env python3
"""The galectin layer -- correction J.

Three classes per islet and per expansion, all of them QuPath's own object-classifier verdicts read
from the annotation export:

    galectin-positive        `Num M2Galectin3`
    galectin+ macrophage     `Num Iba1_CD68: M2Galectin3`
    metabolically-active     `Num M2Galectin3` - `Num Iba1_CD68: M2Galectin3`
                             -- galectin without macrophage markers, the maintainer's definition;
                             null where the difference is negative (see below)
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

**A count that was not exported, or is not a count, is null -- not 0.**  A blank cell in one of
the three count columns stays null (the count is UNAVAILABLE), and `count_Metabolically_active`
is null where `Num Iba1_CD68: M2Galectin3` exceeds `Num M2Galectin3`: the two object-classifier
passes disagree there, and a negative difference is not a count to clip.  `galectin_state` and
`galectin_reason` say which rows and why; the proportions follow their counts.  Until 2026-09-21 a
blank became the integer 0 (`fillna(0)`) and a negative difference was clipped to 0
(`clip(lower=0)`), both indistinguishable from a measured zero.  Measured on the 2026-09-14
layer (15,557 islet-regions): no blank count, and 719 islet-regions where the passes disagree,
every one of them written as 0.  Column names are unchanged; the counts are nullable integers.

**The denominator is a count too.**  A blank `Num Detections` leaves every proportion on the row
without a denominator, so the row is `unavailable` and `galectin_reason` names `n_cells_qupath`
-- a row whose four counts were all read but whose proportions are all null used to say
`measured`.  A `Num Detections` of exactly 0 is present, not missing: the proportion is undefined
(null, as it always was) and the row stays `measured`.

**A value that cannot be a count is refused where the column is read.**  Each of the four count
columns -- the three above and `Num Detections` -- passes through `read_count`, which admits a
blank (the count is UNAVAILABLE, and the row names it) or a finite non-negative integer that a
64-bit count can hold, and refuses everything else on the row by name, with the value: a
negative, a non-finite, a fractional, a non-numeric cell, or a boolean column.  An integral
float such as `10.0` and a numeric string such as `"10"` are admitted as the integer they spell,
because a CSV column with one blank cell reads as float64 and one with one stray character reads
as object.  Until 2026-09-22 a negative or fractional `Num Detections` went straight into the
proportions (measured on this build: -5 gave `measured` with every proportion NaN, `inf` gave
`measured` with every proportion 0.0, 2.5 gave `measured` with a proportion of 2.0), a count of
-1 was cast to -1 and reported a proportion of -0.1 under the classifier-disagreement reason, and
an `inf` or a 2.5 in a count column raised at the `Int64` cast for the whole frame.  Only admitted
values now reach that cast, so it cannot raise.

    python core/multiplex5000/scripts/typing/build_galectin_layer.py
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

#: Beside the four counts on every row.  ``measured`` when each count on the row was read from
#: the export and the difference is a count; otherwise ``unavailable``, with ``galectin_reason``
#: naming which count and why.  A missing measurement is unavailable, not zero.
COUNT_STATE = "galectin_state"
COUNT_REASON = "galectin_reason"
MEASURED = "measured"
UNAVAILABLE = "unavailable"


def parse_annotations(donors: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(ANNOTATIONS)
    frame["case_id"] = frame["Image"].str.extract(r"(\d{3,4})")[0]
    frame = frame[frame["case_id"].isin(donors)].copy()

    parsed = frame["Name"].str.extract(r"^(Islet_\d+)(_exp20um)?$")
    frame["islet_key"] = parsed[0]
    frame["region"] = parsed[1].fillna("").map(REGION_OF)
    frame = frame[frame["islet_key"].notna() & frame["region"].notna()]
    return frame


#: A value at or above this is an integer a 64-bit count cannot hold; refusing it at the read
#: is what keeps the `Int64` cast from ever overflowing.
INT64_CEILING = 2.0 ** 63


def _shown(value) -> str:
    """The cell as it was read, for the reason: a string in quotes, a number as itself."""
    return repr(value) if isinstance(value, str) else str(value)


def read_count(values: pd.Series) -> tuple[pd.Series, dict[int, str]]:
    """The one boundary where a count column of the annotation export is read.

    Returns the column as a nullable integer and, per refused position, why.  A cell is
    admitted when it is blank (null: UNAVAILABLE, and the caller names it) or a finite
    non-negative integer that a 64-bit count can hold.  An integral float such as ``10.0`` and
    a numeric string such as ``"10"`` are that integer: a CSV column with one blank cell reads
    as float64, and one with one stray character reads as object, so the dtype says nothing
    about the cells.  Anything else cannot be a count -- negative, non-finite, fractional,
    non-numeric, or a boolean column, where True is not a count of 1 -- and is refused HERE:
    null in the returned column, and named in ``refused`` with the value as read.  Only
    admitted values reach the ``Int64`` cast, so the cast cannot raise; it used to, on an
    ``inf`` (OverflowError) or a ``2.5`` (TypeError), for the whole frame at once.
    """
    blank = values.isna().to_numpy(dtype=bool)
    if pd.api.types.is_bool_dtype(values.dtype):
        as_float = np.full(len(values), np.nan)
        why_column = "boolean"
    else:
        as_float = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
        why_column = None
    with np.errstate(invalid="ignore"):
        admitted = (np.isfinite(as_float) & (as_float >= 0) & (as_float == np.floor(as_float))
                    & (as_float < INT64_CEILING))

    refused: dict[int, str] = {}
    for i in np.flatnonzero(~blank & ~admitted):
        x = as_float[i]
        why = (why_column if why_column is not None
               else "not a number" if np.isnan(x)
               else "not finite" if not np.isfinite(x)
               else "negative" if x < 0
               else "not integral" if x != np.floor(x)
               else "beyond a 64-bit count")
        refused[int(i)] = f"{_shown(values.iat[i])} ({why})"

    counts = pd.Series(np.where(admitted, as_float, np.nan), index=values.index).astype("Int64")
    return counts, refused


def derive_counts(frame: pd.DataFrame) -> pd.DataFrame:
    """The per-islet-region counts and proportions from the annotation rows.

    ``count_<class>`` is a nullable integer: a blank cell in the export stays null instead of
    becoming 0.  ``count_Metabolically_active`` is the difference of two object-classifier passes
    and is null where the second exceeds the first -- the passes disagree there, and a negative
    number of cells is not a count to clip.  ``prop_<class>`` follows its count: null where the
    count is, and where ``n_cells_qupath`` is 0, exactly as before.  A blank ``Num Detections``
    is a missing denominator, and the row is unavailable for that reason -- a 0 is present and
    is not.  ``galectin_state`` and ``galectin_reason`` say which rows lost a count and why.  A
    row that carries every measurement is numerically what it was before 2026-09-21.  A value
    that cannot be a count -- negative, non-finite, fractional, non-numeric -- is refused where
    its column is read (`read_count`), null from there on, and named with its value.
    """
    out = frame[["case_id", "islet_key", "region"]].copy()
    reasons: list[list[str]] = [[] for _ in range(len(out))]

    # Every count column of the export crosses ONE boundary, `read_count`; a blank or a refused
    # cell is named on its row in the same words, the value in place of `blank`.
    read = {column: read_count(frame[column]) for column in ["Num Detections", *COUNT_COLUMNS]}

    def name(column: str, output: str, tail: str = "") -> None:
        for i in np.flatnonzero(frame[column].isna().to_numpy()):
            reasons[i].append(f"{output}: `{column}` is blank in the annotation export{tail}")
        for i, value in read[column][1].items():
            reasons[i].append(f"{output}: `{column}` is {value} in the annotation export, "
                              f"not a count{tail}")

    out["n_cells_qupath"] = read["Num Detections"][0]
    for column, label in COUNT_COLUMNS.items():
        out[f"count_{label}"] = read[column][0]
        name(column, f"count_{label}")

    # The maintainer's definition: galectin without macrophage markers -- the difference of two
    # classifier passes, which is a count only where the passes agree on its sign.
    positive, macrophage = out["count_Galectin_positive"], out["count_Galectin_macrophage"]
    difference = positive - macrophage
    disagree = (difference < 0).fillna(False).to_numpy(dtype=bool)
    out["count_Metabolically_active"] = difference.mask(disagree)
    for i in np.flatnonzero(disagree):
        reasons[i].append(
            f"count_Metabolically_active: `Num Iba1_CD68: M2Galectin3` ({macrophage.iat[i]}) "
            f"exceeds `Num M2Galectin3` ({positive.iat[i]}) -- the two QuPath object-classifier "
            f"passes disagree, and a negative difference is not a count")
    for i in np.flatnonzero(difference.isna().to_numpy()):
        reasons[i].append("count_Metabolically_active: one of its two inputs is unavailable")

    # The denominator of every proportion.  Blank or refused, it is a missing measurement and
    # the row says so; 0 is a measurement, and the proportion it leaves undefined is null below
    # as before.
    denominator = out["n_cells_qupath"]
    name("Num Detections", "n_cells_qupath",
         tail=" -- no denominator, so no proportion on this row")

    out[COUNT_STATE] = [UNAVAILABLE if r else MEASURED for r in reasons]
    out[COUNT_REASON] = pd.Series(["; ".join(r) if r else None for r in reasons],
                                  index=out.index, dtype="string")

    for label in list(COUNT_COLUMNS.values()) + ["Metabolically_active"]:
        ratio = (out[f"count_{label}"].astype("Float64") / denominator).where(denominator > 0)
        out[f"prop_{label}"] = ratio.astype("float64")
    return out


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

    out = derive_counts(frame)

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
    # The eligible count sits where it always did: before the proportions, after the counts.
    props = [c for c in out.columns if c.startswith("prop_")]
    out = out[[c for c in out.columns if c not in props] + props]

    target = APP / "islet_galectin.parquet"
    out.to_parquet(target, index=False)

    # Unavailable is its own quantity, printed even when it is zero.
    print(f"\n{int(out[COUNT_STATE].eq(UNAVAILABLE).sum()):,} of {len(out):,} islet-regions carry "
          f"an unavailable count (`{COUNT_REASON}` names which and why):")
    for label in list(COUNT_COLUMNS.values()) + ["Metabolically_active"]:
        print(f"  count_{label}: {int(out[f'count_{label}'].isna().sum()):,} unavailable")
    print(f"  n_cells_qupath: {int(out['n_cells_qupath'].isna().sum()):,} unavailable "
          f"(`Num Detections` blank -- no proportion on those rows)")

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
