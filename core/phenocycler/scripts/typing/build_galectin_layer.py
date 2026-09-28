#!/usr/bin/env python3
"""The galectin layer -- v3: rebuilt from typed cells, on the UniFORM cohort scale.

**The maintainer's ruling, 2026-09-28.**  Correction J (below, historical) read three counts
straight off QuPath's own object-classifier verdicts in the annotation export.  The maintainer
reversed that choice for the v3 re-run: UniFORM is now the only normalization of record for this
cohort, and the QuPath object classifiers are retired for it.  Positivity and lineage now come
from the same typing, and the same reference-level machinery every other marker in this lane
already uses, so a reader comparing this layer to any other one in the app is comparing one
normalization throughout, not two.

Per cell, from ``data/development/typing/cells`` (the declared pancreas x PhenoCycler producer,
``docs/PRODUCERS.md``) and the UniFORM-aligned values (``subpanel_io.ALIGNED``):

    macrophage (Myeloid)   ``cell_class == "Immune"`` and ``immune_definition`` in
                            ``{"Imm-macrophage", "Imm-dendritic-cell"}`` -- the producer's OWN
                            per-cell phenotype column, written under ``--emit-detail`` (the
                            declared production command, ``docs/PRODUCERS.md``), grouped the way
                            ``scripts/senior/build_grid_composition.py`` / ``build_compartment_
                            gradient.py`` / ``modality_concordance.py`` already group it into the
                            common lineage "Myeloid" (``PANC_DEFINITION_TO_COMMON`` /
                            ``PANC_PHENO_IMMDEF_TO_COMMON``).  This is NOT ``lineages.py``'s own,
                            broader Myeloid marker list -- that module states what an antibody
                            *means*, before any result is seen; it never decides a call, and
                            ``cell_class`` carries no "Myeloid" value of its own (only the broad
                            "Immune").
    galectin-positive       M2Gal3, on the UniFORM cohort scale (``subpanel_io.apply_scale``),
                            at or above the donor's OWN M2Gal3 reference level -- the p99 of
                            M2Gal3 over that donor's reference population (cells in no
                            pixel-classified object and outside every islet, ``context ==
                            "tissue"`` in the typed-cells table), fitted with
                            ``subpanel_io.reference_level``, the SAME function every structure and
                            the hormone level are fitted with (``subpanel.fit_levels``).  No new
                            threshold is invented here.  The read-and-scale step
                            (``read_aligned_marker`` below) restates ``subpanel.scaled_markers``
                            for one marker rather than importing ``subpanel`` itself: that module
                            also imports ``context``/``islet_membership`` for its OWN geometry
                            work, which pulls in ``shapely`` -- a dependency this layer's lean
                            ``tests/current`` CI job does not install (`pyproject.toml`'s ``test``
                            extra, unlike ``integration``).  The same reasoning `apply_scale`'s own
                            docstring gives for duplicating rather than importing
                            ``celltype/score.apply_scale``.
    metabolically-active    galectin-positive and NOT Myeloid -- the maintainer's definition,
                            galectin without macrophage markers, unchanged from correction J.
    macrophage (denominator) the Myeloid count above, independent of M2Gal3.

A missing M2Gal3 value on a cell, or a donor whose M2Gal3 reference level cannot be fitted (the
donor's panel never carried it, or its reference population is empty), leaves that cell's
galectin state UNAVAILABLE -- never negative, never zero, the same rule every other missing
measurement in this lane follows.  ``n_cells_galectin_evaluable`` sits beside ``n_cells_eligible``
so the gap between "cells in this islet-region" and "cells a galectin call could be made on" stays
visible, the same "two denominators" idiom correction J used for ``n_cells_qupath`` /
``n_cells_eligible``.  The islet keys are the ones ``typing/cells`` already assigns
(``islet_key``, ``cell_region_20um`` mapped to core/peri) -- the same geometric islet outline
membership the rest of the app reads, restricted to cell_class != "Excluded" cells, which is the
current lane's reading of "class Islet only".

Column names, nullable-integer semantics and the app's contract are unchanged from correction J:
every ``count_<class>``/``prop_<class>`` pair keeps its name, and a nullable integer -- not a
0 -- is what a count that could not be made looks like.  ``method`` records
``"typed-cells-uniform-v3"`` on every row, and the same string plus the ruling date are written
into the Parquet file's own schema metadata.

**The retired QuPath-count path.**  ``--source qupath-counts`` runs the untouched pre-2026-09-28
method below (correction J): it reads ``results/AnnotationsMeasurements.csv`` and the two RTrees
classifier passes QuPath fit on ``Cell: M2Gal3: Mean``/``Max``/``Std.Dev.``.  It is kept, behind
an explicit flag, only because it is the table `tests/current/test_silent_failures_galectin_
mixture.py` exercises and because the annotation export cannot be recomputed outside QuPath (see
below) -- not because it is a candidate default.  ``--source typed-cells`` (the default) is the v3
method described above.

    python core/phenocycler/scripts/typing/build_galectin_layer.py
    python core/phenocycler/scripts/typing/build_galectin_layer.py --source qupath-counts

---

**Historical -- correction J, 2026-09-21.  Superseded as the default by the 2026-09-28 ruling
above; kept, verbatim, behind ``--source qupath-counts``.**

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
"""
from __future__ import annotations

import argparse
import re
import sys

import datetime as dt

import numpy as np
import pandas as pd
import pyarrow as pa
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


def _main_qupath_counts(donors: list[str]) -> int:
    """The pre-2026-09-28 method (correction J), untouched: ``--source qupath-counts``."""
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


# =============================================================================================
# v3: rebuilt from typed cells, on the UniFORM cohort scale (the maintainer's ruling, 2026-09-28)
# =============================================================================================

#: The marker this layer classifies.  Spelled the way the panel registry and `lineages.py` spell
#: it (`lineages.LINEAGE_MARKERS["Myeloid"]` carries it, and `lineages.SUPPORTING_ONLY` too -- it
#: may not carry a lineage call on its own, which is exactly why this layer does not ask
#: `lineages.py` for the call and reads the typing's own `immune_definition` instead).
MARKER = "M2Gal3"

#: `immune_definition` values `scripts/senior/build_grid_composition.py` (`PANC_DEFINITION_TO_
#: COMMON`), `build_compartment_gradient.py` and `modality_concordance.py` (`PANC_PHENO_IMMDEF_
#: TO_COMMON`) already group into the common lineage "Myeloid".  Restated here rather than
#: imported -- those modules pull in the whole grid/compartment machinery for one small dict, the
#: same reason `build_cells.py` restates `refine.UNRESOLVED` rather than importing `refine`.  The
#: three call sites are named so a maintainer changing the grouping there is told to change it
#: here too.
MYELOID_DEFINITIONS = ("Imm-macrophage", "Imm-dendritic-cell")

#: `cell_region_20um` (typed cells) -> the app's word for the zone, matching `REGION_OF` above.
TYPED_REGION_OF = {"core": "core", "expansion": "peri"}

#: Columns this method reads off `typing/cells`.  `context` is written on every row regardless of
#: `--emit-detail`; `immune_definition` is not, and its absence means the declared producer's own
#: production command (`docs/PRODUCERS.md`) was not followed for this build.
TYPED_CELLS_COLUMNS = ("object_id", "cell_class", "context", "immune_definition",
                       "islet_key", "cell_region_20um")

METHOD = "typed-cells-uniform-v3"
RULING_DATE = "2026-09-28"


def read_aligned_marker(donor: str, marker: str, scale: pd.DataFrame,
                        object_ids: pd.Series) -> np.ndarray:
    """One marker's UniFORM-aligned value, on the cohort scale, joined to ``object_ids``.

    Restates ``subpanel.scaled_markers`` for a single marker rather than importing ``subpanel``
    itself -- see the module docstring for why: that module's own geometry work pulls in
    ``shapely``, which this layer's lean CI job does not install.  Same contract: a marker absent
    from this donor's panel yields an all-NaN column, by name, never by position -- ``NaN >=
    level`` is False, so a missing marker never fires as a positive.
    """
    path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    have = set(pq.ParquetFile(path).schema_arrow.names)
    if marker not in have:
        return np.full(len(object_ids), np.nan, dtype=float)
    frame = pq.read_table(path, columns=["object_id", marker]).to_pandas()
    frame = object_ids.to_frame("object_id").merge(frame, on="object_id", how="left")
    return S.apply_scale(frame[[marker]].to_numpy(), [marker], scale)[:, 0]


def donor_m2gal3_reference_level(cells: pd.DataFrame, scale: pd.DataFrame) -> tuple[np.ndarray, float | None]:
    """This donor's M2Gal3 values (UniFORM cohort scale, one per row of ``cells``) and its
    reference level -- the SAME machinery every structure and the hormone level are fitted with
    (``read_aligned_marker`` / ``subpanel.scaled_markers``, ``subpanel_io.reference_level``),
    never a new threshold.

    The reference population is ``context == "tissue"``: cells in no pixel-classified object and
    outside every islet, exactly as ``subpanel.fit_levels`` defines it -- read here off the
    typing's own carried ``context`` column rather than recomputed, so this layer cannot drift
    from the classification it is layered on.

    Unlike ``fit_levels``, which fails loudly when a level will not fit, an unfittable donor here
    is UNAVAILABLE rather than a crash: a donor whose panel never carried M2Gal3, or whose tissue
    population is empty, is an expected outcome for one marker on one donor, not a defect in the
    whole cohort build.  Returns ``None`` for the level in that case, and every cell of this
    donor's galectin state is UNAVAILABLE from it.
    """
    m2gal3 = read_aligned_marker(cells.attrs["donor"], MARKER, scale, cells["object_id"])
    tissue = cells["context"].to_numpy() == "tissue"
    level = None
    if tissue.any():
        try:
            level = S.reference_level(m2gal3[tissue])
        except ValueError:
            level = None
    return m2gal3, level


def load_typed_donor(donor: str) -> pd.DataFrame:
    """One donor's typed cells, with the columns this layer reads and no others."""
    path = S.OUT / f"cells/donor_id={donor}/data.parquet"
    have = set(pq.ParquetFile(path).schema_arrow.names)
    missing = [c for c in TYPED_CELLS_COLUMNS if c not in have]
    if missing:
        raise KeyError(
            f"{donor}: typing/cells lacks {missing} -- `immune_definition` is written only under "
            "`subpanel.py --emit-detail`, the declared production command (docs/PRODUCERS.md); "
            "rerun it with that flag.")
    cells = pq.read_table(path, columns=list(TYPED_CELLS_COLUMNS)).to_pandas()
    cells.attrs["donor"] = donor
    return cells


def build_donor_galectin(donor: str, scale: pd.DataFrame) -> tuple[pd.DataFrame, float | None]:
    """One donor's admitted islet-region cells, with a galectin state and a Myeloid flag each.

    "Admitted" is `cell_class != "Excluded"` and `islet_key` not null -- this lane's reading of
    "class Islet only": the geometric islet-outline membership `typing/cells` already carries,
    with the in-islet exocrine misclassification (`Excluded`) removed, same as correction J's own
    `n_cells_eligible`.
    """
    cells = load_typed_donor(donor)
    m2gal3, level = donor_m2gal3_reference_level(cells, scale)
    cells["m2gal3_scaled"] = m2gal3

    admitted = cells[(cells["cell_class"] != "Excluded") & cells["islet_key"].notna()].copy()
    admitted["region"] = admitted["cell_region_20um"].map(TYPED_REGION_OF)
    admitted = admitted.dropna(subset=["region"])

    state = np.full(len(admitted), UNAVAILABLE, dtype=object)
    if level is not None:
        known = admitted["m2gal3_scaled"].notna().to_numpy()
        state[known] = np.where(
            admitted.loc[known, "m2gal3_scaled"].to_numpy() >= level, "positive", "negative")
    admitted["galectin_state"] = state
    admitted["is_myeloid"] = (
        (admitted["cell_class"] == "Immune")
        & admitted["immune_definition"].fillna("").isin(MYELOID_DEFINITIONS))
    admitted["case_id"] = donor
    return admitted, level


def derive_counts_from_typed_cells(donors: list[str], scale: pd.DataFrame | None = None
                                   ) -> pd.DataFrame:
    """Per-islet-region counts and proportions, from typed cells on the UniFORM cohort scale.

    Mirrors ``derive_counts``'s contract exactly: the same four ``count_<class>`` /
    ``prop_<class>`` names, nullable integers, and a ``galectin_state`` / ``galectin_reason`` pair
    naming which rows lost a count and why.  A count that cannot be made -- the donor's M2Gal3
    reference level is unavailable, or not one admitted cell in this islet-region carries a
    measured M2Gal3 value -- is null, never 0; a count that needs no M2Gal3 at all
    (``count_Macrophage``) is unaffected by either.
    """
    scale = S.load_scale() if scale is None else scale
    donor_levels: dict[str, float | None] = {}
    frames = []
    for donor in donors:
        admitted, level = build_donor_galectin(donor, scale)
        donor_levels[donor] = level
        if len(admitted):
            frames.append(admitted)
    if not frames:
        raise ValueError(f"no admitted islet-region cells across {len(donors)} donor(s)")
    cells = pd.concat(frames, ignore_index=True)

    cells["is_known"] = cells["galectin_state"] != UNAVAILABLE
    cells["is_positive"] = cells["galectin_state"] == "positive"
    cells["is_macrophage_positive"] = cells["is_positive"] & cells["is_myeloid"]
    cells["is_active_only"] = cells["is_positive"] & ~cells["is_myeloid"]

    out = cells.groupby(["case_id", "islet_key", "region"], sort=True).agg(
        n_cells_eligible=("object_id", "size"),
        n_cells_galectin_evaluable=("is_known", "sum"),
        count_Macrophage=("is_myeloid", "sum"),
        count_Galectin_positive=("is_positive", "sum"),
        count_Galectin_macrophage=("is_macrophage_positive", "sum"),
        count_Metabolically_active=("is_active_only", "sum"),
    ).reset_index()

    no_coverage = (out["n_cells_galectin_evaluable"] == 0).to_numpy()
    galectin_labels = ["Galectin_positive", "Galectin_macrophage", "Metabolically_active"]
    for label in galectin_labels:
        column = f"count_{label}"
        out[column] = out[column].astype("Int64")
        out.loc[no_coverage, column] = pd.NA
    out["count_Macrophage"] = out["count_Macrophage"].astype("Int64")

    def reason(row) -> str | None:
        if row["n_cells_galectin_evaluable"] > 0:
            return None
        names = ", ".join(f"count_{label}" for label in galectin_labels)
        level = donor_levels.get(row["case_id"])
        if level is None:
            return (f"{names}: donor {row['case_id']}'s M2Gal3 reference level is unavailable "
                    "-- its panel did not carry M2Gal3, or its tissue reference population "
                    "(context == \"tissue\") was empty")
        return f"{names}: no admitted cell in this islet-region carries a measured M2Gal3 value"

    out[COUNT_REASON] = pd.array([reason(r) for _, r in out.iterrows()], dtype="string")
    out[COUNT_STATE] = np.where(no_coverage, UNAVAILABLE, MEASURED)
    out["m2gal3_reference_level"] = out["case_id"].map(donor_levels).astype("float64")
    out["method"] = METHOD

    denominator = out["n_cells_eligible"].astype("Float64")
    for label in list(COUNT_COLUMNS.values()) + ["Metabolically_active"]:
        ratio = (out[f"count_{label}"].astype("Float64") / denominator).where(denominator > 0)
        out[f"prop_{label}"] = ratio.astype("float64")
    return out


def _main_typed_cells(donors: list[str]) -> int:
    """The default method: rebuilt from typed cells, per the maintainer's 2026-09-28 ruling."""
    scale = S.load_scale()
    out = derive_counts_from_typed_cells(donors, scale)

    props = [c for c in out.columns if c.startswith("prop_")]
    out = out[[c for c in out.columns if c not in props] + props]

    target = APP / "islet_galectin.parquet"
    table = pa.Table.from_pandas(out, preserve_index=False)
    table = table.replace_schema_metadata({
        **(table.schema.metadata or {}),
        b"method": METHOD.encode(),
        b"ruling_date": RULING_DATE.encode(),
        b"built_by": b"core/phenocycler/scripts/typing/build_galectin_layer.py",
        b"built_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").encode(),
    })
    pq.write_table(table, target)

    print(f"{len(out):,} islet-regions over {out['case_id'].nunique()} donors "
          f"({(out['region'] == 'core').sum():,} core, {(out['region'] == 'peri').sum():,} peri) "
          f"-- method={METHOD}")
    fittable = sum(1 for v in out.drop_duplicates('case_id')['m2gal3_reference_level'] if pd.notna(v))
    print(f"M2Gal3 reference level fitted for {fittable} of {out['case_id'].nunique()} donors")
    print(f"\n{int((out[COUNT_STATE] == UNAVAILABLE).sum()):,} of {len(out):,} islet-regions carry "
          f"an unavailable count (`{COUNT_REASON}` names which and why):")
    for label in list(COUNT_COLUMNS.values()) + ["Metabolically_active"]:
        print(f"  count_{label}: {int(out[f'count_{label}'].isna().sum()):,} unavailable")
    print("\nmedian proportion per class, by region:")
    for region in ("core", "peri"):
        here = out[out["region"] == region]
        stats = {label: here[f"prop_{label}"].median()
                 for label in list(COUNT_COLUMNS.values()) + ["Metabolically_active"]}
        print(f"  {region:>5}: " + "  ".join(f"{k} {v:.4f}" for k, v in stats.items()))
    print(f"\n{len(out):,} islet-regions -> {target}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--source", choices=("typed-cells", "qupath-counts"),
                        default="typed-cells",
                        help="typed-cells (default): rebuilt from typed cells on the UniFORM "
                             "cohort scale, per the maintainer's 2026-09-28 ruling. "
                             "qupath-counts: the retired pre-2026-09-28 method (correction J), "
                             "QuPath's own object-classifier counts read from the annotation "
                             "export -- kept behind this flag, never the default.")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    if args.source == "qupath-counts":
        return _main_qupath_counts(donors)
    return _main_typed_cells(donors)


if __name__ == "__main__":
    raise SystemExit(main())
