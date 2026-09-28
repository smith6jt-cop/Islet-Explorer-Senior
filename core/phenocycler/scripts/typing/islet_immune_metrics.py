#!/usr/bin/env python3
"""Per-islet and per-donor immune infiltration metrics.

A pure computation module: it takes two already-typed tables and returns metrics, and reads
no cohort file, manifest or hardcoded ``data/`` path of its own. That keeps it usable against
whatever a caller's typing lane currently produces (the v3 same-slide re-run, the 20-donor
islet chain, or a synthetic fixture in a test) without this module having an opinion about
which one is current.

INPUTS

``cells`` -- one row per typed, nucleated cell (a cell with no nucleus is excluded upstream,
per the repository's standing rule; this module does not re-check that):

    donor_id                 str
    object_id                str    the QuPath UUID, carried verbatim
    x_um, y_um                float  centroid, microns
    cell_class                str    the repo's typing label (``cell_class`` in
                                     ``data/development/typing/cells``). CD45 positivity is
                                     read from here: ``cell_class in cd45_positive_classes``
                                     (default ``{"Immune"}``) -- not re-derived from a marker
                                     threshold, per the repository's rule against inferring
                                     positivity outside the typing lane.
    immune_class               str    OPTIONAL. The immune second pass's finer label
                                     (``immune_subtype.py``'s ``immune_class``). T-cell (CD3+)
                                     membership is read from here:
                                     ``immune_class in t_cell_immune_classes`` (default
                                     ``{"Imm-T-cell"}``). Absent entirely, every T-cell-derived
                                     metric is reported ``UNAVAILABLE`` rather than 0.
    cell_region                str    OPTIONAL context column, carried through unexamined --
                                     core/peri membership here is always derived from
                                     ``signed_dist_islet_um`` at the CONFIGURED band width,
                                     which need not match whatever a fixed ``_20um`` column
                                     upstream was cut at.
    islet_id                   the id of the NEAREST ``Islet``-class object to this cell (the
                                same object ``signed_dist_islet_um`` is measured against), or
                                null if the cell has no islet association at all (e.g. beyond
                                every configured cutoff). A cell's association is with exactly
                                one islet.
    signed_dist_islet_um       float  negative inside the islet outline, zero at the border,
                                     positive outside -- the repository's standing sign
                                     convention (``build_cells.py``, ``islet_membership.py``).

``islets`` -- one row per QuPath object the upstream export classified as endocrine tissue:

    donor_id                   str
    islet_id                   matches ``cells.islet_id``
    class                      str    ``"Islet"``, ``"EndocrineCluster"`` or
                                     ``"EndocrineSingle"``. Only ``"Islet"`` rows are islets
                                     for this module's per-islet and per-donor metrics --
                                     Drotar et al. 2025's >=35 um mean caliper diameter is
                                     what the upstream classifier used to draw that line, and
                                     this module trusts the label rather than re-deriving it.
                                     The other two classes are surfaced separately (see
                                     :data:`IsletImmuneMetrics.non_islet_endocrine_objects`)
                                     and never folded into an islet's counts.
    area_um2                   float
    mean_caliper_diameter_um   float

Every count in this module is a count over CLASSIFIED cells. A cell whose relevant label
(``cell_class`` for CD45, ``immune_class`` for T cells) is missing does not count as negative
-- it is simply not counted, and when NO cell in a region has that label at all, the metric
for that region is ``pd.NA`` rather than 0. A count of zero is what this module reports when
it has positive evidence of zero: cells were classified and none of them qualified.
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ contract

ISLET_CLASS = "Islet"
NON_ISLET_ENDOCRINE_CLASSES = ("EndocrineCluster", "EndocrineSingle")

#: Drotar et al. 2025 -- the caliper-diameter floor the upstream classifier drew "Islet" at.
#: Not re-applied here (see module docstring); recorded so a caller can audit against it.
MIN_ISLET_MEAN_CALIPER_DIAMETER_UM = 35.0

CELLS_REQUIRED_COLUMNS = (
    "donor_id", "object_id", "x_um", "y_um", "cell_class",
    "islet_id", "signed_dist_islet_um",
)
CELLS_OPTIONAL_COLUMNS = ("immune_class", "cell_region")

ISLETS_REQUIRED_COLUMNS = (
    "donor_id", "islet_id", "class", "area_um2", "mean_caliper_diameter_um",
)


class IsletImmuneMetricsError(ValueError):
    """A cells or islets table does not satisfy this module's contract."""


def _require_columns(df: pd.DataFrame, columns: tuple, *, name: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise IsletImmuneMetricsError(f"{name} table is missing required column(s): {missing}")


def _validate_cells(cells: pd.DataFrame) -> None:
    _require_columns(cells, CELLS_REQUIRED_COLUMNS, name="cells")


def _validate_islets(islets: pd.DataFrame) -> None:
    _require_columns(islets, ISLETS_REQUIRED_COLUMNS, name="islets")


# ------------------------------------------------------------------ config

@dataclass(frozen=True)
class IsletImmuneMetricsConfig:
    """Every threshold this module applies, in one place and all parameterizable.

    Defaults follow the citations named in the maintainer's request: a 20 um peri-islet band
    (matching the repository's existing ``Islet_N_exp20um`` convention), a >=6 T-cell
    insulitis flag per Bodenmiller/Damond, and the nPOD histological insulitis definition
    (>=15 CD45+ cells in >=3 islets).
    """

    peri_band_um: float = 20.0
    insulitis_t_cell_threshold: int = 6
    infiltration_cutoff_um: float = 100.0
    distance_bins_um: tuple = ((0.0, 20.0), (20.0, 40.0), (40.0, 60.0), (60.0, 100.0))
    npod_min_cd45_per_islet: int = 15
    npod_min_qualifying_islets: int = 3
    cd45_positive_classes: frozenset = field(default_factory=lambda: frozenset({"Immune"}))
    t_cell_immune_classes: frozenset = field(default_factory=lambda: frozenset({"Imm-T-cell"}))


DEFAULT_CONFIG = IsletImmuneMetricsConfig()


# ------------------------------------------------------------------ small helpers

def _classified_positive_count(labels: pd.Series, positive: frozenset):
    """``(n_positive, n_classified)``, both plain ints. ``labels`` may contain null/"" entries.

    A cell with a blank or null label is neither positive nor negative evidence -- it is
    simply not classified, and is excluded from both counts. Callers turn ``n_classified == 0``
    into ``pd.NA`` at the point a rate or flag is computed; the raw counts stay here as plain
    ints so they can still be summed across islets for a donor total.
    """
    text = labels.astype(str).str.strip()
    classified = text.ne("") & labels.notna()
    positive_mask = classified & text.isin(positive)
    return int(positive_mask.sum()), int(classified.sum())


def _count_or_na(n_positive: int, n_classified: int, n_total: int):
    """0 when the region genuinely has no cells at all (a confident zero, nothing to infer);
    ``pd.NA`` when cells are present but none of them carry a determinate label; otherwise the
    classified positive count."""
    if n_total == 0:
        return 0
    if n_classified == 0:
        return pd.NA
    return n_positive


def _safe_ratio(numerator, denominator):
    if numerator is pd.NA or denominator is pd.NA or denominator is None or denominator == 0:
        return pd.NA
    return float(numerator) / float(denominator)


# ------------------------------------------------------------------ per-islet metrics

def compute_per_islet_metrics(
    cells: pd.DataFrame,
    islets: pd.DataFrame,
    config: IsletImmuneMetricsConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """One row per ``class == "Islet"`` object: CD45/T-cell counts, insulitis, infiltration.

    Cells are matched to an islet by ``(donor_id, islet_id)`` -- each cell already carries the
    id of the single nearest ``Islet``-class object, so no distance computation happens
    against any OTHER islet's geometry here. A cell whose ``islet_id`` names an object that is
    not ``class == "Islet"`` (i.e. an EndocrineCluster/EndocrineSingle) contributes to no row
    of this table, which is what keeps those objects out of every islet's counts.
    """
    _validate_cells(cells)
    _validate_islets(islets)

    has_immune_class = "immune_class" in cells.columns
    islet_rows = islets[islets["class"] == ISLET_CLASS].copy()

    cells = cells.copy()
    cells["donor_id"] = cells["donor_id"].astype(str)
    islet_rows["donor_id"] = islet_rows["donor_id"].astype(str)
    cells["islet_id"] = cells["islet_id"].astype(str)
    islet_rows["islet_id"] = islet_rows["islet_id"].astype(str)

    d = cells["signed_dist_islet_um"]
    in_core = d <= 0
    in_peri = (d > 0) & (d <= config.peri_band_um)
    in_core_or_peri = in_core | in_peri

    rows = []
    for islet in islet_rows.to_dict("records"):
        key = cells["islet_id"].to_numpy() == islet["islet_id"]
        donor_key = cells["donor_id"].to_numpy() == str(islet["donor_id"])
        mine = key & donor_key

        core = mine & in_core.to_numpy()
        peri = mine & in_peri.to_numpy()
        core_or_peri = mine & in_core_or_peri.to_numpy()

        n_cells_core = int(core.sum())

        n_cd45_core, n_cd45_core_classified = _classified_positive_count(
            cells.loc[core, "cell_class"], config.cd45_positive_classes)
        n_cd45_peri, n_cd45_peri_classified = _classified_positive_count(
            cells.loc[peri, "cell_class"], config.cd45_positive_classes)

        n_cd45_core_val = _count_or_na(n_cd45_core, n_cd45_core_classified, int(core.sum()))
        n_cd45_peri_val = _count_or_na(n_cd45_peri, n_cd45_peri_classified, int(peri.sum()))
        cd45_per_nucleus_core = _safe_ratio(n_cd45_core_val, n_cells_core)
        has_ge1_cd45_core = (
            pd.NA if n_cd45_core_val is pd.NA else bool(n_cd45_core_val >= 1)
        )

        if has_immune_class:
            n_t_cp, n_t_cp_classified = _classified_positive_count(
                cells.loc[core_or_peri, "immune_class"], config.t_cell_immune_classes)
            n_t_cells_val = _count_or_na(n_t_cp, n_t_cp_classified, int(core_or_peri.sum()))
            insulitic = (
                pd.NA if n_t_cells_val is pd.NA
                else bool(n_t_cells_val >= config.insulitis_t_cell_threshold)
            )
        else:
            n_t_cells_val, insulitic = pd.NA, pd.NA

        # Continuous infiltration score, over CD45+ cells only: 1 inside, 1/d outside up to
        # the cutoff, 0 (excluded) beyond it. 0.0 when the islet has no cells within range at
        # all; NA when cells exist there but none of them carry a determinate cell_class;
        # otherwise the true weighted sum, which is itself a real 0.0 when every classified
        # cell in range turns out CD45-negative.
        within_cutoff = mine & (d.to_numpy() <= config.infiltration_cutoff_um)
        infil_labels = cells.loc[within_cutoff, "cell_class"]
        _, n_infil_classified = _classified_positive_count(infil_labels, config.cd45_positive_classes)
        if int(within_cutoff.sum()) == 0:
            infiltration_score = 0.0
        elif n_infil_classified == 0:
            infiltration_score = pd.NA
        else:
            dd = d.to_numpy()[within_cutoff]
            text = cells.loc[within_cutoff, "cell_class"].astype(str).str.strip()
            classified = (text.ne("")) & cells.loc[within_cutoff, "cell_class"].notna()
            is_cd45 = (classified & text.isin(config.cd45_positive_classes)).to_numpy()
            weight = np.where(dd <= 0, 1.0, np.where(dd > 0, 1.0 / np.maximum(dd, 1e-9), 0.0))
            infiltration_score = float(weight[is_cd45].sum())

        row = {
            "donor_id": str(islet["donor_id"]), "islet_id": islet["islet_id"],
            "class": islet["class"],
            "area_um2": islet["area_um2"],
            "mean_caliper_diameter_um": islet["mean_caliper_diameter_um"],
            "n_cells_core": n_cells_core,
            "n_cd45_positive_core": n_cd45_core_val,
            "n_cd45_positive_peri": n_cd45_peri_val,
            "cd45_per_nucleus_core": cd45_per_nucleus_core,
            "has_ge1_cd45_positive_core": has_ge1_cd45_core,
            "n_t_cells_core_peri": n_t_cells_val,
            "insulitic": insulitic,
            "infiltration_score": infiltration_score,
        }
        rows.append(row)

    columns = [
        "donor_id", "islet_id", "class", "area_um2", "mean_caliper_diameter_um",
        "n_cells_core", "n_cd45_positive_core", "n_cd45_positive_peri",
        "cd45_per_nucleus_core", "has_ge1_cd45_positive_core",
        "n_t_cells_core_peri", "insulitic", "infiltration_score",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns)


def compute_distance_bins(
    cells: pd.DataFrame,
    islets: pd.DataFrame,
    config: IsletImmuneMetricsConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Long-format counts per islet per outward distance bin: total / CD45+ / T-cell cells.

    Bins are measured from the islet boundary (``signed_dist_islet_um``), independent of --
    and may overlap -- the peri band used elsewhere: both are configured widths over the same
    underlying distance, not two different measurements.
    """
    _validate_cells(cells)
    _validate_islets(islets)
    has_immune_class = "immune_class" in cells.columns

    islet_rows = islets[islets["class"] == ISLET_CLASS].copy()
    cells = cells.copy()
    cells["donor_id"] = cells["donor_id"].astype(str)
    islet_rows["donor_id"] = islet_rows["donor_id"].astype(str)
    cells["islet_id"] = cells["islet_id"].astype(str)
    islet_rows["islet_id"] = islet_rows["islet_id"].astype(str)
    d = cells["signed_dist_islet_um"].to_numpy()

    rows = []
    for islet in islet_rows.to_dict("records"):
        mine = (cells["islet_id"].to_numpy() == islet["islet_id"]) & \
               (cells["donor_id"].to_numpy() == str(islet["donor_id"]))
        for lo, hi in config.distance_bins_um:
            in_bin = mine & (d > lo) & (d <= hi)
            n_cells = int(in_bin.sum())
            n_cd45, n_cd45_classified = _classified_positive_count(
                cells.loc[in_bin, "cell_class"], config.cd45_positive_classes)
            n_cd45_val = _count_or_na(n_cd45, n_cd45_classified, n_cells)
            if has_immune_class:
                n_t, n_t_classified = _classified_positive_count(
                    cells.loc[in_bin, "immune_class"], config.t_cell_immune_classes)
                n_t_val = _count_or_na(n_t, n_t_classified, n_cells)
            else:
                n_t_val = pd.NA
            rows.append({
                "donor_id": str(islet["donor_id"]), "islet_id": islet["islet_id"],
                "bin_lo_um": lo, "bin_hi_um": hi,
                "n_cells": n_cells, "n_cd45_positive": n_cd45_val, "n_t_cells": n_t_val,
            })

    columns = ["donor_id", "islet_id", "bin_lo_um", "bin_hi_um",
               "n_cells", "n_cd45_positive", "n_t_cells"]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns)


def non_islet_endocrine_objects(islets: pd.DataFrame) -> pd.DataFrame:
    """EndocrineCluster/EndocrineSingle objects, reported on their own -- never as islets."""
    _validate_islets(islets)
    out = islets[islets["class"].isin(NON_ISLET_ENDOCRINE_CLASSES)].copy()
    return out.reset_index(drop=True)


# ------------------------------------------------------------------ per-donor metrics

def compute_per_donor_metrics(
    per_islet: pd.DataFrame,
    config: IsletImmuneMetricsConfig = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """One row per donor. Every count here is restricted to ``class == "Islet"`` rows --
    ``per_islet`` already carries only those (see :func:`compute_per_islet_metrics`), so no
    further filtering happens here.
    """
    if per_islet.empty:
        return pd.DataFrame(columns=[
            "donor_id", "n_islets", "n_islets_with_cd45_evaluable",
            "pct_islets_with_ge1_cd45_positive", "n_qualifying_islets_npod",
            "insulitis_npod",
        ])

    rows = []
    for donor_id, g in per_islet.groupby("donor_id", sort=True):
        n_islets = len(g)
        evaluable = g["has_ge1_cd45_positive_core"].notna()
        n_evaluable = int(evaluable.sum())
        pct_with_cd45 = (
            float(g.loc[evaluable, "has_ge1_cd45_positive_core"].mean())
            if n_evaluable > 0 else pd.NA
        )

        cd45_counts = g["n_cd45_positive_core"]
        qualifying = cd45_counts.apply(
            lambda v: (v is not pd.NA) and v >= config.npod_min_cd45_per_islet)
        # An islet whose CD45 count is unavailable cannot be confirmed to qualify, but it
        # also cannot be confirmed NOT to -- it is excluded from the numerator, same as every
        # other unavailable-count case in this module, rather than assumed negative.
        n_qualifying = int(qualifying.sum())
        n_cd45_available = int((cd45_counts.apply(lambda v: v is not pd.NA)).sum())
        insulitis = (
            pd.NA if n_cd45_available == 0
            else bool(n_qualifying >= config.npod_min_qualifying_islets)
        )

        rows.append({
            "donor_id": donor_id,
            "n_islets": n_islets,
            "n_islets_with_cd45_evaluable": n_evaluable,
            "pct_islets_with_ge1_cd45_positive": pct_with_cd45,
            "n_qualifying_islets_npod": n_qualifying,
            "insulitis_npod": insulitis,
        })

    return pd.DataFrame(rows)


# ------------------------------------------------------------------ orchestration

@dataclass(frozen=True)
class IsletImmuneMetricsResult:
    per_islet: pd.DataFrame
    distance_bins: pd.DataFrame
    non_islet_endocrine: pd.DataFrame
    per_donor: pd.DataFrame


def compute_islet_immune_metrics(
    cells: pd.DataFrame,
    islets: pd.DataFrame,
    config: IsletImmuneMetricsConfig | None = None,
) -> IsletImmuneMetricsResult:
    cfg = config or DEFAULT_CONFIG
    per_islet = compute_per_islet_metrics(cells, islets, cfg)
    distance_bins = compute_distance_bins(cells, islets, cfg)
    non_islet = non_islet_endocrine_objects(islets)
    per_donor = compute_per_donor_metrics(per_islet, cfg)
    return IsletImmuneMetricsResult(
        per_islet=per_islet, distance_bins=distance_bins,
        non_islet_endocrine=non_islet, per_donor=per_donor,
    )


# ------------------------------------------------------------------ CLI

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", required=True, type=Path, help="typed-cells parquet/CSV")
    ap.add_argument("--islets", required=True, type=Path, help="islet-table parquet/CSV")
    ap.add_argument("--out", required=True, type=Path, help="output directory")
    ap.add_argument("--peri-band-um", type=float, default=DEFAULT_CONFIG.peri_band_um)
    ap.add_argument("--insulitis-t-cell-threshold", type=int,
                    default=DEFAULT_CONFIG.insulitis_t_cell_threshold)
    ap.add_argument("--infiltration-cutoff-um", type=float,
                    default=DEFAULT_CONFIG.infiltration_cutoff_um)
    ap.add_argument("--npod-min-cd45-per-islet", type=int,
                    default=DEFAULT_CONFIG.npod_min_cd45_per_islet)
    ap.add_argument("--npod-min-qualifying-islets", type=int,
                    default=DEFAULT_CONFIG.npod_min_qualifying_islets)
    args = ap.parse_args(argv)

    def _read(path: Path) -> pd.DataFrame:
        return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)

    cells = _read(args.cells)
    islets = _read(args.islets)
    config = replace(
        DEFAULT_CONFIG,
        peri_band_um=args.peri_band_um,
        insulitis_t_cell_threshold=args.insulitis_t_cell_threshold,
        infiltration_cutoff_um=args.infiltration_cutoff_um,
        npod_min_cd45_per_islet=args.npod_min_cd45_per_islet,
        npod_min_qualifying_islets=args.npod_min_qualifying_islets,
    )
    result = compute_islet_immune_metrics(cells, islets, config)

    args.out.mkdir(parents=True, exist_ok=True)
    result.per_islet.to_parquet(args.out / "per_islet.parquet", index=False)
    result.distance_bins.to_parquet(args.out / "distance_bins.parquet", index=False)
    result.non_islet_endocrine.to_parquet(args.out / "non_islet_endocrine.parquet", index=False)
    result.per_donor.to_parquet(args.out / "per_donor.parquet", index=False)
    print(f"[islet_immune_metrics] {len(result.per_islet)} islets, "
          f"{len(result.per_donor)} donors -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
