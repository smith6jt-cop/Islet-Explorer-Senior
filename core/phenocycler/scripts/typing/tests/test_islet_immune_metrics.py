"""Synthetic-islet tests for ``islet_immune_metrics.py``.

Every islet here is hand-built with known cell placements so each metric can be checked
against an arithmetic expectation, not just "did it run". No data/ dependency.

    ../../../../.venv/bin/python -m pytest -q core/phenocycler/scripts/typing/tests
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve()
LANE = HERE.parents[1]
sys.path.insert(0, str(LANE))

import islet_immune_metrics as M  # noqa: E402


def _cell(donor_id, object_id, islet_id, signed_dist_islet_um, cell_class=None,
         immune_class=None, x_um=0.0, y_um=0.0):
    return {
        "donor_id": donor_id, "object_id": object_id, "x_um": x_um, "y_um": y_um,
        "cell_class": cell_class, "immune_class": immune_class,
        "cell_region": "core" if signed_dist_islet_um <= 0 else "peri",
        "islet_id": islet_id, "signed_dist_islet_um": signed_dist_islet_um,
    }


def _islet(donor_id, islet_id, klass, area_um2=1000.0, caliper_um=40.0):
    return {"donor_id": donor_id, "islet_id": islet_id, "class": klass,
            "area_um2": area_um2, "mean_caliper_diameter_um": caliper_um}


# --------------------------------------------------------------------- contract

def test_missing_required_cell_column_raises():
    cells = pd.DataFrame({"donor_id": ["D1"], "object_id": ["o1"]})
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    with pytest.raises(M.IsletImmuneMetricsError, match="cells table"):
        M.compute_per_islet_metrics(cells, islets)


def test_missing_required_islet_column_raises():
    cells = pd.DataFrame([_cell("D1", "o1", "isl-1", -1.0, "Immune")])
    islets = pd.DataFrame({"donor_id": ["D1"], "islet_id": ["isl-1"]})
    with pytest.raises(M.IsletImmuneMetricsError, match="islets table"):
        M.compute_per_islet_metrics(cells, islets)


def test_immune_class_and_cell_region_are_optional():
    """A minimal cells table with no immune_class/cell_region column does not raise."""
    cells = pd.DataFrame([{
        "donor_id": "D1", "object_id": "o1", "x_um": 0.0, "y_um": 0.0,
        "cell_class": "Immune", "islet_id": "isl-1", "signed_dist_islet_um": -1.0,
    }])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    result = M.compute_islet_immune_metrics(cells, islets)
    row = result.per_islet.iloc[0]
    assert row["n_cd45_positive_core"] == 1
    assert row["n_t_cells_core_peri"] is pd.NA          # no immune_class column at all


# --------------------------------------------------------------------- per-islet: CD45/core/peri

def test_cd45_core_and_peri_counts():
    cells = pd.DataFrame([
        _cell("D1", "o1", "isl-1", -5.0, "Immune"),
        _cell("D1", "o2", "isl-1", -2.0, "Immune"),
        _cell("D1", "o3", "isl-1", -1.0, "Beta"),
        _cell("D1", "o4", "isl-1", 10.0, "Immune"),     # peri (<=20um band)
        _cell("D1", "o5", "isl-1", 15.0, "Exocrine"),   # peri
        _cell("D1", "o6", "isl-1", 25.0, "Immune"),     # beyond the 20um band -> not core/peri
    ])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    per_islet = M.compute_per_islet_metrics(cells, islets)
    row = per_islet.iloc[0]

    assert row["n_cells_core"] == 3
    assert row["n_cd45_positive_core"] == 2
    assert row["n_cd45_positive_peri"] == 1
    assert row["cd45_per_nucleus_core"] == pytest.approx(2 / 3)
    assert bool(row["has_ge1_cd45_positive_core"]) is True


def test_peri_band_is_parameterizable():
    cells = pd.DataFrame([_cell("D1", "o1", "isl-1", 25.0, "Immune")])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])

    default_row = M.compute_per_islet_metrics(cells, islets).iloc[0]
    assert default_row["n_cd45_positive_peri"] == 0           # 25um is outside the 20um default

    wide = M.IsletImmuneMetricsConfig(peri_band_um=30.0)
    wide_row = M.compute_per_islet_metrics(cells, islets, wide).iloc[0]
    assert wide_row["n_cd45_positive_peri"] == 1


# --------------------------------------------------------------------- missing != zero

def test_unclassified_cells_report_unavailable_not_zero():
    """Every cell in the core has a null cell_class -- there is no evidence, so CD45 counts
    for this islet must come back NA, not a misleadingly confident 0."""
    cells = pd.DataFrame([
        _cell("D1", "o1", "isl-1", -1.0, None),
        _cell("D1", "o2", "isl-1", -2.0, None),
    ])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    row = M.compute_per_islet_metrics(cells, islets).iloc[0]

    assert row["n_cd45_positive_core"] is pd.NA
    assert row["cd45_per_nucleus_core"] is pd.NA
    assert row["has_ge1_cd45_positive_core"] is pd.NA
    assert row["n_cells_core"] == 2                     # the cells themselves are still counted


def test_a_genuine_zero_is_distinguishable_from_unavailable():
    """Classified cells, none of them CD45+, is a real 0 -- not NA."""
    cells = pd.DataFrame([
        _cell("D1", "o1", "isl-1", -1.0, "Beta"),
        _cell("D1", "o2", "isl-1", -2.0, "Exocrine"),
    ])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    row = M.compute_per_islet_metrics(cells, islets).iloc[0]

    assert row["n_cd45_positive_core"] == 0
    assert row["cd45_per_nucleus_core"] == 0.0
    assert bool(row["has_ge1_cd45_positive_core"]) is False


# --------------------------------------------------------------------- T cells / insulitis

def test_insulitic_flag_uses_the_configurable_threshold():
    t_cells = [_cell("D1", f"t{i}", "isl-1", -1.0, "Immune", "Imm-T-cell") for i in range(6)]
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    cells = pd.DataFrame(t_cells)

    default_row = M.compute_per_islet_metrics(cells, islets).iloc[0]
    assert default_row["n_t_cells_core_peri"] == 6
    assert bool(default_row["insulitic"]) is True

    strict = M.IsletImmuneMetricsConfig(insulitis_t_cell_threshold=7)
    strict_row = M.compute_per_islet_metrics(cells, islets, strict).iloc[0]
    assert bool(strict_row["insulitic"]) is False


def test_t_cells_count_across_core_and_peri():
    cells = pd.DataFrame([
        _cell("D1", "o1", "isl-1", -1.0, "Immune", "Imm-T-cell"),
        _cell("D1", "o2", "isl-1", 5.0, "Immune", "Imm-T-cell"),
        _cell("D1", "o3", "isl-1", -1.0, "Immune", "Imm-B-cell"),   # not a T cell
    ])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    row = M.compute_per_islet_metrics(cells, islets).iloc[0]
    assert row["n_t_cells_core_peri"] == 2


# --------------------------------------------------------------------- infiltration score

def test_infiltration_score_is_one_inside_and_one_over_d_outside():
    cells = pd.DataFrame([
        _cell("D1", "o1", "isl-1", -3.0, "Immune"),     # inside -> weight 1
        _cell("D1", "o2", "isl-1", 25.0, "Immune"),      # outside -> weight 1/25
        _cell("D1", "o3", "isl-1", -1.0, "Beta"),        # not CD45+, contributes nothing
        _cell("D1", "o4", "isl-1", 150.0, "Immune"),     # beyond the 100um default cutoff
    ])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    row = M.compute_per_islet_metrics(cells, islets).iloc[0]
    assert row["infiltration_score"] == pytest.approx(1.0 + 1.0 / 25.0)


def test_infiltration_cutoff_is_parameterizable():
    cells = pd.DataFrame([_cell("D1", "o1", "isl-1", 150.0, "Immune")])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])

    default_row = M.compute_per_islet_metrics(cells, islets).iloc[0]
    # No cell at all falls within the 100um default cutoff -- a real, confident 0, not NA.
    assert default_row["infiltration_score"] == 0.0

    wide = M.IsletImmuneMetricsConfig(infiltration_cutoff_um=200.0)
    wide_row = M.compute_per_islet_metrics(cells, islets, wide).iloc[0]
    assert wide_row["infiltration_score"] == pytest.approx(1.0 / 150.0)


def test_infiltration_score_is_unavailable_when_nothing_in_range_is_classified():
    cells = pd.DataFrame([_cell("D1", "o1", "isl-1", 10.0, None)])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    row = M.compute_per_islet_metrics(cells, islets).iloc[0]
    assert row["infiltration_score"] is pd.NA


# --------------------------------------------------------------------- distance bins

def test_distance_bin_counts():
    cells = pd.DataFrame([
        _cell("D1", "o1", "isl-1", 10.0, "Immune"),   # 0-20
        _cell("D1", "o2", "isl-1", 30.0, "Immune"),   # 20-40
        _cell("D1", "o3", "isl-1", 55.0, "Beta"),     # 40-60, not CD45+
        _cell("D1", "o4", "isl-1", 90.0, "Immune"),   # 60-100
        _cell("D1", "o5", "isl-1", 500.0, "Immune"),  # beyond every bin
    ])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    bins = M.compute_distance_bins(cells, islets)

    assert len(bins) == 4
    by_range = {(r.bin_lo_um, r.bin_hi_um): r for r in bins.itertuples()}
    assert by_range[(0.0, 20.0)].n_cells == 1
    assert by_range[(0.0, 20.0)].n_cd45_positive == 1
    assert by_range[(20.0, 40.0)].n_cd45_positive == 1
    assert by_range[(40.0, 60.0)].n_cells == 1
    assert by_range[(40.0, 60.0)].n_cd45_positive == 0
    assert by_range[(60.0, 100.0)].n_cd45_positive == 1


def test_distance_bins_are_parameterizable():
    cells = pd.DataFrame([_cell("D1", "o1", "isl-1", 5.0, "Immune")])
    islets = pd.DataFrame([_islet("D1", "isl-1", "Islet")])
    config = M.IsletImmuneMetricsConfig(distance_bins_um=((0.0, 10.0),))
    bins = M.compute_distance_bins(cells, islets, config)
    assert len(bins) == 1
    assert bins.iloc[0]["bin_lo_um"] == 0.0 and bins.iloc[0]["bin_hi_um"] == 10.0


# --------------------------------------------------------------------- Islet vs EndocrineCluster/Single

def test_only_class_islet_counts_as_an_islet():
    cells = pd.DataFrame([
        _cell("D1", "o1", "isl-1", -1.0, "Immune"),
        _cell("D1", "o2", "cluster-1", -1.0, "Immune"),
        _cell("D1", "o3", "single-1", -1.0, "Immune"),
    ])
    islets = pd.DataFrame([
        _islet("D1", "isl-1", "Islet"),
        _islet("D1", "cluster-1", "EndocrineCluster", area_um2=100.0, caliper_um=12.0),
        _islet("D1", "single-1", "EndocrineSingle", area_um2=20.0, caliper_um=6.0),
    ])
    result = M.compute_islet_immune_metrics(cells, islets)

    assert list(result.per_islet["islet_id"]) == ["isl-1"]
    assert set(result.non_islet_endocrine["islet_id"]) == {"cluster-1", "single-1"}
    assert set(result.non_islet_endocrine["class"]) == {"EndocrineCluster", "EndocrineSingle"}
    # Never folded into an islet's own counts:
    assert result.per_islet.iloc[0]["n_cells_core"] == 1


# --------------------------------------------------------------------- per-donor

def test_pct_islets_with_cd45_and_npod_insulitis():
    cells = []
    # isl-1: 20 CD45+ cells (qualifies for nPOD's per-islet >=15 threshold)
    cells += [_cell("D1", f"a{i}", "isl-1", -1.0, "Immune") for i in range(20)]
    # isl-2: 16 CD45+ cells (also qualifies)
    cells += [_cell("D1", f"b{i}", "isl-2", -1.0, "Immune") for i in range(16)]
    # isl-3: no CD45+ cells at all, but classified (a real zero, still "evaluable")
    cells += [_cell("D1", "c0", "isl-3", -1.0, "Beta")]
    islets = pd.DataFrame([
        _islet("D1", "isl-1", "Islet"), _islet("D1", "isl-2", "Islet"),
        _islet("D1", "isl-3", "Islet"),
    ])
    cells = pd.DataFrame(cells)

    result = M.compute_islet_immune_metrics(cells, islets)
    donor_row = result.per_donor.iloc[0]

    assert donor_row["n_islets"] == 3
    assert donor_row["pct_islets_with_ge1_cd45_positive"] == pytest.approx(2 / 3)
    # Default nPOD needs >=3 qualifying islets; only 2 qualify here.
    assert donor_row["n_qualifying_islets_npod"] == 2
    assert bool(donor_row["insulitis_npod"]) is False

    # Lower the bar and the same data now qualifies.
    lenient = M.IsletImmuneMetricsConfig(npod_min_cd45_per_islet=15, npod_min_qualifying_islets=2)
    result2 = M.compute_islet_immune_metrics(cells, islets, lenient)
    assert bool(result2.per_donor.iloc[0]["insulitis_npod"]) is True


def test_donor_metrics_restricted_to_islets_never_endocrine_cluster_or_single():
    cells = pd.DataFrame([
        _cell("D1", "o1", "isl-1", -1.0, "Immune"),
        _cell("D1", "o2", "cluster-1", -1.0, "Immune"),
    ])
    islets = pd.DataFrame([
        _islet("D1", "isl-1", "Islet"),
        _islet("D1", "cluster-1", "EndocrineCluster", area_um2=100.0, caliper_um=12.0),
    ])
    result = M.compute_islet_immune_metrics(cells, islets)
    assert result.per_donor.iloc[0]["n_islets"] == 1


def test_empty_islets_table_yields_empty_but_well_formed_outputs():
    cells = pd.DataFrame([_cell("D1", "o1", "isl-1", -1.0, "Immune")])
    islets = pd.DataFrame(columns=list(M.ISLETS_REQUIRED_COLUMNS))
    result = M.compute_islet_immune_metrics(cells, islets)
    assert result.per_islet.empty
    assert result.per_donor.empty
    assert list(result.per_islet.columns) == [
        "donor_id", "islet_id", "class", "area_um2", "mean_caliper_diameter_um",
        "n_cells_core", "n_cd45_positive_core", "n_cd45_positive_peri",
        "cd45_per_nucleus_core", "has_ge1_cd45_positive_core",
        "n_t_cells_core_peri", "insulitic", "infiltration_score",
    ]
