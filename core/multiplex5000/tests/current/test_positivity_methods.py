"""The ported per-image method on synthetic sections: levels, references, coverage and the canary tie rule."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from multiplex5000 import positivity
from multiplex5000.positivity import coverage as AC
from multiplex5000.positivity import identity as ID
from multiplex5000.positivity import levels as LV
from multiplex5000.positivity import registry as RG

PANCREAS = ("codex-uniform-32", "pancreas")


def test_the_source_commit_is_pinned():
    assert positivity.SOURCE_COMMIT == "578767eeaf0f95bb1df7272238f361888a0ae75d"


def test_the_pancreas_registry_ships_as_package_data_unchanged():
    reg = RG.load(*PANCREAS)
    assert len(reg.table) == 35 and len(reg.typing_markers()) == 26
    assert reg.const("POSITIVITY_K") == 3 and reg.const("ISOLATION_UM") == 20
    assert (RG.REGISTRY_DIR / "codex-uniform-32__pancreas__labels.csv").is_file()
    assert reg.spec == (*PANCREAS, str(RG.REGISTRY_DIR))


def test_a_registry_with_a_role_but_no_column_on_an_image_is_refused():
    reg = RG.load(*PANCREAS)
    columns = {reg.column(m) for m in reg.measured_markers()} - {reg.column("CD3e")}
    with pytest.raises(RG.RegistryError, match="CD3e"):
        reg.validate({"section": columns})


def test_a_narrow_background_peak_gives_its_level_and_a_plateau_gives_none():
    rng = np.random.default_rng(0)
    narrow = rng.normal(3.0, 0.2, 20000)
    center, width, level = LV.peak_level(narrow, 3.0, lowest=False)
    assert abs(center - 3.0) < 0.06 and 0.1 < width < 0.3
    assert level == pytest.approx(center + 3.0 * width)
    plateau = rng.uniform(1.0, 9.0, 20000)               # no background peak: wider than MAX_WIDTH
    assert np.isnan(LV.peak_level(plateau, 3.0, lowest=False)[2])
    assert np.isnan(LV.peak_level(narrow[:999], 3.0, lowest=False)[2])   # fewer than 1000 values


def test_negative_quantiles_follow_the_method_name():
    assert LV.negative_quantile("negative_p99") == 99.0
    assert LV.negative_quantile("negative_p999") == pytest.approx(99.9)
    assert LV.negative_quantile("peak_dense") is None


def _frame(values):
    return pd.DataFrame({"Cell__CD20": values})


def test_a_known_negative_level_is_the_reference_percentile_and_needs_enough_cells_above():
    reg = RG.load(*PANCREAS)
    rng = np.random.default_rng(1)
    x = np.concatenate([rng.uniform(1, 100, 5000), rng.uniform(500, 1000, 500)])
    ref = np.zeros(len(x), bool)
    ref[:5000] = True                                      # the reference cannot carry the marker
    lv = LV.level_for(_frame(x), reg, "CD20", "whole_cell", "negative_p99", 3.0, np.nan, ref)
    assert lv["level_log"] == pytest.approx(np.log1p(np.percentile(x[:5000], 99)))
    assert lv["pos"][5000:].all() and lv["pos"][:5000].mean() == pytest.approx(0.01, abs=0.002)
    small = np.zeros(len(x), bool)
    small[:500] = True                                     # 1 % of 500 = 5 cells above < KNOWN_NEGATIVE_MIN_ABOVE
    assert np.isnan(LV.level_for(_frame(x), reg, "CD20", "whole_cell", "negative_p99", 3.0, np.nan,
                                 small)["level_log"])


def test_a_marker_without_measurements_is_unavailable_not_negative_by_level():
    reg = RG.load(*PANCREAS)
    lv = LV.level_for(_frame(np.full(5000, np.nan)), reg, "CD20", "whole_cell", "peak_dense", 3.0, np.nan)
    assert np.isnan(lv["level_log"]) and np.isnan(lv["level_raw"])
    assert not lv["pos"].any()                             # no level, no positive cell: the option fails


def test_an_edge_level_is_floored_at_the_empty_slide_null():
    reg = RG.load(*PANCREAS)
    rng = np.random.default_rng(2)
    x = np.expm1(rng.normal(5.0, 0.3, 20000))
    own = LV.level_for(_frame(x).rename(columns={"Cell__CD20": "Cell__ECAD"}), reg, "ECAD", "whole_cell", "edge",
                       3.0, np.nan)["level_log"]
    floored = LV.level_for(_frame(x).rename(columns={"Cell__CD20": "Cell__ECAD"}), reg, "ECAD", "whole_cell", "edge",
                           3.0, 2000.0)["level_log"]
    assert own < np.log1p(2000.0) and floored == pytest.approx(np.log1p(2000.0))


def test_isolated_rates_leave_out_reference_cells_beside_a_positive_cell():
    p = np.array([True, False, False, True])
    ref = np.array([False, True, True, True])
    xy = np.array([[0.0, 0.0], [5.0, 0.0], [100.0, 0.0], [200.0, 0.0]])
    all_rate, iso_rate, n_iso = LV.isolated_rate(p, ref, xy, 20.0, 1)
    assert all_rate == pytest.approx(1 / 3)
    assert n_iso == 2 and iso_rate == pytest.approx(1 / 2)   # cell 1 is beside positive cell 0


def test_a_canary_tie_is_settled_only_by_brighter_and_more_covering_markers():
    detected = np.array([[True, True]] * 3)
    depth = np.array([[2, 2]] * 3)
    score = np.zeros((3, 2))
    bright = np.array([[2.0, 1.0], [2.0, 1.0], [2.0, 1.0]])
    cover = np.array([[0.9, 0.2], [0.30, 0.20], [0.2, 0.9]])   # clear lead, short lead, brighter but covers less
    status, top, _, _, _, _, settled = ID.group_status(detected, depth, score, 1, 2, bright, cover, lead=0.15)
    assert list(status) == ["resolved", "lineage_ambiguous", "lineage_ambiguous"]
    assert list(settled) == [True, False, False] and top[0] == 0


def test_despeckle_is_the_3x3_median_with_replicated_borders():
    a = np.random.default_rng(3).integers(0, 60000, (17, 23)).astype(">u2")
    padded = np.pad(a.astype(np.uint16), 1, mode="edge")
    windows = np.lib.stride_tricks.sliding_window_view(padded, (3, 3))
    expected = np.median(windows.reshape(17, 23, 9), axis=2).astype(np.uint16)
    assert np.array_equal(AC.despeckle(a), expected)
    assert np.array_equal(AC.despeckle(a.astype(np.int32)), expected.astype(np.int32))     # the SciPy fallback
    cv2 = pytest.importorskip("cv2")
    assert np.array_equal(AC.despeckle(a), cv2.medianBlur(a.astype(np.uint16), 3))


def test_coverage_reads_whole_cell_band_and_interior_pixels():
    reg = RG.load(*PANCREAS)
    square = {"type": "Polygon", "coordinates": [[[2, 2], [12, 2], [12, 12], [2, 12], [2, 2]]]}
    features = [{"id": "c1", "geometry": square}]
    image = np.zeros((16, 16), np.uint16)
    image[2:7, 2:12] = 1000                                # the upper half of the cell is bright
    levels = pd.DataFrame({"image": ["s"], "marker": ["CD3e"], "measurement": ["whole_cell"],
                           "level_raw": [500.0]})
    frame = AC.compute("s", features, image.shape, lambda name: image, ["CD3e"], levels, reg)
    mask, _, _, _ = AC.rasterize(features, image.shape)
    cell = mask == 1
    assert frame.loc[0, "cov__CD3e"] == pytest.approx((AC.despeckle(image)[cell] >= 500).mean())
    interior = frame.loc[0, AC.INTERIOR_COUNT]
    assert 0 < interior < cell.sum()
    assert 0 < frame.loc[0, "covint__CD3e"] < 1
    membrane = levels.assign(measurement="membrane")
    band = AC.compute("s", features, image.shape, lambda name: image, ["CD3e"], membrane, reg)
    assert band.loc[0, "cov__CD3e"] != frame.loc[0, "cov__CD3e"]
    assert band.loc[0, "covint__CD3e"] == frame.loc[0, "covint__CD3e"]   # interior does not depend on measurement


def test_coverage_refuses_a_marker_the_image_does_not_have():
    reg = RG.load(*PANCREAS)
    square = {"type": "Polygon", "coordinates": [[[2, 2], [12, 2], [12, 12], [2, 12], [2, 2]]]}
    levels = pd.DataFrame({"image": ["s"], "marker": ["CD3e"], "measurement": ["whole_cell"], "level_raw": [1.0]})
    with pytest.raises(ValueError, match="no image channel"):
        AC.compute("s", [{"id": "c1", "geometry": square}], (16, 16), lambda n: np.zeros((16, 16), np.uint16),
                   ["DAPI"], levels, reg)
