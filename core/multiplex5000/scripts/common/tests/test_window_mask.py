"""Windowed rasterization must equal the whole-slide rasterization it replaces.

Synthetic geometry only -- no cohort data.  The reference is ``redsea.rasterize_mask``
itself: whatever it draws for the whole slide, a window must reproduce exactly inside its
own bounds.  The invariant that matters most is the label one -- ``label = feature_index+1``
over the WHOLE stream -- because a window that renumbered from 1 would join every cell in
it to the wrong ``object_id`` with every other check still green.

(The same claim is checked against real archived masks by
``data/development/deepcell_types_review/window_mask_verification.json``: 32 windows across
four donors, 93,834 labels, all bit-identical.)
"""
from __future__ import annotations

import json
import pathlib
import sys
import uuid

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from phenocycler.redsea import rasterize_mask  # noqa: E402
from window_mask import Window, rasterize_windows  # noqa: E402

SHAPE = (600, 600)


def _square(cx, cy, r):
    return [[[cx - r, cy - r], [cx + r, cy - r], [cx + r, cy + r], [cx - r, cy + r], [cx - r, cy - r]]]


def _geojson(tmp_path, cells, nuclei=True):
    feats = []
    for cx, cy, r in cells:
        f = {"type": "Feature", "id": str(uuid.uuid4()),
             "geometry": {"type": "Polygon", "coordinates": _square(cx, cy, r)}}
        if nuclei:
            f["nucleusGeometry"] = {"type": "Polygon", "coordinates": _square(cx, cy, max(r // 2, 1))}
        feats.append(f)
    p = tmp_path / "cells.geojson"
    p.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    return p


CELLS = [(80, 80, 12), (140, 90, 10), (300, 300, 15), (310, 330, 9),
         (500, 120, 11), (128, 300, 14), (255, 255, 20)]


def test_window_reproduces_the_whole_slide_raster_exactly(tmp_path):
    gj = _geojson(tmp_path, CELLS)
    full, full_oids = rasterize_mask(gj, SHAPE, 1.0)
    wins = [Window(0, 0, 200), Window(256, 256, 128), Window(64, 256, 200)]
    masks, oids = rasterize_windows(gj, wins)
    assert oids == full_oids
    for w, got in zip(wins, masks):
        np.testing.assert_array_equal(full[w.y0:w.y1, w.x0:w.x1], got)


def test_labels_stay_global_and_index_the_object_ids(tmp_path):
    """A per-window renumber would mislabel every cell with every other gate green."""
    gj = _geojson(tmp_path, CELLS)
    masks, oids = rasterize_windows(gj, [Window(256, 256, 128)])
    labels = np.unique(masks[0][masks[0] > 0])
    assert labels.size >= 1
    assert labels.max() <= len(oids)
    # the 255,255 cell is feature index 6 -> label 7, and it must be that label here
    assert 7 in labels, f"expected the global label 7, got {labels}"
    assert all(isinstance(oids[l - 1], str) for l in labels)


def test_a_cell_straddling_the_edge_is_clipped_not_dropped(tmp_path):
    gj = _geojson(tmp_path, [(200, 200, 30)])
    full, _ = rasterize_mask(gj, SHAPE, 1.0)
    w = Window(180, 180, 64)                      # cuts the cell in both axes
    masks, _ = rasterize_windows(gj, [w])
    np.testing.assert_array_equal(full[w.y0:w.y1, w.x0:w.x1], masks[0])
    assert (masks[0] > 0).any(), "the clipped part of the cell vanished"


def test_nucleus_is_a_strict_subset_of_its_own_cell(tmp_path):
    gj = _geojson(tmp_path, CELLS)
    masks, _, nucs = rasterize_windows(gj, [Window(0, 0, 400)], with_nucleus=True)
    m, n = masks[0], nucs[0]
    assert (n[n > 0] == m[n > 0]).all(), "a nucleus pixel is claimed by a different cell"
    assert (n > 0).sum() < (m > 0).sum()


def test_nucleus_matches_the_whole_slide_nucleus_raster(tmp_path):
    gj = _geojson(tmp_path, CELLS)
    full, _, full_nuc = rasterize_mask(gj, SHAPE, 1.0, with_nucleus=True)
    w = Window(64, 64, 256)
    _, _, nucs = rasterize_windows(gj, [w], with_nucleus=True)
    np.testing.assert_array_equal(full_nuc[w.y0:w.y1, w.x0:w.x1], nucs[0])


def test_many_windows_in_one_pass_agree_with_one_at_a_time(tmp_path):
    gj = _geojson(tmp_path, CELLS)
    wins = [Window(0, 0, 200), Window(256, 256, 128), Window(448, 64, 128)]
    together, _ = rasterize_windows(gj, wins)
    for w, got in zip(wins, together):
        alone, _ = rasterize_windows(gj, [w])
        np.testing.assert_array_equal(alone[0], got)


def test_an_empty_request_raises(tmp_path):
    gj = _geojson(tmp_path, CELLS)
    with pytest.raises(ValueError, match="no windows"):
        rasterize_windows(gj, [])


def test_a_degenerate_window_raises(tmp_path):
    gj = _geojson(tmp_path, CELLS)
    with pytest.raises(ValueError, match="non-positive size"):
        rasterize_windows(gj, [Window(0, 0, 0)])


def test_a_window_over_empty_tissue_is_empty_not_wrong(tmp_path):
    gj = _geojson(tmp_path, CELLS)
    masks, oids = rasterize_windows(gj, [Window(560, 560, 32)])
    assert masks[0].max() == 0
    assert len(oids) == len(CELLS), "the pass still saw every feature"
