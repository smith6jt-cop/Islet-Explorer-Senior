"""Region selection must be blind, reproducible, and unable to quietly return fewer.

The failure this guards is subtle: a draw that silently shrinks when a donor is awkward
would let a thin donor into the gate looking like a full one.  A short draw is a finding
about the donor, so it raises.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from regions import CORE_PX, MIN_CELLS, region_hash, select_regions  # noqa: E402

PX = 0.5077663810243286
SHAPE = (23760, 24960)


def _parquet(tmp_path, n=40000, seed=0, islet_frac=0.05, span=9000):
    """Synthetic centroids in MICRONS, the unit the real expression parquet uses.

    Sized so a 1024-px core holds ~1,200 cells -- the real donors run 784-3,886 -- because
    a fixture below the 500-cell floor tests the floor rather than the draw.
    """
    rng = np.random.default_rng(seed)
    xp = rng.uniform(3000, span, n)
    yp = rng.uniform(3000, span, n)
    region = np.where(rng.random(n) < islet_frac, "core", "tissue")
    df = pd.DataFrame({"X_centroid": xp * PX, "Y_centroid": yp * PX,
                       "cell_region": region, "object_id": np.arange(n).astype(str)})
    f = tmp_path / "data.parquet"
    df.to_parquet(f)
    return f


def test_hash_is_stable_and_donor_specific():
    assert region_hash("6539", 1024, 2048) == region_hash("6539", 1024, 2048)
    assert region_hash("6539", 1024, 2048) != region_hash("6593", 1024, 2048)
    assert region_hash("6539", 1024, 2048) != region_hash("6539", 2048, 1024)


def test_selection_is_reproducible_from_the_donor_id_alone(tmp_path):
    f = _parquet(tmp_path)
    a = select_regions("6539", f, PX, SHAPE)
    b = select_regions("6539", f, PX, SHAPE)
    assert [(r.y0, r.x0, r.draw) for r in a] == [(r.y0, r.x0, r.draw) for r in b]


def test_a_different_donor_id_draws_different_regions(tmp_path):
    f = _parquet(tmp_path)
    a = {(r.y0, r.x0) for r in select_regions("6539", f, PX, SHAPE)}
    b = {(r.y0, r.x0) for r in select_regions("6593", f, PX, SHAPE)}
    assert a != b


def test_eight_distinct_regions_four_of_each_draw(tmp_path):
    regs = select_regions("6539", _parquet(tmp_path), PX, SHAPE)
    assert len(regs) == 8
    assert len({(r.y0, r.x0) for r in regs}) == 8, "a region was drawn twice"
    assert sum(r.draw == "blind" for r in regs) == 4
    assert sum(r.draw == "islet" for r in regs) == 4


def test_islet_draws_actually_contain_an_islet(tmp_path):
    regs = select_regions("6539", _parquet(tmp_path), PX, SHAPE)
    for r in regs:
        if r.draw == "islet":
            assert r.n_islet_core_cells >= 1


def test_every_region_meets_the_cell_floor(tmp_path):
    regs = select_regions("6539", _parquet(tmp_path), PX, SHAPE)
    assert all(r.n_cells >= MIN_CELLS for r in regs)


def test_regions_leave_room_for_the_largest_window(tmp_path):
    """A window-size arm must never silently read a clipped field."""
    regs = select_regions("6539", _parquet(tmp_path), PX, SHAPE)
    margin = (4096 - CORE_PX) // 2
    H, W = SHAPE
    for r in regs:
        assert r.y0 >= margin and r.x0 >= margin
        assert r.y0 + CORE_PX + margin <= H and r.x0 + CORE_PX + margin <= W


def test_a_short_draw_raises_rather_than_returning_fewer(tmp_path):
    f = _parquet(tmp_path, n=900, islet_frac=0.0, span=6000)
    with pytest.raises(ValueError, match="regions available"):
        select_regions("6539", f, PX, SHAPE)


def test_centroids_are_read_as_microns(tmp_path):
    """Pinned: taking them as pixels puts every region in the wrong place.

    Verified against the archived 6539 mask on 2026-08-24 -- 98.2 % of a window's cells
    fall inside it when centroids are divided by the pixel size, 0 % when they are not.
    """
    f = _parquet(tmp_path)
    um = select_regions("6539", f, PX, SHAPE)
    px = select_regions("6539", f, 1.0, SHAPE)      # pretend they were already pixels
    assert {(r.y0, r.x0) for r in um} != {(r.y0, r.x0) for r in px}
    assert max(r.y0 for r in um) > max(r.y0 for r in px)


def test_frozen_artifact_matches_its_own_rule():
    """regions.json must be what the code produces -- not a hand-edited list of origins."""
    f = HERE / "regions.json"
    if not f.exists():
        pytest.skip("regions.json not built in this checkout")
    d = json.loads(f.read_text())
    assert d["core_px"] == CORE_PX and d["min_cells"] == MIN_CELLS
    assert len(d["regions"]) == 8, "the gate cohort is eight donors"
    for donor, regs in d["regions"].items():
        assert len(regs) == 8, donor
        assert len({(r["y0"], r["x0"]) for r in regs}) == 8, donor
        for r in regs:
            assert r["digest"] == region_hash(donor, r["y0"], r["x0"]), (donor, r)
            assert r["n_cells"] >= MIN_CELLS
            if r["draw"] == "islet":
                assert r["n_islet_core_cells"] >= 1
