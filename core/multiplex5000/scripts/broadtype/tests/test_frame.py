"""The frozen frame on a synthetic slide: global labels that join ``object_id``, the core rule.

Synthetic geometry only (the ``common/tests/test_window_mask.py`` idiom: square cells written
to a GeoJSON, an expression parquet with centroids in MICRONS, a whole-slide raster from
``redsea.rasterize_mask`` as the 'archive').  ``pins`` is monkeypatched so ``frame.build_donor``
runs its real path -- the DCT hash draw at the frozen 1024/2048 geometry, the windowed
rasterizer, the archive bit-identity gate -- against that slide.

What has to hold:
  * labels are GLOBAL: ``object_ids[label - 1]`` is the cell's UUID in every window;
  * the core rule: a cell is in the frame iff its label appears inside a 1024-px core, so a
    cell in the 512-px pad only is NOT in the frame, and a cell straddling two adjacent cores
    is scored once, in the first core that holds it;
  * frame cells join the expression table on ``object_id`` (GeoJSON-only cells are flagged);
  * ``tails_in_frame`` counts exactly the Step-1 tail cells the frame holds.
"""
from __future__ import annotations

import json
import uuid

import numpy as np
import pandas as pd
import pytest
import tifffile

import frame
import pins
import window_mask
from criteria import CRITERIA
from phenocycler.redsea import rasterize_mask

DONOR = "9999"
PX = pins.PIXEL_SIZE_UM
CORE, WIN = CRITERIA["core_px"], CRITERIA["window_px"]
PAD = (WIN - CORE) // 2
SHAPE = (6144, 7168)
#: 20 grid cores that satisfy the DCT margin rule (room for a 2048 window) -- 16 get drawn.
GRID = [(y0, x0) for y0 in (1024, 2048, 3072, 4096) for x0 in (1024, 2048, 3072, 4096, 5120)]
PER_CORE = 560                     # > regions.MIN_CELLS = 500


def _square(cx, cy, r):
    return [[[cx - r, cy - r], [cx + r, cy - r], [cx + r, cy + r], [cx - r, cy + r], [cx - r, cy - r]]]


@pytest.fixture(scope="module")
def slide(tmp_path_factory):
    """One synthetic donor: GeoJSON, expression parquet, archive mask, and the planted facts."""
    tmp = tmp_path_factory.mktemp("frame")
    rng = np.random.default_rng(0)
    cells = []                     # (cx, cy, r, kind, islet)
    for (y0, x0) in GRID:
        cx = rng.uniform(x0 + 8, x0 + CORE - 8, PER_CORE)
        cy = rng.uniform(y0 + 8, y0 + CORE - 8, PER_CORE)
        islet = rng.random(PER_CORE) < 0.05
        cells += [(float(a), float(b), 2, "interior", bool(i)) for a, b, i in zip(cx, cy, islet)]
        # straddlers on the shared edge with the right-hand and lower neighbour (when in the grid)
        if (y0, x0 + CORE) in GRID:
            cells.append((float(x0 + CORE), float(y0 + 300), 3, "straddler", False))
        if (y0 + CORE, x0) in GRID:
            cells.append((float(x0 + 300), float(y0 + CORE), 3, "straddler", False))
    # pad-only cells: inside the top-row windows' upper pad, outside EVERY core
    for (y0, x0) in [g for g in GRID if g[0] == 1024]:
        for k in range(5):
            cells.append((float(x0 + 100 + 200 * k), float(y0 - 300), 2, "pad_only", False))
    ids = np.array([str(uuid.UUID(bytes=rng.bytes(16), version=4)) for _ in cells], dtype=object)
    feats = [{"type": "Feature", "id": oid, "properties": {},
              "geometry": {"type": "Polygon", "coordinates": _square(cx, cy, r)}}
             for oid, (cx, cy, r, _, _) in zip(ids, cells)]
    gj = tmp / "cells.geojson"
    gj.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))

    kinds = np.array([c[3] for c in cells])
    geojson_only = np.zeros(len(cells), dtype=bool)
    geojson_only[np.flatnonzero(kinds == "interior")[:7]] = True         # 7 cells the run never scored
    expr = pd.DataFrame({
        "object_id": ids[~geojson_only],
        "X_centroid": np.array([c[0] for c in cells])[~geojson_only] * PX,   # MICRONS, as the run stores them
        "Y_centroid": np.array([c[1] for c in cells])[~geojson_only] * PX,
        "cell_region": np.where(np.array([c[4] for c in cells])[~geojson_only], "core", "tissue"),
    })
    ep = tmp / "data.parquet"
    expr.to_parquet(ep, index=False)

    full, oids = rasterize_mask(gj, SHAPE, 1.0)
    assert list(oids) == list(ids)
    archive = tmp / "archive.tif"
    tifffile.imwrite(archive, full.astype(np.int32), compression="deflate", rowsperstrip=2)
    return {"tmp": tmp, "geojson": gj, "expression": ep, "archive": archive, "full": full, "ids": ids,
            "cells": cells, "kinds": kinds, "geojson_only": geojson_only, "expr": expr}


@pytest.fixture(scope="module")
def built(slide, tmp_path_factory):
    """``frame.build_donor`` on the synthetic slide with every pin redirected."""
    out = tmp_path_factory.mktemp("out") / "frame"
    mp = pytest.MonkeyPatch()
    mp.setattr(pins, "archive_mask", lambda d: slide["archive"])
    mp.setattr(pins, "expression_parquet", lambda d: slide["expression"])
    mp.setattr(pins, "geojson_path", lambda d: slide["geojson"])
    mp.setattr(frame, "FRAME", out)
    try:
        res = frame.build_donor(DONOR, verify=True)
    finally:
        mp.undo()
    cells = pd.read_parquet(out / "cells" / f"donor_id={DONOR}.parquet")
    return {"res": res, "out": out, "cells": cells}


def _expected_frame(slide, regions):
    """The core rule computed independently from the whole-slide raster: labels inside each
    1024 core, in region order, a straddler kept by its first core."""
    full, ids = slide["full"], slide["ids"]
    seen, rows = set(), []
    for i, r in enumerate(regions):
        core = full[r["y0"]:r["y0"] + CORE, r["x0"]:r["x0"] + CORE]
        for lab in np.unique(core[core > 0]):
            oid = ids[int(lab) - 1]
            if oid in seen:
                continue
            seen.add(oid)
            rows.append((oid, i))
    return pd.DataFrame(rows, columns=["object_id", "region_idx"])


# --------------------------------------------------------------------------- the draw and the gate
def test_the_draw_is_the_frozen_geometry_and_the_archive_gate_passes(built):
    res = built["res"]
    regs = res["regions"]
    assert res["image_shape"] == list(SHAPE)
    assert len(regs) == CRITERIA["n_blind"] + CRITERIA["n_islet"] == 16
    assert sum(r["draw"] == "blind" for r in regs) == 8 and sum(r["draw"] == "islet" for r in regs) == 8
    assert all(r["core_px"] == CORE and r["n_cells"] >= 500 for r in regs)
    assert {(r["y0"], r["x0"]) for r in regs} <= set(GRID)
    for r in regs:
        assert r["window"] == {"y0": r["y0"] - PAD, "x0": r["x0"] - PAD, "size": WIN}
    assert res["verify"]["all_identical"] is True
    assert len(res["verify"]["windows"]) == 16 and all(w["identical"] for w in res["verify"]["windows"])
    assert all(w["labels_fresh"] == w["labels_archive"] > 0 for w in res["verify"]["windows"])


def test_a_disagreeing_archive_stops_the_lane(slide, tmp_path):
    """The bit-identity gate is a stop, not a tolerance."""
    bad = slide["full"].copy()
    ys, xs = np.nonzero(bad)
    bad[ys[0], xs[0]] = 0                                   # one pixel
    bad_tif = tmp_path / "bad_archive.tif"
    tifffile.imwrite(bad_tif, bad.astype(np.int32), compression="deflate", rowsperstrip=2)
    mp = pytest.MonkeyPatch()
    mp.setattr(pins, "archive_mask", lambda d: bad_tif)
    mp.setattr(pins, "expression_parquet", lambda d: slide["expression"])
    mp.setattr(pins, "geojson_path", lambda d: slide["geojson"])
    mp.setattr(frame, "FRAME", tmp_path / "frame")
    try:
        with pytest.raises(RuntimeError, match="differs from the archived mask -- the lane stops here"):
            frame.build_donor(DONOR, verify=True)
    finally:
        mp.undo()


# --------------------------------------------------------------------------- global labels
def test_window_masks_carry_global_labels_that_index_the_object_ids(built, slide):
    ids = slide["ids"]
    cells = built["cells"]
    for i, r in enumerate(built["res"]["regions"]):
        m = np.load(built["out"] / "masks" / f"donor_id={DONOR}" / f"w{i:02d}.npy")
        assert m.dtype == np.int32 and m.shape == (WIN, WIN)
        labels = np.unique(m[m > 0])
        assert labels.max() <= len(ids)
        assert labels.max() > len(labels), "labels are global and sparse, not renumbered from 1"
        # the window is the archive's window, so every label is a real feature index + 1
        w = r["window"]
        np.testing.assert_array_equal(m, slide["full"][w["y0"]:w["y0"] + WIN, w["x0"]:w["x0"] + WIN])
        core = m[PAD:PAD + CORE, PAD:PAD + CORE]
        in_core = {ids[int(l) - 1] for l in np.unique(core[core > 0])}
        mine = set(cells.loc[cells["region_idx"] == i, "object_id"])
        assert mine <= in_core, "every frame cell of this region is drawn in its core"
        # anything in the core but not this region's must have gone to an EARLIER region
        earlier = set(cells.loc[cells["region_idx"] < i, "object_id"])
        assert in_core - mine <= earlier


# --------------------------------------------------------------------------- the core rule
def test_frame_cells_are_exactly_the_core_rule_applied_to_the_whole_slide(built, slide):
    want = _expected_frame(slide, built["res"]["regions"])
    got = built["cells"][["object_id", "region_idx"]]
    assert not got["object_id"].duplicated().any()
    pd.testing.assert_frame_equal(got.sort_values("object_id").reset_index(drop=True),
                                  want.sort_values("object_id").reset_index(drop=True), check_dtype=False)
    assert built["res"]["n_cells_frame"] == len(got) == len(want)


def test_pad_only_cells_are_in_a_window_but_not_in_the_frame(built, slide):
    kinds, ids = slide["kinds"], slide["ids"]
    pad_ids = set(ids[kinds == "pad_only"])
    frame_ids = set(built["cells"]["object_id"])
    assert not pad_ids & frame_ids
    # and they really are drawn in some window's pad, so the rule was exercised
    drawn = set()
    for i in range(16):
        m = np.load(built["out"] / "masks" / f"donor_id={DONOR}" / f"w{i:02d}.npy")
        drawn |= {ids[int(l) - 1] for l in np.unique(m[m > 0])}
    assert pad_ids & drawn, "no pad-only cell fell inside a drawn window; the fixture is not testing the rule"


def test_a_straddler_is_scored_once_in_its_first_core(built, slide):
    kinds, ids, full = slide["kinds"], slide["ids"], slide["full"]
    regs = built["res"]["regions"]
    cells = built["cells"].set_index("object_id")
    exercised = 0
    for oid in ids[kinds == "straddler"]:
        lab = int(np.flatnonzero(ids == oid)[0]) + 1
        holders = [i for i, r in enumerate(regs)
                   if (full[r["y0"]:r["y0"] + CORE, r["x0"]:r["x0"] + CORE] == lab).any()]
        if len(holders) >= 2:
            exercised += 1
            assert oid in cells.index
            assert cells.loc[oid, "region_idx"] == min(holders)
    assert exercised >= 1, "no straddler sat in two drawn cores; the dedup rule was not exercised"
    assert not built["cells"].index.duplicated().any()


def test_no_frame_cell_lies_outside_its_region_core_beyond_its_own_radius(built, slide):
    """Centroids (from the expression table, in microns) land inside the region's core, up to a
    straddler's 3-px half-width -- the join is by object_id, so this is a check on the geometry."""
    regs = built["res"]["regions"]
    joined = built["cells"].merge(slide["expr"], on="object_id", how="inner")
    assert len(joined) == built["res"]["n_cells_in_expression"]
    x = joined["X_centroid"].to_numpy() / PX
    y = joined["Y_centroid"].to_numpy() / PX
    for i, r in enumerate(regs):
        sel = joined["region_idx"].to_numpy() == i
        assert sel.any()
        assert (x[sel] >= r["x0"] - 3).all() and (x[sel] <= r["x0"] + CORE + 3).all()
        assert (y[sel] >= r["y0"] - 3).all() and (y[sel] <= r["y0"] + CORE + 3).all()


# --------------------------------------------------------------------------- the join
def test_frame_cells_join_the_expression_table_on_object_id(built, slide):
    cells = built["cells"]
    expr_ids = set(slide["expr"]["object_id"])
    assert set(cells.columns) == {"object_id", "region_idx", "draw", "in_expression"}
    assert cells["in_expression"].dtype == bool
    assert (cells["in_expression"] == cells["object_id"].isin(expr_ids)).all()
    res = built["res"]
    assert res["n_cells_in_expression"] + res["n_geojson_only"] == res["n_cells_frame"]
    assert res["n_cells_in_expression"] == int(cells["in_expression"].sum())
    only = set(slide["ids"][slide["geojson_only"]]) & set(cells["object_id"])
    assert only, "no GeoJSON-only cell landed in the frame; the flag was not exercised"
    assert not cells.loc[cells["object_id"].isin(only), "in_expression"].any()
    assert cells["draw"].isin(["blind", "islet"]).all()
    assert set(cells["region_idx"]) == set(range(16))
    # UUIDs, never a renumbering
    assert all(isinstance(o, str) and len(o) == 36 for o in cells["object_id"])


# --------------------------------------------------------------------------- tails_in_frame
def test_tails_in_frame_counts_exactly_the_step1_tail_cells_the_frame_holds(built, slide, tmp_path):
    rng = np.random.default_rng(3)
    expr = slide["expr"]
    frame_ids = set(built["cells"]["object_id"])
    ids = expr["object_id"].to_numpy()
    n = len(ids)
    in_frame = np.array([o in frame_ids for o in ids])
    assert in_frame.any() and (~in_frame).any()
    tab = pd.DataFrame({"donor_id": DONOR, "object_id": ids, "qc_analysis_eligible": True,
                        "X_centroid": expr["X_centroid"], "Y_centroid": expr["Y_centroid"],
                        "cell_region": expr["cell_region"], "islet_num": -1})
    # CD3e: a continuous marker whose 20 loudest cells are 12 in-frame + 8 out-of-frame, planted
    v = rng.normal(size=n)
    loud = np.concatenate([rng.choice(np.flatnonzero(in_frame), 12, replace=False),
                           rng.choice(np.flatnonzero(~in_frame), 8, replace=False)])
    v[loud] += 50 + rng.random(20)
    tab["CD3e"] = v
    # INS: 30 % NaN, tail planted entirely in the frame
    w = rng.normal(size=n)
    w[rng.choice(n, int(0.3 * n), replace=False)] = np.nan
    fin_frame = np.flatnonzero(in_frame & np.isfinite(w))
    w[rng.choice(fin_frame, 25, replace=False)] += 50
    tab["INS"] = w
    tab["SST"] = np.nan                                            # never acquired
    up = tmp_path / "uniform.parquet"
    tab.to_parquet(up, index=False)

    mp = pytest.MonkeyPatch()
    mp.setattr(pins, "uniform_parquet", lambda d: up)
    mp.setattr(frame, "FRAME", built["out"])
    try:
        df = frame.tails_in_frame([DONOR])
    finally:
        mp.undo()
    assert list(df.columns) == ["donor_id", "marker", "n_finite", "k_tail", "n_tail_in_frame", "judgeable_in_frame"]
    assert df["marker"].tolist() == ["CD3e", "INS", "SST"], "context columns are never markers"
    by = df.set_index("marker")
    k = max(CRITERIA["tail_min_k"], int(CRITERIA["tail_q"] * n))
    assert k == 20
    # independent expectation
    for m in ("CD3e", "INS"):
        vals = tab[m].to_numpy(float)
        fin = np.flatnonzero(np.isfinite(vals))
        top = fin[np.argsort(-vals[fin])[:k]]
        expect = int(in_frame[top].sum())
        assert by.loc[m, "n_finite"] == len(fin) and by.loc[m, "k_tail"] == k
        assert by.loc[m, "n_tail_in_frame"] == expect
        assert bool(by.loc[m, "judgeable_in_frame"]) == (expect >= CRITERIA["min_tail_in_pool"])
    assert by.loc["CD3e", "n_tail_in_frame"] == 12 and not by.loc["CD3e", "judgeable_in_frame"]
    assert by.loc["INS", "n_tail_in_frame"] == 20 and by.loc["INS", "judgeable_in_frame"]
    assert by.loc["SST"].tolist()[1:] == [0, 0, 0, False]
    written = pd.read_csv(built["out"] / "tails_in_frame.csv", dtype={"donor_id": str})
    pd.testing.assert_frame_equal(written, df.reset_index(drop=True), check_dtype=False)


# --------------------------------------------------------------------------- windows_for
def test_windows_for_pads_the_core_symmetrically():
    from types import SimpleNamespace
    regs = [SimpleNamespace(y0=4096, x0=2048)]
    (w,) = frame.windows_for(regs)
    assert isinstance(w, window_mask.Window)
    assert (w.y0, w.x0, w.size) == (4096 - PAD, 2048 - PAD, WIN)
    assert w.y1 - w.y0 == WIN and w.y0 + PAD == 4096
