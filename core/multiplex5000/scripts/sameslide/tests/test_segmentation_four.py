"""segmentation_four end to end on a synthetic preparation and Xenium bundle.

The bundle is written in the vendor's long vertex format, including a cell with two nuclei
(two ``label_id`` rings under one ``cell_id``), so the ring split, the warp into the PhenoCycler
frame, tissue routing and the shared area are all exercised -- not just the metric library.
"""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import pytest
import shapely

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import segmentation_four as sf  # noqa: E402

SHIFT = np.array([1000.0, -500.0])     # Xenium frame -> PhenoCycler frame, a pure translation


class Shift:
    def apply(self, xy):
        return np.asarray(xy, float) + SHIFT


def long_table(ids, polys, labels):
    rows = []
    for cid, p, lab in zip(ids, polys, labels):
        xy = np.asarray(p.exterior.coords)
        rows.append(pd.DataFrame({"cell_id": cid, "vertex_x": xy[:, 0], "vertex_y": xy[:, 1], "label_id": lab}))
    return pd.concat(rows, ignore_index=True)


@pytest.fixture()
def section(tmp_path):
    rng = np.random.default_rng(0)
    centres = np.stack(np.meshgrid(np.arange(30) * 14.0, np.arange(30) * 14.0), -1).reshape(-1, 2)
    n = len(centres)
    # InstanSeg, PhenoCycler frame: cells r=6, nuclei r=3.5. Left half pancreas, right half pLN.
    pc = centres + SHIFT
    cells = [shapely.Point(x, y).buffer(6.0) for x, y in pc]
    nuclei = [shapely.Point(x, y).buffer(3.5) for x, y in pc]
    ids = [f"uuid-{i:04d}" for i in range(n)]
    nuc_wkb = [shapely.to_wkb(g) for g in nuclei]
    nuc_wkb[5] = None                                  # an InstanSeg cell without a nucleus
    prep = tmp_path / "prep"
    (prep / "panc").mkdir(parents=True)
    (prep / "pln_1").mkdir(parents=True)
    pd.DataFrame({"object_id": ids, "wkb": [shapely.to_wkb(g) for g in cells], "nucleus_wkb": nuc_wkb}
                 ).to_parquet(prep / "geometry.parquet")
    left = centres[:, 0] < 200
    pd.DataFrame({"object_id": np.asarray(ids)[left]}).to_parquet(prep / "panc" / "cells.parquet")
    pd.DataFrame({"object_id": np.asarray(ids)[~left]}).to_parquet(prep / "pln_1" / "cells.parquet")

    # Xenium, its own frame, jittered 0.5 um; cell 0 is missing (InstanSeg-only); cell 7 has two
    # nuclei; one extra Xenium-only cell far from any InstanSeg cell but inside the grid.
    xc = centres + rng.normal(0, 0.5, centres.shape)
    keep = np.arange(n) != 0
    xids = [f"x{i}" for i in range(n) if keep[i]]
    xcells = [shapely.Point(x, y).buffer(6.0) for x, y in xc[keep]]
    xnuc = [shapely.Point(x, y).buffer(3.5) for x, y in xc[keep]]
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    nuc_ids, nuc_polys = list(xids), list(xnuc)
    i7 = xids.index("x7")
    c7 = np.asarray(xcells[i7].centroid.coords)[0]
    nuc_polys[i7] = shapely.Point(c7[0] - 2.2, c7[1]).buffer(2.0)
    nuc_ids.insert(i7 + 1, "x7")
    nuc_polys.insert(i7 + 1, shapely.Point(c7[0] + 2.2, c7[1]).buffer(2.0))
    long_table(xids, xcells, range(1, len(xids) + 1)).to_parquet(bundle / "cell_boundaries.parquet")
    long_table(nuc_ids, nuc_polys, range(1, len(nuc_ids) + 1)).to_parquet(bundle / "nucleus_boundaries.parquet")
    methods = np.where(np.arange(len(xids)) % 3 == 0, "Segmented by nucleus expansion of 5.0µm",
                       "Segmented by boundary stain")
    pd.DataFrame({"cell_id": xids, "segmentation_method": methods}).to_parquet(bundle / "cells.parquet")
    sec = {"section_id": "test_panc", "donor_id": "000", "roi": "panc", "xenium_serial": "0000000",
           "bundle": str(bundle)}
    return sec, prep, tmp_path


def test_compare_section_end_to_end(section):
    sec, prep, tmp = section
    out = tmp / "out"
    out.mkdir()
    summary, checks, thresholds = sf.compare_section(sec, prep, Shift(), 0.7, out)
    assert all(all(c.values()) for c in checks.values()), checks
    o = summary["objects"]
    assert o["instanseg_cells"] == 900 and o["instanseg_nuclei"] == 899
    assert o["xenium_cells"] == 899 and o["xenium_nuclei"] == 900
    assert o["xenium_cells_by_nucleus_count"] == {"1": 898, "2": 1}
    # The warp put Xenium back on top of InstanSeg: matched pairs are ~0.5 um apart, not 1118 um.
    cell = summary["levels"]["cell"]["strata"]["all"]
    assert cell["pair_centroid_offset_um"]["distance"]["p50"] < 1.5
    assert cell["headline"]["matched_frac_b"] == pytest.approx(1.0)
    # The two-nucleus Xenium cell splits its InstanSeg nucleus.
    nuc = pd.read_parquet(out / "objects_instanseg_nucleus.parquet")
    assert nuc.loc[nuc.object_id == "uuid-0007", "reason"].item() == "split"
    # Tissue strata come from the preparation's routing, on both sides.
    st = summary["levels"]["cell"]["strata"]
    assert {"pancreas", "pLN"} <= set(st)
    assert st["pancreas"]["n_a"] + st["pLN"]["n_a"] == st["all"]["n_a"]
    assert set(summary["xenium_cells_by_segmentation_method"]) == {
        "Segmented by boundary stain", "Segmented by nucleus expansion of 5.0µm"}
    h = summary["nucleus_in_matched_cells"]
    assert h["same_cell_pair"] == h["matched_nuclei"] - h["no_parent_cell"] - h["a_cell_unmatched"] - h["cells_matched_elsewhere"]
    assert set(thresholds.level) == {"nucleus", "cell"}
    for name in ("objects_instanseg_cell", "objects_xenium_cell", "edges_cell",
                 "objects_instanseg_nucleus", "objects_xenium_nucleus", "edges_nucleus"):
        assert (out / f"{name}.parquet").stat().st_size > 0
    json.dumps(summary, allow_nan=False, default=str)   # the receipt writer refuses NaN


def test_figure_renders(section):
    sec, prep, tmp = section
    out = tmp / "out2"
    out.mkdir()
    summary, _, thresholds = sf.compare_section(sec, prep, Shift(), 0.7, out)
    sf.figure(summary, thresholds, out / "segmentation.png", "test")
    assert (out / "segmentation.png").stat().st_size > 10_000


def test_missing_ceiling_is_null_not_nan(section):
    sec, prep, tmp = section
    out = tmp / "out3"
    out.mkdir()
    summary, checks, _ = sf.compare_section(sec, prep, Shift(), float("nan"), out)
    assert summary["ceiling_shift_um"] is None
    assert all(all(c.values()) for c in checks.values())
    json.dumps(summary, allow_nan=False, default=str)
