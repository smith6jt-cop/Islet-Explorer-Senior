"""Synthetic checks for segcompare: every expected number is derived by hand, not by the module."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import shapely

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import segcompare as sc  # noqa: E402


def square(x, y, s=10.0):
    return shapely.box(x, y, x + s, y + s)


def grid(n, step=20.0, s=10.0, dx=0.0, dy=0.0):
    return np.array([square(i * step + dx, j * step + dy, s) for i in range(n) for j in range(n)],
                    dtype=object)


def disk_iou(r, d):
    """IoU of two disks of radius r whose centres are d apart (exact)."""
    inter = 2 * r * r * np.arccos(d / (2 * r)) - (d / 2) * np.sqrt(4 * r * r - d * d)
    return inter / (2 * np.pi * r * r - inter)


def test_identical_segmentations_match_everything():
    a = grid(4)
    out = sc.compare(a, a.copy())
    s = out["summary"]["strata"]["all"]
    assert s["n_a"] == s["n_b"] == 16
    assert out["thresholds"]["f1"].eq(1.0).all()
    assert s["headline"]["SQ"] == pytest.approx(1.0)
    assert s["headline"]["PQ"] == pytest.approx(1.0)
    assert s["reason_a"] == {"matched": 16} and s["reason_b"] == {"matched": 16}
    assert s["a_split_by_b"]["objects"] == 0 and s["b_split_by_a"]["objects"] == 0


def test_translation_gives_the_analytic_iou_and_threshold_cut():
    # Squares of side 10 shifted 1 along x: IoU = (s - d) / (s + d) = 9 / 11.
    a = grid(3)
    b = grid(3, dx=1.0)
    out = sc.compare(a, b)
    assert np.allclose(out["objects_a"]["pair_iou"], 9 / 11)
    t = out["thresholds"].set_index("iou_threshold")
    assert t.loc[0.8, "matched_a"] == 9 and t.loc[0.9, "matched_a"] == 0
    off = out["summary"]["strata"]["all"]["pair_centroid_offset_um"]
    assert off["mean_dx"] == pytest.approx(1.0) and off["mean_dy"] == pytest.approx(0.0, abs=1e-12)


def test_ceiling_matches_the_disk_formula_in_any_direction():
    r, d = 3.5, 0.7
    a = np.array([shapely.Point(x, 0).buffer(r, quad_segs=256) for x in range(0, 200, 20)],
                 dtype=object)
    got = sc.shift_ceiling(a, d, seed=3)
    assert np.allclose(got, disk_iou(r, d), atol=2e-3)


def test_split_and_merge_are_named_from_both_sides():
    # One 30x10 InstanSeg object; Xenium cuts it into three 10x10 objects: IoU 1/3 each.
    a = np.array([shapely.box(0, 0, 30, 10)], dtype=object)
    b = np.array([square(0, 0), square(10, 0), square(20, 0)], dtype=object)
    out = sc.compare(a, b)
    oa, ob = out["objects_a"], out["objects_b"]
    assert oa.loc[0, "reason"] == "split" and oa.loc[0, "holds_n"] == 3
    assert (ob["reason"] == "merged").all() and (ob["held_by"] == 0).all()
    s = out["summary"]["strata"]["all"]
    assert s["a_split_by_b"] == {"objects": 1, "b_objects_involved": 3}
    assert np.allclose(out["edges"]["iou"], 1 / 3)


def test_unmatched_without_any_overlap_is_a_detection_difference():
    a = grid(2)
    b = np.concatenate([grid(2), [square(500, 500)]])
    out = sc.compare(a, b)
    ob = out["objects_b"]
    assert ob["reason"].tolist() == ["matched"] * 4 + ["no_overlap"]
    h = out["summary"]["strata"]["all"]["headline"]
    assert h["matched_frac_a"] == 1.0 and h["matched_frac_b"] == pytest.approx(4 / 5)
    assert h["RQ"] == pytest.approx(8 / 9)


def test_half_offset_is_a_boundary_disagreement():
    # Shift by half a side: IoU 5/15, each covers exactly half of the other -- not a majority.
    out = sc.compare(np.array([square(0, 0)], dtype=object), np.array([square(5, 0)], dtype=object))
    assert out["objects_a"].loc[0, "reason"] == "boundary"
    assert out["objects_b"].loc[0, "reason"] == "boundary"
    assert out["objects_a"].loc[0, "held_by"] == -1


def test_touching_objects_are_not_edges():
    g = sc.overlap_graph(np.array([square(0, 0)], dtype=object), np.array([square(10, 0)], dtype=object))
    assert g.empty


def test_match_takes_the_most_pairs_not_the_greedy_best():
    # Greedy would take a0-b0 (0.6) and strand a1; the optimum pairs a0-b1 and a1-b0.
    edges = pd.DataFrame({"a": [0, 0, 1], "b": [0, 1, 0], "iou": [0.6, 0.55, 0.52]})
    pa, pb = sc.match(edges, 2, 2, 0.5)
    assert pa.tolist() == [1, 0] and pb.tolist() == [1, 0]


def test_match_breaks_cardinality_ties_by_total_iou():
    # Both perfect matchings have two pairs; a0-b0 + a1-b1 sums 1.5, the cross pairs sum 1.3.
    edges = pd.DataFrame({"a": [0, 0, 1, 1], "b": [0, 1, 0, 1], "iou": [0.8, 0.6, 0.7, 0.7]})
    pa, _ = sc.match(edges, 2, 2, 0.5)
    assert pa.tolist() == [0, 1]


def test_match_is_one_to_one_on_random_overlaps():
    rng = np.random.default_rng(0)
    a = np.array([shapely.Point(rng.uniform(0, 300, 2)).buffer(6) for _ in range(400)], dtype=object)
    b = np.array([shapely.Point(rng.uniform(0, 300, 2)).buffer(6) for _ in range(400)], dtype=object)
    edges = sc.overlap_graph(a, b)
    for t in (0.1, 0.3):
        pa, pb = sc.match(edges, len(a), len(b), t)
        i = np.flatnonzero(pa >= 0)
        assert (pb[pa[i]] == i).all()
        assert len(set(pa[i])) == len(i)
        e = edges.set_index(["a", "b"])["iou"]
        assert (e.loc[list(zip(i, pa[i]))].to_numpy() >= t).all()


def test_edge_iou_equals_direct_shapely():
    rng = np.random.default_rng(1)
    a = np.array([shapely.Point(rng.uniform(0, 100, 2)).buffer(rng.uniform(3, 8)) for _ in range(60)],
                 dtype=object)
    b = np.array([shapely.Point(rng.uniform(0, 100, 2)).buffer(rng.uniform(3, 8)) for _ in range(60)],
                 dtype=object)
    edges = sc.overlap_graph(a, b)
    assert len(edges)
    for r in edges.itertuples():
        inter = a[r.a].intersection(b[r.b]).area
        assert r.iou == pytest.approx(inter / a[r.a].union(b[r.b]).area)


def test_domain_limits_denominators_not_matching():
    a = grid(2)
    b = grid(2)
    b_dom = np.array([True, True, True, False])
    out = sc.compare(a, b, b_in_domain=b_dom)
    s = out["summary"]["strata"]["all"]
    assert s["n_a"] == 4 and s["n_b"] == 3
    assert s["headline"]["matched_a"] == 4 and s["headline"]["matched_b"] == 3
    assert out["objects_a"]["partner"].tolist() == [0, 1, 2, 3]


def test_strata_partition_the_counts():
    a = grid(3)
    lab = np.array(["pancreas"] * 5 + ["pLN"] * 4, dtype=object)
    out = sc.compare(a, a.copy(), a_strata=lab, b_strata=lab)
    st = out["summary"]["strata"]
    assert st["pancreas"]["n_a"] + st["pLN"]["n_a"] == st["all"]["n_a"] == 9


def test_invalid_geometry_is_repaired_and_counted():
    bowtie = shapely.Polygon([(0, 0), (10, 10), (10, 0), (0, 10)])
    out = sc.compare(np.array([bowtie], dtype=object), np.array([square(0, 0)], dtype=object))
    assert out["summary"]["repaired_geometries"]["a"] == 1
    assert out["objects_a"].loc[0, "area_um2"] == pytest.approx(50.0)


def test_degenerate_polygon_is_emptied_and_stays_an_unmatched_object():
    # Real InstanSeg output (HDL059, 2026-11-12): a nucleus ring with every vertex identical.
    # Repair empties it; its centroid is then an empty point, which GEOS 3.12+ will not read.
    dot = shapely.Polygon([(5, 5)] * 5)
    out = sc.compare(np.array([dot, square(20, 0)], dtype=object),
                     np.array([square(0, 0), square(20, 0)], dtype=object))
    oa = out["objects_a"]
    assert out["summary"]["repaired_geometries"]["a"] == 1
    assert np.isnan(oa.loc[0, "x_um"]) and np.isnan(oa.loc[0, "y_um"])
    assert oa.loc[0, "reason"] == "no_overlap" and oa.loc[0, "area_um2"] == 0
    assert oa.loc[1, "reason"] == "matched"
    x, y = sc.centroid_xy(np.array([dot, None, square(0, 0)], dtype=object))
    assert x[2] == pytest.approx(5.0) and y[2] == pytest.approx(5.0)
    assert np.isnan(x[1]) and not np.isnan(x[0])   # unrepaired, the degenerate ring has a centroid


def test_shared_area_counts_only_where_both_segmented_and_reports_the_rest():
    # a covers x 0-600, b covers x 300-900, both y 0-300, one centroid every 8 um.
    g = np.stack(np.meshgrid(np.arange(0, 900, 8.0), np.arange(0, 300, 8.0)), -1).reshape(-1, 2)
    xy_a, xy_b = g[g[:, 0] < 600], g[g[:, 0] >= 300]
    d = sc.shared_area(xy_a, xy_b)
    xa = xy_a[d["a"], 0]
    # counted: the 300-600 overlap, minus the 20 um margin (one bin of slack for gridding)
    assert xa.min() >= 300 + 20 - 10 and xa.max() <= 600 - 20 + 10
    assert not d["a"][xy_a[:, 0] < 290].any() and not d["b"][xy_b[:, 0] > 610].any()
    assert d["a"][(xy_a[:, 0] > 340) & (xy_a[:, 0] < 560) & (xy_a[:, 1] > 40) & (xy_a[:, 1] < 260)].all()
    # a-only and b-only strips are each ~300 x 300 um = 0.09 mm2
    assert d["area_mm2"]["a_only"] == pytest.approx(0.09, rel=0.1)
    assert d["area_mm2"]["b_only"] == pytest.approx(0.09, rel=0.1)


def test_shared_area_closes_small_gaps_but_not_large_ones():
    g = np.stack(np.meshgrid(np.arange(0, 600, 8.0), np.arange(0, 600, 8.0)), -1).reshape(-1, 2)
    small_hole = np.hypot(g[:, 0] - 150, g[:, 1] - 300) < 20      # 40 um lumen in a only
    big_hole = np.hypot(g[:, 0] - 420, g[:, 1] - 300) < 100       # 200 um patch a missed
    d = sc.shared_area(g[~small_hole & ~big_hole], g)
    inside_small = np.hypot(g[:, 0] - 150, g[:, 1] - 300) < 15
    inside_big = np.hypot(g[:, 0] - 420, g[:, 1] - 300) < 60
    assert d["b"][inside_small].all()          # b objects inside a's small lumen still count
    assert not d["b"][inside_big].any()        # a's large gap is not counted ...
    assert d["area_mm2"]["b_only"] > 0.02      # ... and is reported (~pi*0.1^2 = 0.031 mm2)


def test_hierarchy_reads_nucleus_pairs_against_cell_pairs():
    # nuclei a0<->b0 (cells a0, b0 matched), a1<->b1 (cell a1 matched to b5, not b1's cell b1),
    # a2<->b2 (cell a2 unmatched), a3 unmatched.
    h = sc.hierarchy(nuc_partner_a=np.array([0, 1, 2, -1]),
                     nuc_parent_a=np.array([0, 1, 2, 3]),
                     nuc_parent_b=np.array([0, 1, 2]),
                     cell_partner_a=np.array([0, 5, -1, -1]))
    assert h == {"matched_nuclei": 3, "same_cell_pair": 1, "cells_matched_elsewhere": 1,
                 "a_cell_unmatched": 1, "no_parent_cell": 0, "same_cell_pair_frac": pytest.approx(1 / 3)}


def test_valid_non_polygonal_geometry_is_reduced_to_its_polygons():
    gc = shapely.GeometryCollection([square(0, 0), shapely.LineString([(20, 0), (30, 0)])])
    g, n = sc.polygonal(np.array([gc, square(50, 50), None], dtype=object))
    assert n == 1
    assert shapely.get_type_id(g[0]) == 3 and g[0].area == pytest.approx(100.0)
    assert g[2] is None


def test_close_thresholds_get_distinct_columns():
    out = sc.compare(grid(2), grid(2, dx=1.0), thresholds=(0.2, 0.25, 0.5))
    assert {"matched_0.2", "matched_0.25", "matched_0.5"} <= set(out["objects_a"].columns)


def test_hierarchy_can_be_limited_to_counted_nuclei():
    h = sc.hierarchy(np.array([0, 1]), np.array([0, 1]), np.array([0, 1]), np.array([0, -1]),
                     counted_a=np.array([True, False]))
    assert h["matched_nuclei"] == 1 and h["same_cell_pair"] == 1
