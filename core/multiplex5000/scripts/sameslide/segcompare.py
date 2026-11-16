#!/usr/bin/env python3
"""Two segmentations of one physical section, compared object by object.

On a same-slide section the PhenoCycler image and the Xenium image show the SAME nuclei and the
SAME cells, so their two segmentations -- InstanSeg (``fluorescence_nuclei_and_cells``) in QuPath
on the qptiff, and the Xenium Analyzer's own on its morphology image -- can be compared directly,
object against object, once one is carried into the other's frame. This module is that
comparison. It is geometry only: no marker, transcript, call or label enters it.

**Neither segmentation is ground truth.** Every count is reported from both sides ("InstanSeg
objects with a Xenium partner", "Xenium objects with an InstanSeg partner"), never as precision
and recall of one against the other. Where a single symmetric number is wanted it is
``f1 = (matched_a + matched_b) / (n_a + n_b)``, which is F1 when the two match counts are equal
and stays symmetric when a domain edge makes them differ.

**The frame is the fixed one.** Registration maps Xenium (moving) forward into the PhenoCycler
frame; the inverse of a transform with a displacement field is a different transform (see
``register_points.py``), so side ``b`` is always the warped one and side ``a`` is untouched.

What is measured, per level (nucleus, cell) and per stratum:

* **Overlap graph** -- every positive-area intersection between an ``a`` object and a ``b``
  object, with IoU and the share of each object the intersection covers.
* **One-to-one matching at each IoU threshold** -- the largest number of pairs with IoU at or
  above the threshold, and among those the largest total IoU (the criterion of Schmidt et al.
  2018, StarDist ``matching``). Solved exactly with a sparse assignment, not greedily: greedy
  matching loses pairs whenever one object's best partner is another object's only partner.
* **Panoptic quality** at IoU 0.5 (Kirillov et al. 2019): ``SQ`` the mean IoU of matched pairs,
  ``RQ`` the symmetric F1 above, ``PQ = SQ * RQ``. 0.5 is the conventional criterion and the
  lowest IoU at which two segmentations whose own objects do not overlap can match each object
  to at most one partner.
* **Splits and merges** from containment: an object is *held* by a partner that covers more than
  half its area. More than half is a majority, so at most one partner can hold it; that is the
  definition, not a tuned threshold. A ``b`` object holding two or more ``a`` objects merges
  them; an ``a`` object holding two or more ``b`` objects is split by ``b``, and vice versa.
* **Why an object is unmatched** at IoU 0.5, in precedence order: no overlap with any partner
  (a detection the other side does not have), split (it holds two or more partners), merged
  (its holder holds two or more objects of its own side), boundary (it overlaps, but no
  partner reaches IoU 0.5 and no split or merge explains it).
* **Matched-pair geometry** -- log2 area ratio ``b/a`` and the centroid offset ``b - a``. A mean
  offset away from zero is residual registration, not segmentation.
* **The registration ceiling** -- each ``a`` object's IoU with itself translated by a given
  distance in a seeded random direction. Two identical segmentations misregistered by that
  distance cannot do better, so observed IoU is read against it rather than against 1.

Domain: an object counts in a denominator only if ``in_domain`` says so (the stage restricts to
area both platforms imaged and segmented). Matching itself runs on every object, so an object
just inside the domain keeps the partner whose centroid lies just outside it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import shapely
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching

#: IoU thresholds at which one-to-one matches are counted.
THRESHOLDS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
#: The criterion for the headline numbers (PQ, unmatched reasons, matched-pair geometry).
HEADLINE_IOU = 0.5
#: An object is held by a partner covering MORE than this share of its area: a strict majority.
HOLD_FRACTION = 0.5
#: Quantiles reported for every distribution.
QUANTILES = (0.1, 0.25, 0.5, 0.75, 0.9)

_POLYGON, _MULTIPOLYGON = 3, 6


def polygonal(geoms) -> tuple[np.ndarray, int]:
    """Valid, polygonal geometries: repaired where invalid, reduced to their polygon parts.

    ``make_valid`` can return a GeometryCollection holding a stray line or point, and a valid
    input can be one too; an area-based comparison keeps only the polygonal part of either.
    Returns the array and how many geometries needed repair or reduction.
    """
    g = np.asarray(geoms, dtype=object).copy()
    present = ~shapely.is_missing(g)
    bad = present & ~shapely.is_valid(g)
    if bad.any():
        g[bad] = shapely.make_valid(g[bad])
    kind = shapely.get_type_id(g)
    reduce = present & ~shapely.is_empty(g) & ~np.isin(kind, (_POLYGON, _MULTIPOLYGON))
    for i in np.flatnonzero(reduce):
        parts = [p for p in shapely.get_parts(g[i])
                 if shapely.get_type_id(p) in (_POLYGON, _MULTIPOLYGON)]
        g[i] = shapely.union_all(parts) if parts else shapely.Polygon()
    return g, int((bad | reduce).sum())


def overlap_graph(a, b) -> pd.DataFrame:
    """Every positive-area intersection between an object of ``a`` and an object of ``b``.

    Columns: ``a``, ``b`` (positions), ``inter_um2``, ``iou``, ``frac_of_a``, ``frac_of_b``.
    Objects that only touch share no area and are not edges.
    """
    a = np.asarray(a, dtype=object)
    b = np.asarray(b, dtype=object)
    cols = {"a": np.empty(0, np.int64), "b": np.empty(0, np.int64), "inter_um2": np.empty(0),
            "iou": np.empty(0), "frac_of_a": np.empty(0), "frac_of_b": np.empty(0)}
    if len(a) == 0 or len(b) == 0:
        return pd.DataFrame(cols)
    ai, bi = shapely.STRtree(b).query(a, predicate="intersects")
    inter = shapely.area(shapely.intersection(a[ai], b[bi]))
    keep = np.isfinite(inter) & (inter > 0)
    ai, bi, inter = ai[keep].astype(np.int64), bi[keep].astype(np.int64), inter[keep]
    area_a, area_b = shapely.area(a), shapely.area(b)
    union = area_a[ai] + area_b[bi] - inter
    out = pd.DataFrame({"a": ai, "b": bi, "inter_um2": inter, "iou": inter / union,
                        "frac_of_a": inter / area_a[ai], "frac_of_b": inter / area_b[bi]})
    return out.sort_values(["a", "b"], kind="stable").reset_index(drop=True)


def match(edges: pd.DataFrame, n_a: int, n_b: int, iou_min: float) -> tuple[np.ndarray, np.ndarray]:
    """One-to-one pairs with IoU >= ``iou_min``: the most pairs, then the largest total IoU.

    Returns ``(partner_of_a, partner_of_b)``, positions into the other side, -1 where unmatched.

    An edge whose two ends have no other edge above the threshold is its own component and is
    matched outright. Everything else goes to one sparse assignment in which each ``a`` object
    also has a private "unmatched" column. With cost ``1 - iou / (2m)`` for a real pair and 2 for
    unmatched, the total cost is ``2 n - sum(1 + iou / (2m))``; since ``m`` is at least the
    number of pairs, the IoU term sums to less than one pair, so the minimum has the most pairs
    first and the largest IoU second. Components are independent, so one global solve is exact.
    """
    pa = np.full(n_a, -1, np.int64)
    pb = np.full(n_b, -1, np.int64)
    e = edges[edges["iou"].to_numpy() >= iou_min]
    if e.empty:
        return pa, pb
    ai, bi, w = e["a"].to_numpy(np.int64), e["b"].to_numpy(np.int64), e["iou"].to_numpy(float)
    deg_a = np.bincount(ai, minlength=n_a)
    deg_b = np.bincount(bi, minlength=n_b)
    simple = (deg_a[ai] == 1) & (deg_b[bi] == 1)
    pa[ai[simple]] = bi[simple]
    pb[bi[simple]] = ai[simple]
    if simple.all():
        return pa, pb

    ra, rb, rw = ai[~simple], bi[~simple], w[~simple]
    ua, ia = np.unique(ra, return_inverse=True)
    ub, ib = np.unique(rb, return_inverse=True)
    n_rows, n_cols = len(ua), len(ub)
    m = float(max(min(n_rows, n_cols), 1))
    rows = np.concatenate([ia, np.arange(n_rows)])
    cols = np.concatenate([ib, n_cols + np.arange(n_rows)])
    cost = np.concatenate([1.0 - rw / (2.0 * m), np.full(n_rows, 2.0)])
    graph = csr_matrix((cost, (rows, cols)), shape=(n_rows, n_cols + n_rows))
    r, c = min_weight_full_bipartite_matching(graph)
    real = c < n_cols
    pa[ua[r[real]]] = ub[c[real]]
    pb[ub[c[real]]] = ua[r[real]]
    return pa, pb


def holders(edges: pd.DataFrame, n_a: int, n_b: int) -> dict:
    """Containment both ways: who holds more than half of whom.

    ``held_by_b[i]`` is the ``b`` object covering more than half of ``a`` object ``i`` (else -1);
    ``holds_n_a[j]`` how many ``a`` objects ``b`` object ``j`` holds; and the mirror for ``a``.
    Two partners can both cover more than half of one object only if they overlap each other;
    that is counted in ``conflicts_*`` and the larger share is kept.
    """
    out = {}
    for side, other, n_side, n_other, frac in (("a", "b", n_a, n_b, "frac_of_a"),
                                               ("b", "a", n_b, n_a, "frac_of_b")):
        held = np.full(n_side, -1, np.int64)
        e = edges[edges[frac].to_numpy() > HOLD_FRACTION]
        e = e.sort_values([side, frac], ascending=[True, False], kind="stable")
        conflicts = int(e.duplicated(side).sum())
        e = e.drop_duplicates(side)
        held[e[side].to_numpy(np.int64)] = e[other].to_numpy(np.int64)
        out[f"held_by_{other}"] = held
        out[f"holds_n_{side}"] = np.bincount(held[held >= 0], minlength=n_other)
        out[f"conflicts_{side}"] = conflicts
    return out


def unmatched_reason(partner, n_partners, holds_n_other, held_by_other, holds_n_same):
    """Why each object of one side has no partner at the headline threshold.

    ``holds_n_other``: how many partner-side objects this object holds. ``held_by_other``: the
    partner holding this object. ``holds_n_same``: per partner-side object, how many objects of
    THIS side it holds. Matched objects get ``matched``.
    """
    n = len(partner)
    reason = np.full(n, "boundary", dtype=object)
    holder = held_by_other
    merged = np.zeros(n, bool)
    has = holder >= 0
    merged[has] = holds_n_same[holder[has]] >= 2
    reason[merged] = "merged"
    reason[holds_n_other >= 2] = "split"
    reason[n_partners == 0] = "no_overlap"
    reason[partner >= 0] = "matched"
    return reason


def shift_ceiling(geoms, shift_um: float, seed: int = 0) -> np.ndarray:
    """IoU of each geometry with itself moved ``shift_um`` in a seeded random direction."""
    g = np.asarray(geoms, dtype=object)
    if len(g) == 0 or not np.isfinite(shift_um):
        return np.full(len(g), np.nan)
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0.0, 2.0 * np.pi, len(g))
    offset = shift_um * np.column_stack([np.cos(theta), np.sin(theta)])
    # transform() hands the callback every coordinate in get_coordinates order, so the index
    # from get_coordinates says which geometry -- and which offset -- each vertex belongs to.
    _, idx = shapely.get_coordinates(g, return_index=True)
    moved = shapely.transform(g, lambda c: c + offset[idx])
    inter = shapely.area(shapely.intersection(g, moved))
    union = shapely.area(shapely.union(g, moved))
    return np.divide(inter, union, out=np.full(len(g), np.nan), where=union > 0)


def quantiles(x) -> dict:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return {f"p{int(q * 100)}": None for q in QUANTILES} | {"n": 0}
    return {f"p{int(q * 100)}": float(np.quantile(x, q)) for q in QUANTILES} | {"n": int(len(x))}


def compare(a, b, *, a_in_domain=None, b_in_domain=None, a_strata=None, b_strata=None,
            thresholds=THRESHOLDS, ceiling_shift_um: float = np.nan, seed: int = 0) -> dict:
    """Compare segmentation ``a`` (reference frame) with ``b`` (already in ``a``'s frame).

    ``a_strata`` / ``b_strata``: one label per object (``None`` -> a single stratum ``all``).
    Returns ``{'edges', 'objects_a', 'objects_b', 'thresholds', 'summary', 'a', 'b'}``: the
    overlap graph, one row per object per side, a tidy per-threshold table, a JSON-ready
    summary, and the (repaired) geometries the graph was built on.
    """
    a, repaired_a = polygonal(a)
    b, repaired_b = polygonal(b)
    n_a, n_b = len(a), len(b)
    a_dom = np.ones(n_a, bool) if a_in_domain is None else np.asarray(a_in_domain, bool)
    b_dom = np.ones(n_b, bool) if b_in_domain is None else np.asarray(b_in_domain, bool)
    a_str = np.full(n_a, "all", object) if a_strata is None else np.asarray(a_strata, object)
    b_str = np.full(n_b, "all", object) if b_strata is None else np.asarray(b_strata, object)
    for name, arr, n in (("a_in_domain", a_dom, n_a), ("b_in_domain", b_dom, n_b),
                         ("a_strata", a_str, n_a), ("b_strata", b_str, n_b)):
        if len(arr) != n:
            raise ValueError(f"{name} has {len(arr)} entries for {n} objects")

    edges = overlap_graph(a, b)
    area_a, area_b = shapely.area(a), shapely.area(b)
    ca, cb = shapely.centroid(a), shapely.centroid(b)
    xa, ya = shapely.get_x(ca), shapely.get_y(ca)
    xb, yb = shapely.get_x(cb), shapely.get_y(cb)

    hold = holders(edges, n_a, n_b)
    deg_a = np.bincount(edges["a"].to_numpy(np.int64), minlength=n_a)
    deg_b = np.bincount(edges["b"].to_numpy(np.int64), minlength=n_b)
    best_a = np.zeros(n_a)
    best_b = np.zeros(n_b)
    if len(edges):
        np.maximum.at(best_a, edges["a"].to_numpy(np.int64), edges["iou"].to_numpy(float))
        np.maximum.at(best_b, edges["b"].to_numpy(np.int64), edges["iou"].to_numpy(float))

    matches = {t: match(edges, n_a, n_b, t) for t in thresholds}
    if HEADLINE_IOU not in matches:
        matches[HEADLINE_IOU] = match(edges, n_a, n_b, HEADLINE_IOU)
    pa, pb = matches[HEADLINE_IOU]

    # IoU of each headline pair, read from the graph rather than recomputed.
    iou_lookup = pd.Series(edges["iou"].to_numpy(),
                           index=pd.MultiIndex.from_arrays([edges["a"], edges["b"]]))
    has = pa >= 0
    pair_iou_a = np.full(n_a, np.nan)
    if has.any():
        pair_iou_a[has] = iou_lookup.reindex(
            pd.MultiIndex.from_arrays([np.flatnonzero(has), pa[has]])).to_numpy()
    pair_iou_b = np.full(n_b, np.nan)
    hb = pb >= 0
    pair_iou_b[hb] = pair_iou_a[pb[hb]]

    # holds_n_b is per a object (how many b it holds); holds_n_a is per b object.
    reason_a = unmatched_reason(pa, deg_a, hold["holds_n_b"], hold["held_by_b"], hold["holds_n_a"])
    reason_b = unmatched_reason(pb, deg_b, hold["holds_n_a"], hold["held_by_a"], hold["holds_n_b"])

    ceiling = shift_ceiling(a, ceiling_shift_um, seed)

    objects_a = pd.DataFrame({
        "area_um2": area_a, "x_um": xa, "y_um": ya, "in_domain": a_dom, "stratum": a_str,
        "n_partners": deg_a, "best_iou": best_a, "partner": pa, "pair_iou": pair_iou_a,
        "held_by": hold["held_by_b"], "holds_n": hold["holds_n_b"], "reason": reason_a,
        "ceiling_iou": ceiling})
    objects_b = pd.DataFrame({
        "area_um2": area_b, "x_um": xb, "y_um": yb, "in_domain": b_dom, "stratum": b_str,
        "n_partners": deg_b, "best_iou": best_b, "partner": pb, "pair_iou": pair_iou_b,
        "held_by": hold["held_by_a"], "holds_n": hold["holds_n_a"], "reason": reason_b})
    for t, (ta, tb) in matches.items():
        objects_a[f"matched_{t:g}"] = ta >= 0
        objects_b[f"matched_{t:g}"] = tb >= 0

    rows, strata = [], {}
    labels = ["all"] + sorted({str(s) for s in np.concatenate([a_str, b_str])} - {"all"})
    for s in labels:
        sa = a_dom & (np.ones(n_a, bool) if s == "all" else (a_str.astype(str) == s))
        sb = b_dom & (np.ones(n_b, bool) if s == "all" else (b_str.astype(str) == s))
        strata[s] = _summarise(s, sa, sb, matches, objects_a, objects_b, rows, thresholds)
    summary = {
        "n_a": n_a, "n_b": n_b, "n_edges": int(len(edges)),
        "repaired_geometries": {"a": repaired_a, "b": repaired_b},
        "hold_conflicts": {"a": hold["conflicts_a"], "b": hold["conflicts_b"]},
        "ceiling_shift_um": None if not np.isfinite(ceiling_shift_um) else float(ceiling_shift_um),
        "thresholds": list(thresholds), "headline_iou": HEADLINE_IOU,
        "strata": strata,
    }
    return {"edges": edges, "objects_a": objects_a, "objects_b": objects_b,
            "thresholds": pd.DataFrame(rows), "summary": summary, "a": a, "b": b}


def _share_at_least(x, t):
    """Share of the finite values at or above ``t``; NaN when none is finite."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.mean(x >= t)) if len(x) else np.nan


def _summarise(label, sa, sb, matches, oa, ob, rows, thresholds) -> dict:
    """One stratum: counts both ways at every threshold, then headline geometry."""
    n_a, n_b = int(sa.sum()), int(sb.sum())
    for t in thresholds:
        ta, tb = matches[t]
        m_a, m_b = int((ta[sa] >= 0).sum()), int((tb[sb] >= 0).sum())
        rows.append({"stratum": label, "iou_threshold": t, "n_a": n_a, "n_b": n_b,
                     "matched_a": m_a, "matched_b": m_b,
                     "matched_frac_a": m_a / n_a if n_a else np.nan,
                     "matched_frac_b": m_b / n_b if n_b else np.nan,
                     "f1": (m_a + m_b) / (n_a + n_b) if n_a + n_b else np.nan,
                     "ceiling_frac_a": _share_at_least(oa["ceiling_iou"].to_numpy()[sa], t)})

    pa = oa["partner"].to_numpy()
    pb = ob["partner"].to_numpy()
    # A pair belongs to the stratum if either end is in it; each pair counted once, from a.
    via_b = np.zeros(len(oa), bool)
    via_b[pb[sb & (pb >= 0)]] = True
    idx_a = np.flatnonzero((pa >= 0) & (sa | via_b))
    idx_b = pa[idx_a]
    iou = oa["pair_iou"].to_numpy()[idx_a]
    m_a = int((pa[sa] >= 0).sum())
    m_b = int((pb[sb] >= 0).sum())
    rq = (m_a + m_b) / (n_a + n_b) if n_a + n_b else np.nan
    sq = float(np.mean(iou)) if len(iou) else np.nan
    area_a = oa["area_um2"].to_numpy()[idx_a]
    area_b = ob["area_um2"].to_numpy()[idx_b]
    dx = ob["x_um"].to_numpy()[idx_b] - oa["x_um"].to_numpy()[idx_a]
    dy = ob["y_um"].to_numpy()[idx_b] - oa["y_um"].to_numpy()[idx_a]

    def reasons(o, s):
        return {k: int(v) for k, v in o.loc[s, "reason"].value_counts().items()}

    holds_a = oa["holds_n"].to_numpy()
    holds_b = ob["holds_n"].to_numpy()
    return {
        "n_a": n_a, "n_b": n_b,
        "headline": {"iou_threshold": HEADLINE_IOU, "matched_a": m_a, "matched_b": m_b,
                     "matched_frac_a": m_a / n_a if n_a else None,
                     "matched_frac_b": m_b / n_b if n_b else None,
                     "SQ": sq if np.isfinite(sq) else None,
                     "RQ": rq if np.isfinite(rq) else None,
                     "PQ": sq * rq if np.isfinite(sq) and np.isfinite(rq) else None},
        "reason_a": reasons(oa, sa), "reason_b": reasons(ob, sb),
        "a_split_by_b": {"objects": int((holds_a[sa] >= 2).sum()),
                         "b_objects_involved": int(holds_a[sa][holds_a[sa] >= 2].sum())},
        "b_split_by_a": {"objects": int((holds_b[sb] >= 2).sum()),
                         "a_objects_involved": int(holds_b[sb][holds_b[sb] >= 2].sum())},
        "best_iou_a": quantiles(oa["best_iou"].to_numpy()[sa]),
        "best_iou_b": quantiles(ob["best_iou"].to_numpy()[sb]),
        "pair_iou": quantiles(iou),
        "ceiling_iou_a": quantiles(oa["ceiling_iou"].to_numpy()[sa]),
        "area_um2_a": quantiles(oa["area_um2"].to_numpy()[sa]),
        "area_um2_b": quantiles(ob["area_um2"].to_numpy()[sb]),
        "pair_log2_area_ratio_b_over_a": quantiles(np.log2(area_b / area_a)) if len(idx_a) else quantiles([]),
        "pair_centroid_offset_um": {
            "mean_dx": float(np.mean(dx)) if len(dx) else None,
            "mean_dy": float(np.mean(dy)) if len(dy) else None,
            "distance": quantiles(np.hypot(dx, dy))},
    }


#: Grid pitch of the imaged-area masks: about one cell.
DOMAIN_PITCH_UM = 10.0
#: Gaps in one platform's objects narrower than about twice this are closed (a lumen, a vessel).
DOMAIN_CLOSE_UM = 30.0
#: Objects closer than this to the edge of the shared area are not counted: their partner may lie
#: just past the other platform's edge. Two cell diameters.
DOMAIN_MARGIN_UM = 20.0


def _disk(r: int) -> np.ndarray:
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def shared_area(xy_a, xy_b, *, pitch=DOMAIN_PITCH_UM, close=DOMAIN_CLOSE_UM,
                margin=DOMAIN_MARGIN_UM) -> dict:
    """Where both platforms segmented: per-side occupancy, closed, intersected, eroded.

    Each side's footprint is the grid of bins holding one of its object centroids, closed over
    gaps of about ``2 * close``. The counted area is the two footprints' intersection eroded by
    ``margin``. A patch only one platform segmented is not counted -- and is not hidden: its area
    is reported as ``area_mm2.a_only`` / ``b_only``, so a wholesale failure of one segmentation is a
    number in the summary rather than an absence.

    Returns ``{'a': bool per a point, 'b': bool per b point, 'area_mm2': {...}, 'grid': {...}}``.
    """
    from scipy.ndimage import binary_closing, binary_erosion

    xy_a = np.asarray(xy_a, float).reshape(-1, 2)
    xy_b = np.asarray(xy_b, float).reshape(-1, 2)
    finite_a = np.isfinite(xy_a).all(1)
    finite_b = np.isfinite(xy_b).all(1)
    pts = np.vstack([xy_a[finite_a], xy_b[finite_b]])
    if not len(pts):
        raise ValueError("no finite object centroids")
    rc, rm = int(round(close / pitch)), int(round(margin / pitch))
    pad = rc + rm + 1
    lo = pts.min(0) - pad * pitch
    shape = tuple((np.floor((pts.max(0) + pad * pitch - lo) / pitch).astype(int) + 1)[::-1])

    def bins(xy):
        ij = np.floor((xy - lo) / pitch).astype(int)
        return ij[:, 1], ij[:, 0]

    def footprint(xy):
        occ = np.zeros(shape, bool)
        occ[bins(xy)] = True
        return binary_closing(occ, structure=_disk(rc)) if rc else occ

    fa, fb = footprint(xy_a[finite_a]), footprint(xy_b[finite_b])
    both = fa & fb
    counted = binary_erosion(both, structure=_disk(rm)) if rm else both

    mm2 = pitch * pitch / 1e6
    area = {"mask": counted, "origin": lo, "pitch": pitch,
            "area_mm2": {"a_footprint": float(fa.sum() * mm2), "b_footprint": float(fb.sum() * mm2),
                         "a_only": float((fa & ~fb).sum() * mm2), "b_only": float((fb & ~fa).sum() * mm2),
                         "both": float(both.sum() * mm2), "counted": float(counted.sum() * mm2)},
            "grid": {"pitch_um": pitch, "close_um": close, "margin_um": margin}}
    area["a"], area["b"] = inside(area, xy_a), inside(area, xy_b)
    return area


def inside(area: dict, xy) -> np.ndarray:
    """Whether each point falls in a ``shared_area``'s counted mask (False off the grid or NaN)."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    out = np.zeros(len(xy), bool)
    ok = np.isfinite(xy).all(1)
    ij = np.floor((xy[ok] - area["origin"]) / area["pitch"]).astype(int)
    h, w = area["mask"].shape
    on = (ij[:, 0] >= 0) & (ij[:, 0] < w) & (ij[:, 1] >= 0) & (ij[:, 1] < h)
    hit = np.zeros(len(ij), bool)
    hit[on] = area["mask"][ij[on, 1], ij[on, 0]]
    out[np.flatnonzero(ok)] = hit
    return out


def one_side(objects: pd.DataFrame, mask, thresholds=THRESHOLDS) -> dict:
    """One side's objects under a mask only that side has a label for (e.g. Xenium's
    ``segmentation_method``): matched share at each threshold, unmatched reasons, IoU."""
    o = objects.loc[np.asarray(mask, bool)]
    n = int(len(o))
    return {"n": n,
            "matched_frac": {f"{t:g}": float(o[f"matched_{t:g}"].mean()) if n else None
                             for t in thresholds},
            "reason": {k: int(v) for k, v in o["reason"].value_counts().items()},
            "best_iou": quantiles(o["best_iou"]), "pair_iou": quantiles(o["pair_iou"]),
            "area_um2": quantiles(o["area_um2"])}


def hierarchy(nuc_partner_a, nuc_parent_a, nuc_parent_b, cell_partner_a, counted_a=None) -> dict:
    """Do nuclei matched across platforms sit in cells matched to each other?

    ``nuc_partner_a``: per ``a`` nucleus, its matched ``b`` nucleus (-1 unmatched).
    ``nuc_parent_a`` / ``nuc_parent_b``: per nucleus, the position of its cell on its own side
    (-1 if it has none). ``cell_partner_a``: per ``a`` cell, its matched ``b`` cell (-1).
    ``counted_a``: restrict to these ``a`` nuclei (the counted area), like every other number.
    """
    keep = np.asarray(nuc_partner_a) >= 0
    if counted_a is not None:
        keep &= np.asarray(counted_a, bool)
    i = np.flatnonzero(keep)
    j = np.asarray(nuc_partner_a)[i]
    ca = np.asarray(nuc_parent_a)[i]
    cb = np.asarray(nuc_parent_b)[j]
    ok = (ca >= 0) & (cb >= 0)
    partner = np.full(len(i), -2, np.int64)
    partner[ok] = np.asarray(cell_partner_a)[ca[ok]]
    same = ok & (partner == cb)
    other = ok & (partner >= 0) & (partner != cb)
    a_cell_unmatched = ok & (partner == -1)
    n = len(i)
    return {"matched_nuclei": int(n),
            "same_cell_pair": int(same.sum()),
            "cells_matched_elsewhere": int(other.sum()),
            "a_cell_unmatched": int(a_cell_unmatched.sum()),
            "no_parent_cell": int((~ok).sum()),
            "same_cell_pair_frac": float(same.sum() / n) if n else None}
