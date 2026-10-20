"""H1-H4: protein spatial endpoints per donor, under the frozen spatial_v2 protocol as amended.

Ported from the archived evaluate_spatial_H1_H2.py and evaluate_spatial_H3_H4.py (2026-10-06) with
the maintainer's 2026-10-17 rulings applied: no support cutoff anywhere, only computability; the
refined hierarchical identity is the identity, not a sensitivity. An endpoint that cannot be computed
is UNAVAILABLE with its reason, never a measured zero.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from shapely.strtree import STRtree

from . import geometry
from .common import Context, binding, parallel_map, peak_rss_gib, read_json, require, write_json
from .sources import audit_binding

ENDOCRINE = ("Beta", "Alpha", "Delta")
# H2 graph and H1 kernels exactly as the archived spatial_H1_H2_plan.json froze them.
H2_GRAPH_K = 6
H1_BINS = (64, 4096)
H3_CONFIGS = (("primary", 50.0, 4), ("graph_k2", 50.0, 2), ("graph_k8", 50.0, 8),
              ("radius25", 25.0, 4), ("radius75", 75.0, 4))
H4_CONFIGS = (("primary", 50.0), ("radius25", 25.0), ("radius75", 75.0))
NUISANCE_MARKERS = ("HLA_A", "PD_L1", "PD_1", "Granzyme_B", "CD107a", "CD31", "B3TUBB",
                    "Collagen_IV", "E_cadherin", "EpCAM")
CELL_COLUMNS = ["object_id", "X_centroid", "Y_centroid", "nucleus_area_um2", "fundamentals_cell_area_um2",
                "endocrine_object_id", "endocrine_object_class", "immune_class",
                "hormone_measurement_state", "annotation_name", "broad_endocrine", "broad_immune",
                "refinement_keeps_broad_label", "acq_run", "stain_group", "disease_status"]
# Per-task memory budgets (GiB): the largest peaks in each stage's resource_usage.csv on 2026-10-18
# (H1 4.7, H2 2.7, H3 2.9, H4 3.2) plus about a quarter.
GIB_PER_TASK = {"h1": 6, "h2": 4, "h3": 4, "h4": 4}


# ---------------------------------------------------------------- shared inputs

def frozen(context: Context) -> dict:
    return read_json(context.task("protocol") / "frozen_protocol.json")


def donors(context: Context) -> list[dict]:
    qualified = context.task("qualify")
    table = pd.read_csv(qualified / "protein_scope_qualification.csv", dtype={"donor_id": str})
    records = read_json(context.task("reconcile") / "selected_inputs.json")["protein_records"]
    native = {row["donor_id"]: row["native_partition_id"] for row in records}
    return [{"donor": row.donor_id, "native": native[row.donor_id], "disease": row.disease_status,
             "rows": int(row.biological_rows),
             "cells": str(qualified / "protein" / f"donor_id={row.donor_id}" / "cells.parquet")}
            for row in table.itertuples()]


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation in GIB_PER_TASK, f"Unknown protein endpoint {operation}")
    natives = [row["native"] for row in donors(context)]
    return [binding(Path(geometry.__file__)), binding(context.task("qualify") / "completion.json"),
            audit_binding(context, geometry.QUALIFIED_OBJECTS), *geometry.bindings(context, natives)]


def seed(*parts) -> int:
    return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:8], 16)


def identity(cells: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Hierarchical refined identity: broad label kept by refinement, then its subclass."""
    kept = cells.refinement_keeps_broad_label.fillna(False).to_numpy(bool)
    name = cells.annotation_name.astype("string").fillna("").to_numpy(str)
    label = np.full(len(cells), "", dtype=object)
    endocrine = cells.broad_endocrine.fillna(False).to_numpy(bool) & kept
    label[endocrine] = "Endocrine (subtype unresolved)"
    subtype = endocrine & np.isin(name, ENDOCRINE)
    label[subtype] = name[subtype]
    immune = cells.broad_immune.fillna(False).to_numpy(bool) & kept
    label[immune] = "Immune"
    t_cell = immune & cells.immune_class.astype("string").eq("Imm-T-cell").fillna(False).to_numpy(bool)
    return label.astype(str), t_cell


def load(task: dict, markers=()) -> dict:
    cells = pd.read_parquet(task["cells"], columns=CELL_COLUMNS + [f"redsea_aligned_{m}" for m in markers])
    require(cells.object_id.is_unique, f"Duplicate cell UUIDs: {task['donor']}")
    xy = cells[["X_centroid", "Y_centroid"]].to_numpy(float)
    require(np.isfinite(xy).all(), f"Nonfinite coordinates in qualified cells: {task['donor']}")
    nucleus = cells.nucleus_area_um2.to_numpy(float)
    require((np.isfinite(nucleus) & (nucleus > 0)).all(), f"A cell without a nucleus: {task['donor']}")
    area = cells.fundamentals_cell_area_um2.to_numpy(float)
    good = np.isfinite(area) & (area > 0)
    require(good.any(), f"No positive cell area: {task['donor']}")
    label, t_cell = identity(cells)
    islets = geometry.qualified_islets(Path(task["objects"]), task["donor"])
    outlines, audit = geometry.true_islets(Path(task["files"]["islets"]), task["pixel"],
                                           set(islets.endocrine_object_id))
    parent = cells.endocrine_object_id.astype("string").fillna("").to_numpy(str)
    core = cells.endocrine_object_class.eq("Islet").fillna(False).to_numpy(bool) & np.isin(parent, list(outlines))
    members = {}
    for uid, index in pd.Series(np.flatnonzero(core)).groupby(parent[core], sort=False):
        index = index.to_numpy(int)
        inside = shapely.covers(outlines[uid], shapely.points(xy[index]))
        require(inside.all(), f"Islet member outside its own outline: {task['donor']} {uid}")
        members[uid] = index
    audit.update(cells=len(cells), core_members=int(core.sum()),
                 mean_cell_diameter_um=float(np.mean(2 * np.sqrt(area[good] / np.pi))))
    return {"cells": cells, "xy": xy, "label": label, "t_cell": t_cell, "outlines": outlines,
            "object_id": cells.object_id.to_numpy().astype(str),
            "members": members, "audit": audit, "ell": audit["mean_cell_diameter_um"],
            "acq_run": "|".join(sorted(cells.acq_run.dropna().astype(str).unique())),
            "stain_group": "|".join(sorted(cells.stain_group.dropna().astype(str).unique()))}


def base_row(task: dict, data: dict, uid: str, geom) -> dict:
    core = data["members"].get(uid, np.array([], int))
    lab = data["label"][core]
    return {"donor_id": task["donor"], "disease_status": task["disease"], "islet_uuid": uid,
            "islet_area_um2": float(geom.area), "geometry_type": geom.geom_type,
            "n_core_cells": len(core), "n_core_endocrine": int(np.isin(lab, ENDOCRINE).sum()),
            "n_core_immune": int((lab == "Immune").sum())}


def median_area_example(rows: list[dict]) -> str | None:
    """The available islet nearest the donor's median available area; chosen by size, never by effect."""
    available = [row for row in rows if row["status"] == "available"]
    if not available:
        return None
    median = float(np.median([row["islet_area_um2"] for row in available]))
    return min(available, key=lambda row: (abs(row["islet_area_um2"] - median), row["islet_uuid"]))["islet_uuid"]


def unavailable(row: dict, reason: str) -> dict:
    return dict(row, status="unavailable", reason=reason)


# ---------------------------------------------------------------- H1

def ring_segments(geom, exterior_only=False):
    starts, ends = [], []
    for part in (list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]):
        for ring in ([part.exterior] if exterior_only else [part.exterior, *part.interiors]):
            coords = np.asarray(ring.coords, float)
            starts.append(coords[:-1])
            ends.append(coords[1:])
    start, end = np.concatenate(starts), np.concatenate(ends)
    valid = np.linalg.norm(end - start, axis=1) > 0
    return start[valid], end[valid]


def exterior_projection(geom, points, max_distance):
    """Exact distance to, and arclength along, the outer ring; equal nearest segments resolve to the first."""
    start, end = ring_segments(geom, True)
    vector = end - start
    length = np.linalg.norm(vector, axis=1)
    cumulative = np.r_[0, np.cumsum(length)[:-1]]
    tree = STRtree(shapely.linestrings(np.stack([start, end], axis=1)))
    pairs, distance = tree.query_nearest(points, all_matches=True, return_distance=True, max_distance=max_distance)
    order = np.lexsort((pairs[1], pairs[0]))
    point, segment = pairs[:, order]
    distance = distance[order]
    _, first = np.unique(point, return_index=True)
    point, segment, distance = point[first], segment[first], distance[first]
    xy = np.column_stack([shapely.get_x(points[point]), shapely.get_y(points[point])])
    fraction = np.clip(np.sum((xy - start[segment]) * vector[segment], axis=1) / length[segment] ** 2, 0, 1)
    return point, np.asarray(distance, float), cumulative[segment] + fraction * length[segment]


def boundary_groups(data: dict) -> tuple[dict, dict]:
    """Every cell within the kernel's numeric reach of an outline, grouped by its own nearest islet.

    An islet member always belongs to its own islet; a cell inside an outline belongs to it; any other
    cell belongs to the islet whose boundary is nearest. No band or expansion is drawn.
    """
    ell, outlines, xy = data["ell"], data["outlines"], data["xy"]
    ids = list(outlines)
    polygons = np.asarray([outlines[uid] for uid in ids], dtype=object)
    points = shapely.points(xy)
    radius = -2 * ell * np.log(1e-8)
    segments, owners = [], []
    for j, geom in enumerate(polygons):
        start, end = ring_segments(geom)
        segments.append(np.stack([start, end], axis=1))
        owners.append(np.full(len(start), j, int))
    owner = np.concatenate(owners)
    tree = STRtree(shapely.linestrings(np.concatenate(segments)))
    pairs, distance = tree.query_nearest(points, all_matches=False, return_distance=True, max_distance=radius)
    near = np.full(len(xy), -1, int)
    near[pairs[0]] = owner[pairs[1]]
    inside = STRtree(polygons).query(points, predicate="covered_by")
    overlapping = 0
    if inside.shape[1]:
        order = np.lexsort((inside[1], inside[0]))
        point, polygon = inside[:, order]
        unique, first, counts = np.unique(point, return_index=True, return_counts=True)
        near[unique] = polygon[first]
        overlapping = int((counts > 1).sum())
    position = {uid: j for j, uid in enumerate(ids)}
    for uid, index in data["members"].items():
        near[index] = position[uid]
    groups = {}
    for j, index in pd.Series(np.flatnonzero(near >= 0)).groupby(near[near >= 0], sort=False):
        geom = polygons[int(j)]
        if geom.geom_type != "Polygon":
            continue
        index = index.to_numpy(int)
        selected, radial, arc = exterior_projection(geom, points[index], radius)
        keep = np.exp(-radial / (2 * ell)) >= 1e-8
        groups[ids[int(j)]] = (index[selected[keep]], radial[keep], arc[keep])
    return groups, {"kernel_reach_um": float(radius), "cells_inside_overlapping_outlines": overlapping,
                    "kernel": "exp(-distance to the outer boundary / donor mean cell diameter)"}


def circular(residual, profile):
    """Every nonzero circular shift of the profile: the null tests orientation, not abundance."""
    corr = np.fft.ifft(np.conj(np.fft.fft(residual)) * np.fft.fft(profile)).real
    observed = float(np.dot(residual, profile))
    require(np.isclose(corr[0], observed, rtol=1e-7, atol=1e-9), "Circular score identity failed")
    return observed, corr[1:]


def nearest_distance(tree, points):
    index, distance = tree.query_nearest(points, all_matches=False, return_distance=True)
    out = np.full(len(points), np.nan)
    out[index[0]] = distance
    return out


def h1_islet(row, geom, group, label, ell, vascular, neural):
    if geom.geom_type != "Polygon":
        return unavailable(row, "multipart outline has no single circular arclength"), None
    if vascular is None or neural is None:
        return unavailable(row, "no vascular or no neural objects in this section"), None
    if group is None:
        return unavailable(row, "no admitted cell in the boundary kernel"), None
    index, radial, arc = group
    weight = np.exp(-radial / ell)
    immune = label[index] == "Immune"
    mass = float(weight[immune].sum())
    row.update(n_boundary_cells=len(index), n_boundary_immune=int(immune.sum()), immune_kernel_mass=mass,
               allcell_kernel_mass=float(weight.sum()))
    if not immune.any() or mass <= 0:
        return unavailable(row, "no refined immune cell in the boundary kernel"), None
    perimeter = float(geom.exterior.length)
    n = int(min(max(H1_BINS[0], np.ceil(perimeter / ell)), H1_BINS[1]))
    bins = np.minimum((arc / perimeter * n).astype(int), n - 1)
    sigma = max(ell * n / perimeter, 0.5)
    background = gaussian_filter1d(np.bincount(bins, weights=weight, minlength=n), sigma, mode="wrap")
    exposure = gaussian_filter1d(np.bincount(bins[immune], weights=weight[immune], minlength=n), sigma, mode="wrap")
    residual = exposure / exposure.sum() - background / background.sum()
    arclength = (np.arange(n) + 0.5) * perimeter / n
    samples = shapely.line_interpolate_point(geom.exterior, arclength)
    v = np.exp(-nearest_distance(vascular, samples) / ell)
    m = np.exp(-nearest_distance(neural, samples) / ell)
    row.update(n_boundary_bins=n, boundary_perimeter_um=perimeter, vascular_profile_mean=float(v.mean()),
               neural_profile_mean=float(m.mean()))
    if not (np.isfinite(v).all() and np.isfinite(m).all()):
        return unavailable(row, "nonfinite structural distance"), None
    if np.std(v - m) < 1e-8:
        return unavailable(row, "vascular-minus-neural profile does not vary around the boundary"), None
    observed, null = circular(residual, v - m)
    sd = float(null.std(ddof=1))
    if sd < 1e-10:
        return unavailable(row, "circular-shift null has no variance"), None
    vascular_observed, vascular_null = circular(residual, v)
    neural_observed, neural_null = circular(residual, m)
    row.update(status="available", reason="", endpoint=observed - float(null.mean()), observed_score=observed,
               null_mean=float(null.mean()), null_sd=sd, n_null=len(null),
               conditional_rank_p=float((1 + (null >= observed - 1e-12).sum()) / (len(null) + 1)),
               vascular_excess=vascular_observed - float(vascular_null.mean()),
               neural_excess=neural_observed - float(neural_null.mean()),
               null_standardized_excess=(observed - float(null.mean())) / sd)
    example = {"arc_fraction": arclength / perimeter, "residual": residual, "vascular_profile": v,
               "neural_profile": m, "immune_profile": exposure / exposure.sum(),
               "background_profile": background / background.sum(), "null": null, "observed": observed}
    return row, example


def h1_donor(task, data):
    vascular_objects, vascular_audit = geometry.objects(Path(task["files"]["Blood vessel"]), task["pixel"])
    neural_objects, neural_audit = geometry.objects(Path(task["files"]["Nerve"]), task["pixel"])
    vascular = STRtree(vascular_objects) if vascular_objects else None
    neural = STRtree(neural_objects) if neural_objects else None
    groups, kernel = boundary_groups(data)
    rows, examples = [], {}
    for uid, geom in data["outlines"].items():
        row = dict(base_row(task, data, uid, geom), hypothesis="H1", configuration="primary",
                   kernel_scale_um=data["ell"], endpoint=np.nan, status="", reason="")
        row, example = h1_islet(row, geom, groups.get(uid), data["label"], data["ell"], vascular, neural)
        rows.append(row)
        if example is not None:
            examples[uid] = example
    chosen = median_area_example(rows)
    example = None
    if chosen:
        geom = data["outlines"][chosen]
        window = geom.buffer(4 * data["ell"])
        nearby = np.flatnonzero(shapely.covers(window, shapely.points(data["xy"])))
        example = {"islet_uuid": chosen, "outline_wkt": geom.wkt, **examples[chosen],
                   "cell_xy": data["xy"][nearby], "cell_label": data["label"][nearby],
                   "vascular_wkt": [vascular_objects[i].intersection(window).wkt
                                    for i in vascular.query(window, predicate="intersects")],
                   "neural_wkt": [neural_objects[i].intersection(window).wkt
                                  for i in neural.query(window, predicate="intersects")]}
    audit = {"kernel": kernel, "vascular_objects": vascular_audit, "neural_objects": neural_audit}
    return rows, [summarize(task, data, rows, "primary", ("vascular_excess", "neural_excess",
                                                          "null_standardized_excess"))], example, audit


# ---------------------------------------------------------------- H2

def graph_edges(xy, geom, k, ell):
    n = len(xy)
    if n < 2:
        return np.array([], int), np.array([], int), np.array([], float), 0
    _, neighbour = cKDTree(xy).query(xy, k=min(k + 1, n))
    neighbour = neighbour.reshape(n, -1)
    pairs = np.column_stack([np.repeat(np.arange(n), neighbour.shape[1] - 1), neighbour[:, 1:].reshape(-1)])
    pairs = np.unique(np.sort(pairs, axis=1), axis=0)
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    within = shapely.covers(geom, shapely.linestrings(xy[pairs]))
    pairs = pairs[within]
    weight = np.exp(-np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1) / ell)
    return pairs[:, 0], pairs[:, 1], weight, int((~within).sum())


def fragmentation(labels, edges):
    i, j, w, _ = edges
    beta = labels == "Beta"
    denominator = float((w * (beta[i].astype(float) + beta[j].astype(float))).sum())
    return np.nan if denominator <= 0 else 1 - float(2 * w[beta[i] & beta[j]].sum()) / denominator


def interface(labels, edges):
    i, j, w, _ = edges
    beta, immune = labels == "Beta", labels == "Immune"
    denominator = float((w * (beta[i].astype(float) + beta[j].astype(float))).sum())
    return np.nan if denominator <= 0 else float(w[(beta[i] & immune[j]) | (beta[j] & immune[i])].sum()) / denominator


def h2_islet(task, data, uid, geom, null_n):
    row = dict(base_row(task, data, uid, geom), hypothesis="H2", configuration="primary",
               graph_k=H2_GRAPH_K, endpoint=np.nan, status="", reason="",
               fragmentation_status="unavailable", interface_status="unavailable")
    index = data["members"].get(uid, np.array([], int))
    if len(index) < 2:
        return unavailable(row, "fewer than two admitted cells in the outline: no graph"), None
    labels = data["label"][index]
    xy = data["xy"][index]
    endocrine = np.isin(labels, ENDOCRINE)
    sub = labels[endocrine]
    n_beta, n_other = int((sub == "Beta").sum()), int((sub != "Beta").sum())
    row.update(n_beta=n_beta, n_alpha=int((sub == "Alpha").sum()), n_delta=int((sub == "Delta").sum()),
               n_immune=int((labels == "Immune").sum()), n_endocrine=len(sub),
               beta_fraction=n_beta / len(sub) if len(sub) else np.nan)
    if not n_beta:
        return unavailable(row, "no refined Beta cell"), None
    if not n_other:
        return unavailable(row, "no other refined endocrine subtype: the composition-preserving shuffle is the identity"), None
    endocrine_edges = graph_edges(xy[endocrine], geom, H2_GRAPH_K, data["ell"])
    all_edges = graph_edges(xy, geom, H2_GRAPH_K, data["ell"])
    row.update(n_endocrine_edges=len(endocrine_edges[0]), n_core_edges=len(all_edges[0]))
    depth = np.asarray(shapely.distance(shapely.points(xy[endocrine]), geom.boundary), float)
    strata_n = min(5, max(1, len(sub) // 10))
    order = np.lexsort((data["object_id"][index[endocrine]], depth))
    strata = np.empty(len(sub), int)
    for z, members in enumerate(np.array_split(order, strata_n)):
        strata[members] = z
    blocks = [np.flatnonzero(strata == z) for z in range(strata_n)]
    movable = sum(len(b) for b in blocks if len(np.unique(sub[b])) > 1)
    row.update(n_radial_rank_strata=strata_n, n_exchangeable_endocrine_positions=movable)
    if movable < 2:
        return unavailable(row, "no radial stratum mixes subtypes: no exchangeable positions"), None
    observed = np.array([fragmentation(sub, endocrine_edges), interface(labels, all_edges)])
    rng = np.random.default_rng(seed(task["seed"], task["donor"], uid, "H2"))
    null = np.empty((null_n, 2))
    for b in range(null_n):
        shuffled = sub.copy()
        for block in blocks:
            shuffled[block] = rng.permutation(shuffled[block])
        relabelled = labels.copy()
        relabelled[endocrine] = shuffled
        null[b] = [fragmentation(shuffled, endocrine_edges), interface(relabelled, all_edges)]
    mean, sd = null.mean(axis=0), null.std(axis=0, ddof=1)
    z = (observed - mean) / np.where(sd > 1e-10, sd, np.nan)
    usable = np.isfinite(observed) & np.isfinite(sd) & (sd > 1e-10)
    reasons = ("no Beta graph edge, or the shuffle never changes fragmentation",
               "no immune cell or no Beta edge: the immune interface cannot vary")
    for k, name in enumerate(("fragmentation", "interface")):
        row.update({f"observed_{name}": float(observed[k]), f"null_{name}_mean": float(mean[k]),
                    f"null_{name}_sd": float(sd[k]), f"{name}_status": "available" if usable[k] else "unavailable",
                    f"{name}_reason": "" if usable[k] else reasons[k]})
        if usable[k]:
            row.update({f"{name}_z": float(z[k]), f"{name}_excess": float(observed[k] - mean[k])})
    if not usable.all():
        return unavailable(row, "joint needs both components: " + "; ".join(
            reasons[k] for k in range(2) if not usable[k])), None
    joint = ((null - mean) / sd).mean(axis=1)
    endpoint = float(z.mean())
    row.update(status="available", reason="", endpoint=endpoint, n_null=null_n,
               null_q025=float(np.quantile(joint, 0.025)), null_q975=float(np.quantile(joint, 0.975)),
               conditional_rank_p=float((1 + (joint >= endpoint - 1e-12).sum()) / (null_n + 1)),
               both_components_positive=bool((z > 0).all()))
    example = {"xy": xy, "labels": labels, "core_edges": np.column_stack(all_edges[:2]),
               "null_joint": joint, "observed": endpoint}
    return row, example


def h2_donor(task, data):
    rows, examples = [], {}
    for uid, geom in data["outlines"].items():
        row, example = h2_islet(task, data, uid, geom, task["null"])
        rows.append(row)
        if example is not None:
            examples[uid] = example
    summary = summarize(task, data, rows, "primary", ("fragmentation_z", "interface_z", "fragmentation_excess",
                                                       "interface_excess", "beta_fraction"))
    for name in ("fragmentation", "interface"):
        component = [row for row in rows if row.get(f"{name}_status") == "available"]
        summary[f"{name}_z_any_available_islet"] = float(np.mean([row[f"{name}_z"] for row in component])) if component else np.nan
        summary[f"n_islets_{name}_available"] = len(component)
    chosen = median_area_example(rows)
    example = dict(islet_uuid=chosen, outline_wkt=data["outlines"][chosen].wkt, **examples[chosen]) if chosen else None
    return rows, [summary], example, {}


# ---------------------------------------------------------------- H3 / H4

def rankz(values):
    values = np.asarray(values, float)
    out = np.full(len(values), np.nan)
    ok = np.isfinite(values)
    if ok.sum() > 1:
        q = (rankdata(values[ok], method="average") - 0.5) / ok.sum()
        sd = q.std()
        out[ok] = (q - q.mean()) / sd if sd > 1e-10 else 0
    return out


def center(values, groups):
    values = np.array(values, float)
    one = values.ndim == 1
    values = values[:, None] if one else values
    counts = np.bincount(groups)
    for j in range(values.shape[1]):
        values[:, j] = values[:, j] - (np.bincount(groups, weights=values[:, j]) / counts)[groups]
    return values[:, 0] if one else values


def design(values, groups=None):
    x = np.asarray(values, float)
    if groups is not None:
        x = center(x, groups)
    sd = x.std(axis=0)
    keep = sd > 1e-8
    return x[:, keep] / sd[keep], keep


def solver(x, w):
    root = np.sqrt(w)
    inverse = np.linalg.pinv(x * root[:, None], rcond=1e-9)
    return lambda y: inverse @ (y * root)


def neighbour_features(xy, t_mask, immune_mask, measured, endocrine, radius, k=12):
    """Local exposure of each endocrine cell: T PD-1 and cytotoxic state, immune and all-cell density,
    and marker-defined structural context. A cell with no T neighbour has zero T exposure."""
    e = xy[endocrine]
    tree = cKDTree(xy)
    distance, index = tree.query(e, k=min(k + 1, len(xy)), workers=2)
    index, distance = index[:, 1:], distance[:, 1:]
    local = distance <= radius
    density = tree.query_ball_point(e, radius, return_length=True, workers=2) - 1
    immune = (cKDTree(xy[immune_mask]).query_ball_point(e, radius, return_length=True, workers=2)
              if immune_mask.any() else np.zeros(len(e)))
    t_index = np.flatnonzero(t_mask)
    neighbours = cKDTree(xy[t_mask]).query_ball_point(e, radius, workers=2)
    pd1 = rankz(measured["PD_1"][t_mask])
    cyto = (rankz(measured["Granzyme_B"][t_mask]) + rankz(measured["CD107a"][t_mask])) / 2
    count = np.asarray([len(v) for v in neighbours])
    t_pd1, t_cyto = np.full(len(e), np.nan), np.full(len(e), np.nan)
    for i, ids in enumerate(neighbours):
        ids = np.asarray(ids, int)
        if not len(ids):
            t_pd1[i] = t_cyto[i] = 0
            continue
        w = np.exp(-0.5 * (np.linalg.norm(xy[t_index[ids]] - e[i], axis=1) / (radius / 2)) ** 2)
        good, goodc = np.isfinite(pd1[ids]), np.isfinite(cyto[ids])
        if good.any():
            t_pd1[i] = np.average(pd1[ids][good], weights=w[good])
        if goodc.any():
            t_cyto[i] = np.average(cyto[ids][goodc], weights=w[goodc])
    structure = []
    for marker in ("CD31", "B3TUBB", "Collagen_IV"):
        v = rankz(measured[marker])
        valid = local & np.isfinite(v[index])
        n = valid.sum(axis=1)
        structure.append(np.divide(np.where(valid, v[index], 0).sum(axis=1), n,
                                   out=np.full(len(e), np.nan), where=n > 0))
    return {"tcount": count, "pd1": t_pd1, "cyto": t_cyto, "immune": immune, "density": density,
            "structure": np.asarray(structure).T}


def endocrine_frame(data, measured, radius):
    """One row per refined Beta/Alpha/Delta cell inside its own true Islet with a read hormone."""
    cells, xy, label = data["cells"], data["xy"], data["label"]
    member = np.zeros(len(cells), bool)
    for index in data["members"].values():
        member[index] = True
    usable = cells.hormone_measurement_state.eq("usable_for_typing").fillna(False).to_numpy(bool)
    endocrine = np.isin(label, ENDOCRINE) & member & usable
    ei = np.flatnonzero(endocrine)
    parent = cells.endocrine_object_id.astype(str).to_numpy()[ei]
    subtype = label[ei]
    hla, pdl1 = np.full(len(ei), np.nan), np.full(len(ei), np.nan)
    for name in ENDOCRINE:
        target = np.flatnonzero(subtype == name)
        hla[target] = rankz(measured["HLA_A"][ei[target]])
        pdl1[target] = rankz(measured["PD_L1"][ei[target]])
    depth, size = np.empty(len(ei)), np.empty(len(ei))
    for uid in np.unique(parent):
        chosen = np.flatnonzero(parent == uid)
        geom = data["outlines"][uid]
        depth[chosen] = shapely.distance(shapely.points(xy[ei[chosen]]), geom.boundary)
        size[chosen] = geom.area
    features = neighbour_features(xy, data["t_cell"], label == "Immune", measured, endocrine, radius)
    nucleus = cells.nucleus_area_um2.to_numpy(float)[ei]
    area = cells.fundamentals_cell_area_um2.to_numpy(float)[ei]
    frame = pd.DataFrame({
        "cell_id": cells.object_id.to_numpy()[ei], "islet": parent, "subtype": subtype,
        "x": xy[ei, 0], "y": xy[ei, 1], "hla": hla, "pdl1": pdl1, "pos": depth / np.sqrt(size / np.pi),
        "nuc": np.log1p(nucleus), "area": np.log1p(area), "tcount": np.log1p(features["tcount"]),
        "immune": np.log1p(features["immune"]), "density": np.log1p(features["density"]),
        "pd1": features["pd1"], "cyto": features["cyto"], "vascular": features["structure"][:, 0],
        "neural": features["structure"][:, 1], "collagen": features["structure"][:, 2],
        "ecad": rankz(measured["E_cadherin"][ei]), "epcam": rankz(measured["EpCAM"][ei]),
        "T_neighbour_count": features["tcount"]})
    frame["pos2"] = frame.pos ** 2
    return frame


H3_VARIABLES = ["hla", "pos", "pos2", "nuc", "area", "tcount", "immune", "density", "pd1", "cyto",
                "vascular", "neural", "collagen", "ecad", "epcam"]


def h3_graph(frame, outlines, degree, cap):
    groups, edge_i, edge_j, edge_group = [], [], [], []
    cell_group = np.full(len(frame), -1, int)
    for (islet, subtype), members in frame.groupby(["islet", "subtype"], sort=True):
        if len(members) < 2:
            continue
        ids = members.index.to_numpy(int)
        xy = members[["x", "y"]].to_numpy(float)
        distance, neighbour = cKDTree(xy).query(xy, k=min(degree + 1, len(xy)))
        pairs = np.column_stack([np.repeat(np.arange(len(xy)), neighbour.shape[1] - 1), neighbour[:, 1:].ravel()])
        pairs = np.unique(np.sort(pairs[distance[:, 1:].ravel() <= cap], axis=1), axis=0)
        if not len(pairs):
            continue
        pairs = pairs[shapely.covers(outlines[islet], shapely.linestrings(xy[pairs]))]
        if not len(pairs):
            continue
        g = len(groups)
        groups.append({"islet": islet, "subtype": subtype, "n": len(ids), "edges": len(pairs)})
        cell_group[ids] = g
        edge_i.extend(ids[pairs[:, 0]])
        edge_j.extend(ids[pairs[:, 1]])
        edge_group.extend([g] * len(pairs))
    return groups, np.asarray(edge_i, int), np.asarray(edge_j, int), np.asarray(edge_group, int), cell_group


def h3_config(frame, outlines, degree, cap, null_n, rng):
    out = {"donor_effect": np.nan, "status": "unavailable", "n_cells": 0, "n_islets": 0, "n_edges": 0}
    groups, ei, ej, eg, cg = h3_graph(frame, outlines, degree, cap)
    if not groups:
        return dict(out, reason="no homotypic graph edge inside any outline"), [], None
    retained = np.flatnonzero(cg >= 0)
    sub = frame.iloc[retained].reset_index(drop=True)
    remap = np.full(len(frame), -1, int)
    remap[retained] = np.arange(len(retained))
    ei, ej, cg = remap[ei], remap[ej], cg[retained]
    y = center(sub.hla.to_numpy(), cg)
    x, _ = design(sub[H3_VARIABLES[1:]].to_numpy(), cg)
    islet_index = pd.factorize(sub.islet, sort=True)[0]
    w = 1 / np.bincount(islet_index)[islet_index]
    fit = solver(x, w)
    fitted = x @ fit(y)
    residual = y - fitted
    out.update(n_cells=len(sub), n_edges=len(ei), nuisance_rank=int(np.linalg.matrix_rank(x)),
               residual_variance=float(residual.var()))
    if residual.var() < 1e-12:
        return dict(out, reason="no residual HLA-A variance after the nuisance fit"), [], None
    members = [np.flatnonzero(cg == g) for g in range(len(groups))]
    size = np.asarray([len(m) for m in members])
    n_edges = np.bincount(eg, minlength=len(groups))
    group_islet = np.asarray([g["islet"] for g in groups])

    def metric(r):
        ss = np.bincount(cg, weights=r * r, minlength=len(groups))
        product = np.bincount(eg, weights=r[ei] * r[ej], minlength=len(groups))
        return np.divide(size * product, n_edges * ss, out=np.full(len(groups), np.nan),
                         where=(ss > 1e-8) & (n_edges > 0))

    observed = metric(residual)
    null = np.full((null_n, len(groups)), np.nan)
    for b in range(null_n):
        shuffled = residual.copy()
        for m in members:
            shuffled[m] = rng.permutation(residual[m])
        refit = fitted + shuffled
        null[b] = metric(refit - x @ fit(refit))
    mean, sd = np.nanmean(null, axis=0), np.nanstd(null, axis=0)
    exchangeable = np.isfinite(observed) & np.isfinite(mean) & (sd > 1e-8)
    rows = []
    for islet in np.unique(group_islet[exchangeable]):
        keep = (group_islet == islet) & exchangeable
        ob = float(observed[keep].mean())
        nu = null[:, keep].mean(axis=1)
        rows.append({"islet_uuid": islet, "observed": ob, "null_mean": float(nu.mean()), "null_sd": float(nu.std()),
                     "endpoint": ob - float(nu.mean()), "n_cells": int(size[keep].sum()),
                     "n_edges": int(n_edges[keep].sum()), "n_subtypes": int(keep.sum())})
    out.update(n_islets=len(rows), n_groups=int(exchangeable.sum()), null_replicates=null_n)
    if not rows:
        return dict(out, reason="no islet whose residual graph null varies"), [], None
    islets = np.unique(group_islet[exchangeable])
    null_mean_by_replicate = np.asarray([np.mean([null[b, (group_islet == i) & exchangeable].mean() for i in islets])
                                         for b in range(null_n)])
    effect = float(np.mean([r["endpoint"] for r in rows]))
    observed_mean = float(np.mean([r["observed"] for r in rows]))
    centre = null_mean_by_replicate.mean()
    out.update(donor_effect=effect, status="available", reason="",
               conditional_null_p=float((1 + np.sum(np.abs(null_mean_by_replicate - centre) >=
                                                     abs(observed_mean - centre) - 1e-12)) / (null_n + 1)))
    display = sub.assign(hla_residual=residual)
    return out, rows, display


def h4_config(frame, null_n, rng):
    keep = (frame.T_neighbour_count > 0) & np.isfinite(frame[["pdl1", "hla", "pd1", "cyto"]]).all(axis=1)
    z = frame.loc[keep].reset_index(drop=True)
    block = pd.factorize(z.islet, sort=True)[0]
    out = {"donor_effect": np.nan, "status": "unavailable", "n_cells": len(z), "n_islets": int(len(np.unique(block)))}
    if not len(z):
        return dict(out, reason="no refined endocrine cell with a measured T neighbour within the radius"), []
    if z.pdl1.nunique() < 2:
        return dict(out, reason="endocrine PD-L1 does not vary"), []
    h, p, c, y = rankz(z.hla), rankz(z.pd1), rankz(z.cyto), rankz(z.pdl1)
    nuisance = z[["pos", "pos2", "nuc", "area", "tcount", "immune", "density", "vascular", "neural",
                  "collagen", "ecad", "epcam"]].to_numpy(float)
    subtype = pd.get_dummies(z.subtype, drop_first=True, dtype=float).to_numpy()
    reduced = np.column_stack([h, p, c, h * p, h * c, p * c, nuisance, subtype])
    xr, kept = design(reduced, block)
    yy = center(y, block)
    target = center(h * p * c, block)
    sd = target.std()
    if not np.isfinite(sd) or sd < 1e-8:
        return dict(out, reason="the triple interaction does not vary within islets"), []
    target = target / sd
    w = 1 / np.bincount(block)[block]
    reduced_fit = solver(xr, w)
    t_residual = target - xr @ reduced_fit(target)
    partial = float(np.average(t_residual * t_residual, weights=w))
    if partial < 1e-8:
        return dict(out, reason="the interaction is aliased with lower-order and nuisance terms"), []
    xf = np.column_stack([xr, target])
    if np.linalg.matrix_rank(xf) < xf.shape[1]:
        return dict(out, reason="rank-deficient full model"), []
    full_fit = solver(xf, w)
    coefficient = full_fit(yy)
    predicted = xr @ reduced_fit(yy)
    residual = yy - predicted
    nb = int(block.max()) + 1
    null = np.asarray([float(full_fit(center(predicted + residual * rng.choice([-1, 1], nb)[block], block))[-1])
                       for _ in range(null_n)])
    observed, baseline = float(coefficient[-1]), float(null.mean())
    out.update(donor_effect=observed - baseline, observed_coefficient=observed, null_mean=baseline,
               null_sd=float(null.std()), status="available", reason="", null_replicates=null_n,
               conditional_null_p=float((1 + np.sum(np.abs(null - baseline) >= abs(observed - baseline) - 1e-12)) / (null_n + 1)),
               nuisance_rank=int(np.linalg.matrix_rank(xr)), full_rank=int(np.linalg.matrix_rank(xf)),
               interaction_partial_variance=partial)
    original_sd = center(reduced, block).std(axis=0)
    raw = np.zeros(reduced.shape[1])
    raw[kept] = coefficient[:-1] / original_sd[kept]
    triple = coefficient[-1] / sd
    for hh, cc, tag in ((-1, -1, "lowHLA_lowCyto"), (1, -1, "highHLA_lowCyto"),
                        (-1, 1, "lowHLA_highCyto"), (1, 1, "highHLA_highCyto")):
        out[f"PD1_simple_slope_{tag}"] = float(raw[1] + raw[3] * hh + raw[5] * cc + triple * hh * cc)
        out[f"PD1_pairwise_slope_{tag}"] = float(raw[1] + raw[3] * hh + raw[5] * cc)
    out["PD1_slope_high_high_minus_low_low"] = out["PD1_simple_slope_highHLA_highCyto"] - out["PD1_simple_slope_lowHLA_lowCyto"]
    blocks = z.groupby("islet").size()
    return out, [{"islet_uuid": i, "n_cells": int(n)} for i, n in blocks.items()]


def t_scope(data) -> bool:
    return bool(data["t_cell"].any())


REQUIRED = {"h3": ("HLA_A",), "h4": ("HLA_A", "PD_L1", "PD_1")}


def scope_reason(task, data, measured) -> str:
    """Why the donor has no computable H3/H4 at all, or empty when it has one."""
    if not t_scope(data):
        return "no refined T-lineage cell: the T exposure scope is absent"
    missing = [m for m in REQUIRED[task["operation"]] if not np.isfinite(measured[m]).any()]
    return f"required marker not measured: {', '.join(missing)}" if missing else ""


def h3_donor(task, data):
    measured = {m: data["cells"][f"redsea_aligned_{m}"].to_numpy(float) for m in NUISANCE_MARKERS}
    rows, islets, example = [], [], None
    frames = {}
    for name, radius, degree in H3_CONFIGS:
        row = {"hypothesis": "H3", "configuration": name, "neighbourhood_um": radius, "graph_k": degree,
               "edge_cap_um": 2 * data["ell"]}
        reason = scope_reason(task, data, measured)
        if reason:
            rows.append(dict(row, donor_effect=np.nan, status="unavailable", reason=reason))
            continue
        if radius not in frames:
            frame = endocrine_frame(data, measured, radius)
            frames[radius] = frame.loc[np.isfinite(frame[H3_VARIABLES].to_numpy(float)).all(axis=1)].reset_index(drop=True)
        rng = np.random.default_rng(seed(task["seed"], task["donor"], "H3", name))
        result, islet_rows, display = h3_config(frames[radius], data["outlines"], degree, 2 * data["ell"], task["null"], rng)
        rows.append(dict(row, **result))
        if name == "primary":
            islets = [dict(r, donor_id=task["donor"], disease_status=task["disease"], hypothesis="H3",
                           configuration=name) for r in islet_rows]
            if display is not None and islet_rows:
                areas = {uid: data["outlines"][uid].area for uid in display.islet.unique()}
                available = [r["islet_uuid"] for r in islet_rows]
                median = np.median([areas[u] for u in available])
                chosen = min(available, key=lambda u: (abs(areas[u] - median), u))
                cells = display.loc[display.islet.eq(chosen)]
                example = {"islet_uuid": chosen, "outline_wkt": data["outlines"][chosen].wkt,
                           "xy": cells[["x", "y"]].to_numpy(), "subtype": cells.subtype.to_numpy(str),
                           "hla_residual": cells.hla_residual.to_numpy()}
    return islets, [donor_row(task, data, row) for row in rows], example, {}


def h4_donor(task, data):
    measured = {m: data["cells"][f"redsea_aligned_{m}"].to_numpy(float) for m in NUISANCE_MARKERS}
    rows, islets = [], []
    for name, radius in H4_CONFIGS:
        row = {"hypothesis": "H4", "configuration": name, "neighbourhood_um": radius}
        reason = scope_reason(task, data, measured)
        if reason:
            rows.append(dict(row, donor_effect=np.nan, status="unavailable", reason=reason))
            continue
        frame = endocrine_frame(data, measured, radius)
        frame = frame.loc[np.isfinite(frame[H3_VARIABLES].to_numpy(float)).all(axis=1)].reset_index(drop=True)
        rng = np.random.default_rng(seed(task["seed"], task["donor"], "H4", name))
        result, blocks = h4_config(frame, task["null"], rng)
        rows.append(dict(row, **result))
        if name == "primary":
            islets = [dict(b, donor_id=task["donor"], disease_status=task["disease"], hypothesis="H4",
                           configuration=name) for b in blocks]
    return islets, [donor_row(task, data, row) for row in rows], None, {}


# ---------------------------------------------------------------- per donor

def donor_row(task, data, row):
    return {"donor_id": task["donor"], "disease_status": task["disease"], "acq_run": data["acq_run"],
            "stain_group": data["stain_group"], "n_true_islets": len(data["outlines"]),
            "mean_cell_diameter_um": data["ell"], **row}


def summarize(task, data, rows, configuration, fields):
    available = [row for row in rows if row["status"] == "available"]
    reasons = pd.Series([row["reason"] for row in rows if row["status"] != "available"]).value_counts()
    out = {"hypothesis": task["operation"].upper(), "configuration": configuration,
           "n_islets": len(rows), "n_islets_available": len(available),
           "donor_effect": float(np.mean([row["endpoint"] for row in available])) if available else np.nan,
           "status": "available" if available else "unavailable",
           "reason": "" if available else ("; ".join(f"{k} ({v})" for k, v in reasons.items()) or "no true islet"),
           "unavailable_islet_reasons": json.dumps(reasons.to_dict())}
    for field in fields:
        values = [row[field] for row in available if np.isfinite(row.get(field, np.nan))]
        out[field] = float(np.mean(values)) if values else np.nan
    return donor_row(task, data, out)


OPERATIONS = {"h1": (h1_donor, ()), "h2": (h2_donor, ()), "h3": (h3_donor, NUISANCE_MARKERS),
              "h4": (h4_donor, NUISANCE_MARKERS)}


def donor_task(task: dict) -> dict:
    started = time.perf_counter()
    function, markers = OPERATIONS[task["operation"]]
    data = load(task, markers)
    islets, donor_rows, example, audit = function(task, data)
    folder = Path(task["stage"])
    if islets:
        pd.DataFrame(islets).to_parquet(folder / "islets" / f"donor_id={task['donor']}.parquet", index=False)
    if example is not None:
        (folder / "examples" / f"donor_id={task['donor']}.json").write_text(json.dumps(jsonable(example)))
    print(f"{task['operation']} {task['donor']}: " + "; ".join(
        f"{r['configuration']} {r['status']} {r.get('donor_effect')}" for r in donor_rows), flush=True)
    return {"donor_rows": donor_rows,
            "audit": {"donor_id": task["donor"], "cells_and_outlines": data["audit"], **audit},
            "resource": {"donor_id": task["donor"], "rows": task["rows"], "peak_rss_gib": peak_rss_gib(),
                         "seconds": time.perf_counter() - started}}


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def run(context: Context, stage: Path, operation: str) -> dict:
    protocol = frozen(context)
    statistics = protocol["statistics"]
    images = geometry.image_records(context)
    objects = audit_binding(context, geometry.QUALIFIED_OBJECTS)["path"]
    (stage / "islets").mkdir()
    (stage / "examples").mkdir()
    tasks = [dict(row, operation=operation, stage=str(stage), objects=objects, pixel=images[row["native"]]["pixel_um"],
                  files={k: str(v) for k, v in geometry.donor_files(context, row["native"]).items()},
                  seed=statistics["seed"], null=statistics["conditional_null_replicates"])
             for row in donors(context)]
    results = parallel_map(donor_task, tasks, GIB_PER_TASK[operation], operation,
                           cost=[task["rows"] for task in tasks])
    endpoints = pd.DataFrame([row for result in results for row in result["donor_rows"]])
    endpoints.to_csv(stage / "donor_endpoints.csv", index=False)
    pd.DataFrame([result["resource"] for result in results]).to_csv(stage / "resource_usage.csv", index=False)
    write_json(stage / "audit.json", jsonable([result["audit"] for result in results]))
    primary = endpoints.loc[endpoints.configuration.eq("primary")]
    return {"hypothesis": operation.upper(), "donors": len(primary),
            "donors_available": int(primary.status.eq("available").sum()),
            "unavailable": {row.donor_id: row.reason for row in primary.itertuples() if row.status != "available"}}


def h1(context: Context, stage: Path) -> dict:
    return run(context, stage, "h1")


def h2(context: Context, stage: Path) -> dict:
    return run(context, stage, "h2")


def h3(context: Context, stage: Path) -> dict:
    return run(context, stage, "h3")


def h4(context: Context, stage: Path) -> dict:
    return run(context, stage, "h4")
