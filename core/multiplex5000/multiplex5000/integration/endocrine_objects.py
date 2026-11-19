"""All-size endocrine candidates using a pinned PubApp whole-cell constructor.

RNA supplies proposals from exact named cells; protein supplies image classifier
proposals. The *same* QuPath constructor attaches whole nucleated segmentations,
resolves ownership, splits components, and designates Single/Cluster/Islet.
This module never treats an anatomical category as a cell identity.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import uuid
import os
import subprocess

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Polygon, mapping

from .structures import dbscan_labels

CLASSES = ("Islet", "EndocrineCluster", "EndocrineSingle")


def area_only(g):
    """Preserve polygonal parts after validity repair; discard dangling lines."""
    if g.geom_type in ("Polygon", "MultiPolygon"):
        return g
    parts = [area_only(p) for p in getattr(g, "geoms", [])]
    parts = [p for p in parts if not p.is_empty and p.area > 0]
    return shapely.union_all(parts) if parts else Polygon()


def digest(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def verify_source(path, expected):
    if digest(path) != expected:
        raise ValueError(f"Pinned source changed: {path}")


def execute_constructor(*, config_path, project_path, qupath, source_hashes, log_path):
    """Execute pinned source bytes through QuPath; no writes to source projects.

    ``project_path`` must be a new scratch project whose data directory is absent.
    All declared code/input hashes are checked before execution. Outputs have to
    be absent, and the constructor receipt as well as process status must pass.
    """
    cfg = json.loads(Path(config_path).read_text())
    for path, expected in source_hashes.items():
        verify_source(path, expected)
    wrapper = Path(__file__).with_name("construct_whole_cells.groovy")
    for key in ("constructor_script", "cells_geojson", "proposals_geojson"):
        if cfg[key] not in source_hashes:
            raise ValueError(f"Unpinned constructor input: {key}")
    if str(wrapper) not in source_hashes:
        raise ValueError("Unpinned constructor wrapper")
    if (Path(project_path).parent / "data").exists():
        raise ValueError("Constructor requires a fresh project without source qpdata")
    for key in ("output_geojson", "output_precise", "output_receipt"):
        if Path(cfg[key]).exists():
            raise FileExistsError(f"Refusing to overwrite {cfg[key]}")
    env = {**os.environ, "_JAVA_OPTIONS": "-Xmx24g -XX:ActiveProcessorCount=2"}
    with Path(log_path).open("w") as log:
        proc = subprocess.run([str(qupath), "script", "-p", str(project_path),
            "-a", str(config_path), str(wrapper)], env=env, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"QuPath constructor failed; inspect {log_path}")
    if not Path(cfg["output_receipt"]).exists():
        raise RuntimeError(f"QuPath produced no successful receipt; inspect {log_path}")
    receipt = json.loads(Path(cfg["output_receipt"]).read_text())
    if not receipt.get("ok"):
        raise RuntimeError("Constructor did not report successful completion")
    return {"exit_code": proc.returncode, "source_hashes": source_hashes,
            "outputs": {k: digest(cfg[k]) for k in ("output_geojson", "output_precise", "output_receipt")}}


def native_uuid(acquisition_id, native_id):
    """Reversible through the mandatory sidecar, never a replacement native key."""
    return str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([acquisition_id, str(native_id)])))


def vertex_polygons(vertices):
    """Vendor vertices in recorded order; no centroid-derived cell outlines."""
    required = ["cell_id", "vertex_x", "vertex_y"]
    v = vertices[required]
    if v.cell_id.isna().any() or not np.isfinite(v[["vertex_x", "vertex_y"]]).all().all():
        raise ValueError("Missing native ID or nonfinite cell vertex")
    ids = v.cell_id.astype(str).to_numpy()
    starts = np.r_[0, np.flatnonzero(ids[1:] != ids[:-1]) + 1]
    if len(np.unique(ids[starts])) != len(starts):
        raise ValueError("Cell vertex records must be contiguous and uniquely keyed")
    xy = v[["vertex_x", "vertex_y"]].to_numpy(float)
    geometries = []
    for lo, hi in zip(starts, np.r_[starts[1:], len(v)]):
        g = Polygon(xy[lo:hi]) if hi - lo >= 3 else Polygon()
        if not g.is_valid:
            g = area_only(shapely.make_valid(g))
        if g.is_empty or g.area <= 0 or g.geom_type not in ("Polygon", "MultiPolygon"):
            raise ValueError(f"Unavailable whole-cell area: {ids[lo]}")
        geometries.append(g)
    return pd.DataFrame({"cell_id": ids[starts], "geometry": geometries})


def rna_proposals(cells, geometries, *, acquisition_id, eps_um=50.):
    """Keep singleton and sparse groups; hulls enclose measured seed-cell areas.

    The existing 50-um neighbor radius is retained. min_samples=1 ensures no
    eligible seed disappears as density noise. The constructor subsequently
    applies whole-cell ownership and connected-piece rules. Proposal groups are
    estimates, not histological contours. Every original native key is retained.
    """
    if not np.isfinite(eps_um) or eps_um <= 0:
        raise ValueError("Positive finite radius required")
    if cells.cell_id.isna().any() or cells.cell_id.astype(str).duplicated().any():
        raise ValueError("Unique native cell IDs required")
    if not pd.api.types.is_bool_dtype(cells.seed_evidence) or cells.seed_evidence.isna().any():
        raise ValueError("Explicit Boolean RNA seed evidence required")
    if geometries.cell_id.duplicated().any():
        raise ValueError("Duplicate native geometry IDs")
    admitted = cells.nucleus_count.gt(0) & cells.depth.gt(0)
    seeds = cells.loc[cells.seed_evidence & admitted].copy().sort_values("cell_id")
    if not np.isfinite(seeds[["x_um", "y_um"]]).all().all():
        raise ValueError("Finite seed coordinates required")
    seeds = seeds.merge(geometries, on="cell_id", how="left", validate="one_to_one")
    if seeds.geometry.isna().any():
        raise ValueError("Seed cell segmentation missing; no inferred circle fallback")
    labels = dbscan_labels(seeds[["x_um", "y_um"]].to_numpy(float), eps_um, 1) if len(seeds) else []
    seeds["proposal_group"] = labels
    rows = []
    for _, group in seeds.groupby("proposal_group", sort=True):
        members = sorted(group.cell_id.astype(str))
        uid = native_uuid(acquisition_id, json.dumps(members))
        g = shapely.union_all(group.geometry.to_list()).convex_hull
        # A single must remain its actual cell, including nonconvexity and holes.
        if len(group) == 1:
            g = group.geometry.iloc[0]
        rows.append(dict(proposal_id=uid, n_seed_cells=len(group), geometry=g,
                         seed_cell_ids=members))
    return pd.DataFrame(rows, columns=["proposal_id", "n_seed_cells", "geometry", "seed_cell_ids"]), seeds


def mean_caliper_um(g):
    """Same minimum-circle / minimum-caliper mean used by the PubApp constructor."""
    hull = g.convex_hull
    c = np.asarray(hull.exterior.coords)[:-1]
    edges = np.roll(c, -1, axis=0) - c
    normals = np.column_stack([-edges[:, 1], edges[:, 0]])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    width = np.ptp(c @ normals.T, axis=0).min()
    return float((2 * shapely.minimum_bounding_radius(g) + width) / 2)


def designate(n_whole_cells, mean_caliper):
    if n_whole_cells < 1 or int(n_whole_cells) != n_whole_cells:
        raise ValueError("Positive integral construction cell count required")
    if not np.isfinite(mean_caliper) or mean_caliper <= 0:
        raise ValueError("Measured positive diameter required")
    return "EndocrineSingle" if n_whole_cells == 1 else ("Islet" if mean_caliper >= 35 else "EndocrineCluster")


def proposal_features(proposals, scale=1.):
    from shapely.affinity import scale as scale_geometry
    return [dict(type="Feature", id=r.proposal_id,
                 geometry=mapping(scale_geometry(r.geometry, xfact=scale, yfact=scale, origin=(0, 0))),
                 properties=dict(objectType="annotation", classification=dict(name="Islet")))
            for r in proposals.itertuples()]


def validate_constructed(objects, cell_geometries, *, pixel_size_um):
    """Independent native-ID construction checks; return exact member links.

    This uses whole-area support, never nearest-centroid assignment. Object
    counts are construction counts and stay distinct from RNA endocrine counts.
    """
    from shapely.geometry import shape
    from shapely.affinity import scale
    features = objects["features"] if isinstance(objects, dict) else objects
    features = [f for f in features if f["properties"]["classification"]["name"] in CLASSES]
    geos = np.array([shapely.make_valid(scale(shape(f["geometry"]), xfact=pixel_size_um,
                     yfact=pixel_size_um, origin=(0, 0))) for f in features], dtype=object)
    cg = np.array(cell_geometries.geometry.to_list(), dtype=object)
    oi, ci = shapely.STRtree(cg).query(geos, predicate="intersects") if len(geos) and len(cg) else (np.array([],int), np.array([],int))
    shares = shapely.area(shapely.intersection(geos[oi], cg[ci])) / shapely.area(cg[ci])
    members = shares >= .95
    links = pd.DataFrame({"object_id": [str(features[i]["id"]) for i in oi[members]],
                          "cell_id": cell_geometries.cell_id.iloc[ci[members]].to_numpy(),
                          "cell_area_fraction": shares[members]})
    count = links.groupby("object_id").size()
    rows = []
    for f, g in zip(features, geos):
        uid = str(f["id"])
        props = f["properties"]
        measurements = props.get("measurements", {})
        if isinstance(measurements, list):
            measurements = {r["name"]: r["value"] for r in measurements}
        n = int(measurements["Whole cells"])
        diameter = mean_caliper_um(g)
        cls = props["classification"]["name"]
        if designate(n, diameter) != cls:
            raise ValueError(f"Constructor category disagrees with independent diameter: {uid}")
        measured = int(count.get(uid, 0))
        single_difference = None
        single_boundary_error = None
        serialization_tolerance = None
        if cls == "EndocrineSingle" and measured == 1:
            cell_id = links.loc[links.object_id.eq(uid), "cell_id"].iloc[0]
            cell = cell_geometries.loc[cell_geometries.cell_id.eq(cell_id), "geometry"].iloc[0]
            single_difference = float(g.symmetric_difference(cell).area / cell.area)
            # QuPath 0.7 GeometryTools uses a .01-unit precision model and
            # float32 coordinate storage. Preserve the raw discrepancy
            # and use the bound implied by those operations, not an area-percent
            # threshold fitted to these results. Both boundary directions count.
            extent = max(abs(v / pixel_size_um) for v in cell.bounds)
            serialization_tolerance = float(np.sqrt(2) * (.005 + np.spacing(np.float32(extent)) / 2) * pixel_size_um + 1e-8)
            single_boundary_error = float(g.boundary.hausdorff_distance(cell.boundary))
        rows.append(dict(structure_id=uid, object_class=cls, n_whole_cells=n,
                         n_area_verified_members=measured, construction_count_matches=measured == n,
                         single_relative_area_difference=single_difference, mean_caliper_um=diameter,
                         single_boundary_error_um=single_boundary_error,
                         single_serialization_tolerance_um=serialization_tolerance,
                         area_um2=float(g.area), geometry_wkb=g.wkb))
    return pd.DataFrame(rows), links, dict(n_objects=len(rows),
        n_multiply_owned_cells=int(links.cell_id.value_counts().gt(1).sum()),
        n_cut_cells=int(len(np.unique(ci[(shares > .05) & (shares < .95)]))),
        n_construction_count_disagreements=sum(not r["construction_count_matches"] for r in rows),
        n_single_geometry_disagreements=sum(r["single_boundary_error_um"] is None or
            r["single_boundary_error_um"] > r["single_serialization_tolerance_um"] for r in rows if r["object_class"] == "EndocrineSingle"))


def proposal_lineage(objects, proposals):
    """Retain every positive-area original proposal relationship after splitting.

    Whole-cell pickup may extend outside a classifier region. This is provenance,
    not a requirement for high projected overlap or a claim of endocrine identity.
    """
    og=np.asarray(objects.geometry.to_list(),dtype=object)
    pg=np.asarray(proposals.geometry.to_list(),dtype=object)
    oi,pi=shapely.STRtree(pg).query(og,predicate='intersects') if len(og) and len(pg) else (np.array([],int),np.array([],int))
    areas=shapely.area(shapely.intersection(og[oi],pg[pi]));positive=areas>0
    return pd.DataFrame(dict(structure_id=objects.structure_id.to_numpy()[oi[positive]],
        original_proposal_id=proposals.proposal_id.to_numpy()[pi[positive]],intersection_area_um2=areas[positive],
        fraction_of_object=areas[positive]/shapely.area(og[oi[positive]])))
