"""Serial Integration, MALDI route: spots placed on the neighbouring slide's objects, never on its
cells (2026-10-31).

A MALDI section is a serial section of the slide it neighbours: different cells. Its spots are
carried into the neighbour's frame by the accepted Serial Registration and related to the
neighbour's QuPath objects only -- by where each spot's centre falls, and by its continuous signed
distance to the nearest islet (no bands, no expansion rings).

For each registered section, on a PhenoCycler pancreas neighbour (the site's
``typing_export_root``: ``typing_islets/Islets__<image>.geojson`` and
``structures/<Blood_vessel|Lymphatic|Nerve>/<class>__<image>.geojson``; and the
``pancreas_annotations`` Tissue and pLN outlines):

* **per spot** (``data/spots/tissue_id=*/data.parquet``) -- its position in the neighbour's frame
  (um), ``signed_dist_islet_um`` (negative inside an islet object, as QuPath's signed distance is),
  the islet and the structure object its centre falls in, the annotation region, and the
  registration's outline RMSE as its uncertainty;
* **per object** (``data/objects/tissue_id=*/data.parquet``) -- the spots whose centre falls in it,
  and, where it is ``evaluable``, each feature's detection fraction (a 0 is not detected and
  counts) and its median among the spots that detected it. An object is evaluable only if it
  still holds a spot centre after shrinking it by the registration's outline RMSE plus half a
  spot pitch: an object smaller than that cannot be told from its surroundings at this
  registration, and is reported ``unavailable`` with that reason -- never dropped.

A section next to a dedicated pLN slide has no QuPath objects here: its spots are carried and its
object table is empty, with the reason. Every table carries ``relationship = serial_section``; no
column names a cell. The gate is automatic: integrity checks, and a figure per section.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .. import gates, maldi, profiles
from ..integration.tissues import PANCREAS
from ..integration.transform import Transform
from ..tasks import Input, Task

STAGE = "08_integration"
STRUCTURES = ("Blood_vessel", "Lymphatic", "Nerve")
REGIONS = ("Tissue", "pLN")
RELATIONSHIP = "serial_section"


class MaldiIntegrationError(ValueError):
    """A registered MALDI section cannot be placed on its neighbour's objects as recorded."""


def _release(project) -> Path:
    return Path(project.release_dir)


def _rows(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _tag(image_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", image_id)


def _registered(project) -> list[dict]:
    rows = _rows(_release(project) / "07_registration/registration.csv")
    return [r for r in rows if r["status"] == "registered"]


def _object_inputs(project, row: dict) -> list[Input]:
    if row["neighbour_platform"] != "phenocycler" or row["neighbour_tissue"] != PANCREAS:
        return []
    root, annotations = project.source("typing_export_root"), project.source("pancreas_annotations")
    if root is None or annotations is None:
        raise MaldiIntegrationError("the typing_export_root and pancreas_annotations sources "
                                    "must be named")
    tag = _tag(row["neighbour_acquisition"].split(":", 1)[1])
    tid = row["tissue_id"]
    found = [Input(f"islets:{tid}", Path(root) / "typing_islets" / f"Islets__{tag}.geojson"),
             Input(f"annotations:{tid}", Path(annotations) / f"{tag}.geojson")]
    found += [Input(f"{s}:{tid}", Path(root) / "structures" / s / f"{s}__{tag}.geojson")
              for s in STRUCTURES]
    return found


def settings(project) -> dict:
    return {"route": "maldi_serial", "relationship": RELATIONSHIP, "structures": list(STRUCTURES),
            "evaluable": "a spot centre remains after shrinking the object by the outline RMSE "
                         "plus half a spot pitch"}


def inputs(project) -> list[Input]:
    release = _release(project)
    found = []
    for stage in ("04b_maldi_ingest", "07_registration"):
        found += [Input(f"{stage}_receipt", release / stage / "receipt.json"),
                  Input(f"{stage}_gate", release / stage / "gate" / "gate.json"),
                  Input(f"{stage}_acceptance", release / stage / "gate" / "accepted.json")]
    found += [Input("registration", release / "07_registration/registration.csv"),
              Input("features", release / "04b_maldi_ingest/features.csv")]
    for row in _registered(project):
        tid = row["tissue_id"]
        found += [Input(f"transform:{tid}", release / f"07_registration/transforms/{tid}/transform.json"),
                  Input(f"spots:{tid}",
                        release / f"04b_maldi_ingest/data/spots/tissue_id={tid}/data.parquet")]
        found += _object_inputs(project, row)
    return found


def _objects(path: Path, pixel_um: float, wanted=None) -> list[dict]:
    """Each polygon feature of a QuPath GeoJSON export, in um: id, name, class and geometry."""
    import shapely

    data = json.loads(Path(path).read_text())
    features = data.get("features", []) if isinstance(data, dict) else data
    out = []
    for f in features:
        props = f.get("properties") or {}
        name = (props.get("classification") or {}).get("name", "")
        if wanted is not None and name not in wanted:
            continue
        geometry = shapely.geometry.shape(f["geometry"])
        if geometry.is_empty or geometry.geom_type not in ("Polygon", "MultiPolygon"):
            continue
        geometry = shapely.affinity.scale(geometry, pixel_um, pixel_um, origin=(0, 0))
        if not geometry.is_valid:
            geometry = shapely.make_valid(geometry)
        out.append(dict(object_id=str(f.get("id") or ""), object_name=props.get("name") or "",
                        object_class=name, size_class=props.get("size_class") or "",
                        geometry=geometry))
    return out


def _signed_distance(points, islets: list[dict]) -> np.ndarray:
    import shapely

    if not islets:
        return np.full(len(points), np.nan)
    union = shapely.union_all([o["geometry"] for o in islets])
    distance = shapely.distance(points, union.boundary)
    inside = shapely.contains_xy(union, *np.asarray([shapely.get_x(points), shapely.get_y(points)]))
    return np.where(inside, -distance, distance)


def _containing(points, objects: list[dict]) -> np.ndarray:
    """The index of the object each point's centre falls in (the first, by input order), else -1."""
    import shapely

    found = np.full(len(points), -1)
    if not objects:
        return found
    tree = shapely.STRtree([o["geometry"] for o in objects])
    point_index, object_index = tree.query(points, predicate="within")
    order = np.lexsort((object_index, point_index))
    point_index, object_index = point_index[order], object_index[order]
    first = np.unique(point_index, return_index=True)[1]
    found[point_index[first]] = object_index[first]
    return found


def _object_table(objects: list[dict], points_xy: np.ndarray, values: np.ndarray,
                  names: list[str], shrink_um: float, base: dict) -> list[dict]:
    import shapely

    out = []
    points = shapely.points(points_xy) if len(points_xy) else np.array([])
    tree = shapely.STRtree(points) if len(points_xy) else None
    for o in objects:
        inside = (tree.query(o["geometry"], predicate="contains") if tree is not None
                  else np.array([], int))
        shrunk = o["geometry"].buffer(-shrink_um)
        core = (tree.query(shrunk, predicate="contains")
                if tree is not None and not shrunk.is_empty else np.array([], int))
        row = dict(base, object_class=o["object_class"], object_id=o["object_id"],
                   object_name=o["object_name"], size_class=o["size_class"],
                   area_um2=float(o["geometry"].area), n_spots=int(len(inside)),
                   n_spots_after_rmse=int(len(core)), evaluable=bool(len(core)))
        row["reason"] = ("" if row["evaluable"] else
                         "no spot centre inside" if not len(inside) else
                         "smaller than a spot after the registration RMSE")
        chosen = values[inside] if row["evaluable"] else np.empty((0, len(names)))
        for k, name in enumerate(names):
            column = chosen[:, k]
            column = column[np.isfinite(column)]
            detected = column[column > 0]
            row[f"detected_fraction__{name}"] = (float((column > 0).mean())
                                                 if column.size else np.nan)
            row[f"median_detected__{name}"] = float(np.median(detected)) if detected.size else np.nan
        out.append(row)
    return out


def run(task: Task) -> dict:
    import shapely

    pixel_um = float(profiles.load_all()[("platform", "phenocycler")].pixel_size_um)
    names = [r["feature"] for r in _rows(task.input("features"))]
    shown, sections = [], []
    for row in _registered(task.project):
        tid = row["tissue_id"]
        transform = Transform.load(task.input(f"transform:{tid}").parent)
        spots = pq.read_table(task.input(f"spots:{tid}")).to_pandas()
        on = spots["on_tissue"].to_numpy(bool)
        xy = transform.apply(spots[["x_um", "y_um"]].to_numpy(float))
        points = shapely.points(xy)
        rmse = float(row["outline_rmse_um"])
        base = dict(tissue_id=tid, donor_id=row["donor_id"],
                    neighbour_acquisition=row["neighbour_acquisition"], relationship=RELATIONSHIP)
        out = pd.DataFrame(dict(spot_uid=spots["spot_uid"], x_neighbour_um=xy[:, 0],
                                y_neighbour_um=xy[:, 1], on_tissue=on,
                                registration_rmse_um=rmse, **{k: v for k, v in base.items()
                                                              if k != "tissue_id"}))
        objects, islets = [], []
        reason = ""
        if f"islets:{tid}" in task.inputs:
            islets = _objects(task.input(f"islets:{tid}"), pixel_um, {"Islet"})
            structures = [o for s in STRUCTURES
                          for o in _objects(task.input(f"{s}:{tid}"), pixel_um)]
            regions = _objects(task.input(f"annotations:{tid}"), pixel_um, set(REGIONS))
            out["signed_dist_islet_um"] = _signed_distance(points, islets)
            index = _containing(points, islets)
            out["islet_object_id"] = [islets[i]["object_id"] if i >= 0 else "" for i in index]
            index = _containing(points, structures)
            out["structure_class"] = [structures[i]["object_class"] if i >= 0 else ""
                                      for i in index]
            out["structure_object_id"] = [structures[i]["object_id"] if i >= 0 else ""
                                          for i in index]
            index = _containing(points, regions)
            out["region"] = [regions[i]["object_class"] if i >= 0 else "" for i in index]
            values = spots[names].to_numpy(float)
            objects = _object_table(islets + structures, xy[on], values[on], names,
                                    rmse + maldi.SPOT_PITCH_UM / 2, base)
        else:
            reason = "the neighbour is not a PhenoCycler pancreas scan; no QuPath objects here"
            for column in ("islet_object_id", "structure_class", "structure_object_id", "region"):
                out[column] = ""
            out["signed_dist_islet_um"] = np.nan
        pq.write_table(pa.Table.from_pandas(out, preserve_index=False),
                       task.output(f"data/spots/tissue_id={tid}/data.parquet"))
        table = pd.DataFrame(objects) if objects else pd.DataFrame(
            columns=[*base, "object_class", "object_id", "evaluable", "reason"])
        pq.write_table(pa.Table.from_pandas(table, preserve_index=False),
                       task.output(f"data/objects/tissue_id={tid}/data.parquet"))
        islet_rows = [o for o in objects if o["object_class"] == "Islet"]
        sections.append(dict(
            tissue_id=tid, donor_id=row["donor_id"], neighbour=row["neighbour_acquisition"],
            spots=len(out), spots_on_tissue=int(on.sum()),
            spots_in_islets=int((out["islet_object_id"][on] != "").sum()),
            islets=len(islet_rows), islets_evaluable=sum(o["evaluable"] for o in islet_rows),
            structures=len(objects) - len(islet_rows),
            structures_evaluable=sum(o["evaluable"] for o in objects
                                     if o["object_class"] != "Islet"),
            registration_rmse_um=rmse, no_objects_reason=reason))
        shown.append((tid, xy[on], out["signed_dist_islet_um"].to_numpy()[on], islets))
    _write(task.output("sections.csv"), sections)
    _figure(task, shown)
    return {"sections": len(sections),
            "spots_on_tissue": sum(s["spots_on_tissue"] for s in sections),
            "by_section": {s["tissue_id"]: f"{s['islets_evaluable']}/{s['islets']} islets, "
                                           f"{s['structures_evaluable']}/{s['structures']} "
                                           f"structures evaluable; RMSE "
                                           f"{s['registration_rmse_um']:.0f} um"
                           for s in sections},
            "without_objects": sorted(f"{s['tissue_id']}: {s['no_objects_reason']}"
                                      for s in sections if s["no_objects_reason"])}


def _write(path: Path, rows: list[dict]) -> None:
    columns: list[str] = []
    for row in rows:
        columns += [k for k in row if k not in columns]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or ["empty"])
        writer.writeheader()
        writer.writerows(rows)


def _figure(task: Task, shown: list) -> None:
    """Each section's on-tissue spots in the neighbour's frame, coloured by signed distance to the
    nearest islet (clipped at +-500 um for display only), with the islet outlines."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = max(len(shown), 1)
    per_row = min(n, 4)
    rows = -(-n // per_row)
    fig, axes = plt.subplots(rows, per_row, figsize=(4.2 * per_row, 4.4 * rows), squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for ax, (tid, xy, distance, islets) in zip(axes.ravel(), shown):
        finite = np.isfinite(distance)
        ax.scatter(xy[~finite, 0], xy[~finite, 1], s=0.4, c="0.7", linewidths=0)
        art = ax.scatter(xy[finite, 0], xy[finite, 1], s=0.4, linewidths=0,
                         c=np.clip(distance[finite], -500, 500), cmap="RdBu_r", vmin=-500, vmax=500)
        for o in islets:
            geometry = o["geometry"]
            for polygon in getattr(geometry, "geoms", [geometry]):
                ax.plot(*polygon.exterior.xy, color="black", linewidth=0.3)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(tid, fontsize=8)
        if finite.any():
            fig.colorbar(art, ax=ax, shrink=0.6, label="signed distance to islet (um)")
    fig.suptitle("MALDI spots on the neighbouring slide's islet objects (serial sections: "
                 "different cells)", fontsize=10)
    fig.tight_layout()
    fig.savefig(gates.figure_path(task, "maldi_integration.png"), dpi=110)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    C = gates.Check
    checks = [
        C("every registered MALDI section is placed on its neighbour", summary["sections"],
          "reported", summary["sections"] > 0,
          detail="spots carried by the accepted transform; objects by where a spot's centre falls"),
        C("objects evaluable, by section", summary["by_section"], "reported", True,
          blocking=False, detail="an object is evaluable only if a spot centre remains after "
                                 "shrinking it by the outline RMSE plus half a spot pitch"),
        C("sections with no QuPath objects", summary["without_objects"], "reported", True,
          blocking=False),
    ]
    return gates.Gate(STAGE, "MALDI spots on the neighbouring slide's objects and their signed "
                      "distance to the nearest islet; serial sections, never the same cells.",
                      checks, ["maldi_integration.png"])
