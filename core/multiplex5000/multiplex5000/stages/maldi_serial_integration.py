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

Dedicated pLN and Xenium neighbours may supply explicitly declared, receipt-bound area objects
through ``sources.maldi_object_sources``. Each acquisition names its frame, units and producer;
an annotation region is not promoted to an islet or a cell identity. Without declared objects,
spots are still carried and object evidence remains unavailable with its reason. Every table
carries ``relationship = serial_section``; no column names a cell.
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
from ..artifacts import sha256_file
from ..integration.tissues import PANCREAS
from ..integration.transform import Transform
from ..resolved_parameters import write_resolved_parameters
from ..tasks import Input, Task

STAGE = "08_integration"
STRUCTURES = ("Blood_vessel", "Lymphatic", "Nerve")
REGIONS = ("Tissue", "pLN")
RELATIONSHIP = "serial_section"
OBJECT_SCHEMA = "multiplex5000-maldi-object-sources/1"
OBJECT_ROLES = {"islets", "structures", "regions"}


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
    if len({r["tissue_id"] for r in rows}) != len(rows):
        raise MaldiIntegrationError("duplicate registered section identifiers")
    return [r for r in rows if r["status"] == "registered"]


def _object_declarations(project) -> dict[str, dict]:
    """Explicit acquisition sources; no donor, tissue token or filename matching."""
    path = project.source("maldi_object_sources")
    if path is None:
        return {}
    data = json.loads(Path(path).read_text())
    if data.get("schema") != OBJECT_SCHEMA or not isinstance(data.get("acquisitions"), list):
        raise MaldiIntegrationError("invalid MALDI object-source manifest schema")
    found = {}
    for entry in data["acquisitions"]:
        aid = entry.get("acquisition_id")
        if not isinstance(aid, str) or not aid or aid in found:
            raise MaldiIntegrationError("object sources require unique exact acquisition IDs")
        if entry.get("frame_id") != f"{aid}:um" or not isinstance(entry.get("sources"), list):
            raise MaldiIntegrationError(f"{aid}: object source needs its exact acquisition um frame")
        for source in entry["sources"]:
            if source.get("binding_kind") not in {"producer_output", "observed_export"}:
                raise MaldiIntegrationError(f"{aid}: source binding kind must be explicit")
            if source.get("role") not in OBJECT_ROLES or source.get("units") not in {"um", "pixel"}:
                raise MaldiIntegrationError(f"{aid}: unknown object role or coordinate units")
            if source["units"] == "pixel":
                pixel = source.get("pixel_size_um")
                if type(pixel) not in (int, float) or not np.isfinite(pixel) or pixel <= 0:
                    raise MaldiIntegrationError(f"{aid}: pixel objects need a finite positive calibration")
            classes = source.get("classes")
            if not isinstance(classes, list) or not classes or any(not isinstance(c, str) or not c for c in classes):
                raise MaldiIntegrationError(f"{aid}: object classes must be explicitly selected")
            if source["role"] == "islets" and set(classes) != {"Islet"}:
                raise MaldiIntegrationError(f"{aid}: only declared Islet objects supply islet distances")
            for key in ("path", "source_receipt"):
                if not isinstance(source.get(key), str) or not Path(source[key]).is_absolute():
                    raise MaldiIntegrationError(f"{aid}: {key} must be an absolute source path")
            for key in ("sha256", "source_receipt_sha256"):
                value = source.get(key)
                if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
                    raise MaldiIntegrationError(f"{aid}: {key} must bind source content")
        found[aid] = entry
    return found


def _object_inputs(project, row: dict, declarations=None) -> list[Input]:
    declarations = _object_declarations(project) if declarations is None else declarations
    explicit = declarations.get(row["neighbour_acquisition"])
    if explicit is not None:
        tid = row["tissue_id"]
        return [item for i, source in enumerate(explicit["sources"]) for item in (
            Input(f"object_source:{tid}:{i}", Path(source["path"])),
            Input(f"object_receipt:{tid}:{i}", Path(source["source_receipt"])))]
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
            "object_sources_schema": OBJECT_SCHEMA,
            "overlap_policy": "ambiguous unique assignment; nonexclusive object summaries",
            "denominators": "finite signal-mask spots plus explicit all-acquired-raster comparison; zeros retained",
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
    declarations = _object_declarations(project)
    if project.source("maldi_object_sources") is not None:
        found.append(Input("object_sources_manifest", Path(project.source("maldi_object_sources"))))
    for row in _registered(project):
        tid = row["tissue_id"]
        found += [Input(f"transform:{tid}", release / f"07_registration/transforms/{tid}/transform.json"),
                  Input(f"displacement:{tid}", release / f"07_registration/transforms/{tid}/transform_displacement.npz",
                        required=False),
                  Input(f"spots:{tid}",
                        release / f"04b_maldi_ingest/data/spots/tissue_id={tid}/data.parquet")]
        found += _object_inputs(project, row, declarations)
    return found


def _objects(path: Path, pixel_um: float, wanted=None, *, strict: bool = False) -> list[dict]:
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
        if strict and (not f.get("id") or props.get("objectType") in {"cell", "nucleus"}):
            raise MaldiIntegrationError("declared area objects need native IDs and cannot be cells/nuclei")
        geometry = shapely.geometry.shape(f["geometry"])
        reason = ""
        if geometry.is_empty or geometry.geom_type not in ("Polygon", "MultiPolygon"):
            if strict:
                reason = f"unavailable area geometry: {geometry.geom_type}"
                geometry = None
            else:
                continue
        if geometry is not None:
            geometry = shapely.affinity.scale(geometry, pixel_um, pixel_um, origin=(0, 0))
            if not geometry.is_valid:
                if strict:
                    reason, geometry = "unavailable invalid area geometry", None
                else:
                    geometry = shapely.make_valid(geometry)
        out.append(dict(object_id=str(f.get("id") or ""), object_name=props.get("name") or "",
                        object_class=name, size_class=props.get("size_class") or "",
                        geometry=geometry, unavailable_reason=reason))
    if strict and len({o["object_id"] for o in out}) != len(out):
        raise MaldiIntegrationError("duplicate native object IDs in an export")
    return out


def _accepted_stage(task: Task, stage: str) -> dict:
    """Acceptance must bind the current receipt and gate, not just exist on disk."""
    receipt_path = task.input(f"{stage}_receipt")
    receipt = json.loads(receipt_path.read_text())
    acceptance = json.loads(task.input(f"{stage}_acceptance").read_text())
    gate = json.loads(task.input(f"{stage}_gate").read_text())
    if (receipt.get("status") != "complete" or receipt.get("stage") != stage
            or acceptance.get("stage") != stage or gate.get("stage") != stage
            or acceptance.get("release") != receipt.get("release") or gate.get("release") != receipt.get("release")
            or receipt.get("release") != task.release or gate.get("passed") is not True
            or any(not c.get("passed") for c in gate.get("checks", []) if c.get("blocking", True))
            or acceptance.get("receipt_sha256") != sha256_file(receipt_path)
            or acceptance.get("gate_sha256") != sha256_file(task.input(f"{stage}_gate"))):
        raise MaldiIntegrationError(f"{stage} acceptance does not bind its current receipt and gate")
    return receipt


def _bound_output(task: Task, receipt: dict, stage: str, name: str) -> None:
    path = task.input(name)
    relative = str(path.relative_to(_release(task.project) / stage))
    outputs = {r["path"]: r.get("sha256") for r in receipt["outputs"]}
    if len(outputs) != len(receipt["outputs"]):
        raise MaldiIntegrationError(f"accepted {stage} receipt has duplicate output paths")
    if not path.is_file():
        if name.startswith("displacement:") and relative not in outputs:
            return
        raise MaldiIntegrationError(f"accepted {stage} output is missing: {relative}")
    if outputs.get(relative) != sha256_file(path):
        raise MaldiIntegrationError(f"accepted {stage} output changed: {relative}")


def _declared_objects(task: Task, row: dict, transform: Transform, entry: dict) -> dict[str, list]:
    """Load source-bound objects directly in the registered neighbour's calibrated frame."""
    aid, tid = row["neighbour_acquisition"], row["tissue_id"]
    if (entry["acquisition_id"] != aid or transform.meta.get("fixed_acquisition_id") != aid
            or transform.meta.get("fixed_frame") != entry["frame_id"]
            or row.get("fixed_frame") != entry["frame_id"]):
        raise MaldiIntegrationError(f"{tid}: object and registration acquisition/frame differ")
    loaded = {role: [] for role in OBJECT_ROLES}
    seen = set()
    for i, source in enumerate(entry["sources"]):
        path, receipt_path = task.input(f"object_source:{tid}:{i}"), task.input(f"object_receipt:{tid}:{i}")
        if sha256_file(path) != source["sha256"] or sha256_file(receipt_path) != source["source_receipt_sha256"]:
            raise MaldiIntegrationError(f"{tid}: declared object source or source receipt changed")
        receipt = json.loads(receipt_path.read_text())
        if source["binding_kind"] == "observed_export":
            if (receipt.get("schema") != "multiplex5000-object-source-observation/1"
                    or receipt.get("acquisition_id") != aid
                    or receipt.get("frame_id") != entry["frame_id"]
                    or receipt.get("geometry_kind") != "annotation"
                    or set(receipt.get("classes", [])) != set(source["classes"])
                    or (source["units"] == "pixel" and receipt.get("frame_pixel_size_um") != source["pixel_size_um"])):
                raise MaldiIntegrationError(f"{tid}: observation receipt does not name this acquisition/frame")
        bindings = []
        for binding in receipt.get("outputs", []):
            recorded = Path(binding["path"])
            recorded = recorded if recorded.is_absolute() else receipt_path.parent / recorded
            if recorded.resolve() == path.resolve():
                bindings.append(binding.get("sha256"))
        if receipt.get("status") != "complete" or bindings != [source["sha256"]]:
            raise MaldiIntegrationError(f"{tid}: source receipt does not uniquely bind the object export")
        pixel = source["pixel_size_um"] if source["units"] == "pixel" else 1.0
        if source["units"] == "pixel":
            frame = transform.meta.get("fixed_geometry", {})
            calibration = frame.get("pixel_um", frame.get("base_pixel_um"))
            if calibration != pixel:
                raise MaldiIntegrationError(f"{tid}: object calibration differs from registration frame")
        objects = _objects(path, pixel, set(source["classes"]), strict=True)
        if (source["binding_kind"] == "observed_export"
                and {o["object_id"] for o in objects} != set(receipt.get("native_object_ids", []))):
            raise MaldiIntegrationError(f"{tid}: observed source native object IDs changed")
        for obj in objects:
            if obj["object_id"] in seen:
                raise MaldiIntegrationError(f"{tid}: duplicate object ID across declared sources")
            seen.add(obj["object_id"])
            obj["object_role"] = source["role"]
        loaded[source["role"]].extend(objects)
    return loaded


def _signed_distance(points, islets: list[dict]) -> np.ndarray:
    import shapely

    if not islets or any(o["geometry"] is None for o in islets):
        return np.full(len(points), np.nan)
    union = shapely.union_all([o["geometry"] for o in islets])
    distance = shapely.distance(points, union.boundary)
    inside = shapely.contains_xy(union, *np.asarray([shapely.get_x(points), shapely.get_y(points)]))
    return np.where(inside, -distance, distance)


def _containing(points, objects: list[dict]) -> np.ndarray:
    """Unique object index, -1 outside, or -2 for overlapping/ambiguous membership."""
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
    count = np.bincount(point_index, minlength=len(points))
    found[count > 1] = -2
    return found


def _object_table(objects: list[dict], points_xy: np.ndarray, values: np.ndarray,
                  names: list[str], shrink_um: float, base: dict) -> list[dict]:
    import shapely

    if not np.isfinite(shrink_um) or shrink_um < 0:
        raise MaldiIntegrationError("object uncertainty radius must be finite and nonnegative")
    out = []
    points = shapely.points(points_xy) if len(points_xy) else np.array([])
    tree = shapely.STRtree(points) if len(points_xy) else None
    for o in objects:
        valid = o["geometry"] is not None
        inside = (tree.query(o["geometry"], predicate="contains") if tree is not None and valid
                  else np.array([], int))
        shrunk = o["geometry"].buffer(-shrink_um) if valid else None
        core = (tree.query(shrunk, predicate="contains")
                if tree is not None and shrunk is not None and not shrunk.is_empty else np.array([], int))
        row = dict(base, object_class=o["object_class"], object_id=o["object_id"],
                   object_name=o["object_name"], size_class=o["size_class"],
                   object_role=o.get("object_role", "islets" if o["object_class"] == "Islet" else "structures"),
                   area_um2=float(o["geometry"].area) if valid else np.nan, n_spots=int(len(inside)),
                   n_spots_after_rmse=int(len(core)), evaluable=bool(len(core)),
                   uncertainty_radius_um=shrink_um,
                   membership_policy="nonexclusive geometric membership; ambiguous spots have no unique object ID")
        row["reason"] = (o.get("unavailable_reason", "") if not valid else "" if row["evaluable"] else
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
            row[f"finite_spots__{name}"] = int(column.size)
            row[f"detected_spots__{name}"] = int(detected.size)
            stable = values[core, k] if row["evaluable"] else np.array([])
            stable = stable[np.isfinite(stable)]
            row[f"core_detected_fraction__{name}"] = float((stable > 0).mean()) if stable.size else np.nan
            row[f"core_finite_spots__{name}"] = int(stable.size)
        out.append(row)
    return out


def run(task: Task) -> dict:
    import shapely

    pixel_um = float(profiles.load_all()[("platform", "phenocycler")].pixel_size_um)
    accepted = {stage: _accepted_stage(task, stage) for stage in ("04b_maldi_ingest", "07_registration")}
    _bound_output(task, accepted["07_registration"], "07_registration", "registration")
    _bound_output(task, accepted["04b_maldi_ingest"], "04b_maldi_ingest", "features")
    declarations = _object_declarations(task.project)
    names = [r["feature"] for r in _rows(task.input("features"))]
    if not names or len(set(names)) != len(names):
        raise MaldiIntegrationError("feature identifiers must be nonempty and unique")
    shown, sections = [], []
    for row in _registered(task.project):
        tid = row["tissue_id"]
        for name in (f"transform:{tid}", f"displacement:{tid}"):
            _bound_output(task, accepted["07_registration"], "07_registration", name)
        _bound_output(task, accepted["04b_maldi_ingest"], "04b_maldi_ingest", f"spots:{tid}")
        transform = Transform.load(task.input(f"transform:{tid}").parent)
        aid = row["neighbour_acquisition"]
        if (transform.meta.get("fixed_acquisition_id") != aid
                or transform.meta.get("fixed_frame") != f"{aid}:um"
                or row.get("fixed_frame") != f"{aid}:um"
                or transform.meta.get("moving_acquisition_id") != f"maldi:{tid}"
                or transform.meta.get("moving_frame") != f"maldi:{tid}:stage_um"
                or row.get("moving_frame") != f"maldi:{tid}:stage_um"
                or transform.meta.get("coordinate_units") != "um"):
            raise MaldiIntegrationError(f"{tid}: transform is not bound to this exact acquisition/frame")
        spots = pq.read_table(task.input(f"spots:{tid}")).to_pandas()
        if spots.spot_uid.isna().any() or not spots.spot_uid.is_unique:
            raise MaldiIntegrationError(f"{tid}: MALDI spot identifiers are missing or repeated")
        on = spots["on_tissue"].to_numpy(bool)
        xy = transform.apply(spots[["x_um", "y_um"]].to_numpy(float))
        points = shapely.points(xy)
        rmse = float(row["outline_rmse_um"])
        if not np.isfinite(rmse) or rmse < 0 or not np.isfinite(xy).all():
            raise MaldiIntegrationError(f"{tid}: registration coordinates/RMSE are not finite and physical")
        base = dict(tissue_id=tid, donor_id=row["donor_id"],
                    neighbour_acquisition=row["neighbour_acquisition"], relationship=RELATIONSHIP)
        out = pd.DataFrame(dict(spot_uid=spots["spot_uid"], x_neighbour_um=xy[:, 0],
                                y_neighbour_um=xy[:, 1], on_tissue=on,
                                registration_rmse_um=rmse, **{k: v for k, v in base.items()
                                                              if k != "tissue_id"}))
        objects, islets, structures, regions = [], [], [], []
        reason = ""
        explicit = declarations.get(row["neighbour_acquisition"])
        if explicit is not None:
            loaded = _declared_objects(task, row, transform, explicit)
            islets, structures, regions = loaded["islets"], loaded["structures"], loaded["regions"]
        elif f"islets:{tid}" in task.inputs:
            if transform.meta.get("fixed_geometry", {}).get("pixel_um") != pixel_um:
                raise MaldiIntegrationError(f"{tid}: legacy object calibration differs from registration frame")
            islets = _objects(task.input(f"islets:{tid}"), pixel_um, {"Islet"})
            structures = [o for s in STRUCTURES
                          for o in _objects(task.input(f"{s}:{tid}"), pixel_um)]
            regions = _objects(task.input(f"annotations:{tid}"), pixel_um, set(REGIONS))
        else:
            reason = "no object sources declared for this exact neighbour acquisition"
        if not islets and not structures and not regions and not reason:
            reason = "declared exports contain no selected area objects"
        out["signed_dist_islet_um"] = _signed_distance(points, islets)
        out["signed_dist_islet_status"] = ("unavailable_no_islet_objects" if not islets else
                                          "unavailable_invalid_islet_geometry" if any(o["geometry"] is None for o in islets) else "available")
        out["object_evidence_state"] = ("unavailable" if reason else "partial_unavailable_objects"
                                        if any(o["geometry"] is None for o in islets + structures + regions) else "available")
        index = _containing(points, islets)
        out["islet_object_id"] = [islets[i]["object_id"] if i >= 0 else "" for i in index]
        out["islet_membership_state"] = np.where(index == -2, "ambiguous_overlap", np.where(index >= 0, "unique", "outside_or_unavailable"))
        index = _containing(points, structures)
        out["structure_class"] = [structures[i]["object_class"] if i >= 0 else "" for i in index]
        out["structure_object_id"] = [structures[i]["object_id"] if i >= 0 else "" for i in index]
        out["structure_membership_state"] = np.where(index == -2, "ambiguous_overlap", np.where(index >= 0, "unique", "outside_or_unavailable"))
        index = _containing(points, regions)
        out["region"] = [regions[i]["object_class"] if i >= 0 else "" for i in index]
        out["region_object_id"] = [regions[i]["object_id"] if i >= 0 else "" for i in index]
        out["region_membership_state"] = np.where(index == -2, "ambiguous_overlap", np.where(index >= 0, "unique", "outside_or_unavailable"))
        values = spots[names].to_numpy(float)
        # Existing pancreas summaries remain islets/structures. Explicit region exports may
        # additionally summarize their declared ROIs without becoming biological structures.
        summarized = islets + structures + (regions if explicit is not None else [])
        objects = _object_table(summarized, xy[on], values[on], names,
                                rmse + maldi.SPOT_PITCH_UM / 2, base)
        all_raster = _object_table(summarized, xy, values, names,
                                  rmse + maldi.SPOT_PITCH_UM / 2, base)
        for output, comparison in zip(objects, all_raster, strict=True):
            output["denominator_policy"] = "finite signal-mask spots; all-raster fields are explicit sensitivity evidence"
            for key, value in comparison.items():
                if key.startswith(("detected_fraction__", "finite_spots__", "detected_spots__", "core_")) or key in {"n_spots", "n_spots_after_rmse", "evaluable"}:
                    output["all_raster_" + key] = value
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
            structures=sum(o["object_role"] == "structures" for o in objects),
            structures_evaluable=sum(o["evaluable"] for o in objects
                                     if o["object_role"] == "structures"),
            regions=sum(o["object_role"] == "regions" for o in objects),
            regions_evaluable=sum(o["evaluable"] for o in objects if o["object_role"] == "regions"),
            object_source_mode="explicit_acquisition_manifest" if explicit is not None else "legacy_pancreas" if islets or structures or regions else "unavailable",
            fixed_frame=transform.meta.get("fixed_frame", ""),
            registration_rmse_um=rmse, no_objects_reason=reason))
        shown.append((tid, xy[on], out["signed_dist_islet_um"].to_numpy()[on],
                      [o for o in islets + structures + regions if o["geometry"] is not None]))
    unavailable = [r for r in _rows(task.input("registration")) if r["status"] != "registered"]
    _write(task.output("registration_unavailable.csv"), unavailable)
    _write(task.output("sections.csv"), sections)
    _figure(task, shown)
    summary = {"sections": len(sections), "declared_sections": len(sections) + len(unavailable),
            "registration_unavailable": [{"tissue_id": r["tissue_id"], "reason": r.get("reason", "")} for r in unavailable],
            "spots_on_tissue": sum(s["spots_on_tissue"] for s in sections),
            "by_section": {s["tissue_id"]: f"{s['islets_evaluable']}/{s['islets']} islets, "
                                           f"{s['structures_evaluable']}/{s['structures']} "
                                           f"structures; {s['regions_evaluable']}/{s['regions']} "
                                           f"declared regions evaluable; RMSE "
                                           f"{s['registration_rmse_um']:.0f} um"
                           for s in sections},
            "without_objects": sorted(f"{s['tissue_id']}: {s['no_objects_reason']}"
                                      for s in sections if s["no_objects_reason"])}
    summary["parameter_recording"] = write_resolved_parameters(
        task, operation="serial_object_summary",
        context={"modality": "maldi", "relationship_type": RELATIONSHIP, "units": "um and exported feature intensity"},
        parameters=[{"name": "object_uncertainty", "role": "selected_method", "value": "outline RMSE + half spot pitch", "units": "um"},
                    {"name": "overlap_policy", "role": "selected_method", "value": "no unique assignment for overlapping objects; nonexclusive object summaries", "units": None},
                    {"name": "denominator_policy", "role": "selected_method", "value": settings(task.project)["denominators"], "units": "finite spots"}],
        acquisitions=sections + [{"tissue_id": r["tissue_id"], "unavailable_reason": r.get("reason", "")} for r in unavailable],
        input_names=tuple(name for name in ("04b_maldi_ingest_receipt", "07_registration_receipt", "object_sources_manifest") if name in task.inputs),
        producer_paths=(Path(__file__),))
    return summary


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
    for ax, (tid, xy, distance, islets) in zip(axes.ravel(), shown, strict=False):
        ax.set_axis_on()
        ax.set_xlabel("Neighbour x (µm)", fontsize=8)
        ax.set_ylabel("Neighbour y (µm)", fontsize=8)
        ax.tick_params(labelsize=7)
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
    fig.suptitle("MALDI spots and declared objects\nSerial sections: different cells", fontsize=10)
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
