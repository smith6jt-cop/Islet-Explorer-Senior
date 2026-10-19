"""Original true-Islet outlines and pixel-classifier objects, in micrometres, bound by hash.

The outlines are the QuPath islet annotation exports the current composition was built from; their
feature id is the cells' `endocrine_object_id`. Objects (Blood vessel, Nerve, Lymphatic) are
structural regions, never cell identities. Nothing here expands, simplifies or repairs a boundary.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import shapely
from shapely.affinity import scale
from shapely.geometry import shape

from .common import Context, binding, read_json, require, threaded

# Where v3's `subpanel_io.EXPORT` and `GEOJSON` point; the islet resolver searches it recursively.
EXPORT = Path.home() / "IO60panc2nd/results/v3/export"
ISLETS = EXPORT / "typing_islets"
OBJECTS = {"Blood vessel": EXPORT / "structures/Blood_vessel", "Nerve": EXPORT / "structures/Nerve",
           "Lymphatic": EXPORT / "structures/Lymphatic"}
QUALIFIED_OBJECTS = "protein_composition/qualified_broad_composition/qualified_broad_objects.parquet"


def qupath_manifest(context: Context) -> Path:
    return context.original_repo / "core/phenocycler/data/qupath_manifest.json"


def image_records(context: Context) -> dict[str, dict]:
    """Native partition id -> pixel size and image stem, as v3's build_islet_geometry._manifest reads it."""
    raw = read_json(qupath_manifest(context))
    out = {}
    for row in raw["images"]:
        require(row["pixel_size_um_x"] == row["pixel_size_um_y"], f"Anisotropic pixels: {row['donor_id']}")
        out[row["donor_id"]] = {"pixel_um": float(row["pixel_size_um_x"]),
                                "stem": row["image_id"].split(" - ")[0].replace(".qptiff", "")}
    return out


def only(paths: list[Path], what: str) -> Path:
    require(len(paths) == 1, f"{what}: expected exactly one export, found {len(paths)}")
    return paths[0]


def donor_files(context: Context, native: str) -> dict[str, Path]:
    stem = image_records(context)[native]["stem"]
    files = {"islets": only([p for p in ISLETS.rglob("*.geojson") if stem in p.name], f"{native} islets")}
    for name, folder in OBJECTS.items():
        files[name] = only([p for p in folder.glob("*.geojson") if stem in p.name], f"{native} {name}")
    return files


def bindings(context: Context, natives: list[str]) -> list[dict]:
    paths = [qupath_manifest(context)]
    for native in natives:
        paths.extend(donor_files(context, native).values())
    return threaded(binding, paths)


def read_features(path: Path, pixel: float) -> list[tuple[str, dict, object]]:
    out = []
    for feature in json.loads(Path(path).read_text())["features"]:
        if not feature.get("geometry"):
            continue
        out.append((str(feature.get("id", "")), feature.get("properties") or {},
                    scale(shape(feature["geometry"]), xfact=pixel, yfact=pixel, origin=(0, 0))))
    return out


def true_islets(path: Path, pixel: float, wanted: set[str]) -> tuple[dict, dict]:
    """The donor's qualified true-Islet outlines by UUID; invalid outlines are reported, not repaired."""
    outlines, audit = {}, {"source_features": 0, "invalid_or_empty": 0}
    for uid, properties, geometry in read_features(path, pixel):
        audit["source_features"] += 1
        if uid not in wanted:
            continue
        require(properties.get("size_class") == "Islet", f"Qualified Islet {uid} is not an Islet outline")
        if geometry.is_empty or not geometry.is_valid:
            audit["invalid_or_empty"] += 1
            continue
        require(uid not in outlines, f"Duplicate islet UUID {uid}")
        shapely.prepare(geometry)
        outlines[uid] = geometry
    audit["qualified_true_islets"] = len(wanted)
    audit["valid_outlines"] = len(outlines)
    audit["qualified_without_outline"] = len(wanted - set(outlines))
    return outlines, audit


def objects(path: Path, pixel: float) -> tuple[list, dict]:
    out, audit = [], {"source_features": 0, "invalid_or_empty": 0, "not_areal": 0}
    for _, _, geometry in read_features(path, pixel):
        audit["source_features"] += 1
        if geometry.is_empty or not geometry.is_valid:
            audit["invalid_or_empty"] += 1
        elif geometry.geom_type not in ("Polygon", "MultiPolygon"):
            audit["not_areal"] += 1
        else:
            out.append(geometry)
    audit["valid_objects"] = len(out)
    return out, audit


def qualified_islets(objects_path: Path, donor: str) -> pd.DataFrame:
    table = pd.read_parquet(objects_path)
    table = table.loc[table.donor_id.astype(str).eq(donor) & table.original_object_class.eq("Islet")]
    require(table.endocrine_object_id.is_unique, f"Duplicate qualified islet UUIDs: {donor}")
    return table.reset_index(drop=True)
