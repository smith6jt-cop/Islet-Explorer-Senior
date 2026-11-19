#!/usr/bin/env python3
"""Inspect existing pLN annotation sources without selecting any MALDI neighbour."""
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote, urlparse

import shapely
import yaml

BASE = Path(__file__).resolve().parents[2]
PLN = Path("/home/smith6jt/IO60pLN")
OUT = BASE / "outputs/maldi_completion_20261102"


def ref(path):
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"path": str(path), "sha256": digest}


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def main():
    project_path = PLN / "project.qpproj"
    annotations_path = BASE / "core/phenocycler/data/integration/pln_node_annotations.csv"
    inventory_path = BASE / "docs/evaluations/dynamic_portability_20261102/inventories/io60_acquisitions.csv"
    audit_path = PLN / "results/v3/export/node_geojson/export_audit.tsv"
    profile_path = BASE / "core/multiplex5000/profiles/platforms/phenocycler.yaml"
    registration_pixel = yaml.safe_load(profile_path.read_text())["pixel_size_um"]
    images = json.loads(project_path.read_text())["images"]
    by_entry = {str(image["entryID"]): image for image in images}
    assert len(by_entry) == len(images)
    inventory = [r for r in rows(inventory_path) if r["platform"] == "phenocycler" and r["tissue"] == "pancreatic_lymph_node"]
    by_file = {r["qptiff"]: r for r in inventory}
    assert len(by_file) == len(inventory)
    output = []
    for declaration in rows(annotations_path):
        entry = declaration["qupath_entry"]
        image = by_entry[entry]
        builder = image["serverBuilder"]
        filename = Path(unquote(urlparse(builder["uri"]).path)).name
        assert filename == declaration["image_file"]
        acquisition = by_file[filename]
        source = PLN / f"results/v3/export/node_geojson/entry{entry}.geojson"
        collection = json.loads(source.read_text())
        selected_classes = declaration["node_classes"].split(";")
        features = [f for f in collection["features"]
                    if f.get("properties", {}).get("classification", {}).get("name") in selected_classes]
        native_ids = [f.get("id") for f in features]
        assert len(features) == int(declaration["n_nodes"])
        assert all(native_ids) and len(set(native_ids)) == len(native_ids)
        geometry = [shapely.geometry.shape(f["geometry"]) for f in features]
        calibration = builder["metadata"]["pixelCalibration"]
        output.append({"acquisition_id": acquisition["acquisition_id"], "image_file": filename,
                       "qupath_entry": entry, "donor_id": declaration["donor_id"],
                       "source": ref(source), "native_object_ids": native_ids,
                       "selected_classes": selected_classes, "objects": len(features),
                       "geometry_units": "level-zero image pixels", "project_pixel_calibration": calibration,
                       "registration_profile_um_per_pixel": registration_pixel,
                       "calibration_note": "Project metadata and active registration profile use different documented calibrations; a future consumer must bind one target coordinate convention explicitly.",
                       "declared_by": declaration["recorded_by"], "declared_on": declaration["recorded_on"],
                       "areas_px2": [float(g.area) for g in geometry],
                       "geometry_validity": [{"object_id": oid, "valid": bool(g.is_valid),
                                               "empty": bool(g.is_empty), "type": g.geom_type,
                                               "area_object_available": bool(g.is_valid and not g.is_empty and g.geom_type in {"Polygon", "MultiPolygon"}),
                                               "reason": str(shapely.is_valid_reason(g))}
                                              for oid, g in zip(native_ids, geometry, strict=True)],
                       "axis_aligned_rectangles": [bool(g.equals(g.envelope)) for g in geometry],
                       "source_identity_check": "pass",
                       "interpretation": "Existing declared node annotation ROIs; bounding rectangles are not tissue outlines or inferred follicle/islet objects.",
                       "maldi_neighbour_assignment": None, "registration_acceptance": None})
    assert {r["acquisition_id"] for r in output} == {r["acquisition_id"] for r in inventory}
    report = {"schema": "maldi-object-availability-audit/1", "scope": "Read-only source availability and exact QuPath-entry/image/acquisition reconciliation; no active object-source routing, MALDI pairing or acceptance",
              "acquisitions": len(output), "objects": sum(r["objects"] for r in output),
              "area_objects": sum(v["area_object_available"] for r in output for v in r["geometry_validity"]),
              "project_entries": sorted(int(x) for x in by_entry),
              "entry_count_note": "There are19 current entries/exports; the largest entry ID is21 and IDs3/8 are absent. Entry ID is not a count.",
              "sources": [ref(p) for p in (project_path, annotations_path, inventory_path, audit_path,
                                            profile_path)],
              "rows": output}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pln_object_availability.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"acquisitions": report["acquisitions"], "objects": report["objects"],
                      "invalid_geometry": [v for r in output for v in r["geometry_validity"] if not v["valid"]],
                      "unavailable_area_objects": [v for r in output for v in r["geometry_validity"] if not v["area_object_available"]],
                      "source_identity_checks": "all pass", "maldi_neighbour_assignments": 0}, indent=2))


if __name__ == "__main__":
    main()
