#!/usr/bin/env python3
"""Bind existing, exactly inventoried pLN annotation exports for the new object reader.

This observes existing exports; it does not claim an original producer receipt,
select a MALDI neighbour, approve registration, or turn annotations into islets.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from multiplex5000.stages import maldi_serial_integration as M
from multiplex5000.tasks import sha256_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    audit = json.loads(args.audit.read_text())
    for source in audit["sources"] + [r["source"] for r in audit["rows"]]:
        if sha256_file(Path(source["path"])) != source["sha256"]:
            raise ValueError(f"source changed since exact acquisition audit: {source['path']}")
    args.out.mkdir(parents=True, exist_ok=False)
    acquisitions, checked = [], []
    for row in audit["rows"]:
        aid, source = row["acquisition_id"], row["source"]
        frame = f"{aid}:um"
        pixel = row["registration_profile_um_per_pixel"]
        receipt = {"schema": "multiplex5000-object-source-observation/1", "status": "complete",
                   "acquisition_id": aid, "frame_id": frame, "geometry_kind": "annotation",
                   "historical_producer": None,
                   "scope": "existing source observed and hashed; no original producer or scientific acceptance inferred",
                   "qupath_entry": row["qupath_entry"], "image_file": row["image_file"],
                   "native_object_ids": row["native_object_ids"], "classes": row["selected_classes"],
                   "source_project_calibration": row["project_pixel_calibration"],
                   "frame_pixel_size_um": pixel, "source_units": "level-zero image pixels",
                   "inputs": [{"path": str(args.audit.resolve()), "sha256": sha256_file(args.audit)}, *audit["sources"]],
                   "outputs": [{"path": source["path"], "sha256": source["sha256"]}]}
        path = args.out / f"entry{row['qupath_entry']}.source_observation.json"
        path.write_text(json.dumps(receipt, indent=2) + "\n")
        acquisitions.append({"acquisition_id": aid, "frame_id": frame, "sources": [{
            "binding_kind": "observed_export", "role": "regions", "classes": row["selected_classes"],
            "units": "pixel", "pixel_size_um": pixel, "path": source["path"], "sha256": source["sha256"],
            "source_receipt": str(path.resolve()), "source_receipt_sha256": sha256_file(path)}]})
        objects = M._objects(Path(source["path"]), pixel, set(row["selected_classes"]), strict=True)
        if {o["object_id"] for o in objects} != set(row["native_object_ids"]):
            raise ValueError(f"native object universe changed: {aid}")
        checked.append({"acquisition_id": aid, "native_objects": len(objects),
                        "area_objects": sum(o["geometry"] is not None for o in objects),
                        "unavailable": [{"object_id": o["object_id"], "reason": o["unavailable_reason"]}
                                        for o in objects if o["geometry"] is None]})
    manifest = {"schema": M.OBJECT_SCHEMA, "acquisitions": acquisitions}
    (args.out / "object_sources.json").write_text(json.dumps(manifest, indent=2) + "\n")
    validation = {"schema": "maldi-object-source-validation/1", "status": "pass", "acquisitions": checked,
                  "scope": "source identity, calibration and geometry availability; no MALDI relationship or acceptance"}
    (args.out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps({"acquisitions": len(checked), "objects": sum(r["native_objects"] for r in checked),
                      "area_objects": sum(r["area_objects"] for r in checked)}, indent=2))


if __name__ == "__main__":
    main()
