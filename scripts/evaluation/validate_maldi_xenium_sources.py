#!/usr/bin/env python3
"""Validate every declared Xenium morphology source without choosing MALDI neighbours.

This writes isolated research diagnostics, not release gate acceptance. Exact
region input files are content-hashed; no raw image is modified or substituted.
"""
from __future__ import annotations

import argparse
import csv
import json
import platform
import shutil
import sys
from importlib.metadata import version
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from multiplex5000.stages import maldi_serial_registration as R
from multiplex5000.tasks import Task, sha256_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory-release", type=Path, required=True)
    parser.add_argument("--source-index", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    project = SimpleNamespace(release_dir=args.inventory_release.resolve(),
                              source=lambda key: args.source_index if key == "xenium_source_index" else None)
    inventory = project.release_dir / "00_inventory/acquisitions.csv"
    with inventory.open(newline="") as handle:
        neighbours = [r for r in csv.DictReader(handle) if r["platform"] == "xenium"]
    ids = [n["acquisition_id"] for n in neighbours]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate inventory acquisition IDs")
    snapshot = args.out / "producer_maldi_serial_registration.py"
    shutil.copyfile(R.__file__, snapshot)
    script_snapshot = args.out / "producer_validate_maldi_xenium_sources.py"
    shutil.copyfile(__file__, script_snapshot)
    report = {"schema": "maldi-xenium-source-validation/1", "scope": "exact-region morphology only; no neighbour selection or accepted registration",
              "producer_source": str(snapshot.resolve()), "producer_sha256": sha256_file(snapshot),
              "validation_script": str(script_snapshot.resolve()), "validation_script_sha256": sha256_file(script_snapshot),
              "environment": {"python": sys.version, "platform": platform.platform(),
                              "packages": {name: version(name) for name in ("numpy", "pandas", "tifffile", "scikit-image", "scipy")}},
              "source_index": ({"path": str(args.source_index.resolve()), "sha256": sha256_file(args.source_index)}
                               if args.source_index is not None else None),
              "inventory": str(inventory), "inventory_sha256": sha256_file(inventory),
              "declared_acquisitions": len(ids), "complete": False, "acquisitions": []}
    panels = []
    for n in neighbours:
        record = {"acquisition_id": n["acquisition_id"], "condition": n.get("condition"),
                  "accepted_registration": False, "neighbour_assignment": None}
        try:
            source = R._xenium_source(project, n)
            record.update(status=source["status"], reason=source["reason"])
            if source["status"] == "ready":
                token = n["acquisition_id"].replace(":", "_")
                row = {"tissue_id": token, "neighbour": n}
                declared = R._neighbour_inputs(project, row)
                task = Task("07_registration", "source-validation", args.out,
                            {i.name: i for i in declared}, project=project)
                before = source["image"].stat()
                fixed = R._fixed_xenium(task, row)
                digest = sha256_file(source["image"])
                after = source["image"].stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise ValueError("morphology file changed during validation")
                destination = args.out / "masks" / f"{token}.npz"
                destination.parent.mkdir(exist_ok=True)
                np.savez_compressed(destination, dapi=fixed["dapi"], mask=fixed["mask"],
                                    grid=json.dumps(fixed["mask_grid"].as_dict()))
                record.update(status="validated", image=str(source["image"]), image_sha256=digest,
                              image_bytes=after.st_size, source_geometry=fixed["frame_metadata"],
                              mask_parameters=fixed["frame_metadata"].get("mask_parameters"),
                              grid=fixed["mask_grid"].as_dict(), mask_pixels=int(fixed["mask"].sum()),
                              mask_artifact=str(destination), mask_artifact_sha256=sha256_file(destination))
                panels.append((n["acquisition_id"], fixed["dapi"], fixed["mask"]))
        except (ValueError, KeyError, OSError) as exc:
            record.update(status="unavailable", reason=f"{type(exc).__name__}: {exc}")
        report["acquisitions"].append(record)
        (args.out / "source_validation.json").write_text(json.dumps(report, indent=2) + "\n")
        print(n["acquisition_id"], record["status"], flush=True)
    report["complete"] = len(report["acquisitions"]) == len(ids)
    report["validated"] = sum(r["status"] == "validated" for r in report["acquisitions"])
    report["unavailable"] = len(ids) - report["validated"]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for start in range(0, len(panels), 8):
        fig, axes = plt.subplots(2, 4, figsize=(16, 9))
        for ax in axes.flat:
            ax.set_axis_off()
        for ax, (name, image, mask) in zip(axes.flat, panels[start:start + 8], strict=False):
            ax.imshow(np.log1p(image), cmap="gray")
            ax.contour(mask, levels=[0.5], colors=["#46c97c"], linewidths=0.6)
            ax.set_title(name, fontsize=9)
        fig.suptitle("Exact Xenium-region DAPI / existing mask at 50 µm; no MALDI relationship or registration accepted")
        fig.tight_layout()
        fig.savefig(args.out / f"source_masks_{start // 8 + 1}.png", dpi=100)
        plt.close(fig)
    (args.out / "source_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("complete", "declared_acquisitions", "validated", "unavailable")}, indent=2))


if __name__ == "__main__":
    main()
