"""MALDI Ingest: each inventoried MALDI section's spots, as measured, keyed and placed on its tissue
(2026-10-31).

Reads the MALDI export the accepted Source Inventory bound (the ``maldi_spots`` source: the two
feature files and their ``SHA256SUMS``), for the MALDI sections it lists (``maldi_acquisitions.csv``:
the IO60 donors' sections only). For each section it writes one Parquet partition,
``data/spots/tissue_id=<Tissue_id>/data.parquet``:

* ``spot_uid`` -- ``identity.cell_uid("maldi:<Tissue_id>", "spot50", <spotId>)``, checked unique;
* ``acquisition_id``, ``tissue_id``, ``spot_id``, ``raster``, ``x_um``, ``y_um`` (the MALDI stage
  frame, um);
* every feature as measured: 0 is *not detected*, kept as a measured value (maintainer,
  2026-10-31); a value missing from the file stays missing and is counted;
* ``on_tissue`` -- ``maldi.ON_TISSUE_METHOD``: inside the raster drawn to the tissue's outline and
  at least one ion detected; never a cutoff on a feature's value.

and ``images/<Tissue_id>.npz``: the summed-ion image and the on-tissue mask on the section's own
50 um spot grid, with that grid -- what Serial Registration reads. No intensity is normalized
here; MALDI normalization is a separate, later decision.

Its gate is automatic, as the measurement stages' are: integrity checks block, and the
indicators (spots on and off the tissue, detection per feature, missing values) are reported.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .. import gates, identity, maldi
from ..tasks import Input, Task

STAGE = "04b_maldi_ingest"
SEGMENTATION = "spot50"


class MaldiIngestError(ValueError):
    """The MALDI export cannot be ingested as the Source Inventory recorded it."""


def _inventory(project) -> Path:
    return Path(project.release_dir) / "00_inventory"


def _folder(project) -> Path:
    folder = project.source("maldi_spots")
    if folder is None:
        raise MaldiIngestError("MALDI Ingest needs the maldi_spots source (site or project)")
    return Path(folder)


def _selected(project) -> list[str] | None:
    chosen = project.settings.get("maldi_tissues")
    if chosen is None:
        return None
    if (not isinstance(chosen, list) or not chosen or len(set(chosen)) != len(chosen)
            or not all(isinstance(t, str) and t for t in chosen)):
        raise MaldiIngestError("settings.maldi_tissues needs distinct Tissue_ids")
    return chosen


def settings(project) -> dict:
    chosen = _selected(project)
    return {"maldi_tissues": chosen if chosen is not None else "all_inventoried",
            "segmentation_id": SEGMENTATION, "on_tissue_method": maldi.ON_TISSUE_METHOD,
            "normalization": "none"}


def inputs(project) -> list[Input]:
    inventory = _inventory(project)
    found = [Input("inventory_receipt", inventory / "receipt.json"),
             Input("inventory_gate", inventory / "gate" / "gate.json"),
             Input("inventory_acceptance", inventory / "gate" / "accepted.json"),
             Input("maldi_acquisitions", inventory / "maldi_acquisitions.csv"),
             Input("maldi_checksums", _folder(project) / maldi.CHECKSUMS)]
    for f in maldi.source_files(_folder(project)):
        found.append(Input(f"maldi_{f.role}", f.path, hash=False, recorded_sha256=f.sha256))
    return found


def _rows(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _bound_by_inventory(task: Task) -> None:
    """The export files must be the ones the accepted Source Inventory read: same path, size,
    modification time and recorded digest."""
    receipt = json.loads(task.input("inventory_receipt").read_text())
    recorded = {i["name"]: i for i in receipt.get("inputs", [])}
    for role, _ in maldi.ROLES:
        name = f"maldi_{role}"
        before, path = recorded.get(name), task.input(name)
        if before is None:
            raise MaldiIngestError(f"the accepted Source Inventory did not read {name}")
        stat = path.stat()
        if (before["path"] != str(path) or before.get("bytes") != stat.st_size
                or before.get("mtime_ns") != stat.st_mtime_ns
                or before.get("recorded_sha256") != task.inputs[name].recorded_sha256):
            raise MaldiIngestError(f"{name} is not the file the accepted Source Inventory read")


def run(task: Task) -> dict:
    _bound_by_inventory(task)
    files = [maldi.SourceFile(role, task.input(f"maldi_{role}"),
                              task.inputs[f"maldi_{role}"].recorded_sha256)
             for role, _ in maldi.ROLES]
    listed = {r["tissue_id"]: r for r in _rows(task.input("maldi_acquisitions"))}
    chosen = task.settings["maldi_tissues"]
    if isinstance(chosen, list):
        unknown = sorted(set(chosen) - set(listed))
        if unknown:
            raise MaldiIngestError(f"settings.maldi_tissues names sections the inventory does "
                                   f"not list: {unknown}")
        listed = {t: listed[t] for t in chosen}
    con = maldi.connect()
    feats = maldi.features(con, files)
    names = [n for _, n in feats]
    sections, detection = [], []
    for tissue_id in sorted(listed):
        maldi.safe_key(tissue_id)
        table = con.execute(maldi.spots_sql(files, feats, [tissue_id])).to_arrow_table()
        frame = table.to_pandas()
        if len(frame) != int(listed[tissue_id]["n_spots"]):
            raise MaldiIngestError(f"{tissue_id}: {len(frame)} spots read, the inventory "
                                   f"recorded {listed[tissue_id]['n_spots']}")
        if frame["mismatch"].any():
            raise MaldiIngestError(f"{tissue_id}: {int(frame['mismatch'].sum())} spots differ "
                                   "between the two export files")
        acquisition = f"maldi:{tissue_id}"
        frame.insert(0, "acquisition_id", acquisition)
        frame["segmentation_id"] = SEGMENTATION
        frame["native_cell_id"] = frame["spot_id"].astype("int64").astype(str)
        frame.insert(0, "spot_uid", identity.cell_uids(frame).to_numpy())
        identity.check_cell_keys(frame.rename(columns={"spot_uid": "cell_uid"}))
        values = frame[names].to_numpy(dtype=float)
        summed = np.nansum(values, axis=1)
        frame["on_tissue"] = maldi.on_tissue(summed)
        frame = frame.drop(columns=["mismatch", "segmentation_id", "native_cell_id"])
        out = task.output(f"data/spots/tissue_id={tissue_id}/data.parquet")
        pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), out)
        grid = maldi.spot_grid(frame["x_um"].to_numpy(), frame["y_um"].to_numpy())
        np.savez_compressed(
            task.output(f"images/{tissue_id}.npz"),
            summed=maldi.spot_image(frame["x_um"], frame["y_um"], summed, grid),
            on_tissue=maldi.spot_image(frame["x_um"], frame["y_um"],
                                       frame["on_tissue"].to_numpy(float), grid, fill=0.0) > 0,
            grid=json.dumps(grid.as_dict()))
        on = frame["on_tissue"].to_numpy()
        missing = int(np.isnan(values).sum())
        sections.append(dict(acquisition_id=acquisition, tissue_id=tissue_id,
                             donor_id=listed[tissue_id]["donor_id"], spots=len(frame),
                             on_tissue=int(on.sum()), off_tissue=int((~on).sum()),
                             missing_values=missing, rasters=int(frame["raster"].nunique())))
        for (role, name), column in zip(feats, values.T, strict=True):
            measured = column[on & np.isfinite(column)]
            detection.append(dict(tissue_id=tissue_id, role=role, feature=name,
                                  spots_on_tissue=int(on.sum()),
                                  detected=int((measured > 0).sum()),
                                  detected_fraction=(float((measured > 0).mean())
                                                     if measured.size else float("nan")),
                                  missing=int((~np.isfinite(column[on])).sum())))
    _write(task.output("sections.csv"), sections)
    _write(task.output("detection.csv"), detection)
    _write(task.output("features.csv"), [dict(role=r, feature=n) for r, n in feats])
    _figure(task, sections)
    return {"sections": len(sections), "spots": sum(s["spots"] for s in sections),
            "spots_on_tissue": sum(s["on_tissue"] for s in sections),
            "off_tissue_by_section": {s["tissue_id"]: s["off_tissue"] for s in sections
                                      if s["off_tissue"]},
            "missing_values": sum(s["missing_values"] for s in sections),
            "features": len(feats), "inventoried": len(listed),
            "lowest_detected_fraction": sorted(
                (f"{d['tissue_id']} {d['feature']}: {d['detected_fraction']:.3f}"
                 for d in detection if np.isfinite(d["detected_fraction"])
                 and d["detected_fraction"] < 0.05))[:20]}


def _write(path: Path, rows: list[dict]) -> None:
    columns: list[str] = []
    for row in rows:
        columns += [k for k in row if k not in columns]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or ["empty"])
        writer.writeheader()
        writer.writerows(rows)


def _figure(task: Task, sections: list[dict]) -> None:
    """Each section's summed-ion image with its off-tissue spots outlined: what the tissue mask
    took out, on the spots themselves."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = max(len(sections), 1)
    per_row = min(n, 6)
    rows = -(-n // per_row)
    fig, axes = plt.subplots(rows, per_row, figsize=(2.6 * per_row, 2.8 * rows), squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for ax, s in zip(axes.ravel(), sections):
        with np.load(task.output(f"images/{s['tissue_id']}.npz")) as data:
            summed, on = data["summed"], data["on_tissue"]
        shown = np.log1p(np.maximum(np.nan_to_num(summed, nan=0), 0))
        ax.imshow(np.where(np.isfinite(summed), shown, np.nan), cmap="magma",
                  interpolation="nearest")
        off = np.isfinite(summed) & ~on
        if off.any():
            rr, cc = np.nonzero(off)
            ax.scatter(cc, rr, s=2, c="#00b4d8", marker="s", linewidths=0)
        ax.set_title(f"{s['tissue_id']}\n{s['on_tissue']} on, {s['off_tissue']} off tissue",
                     fontsize=7)
    fig.suptitle("MALDI sections: summed ions; off-tissue spots (no ion detected) in blue",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(gates.figure_path(task, "maldi_spots.png"), dpi=110)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    C = gates.Check
    checks = [
        C("every inventoried MALDI section is ingested", summary["sections"],
          f"{summary['inventoried']} sections", summary["sections"] == summary["inventoried"],
          detail="spot counts equal the Source Inventory's, and the two files agree on every "
                 "spot (checked per section; a difference stops the stage)"),
        C("spots on the tissue", summary["spots_on_tissue"], "reported", True, blocking=False,
          detail=f"{summary['spots']} spots; {maldi.ON_TISSUE_METHOD}. Off-tissue by section: "
                 + (", ".join(f"{k} {v}" for k, v in summary["off_tissue_by_section"].items())
                    or "none")),
        C("values missing from the export", summary["missing_values"], "reported", True,
          blocking=False, detail="unavailable, never filled; a 0 is not detected and is kept"),
        C("features detected in under 5% of a section's on-tissue spots",
          summary["lowest_detected_fraction"], "reported", True, blocking=False,
          detail="detection.csv has every feature in every section; nothing is dropped for it"),
    ]
    return gates.Gate(STAGE, "Each MALDI section's spots as measured, keyed by spot and placed on "
                      "its tissue; no intensity normalized.", checks, ["maldi_spots.png"])
