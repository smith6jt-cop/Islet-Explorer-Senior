"""Stage 03 adapter for dedicated pLN v3 QuPath exports.

This bounded pilot preserves native object identity and raw compartment means
for review. It does not substitute for the approved own-compartment REDSEA
chain when a release is built. Section-local marker calls require Stage 05.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .. import gates, roots
from ..pln_sources import PLNExport, PLNSourceError, load_exports
from ..tasks import Input, Task

STAGE = "03_protein_measurement"
ROSTER = roots.INTEGRATION / "pln_node_annotations.csv"
INDEX = roots.REPO / "core/multiplex5000/sites/pln_qupath_exports_v3.csv"
NUCLEUS = "Nucleus: Area µm^2"
X = "Centroid X µm"
Y = "Centroid Y µm"
SUMMARY_COLUMNS = ("acquisition_id", "protein_entry", "donor_id", "roi", "tissue",
                   "cells", "eligible", "no_nucleus", "missing_nucleus", "markers",
                   "segmentation_id", "source_path")


class PLNProteinPilotError(ValueError):
    """A bounded pLN export cannot be admitted as review-only evidence."""


def _chosen(project) -> tuple[PLNExport, ...]:
    names = project.settings.get("pln_export_images")
    if (not isinstance(names, list) or not names
            or not all(isinstance(name, str) and name for name in names)
            or len(set(names)) != len(names)):
        raise PLNProteinPilotError("settings.pln_export_images needs distinct roster image names")
    root = project.source("pln_v3_export_cells")
    if root is None:
        raise PLNProteinPilotError("site sources.pln_v3_export_cells is required")
    roster = project.source("pln_export_roster") or ROSTER
    index = project.source("pln_export_index") or INDEX
    return load_exports(roster, index, source_root=root, selected_images=set(names))


def settings(project) -> dict:
    chosen = _chosen(project)
    return {"pln_export_images": [item.image_file for item in chosen],
            "scope": "bounded_raw_qupath_candidate"}


def inputs(project) -> list[Input]:
    roster = project.source("pln_export_roster") or ROSTER
    index = project.source("pln_export_index") or INDEX
    found = [Input("pln_roster", roster), Input("pln_export_index", index)]
    found.extend(Input(f"pln_export:{item.image_file}", item.source_path)
                 for item in _chosen(project))
    return found


def _entry(item: PLNExport) -> str:
    return f"{item.donor_id}__{item.roi}"


def _markers(source: pq.ParquetFile) -> dict[str, str]:
    names = source.schema_arrow.names
    if len(names) != len(set(names)):
        raise PLNProteinPilotError("QuPath export repeats a column name")
    found = {name[len("Cell: "):-len(": Mean")]: name for name in names
             if name.startswith("Cell: ") and name.endswith(": Mean")}
    found.pop("DAPI", None)
    if not found or len(found) != sum(name.startswith("Cell: ") and name.endswith(": Mean")
                                     and name != "Cell: DAPI: Mean" for name in names):
        raise PLNProteinPilotError("QuPath export has no unique non-DAPI Cell Mean markers")
    return found


def _schema(markers: dict[str, str]) -> pa.Schema:
    return pa.schema([
        pa.field("cell_uid", pa.string()), pa.field("object_id", pa.string()),
        pa.field("acquisition_id", pa.string()), pa.field("donor_id", pa.string()),
        pa.field("tissue", pa.string()), pa.field("x_um", pa.float64()),
        pa.field("y_um", pa.float64()), pa.field("nucleus_area_um2", pa.float64()),
        pa.field("protein_eligible", pa.bool_()),
        *(pa.field(f"marker:{marker}", pa.float64()) for marker in sorted(markers)),
    ])


def run(task: Task) -> dict:
    selected = _chosen(task.project)
    if {item.image_file for item in selected} != set(task.settings["pln_export_images"]):
        raise PLNProteinPilotError("selected pLN images changed after input declaration")
    rows, entries, marker_sets = [], set(), set()
    for item in selected:
        entry = _entry(item)
        if entry in entries:
            raise PLNProteinPilotError(f"more than one image names section {entry}")
        entries.add(entry)
        path = task.input(f"pln_export:{item.image_file}")
        if path != item.source_path:
            raise PLNProteinPilotError(f"{item.image_file}: export changed after declaration")
        source = pq.ParquetFile(path)
        markers = _markers(source)
        marker_sets.add(frozenset(markers))
        schema = _schema(markers)
        columns = ["Image", "Object ID", X, Y, NUCLEUS, *markers.values()]
        acquisition = f"phenocycler:{item.image_file}"
        seen: set[str] = set()
        cells = eligible = no_nucleus = missing_nucleus = 0
        with pq.ParquetWriter(task.output(f"cells/{entry}.parquet"), schema,
                              compression="zstd") as writer:
            for batch in source.iter_batches(batch_size=50_000, columns=columns):
                chunk = batch.to_pydict()
                declared = set(chunk["Image"])
                if declared not in ({item.image_file},
                                    {f"{item.image_file} - resolution #1"}):
                    raise PLNProteinPilotError(f"{item.image_file}: mixed or wrong Image values")
                object_ids = chunk["Object ID"]
                if (any(not isinstance(oid, str) or not oid for oid in object_ids)
                        or len(set(object_ids)) != len(object_ids)
                        or seen.intersection(object_ids)):
                    raise PLNProteinPilotError(f"{item.image_file}: blank or repeated Object ID")
                seen.update(object_ids)
                numeric = {name: np.asarray(chunk[name], dtype=float)
                           for name in (X, Y, NUCLEUS, *markers.values())}
                if not np.isfinite(numeric[X]).all() or not np.isfinite(numeric[Y]).all():
                    raise PLNProteinPilotError(f"{item.image_file}: invalid cell coordinates")
                if any(np.isinf(values).any() for name, values in numeric.items()
                       if name not in (X, Y)):
                    raise PLNProteinPilotError(f"{item.image_file}: infinite measurement")
                nucleus = numeric[NUCLEUS]
                if (nucleus[np.isfinite(nucleus)] < 0).any():
                    raise PLNProteinPilotError(f"{item.image_file}: negative nucleus area")
                admitted = np.isfinite(nucleus) & (nucleus > 0)
                cells += len(object_ids)
                eligible += int(admitted.sum())
                no_nucleus += int(np.count_nonzero(nucleus == 0))
                missing_nucleus += int(np.count_nonzero(np.isnan(nucleus)))
                payload = {
                    "cell_uid": [f"{acquisition}:{item.segmentation_id}:{oid}"
                                 for oid in object_ids],
                    "object_id": object_ids, "acquisition_id": [acquisition] * len(object_ids),
                    "donor_id": [item.donor_id] * len(object_ids),
                    "tissue": ["pancreatic_lymph_node"] * len(object_ids),
                    "x_um": numeric[X], "y_um": numeric[Y],
                    "nucleus_area_um2": nucleus, "protein_eligible": admitted,
                    **{f"marker:{marker}": numeric[column]
                       for marker, column in markers.items()},
                }
                writer.write_table(pa.Table.from_pydict(payload, schema=schema))
        if cells != source.metadata.num_rows or not cells:
            raise PLNProteinPilotError(f"{item.image_file}: Parquet row count changed or is empty")
        rows.append({"acquisition_id": acquisition, "protein_entry": entry,
                     "donor_id": item.donor_id, "roi": item.roi,
                     "tissue": "pancreatic_lymph_node", "cells": cells,
                     "eligible": eligible, "no_nucleus": no_nucleus,
                     "missing_nucleus": missing_nucleus, "markers": len(markers),
                     "segmentation_id": item.segmentation_id, "source_path": str(path)})
    with task.output("sections.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    _figure(task, rows)
    return {"sections": len(rows), "cells": sum(row["cells"] for row in rows),
            "eligible": sum(row["eligible"] for row in rows),
            "no_nucleus": sum(row["no_nucleus"] for row in rows),
            "missing_nucleus": sum(row["missing_nucleus"] for row in rows),
            "marker_sets": len(marker_sets), "candidate_only": True}


def _figure(task: Task, rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, max(3, len(rows) * 0.6 + 1)), layout="constrained")
    y = range(len(rows))
    ax.barh(y, [row["eligible"] for row in rows], label="Measured nucleus")
    ax.barh(y, [row["no_nucleus"] for row in rows],
            left=[row["eligible"] for row in rows], label="No nucleus")
    ax.barh(y, [row["missing_nucleus"] for row in rows],
            left=[row["eligible"] + row["no_nucleus"] for row in rows],
            label="Nucleus missing")
    ax.set(yticks=list(y), yticklabels=[row["protein_entry"] for row in rows],
           xlabel="Dedicated pLN QuPath cells (n)",
           title=f"Raw pLN protein exports: {len(rows)} physical sections")
    ax.invert_yaxis()
    ax.legend(fontsize=8)
    fig.savefig(gates.figure_path(task, "protein_cells.png"), dpi=150)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [
        gates.Check("selected pLN exports ingested", summary["sections"],
                    str(len(task.settings["pln_export_images"])),
                    summary["sections"] == len(task.settings["pln_export_images"])),
        gates.Check("nucleus accounting", summary["cells"],
                    "eligible + no nucleus + missing nucleus",
                    summary["cells"] == summary["eligible"] + summary["no_nucleus"]
                    + summary["missing_nucleus"]),
        gates.Check("common measured marker set", summary["marker_sets"], "1",
                    summary["marker_sets"] == 1),
        gates.Check("REDSEA compensation", "not supplied", "required for canonical release",
                    False, blocking=False),
    ]
    return gates.Gate(STAGE, "Dedicated pLN raw QuPath means for review only; native Object ID "
                      "and nucleus admission preserved. No REDSEA, normalization or typing.",
                      checks, ["protein_cells.png"])
