"""Stage 03 adapter for recorded same-slide QuPath cell-measurement exports.

This reads the existing full-cell exports as measured evidence for a bounded
pilot. It preserves native object IDs and raw marker means. It does not claim
REDSEA compensation, UniFORM normalization, or marker positivity.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from .. import gates
from ..artifacts import sha256_file
from ..same_slide_sources import load_pairs
from ..tasks import Input, Task
from .pilot_inventory import _manifest

STAGE = "03_protein_measurement"
_MARKER = re.compile(r"^Cell: (.+): Mean$")
_NUCLEUS = "Nucleus: Area µm^2"
_COLUMNS = ("acquisition_id", "protein_entry", "donor_id", "tissue", "cells",
            "eligible", "no_nucleus", "missing_nucleus", "markers", "source_sha256")


class PilotProteinError(ValueError):
    """A recorded pilot export cannot be admitted without changing its meaning."""


def _selected(project) -> list:
    pairs = load_pairs(_manifest(project))
    names = project.settings.get("protein_entries")
    if names is None:
        return list(pairs)
    if (not isinstance(names, list) or not names
            or not all(isinstance(name, str) and name for name in names)
            or len(set(names)) != len(names)):
        raise PilotProteinError("settings.protein_entries needs distinct section entry names")
    chosen = [pair for pair in pairs if pair.protein_entry in names]
    if len(chosen) != len(names):
        raise PilotProteinError("protein_entries names a section outside Stage 00's source pairs")
    return chosen


def settings(project) -> dict:
    return {"protein_entries": [pair.protein_entry for pair in _selected(project)],
            "scope": "selected" if project.settings.get("protein_entries") else "all_accepted"}


def inputs(project) -> list[Input]:
    inventory = project.release_dir / "00_inventory"
    found = [Input("inventory_acquisitions", inventory / "acquisitions.csv"),
             Input("inventory_receipt", inventory / "receipt.json"),
             Input("inventory_acceptance", inventory / "gate/accepted.json"),
             Input("same_slide_source_manifest", _manifest(project))]
    for pair in _selected(project):
        found.append(Input(f"protein_export:{pair.protein_entry}", pair.protein_export))
    return found


def _header(path: Path, expected_sha: str) -> tuple[list[str], dict[str, str]]:
    with gzip.open(path, "rb") as handle:
        raw = handle.readline()
    if hashlib.sha256(raw).hexdigest() != expected_sha:
        raise PilotProteinError(f"{path}: full-cell export header differs from the recorded panel")
    columns = raw.decode("utf-8-sig").rstrip("\r\n").split("\t")
    if len(columns) != len(set(columns)) or not {"object_id", "x_um", "y_um", _NUCLEUS} <= set(
            columns):
        raise PilotProteinError(f"{path}: missing or repeated identity/geometry columns")
    markers = {match.group(1): column for column in columns
               if (match := _MARKER.fullmatch(column))}
    if not markers or len(markers) != sum(_MARKER.fullmatch(c) is not None for c in columns):
        raise PilotProteinError(f"{path}: no unique cell-mean marker columns")
    return columns, markers


def _inventory_identity(task: Task) -> tuple[dict[str, dict], dict[str, str]]:
    with task.input("inventory_acquisitions").open(newline="") as handle:
        rows = {row["acquisition_id"]: row for row in csv.DictReader(handle)}
    receipt = json.loads(task.input("inventory_receipt").read_text())
    digests = {item["name"].removeprefix("protein_export:"): item.get("sha256")
               for item in receipt["inputs"] if item["name"].startswith("protein_export:")}
    return rows, digests


def run(task: Task) -> dict:
    pairs = load_pairs(task.input("same_slide_source_manifest"))
    chosen = set(task.settings["protein_entries"])
    selected = [pair for pair in pairs if pair.protein_entry in chosen]
    if len(selected) != len(chosen):
        raise PilotProteinError("selected protein sections changed after input declaration")
    source_spec = yaml.safe_load(task.input("same_slide_source_manifest").read_text())
    expected_header = source_spec["protein_header_sha256"]
    acquisitions, digests = _inventory_identity(task)
    rows = []
    common_markers: set[str] | None = None
    for pair in selected:
        acquisition = f"phenocycler:{pair.protein_entry}"
        recorded = acquisitions.get(acquisition)
        digest = digests.get(pair.protein_entry)
        if (not recorded or recorded["platform"] != "phenocycler"
                or recorded["donor_id"] != pair.donor_id
                or recorded["tissue"] != pair.tissue or not digest):
            raise PilotProteinError(f"{pair.protein_entry}: Stage 00 identity or digest missing")
        source = task.input(f"protein_export:{pair.protein_entry}")
        if source != pair.protein_export:
            raise PilotProteinError(f"{pair.protein_entry}: export changed after declaration")
        if sha256_file(source) != digest:
            raise PilotProteinError(f"{pair.protein_entry}: export differs from accepted Stage 00")
        _, marker_columns = _header(source, expected_header)
        markers = set(marker_columns)
        if common_markers is not None and markers != common_markers:
            raise PilotProteinError("selected protein sections have different marker columns")
        common_markers = markers
        schema = pa.schema([
            pa.field("cell_uid", pa.string()), pa.field("object_id", pa.string()),
            pa.field("acquisition_id", pa.string()), pa.field("donor_id", pa.string()),
            pa.field("tissue", pa.string()), pa.field("x_um", pa.float64()),
            pa.field("y_um", pa.float64()), pa.field("nucleus_area_um2", pa.float64()),
            pa.field("protein_eligible", pa.bool_()),
            *(pa.field(f"marker:{m}", pa.float64()) for m in sorted(markers)),
        ])
        output = task.output(f"cells/{pair.protein_entry}.parquet")
        seen: set[str] = set()
        cells = eligible = no_nucleus = missing_nucleus = 0
        with pq.ParquetWriter(output, schema, compression="zstd") as writer:
            for chunk in pd.read_csv(source, sep="\t", chunksize=50_000, dtype={"object_id": str}):
                object_ids = chunk["object_id"].fillna("").astype(str)
                if (object_ids.eq("").any() or object_ids.duplicated().any()
                        or seen.intersection(object_ids)):
                    raise PilotProteinError(f"{pair.protein_entry}: blank or repeated object_id")
                seen.update(object_ids)
                numeric = {name: pd.to_numeric(chunk[name], errors="raise").to_numpy(float)
                           for name in ("x_um", "y_um", _NUCLEUS, *marker_columns.values())}
                if not np.isfinite(numeric["x_um"]).all() or not np.isfinite(
                        numeric["y_um"]).all():
                    raise PilotProteinError(f"{pair.protein_entry}: invalid cell coordinates")
                if any(np.isinf(numeric[name]).any() for name in (_NUCLEUS,
                                                                  *marker_columns.values())):
                    raise PilotProteinError(f"{pair.protein_entry}: infinite cell measurement")
                nucleus = numeric[_NUCLEUS]
                if (nucleus[np.isfinite(nucleus)] < 0).any():
                    raise PilotProteinError(f"{pair.protein_entry}: negative nucleus area")
                admitted = np.isfinite(nucleus) & (nucleus > 0)
                cells += len(chunk)
                eligible += int(admitted.sum())
                no_nucleus += int(np.count_nonzero(nucleus == 0))
                missing_nucleus += int(np.count_nonzero(np.isnan(nucleus)))
                payload = {
                    "cell_uid": [f"{acquisition}:{digest}:{item}" for item in object_ids],
                    "object_id": object_ids.tolist(),
                    "acquisition_id": [acquisition] * len(chunk),
                    "donor_id": [pair.donor_id] * len(chunk),
                    "tissue": [pair.tissue] * len(chunk),
                    "x_um": numeric["x_um"], "y_um": numeric["y_um"],
                    "nucleus_area_um2": nucleus, "protein_eligible": admitted,
                    **{f"marker:{marker}": numeric[column]
                       for marker, column in marker_columns.items()},
                }
                writer.write_table(pa.Table.from_pydict(payload, schema=schema))
        if cells == 0:
            raise PilotProteinError(f"{pair.protein_entry}: empty cell export")
        rows.append({"acquisition_id": acquisition, "protein_entry": pair.protein_entry,
                     "donor_id": pair.donor_id, "tissue": pair.tissue, "cells": cells,
                     "eligible": eligible, "no_nucleus": no_nucleus,
                     "missing_nucleus": missing_nucleus, "markers": len(markers),
                     "source_sha256": digest})
    with task.output("sections.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    _figure(task, rows)
    return {"sections": len(rows), "accepted_sections": len(pairs),
            "scope": task.settings["scope"], "cells": sum(r["cells"] for r in rows),
            "eligible": sum(r["eligible"] for r in rows),
            "no_nucleus": sum(r["no_nucleus"] for r in rows),
            "missing_nucleus": sum(r["missing_nucleus"] for r in rows),
            "markers": len(common_markers or ())}


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
           xlabel="QuPath cells in the export (n)",
           title=f"Spleen protein exports: {len(rows)} sections, "
                 f"{rows[0]['markers']} measured markers")
    ax.invert_yaxis()
    ax.legend(fontsize=8)
    fig.savefig(gates.figure_path(task, "protein_cells.png"), dpi=150)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [
        gates.Check("selected protein sections were ingested", summary["sections"],
                    str(len(task.settings["protein_entries"])),
                    summary["sections"] == len(task.settings["protein_entries"])),
        gates.Check("cell nucleus accounting", summary["cells"],
                    "eligible + no nucleus + missing nucleus",
                    summary["cells"] == summary["eligible"] + summary["no_nucleus"]
                    + summary["missing_nucleus"]),
        gates.Check("scope of the accepted source inventory", summary["sections"],
                    f"of {summary['accepted_sections']} sections", True, blocking=False),
        gates.Check("REDSEA and normalized marker values", "not yet supplied",
                    "required before typing", False, blocking=False),
    ]
    return gates.Gate(STAGE, "Raw full-cell QuPath measurements and nucleus admission; "
                      "no REDSEA, normalization or typing.", checks, ["protein_cells.png"])
