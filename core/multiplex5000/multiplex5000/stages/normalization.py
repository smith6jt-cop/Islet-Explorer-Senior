"""Stage 05: one UniFORM shift for each measured protein marker and section.

Both same-slide and dedicated pLN pilot adapters emit the same Stage 03 cell
table contract. This stage fits one shared marker basis across selected
sections, without mixing raw and aligned values or choosing positivity cutoffs.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from .. import gates, normalize, roots
from ..same_slide_sources import load_pairs
from ..tasks import Input, Task
from .pilot_inventory import _manifest

STAGE = "05_normalization"
_SCALE_COLUMNS = ("protein_entry", "marker", "reference", "scale", "shift_log",
                  "n_cells", "mode_sd_before", "mode_sd_after")
_QC_COLUMNS = ("marker", "samples", "unavailable_samples", "reference",
               "mode_sd_before", "mode_sd_after", "mode_fold_before", "mode_fold_after")


class NormalizationError(ValueError):
    """Stage 03 evidence or the pinned normalization source is incomplete."""


def _require_pilot(project) -> None:
    if not (project.settings.get("same_slide_source_manifest")
            or project.settings.get("pln_export_images")):
        raise NormalizationError("Stage 05 needs an explicit measurement-source adapter")


def _route(project):
    if getattr(project, "cohort", None) == "io60_v3" and not (
            project.settings.get("same_slide_source_manifest")
            or project.settings.get("pln_export_images")):
        from . import io60_normalization

        return io60_normalization
    return None


def _checkout(project) -> Path:
    named = project.source("uniform_checkout")
    return Path(named) if named else roots.REPO / "external/UniFORM"


def _entries(project) -> list[str]:
    if project.settings.get("pln_export_images"):
        path = project.release_dir / "03_protein_measurement/sections.csv"
        if not path.is_file():
            raise NormalizationError("dedicated pLN Stage 03 sections are not published")
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        entries = [row["protein_entry"] for row in rows]
        if (not entries or len(set(entries)) != len(entries)
                or len(entries) != len(project.settings["pln_export_images"])):
            raise NormalizationError("dedicated pLN Stage 03 section roster is ambiguous")
        return entries
    pairs = load_pairs(_manifest(project))
    selected = project.settings.get("protein_entries")
    if selected is None:
        return [pair.protein_entry for pair in pairs]
    if not set(selected) <= {pair.protein_entry for pair in pairs}:
        raise NormalizationError("Stage 05 names a protein entry outside Stage 00")
    return list(selected)


def settings(project) -> dict:
    route = _route(project)
    if route is not None:
        return route.settings(project)
    _require_pilot(project)
    return {"uniform_commit": normalize.UNIFORM_COMMIT,
            "uniform_checkout": str(_checkout(project)),
            "protein_entries": _entries(project), "num_bins": normalize.NUM_BINS}


def inputs(project) -> list[Input]:
    route = _route(project)
    if route is not None:
        return route.inputs(project)
    _require_pilot(project)
    stage = project.release_dir / "03_protein_measurement"
    found = [Input("protein_receipt", stage / "receipt.json"),
             Input("protein_sections", stage / "sections.csv"),
             Input("protein_acceptance", stage / "gate/accepted.json")]
    if project.settings.get("same_slide_source_manifest"):
        found.append(Input("same_slide_source_manifest", _manifest(project)))
    for entry in _entries(project):
        found.append(Input(f"protein_cells:{entry}", stage / "cells" / f"{entry}.parquet"))
    package = _checkout(project) / "UniFORM"
    if not package.is_dir():
        raise NormalizationError(f"pinned UniFORM package missing at {package}")
    found += [Input(f"uniform_code:{path.name}", path)
              for path in sorted(package.glob("*.py"))]
    return found


def _table_markers(path: Path) -> set[str]:
    names = pq.ParquetFile(path).schema_arrow.names
    if not {"cell_uid", "object_id", "protein_eligible"} <= set(names):
        raise NormalizationError(f"{path}: Stage 03 identity or admission columns missing")
    markers = {name.removeprefix("marker:") for name in names if name.startswith("marker:")}
    if not markers:
        raise NormalizationError(f"{path}: no measured marker columns")
    return markers


def _values(path: Path, marker: str) -> np.ndarray:
    table = pq.read_table(path, columns=["protein_eligible", f"marker:{marker}"])
    eligible = table["protein_eligible"].to_numpy(zero_copy_only=False).astype(bool)
    # Arrow nulls remain NaN: missing measurements are never treated as zero.
    values = table[f"marker:{marker}"].to_pandas().to_numpy(dtype=float)
    return values[eligible]


def run(task: Task) -> dict:
    route = _route(task.project)
    if route is not None:
        return route.run(task)
    import contextlib
    import io

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    entries = task.settings["protein_entries"]
    with task.input("protein_sections").open(newline="") as handle:
        sections = {row["protein_entry"]: row for row in csv.DictReader(handle)}
    if set(sections) != set(entries):
        raise NormalizationError("Stage 03 sections differ from the selected normalization roster")
    receipt = json.loads(task.input("protein_receipt").read_text())
    if receipt.get("status") != "complete" or receipt.get("stage") != "03_protein_measurement":
        raise NormalizationError("Stage 03 has no complete measurement receipt")
    paths = {entry: task.input(f"protein_cells:{entry}") for entry in entries}
    marker_sets = {entry: _table_markers(path) for entry, path in paths.items()}
    if len({frozenset(markers) for markers in marker_sets.values()}) != 1:
        raise NormalizationError("protein sections have different measured marker sets")
    markers = sorted(next(iter(marker_sets.values())))
    checkout = Path(task.settings["uniform_checkout"])
    registration, normalization = normalize.load_uniform(checkout)
    scales, qc = [], []
    for marker in markers:
        values = {entry: _values(path, marker) for entry, path in paths.items()}
        measurable = {entry: data for entry, data in values.items()
                      if np.any(np.isfinite(data) & (data >= 1.0))}
        if len(measurable) < 2:
            qc.append({"marker": marker, "samples": len(measurable),
                       "unavailable_samples": len(entries) - len(measurable),
                       "reference": "", "mode_sd_before": "", "mode_sd_after": "",
                       "mode_fold_before": "", "mode_fold_after": ""})
            continue
        with contextlib.redirect_stdout(io.StringIO()):
            fitted = normalize.fit_marker(marker, measurable, registration, normalization)
        fig = plt.gcf()
        figure_path = task.output(f"gate/markers/{marker.replace('/', '_')}.png")
        fig.savefig(figure_path, dpi=110, bbox_inches="tight")
        plt.close(fig)
        qc.append({"marker": marker, "samples": len(measurable),
                   "unavailable_samples": len(entries) - len(measurable),
                   "reference": fitted.reference,
                   "mode_sd_before": fitted.mode_sd_before,
                   "mode_sd_after": fitted.mode_sd_after,
                   "mode_fold_before": fitted.mode_fold_before,
                   "mode_fold_after": fitted.mode_fold_after})
        for entry in entries:
            scales.append({"protein_entry": entry, "marker": marker,
                           "reference": fitted.reference,
                           "scale": fitted.scale.get(entry, ""),
                           "shift_log": fitted.shift_log.get(entry, ""),
                           "n_cells": fitted.n_cells.get(entry, 0),
                           "mode_sd_before": fitted.mode_sd_before,
                           "mode_sd_after": fitted.mode_sd_after})
    for filename, fields, rows in (("scales.csv", _SCALE_COLUMNS, scales),
                                   ("marker_qc.csv", _QC_COLUMNS, qc)):
        with task.output(filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    _figure(task, qc, len(entries))
    shifted = {row["protein_entry"] for row in scales if row["marker"] == "CD8"
               and row["scale"] != ""}
    return {"markers_measured": len(markers), "markers_assessed": len(qc),
            "markers_fitted": sum(row["reference"] != "" for row in qc),
            "sections": len(entries), "CD8_shifted_sections": len(shifted),
            "over_bar": [row["marker"] for row in qc if row["reference"]
                         and float(row["mode_sd_after"]) > 0.10]}


def _figure(task: Task, rows: list[dict], sections: int) -> None:
    import matplotlib.pyplot as plt

    labels = [row["marker"] for row in rows]
    before = [float(row["mode_sd_before"]) if row["reference"] else np.nan for row in rows]
    after = [float(row["mode_sd_after"]) if row["reference"] else np.nan for row in rows]
    fig, ax = plt.subplots(figsize=(11, max(5, len(labels) * 0.22 + 1)), layout="constrained")
    y = np.arange(len(labels))
    ax.scatter(before, y, s=13, label="Before")
    ax.scatter(after, y, s=13, label="After")
    ax.axvline(0.10, color="crimson", linestyle="--", label="review bar 0.10")
    ax.set(yticks=y, yticklabels=labels,
           xlabel="Background-mode spread across sections (SD, log intensity)",
           title=f"UniFORM shifts: {len(labels)} measured markers across {sections} sections")
    ax.invert_yaxis()
    ax.legend(fontsize=8)
    fig.savefig(gates.figure_path(task, "mode_spread.png"), dpi=140)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    route = _route(task.project)
    if route is not None:
        return route.gate(task, summary)
    checks = [
        gates.Check("all measured markers assessed", summary["markers_measured"],
                    str(summary["markers_assessed"]),
                    summary["markers_measured"] == summary["markers_assessed"]),
        gates.Check("markers with a usable shift", summary["markers_fitted"],
                    f"of {summary['markers_measured']} measured markers; missing stays unavailable",
                    True, blocking=False),
        gates.Check("CD8 has a UniFORM shift for every section",
                    summary["CD8_shifted_sections"], str(summary["sections"]),
                    summary["CD8_shifted_sections"] == summary["sections"]),
        gates.Check("markers over the background-mode review bar", summary["over_bar"],
                    "review their marker figures", True, blocking=False),
    ]
    figures = ["mode_spread.png"] + [f"markers/{row['marker'].replace('/', '_')}.png"
                                   for row in _read_qc(task) if row["reference"]]
    return gates.Gate(STAGE, "UniFORM shifts on raw cell means; no pooled positivity cutoffs "
                      "or tissue calls.", checks, figures)


def _read_qc(task: Task) -> list[dict]:
    with (task.staging / "marker_qc.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))
