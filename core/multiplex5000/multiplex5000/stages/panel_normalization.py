"""Normalization fitted per panel, from the accepted Protein Measurement expression.

Selected by ``settings.normalization: per_panel`` (the spleen primary cohort). Each panel's
sections enter UniFORM together and never with another panel's: a 33-plex and a 67-plex section
share only some markers, and those they share were acquired on different runs at different scales.
UniFORM's automatic registration picks each marker's reference section; the choice is recorded in
``scales.csv`` and can be pinned for a rerun with the ``normalization_references`` source (a
previous ``scales.csv``).

A marker with fewer than two measurable sections in a panel, or whose fit UniFORM refuses, is
``unavailable`` there: its aligned values are NaN, never the raw value. Every cell keeps its row;
only ``qc_analysis_eligible`` cells enter the fit.

Cell Typing does not read these values. The spleen typing route reads each marker's REDSEA value
from Protein Measurement against that section's own level, so this stage is a prerequisite for
typing (its gate is the alignment review), not its input. The aligned values serve region
definitions, classifier rescaling and display.
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from .. import gates, normalize, profiles, roots
from ..tasks import Input, Task

STAGE = "05_normalization"
SETTING = "per_panel"
MODE_SD_BAR = 0.10
#: Protein Measurement expression columns that are not markers; every other column is one.
META = ("donor_id", "object_id", "qc_analysis_eligible", "qc_estimation_eligible", "image",
        "cell_region", "islet_num")
_SCALE_COLUMNS = ("donor_id", "panel_id", "marker", "reference_donor", "scale", "shift_log",
                  "n_fit_cells", "status", "mode_sd_before", "mode_sd_after")


class PanelNormalizationError(ValueError):
    """The per-panel normalization inputs cannot be combined as declared."""


def applies(project) -> bool:
    return project.settings.get("normalization") == SETTING


def _spec(project) -> Path:
    named = project.source("cohort_spec")
    if named:
        return Path(named)
    return profiles.REPO / profiles.load_all()[("cohort", project.cohort)].membership


def _panels(project) -> dict[str, list[str]]:
    """{panel_id: [donor_id, ...]} in the cohort spec's order."""
    spec = json.loads(_spec(project).read_text())
    default = spec.get("defaults", {}).get("panel_id")
    panels: dict[str, list[str]] = {}
    seen = set()
    for image in spec["images"]:
        donor = str(image["donor_id"])
        panel = image.get("panel_id", default)
        if not panel:
            raise PanelNormalizationError(f"{donor}: the cohort spec names no panel")
        if donor in seen:
            raise PanelNormalizationError(f"{donor}: more than one image in the cohort spec")
        seen.add(donor)
        panels.setdefault(str(panel), []).append(donor)
    return panels


def _uniform(project) -> Path:
    return Path(project.source("uniform_checkout") or roots.REPO / "external/UniFORM")


def _expression(project, donor: str) -> Path:
    return (project.release_dir / "03_protein_measurement/data/expression"
            / f"donor_id={donor}" / "data.parquet")


def settings(project) -> dict:
    return {"panels": _panels(project), "uniform_commit": normalize.UNIFORM_COMMIT,
            "uniform_checkout": str(_uniform(project)), "mode_sd_bar": MODE_SD_BAR,
            "num_bins": normalize.NUM_BINS,
            "pinned_references": bool(project.source("normalization_references"))}


def inputs(project) -> list[Input]:
    stage = project.release_dir / "03_protein_measurement"
    found = [Input("protein_receipt", stage / "receipt.json"),
             Input("protein_acceptance", stage / "gate/accepted.json"),
             Input("cohort_spec", _spec(project))]
    for donors in _panels(project).values():
        found.extend(Input(f"expression:{d}", _expression(project, d)) for d in donors)
    pinned = project.source("normalization_references")
    if pinned:
        found.append(Input("normalization_references", Path(pinned)))
    package = _uniform(project) / "UniFORM"
    if not package.is_dir():
        raise PanelNormalizationError(f"pinned UniFORM package missing at {package}")
    found.extend(Input(f"uniform_code:{path.name}", path) for path in sorted(package.glob("*.py")))
    return found


def _pinned(path: Path | None) -> dict[tuple[str, str], str]:
    """(panel, marker) -> reference donor, from a previous run's scales.csv."""
    if path is None:
        return {}
    with path.open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r["status"] == "aligned"]
    out: dict[tuple[str, str], str] = {}
    for row in rows:
        key = (row["panel_id"], row["marker"])
        if out.setdefault(key, row["reference_donor"]) != row["reference_donor"]:
            raise PanelNormalizationError(f"{key}: the pinned record names two references")
    return out


def _markers(paths: dict[str, Path]) -> list[str]:
    schemas = {d: pq.ParquetFile(p).schema_arrow.names for d, p in paths.items()}
    first = next(iter(schemas.values()))
    if any(names != first for names in schemas.values()):
        raise PanelNormalizationError("expression partitions differ in their columns")
    missing = [c for c in ("donor_id", "object_id", "qc_analysis_eligible") if c not in first]
    if missing:
        raise PanelNormalizationError(f"expression is missing {missing}")
    return [c for c in first if c not in META]


def _sample(path: Path, eligible: np.ndarray, marker: str) -> np.ndarray | None:
    """The eligible cells' values, or None when nothing is measurable (all NaN or below 1)."""
    column = pq.ParquetFile(path).read(columns=[marker]).column(0).to_numpy(zero_copy_only=False)
    values = np.asarray(column, dtype="float64")[eligible]
    return values if np.any(np.isfinite(values) & (values >= 1.0)) else None


def _fit_panel(task: Task, panel: str, donors: list[str], markers: list[str],
               pinned: dict, registration, normalization) -> list[dict]:
    import matplotlib.pyplot as plt

    paths = {d: task.input(f"expression:{d}") for d in donors}
    eligible = {d: pq.ParquetFile(p).read(columns=["qc_analysis_eligible"]).column(0)
                .to_numpy(zero_copy_only=False).astype(bool) for d, p in paths.items()}
    records = []
    for marker in markers:
        sample = {}
        for d in donors:
            values = _sample(paths[d], eligible[d], marker)
            if values is not None:
                sample[d] = values
        result, refusal = None, ""
        if len(sample) >= 2:
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    result = normalize.fit_marker(marker, sample, registration, normalization,
                                                  reference=pinned.get((panel, marker)))
            except normalize.NormalizeError as err:
                refusal = str(err)
            figure = plt.gcf()
            if result is not None and figure.axes:
                figure.savefig(task.output(f"gate/markers/{panel}__{marker}.png"),
                               dpi=110, bbox_inches="tight")
            plt.close("all")
        for d in donors:
            fitted = result is not None and d in sample
            records.append({
                "donor_id": d, "panel_id": panel, "marker": marker,
                "reference_donor": result.reference if result is not None else "",
                "scale": result.scale[d] if fitted else None,
                "shift_log": result.shift_log[d] if fitted else None,
                "n_fit_cells": (result.n_cells.get(d, 0) if result is not None else
                                int(np.count_nonzero(np.isfinite(sample[d]) & (sample[d] >= 1.0)))
                                if d in sample else 0),
                "status": ("aligned" if fitted else "refused" if refusal and d in sample
                           else "unavailable"),
                "mode_sd_before": result.mode_sd_before if result is not None else None,
                "mode_sd_after": result.mode_sd_after if result is not None else None,
            })
    return records


def _write_aligned(task: Task, panels: dict[str, list[str]], markers: list[str],
                   records: list[dict]) -> None:
    scale = {(r["donor_id"], r["marker"]): r["scale"] for r in records}
    for donors in panels.values():
        for d in donors:
            frame = pq.ParquetFile(task.input(f"expression:{d}")).read().to_pandas()
            for marker in markers:
                factor = scale[(d, marker)]
                frame[marker] = (frame[marker].to_numpy("float64") * factor
                                 if factor is not None else np.nan)
            frame.to_parquet(task.output(f"data/aligned/donor_id={d}/data.parquet"),
                             index=False, compression="zstd")


def run(task: Task) -> dict:
    import matplotlib

    matplotlib.use("Agg")

    panels = task.settings["panels"]
    donors = [d for ds in panels.values() for d in ds]
    receipt = json.loads(task.input("protein_receipt").read_text())
    if (receipt.get("stage") != "03_protein_measurement" or receipt.get("status") != "complete"
            or sorted(receipt.get("summary", {}).get("donors", [])) != sorted(donors)):
        raise PanelNormalizationError("the accepted Protein Measurement roster differs from the "
                                      "cohort spec")
    if _panels_from(task.input("cohort_spec")) != panels:
        raise PanelNormalizationError("the cohort spec changed after its inputs were declared")
    markers = _markers({d: task.input(f"expression:{d}") for d in donors})
    pinned = _pinned(task.input("normalization_references")
                     if "normalization_references" in task.inputs else None)
    registration, normalization = normalize.load_uniform(Path(task.settings["uniform_checkout"]))
    records = []
    for panel, panel_donors in panels.items():
        records += _fit_panel(task, panel, panel_donors, markers, pinned,
                              registration, normalization)
    with task.output("scales.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_SCALE_COLUMNS)
        writer.writeheader()
        writer.writerows(records)
    _write_aligned(task, panels, markers, records)
    _figure(task, records, panels)
    return _summary(records, panels, markers)


def _panels_from(spec: Path) -> dict[str, list[str]]:
    from types import SimpleNamespace

    return _panels(SimpleNamespace(source=lambda name: spec if name == "cohort_spec" else None,
                                   cohort=None))


def _summary(records: list[dict], panels: dict[str, list[str]], markers: list[str]) -> dict:
    out = {"sections": sum(len(d) for d in panels.values()), "markers": len(markers),
           "panels": {}}
    for panel, donors in panels.items():
        rows = [r for r in records if r["panel_id"] == panel]
        by_marker: dict[str, list[dict]] = {}
        for r in rows:
            by_marker.setdefault(r["marker"], []).append(r)
        aligned = sorted(m for m, rs in by_marker.items() if any(r["status"] == "aligned" for r in rs))
        refused = sorted(m for m, rs in by_marker.items() if any(r["status"] == "refused" for r in rs))
        over = sorted(m for m, rs in by_marker.items()
                      if rs[0]["mode_sd_after"] is not None and rs[0]["mode_sd_after"] > MODE_SD_BAR)
        cd8 = sum(r["marker"] == "CD8" and r["status"] == "aligned" for r in rows)
        out["panels"][panel] = {
            "sections": len(donors), "markers_assessed": len(by_marker),
            "markers_aligned": aligned, "markers_refused": refused,
            "markers_unavailable": sorted(set(by_marker) - set(aligned) - set(refused)),
            "markers_over_bar": over, "CD8_aligned_sections": cd8}
    return out


def _figure(task: Task, records: list[dict], panels: dict[str, list[str]]) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(panels), figsize=(6 * len(panels), 14), squeeze=False,
                             layout="constrained")
    for ax, (panel, donors) in zip(axes[0], panels.items()):
        first = {}
        for r in records:
            if r["panel_id"] == panel and r["mode_sd_before"] is not None:
                first.setdefault(r["marker"], r)
        labels = sorted(first)
        y = np.arange(len(labels))
        ax.scatter([first[m]["mode_sd_before"] for m in labels], y, s=13, label="Before")
        ax.scatter([first[m]["mode_sd_after"] for m in labels], y, s=13, label="After")
        ax.axvline(MODE_SD_BAR, color="crimson", linestyle="--", label=f"review bar {MODE_SD_BAR}")
        ax.set(yticks=y, yticklabels=labels, title=f"{panel}: {len(labels)} markers, "
               f"{len(donors)} sections",
               xlabel="Background-mode spread across sections (SD, log intensity)")
        ax.tick_params(axis="y", labelsize=7)
        ax.invert_yaxis()
        ax.legend(fontsize=8)
    fig.savefig(gates.figure_path(task, "mode_spread.png"), dpi=130)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    panels = summary["panels"]
    checks = [
        gates.Check("every marker assessed in every panel",
                    {p: v["markers_assessed"] for p, v in panels.items()},
                    str(summary["markers"]),
                    all(v["markers_assessed"] == summary["markers"] for v in panels.values())),
        gates.Check("CD8 aligned on every section",
                    {p: v["CD8_aligned_sections"] for p, v in panels.items()},
                    "every section of every panel",
                    all(v["CD8_aligned_sections"] == v["sections"] for v in panels.values())),
        gates.Check("markers aligned, by panel",
                    {p: len(v["markers_aligned"]) for p, v in panels.items()},
                    "reported", True, blocking=False),
        gates.Check("markers unavailable (NaN in aligned cells), by panel",
                    {p: v["markers_unavailable"] for p, v in panels.items()},
                    "fewer than two measurable sections", True, blocking=False),
        gates.Check("markers whose fit UniFORM refused, by panel",
                    {p: v["markers_refused"] for p, v in panels.items()},
                    "unavailable; review", True, blocking=False),
        gates.Check("markers over the background-mode review bar after alignment",
                    {p: v["markers_over_bar"] for p, v in panels.items()},
                    f"mode SD > {MODE_SD_BAR}; review their figures", True, blocking=False),
    ]
    figures = ["mode_spread.png"] + sorted(
        f"markers/{p.name}" for p in (task.staging / gates.GATE_DIR / "markers").glob("*.png"))
    return gates.Gate(STAGE, "UniFORM shifts fitted per panel from Protein Measurement "
                      "expression. Cell Typing reads REDSEA values against each section's own "
                      "level, not these aligned values.", checks, figures)
