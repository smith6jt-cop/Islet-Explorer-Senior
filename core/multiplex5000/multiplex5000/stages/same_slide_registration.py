"""Stage 07 pilot: admit and inspect a section-specific same-slide transform.

The transform is an explicit site input, not inferred from donor names. This
stage binds it to the prior assay receipts, exact cells and source manifest,
then plots the two coordinate clouds for manual registration review.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from .. import gates
from ..integration.paired_typing import validate_pair_cells, validate_pair_transform
from ..integration.sameslide import match_by_centroid
from ..integration.transform import Transform
from ..same_slide_sources import load_pairs
from ..tasks import Input, Task
from .pilot_inventory import _manifest

STAGE = "07_registration"
MAX_DIST_UM = 5.0


class SameSlideRegistrationError(ValueError):
    """This pilot has no exact, reviewable registration source."""


def _selected(project):
    chosen = project.settings.get("rna_acquisitions")
    if not isinstance(chosen, list) or len(chosen) != 1:
        raise SameSlideRegistrationError("the first registration pilot selects one Xenium region")
    matches = [pair for pair in load_pairs(_manifest(project))
               if pair.xenium_acquisition_id == chosen[0]]
    if len(matches) != 1:
        raise SameSlideRegistrationError("Xenium region has no unique recorded source pair")
    return matches[0]


def _sample(project, pair) -> str:
    with (project.release_dir / "04_rna_ingest/regions.csv").open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle)
                if row["acquisition_id"] == pair.xenium_acquisition_id]
    if (len(rows) != 1 or rows[0]["donor_id"] != pair.donor_id
            or rows[0]["tissue"] != pair.tissue
            or Path(rows[0]["source_path"]) != pair.xenium_output):
        raise SameSlideRegistrationError("Stage 04 region differs from source pair")
    return rows[0]["sample_key"]


def _transform_dir(project, pair) -> Path:
    root = project.source("same_slide_registrations")
    if root is None:
        raise SameSlideRegistrationError(
            "site sources.same_slide_registrations must name a reviewed transform root")
    return Path(root) / pair.protein_entry


def settings(project) -> dict:
    pair = _selected(project)
    return {"protein_entry": pair.protein_entry,
            "xenium_acquisition_id": pair.xenium_acquisition_id,
            "max_dist_um": MAX_DIST_UM, "scope": "registration_review_candidate"}


def inputs(project) -> list[Input]:
    pair = _selected(project)
    sample = _sample(project, pair)
    folder = _transform_dir(project, pair)
    release = project.release_dir
    return [Input("same_slide_source_manifest", _manifest(project)),
            Input("rna_regions", release / "04_rna_ingest/regions.csv"),
            Input("typing_receipt", release / "06_typing/receipt.json"),
            Input("typing_acceptance", release / "06_typing/gate/accepted.json"),
            Input("protein_cells", release / "03_protein_measurement/cells"
                  / f"{pair.protein_entry}.parquet"),
            Input("rna_cells", release / "04_rna_ingest/cells" / f"{sample}.parquet"),
            Input("transform_json", folder / "transform.json"),
            Input("transform_displacement", folder / "transform_displacement.npz",
                  required=False)]


def run(task: Task) -> dict:
    pair = _selected(task.project)
    if task.input("transform_json").parent != _transform_dir(task.project, pair):
        raise SameSlideRegistrationError("transform source changed after input declaration")
    receipt = json.loads(task.input("typing_receipt").read_text())
    if receipt.get("stage") != "06_typing" or receipt.get("status") != "complete":
        raise SameSlideRegistrationError("Stage 06 has no complete typing receipt")
    transform = Transform.load(task.input("transform_json").parent)
    validate_pair_transform(pair, transform)
    protein = pq.read_table(task.input("protein_cells"), columns=[
        "cell_uid", "donor_id", "tissue", "acquisition_id", "x_um", "y_um",
        "protein_eligible"]).to_pandas()
    rna = pq.read_table(task.input("rna_cells"), columns=[
        "cell_uid", "donor_id", "tissue", "acquisition_id", "x_um", "y_um",
        "typing_eligible"]).to_pandas()
    protein, rna = validate_pair_cells(pair, protein, rna)
    p = protein.loc[protein["protein_eligible"]].reset_index(drop=True)
    x = rna.loc[rna["typing_eligible"]].reset_index(drop=True)
    fixed = p[["x_um", "y_um"]].to_numpy(dtype=float)
    moving = transform.apply(x[["x_um", "y_um"]].to_numpy(dtype=float))
    if not np.isfinite(moving).all():
        raise SameSlideRegistrationError("transform produced non-finite coordinates")
    matches = match_by_centroid(moving, fixed, max_dist_um=MAX_DIST_UM)
    shutil.copyfile(task.input("transform_json"), task.output("registration/transform.json"))
    displacement = task.input("transform_displacement")
    if displacement.is_file():
        shutil.copyfile(displacement, task.output("registration/transform_displacement.npz"))
    _figure(task, pair, fixed, moving, len(matches))
    return {"donor_id": pair.donor_id, "tissue": pair.tissue,
            "protein_entry": pair.protein_entry,
            "xenium_acquisition_id": pair.xenium_acquisition_id,
            "protein_cells": len(protein), "rna_cells": len(rna),
            "protein_eligible": len(p), "rna_eligible": len(x),
            "mutual_nearest_candidates": len(matches),
            "median_candidate_distance_um": (float(matches["distance_um"].median())
                                             if len(matches) else None),
            "transform_scale": transform.scale,
            "transform_max_displacement_um": transform.max_displacement_um,
            "candidate_only": True}


def _figure(task: Task, pair, fixed: np.ndarray, moving: np.ndarray, matches: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 8), layout="constrained")
    for points, label, color in ((fixed, "PhenoCycler eligible", "#df8f44"),
                                 (moving, "Xenium eligible, transformed", "#377eb8")):
        stride = max(1, len(points) // 3000)
        ax.scatter(points[::stride, 0], points[::stride, 1], s=2, alpha=0.35,
                   label=f"{label} (displayed {len(points[::stride]):,}/{len(points):,})",
                   color=color)
    ax.set(xlabel="Protein crop-local x (µm)", ylabel="Protein crop-local y (µm)",
           title=f"{pair.donor_id} · {pair.protein_entry}\n"
                 f"Candidate registration: {matches:,} mutual-nearest cells within "
                 f"{MAX_DIST_UM:g} µm")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(fontsize=8, markerscale=3)
    fig.savefig(gates.figure_path(task, "registration_overlay.png"), dpi=150)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [
        gates.Check("transform and exact section identities", summary["protein_entry"],
                    task.settings["protein_entry"],
                    summary["protein_entry"] == task.settings["protein_entry"]),
        gates.Check("finite candidate registration geometry",
                    summary["mutual_nearest_candidates"], "at least one candidate",
                    summary["mutual_nearest_candidates"] > 0),
        gates.Check("transform scale in micron frames", summary["transform_scale"],
                    "review overlay and calibration", True, blocking=False),
        gates.Check("non-rigid displacement", summary["transform_max_displacement_um"],
                    "review residuals", True, blocking=False),
    ]
    return gates.Gate(STAGE, "Candidate same-slide registration; inspect overlay before "
                      "accepting. No cell identity or type is verified by this gate.",
                      checks, ["registration_overlay.png"])
