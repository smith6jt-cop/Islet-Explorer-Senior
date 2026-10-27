"""Stage 08 pilot: candidate cell links for an accepted same-slide registration.

This joins the separate Stage 06 assay evidence by exact Stage 03/04 cell keys
and Stage 07 geometry. It never turns a candidate link into ground truth or
publishes a final tissue cell class.
"""

from __future__ import annotations

import json

import pyarrow as pa
import pyarrow.parquet as pq

from .. import gates
from ..artifacts import sha256_file
from ..integration.paired_typing import link_candidate_cells
from ..integration.transform import Transform
from ..tasks import Input, Task
from .pilot_inventory import _manifest
from .same_slide_registration import MAX_DIST_UM, _sample, _selected

STAGE = "08_integration"


class SameSlideIntegrationError(ValueError):
    """The accepted registration or assay evidence is not the declared pair."""


def settings(project) -> dict:
    pair = _selected(project)
    return {"protein_entry": pair.protein_entry,
            "xenium_acquisition_id": pair.xenium_acquisition_id,
            "max_dist_um": MAX_DIST_UM, "scope": "candidate_cell_links_only"}


def inputs(project) -> list[Input]:
    pair = _selected(project)
    sample = _sample(project, pair)
    release = project.release_dir
    registration = release / "07_registration"
    return [Input("same_slide_source_manifest", _manifest(project)),
            Input("rna_regions", release / "04_rna_ingest/regions.csv"),
            Input("registration_receipt", registration / "receipt.json"),
            Input("registration_acceptance", registration / "gate/accepted.json"),
            Input("transform_json", registration / "registration/transform.json"),
            Input("transform_displacement",
                  registration / "registration/transform_displacement.npz", required=False),
            Input("protein_cells", release / "03_protein_measurement/cells"
                  / f"{pair.protein_entry}.parquet"),
            Input("rna_cells", release / "04_rna_ingest/cells" / f"{sample}.parquet"),
            Input("protein_calls", release / "06_typing/protein"
                  / f"{pair.protein_entry}.parquet"),
            Input("rna_calls", release / "06_typing/rna" / f"{sample}.parquet")]


def _accepted_transform(task: Task) -> Transform:
    receipt_path = task.input("registration_receipt")
    receipt = json.loads(receipt_path.read_text())
    acceptance = json.loads(task.input("registration_acceptance").read_text())
    if (receipt.get("status") != "complete" or receipt.get("stage") != "07_registration"
            or acceptance.get("receipt_sha256") != sha256_file(receipt_path)):
        raise SameSlideIntegrationError("Stage 07 acceptance does not bind its current receipt")
    bound = {row["name"]: row.get("sha256") for row in receipt.get("inputs", [])}
    if bound.get("same_slide_source_manifest") != sha256_file(
            task.input("same_slide_source_manifest")):
        raise SameSlideIntegrationError("source manifest changed after registration review")
    outputs = {row["path"]: row["sha256"] for row in receipt["outputs"]}
    for name in ("transform_json", "transform_displacement"):
        path = task.input(name)
        if path.is_file() and outputs.get(f"registration/{path.name}") != sha256_file(path):
            raise SameSlideIntegrationError(f"accepted Stage 07 {name} changed")
    return Transform.load(task.input("transform_json").parent)


def run(task: Task) -> dict:
    pair = _selected(task.project)
    transform = _accepted_transform(task)
    protein = pq.read_table(task.input("protein_cells"), columns=[
        "cell_uid", "object_id", "donor_id", "tissue", "acquisition_id",
        "x_um", "y_um", "protein_eligible"]).to_pandas()
    rna = pq.read_table(task.input("rna_cells"), columns=[
        "cell_uid", "donor_id", "tissue", "acquisition_id",
        "x_um", "y_um", "typing_eligible"]).to_pandas()
    protein_calls = pq.read_table(task.input("protein_calls")).to_pandas()
    rna_calls = pq.read_table(task.input("rna_calls")).to_pandas()
    result = link_candidate_cells(pair, protein, protein_calls, rna, rna_calls,
                                  transform, max_dist_um=MAX_DIST_UM)
    pq.write_table(pa.Table.from_pandas(result.links, preserve_index=False),
                   task.output("candidate_links.parquet"), compression="zstd")
    summary = {"donor_id": pair.donor_id, "tissue": pair.tissue,
               "physical_section": pair.protein_entry,
               "xenium_acquisition_id": pair.xenium_acquisition_id,
               "protein_cells": result.protein_cells, "rna_cells": result.rna_cells,
               "protein_eligible": result.protein_eligible,
               "rna_eligible": result.rna_eligible,
               "candidate_links": len(result.links),
               "unmatched_protein_eligible": result.unmatched_protein_eligible,
               "unmatched_rna_eligible": result.unmatched_rna_eligible,
               "pairing_method": "mutual_nearest_centroid_candidate",
               "max_dist_um": MAX_DIST_UM,
               "final_cell_classes": 0}
    task.output("pairing_counts.json").write_text(json.dumps(summary, indent=2) + "\n")
    _figure(task, summary)
    return summary


def _figure(task: Task, summary: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matched = summary["candidate_links"]
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    labels = ["PhenoCycler eligible", "Xenium eligible"]
    totals = [summary["protein_eligible"], summary["rna_eligible"]]
    ax.barh(labels, [matched, matched], color="#377eb8", label="candidate link")
    ax.barh(labels, [total - matched for total in totals], left=[matched, matched],
            color="#999999", label="unmatched")
    ax.set(xlabel="Cells (n)", title=f"{summary['donor_id']} · {summary['physical_section']}\n"
           f"Mutual-nearest candidates within {MAX_DIST_UM:g} µm; no final type")
    ax.legend()
    fig.savefig(gates.figure_path(task, "candidate_links.png"), dpi=150)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [
        gates.Check("at least one spatial candidate", summary["candidate_links"],
                    "nonzero after reviewed registration",
                    summary["candidate_links"] > 0),
        gates.Check("protein candidate accounting",
                    summary["candidate_links"] + summary["unmatched_protein_eligible"],
                    str(summary["protein_eligible"]),
                    summary["candidate_links"] + summary["unmatched_protein_eligible"]
                    == summary["protein_eligible"]),
        gates.Check("RNA candidate accounting",
                    summary["candidate_links"] + summary["unmatched_rna_eligible"],
                    str(summary["rna_eligible"]),
                    summary["candidate_links"] + summary["unmatched_rna_eligible"]
                    == summary["rna_eligible"]),
        gates.Check("final cell classes", summary["final_cell_classes"],
                    "none before typing review", summary["final_cell_classes"] == 0),
    ]
    return gates.Gate(STAGE, "Same-section spatial candidate links with assay evidence; "
                      "not a verified same-cell or final typing result.", checks,
                      ["candidate_links.png"])
