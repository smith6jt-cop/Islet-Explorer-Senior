"""Stage 08 pilot: candidate cell links for an accepted same-slide registration.

This joins the separate Stage 06 assay evidence by exact Stage 03/04 cell keys
and Stage 07 geometry. It never turns a candidate link into ground truth or
publishes a final tissue cell class.
"""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .. import gates
from ..artifacts import sha256_file
from ..integration import paired_typing
from ..integration.paired_typing import link_candidate_cells
from ..integration.transform import Transform
from ..resolved_parameters import write_resolved_parameters
from ..tasks import Input, Task
from .pilot_inventory import _manifest
from .same_slide_registration import (
    MAX_DIST_UM,
    SameSlideRegistrationError,
    _cell_table,
    _receipt_hashes,
    _sample,
    _selected,
    _typed_cells_current,
)

STAGE = "08_integration"


class SameSlideIntegrationError(ValueError):
    """The accepted registration or assay evidence is not the declared pair."""


def settings(project) -> dict:
    pair = _selected(project)
    return {"protein_entry": pair.protein_entry,
            "xenium_acquisition_id": pair.xenium_acquisition_id,
            "max_dist_um": MAX_DIST_UM, "scope": "candidate_cell_links_only",
            "hierarchy_map": project.settings.get("same_slide_hierarchy_map")}


def inputs(project) -> list[Input]:
    pair = _selected(project)
    sample = _sample(project, pair)
    release = project.release_dir
    registration = release / "07_registration"
    return [Input("same_slide_source_manifest", _manifest(project)),
            Input("rna_regions", release / "04_rna_ingest/regions.csv"),
            Input("registration_receipt", registration / "receipt.json"),
            Input("registration_acceptance", registration / "gate/accepted.json"),
            Input("typing_receipt", release / "06_typing/receipt.json"),
            Input("typing_acceptance", release / "06_typing/gate/accepted.json"),
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
            or acceptance.get("stage") != "07_registration"
            or acceptance.get("receipt_sha256") != sha256_file(receipt_path)):
        raise SameSlideIntegrationError("Stage 07 acceptance does not bind its current receipt")
    try:
        bound = _receipt_hashes(receipt.get("inputs", []), "name")
        outputs = _receipt_hashes(receipt.get("outputs", []), "path")
    except SameSlideRegistrationError as error:
        raise SameSlideIntegrationError(str(error)) from error
    if bound.get("same_slide_source_manifest") != sha256_file(
            task.input("same_slide_source_manifest")):
        raise SameSlideIntegrationError("source manifest changed after registration review")
    for name in ["protein_cells", "rna_cells", "rna_regions", "typing_receipt", "typing_acceptance"]:
        if bound.get(name) != sha256_file(task.input(name)):
            raise SameSlideIntegrationError(f"accepted Stage 07 source {name} changed or is unbound")
    typing_path = task.input("typing_receipt")
    typing = json.loads(typing_path.read_text())
    typing_acceptance = json.loads(task.input("typing_acceptance").read_text())
    if (typing.get("stage") != "06_typing" or typing.get("status") != "complete"
            or typing_acceptance.get("stage") != "06_typing"
            or typing_acceptance.get("receipt_sha256") != sha256_file(typing_path)):
        raise SameSlideIntegrationError("Stage 06 acceptance does not bind the independent calls")
    try:
        _typed_cells_current(task, typing)
        typing_outputs = _receipt_hashes(typing.get("outputs", []), "path")
    except SameSlideRegistrationError as error:
        raise SameSlideIntegrationError(str(error)) from error
    for name, folder in [("protein_calls", "protein"), ("rna_calls", "rna")]:
        path = task.input(name)
        if typing_outputs.get(f"{folder}/{path.name}") != sha256_file(path):
            raise SameSlideIntegrationError(f"accepted Stage 06 {name} changed or is unbound")
    for name in ("transform_json", "transform_displacement"):
        path = task.input(name)
        expected = outputs.get(f"registration/{path.name}")
        if (expected is not None and not path.is_file()) or (path.is_file() and expected != sha256_file(path)):
            raise SameSlideIntegrationError(f"accepted Stage 07 {name} changed")
    return Transform.load(task.input("transform_json").parent)


def run(task: Task) -> dict:
    pair = _selected(task.project)
    transform = _accepted_transform(task)
    protein_columns = [
        "cell_uid", "object_id", "donor_id", "tissue", "acquisition_id",
        "x_um", "y_um", "protein_eligible"]
    rna_columns = [
        "cell_uid", "donor_id", "tissue", "acquisition_id",
        "x_um", "y_um", "typing_eligible"]
    protein = _cell_table(task, "protein_cells", protein_columns)
    rna = _cell_table(task, "rna_cells", rna_columns)
    protein_calls = pq.read_table(task.input("protein_calls")).to_pandas()
    rna_calls = pq.read_table(task.input("rna_calls")).to_pandas()
    result = link_candidate_cells(pair, protein, protein_calls, rna, rna_calls,
                                  transform, max_dist_um=MAX_DIST_UM,
                                  hierarchy_map=task.settings["hierarchy_map"])
    pq.write_table(pa.Table.from_pandas(result.links, preserve_index=False),
                   task.output("candidate_links.parquet"), compression="zstd")
    pq.write_table(pa.Table.from_pandas(result.candidates, preserve_index=False),
                   task.output("spatial_candidates.parquet"), compression="zstd")
    pq.write_table(pa.Table.from_pandas(result.cell_dispositions, preserve_index=False),
                   task.output("cell_link_dispositions.parquet"), compression="zstd")
    task.output("linking_metadata.json").write_text(json.dumps(result.metadata, indent=2) + "\n")
    summary = {"donor_id": pair.donor_id, "tissue": pair.tissue,
               "physical_section": pair.protein_entry,
               "xenium_acquisition_id": pair.xenium_acquisition_id,
               "protein_cells": result.protein_cells, "rna_cells": result.rna_cells,
               "protein_eligible": result.protein_eligible,
               "rna_eligible": result.rna_eligible,
               "candidate_links": len(result.links),
               "unmatched_protein_eligible": result.unmatched_protein_eligible,
               "unmatched_rna_eligible": result.unmatched_rna_eligible,
               "spatial_candidate_edges": len(result.candidates),
               "ambiguous_candidate_edges": int(result.candidates["nominal_state"].eq("ambiguous_nearest").sum()),
               "cell_dispositions": result.cell_dispositions.groupby(["assay", "link_state"]).size().rename("cells").reset_index().to_dict("records"),
               "registration_uncertainty": result.metadata["registration_uncertainty"],
               "hierarchy_comparison": result.metadata["hierarchy_comparison"],
               "pairing_method": result.metadata["method"],
               "max_dist_um": MAX_DIST_UM,
               "final_cell_classes": 0}
    summary["parameter_recording"] = write_resolved_parameters(
        task, operation="same_slide_candidate_matching",
        context={"modalities": ["protein", "RNA"], "tissue": pair.tissue,
                 "relationship_type": pair.relationship, "coordinate_units": "um"},
        parameters=[{"name": "max_dist_um", "role": "fixed", "value": MAX_DIST_UM, "units": "um"},
                    {"name": "matching_method", "role": "selected_method", "value": result.metadata["method"]},
                    {"name": "registration_uncertainty", "role": "fitted",
                     "value": result.metadata["registration_uncertainty"]},
                    {"name": "hierarchy_comparison", "role": "selected_method",
                     "value": result.metadata["hierarchy_comparison"]}],
        acquisitions=[{"protein_entry": pair.protein_entry, "rna_acquisition_id": pair.xenium_acquisition_id}],
        input_names=["registration_receipt", "registration_acceptance", "typing_receipt", "typing_acceptance", "transform_json"],
        producer_paths=[Path(__file__), Path(__file__).with_name("same_slide_registration.py"),
                        Path(paired_typing.__file__), Path(paired_typing.__file__).with_name("transform.py")])
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
