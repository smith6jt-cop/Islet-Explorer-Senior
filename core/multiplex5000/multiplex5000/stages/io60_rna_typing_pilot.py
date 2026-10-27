"""Stage 06 RNA pilot on exact accepted pancreas and pLN Xenium acquisitions.

This route calls the packaged tissue workers. It does not discover h5ads by
donor, use frozen-partition overrides on new clusters, type protein cells, or
claim same-cell links. A recorded damaged region or missing nucleus values are
reported as unavailable for that acquisition without stopping other regions.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd

from .. import gates, roots
from ..artifacts import sha256_file
from ..cohort import excluded_on
from ..integration import tissues
from ..tasks import Input, Task
from ..tissues import lymph_node_xenium_worker as pln_worker
from ..tissues import pancreas_xenium_worker as panc_worker
from ..xenium_typing_input import load_accepted_region, write_worker_h5ad

STAGE = "06_typing"
PANELS = roots.REPO / "core/panels/data"
_COLUMNS = ("acquisition_id", "sample_key", "donor_id", "tissue", "roi", "cells",
            "typing_eligible", "called_cells", "unavailable_cells", "status", "reason")


class IO60RNATypingError(ValueError):
    """A selected RNA acquisition cannot enter the tissue-specific worker."""


def _regions(project) -> list[dict[str, str]]:
    chosen = project.settings.get("rna_acquisitions")
    if (not isinstance(chosen, list) or not chosen
            or not all(isinstance(value, str) and value for value in chosen)
            or len(set(chosen)) != len(chosen)):
        raise IO60RNATypingError("RNA pilot needs distinct selected acquisition IDs")
    path = project.release_dir / "04_rna_ingest/regions.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["acquisition_id"] for row in rows]
    if (set(ids) != set(chosen) or len(ids) != len(chosen)
            or any(row["tissue"] not in (tissues.PANCREAS, tissues.LYMPH_NODE)
                   for row in rows)):
        raise IO60RNATypingError("Stage 04 rows differ from the selected pancreas/pLN pilot")
    return rows


def settings(project) -> dict:
    return {"rna_acquisitions": [row["acquisition_id"] for row in _regions(project)],
            "pancreas_naming_policy": "shared_er_exact_v1",
            "frozen_partition_overrides": "not_applied_to_new_clusters",
            "scope": "rna_pilot_review_only"}


def inputs(project) -> list[Input]:
    release = project.release_dir
    inventory = release / "00_inventory"
    ingest = release / "04_rna_ingest"
    metadata = project.source("donor_metadata")
    if metadata is None:
        raise IO60RNATypingError("site must name the recorded donor_metadata workbook")
    found = [Input("inventory_receipt", inventory / "receipt.json"),
             Input("inventory_acquisitions", inventory / "acquisitions.csv"),
             Input("inventory_gate_summary", inventory / "gate/gate.json"),
             Input("inventory_gate_acceptance", inventory / "gate/accepted.json"),
             Input("donor_metadata", Path(metadata)),
             Input("rna_receipt", ingest / "receipt.json"),
             Input("rna_regions", ingest / "regions.csv"),
             Input("rna_gate_summary", ingest / "gate/gate.json"),
             Input("rna_gate_acceptance", ingest / "gate/accepted.json")]
    for row in _regions(project):
        if row["condition"]:
            # A recorded damaged region is reported, not opened. Its source
            # cannot make another acquisition's Stage 06 input incomplete.
            continue
        sample = row["sample_key"]
        folder = Path(row["source_path"])
        found += [Input(f"rna_cells:{sample}", ingest / "cells" / f"{sample}.parquet"),
                  Input(f"raw_matrix:{sample}", folder / "cell_feature_matrix.h5"),
                  Input(f"raw_cells:{sample}", folder / "cells.parquet")]
    for folder in (PANELS / "identity_core", PANELS / "tissues/pancreas/identity_core"):
        found += [Input(f"panel:{path.relative_to(PANELS)}", path)
                  for path in sorted(folder.glob("*.csv"))]
    return found


def _inventory(task: Task) -> dict[str, dict[str, str]]:
    root = task.input("inventory_receipt").parent
    accepted = gates.require_accepted(root)
    receipt = json.loads(task.input("inventory_receipt").read_text())
    if (receipt.get("stage") != "00_inventory" or receipt.get("status") != "complete"
            or accepted.get("receipt_sha256") != sha256_file(task.input("inventory_receipt"))):
        raise IO60RNATypingError("Stage 00 receipt differs from its acceptance")
    bound = {item["name"]: item.get("sha256") for item in receipt["inputs"]}
    if bound.get("donor_metadata") != sha256_file(task.input("donor_metadata")):
        raise IO60RNATypingError("donor metadata differs from accepted Stage 00")
    output = [item for item in receipt["outputs"] if item["path"] == "acquisitions.csv"]
    if len(output) != 1 or output[0]["sha256"] != sha256_file(
            task.input("inventory_acquisitions")):
        raise IO60RNATypingError("Stage 00 acquisition roster changed")
    with task.input("inventory_acquisitions").open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["platform"] == "xenium"]
    ids = [row["acquisition_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise IO60RNATypingError("repeated Stage 00 Xenium acquisition")
    return dict(zip(ids, rows, strict=True))


def _donors(task: Task) -> dict[str, tuple[str, str]]:
    frame = pd.read_excel(task.input("donor_metadata"), dtype={"donor_id": str})
    if not {"donor_id", "sex", "disease.status"}.issubset(frame.columns):
        raise IO60RNATypingError("donor metadata lacks sex or disease.status")
    if frame.donor_id.duplicated().any():
        raise IO60RNATypingError("donor metadata has repeated donor_id")
    return {str(row["donor_id"]): (str(row["sex"]), str(row["disease.status"]))
            for _, row in frame.iterrows()}


def _worker_result(path: Path, row: dict[str, str], region) -> tuple[int, int]:
    calls = pd.read_parquet(path, columns=["cell_uid", "cell_id", "acquisition_id", "lineage"])
    if (calls.cell_uid.tolist() != region.cell_uid.tolist()
            or calls.cell_id.tolist() != region.cell_id.tolist()
            or not calls.acquisition_id.eq(row["acquisition_id"]).all()):
        raise IO60RNATypingError(f"{row['acquisition_id']}: worker changed cell identity")
    no_nucleus = ~region.eligible
    if not calls.loc[no_nucleus, "lineage"].eq("Unavailable_no_nucleus").all():
        raise IO60RNATypingError(f"{row['acquisition_id']}: a no-nucleus cell got a call")
    unavailable = int(calls.lineage.str.startswith(("Unavailable_", "Uncalled_", "Unresolved_")).sum())
    return len(calls) - unavailable, unavailable


def run(task: Task) -> dict:
    inventory = _inventory(task)
    donors = _donors(task)
    with task.input("rna_regions").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if [row["acquisition_id"] for row in rows] != task.settings["rna_acquisitions"]:
        raise IO60RNATypingError("Stage 04 selection changed after input declaration")
    panc_panels = pln_panels = None
    results = []
    for row in rows:
        acquisition = row["acquisition_id"]
        source = inventory.get(acquisition)
        if (source is None or source["donor_id"] != row["donor_id"]
                or source["tissue"] != row["tissue"] or source["roi"] != row["roi"]):
            raise IO60RNATypingError(f"{acquisition}: Stage 00 identity differs from Stage 04")
        donor, sample = row["donor_id"], row["sample_key"]
        if excluded_on(donor, "xenium"):
            raise IO60RNATypingError(f"{acquisition}: centrally excluded Xenium acquisition")
        record = dict.fromkeys(_COLUMNS, "")
        record.update(acquisition_id=acquisition, sample_key=sample, donor_id=donor,
                      tissue=row["tissue"], roi=row["roi"])
        if row["cells"]:
            record["cells"] = int(row["cells"])
            record["typing_eligible"] = int(row["typing_eligible"])
        if row["condition"]:
            record.update(status="unavailable", reason=f"recorded tissue condition: {row['condition']}",
                          called_cells=0, unavailable_cells=record["cells"])
            results.append(record)
            continue
        if not row["cells"]:
            raise IO60RNATypingError(f"{acquisition}: intact region lacks Stage 04 cell counts")
        region = load_accepted_region(
            row, task.input("rna_receipt"), task.input("rna_regions"),
            task.input(f"rna_cells:{sample}"), task.input(f"raw_matrix:{sample}"),
            task.input(f"raw_cells:{sample}"))
        if region.n_cells != record["cells"] or int(region.eligible.sum()) != record["typing_eligible"]:
            raise IO60RNATypingError(f"{acquisition}: Stage 04 cell counts changed")
        if region.nucleus_missing.any():
            record.update(status="unavailable",
                          reason=f"{int(region.nucleus_missing.sum())} missing nucleus measurements",
                          called_cells=0, unavailable_cells=record["cells"])
            results.append(record)
            continue
        sex, disease = donors.get(donor, ("", ""))
        if not sex or sex == "nan" or not disease or disease == "nan":
            record.update(status="unavailable", reason="recorded sex or disease status unavailable",
                          called_cells=0, unavailable_cells=record["cells"])
            results.append(record)
            continue
        serial = source["xenium_serial"]
        if not serial or serial not in acquisition:
            raise IO60RNATypingError(f"{acquisition}: missing exact Stage 00 serial")
        h5ad = task.output(f"_worker/{sample}.h5ad")
        write_worker_h5ad(region, h5ad, sex=sex, disease_status=disease,
                          batch=source["batch_number"])
        work = task.staging / "rna" / sample
        if row["tissue"] == tissues.PANCREAS:
            if panc_panels is None:
                panc_panels = panc_worker.load_panels()
            report = panc_worker.process_section(
                donor, sample, h5ad, panc_panels, work,
                naming_policy="shared_er_exact_v1", source_serial=serial,
                exact_cells=task.input(f"raw_cells:{sample}"), acquisition_id=acquisition)
        else:
            if pln_panels is None:
                pln_panels = pln_worker.load_panels()
            report = pln_worker.process_exact_section(
                donor, sample, h5ad, task.input(f"raw_cells:{sample}"), pln_panels,
                work, acquisition_id=acquisition, roi=row["roi"])
        h5ad.unlink()
        called, unavailable = _worker_result(work / f"{sample}.calls.parquet", row, region)
        record.update(status="candidate", reason="new section-specific clustering and naming",
                      called_cells=called, unavailable_cells=unavailable)
        results.append(record)
        if report.get("acquisition_id") != acquisition:
            raise IO60RNATypingError(f"{acquisition}: worker report lacks source identity")
    with task.output("rna_acquisitions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_COLUMNS)
        writer.writeheader()
        writer.writerows(results)
    _figure(task, results)
    return {"acquisitions": len(results), "candidate": sum(r["status"] == "candidate"
                                                       for r in results),
            "unavailable": [r["acquisition_id"] for r in results if r["status"] == "unavailable"],
            "cells": sum(r["cells"] for r in results if isinstance(r["cells"], int)),
            "called_cells": sum(r["called_cells"] for r in results),
            "unavailable_cells": sum(r["unavailable_cells"] for r in results
                                     if isinstance(r["unavailable_cells"], int)),
            "unknown_cell_count": [r["acquisition_id"] for r in results if r["cells"] == ""],
            "protein_cells_typed": 0, "paired_cell_links": 0,
            "scope": "selected_xenium_acquisitions_only"}


def _figure(task: Task, rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    measured = [row for row in rows if isinstance(row["cells"], int)]
    names = [f"{row['sample_key']}\n{row['tissue']}" for row in measured]
    called = [row["called_cells"] for row in measured]
    unavailable = [row["unavailable_cells"] for row in measured]
    fig, ax = plt.subplots(figsize=(max(8, len(rows) * 2.5), 5), layout="constrained")
    ax.bar(names, called, label="candidate called")
    ax.bar(names, unavailable, bottom=called, label="unavailable or unresolved")
    ax.set(ylabel="Xenium segmented cells (n)",
           title="Stage 06 tissue-specific RNA pilot; no protein typing or paired-cell claim")
    unknown = [row["acquisition_id"] for row in rows if row["cells"] == ""]
    if unknown:
        ax.text(0.01, 0.98, f"{len(unknown)} damaged source(s): cell counts unknown",
                transform=ax.transAxes, va="top")
    ax.legend()
    fig.savefig(gates.figure_path(task, "rna_typing_scope.png"), dpi=160)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [gates.Check("selected acquisitions accounted", summary["acquisitions"],
                          str(len(task.settings["rna_acquisitions"])),
                          summary["acquisitions"] == len(task.settings["rna_acquisitions"])),
              gates.Check("cell accounting", summary["cells"],
                          "called plus unavailable or unresolved",
                          summary["cells"] == summary["called_cells"]
                          + summary["unavailable_cells"]),
              gates.Check("unavailable acquisitions", summary["unavailable"],
                          "reported; one section does not suppress another", True,
                          blocking=False),
              gates.Check("unknown cell counts", summary["unknown_cell_count"],
                          "reported for damaged sources without raw cells", True,
                          blocking=False),
              gates.Check("paired cell links", summary["paired_cell_links"], "none before Stage 08",
                          summary["paired_cell_links"] == 0)]
    return gates.Gate(STAGE, "New-cluster pancreas/pLN Xenium typing candidate by exact "
                      "accepted acquisition. No protein calls, same-cell links, verified accuracy "
                      "or active-release claim.", checks, ["rna_typing_scope.png"])
