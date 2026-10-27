"""Stage 06 RNA pilot on exact accepted pancreas and pLN Xenium acquisitions.

This route calls the packaged tissue workers. It does not discover h5ads by
donor, use frozen-partition overrides on new clusters, type protein cells, or
claim same-cell links. A recorded damaged region or missing nucleus values are
reported as unavailable for that acquisition without stopping other regions.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .. import gates, profiles, roots
from ..artifacts import sha256_file
from ..cohort import excluded_on
from ..integration import tissues
from ..integration.manifest import parse_bundle_path
from ..tasks import Input, Task
from ..tissues import lymph_node_xenium_worker as pln_worker
from ..tissues import pancreas_xenium_worker as panc_worker
from ..xenium_typing_input import load_accepted_region, write_worker_h5ad

STAGE = "06_typing"
PANELS = roots.REPO / "core/panels/data"
_COLUMNS = ("acquisition_id", "sample_key", "donor_id", "tissue", "roi", "cells",
            "typing_eligible", "called_cells", "broad_cells", "unresolved_cells",
            "uncalled_cells", "unavailable_cells", "status", "reason")
_STATES = ("class", "broad", "unresolved", "uncalled", "unavailable")
_CONTROL_COLUMNS = ("acquisition_id", "tissue", "sex", "cells", "nucleus_excluded",
                    "y_floor_detect_frac_any", "cd3_detect_frac_in_T_NK",
                    "cd3_detect_frac_in_B_Plasma", "cd3_detect_frac_in_Myeloid")
_EVIDENCE = {"cluster_agree_fraction": ("cell_agrees_with_cluster_call", False),
             "naming_gene_fraction": ("n_naming_detected_for_cluster_call", True),
             "cd3_fraction": ("cd3_detected", False),
             "trac_fraction": ("trac_detected", False),
             "cd4_fraction": ("cd4_detected", False),
             "cd8_fraction": ("cd8_detected", False),
             "foxp3_fraction": ("foxp3_detected", False)}
_CLASS_COLUMNS = ("acquisition_id", "tissue", "lineage", "state", "cells", *_EVIDENCE)


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
            "rna_min_median_detected_genes": _min_genes(project),
            "pancreas_naming_policy": "shared_er_exact_v1",
            "frozen_partition_overrides": "not_applied_to_new_clusters",
            "scope": "rna_pilot_review_only"}


def _min_genes(project) -> int:
    value = project.settings.get("rna_min_median_detected_genes", 0)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise IO60RNATypingError("rna_min_median_detected_genes must be a nonnegative integer")
    return value


def _low_signal(row: dict[str, str], minimum: int) -> bool:
    if minimum == 0 or row["condition"]:
        return False
    try:
        median = float(row["median_detected_genes"])
    except (KeyError, TypeError, ValueError) as exc:
        raise IO60RNATypingError(
            f"{row['acquisition_id']}: missing Stage 04 median detected genes") from exc
    if not math.isfinite(median) or median < 0:
        raise IO60RNATypingError(
            f"{row['acquisition_id']}: invalid Stage 04 median detected genes")
    return median < minimum


def inputs(project) -> list[Input]:
    release = project.release_dir
    inventory = release / "00_inventory"
    ingest = release / "04_rna_ingest"
    metadata = project.source("donor_metadata")
    if metadata is None:
        raise IO60RNATypingError("site must name the recorded donor_metadata workbook")
    npod = project.source("npod_cases")
    if npod is None:
        raise IO60RNATypingError("site must name the Stage 00-bound npod_cases workbook")
    found = [Input("inventory_receipt", inventory / "receipt.json"),
             Input("inventory_acquisitions", inventory / "acquisitions.csv"),
             Input("inventory_gate_summary", inventory / "gate/gate.json"),
             Input("inventory_gate_acceptance", inventory / "gate/accepted.json"),
             Input("donor_metadata", Path(metadata)),
             Input("npod_cases", Path(npod)),
             Input("rna_receipt", ingest / "receipt.json"),
             Input("rna_regions", ingest / "regions.csv"),
             Input("rna_gate_summary", ingest / "gate/gate.json"),
             Input("rna_gate_acceptance", ingest / "gate/accepted.json")]
    for row in _regions(project):
        if row["condition"] or _low_signal(row, _min_genes(project)):
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
    if bound.get("npod_cases") != sha256_file(task.input("npod_cases")):
        raise IO60RNATypingError("nPOD cases differ from accepted Stage 00")
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
    cohort = profiles.load_all()[("cohort", task.project.cohort)]
    aliases = cohort.donor_aliases
    statuses = cohort.disease_status["spellings"]
    frame = pd.read_excel(task.input("donor_metadata"), dtype={"donor_id": str})
    if not {"donor_id", "sex", "disease.status"}.issubset(frame.columns):
        raise IO60RNATypingError("donor metadata lacks sex or disease.status")
    if frame.donor_id.duplicated().any():
        raise IO60RNATypingError("donor metadata has repeated donor_id")
    donors: dict[str, tuple[str, str]] = {}
    for _, row in frame.iterrows():
        raw = str(row["donor_id"])
        donor = aliases.get(raw, raw)
        sex = str(row["sex"])
        disease = statuses.get(str(row["disease.status"]), "")
        if donor in donors and donors[donor] != (sex, disease):
            raise IO60RNATypingError(f"{donor}: conflicting donor metadata aliases")
        donors[donor] = (sex, disease)
    cases = pd.read_excel(task.input("npod_cases"), dtype={"case_id": str})
    if not {"case_id", "sex", "donor_type"}.issubset(cases.columns):
        raise IO60RNATypingError("nPOD cases lack case_id, sex or donor_type")
    for _, row in cases.iterrows():
        raw = str(row["case_id"])
        donor = aliases.get(raw, raw)
        if donor in donors:
            continue  # The accepted project workbook is authoritative where it has a row.
        sex = {"Female": "F", "Male": "M"}.get(str(row["sex"]), "")
        disease = statuses.get(str(row["donor_type"]), "")
        if sex and disease:
            donors[donor] = (sex, disease)
    return donors


def _call_state(lineage: str) -> str:
    if lineage.startswith("Unavailable_"):
        return "unavailable"
    if lineage.startswith("Uncalled_") or lineage.endswith("_uncalled"):
        return "uncalled"
    if lineage.startswith("Unresolved_") or lineage.endswith("_unresolved"):
        return "unresolved"
    if lineage == "Non_immune":
        return "broad"
    return "class"


def _worker_result(path: Path, row: dict[str, str], region) -> tuple[dict[str, int], list[dict]]:
    calls = pd.read_parquet(path)
    if (calls.cell_uid.tolist() != region.cell_uid.tolist()
            or calls.cell_id.tolist() != region.cell_id.tolist()
            or not calls.acquisition_id.eq(row["acquisition_id"]).all()):
        raise IO60RNATypingError(f"{row['acquisition_id']}: worker changed cell identity")
    no_nucleus = ~region.eligible
    if not calls.loc[no_nucleus, "lineage"].eq("Unavailable_no_nucleus").all():
        raise IO60RNATypingError(f"{row['acquisition_id']}: a no-nucleus cell got a call")
    states = dict.fromkeys(_STATES, 0)
    classes = []
    for lineage, count in calls.lineage.value_counts().items():
        state = _call_state(lineage)
        states[state] += int(count)
        group = calls.loc[calls.lineage.eq(lineage)]
        evidence = {}
        for name, (column, positive) in _EVIDENCE.items():
            if column not in group:
                evidence[name] = ""
                continue
            values = pd.to_numeric(group[column], errors="coerce").dropna()
            if not len(values):
                evidence[name] = ""
            else:
                evidence[name] = (float((values > 0).mean()) if positive
                                  else float(values.mean()))
        classes.append({"acquisition_id": row["acquisition_id"], "tissue": row["tissue"],
                        "lineage": lineage, "state": state, "cells": int(count), **evidence})
    if sum(states.values()) != region.n_cells:
        raise IO60RNATypingError(f"{row['acquisition_id']}: lineage states lost cells")
    return states, classes


def run(task: Task) -> dict:
    inventory = _inventory(task)
    donors = _donors(task)
    with task.input("rna_regions").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if [row["acquisition_id"] for row in rows] != task.settings["rna_acquisitions"]:
        raise IO60RNATypingError("Stage 04 selection changed after input declaration")
    panc_panels = pln_panels = None
    results = []
    class_rows = []
    controls = []
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
                          called_cells=0, broad_cells=0, unresolved_cells=0,
                          uncalled_cells=0, unavailable_cells=record["cells"])
            results.append(record)
            continue
        if _low_signal(row, task.settings["rna_min_median_detected_genes"]):
            record.update(status="unavailable",
                          reason=(f"low RNA signal: median {row['median_detected_genes']} "
                                  "detected genes/cell below recorded minimum "
                                  f"{task.settings['rna_min_median_detected_genes']}"),
                          called_cells=0, broad_cells=0, unresolved_cells=0,
                          uncalled_cells=0, unavailable_cells=record["cells"])
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
                          called_cells=0, broad_cells=0, unresolved_cells=0,
                          uncalled_cells=0, unavailable_cells=record["cells"])
            results.append(record)
            continue
        sex, disease = donors.get(donor, ("", ""))
        if not sex or sex == "nan" or not disease or disease == "nan":
            record.update(status="unavailable", reason="recorded sex or disease status unavailable",
                          called_cells=0, broad_cells=0, unresolved_cells=0,
                          uncalled_cells=0, unavailable_cells=record["cells"])
            results.append(record)
            continue
        serial, _ = parse_bundle_path(source["output"])
        recorded = source.get("xenium_serial", "")
        if (not serial or acquisition.split(":")[1] != serial
                or (recorded and recorded != serial)):
            raise IO60RNATypingError(f"{acquisition}: Stage 00 serial differs from output")
        h5ad = task.output(f"_worker/{sample}.h5ad")
        batch = source["batch_number"] or source["run"]
        if not batch:
            raise IO60RNATypingError(f"{acquisition}: no accepted batch or run identity")
        write_worker_h5ad(region, h5ad, sex=sex, disease_status=disease, batch=batch)
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
        states, classes = _worker_result(work / f"{sample}.calls.parquet", row, region)
        class_rows.extend(classes)
        record.update(status="candidate", reason="new section-specific clustering and naming",
                      called_cells=states["class"], broad_cells=states["broad"],
                      unresolved_cells=states["unresolved"], uncalled_cells=states["uncalled"],
                      unavailable_cells=states["unavailable"])
        results.append(record)
        if report.get("acquisition_id") != acquisition:
            raise IO60RNATypingError(f"{acquisition}: worker report lacks source identity")
        cd3 = report.get("cd3_gate_control") or {}
        controls.append({"acquisition_id": acquisition, "tissue": row["tissue"],
                         "sex": sex, "cells": region.n_cells,
                         "nucleus_excluded": int((~region.eligible).sum()),
                         "y_floor_detect_frac_any": report.get("y_floor_detect_frac_any", ""),
                         **{key: cd3.get(key, "") for key in _CONTROL_COLUMNS
                            if key.startswith("cd3_")}})
    with task.output("rna_acquisitions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_COLUMNS)
        writer.writeheader()
        writer.writerows(results)
    with task.output("rna_classes.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_CLASS_COLUMNS)
        writer.writeheader()
        writer.writerows(class_rows)
    with task.output("rna_controls.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_CONTROL_COLUMNS)
        writer.writeheader()
        writer.writerows(controls)
    _figure(task, results, class_rows)
    return {"acquisitions": len(results), "candidate": sum(r["status"] == "candidate"
                                                       for r in results),
            "unavailable": [r["acquisition_id"] for r in results if r["status"] == "unavailable"],
            "cells": sum(r["cells"] for r in results if isinstance(r["cells"], int)),
            "called_cells": sum(r["called_cells"] for r in results),
            "broad_cells": sum(r["broad_cells"] for r in results),
            "unresolved_cells": sum(r["unresolved_cells"] for r in results),
            "uncalled_cells": sum(r["uncalled_cells"] for r in results),
            "unavailable_cells": sum(r["unavailable_cells"] for r in results
                                     if isinstance(r["unavailable_cells"], int)),
            "classes": len({(r["tissue"], r["lineage"]) for r in class_rows
                            if r["state"] == "class"}),
            "unknown_cell_count": [r["acquisition_id"] for r in results if r["cells"] == ""],
            "low_signal_regions": [r["acquisition_id"] for r in results
                                   if str(r["reason"]).startswith("low RNA signal:")],
            "protein_cells_typed": 0, "paired_cell_links": 0,
            "scope": "selected_xenium_acquisitions_only"}


def _figure(task: Task, rows: list[dict], class_rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    measured = [row for row in rows if isinstance(row["cells"], int)]
    names = [f"{row['sample_key']} ({'pLN' if row['tissue'] == tissues.LYMPH_NODE else 'panc'})"
             for row in measured]
    fig, ax = plt.subplots(figsize=(13, max(5, 1.8 + 0.34 * len(measured))),
                           layout="constrained")
    y = list(range(len(measured)))
    left = [0] * len(measured)
    colors = {"called_cells": "#21867a", "broad_cells": "#78a6a0",
              "unresolved_cells": "#e0a047", "uncalled_cells": "#bd784b",
              "unavailable_cells": "#9b72ad"}
    for key, label in (("called_cells", "specific lineage"), ("broad_cells", "broad state"),
                       ("unresolved_cells", "unresolved subset"),
                       ("uncalled_cells", "uncalled"), ("unavailable_cells", "unavailable")):
        values = [row[key] for row in measured]
        ax.barh(y, values, left=left, color=colors[key], label=label)
        left = [a + b for a, b in zip(left, values, strict=True)]
    ax.set(yticks=y, yticklabels=names, xlabel="Segmented Xenium cells (n)",
           title="RNA typing outcome by exact acquisition")
    ax.invert_yaxis()
    unknown = [row["acquisition_id"] for row in rows if row["cells"] == ""]
    if unknown:
        ax.text(0.01, 0.01, f"{len(unknown)} damaged source(s): cell counts unknown",
                transform=ax.transAxes, va="bottom")
    ax.legend(loc="lower right", fontsize=8)
    fig.savefig(gates.figure_path(task, "rna_typing_scope.png"), dpi=160)
    plt.close(fig)

    present = sorted({row["tissue"] for row in class_rows})
    n_labels = len({(row["tissue"], row["lineage"]) for row in class_rows})
    fig, axes = plt.subplots(1, max(1, len(present)),
                             figsize=(14, max(6, 2 + 0.23 * n_labels)),
                             layout="constrained")
    axes = [axes] if len(present) <= 1 else axes
    for ax, tissue in zip(axes, present, strict=False):
        counts: dict[str, int] = {}
        states: dict[str, str] = {}
        for row in class_rows:
            if row["tissue"] == tissue:
                counts[row["lineage"]] = counts.get(row["lineage"], 0) + row["cells"]
                states[row["lineage"]] = row["state"]
        ordered = sorted(counts, key=lambda label: counts[label])
        ax.barh(ordered, [counts[label] for label in ordered],
                color=[{"class": "#21867a", "broad": "#78a6a0", "unresolved": "#e0a047",
                        "uncalled": "#bd784b", "unavailable": "#9b72ad"}[states[label]]
                       for label in ordered])
        ax.set(title=tissue.replace("_", " "), xlabel="Cells with lineage label (n)")
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(axis="x", color="#dddddd", linewidth=0.5)
    if not present:
        axes[0].text(0.5, 0.5, "No measured cells to type", ha="center", va="center")
    fig.suptitle("RNA lineage distribution; labels are candidates from section-specific evidence")
    fig.savefig(gates.figure_path(task, "rna_typing_classes.png"), dpi=160)
    plt.close(fig)
    _support_figure(task, class_rows)


def _support_figure(task: Task, class_rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups = ((tissues.PANCREAS,
               ("cluster_agree_fraction", "naming_gene_fraction")),
              (tissues.LYMPH_NODE,
               ("cluster_agree_fraction", "cd3_fraction", "trac_fraction",
                "cd4_fraction", "cd8_fraction", "foxp3_fraction")))
    sizes = [len({row["lineage"] for row in class_rows if row["tissue"] == tissue
                  and row["state"] in ("class", "broad", "unresolved")})
             for tissue, _ in groups]
    fig, axes = plt.subplots(1, 2, figsize=(15, max(7, 2.5 + 0.38 * max(sizes))),
                             layout="constrained")
    cmap = plt.colormaps["viridis"].with_extremes(bad="#e5e5e5")
    for ax, (tissue, metrics) in zip(axes, groups, strict=True):
        chosen = [row for row in class_rows if row["tissue"] == tissue
                  and row["state"] in ("class", "broad", "unresolved")]
        totals = {label: sum(row["cells"] for row in chosen if row["lineage"] == label)
                  for label in {row["lineage"] for row in chosen}}
        labels = sorted(totals, key=lambda label: -totals[label])
        matrix = np.full((max(1, len(labels)), len(metrics)), np.nan)
        for i, label in enumerate(labels):
            matching = [row for row in chosen if row["lineage"] == label]
            for j, metric in enumerate(metrics):
                measured = [row for row in matching if row[metric] != ""]
                if measured:
                    denominator = sum(row["cells"] for row in measured)
                    matrix[i, j] = sum(row[metric] * row["cells"] for row in measured) / denominator
        image = ax.imshow(np.ma.masked_invalid(matrix), aspect="auto", vmin=0, vmax=1,
                          cmap=cmap)
        ax.set(xticks=range(len(metrics)),
               xticklabels=[metric.replace("_fraction", "").replace("_", " ")
                            for metric in metrics],
               yticks=range(len(labels)),
               yticklabels=[f"{label} (n={totals[label]:,})" for label in labels],
               title=tissue.replace("_", " "))
        ax.tick_params(axis="y", labelsize=8)
        plt.setp(ax.get_xticklabels(), rotation=35, ha="right")
        for i in range(len(labels)):
            for j in range(len(metrics)):
                if np.isfinite(matrix[i, j]):
                    value = matrix[i, j]
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center",
                            color="white" if value < 0.55 else "black", fontsize=7)
        fig.colorbar(image, ax=ax, shrink=0.6, label="Fraction of cells")
    fig.suptitle("Internal RNA class support; gray = not measured; "
                 "these fractions are not ground-truth accuracy")
    fig.savefig(gates.figure_path(task, "rna_typing_support.png"), dpi=160)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [gates.Check("selected acquisitions accounted", summary["acquisitions"],
                          str(len(task.settings["rna_acquisitions"])),
                          summary["acquisitions"] == len(task.settings["rna_acquisitions"])),
              gates.Check("cell accounting", summary["cells"],
                          "specific plus broad plus unresolved plus uncalled plus unavailable",
                          summary["cells"] == summary["called_cells"]
                          + summary["broad_cells"] + summary["unresolved_cells"]
                          + summary["uncalled_cells"] + summary["unavailable_cells"]),
              gates.Check("unavailable acquisitions", summary["unavailable"],
                          "reported; one section does not suppress another", True,
                          blocking=False),
              gates.Check("unknown cell counts", summary["unknown_cell_count"],
                          "reported for damaged sources without raw cells", True,
                          blocking=False),
              gates.Check("low RNA signal", summary["low_signal_regions"],
                          "reported; no cell classes issued below the declared median-gene floor",
                          True, blocking=False),
              gates.Check("specific lineage classes", summary["classes"],
                          "reported from exact per-cell worker outputs", True,
                          blocking=False),
              gates.Check("paired cell links", summary["paired_cell_links"], "none before Stage 08",
                          summary["paired_cell_links"] == 0)]
    return gates.Gate(STAGE, "New-cluster pancreas/pLN Xenium typing candidate by exact "
                      "accepted acquisition. No protein calls, same-cell links, verified accuracy "
                      "or active-release claim.", checks,
                      ["rna_typing_scope.png", "rna_typing_classes.png",
                       "rna_typing_support.png"])
