"""Stage 06 review-only protein-rule evidence for dedicated pLN sections.

One selected, full 6551-derived rule source is evaluated by the same section
engine as spleen. Raw QuPath means and UniFORM shifts remain a candidate until
the approved compartment chain and processing review are complete. Rule matches
are not a served final cell class.
"""

from __future__ import annotations

import csv
import json

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .. import gates, roots
from ..tasks import Input, Task
from ..tissues.evaluate import CANDIDATE_STATES, resolve_unique_matched_rules
from ..tissues.lymph_node import load_rule_source
from ..tissues.protein_section import evaluate_protein_section

STAGE = "06_typing"
PANEL = roots.REPO / "core/panels/data/tissues/pancreatic_lymph_node/manifest.yml"
RULES = PANEL.parent / "phenotype_rules.csv"
THRESHOLD_COLUMNS = ("protein_entry", "marker", "floor", "top", "cutoff", "shift_scale",
                     "n_background_support", "unavailable_reason", "detected", "undetected",
                     "unavailable")
RULE_COLUMNS = ("protein_entry", "phenotype", "matched", "not_matched", "unavailable")
CANDIDATE_COLUMNS = ("protein_entry", "status", "cells")


class PLNTypingPilotError(ValueError):
    """A pLN candidate cannot be evaluated from the declared stage outputs."""


def _sections(project) -> list[dict[str, str]]:
    path = project.release_dir / "03_protein_measurement/sections.csv"
    if not path.is_file():
        raise PLNTypingPilotError("dedicated pLN Stage 03 sections are not published")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    entries = [row["protein_entry"] for row in rows]
    if (not rows or len(set(entries)) != len(entries)
            or len(rows) != len(project.settings["pln_export_images"])
            or any(row["tissue"] != "pancreatic_lymph_node" for row in rows)):
        raise PLNTypingPilotError("dedicated pLN Stage 03 section roster is ambiguous")
    return rows


def settings(project) -> dict:
    return {"protein_entries": [row["protein_entry"] for row in _sections(project)],
            "protein_cutoff_fraction": 0.10,
            "scope": "review_candidate_no_final_cell_class"}


def inputs(project) -> list[Input]:
    release = project.release_dir
    found = [Input("protein_sections", release / "03_protein_measurement/sections.csv"),
             Input("protein_receipt", release / "03_protein_measurement/receipt.json"),
             Input("protein_gate", release / "03_protein_measurement/gate/accepted.json"),
             Input("normalization_scales", release / "05_normalization/scales.csv"),
             Input("normalization_receipt", release / "05_normalization/receipt.json"),
             Input("normalization_gate", release / "05_normalization/gate/accepted.json"),
             Input("pln_panel", PANEL), Input("pln_rules", RULES)]
    for row in _sections(project):
        entry = row["protein_entry"]
        found.append(Input(f"protein_cells:{entry}",
                           release / "03_protein_measurement/cells" / f"{entry}.parquet"))
    return found


def _scales(task: Task) -> dict[tuple[str, str], float | None]:
    with task.input("normalization_scales").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    found = {}
    for row in rows:
        key = (row["protein_entry"], row["marker"])
        if key in found:
            raise PLNTypingPilotError(f"duplicate UniFORM scale: {key}")
        found[key] = float(row["scale"]) if row["scale"] else None
    return found


def _complete(task: Task, name: str, stage: str) -> None:
    receipt = json.loads(task.input(name).read_text())
    if receipt.get("status") != "complete" or receipt.get("stage") != stage:
        raise PLNTypingPilotError(f"{stage} has no complete receipt")


def run(task: Task) -> dict:
    _complete(task, "protein_receipt", "03_protein_measurement")
    _complete(task, "normalization_receipt", "05_normalization")
    with task.input("protein_sections").open(newline="") as handle:
        sections = list(csv.DictReader(handle))
    entries = [row["protein_entry"] for row in sections]
    if entries != task.settings["protein_entries"]:
        raise PLNTypingPilotError("pLN sections changed after input declaration")
    source = load_rule_source(task.input("pln_panel"), table_path=task.input("pln_rules"))
    scales = _scales(task)
    threshold_rows, rule_rows, candidate_rows = [], [], []
    total_cells = no_nucleus = 0
    for row in sections:
        entry = row["protein_entry"]
        table = pq.read_table(task.input(f"protein_cells:{entry}"))
        markers = sorted(name.removeprefix("marker:") for name in table.schema.names
                         if name.startswith("marker:"))
        eligible = table["protein_eligible"].to_numpy(zero_copy_only=False).astype(bool)
        if len(eligible) != int(row["cells"]):
            raise PLNTypingPilotError(f"{entry}: Stage 03 cell count changed")
        measurements = {marker: table[f"marker:{marker}"].to_pandas().to_numpy(dtype=float)
                        for marker in markers}
        section = evaluate_protein_section(
            source, entry, measurements, eligible,
            {marker: scales.get((entry, marker)) for marker in markers})
        candidates = resolve_unique_matched_rules(source, section.rules)
        payload = {"cell_uid": table["cell_uid"], "object_id": table["object_id"],
                   "protein_eligible": table["protein_eligible"],
                   "candidate_phenotype": pa.array(candidates.phenotype),
                   "candidate_state": pa.array(candidates.status),
                   "candidate_unavailable_rules": pa.array(candidates.unavailable_rules)}
        for marker, evidence in section.markers.items():
            payload[f"marker_call:{marker}"] = pa.array(evidence.states)
            threshold_rows.append({"protein_entry": entry, "marker": marker,
                                   "floor": evidence.floor, "top": evidence.top,
                                   "cutoff": evidence.cutoff,
                                   "shift_scale": evidence.shift_scale,
                                   "n_background_support": evidence.n_background_support,
                                   "unavailable_reason": evidence.unavailable_reason,
                                   "detected": int(np.count_nonzero(evidence.states == 1)),
                                   "undetected": int(np.count_nonzero(evidence.states == 0)),
                                   "unavailable": int(np.count_nonzero(evidence.states == -1))})
        for phenotype, states in section.rules.states.items():
            payload[f"rule:{phenotype}"] = pa.array(states)
            rule_rows.append({"protein_entry": entry, "phenotype": phenotype,
                              "matched": int(np.count_nonzero(states == 1)),
                              "not_matched": int(np.count_nonzero(states == 0)),
                              "unavailable": int(np.count_nonzero(states == -1))})
        candidate_rows.extend({"protein_entry": entry, "status": status,
                               "cells": int(np.count_nonzero(candidates.status == status))}
                              for status in CANDIDATE_STATES)
        pq.write_table(pa.table(payload), task.output(f"protein/{entry}.parquet"),
                       compression="zstd")
        total_cells += len(eligible)
        no_nucleus += int(np.count_nonzero(~eligible))
    for filename, columns, rows in (("protein_thresholds.csv", THRESHOLD_COLUMNS, threshold_rows),
                                    ("protein_rule_counts.csv", RULE_COLUMNS, rule_rows),
                                    ("protein_candidate_status.csv", CANDIDATE_COLUMNS,
                                     candidate_rows)):
        with task.output(filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
    _marker_figure(task, threshold_rows, sections)
    _rule_figures(task, rule_rows, sections)
    _candidate_figure(task, candidate_rows, sections)
    return {"sections": len(sections), "protein_cells": total_cells,
            "protein_cells_without_nucleus": no_nucleus,
            "protein_markers": len(threshold_rows),
            "protein_markers_with_cutoff": sum(row["cutoff"] is not None
                                               for row in threshold_rows),
            "protein_rules": len(rule_rows),
            "expected_protein_rules": len(sections) * len(source.rules),
            "candidate_status_cells": sum(row["cells"] for row in candidate_rows),
            "unique_rule_candidates": sum(row["cells"] for row in candidate_rows
                                          if row["status"] == "unique_candidate"),
            "candidate_only": True, "final_cell_classes": 0}


def _marker_figure(task: Task, rows: list[dict], sections: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    markers = sorted({row["marker"] for row in rows})
    entries = [row["protein_entry"] for row in sections]
    eligible = {row["protein_entry"]: int(row["eligible"]) for row in sections}
    indexed = {(row["protein_entry"], row["marker"]): row for row in rows}
    values = np.full((len(markers), len(entries)), np.nan)
    for y, marker in enumerate(markers):
        for x, entry in enumerate(entries):
            row = indexed[(entry, marker)]
            if row["cutoff"] is not None and eligible[entry]:
                values[y, x] = 100 * row["detected"] / eligible[entry]
    fig, ax = plt.subplots(figsize=(max(7, len(entries) * 2.2),
                                    max(9, len(markers) * 0.23)), layout="constrained")
    cmap = plt.get_cmap("viridis").with_extremes(bad="#d9d9d9")
    image = ax.imshow(np.ma.masked_invalid(values), aspect="auto", cmap=cmap,
                      vmin=0, vmax=100)
    ax.set(xticks=range(len(entries)),
           xticklabels=[f"{entry}\n{eligible[entry]:,} eligible" for entry in entries],
           yticks=range(len(markers)), yticklabels=markers,
           xlabel="Physical pLN section (n eligible nuclei)", ylabel="Measured protein marker",
           title="Section-local 0.10 candidate positivity; gray = unavailable")
    ax.tick_params(axis="y", labelsize=7)
    fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02,
                 label="Candidate-positive eligible cells (%)")
    fig.savefig(gates.figure_path(task, "pln_marker_support.png"), dpi=140)
    plt.close(fig)


def _rule_figures(task: Task, rows: list[dict], sections: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for section in sections:
        entry = section["protein_entry"]
        grouped = sorted((row for row in rows if row["protein_entry"] == entry),
                         key=lambda row: (-row["matched"], row["phenotype"]))
        names = [row["phenotype"] for row in grouped]
        matched = [row["matched"] for row in grouped]
        unavailable = [row["unavailable"] for row in grouped]
        fig, ax = plt.subplots(figsize=(11, max(10, len(grouped) * 0.29 + 1)),
                               layout="constrained")
        ax.barh(names, matched, label="rule matched")
        ax.barh(names, unavailable, left=matched, label="rule unavailable")
        ax.set(xlabel="Cells in this physical pLN section (n)",
               title=f"{entry}: {int(section['cells']):,} cells; "
                     "hierarchy matches and unavailable evidence")
        ax.invert_yaxis()
        ax.tick_params(axis="y", labelsize=7)
        ax.legend(fontsize=8)
        fig.savefig(gates.figure_path(task, f"pln_rule_evidence_{entry}.png"), dpi=140)
        plt.close(fig)


def _candidate_figure(task: Task, rows: list[dict], sections: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    entries = [row["protein_entry"] for row in sections]
    counts = {(row["protein_entry"], row["status"]): row["cells"] for row in rows}
    fig, ax = plt.subplots(figsize=(max(8, len(entries) * 1.5), 5), layout="constrained")
    bottom = np.zeros(len(entries), dtype=int)
    for status in CANDIDATE_STATES:
        values = np.array([counts[(entry, status)] for entry in entries])
        ax.bar(entries, values, bottom=bottom, label=status.replace("_", " "))
        bottom += values
    ax.set(ylabel="PhenoCycler segmented cells (n)",
           title="pLN hierarchy: deepest unique rule is a review candidate")
    ax.tick_params(axis="x", rotation=60)
    ax.legend()
    fig.savefig(gates.figure_path(task, "pln_candidate_status.png"), dpi=140)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [
        gates.Check("all selected pLN sections evaluated", summary["sections"],
                    str(len(task.settings["protein_entries"])),
                    summary["sections"] == len(task.settings["protein_entries"])),
        gates.Check("full recorded hierarchy evaluated", summary["protein_rules"],
                    str(summary["expected_protein_rules"]),
                    summary["protein_rules"] == summary["expected_protein_rules"]),
        gates.Check("marker cutoffs available", summary["protein_markers_with_cutoff"],
                    f"of {summary['protein_markers']}; unavailable is preserved",
                    True, blocking=False),
        gates.Check("candidate status accounts for every cell", summary["candidate_status_cells"],
                    str(summary["protein_cells"]),
                    summary["candidate_status_cells"] == summary["protein_cells"]),
        gates.Check("unique deepest rule candidates", summary["unique_rule_candidates"],
                    "review only; overlapping matches stay ambiguous", True, blocking=False),
        gates.Check("final cell classes released", 0, "processing approval required",
                    False, blocking=False),
    ]
    return gates.Gate(STAGE, "Dedicated pLN per-section protein rule evidence for review; "
                      "no final class or release admission.", checks,
                      ["pln_marker_support.png", "pln_candidate_status.png"] + [
                          f"pln_rule_evidence_{entry}.png"
                          for entry in task.settings["protein_entries"]])
