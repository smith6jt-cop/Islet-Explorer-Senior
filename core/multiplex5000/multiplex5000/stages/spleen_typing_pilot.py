"""Stage 06 pilot: separate protein-rule and RNA-gene evidence on one spleen section.

The protein export and Xenium region output share a physical section, not cell identities. This
stage does not spatially link their segmented cells or publish a final class.
Its two assay tables and review figures feed the registration/integration work.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .. import gates, roots
from ..rna_evidence import detect_genes
from ..same_slide_sources import load_pairs
from ..tasks import Input, Task, sha256_file
from ..tissues.evaluate import CANDIDATE_STATES, resolve_unique_matched_rules
from ..tissues.protein_section import evaluate_protein_section
from ..tissues.spleen import load_rule_source
from ..xenium_typing_input import load_accepted_region
from .pilot_inventory import _manifest

STAGE = "06_typing"
PANEL = roots.REPO / "core/panels/data/tissues/spleen/manifest.yml"
IDENTITY = roots.REPO / "core/panels/data/identity_core"
GENE_SHEETS = tuple(sorted(IDENTITY.glob("0[3456]_immune_*.csv")))
_THRESHOLD_COLUMNS = ("protein_entry", "marker", "floor", "top", "cutoff", "shift_scale",
                      "n_background_support", "unavailable_reason", "detected", "undetected",
                      "unavailable")


class SpleenTypingPilotError(ValueError):
    """The accepted assay sources cannot support this bounded typing pilot."""


def _selected(project) -> list[tuple[object, dict]]:
    manifest = project.settings.get("same_slide_source_manifest")
    chosen = project.settings.get("rna_acquisitions")
    if not manifest or not isinstance(chosen, list) or not chosen:
        raise SpleenTypingPilotError("Stage 06 needs recorded same-slide RNA selections")
    pairs = load_pairs(_manifest(project))
    with (project.release_dir / "04_rna_ingest/regions.csv").open(newline="") as handle:
        rows = {row["acquisition_id"]: row for row in csv.DictReader(handle)}
    found = []
    for pair in pairs:
        if pair.xenium_acquisition_id in chosen:
            row = rows.get(pair.xenium_acquisition_id)
            if (row is None or row["donor_id"] != pair.donor_id or row["tissue"] != "spleen"
                    or Path(row["source_path"]) != pair.xenium_output):
                raise SpleenTypingPilotError(f"{pair.xenium_acquisition_id}: Stage 04 mismatch")
            found.append((pair, row))
    if {pair.xenium_acquisition_id for pair, _ in found} != set(chosen):
        raise SpleenTypingPilotError("Stage 04 selection is not in the same-slide manifest")
    if len(found) != 1:
        raise SpleenTypingPilotError("the first Stage 06 pilot selects one paired section")
    return found


def settings(project) -> dict:
    return {"rna_acquisitions": [pair.xenium_acquisition_id
                                 for pair, _ in _selected(project)],
            "protein_cutoff_fraction": 0.10,
            "scope": "pilot_candidate_no_final_cell_class"}


def inputs(project) -> list[Input]:
    release = project.release_dir
    found = [Input("same_slide_source_manifest", _manifest(project)),
             Input("rna_regions", release / "04_rna_ingest/regions.csv"),
             Input("rna_receipt", release / "04_rna_ingest/receipt.json"),
             Input("rna_gate_summary", release / "04_rna_ingest/gate/gate.json"),
             Input("rna_gate", release / "04_rna_ingest/gate/accepted.json"),
             Input("protein_receipt", release / "03_protein_measurement/receipt.json"),
             Input("protein_gate", release / "03_protein_measurement/gate/accepted.json"),
             Input("normalization_scales", release / "05_normalization/scales.csv"),
             Input("normalization_receipt", release / "05_normalization/receipt.json"),
             Input("normalization_gate", release / "05_normalization/gate/accepted.json"),
             Input("spleen_panel", PANEL),
             Input("pln_panel", PANEL.parent.parent / "pancreatic_lymph_node/manifest.yml"),
             Input("pln_rules", PANEL.parent.parent / "pancreatic_lymph_node/phenotype_rules.csv"),
             Input("pc67_header", PANEL.parent / "pc67_export_header.tsv")]
    found += [Input(f"gene_sheet:{path.stem}", path) for path in GENE_SHEETS]
    for pair, row in _selected(project):
        entry, sample = pair.protein_entry, row["sample_key"]
        folder = Path(row["source_path"])
        found += [Input(f"protein_cells:{entry}", release / "03_protein_measurement/cells"
                        / f"{entry}.parquet"),
                  Input(f"rna_cells:{sample}", release / "04_rna_ingest/cells"
                        / f"{sample}.parquet"),
                  Input(f"raw_matrix:{sample}", folder / "cell_feature_matrix.h5"),
                  Input(f"raw_cells:{sample}", folder / "cells.parquet")]
    return found


def _genes(task: Task) -> tuple[str, ...]:
    names = set()
    for path in GENE_SHEETS:
        with task.input(f"gene_sheet:{path.stem}").open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("gene"):
                    names.add(row["gene"])
    if not names:
        raise SpleenTypingPilotError("identity-core gene panels are empty")
    return tuple(sorted(names))


def _scale_rows(task: Task) -> dict[tuple[str, str], float | None]:
    with task.input("normalization_scales").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    found = {}
    for row in rows:
        key = (row["protein_entry"], row["marker"])
        if key in found:
            raise SpleenTypingPilotError(f"duplicate UniFORM scale: {key}")
        found[key] = float(row["scale"]) if row["scale"] else None
    return found


def run(task: Task) -> dict:
    protein_receipt = json.loads(task.input("protein_receipt").read_text())
    if (protein_receipt.get("status") != "complete"
            or protein_receipt.get("stage") != "03_protein_measurement"):
        raise SpleenTypingPilotError("Stage 03 has no complete protein receipt")
    inputs = {item["name"]: item.get("sha256") for item in protein_receipt["inputs"]}
    if inputs.get("same_slide_source_manifest") != sha256_file(
            task.input("same_slide_source_manifest")):
        raise SpleenTypingPilotError("paired-section manifest changed after Stage 03")
    normalization_receipt = json.loads(task.input("normalization_receipt").read_text())
    if (normalization_receipt.get("status") != "complete"
            or normalization_receipt.get("stage") != "05_normalization"):
        raise SpleenTypingPilotError("Stage 05 has no complete normalization receipt")
    selected = _selected(task.project)
    source = load_rule_source(task.input("spleen_panel"))
    genes = _genes(task)
    shifts = _scale_rows(task)
    threshold_rows, rule_rows, gene_rows, candidate_rows = [], [], [], []
    protein_cells = rna_cells = 0
    for pair, row in selected:
        entry, sample = pair.protein_entry, row["sample_key"]
        region = load_accepted_region(
            row, task.input("rna_receipt"), task.input("rna_regions"),
            task.input(f"rna_cells:{sample}"), task.input(f"raw_matrix:{sample}"),
            task.input(f"raw_cells:{sample}"))
        rna = detect_genes(region.counts, region.genes, genes, region.eligible)
        rna_payload = {"cell_uid": pa.array(region.cell_uid),
                       "typing_eligible": pa.array(region.eligible),
                       **{f"gene:{gene}": pa.array(states)
                          for gene, states in rna.states.items()}}
        pq.write_table(pa.table(rna_payload), task.output(f"rna/{sample}.parquet"),
                       compression="zstd")
        for gene, states in rna.states.items():
            gene_rows.append({"sample_key": sample, "gene": gene,
                              "measured": gene not in rna.missing_genes,
                              "detected": int(np.count_nonzero(states == 1)),
                              "undetected": int(np.count_nonzero(states == 0)),
                              "unavailable": int(np.count_nonzero(states == -1))})
        rna_cells += region.n_cells
        table = pq.read_table(task.input(f"protein_cells:{entry}"))
        names = table.schema.names
        measured = sorted(name.removeprefix("marker:") for name in names
                          if name.startswith("marker:"))
        eligible = table["protein_eligible"].to_numpy(zero_copy_only=False).astype(bool)
        protein_payload = {"cell_uid": table["cell_uid"], "object_id": table["object_id"],
                           "protein_eligible": table["protein_eligible"]}
        measurements = {marker: table[f"marker:{marker}"].to_pandas().to_numpy(dtype=float)
                        for marker in measured}
        section_evidence = evaluate_protein_section(
            source, entry, measurements, eligible,
            {marker: shifts.get((entry, marker)) for marker in measured})
        candidates = resolve_unique_matched_rules(source, section_evidence.rules)
        protein_payload["candidate_phenotype"] = pa.array(candidates.phenotype)
        protein_payload["candidate_state"] = pa.array(candidates.status)
        protein_payload["candidate_unavailable_lineage_rules"] = pa.array(
            candidates.unavailable_lineage_rules)
        for marker, evidence in section_evidence.markers.items():
            protein_payload[f"marker_call:{marker}"] = pa.array(evidence.states)
            threshold_rows.append({"protein_entry": entry, "marker": marker,
                                   "floor": evidence.floor, "top": evidence.top,
                                   "cutoff": evidence.cutoff,
                                   "shift_scale": evidence.shift_scale,
                                   "n_background_support": evidence.n_background_support,
                                   "unavailable_reason": evidence.unavailable_reason,
                                   "detected": int(np.count_nonzero(evidence.states == 1)),
                                   "undetected": int(np.count_nonzero(evidence.states == 0)),
                                   "unavailable": int(np.count_nonzero(evidence.states == -1))})
        for phenotype, states in section_evidence.rules.states.items():
            protein_payload[f"rule:{phenotype}"] = pa.array(states)
            rule_rows.append({"protein_entry": entry, "phenotype": phenotype,
                              "matched": int(np.count_nonzero(states == 1)),
                              "not_matched": int(np.count_nonzero(states == 0)),
                              "unavailable": section_evidence.rules.unavailable_counts[phenotype]})
        candidate_rows.extend({"protein_entry": entry, "status": status,
                               "cells": int(np.count_nonzero(candidates.status == status))}
                              for status in CANDIDATE_STATES)
        pq.write_table(pa.table(protein_payload), task.output(f"protein/{entry}.parquet"),
                       compression="zstd")
        protein_cells += len(table)
    for filename, rows in (("protein_thresholds.csv", threshold_rows),
                           ("protein_rule_counts.csv", rule_rows),
                           ("protein_candidate_status.csv", candidate_rows),
                           ("rna_gene_counts.csv", gene_rows)):
        with task.output(filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    _figure(task, threshold_rows, rule_rows, gene_rows, protein_cells, rna_cells)
    _candidate_figure(task, candidate_rows, protein_cells)
    return {"pairs": len(selected), "protein_cells": protein_cells, "rna_cells": rna_cells,
            "protein_markers": len(threshold_rows), "protein_rules": len(rule_rows),
            "protein_markers_with_cutoff": sum(row["cutoff"] is not None
                                               for row in threshold_rows),
            "expected_protein_rules": len(selected) * len(source.rules),
            "rna_genes_requested": len(genes) * len(selected),
            "rna_genes_measured": sum(row["measured"] for row in gene_rows),
            "candidate_status_cells": sum(row["cells"] for row in candidate_rows),
            "unique_rule_candidates": sum(row["cells"] for row in candidate_rows
                                          if row["status"] == "unique_candidate"),
            "candidate_only": True, "paired_cell_links": 0}


def _figure(task: Task, protein: list[dict], rules: list[dict], rna: list[dict],
            protein_cells: int, rna_cells: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (left, middle, right) = plt.subplots(1, 3, figsize=(24, 14),
                                             layout="constrained")
    p = sorted(protein, key=lambda row: row["detected"], reverse=True)
    left.barh([row["marker"] for row in p], [row["detected"] for row in p],
              label="marker detected")
    left.set(xlabel="Eligible protein cells above section-local candidate cutoff (n)",
             title=f"PhenoCycler raw-mean candidate: {len(p)} markers, "
                   f"{protein_cells:,} segmented cells")
    left.invert_yaxis()
    grouped = sorted(rules, key=lambda row: row["matched"], reverse=True)
    y = [row["phenotype"] for row in grouped]
    matched = [row["matched"] for row in grouped]
    unavailable = [row["unavailable"] for row in grouped]
    middle.barh(y, matched, label="rule matched")
    middle.barh(y, unavailable, left=matched, label="rule unavailable")
    middle.set(xlabel="PhenoCycler cells (n)",
               title=f"Spleen rule-match candidate: {len(grouped)} rules, "
                     f"{protein_cells:,} segmented cells")
    middle.invert_yaxis()
    middle.legend(fontsize=8)
    top = sorted((row for row in rna if row["measured"]),
                 key=lambda row: row["detected"], reverse=True)[:30]
    right.barh([row["gene"] for row in top], [row["detected"] for row in top])
    right.set(xlabel="Nucleus-eligible Xenium cells meeting fixed-depth gene predicate (n)",
              title=f"Xenium gene evidence: top {len(top)} measured identity genes, "
                    f"{rna_cells:,} segmented cells")
    right.invert_yaxis()
    fig.savefig(gates.figure_path(task, "assay_evidence.png"), dpi=140)
    plt.close(fig)


def _candidate_figure(task: Task, rows: list[dict], protein_cells: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    counts = {row["status"]: row["cells"] for row in rows}
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    bottom = 0
    for status in CANDIDATE_STATES:
        value = counts[status]
        ax.bar(["selected spleen section"], [value], bottom=[bottom],
               label=status.replace("_", " "))
        bottom += value
    ax.set(ylabel="PhenoCycler segmented cells (n)",
           title=f"Spleen lineage candidates among {protein_cells:,} cells; review only")
    ax.legend()
    fig.savefig(gates.figure_path(task, "spleen_candidate_status.png"), dpi=140)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [
        gates.Check("selected physical-section pairs processed", summary["pairs"],
                    str(len(task.settings["rna_acquisitions"])),
                    summary["pairs"] == len(task.settings["rna_acquisitions"])),
        gates.Check("protein and RNA cells retained separately",
                    [summary["protein_cells"], summary["rna_cells"]],
                    "two assay identities; zero inferred cell links",
                    summary["protein_cells"] > 0 and summary["rna_cells"] > 0
                    and summary["paired_cell_links"] == 0),
        gates.Check("candidate typing is not a verified cell class", True,
                    "rule matches and gene evidence only", summary["candidate_only"]),
        gates.Check("all protein hierarchy rules evaluated", summary["protein_rules"],
                    str(summary["expected_protein_rules"]),
                    summary["protein_rules"] == summary["expected_protein_rules"]),
        gates.Check("candidate status accounts for every protein cell",
                    summary["candidate_status_cells"], str(summary["protein_cells"]),
                    summary["candidate_status_cells"] == summary["protein_cells"]),
        gates.Check("unique deepest rule candidates", summary["unique_rule_candidates"],
                    "lineage rules only; overlapping lineages stay ambiguous",
                    True, blocking=False),
        gates.Check("protein candidate cutoffs available",
                    summary["protein_markers_with_cutoff"],
                    f"of {summary['protein_markers']} measured markers",
                    summary["protein_markers_with_cutoff"] > 0),
        gates.Check("assayed identity genes available", summary["rna_genes_measured"],
                    f"of {summary['rna_genes_requested']} requested gene-section pairs",
                    summary["rna_genes_measured"] > 0),
        gates.Check("unmeasured identity genes reported",
                    summary["rna_genes_requested"] - summary["rna_genes_measured"],
                    "unavailable, not negative", True, blocking=False),
    ]
    return gates.Gate(STAGE, "Review-only per-assay typing evidence. Protein rules use "
                      "section-local 0.10 candidate cutoffs; RNA uses fixed-depth gene "
                      "detection. No final class or paired-cell result.", checks,
                      ["assay_evidence.png", "spleen_candidate_status.png"])
