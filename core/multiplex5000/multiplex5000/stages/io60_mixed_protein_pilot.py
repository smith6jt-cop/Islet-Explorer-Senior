"""Stage 06 pilot: selected pancreas hierarchy plus same-slide pLN sheet evidence.

This project types only its named donors. It reads one accepted IO60 Stage 00/05
release and never publishes final pLN classes or activates a release. The old
pancreas method runs from this package's tracked source (or wheel snapshot).
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from .. import cohort_scale, gates, roots
from ..identity import check_cell_keys
from ..tasks import Input, Task
from ..tissues.evaluate import resolve_unique_matched_rules
from ..tissues.lymph_node import load_rule_source
from ..tissues.pancreas_mixed import assemble
from ..tissues.protein_section import evaluate_protein_section

STAGE = "06_typing"
PANEL = roots.REPO / "core/panels/data/tissues/pancreatic_lymph_node/manifest.yml"
RULES = PANEL.parent / "phenotype_rules.csv"
COHORT = roots.REPO / "core/multiplex5000/cohort-spec-v3.json"
RUNNER = roots.REPO / "core/multiplex5000/scripts/typing/run_calibration_arm.py"
STRUCTURE_FOLDERS = ("Blood_vessel", "Lymphatic", "Nerve")
META = {"donor_id", "object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
        "cell_region", "islet_num", "in_alignment_block"}


class IO60MixedPilotError(ValueError):
    """The selected mixed-slide pilot has missing, ambiguous or cross-tissue inputs."""


def _source(project, name: str) -> Path:
    value = project.source(name)
    if value is None:
        raise IO60MixedPilotError(f"Stage 06 needs sources.{name} in the site or project file")
    return Path(value)


def _cohort_rows() -> list[dict]:
    rows = json.loads(COHORT.read_text())["images"]
    donors = [str(row["donor_id"]) for row in rows]
    if len(donors) != len(set(donors)):
        raise IO60MixedPilotError("the v3 cohort has repeated donors")
    return rows


def _donors(project) -> list[str]:
    value = project.settings.get("typing_donors")
    if (not isinstance(value, list) or not value
            or not all(isinstance(donor, str) and donor for donor in value)
            or len(value) != len(set(value))):
        raise IO60MixedPilotError("typing_donors must be a nonempty unique list of text IDs")
    known = {str(row["donor_id"]) for row in _cohort_rows()}
    if not set(value) <= known or "6457" in value:
        raise IO60MixedPilotError("typing_donors names an excluded or unknown IO60 donor")
    return value


def _upstream(project) -> Path:
    if project.settings.get("typing_route") == "io60_joint_typing":
        return Path(project.release_dir)
    return _source(project, "io60_upstream_release")


def _aligned(project, donor: str) -> Path:
    return _upstream(project) / f"05_normalization/data/aligned/donor_id={donor}/data.parquet"


def _exports(project, donor: str) -> list[Input]:
    manifest = json.loads(_source(project, "qupath_manifest").read_text())
    rows = [row for row in manifest["images"] if str(row["donor_id"]) == donor]
    if len(rows) != 1:
        raise IO60MixedPilotError(f"{donor}: expected one recorded QuPath image")
    stem = rows[0]["image_id"].split(" - ")[0].replace(".qptiff", "")
    root = _source(project, "typing_export_root")
    found = []
    for folder in STRUCTURE_FOLDERS:
        matches = sorted(path for path in (root / "structures" / folder).glob("*.geojson")
                         if stem in path.name)
        if len(matches) != 1:
            raise IO60MixedPilotError(f"{donor}: {folder} has {len(matches)} matching exports")
        found.append(Input(f"structure:{donor}:{folder}", matches[0]))
    islets = sorted(path for path in (root / "typing_islets").rglob("*.geojson")
                    if stem in path.name)
    if not islets:
        raise IO60MixedPilotError(f"{donor}: no islet annotation export")
    # The selected resolver checks every matching file by newest mtime and
    # refuses detection exports. Declare all candidates, including rejected ones.
    found.extend(Input(f"islet_export:{donor}:{index}", path,
                       hash=path.stat().st_size < 1024 ** 3)
                 for index, path in enumerate(islets))
    return found


def settings(project) -> dict:
    return {"donors": _donors(project), "scope": "bounded_review_pilot",
            "upstream_release": str(_upstream(project)),
            "pancreas_method": "selected_hierarchy_with_single_donor_representation_fit",
            "pln_method": "full_6551_sheet_section_local_0.10_candidate"}


def inputs(project) -> list[Input]:
    upstream = _upstream(project)
    geometry = _source(project, "typing_geometry_root")
    found = [Input("cohort_spec", COHORT), Input("pancreas_runner", RUNNER),
             Input("pln_panel", PANEL), Input("pln_rules", RULES),
             Input("inventory_acquisitions", upstream / "00_inventory/acquisitions.csv"),
             Input("inventory_receipt", upstream / "00_inventory/receipt.json"),
             Input("inventory_gate", upstream / "00_inventory/gate/gate.json"),
             Input("inventory_acceptance", upstream / "00_inventory/gate/accepted.json"),
             Input("normalization_receipt", upstream / "05_normalization/receipt.json"),
             Input("normalization_gate", upstream / "05_normalization/gate/gate.json"),
             Input("normalization_acceptance", upstream / "05_normalization/gate/accepted.json"),
             Input("normalization_scales", upstream / "05_normalization/scales.csv"),
             Input("qupath_manifest", _source(project, "qupath_manifest")),
             Input("donor_acquisition_runs", _source(project, "donor_acquisition_runs")),
             Input("typing_registry_log", _source(project, "typing_registry_log")),
             Input("typing_shifts", _source(project, "typing_shifts"))]
    for row in _cohort_rows():
        donor = str(row["donor_id"])
        found.append(Input(f"aligned:{donor}", _aligned(project, donor)))
    for donor in _donors(project):
        found.extend((Input(f"revised:{donor}",
                            geometry / f"revised_cells/donor_id={donor}/data.parquet"),
                      Input(f"crossed:{donor}", geometry / f"crossed/crossed_{donor}.parquet"),
                      Input(f"fundamentals:{donor}",
                            geometry / f"qupath_fundamentals/donor_id={donor}/data.parquet"),
                      Input(f"pln_blocks:{donor}",
                            upstream / f"05_normalization/alignment_block/donor_id={donor}.parquet",
                            required=False)))
        found.extend(_exports(project, donor))
    return found


def _accepted(task: Task, name: str, stage: str) -> None:
    receipt = json.loads(task.input(name).read_text())
    if receipt.get("stage") != stage or receipt.get("status") != "complete":
        raise IO60MixedPilotError(f"upstream {stage} has no complete receipt")


def _marker_roster(task: Task) -> list[str]:
    with task.input("normalization_scales").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    # A marker refused for every donor (CD57 in this selected Stage 05 run)
    # has no aligned observation. It cannot enter a comparison as zero.
    markers = sorted({row["marker"] for row in rows if row["status"] == "aligned"})
    donors = {row["donor_id"] for row in rows}
    if not markers or donors != {str(row["donor_id"]) for row in _cohort_rows()}:
        raise IO60MixedPilotError("Stage 05 scale roster differs from the v3 cohort")
    return markers


def _acquisition(task: Task, donor: str) -> tuple[str, str]:
    table = pd.read_csv(task.input("inventory_acquisitions"), dtype=str)
    rows = table.loc[table.donor_id.eq(donor) & table.source.eq("qupath_manifest")
                     & table.tissue.eq("pancreas")]
    if len(rows) != 1 or rows.acquisition_id.isna().any() or rows.segmentation_version.isna().any():
        raise IO60MixedPilotError(f"{donor}: Stage 00 has no unique pancreas acquisition key")
    row = rows.iloc[0]
    return row.acquisition_id, row.segmentation_version


def _protein_candidate(task: Task, donor: str, aligned: pd.DataFrame):
    node = aligned.loc[aligned.cell_region.eq("pLN")]
    if node.empty:
        return pd.DataFrame({"object_id": [], "candidate_phenotype": [],
                             "candidate_state": []}), {}
    source = load_rule_source(task.input("pln_panel"), table_path=task.input("pln_rules"))
    measurements = {marker: node[marker].to_numpy("float64")
                    for marker in node.columns if marker not in META}
    # The selected sheet uses PD-L1; dedicated pLN exports call it PDL1,
    # while this recorded IO60 panel calls the same channel PD_L1.
    if "PD_L1" in measurements:
        if "PDL1" in measurements:
            raise IO60MixedPilotError(f"{donor}: both PD-L1 export spellings are present")
        measurements["PDL1"] = measurements.pop("PD_L1")
    evidence = evaluate_protein_section(
        source, donor, measurements, node.qc_analysis_eligible.to_numpy(bool),
        # Stage 05 has already multiplied every marker by its UniFORM shift.
        {marker: 1.0 for marker in measurements})
    calls = resolve_unique_matched_rules(source, evidence.rules)
    candidate = pd.DataFrame({"object_id": node.object_id.to_numpy(),
                              "candidate_phenotype": calls.phenotype,
                              "candidate_state": calls.status})
    marker_counts = {marker: (int(np.count_nonzero(item.states == 1))
                              if item.cutoff is not None else None)
                     for marker, item in evidence.markers.items()}
    for rule in source.rules:
        for marker in rule.markers:
            marker_counts.setdefault(source.panel_aliases.get(marker, marker), None)
    return candidate, marker_counts


def _published_inner_manifest(task: Task, record: dict) -> dict:
    """Point recorded outputs at the path the task contract will publish."""
    published = task.project.release_dir / STAGE
    record["typing_root"] = str(published / "pancreas")
    record["scale"] = str(published / "pancreas_scale.parquet")
    for step in record.get("stages", []):
        step["log"] = str(published / "pancreas" / Path(step["log"]).name)
    # Commands preserve the actual staging path used during execution, while
    # output fields resolve after the atomic rename into the published stage.
    record["publication_note"] = "Commands ran in task staging; output paths name the published stage"
    return record


def run(task: Task) -> dict:
    _accepted(task, "inventory_receipt", "00_inventory")
    _accepted(task, "normalization_receipt", "05_normalization")
    gates.require_accepted(_upstream(task.project) / "00_inventory")
    gates.require_accepted(_upstream(task.project) / "05_normalization")
    markers = _marker_roster(task)
    aligned_inputs = {str(row["donor_id"]): task.input(f"aligned:{row['donor_id']}")
                      for row in _cohort_rows()}
    scale = cohort_scale.build_scale(aligned_inputs, markers)
    scale.reset_index().to_parquet(task.output("pancreas_scale.parquet"), index=False)
    runner = [sys.executable, str(task.input("pancreas_runner")),
              "--aligned", str(_upstream(task.project) / "05_normalization/data/aligned"),
              "--scale", str(task.staging / "pancreas_scale.parquet"),
              "--geometry-root", str(_source(task.project, "typing_geometry_root")),
              "--registry-log", str(task.input("typing_registry_log")),
              "--shifts", str(task.input("typing_shifts")),
              "--donors", *task.settings["donors"],
              "--representation-donors", *task.settings["donors"],
              "--out", str(task.staging / "pancreas")]
    env = os.environ.copy()
    env.update(IO60_EXPORT_ROOT=str(_source(task.project, "typing_export_root")),
               IO60_QUPATH_MANIFEST=str(task.input("qupath_manifest")),
               IO60_DONOR_ACQUISITION_RUNS=str(task.input("donor_acquisition_runs")))
    result = subprocess.run(runner, env=env, cwd=roots.REPO, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    task.output("pancreas_runner.log").write_text(result.stdout)
    if result.returncode:
        raise IO60MixedPilotError(f"pancreas hierarchy exited {result.returncode}: "
                                  f"{result.stdout[-1500:]}")
    inner = json.loads((task.staging / "pancreas/run.json").read_text())
    if inner.get("status") != "complete" or inner.get("donors") != task.settings["donors"]:
        raise IO60MixedPilotError("pancreas hierarchy did not record complete pilot scope")
    inner = _published_inner_manifest(task, inner)
    (task.staging / "pancreas/run.json").write_text(json.dumps(inner, indent=2, sort_keys=True)
                                                        + "\n")

    totals = {"segmented_cells": 0, "eligible_pancreas": 0, "pln_cells": 0,
              "pln_eligible": 0, "unique_cell_uids": 0}
    states: dict[str, int] = {}
    pln_states: dict[str, int] = {}
    classes: dict[str, int] = {}
    marker_support: dict[str, dict[str, int]] = {}
    for donor in task.settings["donors"]:
        aligned = pq.read_table(task.input(f"aligned:{donor}")).to_pandas()
        pancreas = pq.read_table(task.staging / "pancreas/cells"
                                 / f"donor_id={donor}/data.parquet",
                                 columns=["object_id", "cell_class", "refined_class"]).to_pandas()
        pln, marker_support[donor] = _protein_candidate(task, donor, aligned)
        block_path = task.input(f"pln_blocks:{donor}")
        blocks = (pd.read_parquet(block_path, columns=["object_id", "block"])
                  if block_path.exists() else None)
        acquisition, segmentation = _acquisition(task, donor)
        joined = assemble(aligned[["object_id", "qc_analysis_eligible", "cell_region",
                                   "in_alignment_block"]], pancreas, pln, donor_id=donor,
                          acquisition_id=acquisition, segmentation_id=segmentation,
                          pln_blocks=blocks)
        check_cell_keys(joined)
        joined.to_parquet(task.output(f"mixed/donor_id={donor}/data.parquet"),
                          index=False, compression="zstd")
        pln.to_parquet(task.output(f"pln_candidates/donor_id={donor}/data.parquet"),
                       index=False, compression="zstd")
        totals["segmented_cells"] += len(joined)
        totals["eligible_pancreas"] += int((joined.typing_state == "pancreas_hierarchy").sum())
        totals["pln_cells"] += int(joined.tissue.eq("pancreatic_lymph_node").sum())
        totals["pln_eligible"] += int((joined.tissue.eq("pancreatic_lymph_node")
                                       & joined.qc_analysis_eligible).sum())
        totals["unique_cell_uids"] += joined.cell_uid.nunique()
        for name, n in joined.typing_state.value_counts().items():
            states[name] = states.get(name, 0) + int(n)
        for name, n in joined.loc[joined.tissue.eq("pancreatic_lymph_node"),
                                  "typing_state"].value_counts().items():
            pln_states[name] = pln_states.get(name, 0) + int(n)
        for name, n in pancreas.refined_class.fillna("unavailable").value_counts().items():
            classes[name] = classes.get(name, 0) + int(n)
    _figures(task, totals, states, pln_states, classes, marker_support)
    legacy = sorted({str(item.path.resolve()) for item in task.inputs.values()
                     if "Islet-Explorer-Senior-v3" in str(item.path.resolve())})
    task.output("typing_counts.json").write_text(json.dumps(
        {"totals": totals, "states": states, "pln_states": pln_states,
         "pancreas_refined_classes": classes,
         "pln_marker_positive": marker_support,
         "legacy_checkout_inputs": legacy}, indent=2, sort_keys=True) + "\n")
    return {**totals, "donors": task.settings["donors"], "typing_states": states,
            "pln_states": pln_states,
            "pancreas_refined_classes": classes, "candidate_only": True,
            "final_pln_cell_classes": 0, "representation_donors": task.settings["donors"],
            "legacy_checkout_inputs": legacy}


def _figures(task: Task, totals: dict, states: dict, pln_states: dict, classes: dict,
             marker_support: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["pancreas hierarchy", "pLN sheet candidate", "pancreas QC unavailable"]
    values = [totals["eligible_pancreas"], totals["pln_cells"],
              states.get("unavailable_qc", 0) -
              (totals["pln_cells"] - totals["pln_eligible"])]
    fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    ax.bar(labels, values)
    ax.set(ylabel="Segmented cells (n)", title=f"Mixed IO60 slide route; "
           f"{totals['segmented_cells']:,} cells across {len(task.settings['donors'])} donor(s)")
    for x, value in enumerate(values):
        ax.text(x, value, f"{value:,}", ha="center", va="bottom")
    fig.savefig(gates.figure_path(task, "mixed_tissue_coverage.png"), dpi=140)
    plt.close(fig)

    candidate_states = ["unique_candidate", "ambiguous", "no_match", "unavailable",
                        "unavailable_qc"]
    numbers = [pln_states.get(name, 0) for name in candidate_states]
    fig, ax = plt.subplots(figsize=(9, 4), layout="constrained")
    ax.bar(candidate_states, numbers)
    ax.set(ylabel="pLN cells (n)", title="6551 sheet: review candidates, not accepted classes")
    ax.tick_params(axis="x", rotation=22)
    fig.savefig(gates.figure_path(task, "pln_candidate_states.png"), dpi=140)
    plt.close(fig)

    grouped = sorted(classes.items(), key=lambda item: (-item[1], item[0]))[:20]
    fig, ax = plt.subplots(figsize=(9, max(5, len(grouped) * 0.28)), layout="constrained")
    ax.barh([name for name, _ in reversed(grouped)], [n for _, n in reversed(grouped)])
    ax.set(xlabel="Eligible pancreas cells (n)",
           title="Selected pancreas hierarchy, final refinement step; pilot fit")
    fig.savefig(gates.figure_path(task, "pancreas_refined_classes.png"), dpi=140)
    plt.close(fig)

    markers = sorted({marker for support in marker_support.values() for marker in support})
    donors = task.settings["donors"]
    matrix = np.array([[marker_support[donor].get(marker) if
                        marker_support[donor].get(marker) is not None else np.nan
                        for donor in donors]
                       for marker in markers], dtype=float)
    fig, ax = plt.subplots(figsize=(max(7, len(donors) * 1.8),
                                    max(9, len(markers) * 0.23)), layout="constrained")
    image = ax.imshow(np.ma.masked_invalid(matrix), aspect="auto",
                      cmap=plt.get_cmap("viridis").with_extremes(bad="#d9d9d9"))
    ax.set(xticks=range(len(donors)), xticklabels=donors,
           yticks=range(len(markers)), yticklabels=markers,
           xlabel=f"Mixed slide donor (n={totals['pln_eligible']:,} eligible pLN cells total)",
           ylabel="Measured protein marker",
           title="Section-local 0.10 pLN marker support (n); gray = unavailable")
    ax.tick_params(axis="y", labelsize=7)
    fig.colorbar(image, ax=ax, label="Candidate-positive eligible pLN cells (n)")
    fig.savefig(gates.figure_path(task, "pln_marker_support.png"), dpi=140)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    total = summary["segmented_cells"]
    states = summary["typing_states"]
    return gates.Gate(STAGE, "Bounded mixed pancreas/pLN typing pilot; review evidence only. "
                      "The pancreas refinement fit uses only the named pilot donors.", [
        gates.Check("all segmented cells retain unique acquisition-scoped keys",
                    summary["unique_cell_uids"], str(total),
                    summary["unique_cell_uids"] == total),
        gates.Check("every segmented cell has a tissue typing state", sum(states.values()),
                    str(total), sum(states.values()) == total),
        gates.Check("pancreas hierarchy covers eligible non-pLN cells",
                    summary["eligible_pancreas"], ">0", summary["eligible_pancreas"] > 0),
        gates.Check("pLN sheet candidate covers recorded pLN cells", summary["pln_cells"],
                    ">0", summary["pln_cells"] > 0),
        gates.Check("pLN candidates are not final classes", summary["final_pln_cell_classes"],
                    "0", summary["final_pln_cell_classes"] == 0),
        gates.Check("biology and control review", "pending", "maintainer review", False,
                    blocking=False,
                    detail="FOXP3 scale and tissue-specific marker controls require review"),
        gates.Check("inputs still resolving into the v3 checkout",
                    summary["legacy_checkout_inputs"], "none before promotion",
                    not summary["legacy_checkout_inputs"], blocking=False,
                    detail="Copy and verify retained inputs before any active release"),
    ], ["mixed_tissue_coverage.png", "pln_candidate_states.png",
        "pancreas_refined_classes.png", "pln_marker_support.png"])
