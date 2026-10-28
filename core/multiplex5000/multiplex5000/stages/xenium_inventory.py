"""Stage 00 for an explicitly selected Xenium-only RNA typing project.

The full IO60 inventory still owns mixed-modality coverage. This route admits
only the Xenium regions named by the project, while retaining the recorded
donor and pLN context. It never asks for a QuPath image or infers a same-cell
pair from a donor or slide number.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from .. import acquisition_records, gates, panel_names, profiles, relationships, roots
from ..integration import tissues
from ..integration.manifest import parse_bundle_path
from ..tasks import Input, Task
from ..xenium_sources import INSTRUMENT_FILES

STAGE = "00_inventory"
_COLUMNS = ("platform", "acquisition_id", "donor_id", "roi", "tissue", "output", "path",
            "condition", "panel", "run", "run_date", "batch_number", "xenium_serial",
            "recorded_region", "section_mode", "same_block_pln_on_pheno_slide",
            "different_block_pln_on_pheno_slide", "lymphoid_in_xenium_pancreas",
            "note")
_INDEX_COLUMNS = ("acquisition_id", "output", "source_path")


class XeniumInventoryError(ValueError):
    """The selected region cannot be tied to its recorded source."""


def _rows(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def _source(project, name: str, default: Path | None = None) -> Path:
    path = project.source(name) or default
    if path is None:
        raise XeniumInventoryError(f"Xenium-only Stage 00 needs sources.{name}")
    return Path(path)


def _selection(project) -> list[str]:
    values = project.settings.get("inventory_acquisitions")
    if (not isinstance(values, list) or not values
            or any(not isinstance(value, str) or not value.startswith("xenium:")
                   for value in values)
            or len(values) != len(set(values))):
        raise XeniumInventoryError("inventory_acquisitions needs distinct Xenium IDs")
    typed = project.settings.get("rna_acquisitions")
    if typed is not None and not set(typed) <= set(values):
        raise XeniumInventoryError("RNA selection is outside the Xenium inventory")
    return values


def _catalog(paths: Path, index: Path, selected: list[str]) -> list[dict[str, str]]:
    listed = {}
    for row in _rows(paths):
        output = Path(row["source_path"].strip()).name
        if not output or output in listed:
            raise XeniumInventoryError(f"repeated or empty Xenium output in {paths}: {output}")
        listed[output] = row
    with Path(index).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != _INDEX_COLUMNS:
            raise XeniumInventoryError("Xenium source index needs acquisition_id,output,source_path")
        indexed = list(reader)
    if len({row["acquisition_id"] for row in indexed}) != len(indexed):
        raise XeniumInventoryError("repeated acquisition ID in Xenium source index")
    if len({row["output"] for row in indexed}) != len(indexed):
        raise XeniumInventoryError("repeated output in Xenium source index")
    if not {row["output"] for row in indexed} <= set(listed):
        raise XeniumInventoryError("site source index names an untracked Xenium output")
    by_id = {row["acquisition_id"]: row for row in indexed}
    if not set(selected) <= set(by_id):
        raise XeniumInventoryError(f"selected Xenium IDs absent from site index: "
                                   f"{sorted(set(selected) - set(by_id))}")
    chosen = []
    for acquisition in selected:
        row = by_id[acquisition]
        source = Path(row["source_path"])
        serial, _ = parse_bundle_path(row["output"])
        if (not source.is_absolute() or source.name != row["output"]
                or not serial or acquisition.split(":")[1] != serial):
            raise XeniumInventoryError(f"{acquisition}: invalid exact Xenium source")
        recorded = listed[row["output"]]
        tissue = tissues.for_roi(recorded["roi"])
        if not tissue:
            raise XeniumInventoryError(f"{acquisition}: unrecognized ROI {recorded['roi']}")
        chosen.append({**row, "donor_id": recorded["donor_id"],
                       "roi": recorded["roi"], "tissue": tissue,
                       "batch": recorded.get("batch", "")})
    if len({row["source_path"] for row in chosen}) != len(chosen):
        raise XeniumInventoryError("selected Xenium regions reuse a source path")
    return chosen


def settings(project) -> dict:
    return {"inventory_scope": "xenium_only", "inventory_acquisitions": _selection(project)}


def inputs(project) -> list[Input]:
    paths = _source(project, "xenium_paths", roots.INTEGRATION / "xenium_paths.csv")
    index = _source(project, "xenium_source_index")
    condition = _source(project, "tissue_condition", roots.INTEGRATION / "tissue_condition.csv")
    selected = _catalog(paths, index, _selection(project))
    damaged = {row["acquisition_id"] for row in _rows(condition)
               if row["acquisition_id"].startswith("xenium:")}
    found = [Input("xenium_paths", paths), Input("xenium_source_index", index),
             Input("tissue_condition", condition),
             Input("section_modes", _source(project, "section_modes",
                                           roots.INTEGRATION / "section_modes.csv")),
             Input("xenium_lymphoid_inspection",
                   _source(project, "xenium_lymphoid_inspection",
                           roots.INTEGRATION / "xenium_lymphoid_inspection.csv")),
             Input("donor_notes", _source(project, "donor_notes")),
             Input("donor_notes_copy", _source(project, "donor_notes_copy")),
             Input("pln_blocks", _source(project, "pln_blocks")),
             Input("donor_metadata", _source(project, "donor_metadata")),
             Input("npod_cases", _source(project, "npod_cases"))]
    for row in selected:
        for filename in INSTRUMENT_FILES:
            found.append(Input(f"xenium:{row['output']}:{filename}",
                               Path(row["source_path"]) / filename,
                               required=row["acquisition_id"] not in damaged))
    return found


def _write(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run(task: Task) -> dict:
    from . import inventory

    selected = _catalog(task.input("xenium_paths"), task.input("xenium_source_index"),
                        task.settings["inventory_acquisitions"])
    cohort = profiles.load_all()[("cohort", task.project.cohort)]
    aliases = cohort.donor_aliases
    reading = relationships.load(
        notes=task.input("donor_notes"), notes_copy=task.input("donor_notes_copy"),
        pln_blocks=task.input("pln_blocks"), section_modes=task.input("section_modes"),
        aliases=aliases, cohort_donors=[])
    conditions = {row["acquisition_id"]: row for row in _rows(task.input("tissue_condition"))
                  if row["acquisition_id"].startswith("xenium:")}
    inspections = {row["slide"]: row for row in _rows(task.input("xenium_lymphoid_inspection"))}
    rows, issues, missing = [], [], []
    panel_rule = profiles.load_all()[("platform", "xenium")].panel_names
    for entry in selected:
        acquisition, output = entry["acquisition_id"], entry["output"]
        donor = aliases.get(entry["donor_id"], entry["donor_id"])
        condition = conditions.get(acquisition)
        if condition and (aliases.get(condition["donor_id"], condition["donor_id"]) != donor
                          or entry["tissue"] not in
                          {tissues.normalize(t) for t in condition["tissue"].split(";")}):
            issues.append(f"{acquisition}: damage table names another donor or tissue")
        files = {name: task.input(f"xenium:{output}:{name}") for name in INSTRUMENT_FILES}
        absent = [name for name, path in files.items() if not path.is_file()]
        if absent and not condition:
            issues.append(f"{acquisition}: missing instrument record(s) {absent}")
        if absent:
            missing.append(acquisition)
            record = None
        else:
            try:
                record = acquisition_records.read_xenium_records(
                    *(files[name] for name in INSTRUMENT_FILES))
            except (OSError, ValueError, KeyError) as err:
                issues.append(f"{acquisition}: {err}")
                record = None
        serial = acquisition.split(":")[1]
        panel = ""
        if record:
            if (record.slide_id != serial or
                    record.region != acquisition.split(":", 2)[2] or record.differences):
                issues.append(f"{acquisition}: instrument identity differs: "
                              f"{record.slide_id}/{record.region} {record.differences}")
            try:
                panel = panel_names.xenium(record.experiment, panel_rule)
            except panel_names.PanelNameError as err:
                issues.append(f"{acquisition}: {err}")
        note = reading.donors.get(donor)
        modes = [mode for mode in reading.modes
                 if mode.donor_id == donor and mode.roi == entry["roi"]
                 and mode.xenium_serial == serial]
        inspection = inspections.get(serial) if entry["tissue"] == tissues.PANCREAS else None
        rows.append({
            "platform": "xenium", "acquisition_id": acquisition, "donor_id": donor,
            "roi": entry["roi"], "tissue": entry["tissue"], "output": output,
            "path": entry["source_path"], "condition": condition["condition"] if condition else "",
            "panel": panel, "run": record.run_name if record else "",
            "run_date": record.run_date if record else "", "batch_number": entry["batch"],
            "xenium_serial": serial, "recorded_region": record.region if record else "",
            "section_mode": modes[0].mode if len(modes) == 1 else "unknown",
            "same_block_pln_on_pheno_slide": (note.same_block_pln if note and
                                                note.same_block_pln is not None else ""),
            "different_block_pln_on_pheno_slide": (note.different_block_pln if note and
                                                     note.different_block_pln is not None else ""),
            "lymphoid_in_xenium_pancreas": (inspection["lymphoid_in_pancreas"]
                                            if inspection else ""),
            "note": "Xenium only; no protein-cell pairing claimed",
        })
    donors = {row["donor_id"] for row in rows}
    statuses, status_disagree, status_unknown = inventory._statuses(task, cohort, donors)
    _write(task.output("acquisitions.csv"), _COLUMNS, rows)
    _write(task.output("source_index.csv"), _INDEX_COLUMNS,
           [{key: entry[key] for key in _INDEX_COLUMNS} for entry in selected])
    _figure(task, rows)
    return {
        "acquisitions": len(rows), "selected": len(task.settings["inventory_acquisitions"]),
        "tissues": dict(Counter(row["tissue"] for row in rows)),
        "issues": issues, "missing_damaged_records": missing,
        "donors_without_status": sorted(donors - set(statuses)),
        "status_disagreements": status_disagree + status_unknown,
        "relationship_conflicts": [f"{c.donor_id}: {c.what}" for c in reading.conflicts],
        "unknown_section_modes": [row["acquisition_id"] for row in rows
                                  if row["section_mode"] == "unknown"],
        "scope": "selected Xenium regions only; no PhenoCycler image inputs",
    }


def _figure(task: Task, rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    donors = sorted({row["donor_id"] for row in rows})
    panc = [sum(row["donor_id"] == donor and row["tissue"] == tissues.PANCREAS
                for row in rows) for donor in donors]
    pln = [sum(row["donor_id"] == donor and row["tissue"] == tissues.LYMPH_NODE
               for row in rows) for donor in donors]
    fig, ax = plt.subplots(figsize=(11, max(3, 1.4 + .38 * len(donors))),
                           layout="constrained")
    y = list(range(len(donors)))
    ax.barh(y, panc, color="#287ba8", label="pancreas")
    ax.barh(y, pln, left=panc, color="#ce6b27", label="pLN")
    ax.set(yticks=y, yticklabels=donors, xlabel="Recorded Xenium regions (n)",
           title=f"Xenium-only source inventory: {len(rows)} selected regions")
    ax.invert_yaxis()
    ax.legend(loc="lower right")
    fig.savefig(gates.figure_path(task, "xenium_source_coverage.png"), dpi=160)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [
        gates.Check("selected Xenium acquisitions", summary["acquisitions"],
                    str(summary["selected"]), summary["acquisitions"] == summary["selected"]),
        gates.Check("instrument identity, panel and damage agreement", summary["issues"],
                    "none", not summary["issues"]),
        gates.Check("donor disease status", summary["donors_without_status"],
                    "none", not summary["donors_without_status"]),
        gates.Check("status sources agree", summary["status_disagreements"],
                    "none", not summary["status_disagreements"]),
        gates.Check("damaged instrument records unavailable",
                    summary["missing_damaged_records"], "reported", True, blocking=False),
        gates.Check("recorded section relationships unresolved",
                    summary["unknown_section_modes"], "reported; no same-cell claim",
                    True, blocking=False),
        gates.Check("relationship source disagreements", summary["relationship_conflicts"],
                    "reported; no relationship inferred", True, blocking=False),
        gates.Check("source scope", summary["scope"], "selected Xenium regions only", True,
                    blocking=False),
    ]
    return gates.Gate(STAGE, "Exact Xenium source inventory for RNA-only typing; "
                      "no protein measurement or paired-cell claim.", checks,
                      ["xenium_source_coverage.png"])
