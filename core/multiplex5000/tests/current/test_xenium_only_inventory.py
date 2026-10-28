"""A Xenium-only project can start on another site without QuPath images."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from multiplex5000 import gates, relationships, tasks
from multiplex5000.stages import inventory, rna_ingest
from multiplex5000.xenium_sources import XeniumSourceError, build_roster


def _write(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _bundle(root: Path, serial: str, *, present: bool = True) -> Path:
    path = root / f"output-XETG00298__{serial}__Region_1__20260416__180001"
    if not present:
        return path
    path.mkdir()
    (path / "experiment.xenium").write_text(json.dumps({
        "run_name": "MM050_XeniumBatch15", "run_start_time": "2026-04-16T18:00:01Z",
        "slide_id": serial, "region_name": "Region_1", "panel_design_id": "base",
        "panel_num_targets_predesigned": 5001, "panel_num_targets_custom": 0,
    }))
    (path / "gene_panel.json").write_text(json.dumps({
        "payload": {"panel": {"identity": {"design_id": "base"}}}}))
    _write(path / "metrics_summary.csv", ("run_name", "panel_design_id", "region_name"),
           [{"run_name": "MM050_XeniumBatch15", "panel_design_id": "base",
             "region_name": "Region_1"}])
    return path


def _project(tmp_path: Path):
    selected = [f"xenium:{serial}:Region_1" for serial in ("0041323", "0041326")]
    bundles = [_bundle(tmp_path, serial, present=serial != "0041400")
               for serial in ("0041323", "0041326", "0041400")]
    paths = tmp_path / "xenium_paths.csv"
    _write(paths, ("donor_id", "roi", "batch", "source_path"), [
        {"donor_id": donor, "roi": "panc", "batch": "XeniumOnly",
         "source_path": str(bundle)}
        for donor, bundle in zip(("6512", "6618", "other"), bundles, strict=True)])
    index = tmp_path / "source_index.csv"
    _write(index, ("acquisition_id", "output", "source_path"), [
        {"acquisition_id": f"xenium:{serial}:Region_1",
         "output": bundle.name, "source_path": str(bundle)}
        for serial, bundle in zip(("0041323", "0041326", "0041400"), bundles,
                                  strict=True)])
    conditions = tmp_path / "tissue_condition.csv"
    _write(conditions, ("acquisition_id", "donor_id", "tissue", "condition"), [])
    inspection = tmp_path / "inspection.csv"
    _write(inspection, ("slide", "lymphoid_in_pancreas"), [])
    metadata = tmp_path / "donors.xlsx"
    pd.DataFrame({"donor_id": ["6512", "6618"], "disease.status": ["ND", "T1D"],
                  "sex": ["F", "M"]}).to_excel(metadata, index=False)
    npod = tmp_path / "npod.xlsx"
    pd.DataFrame({"case_id": ["6512", "6618"], "donor_type":
                  ["No Diabetes", "Type 1 Diabetes"], "sex": ["Female", "Male"]}
                 ).to_excel(npod, index=False)
    sources = {"xenium_paths": paths, "xenium_source_index": index,
               "tissue_condition": conditions, "xenium_lymphoid_inspection": inspection,
               "donor_metadata": metadata, "npod_cases": npod}
    for name in ("section_modes", "donor_notes", "donor_notes_copy", "pln_blocks"):
        file = tmp_path / f"{name}.txt"
        file.write_text("")
        sources[name] = file
    return SimpleNamespace(release_dir=tmp_path / "release", release="xenium-only",
                           cohort="io60_v3", path=tmp_path / "project.yaml",
                           settings={"inventory_scope": "xenium_only",
                                     "inventory_acquisitions": selected,
                                     "rna_acquisitions": selected,
                                     "typing_route": "io60_xenium_pilot"},
                           source=lambda name: sources.get(name), previous=lambda stage: None,
                           automatic_gates=())


def test_selected_xenium_regions_run_without_qupath_and_bind_stage04(tmp_path, monkeypatch):
    monkeypatch.setattr(relationships, "load", lambda **kwargs: relationships.Reading())
    project = _project(tmp_path)
    declared = inventory.inputs(project)
    assert not any("qptiff" in item.name or "protein" in item.name for item in declared)
    assert len([item for item in declared if item.name.startswith("xenium:")]) == 6
    result = tasks.run_stage(inventory, project)
    assert result.gate_passed
    stage = project.release_dir / "00_inventory"
    with (stage / "acquisitions.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["acquisition_id"] for row in rows] == project.settings["inventory_acquisitions"]
    assert {row["tissue"] for row in rows} == {"pancreas"}
    assert all(row["same_block_pln_on_pheno_slide"] == "" for row in rows)
    gates.accept(stage, by="test")
    assert rna_ingest._locations(project)[0] == stage / "source_index.csv"
    roster = build_roster(stage, site_index=stage / "source_index.csv",
                          acquisition_ids=set(project.settings["rna_acquisitions"]))
    assert len(roster) == 2
    assert all(row.status == "incomplete" for row in roster)  # no cell data in this fixture
    with (stage / "source_index.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(XeniumSourceError, match="selected source index differs"):
        build_roster(stage, site_index=stage / "source_index.csv")
