"""A single project routes accepted pancreas and pLN RNA acquisitions by tissue."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import h5py
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from multiplex5000 import gates, tasks
from multiplex5000.stages import io60_rna_typing_pilot as pilot
from multiplex5000.stages import rna_ingest, typing
from multiplex5000.xenium_sources import INSTRUMENT_FILES


def _bundle(tmp_path: Path, roi: str, serial: str = "0040001") -> Path:
    path = tmp_path / f"output-XETG00298__{serial}__{roi}__20261026__120000"
    path.mkdir()
    for name in INSTRUMENT_FILES:
        (path / name).write_text(f"{serial} {roi} {name}")
    (path / "cell_boundaries.parquet").write_text("fixture")
    (path / "nucleus_boundaries.parquet").write_text("fixture")
    with h5py.File(path / "cell_feature_matrix.h5", "w") as file:
        matrix = file.create_group("matrix")
        matrix.create_dataset("shape", data=[2, 2])
        matrix.create_dataset("barcodes", data=[b"same-1", b"same-2"])
        features = matrix.create_group("features")
        features.create_dataset("name", data=[b"CD3E", b"TRDC"])
        features.create_dataset("feature_type", data=[b"Gene Expression", b"Gene Expression"])
        matrix.create_dataset("data", data=[2, 1, 3])
        matrix.create_dataset("indices", data=[0, 1, 1])
        matrix.create_dataset("indptr", data=[0, 2, 3])
    pq.write_table(pa.table({"cell_id": ["same-1", "same-2"],
                             "x_centroid": [1.0, 2.0], "y_centroid": [3.0, 4.0],
                             "nucleus_count": [1, 0]}), path / "cells.parquet")
    return path


def _project(tmp_path: Path):
    release = tmp_path / "release"
    inventory = release / "00_inventory"
    inventory.mkdir(parents=True)
    panc, pln = _bundle(tmp_path, "Panc"), _bundle(tmp_path, "PanT")
    damaged = _bundle(tmp_path, "PanT", "0040002")
    metadata = tmp_path / "donors.xlsx"
    pd.DataFrame({"donor_id": ["6539"], "sex": ["M"],
                  "disease.status": ["ND"]}).to_excel(metadata, index=False)
    npod = tmp_path / "npod.xlsx"
    pd.DataFrame({"case_id": ["6539"], "sex": ["Male"],
                  "donor_type": ["No Diabetes"]}).to_excel(npod, index=False)
    fields = ("platform", "acquisition_id", "donor_id", "roi", "tissue", "output", "path",
              "condition", "batch_number", "xenium_serial")
    acquisitions = inventory / "acquisitions.csv"
    with acquisitions.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for path, serial, donor, roi, tissue, condition in (
                (panc, "0040001", "6539", "Panc", "pancreas", ""),
                (pln, "0040001", "6539", "PanT", "pancreatic_lymph_node", ""),
                (damaged, "0040002", "6533", "PanT", "pancreatic_lymph_node", "damaged")):
            writer.writerow({"platform": "xenium", "acquisition_id": f"xenium:{serial}:{roi}",
                             "donor_id": donor, "roi": roi, "tissue": tissue,
                             "output": path.name, "path": str(path), "condition": condition,
                             "batch_number": "10",
                             "xenium_serial": "" if roi == "Panc" else serial})
    input_hashes = [{"name": "donor_metadata", "path": str(metadata),
                     "sha256": tasks.sha256_file(metadata)},
                    {"name": "npod_cases", "path": str(npod),
                     "sha256": tasks.sha256_file(npod)}]
    input_hashes += [{"name": f"xenium:{path.name}:{name}",
                      "sha256": tasks.sha256_file(path / name)}
                     for path in (panc, pln, damaged) for name in INSTRUMENT_FILES]
    (inventory / "receipt.json").write_text(json.dumps({
        "stage": "00_inventory", "status": "complete", "inputs": input_hashes,
        "outputs": [{"path": "acquisitions.csv", "sha256": tasks.sha256_file(acquisitions)}]}))
    gate = inventory / "gate"
    gate.mkdir()
    (gate / "gate.json").write_text(json.dumps({"stage": "00_inventory", "release": "pilot",
                                                "passed": True, "checks": []}))
    gates.accept(inventory, by="test")
    selected = ["xenium:0040001:Panc", "xenium:0040001:PanT", "xenium:0040002:PanT"]
    return SimpleNamespace(release_dir=release, release="pilot", cohort="io60_v3",
                           path=tmp_path / "pilot.yaml",
                           settings={"typing_route": "io60_xenium_pilot",
                                     "rna_acquisitions": selected},
                           source=lambda name: {"donor_metadata": metadata,
                                                "npod_cases": npod}.get(name),
                           previous=lambda stage: "00_inventory" if stage == "04_rna_ingest"
                           else "04_rna_ingest", automatic_gates=("04_rna_ingest",))


@pytest.mark.parametrize("missing_before_ingest", [False, True])
def test_stage06_rna_routes_both_tissues_and_keeps_cell_ids_distinct(
        tmp_path, monkeypatch, missing_before_ingest):
    monkeypatch.setattr(tasks.roots, "tmp_dir", lambda: None)
    project = _project(tmp_path)
    damaged = next(tmp_path.glob("output-*__0040002__PanT__*"))
    if missing_before_ingest:
        (damaged / "cell_feature_matrix.h5").unlink()
        (damaged / "cells.parquet").unlink()
    stage04 = tasks.run_stage(rna_ingest, project)
    assert stage04.gate_passed and stage04.accepted_automatically
    # This selected region is recorded as damaged. Its raw files can become
    # unavailable without holding up the two usable acquisitions.
    if not missing_before_ingest:
        (damaged / "cell_feature_matrix.h5").unlink()
        (damaged / "cells.parquet").unlink()
    seen = []

    def write_calls(donor, section, h5ad, work, acquisition_id, exact_cells):
        data = ad.read_h5ad(h5ad)
        assert data.obs_names.tolist() == ["same-1", "same-2"]
        assert exact_cells == next(path for path in tmp_path.glob("output-*/cells.parquet")
                                   if acquisition_id.split(":")[-1] in path.parent.name)
        work.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"cell_uid": [f"{acquisition_id}:cells:{cell}" for cell in data.obs_names],
                      "cell_id": data.obs_names.tolist(), "acquisition_id": acquisition_id,
                      "lineage": ["Immune_T", "Unavailable_no_nucleus"],
                      "cell_agrees_with_cluster_call": [True, False],
                      "n_naming_detected_for_cluster_call": [2, None],
                      "cd3_detected": [True, False]}).to_parquet(
                          work / f"{section}.calls.parquet", index=False)
        seen.append(acquisition_id)
        return {"acquisition_id": acquisition_id}

    monkeypatch.setattr(pilot.panc_worker, "load_panels", lambda: {})
    monkeypatch.setattr(pilot.pln_worker, "load_panels", lambda: {})

    def panc(donor, section, h5ad, panels, work, **kwargs):
        assert kwargs["naming_policy"] == "shared_er_exact_v1"
        assert kwargs["source_serial"] == "0040001"
        return write_calls(donor, section, h5ad, work, kwargs["acquisition_id"],
                           kwargs["exact_cells"])

    def pln(donor, section, h5ad, cells, panels, work, **kwargs):
        assert kwargs["roi"] == "PanT"
        return write_calls(donor, section, h5ad, work, kwargs["acquisition_id"], cells)

    monkeypatch.setattr(pilot.panc_worker, "process_section", panc)
    monkeypatch.setattr(pilot.pln_worker, "process_exact_section", pln)
    result = tasks.run_stage(typing, project)
    assert result.gate_passed
    assert seen == project.settings["rna_acquisitions"][:2]
    summary = json.loads((result.path / "receipt.json").read_text())["summary"]
    assert summary["acquisitions"] == 3
    assert summary["called_cells"] == 2
    assert summary["specific_cells_by_tissue"] == {
        "pancreas": 1, "pancreatic_lymph_node": 1}
    assert summary["broad_cells"] == 0
    assert summary["unresolved_cells"] == 0
    assert summary["uncalled_cells"] == 0
    assert summary["unavailable_cells"] == (2 if missing_before_ingest else 4)
    assert summary["unavailable"] == ["xenium:0040002:PanT"]
    assert summary["unknown_cell_count"] == (["xenium:0040002:PanT"]
                                              if missing_before_ingest else [])
    assert summary["paired_cell_links"] == 0
    with (result.path / "rna_acquisitions.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["tissue"] for row in rows} == {"pancreas", "pancreatic_lymph_node"}
    assert (result.path / "gate/rna_typing_scope.png").is_file()
    assert (result.path / "gate/rna_typing_classes.png").is_file()
    assert (result.path / "gate/rna_typing_support.png").is_file()
    with (result.path / "rna_classes.csv").open(newline="") as handle:
        evidence = [row for row in csv.DictReader(handle) if row["lineage"] == "Immune_T"]
    assert len(evidence) == 2
    assert all(float(row["cluster_agree_fraction"]) == 1.0 for row in evidence)
    assert all(float(row["naming_gene_fraction"]) == 1.0 for row in evidence)


def test_unresolved_t_subsets_and_broad_nonimmune_are_not_counted_as_specific_classes():
    assert pilot._call_state("T_NK_unresolved") == "unresolved"
    assert pilot._call_state("T_subset_unresolved") == "unresolved"
    assert pilot._call_state("Non_immune") == "broad"
    assert pilot._call_state("Uncalled_ambiguous") == "uncalled"
    assert pilot._call_state("Unavailable_no_nucleus") == "unavailable"
    assert pilot._call_state("T_CD8") == "class"


def test_low_rna_signal_keeps_cells_but_does_not_issue_classes(tmp_path, monkeypatch):
    monkeypatch.setattr(tasks.roots, "tmp_dir", lambda: None)
    project = _project(tmp_path)
    tasks.run_stage(rna_ingest, project)
    project.settings["rna_min_median_detected_genes"] = 20

    def fail_if_called(*args, **kwargs):
        pytest.fail("low-RNA section must not enter a typing worker")

    monkeypatch.setattr(pilot.panc_worker, "process_section", fail_if_called)
    monkeypatch.setattr(pilot.pln_worker, "process_exact_section", fail_if_called)
    result = tasks.run_stage(typing, project)
    assert result.gate_passed
    summary = json.loads((result.path / "receipt.json").read_text())["summary"]
    assert summary["called_cells"] == 0
    assert summary["candidate_tissues"] == []
    assert summary["unavailable_cells"] == summary["cells"] == 6
    assert summary["low_signal_regions"] == project.settings["rna_acquisitions"][:2]
    with (result.path / "rna_acquisitions.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert all(row["status"] == "unavailable" for row in rows)
    assert all("median" in row["reason"] for row in rows[:2])


def test_rna_gate_refuses_processed_tissue_with_no_specific_calls():
    task = SimpleNamespace(settings={"rna_acquisitions": ["xenium:panc", "xenium:pln"]})
    summary = {"acquisitions": 2, "cells": 4, "called_cells": 1, "broad_cells": 0,
               "unresolved_cells": 0, "uncalled_cells": 3, "unavailable_cells": 0,
               "unavailable": [], "unknown_cell_count": [], "low_signal_regions": [],
               "classes": 1, "paired_cell_links": 0,
               "candidate_tissues": ["pancreas", "pancreatic_lymph_node"],
               "specific_cells_by_tissue": {"pancreas": 1, "pancreatic_lymph_node": 0}}
    result = pilot.gate(task, summary)
    assert not result.passed
    assert [check.name for check in result.failing] == [
        "specific cells in pancreatic_lymph_node"]

    # A recorded damaged region has no processed candidate and cannot stop
    # an intact tissue's pilot merely because its cells were unavailable.
    summary["candidate_tissues"] = ["pancreas"]
    summary["unavailable"] = ["xenium:damaged:pln"]
    summary["uncalled_cells"] = 0
    summary["unavailable_cells"] = 3
    assert pilot.gate(task, summary).passed


def test_low_signal_rule_refuses_missing_qc_measurement():
    with pytest.raises(pilot.IO60RNATypingError, match="median detected genes"):
        pilot._low_signal({"acquisition_id": "xenium:test:region", "condition": "",
                           "median_detected_genes": "nan"}, 20)


def test_donor_metadata_uses_recorded_alias_and_npod_fallback(tmp_path):
    metadata = tmp_path / "donors.xlsx"
    npod = tmp_path / "npod.xlsx"
    pd.DataFrame({"donor_id": ["112"], "sex": ["F"],
                  "disease.status": ["ND"]}).to_excel(metadata, index=False)
    pd.DataFrame({"case_id": ["6512", "6618"], "sex": ["Female", "Male"],
                  "donor_type": ["Autoantibody Positive", "Type 1 Diabetes"]}).to_excel(
                      npod, index=False)
    task = SimpleNamespace(project=SimpleNamespace(cohort="io60_v3"),
                           input=lambda name: {"donor_metadata": metadata,
                                               "npod_cases": npod}[name])
    donors = pilot._donors(task)
    assert donors["0112"] == ("F", "ND")
    assert donors["6512"] == ("F", "Aab+")
    assert donors["6618"] == ("M", "T1D")
