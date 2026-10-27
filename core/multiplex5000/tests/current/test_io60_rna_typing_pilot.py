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
                     "sha256": tasks.sha256_file(metadata)}]
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
    return SimpleNamespace(release_dir=release, release="pilot", path=tmp_path / "pilot.yaml",
                           settings={"typing_route": "io60_xenium_pilot",
                                     "rna_acquisitions": selected},
                           source=lambda name: metadata if name == "donor_metadata" else None,
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
                      "lineage": ["Immune_T", "Unavailable_no_nucleus"]}).to_parquet(
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
    assert summary["unavailable_cells"] == (2 if missing_before_ingest else 4)
    assert summary["unavailable"] == ["xenium:0040002:PanT"]
    assert summary["unknown_cell_count"] == (["xenium:0040002:PanT"]
                                              if missing_before_ingest else [])
    assert summary["paired_cell_links"] == 0
    with (result.path / "rna_acquisitions.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["tissue"] for row in rows} == {"pancreas", "pancreatic_lymph_node"}
    assert (result.path / "gate/rna_typing_scope.png").is_file()
