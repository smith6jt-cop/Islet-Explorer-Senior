"""Stage 04 ingests both tissues from their accepted Xenium regions."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import h5py
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from multiplex5000 import gates, tasks
from multiplex5000.artifacts import sha256_file
from multiplex5000.stages import rna_ingest
from multiplex5000.xenium_sources import CELL_FILES, INSTRUMENT_FILES
from multiplex5000.xenium_typing_input import (
    XeniumTypingInputError,
    load_accepted_region,
    write_worker_h5ad,
)


def _bundle(root: Path, serial: str, roi: str, nuclei: list[int | None]) -> Path:
    path = root / f"output-XETG00298__{serial}__{roi}__20261026__120000"
    path.mkdir(parents=True)
    for name in (*INSTRUMENT_FILES, "cell_boundaries.parquet", "nucleus_boundaries.parquet"):
        (path / name).write_text(f"{serial} {roi} {name}")
    with h5py.File(path / "cell_feature_matrix.h5", "w") as file:
        matrix = file.create_group("matrix")
        matrix.create_dataset("shape", data=[2, 2])
        matrix.create_dataset("barcodes", data=[b"shared-1", b"other-1"])
        features = matrix.create_group("features")
        features.create_dataset("name", data=[b"CD3E", b"TRDC"])
        features.create_dataset("feature_type", data=[b"Gene Expression", b"Gene Expression"])
        matrix.create_dataset("data", data=[2, 1, 3])
        matrix.create_dataset("indices", data=[0, 1, 1])
        matrix.create_dataset("indptr", data=[0, 2, 3])
    pq.write_table(pa.table({"cell_id": ["shared-1", "other-1"],
                             "x_centroid": [10.0, 20.0], "y_centroid": [30.0, 40.0],
                             "nucleus_count": nuclei}), path / "cells.parquet")
    assert all((path / name).is_file() for name in CELL_FILES)
    return path


def _project(tmp_path: Path, sources: list[Path]):
    release = tmp_path / "release"
    stage = release / "00_inventory"
    stage.mkdir(parents=True)
    index = tmp_path / "site_index.csv"
    fields = ("acquisition_id", "output", "source_path")
    with index.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for source, roi in zip(sources, ("Panc", "PanT"), strict=True):
            writer.writerow({"acquisition_id": f"xenium:0040001:{roi}",
                             "output": source.name, "source_path": str(source)})
    table = stage / "acquisitions.csv"
    columns = ("platform", "acquisition_id", "donor_id", "roi", "tissue", "output", "path",
               "condition")
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for source, roi in zip(sources, ("Panc", "PanT"), strict=True):
            writer.writerow({"platform": "xenium", "acquisition_id": f"xenium:0040001:{roi}",
                             "donor_id": "6539", "roi": roi,
                             "tissue": "pancreas" if roi == "Panc" else "pancreatic_lymph_node",
                             "output": source.name, "path": str(source),
                             "condition": "damaged" if roi == "PanT" else ""})
    input_hashes = [{"name": f"xenium:{source.name}:{name}",
                     "sha256": sha256_file(source / name)}
                    for source in sources for name in INSTRUMENT_FILES]
    (stage / "receipt.json").write_text(json.dumps({
        "stage": "00_inventory", "status": "complete", "inputs": input_hashes,
        "outputs": [{"path": "acquisitions.csv", "sha256": sha256_file(table)}]}))
    gate_dir = stage / "gate"
    gate_dir.mkdir()
    (gate_dir / "gate.json").write_text(json.dumps({"stage": "00_inventory", "release": "test",
                                                   "checks": [], "passed": True}))
    gates.accept(stage, by="test")
    return SimpleNamespace(release_dir=release, release="test", path=tmp_path / "project.yaml",
                           source=lambda name: index if name == "xenium_source_index" else None,
                           previous=lambda stage_id: "00_inventory",
                           automatic_gates=("04_rna_ingest",))


def test_stage_04_same_native_ids_and_damage_do_not_lose_other_region(tmp_path, monkeypatch):
    monkeypatch.setattr(tasks.roots, "tmp_dir", lambda: None)
    a = _bundle(tmp_path, "0040001", "Panc", [1, 0])
    b = _bundle(tmp_path, "0040001", "PanT", [1, None])
    project = _project(tmp_path, [a, b])
    result = tasks.run_stage(rna_ingest, project)
    assert result.gate_passed and result.accepted_automatically
    assert result.path == project.release_dir / "04_rna_ingest"
    receipt = json.loads((result.path / "receipt.json").read_text())
    assert receipt["summary"]["regions"] == 2
    assert receipt["summary"]["cells"] == 4
    assert receipt["summary"]["typing_eligible"] == 2
    assert receipt["summary"]["no_nucleus"] == 1
    assert receipt["summary"]["nucleus_missing"] == 1
    assert receipt["summary"]["damaged_regions"] == ["xenium:0040001:PanT"]
    assert any(item["name"] == "xenium_source_index" and item["sha256"]
               for item in receipt["inputs"])
    assert sum(item["name"].endswith("cell_feature_matrix.h5") for item in receipt["inputs"]) == 2
    panc = pq.read_table(result.path / "cells/0040001__Panc.parquet").to_pydict()
    pln = pq.read_table(result.path / "cells/0040001__PanT.parquet").to_pydict()
    assert panc["cell_id"][0] == pln["cell_id"][0]
    assert panc["cell_uid"][0] != pln["cell_uid"][0]
    assert panc["nucleus_count"] == [1, 0]
    assert pln["nucleus_count"] == [1, None]
    assert panc["typing_eligible"] == [True, False]
    assert pln["typing_eligible"] == [True, False]


def test_typing_input_uses_exact_stage04_acquisition_and_raw_counts(tmp_path, monkeypatch):
    import anndata as ad

    monkeypatch.setattr(tasks.roots, "tmp_dir", lambda: None)
    a = _bundle(tmp_path, "0040001", "Panc", [1, 0])
    b = _bundle(tmp_path, "0040001", "PanT", [1, None])
    project = _project(tmp_path, [a, b])
    stage = tasks.run_stage(rna_ingest, project).path
    with (stage / "regions.csv").open(newline="") as handle:
        rows = {row["roi"]: row for row in csv.DictReader(handle)}
    row = rows["PanT"]
    region = load_accepted_region(
        row, stage / "receipt.json", stage / "regions.csv",
        stage / f"cells/{row['sample_key']}.parquet",
        b / "cell_feature_matrix.h5", b / "cells.parquet")
    assert region.cell_uid[0].startswith("xenium:0040001:PanT:")
    assert region.nucleus_missing.tolist() == [False, True]
    assert region.eligible.tolist() == [True, False]
    h5ad = tmp_path / "worker-input.h5ad"
    write_worker_h5ad(region, h5ad, sex="F", disease_status="T1D", batch="0040001")
    worker = ad.read_h5ad(h5ad)
    assert worker.obs_names.tolist() == region.cell_id.tolist()
    assert worker.var_names.tolist() == list(region.genes)
    assert worker.X.toarray().tolist() == region.counts.toarray().tolist()
    assert worker.obs["disease.status"].tolist() == ["T1D", "T1D"]
    with pytest.raises(XeniumTypingInputError, match="raw inputs point elsewhere"):
        load_accepted_region(row, stage / "receipt.json", stage / "regions.csv",
                             stage / f"cells/{row['sample_key']}.parquet",
                             a / "cell_feature_matrix.h5", a / "cells.parquet")
    with (b / "cells.parquet").open("ab") as handle:
        handle.write(b"changed after Stage 04")
    with pytest.raises(XeniumTypingInputError, match="changed after ingest"):
        load_accepted_region(row, stage / "receipt.json", stage / "regions.csv",
                             stage / f"cells/{row['sample_key']}.parquet",
                             b / "cell_feature_matrix.h5", b / "cells.parquet")


def test_stage_04_refuses_missing_indexed_source(tmp_path):
    a = _bundle(tmp_path, "0040001", "Panc", [1, 1])
    b = _bundle(tmp_path, "0040001", "PanT", [1, 1])
    project = _project(tmp_path, [a, b])
    (b / "cell_feature_matrix.h5").unlink()
    with pytest.raises(rna_ingest.RNAIngestError, match="every selected source ready"):
        rna_ingest.inputs(project)


def test_selected_pilot_does_not_require_or_claim_unselected_region(tmp_path, monkeypatch):
    monkeypatch.setattr(tasks.roots, "tmp_dir", lambda: None)
    a = _bundle(tmp_path, "0040001", "Panc", [1, 1])
    b = _bundle(tmp_path, "0040001", "PanT", [1, 1])
    project = _project(tmp_path, [a, b])
    project.settings = {"rna_acquisitions": ["xenium:0040001:Panc"]}
    (b / "cell_feature_matrix.h5").unlink()
    result = tasks.run_stage(rna_ingest, project)
    assert result.gate_passed
    receipt = json.loads((result.path / "receipt.json").read_text())
    assert receipt["summary"]["regions"] == 1
    assert receipt["summary"]["accepted_regions"] == 2
    assert receipt["summary"]["scope"] == "selected"
    assert receipt["settings"]["rna_acquisitions"] == ["xenium:0040001:Panc"]
    assert not (result.path / "cells/0040001__PanT.parquet").exists()
    assert all("PanT" not in item["name"] for item in receipt["inputs"])


def test_local_pilot_uses_accepted_paths_and_records_manifest_without_site_index(
        tmp_path, monkeypatch):
    monkeypatch.setattr(tasks.roots, "tmp_dir", lambda: None)
    a = _bundle(tmp_path, "0040001", "Panc", [1, 1])
    b = _bundle(tmp_path, "0040001", "PanT", [1, 1])
    project = _project(tmp_path, [a, b])
    project.settings = {"rna_acquisitions": ["xenium:0040001:Panc"]}
    records = tmp_path / "records"
    records.mkdir()
    manifest = records / "MANIFEST.tsv"
    manifest.write_text("output\tfile\tsource\tsha256\tbytes\tcopied_utc\n")
    stage = project.release_dir / "00_inventory"
    receipt = json.loads((stage / "receipt.json").read_text())
    receipt["inputs"].append({"name": "xenium_records_manifest",
                              "sha256": sha256_file(manifest)})
    (stage / "receipt.json").write_text(json.dumps(receipt))
    gates.accept(stage, by="test")
    project.source = lambda name: records if name == "xenium_records" else None
    (b / "cell_feature_matrix.h5").unlink()
    result = tasks.run_stage(rna_ingest, project)
    assert result.gate_passed
    bound = json.loads((result.path / "receipt.json").read_text())["inputs"]
    assert any(item["name"] == "xenium_records_manifest" for item in bound)
    assert not any(item["name"] == "xenium_source_index" for item in bound)
