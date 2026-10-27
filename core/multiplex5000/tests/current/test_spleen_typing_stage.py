"""A bounded Stage 06 pilot retains separate assay identities and review evidence."""

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from scipy import sparse

from multiplex5000 import gates, tasks
from multiplex5000.stages import spleen_typing_pilot as pilot


def _accepted(folder: Path, stage: str) -> None:
    folder.mkdir(parents=True)
    (folder / "receipt.json").write_text(json.dumps({"stage": stage, "status": "complete"}))
    gate = folder / "gate"
    gate.mkdir()
    (gate / "gate.json").write_text(json.dumps({"stage": stage, "release": "pilot",
                                               "checks": [], "passed": True}))
    gates.accept(folder, by="test")


def test_stage06_writes_candidate_evidence_without_joining_cells(tmp_path, monkeypatch):
    release = tmp_path / "release"
    for stage in ("03_protein_measurement", "04_rna_ingest", "05_normalization"):
        _accepted(release / stage, stage)
    entry, acquisition, sample = "HDL090_PC67_0110618-S1", "xenium:0110618:090", "0110618__090"
    manifest = tmp_path / "pairs.yaml"
    manifest.write_text("schema: test\n")
    raw = tmp_path / "raw-output"
    pair = SimpleNamespace(protein_entry=entry, xenium_acquisition_id=acquisition,
                           donor_id="HDL090", xenium_output=raw)
    monkeypatch.setattr(pilot, "load_pairs", lambda _: (pair,))
    (release / "03_protein_measurement/receipt.json").write_text(json.dumps({
        "stage": "03_protein_measurement", "status": "complete",
        "inputs": [{"name": "same_slide_source_manifest",
                    "sha256": tasks.sha256_file(manifest)}]}))
    gates.accept(release / "03_protein_measurement", by="test")
    cells = 1220
    protein = release / "03_protein_measurement/cells"
    protein.mkdir()
    eligible = np.ones(cells, dtype=bool)
    eligible[-1] = False
    values = np.r_[np.linspace(5, 6, 1200), np.full(20, 100.0)]
    pq.write_table(pa.table({"cell_uid": [f"protein:{n}" for n in range(cells)],
                             "object_id": [str(n) for n in range(cells)],
                             "protein_eligible": eligible,
                             "marker:CD8": values}), protein / f"{entry}.parquet")
    raw.mkdir()
    for name in ("cell_feature_matrix.h5", "cells.parquet"):
        (raw / name).write_bytes(name.encode())
    receipt = release / "04_rna_ingest/receipt.json"
    receipt.write_text(json.dumps({"stage": "04_rna_ingest", "status": "complete",
                                   "inputs": [{"name": f"raw:{acquisition}:{name}",
                                               "sha256": tasks.sha256_file(raw / name)}
                                              for name in ("cell_feature_matrix.h5",
                                                           "cells.parquet")]}))
    # Acceptance binds the rewritten Stage 04 receipt.
    gates.accept(release / "04_rna_ingest", by="test")
    with (release / "04_rna_ingest/regions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["acquisition_id", "sample_key", "donor_id",
                                                   "tissue", "roi", "source_path", "condition"])
        writer.writeheader()
        writer.writerow({"acquisition_id": acquisition, "sample_key": sample,
                         "donor_id": "HDL090", "tissue": "spleen", "roi": "090",
                         "source_path": str(raw), "condition": ""})
    rna_folder = release / "04_rna_ingest/cells"
    rna_folder.mkdir()
    rna_uid = [f"{acquisition}:cells:{n}" for n in range(cells)]
    pq.write_table(pa.table({"cell_uid": rna_uid, "typing_eligible": eligible}),
                   rna_folder / f"{sample}.parquet")
    with (release / "05_normalization/scales.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["protein_entry", "marker", "scale"])
        writer.writeheader()
        writer.writerow({"protein_entry": entry, "marker": "CD8", "scale": 1.0})
    raw_counts = sparse.csr_matrix(np.column_stack([
        np.ones(cells, dtype=np.int64), np.full(cells, 9, dtype=np.int64)]))
    monkeypatch.setattr(pilot, "load_accepted_region", lambda *args: SimpleNamespace(
        counts=raw_counts, genes=("CD3E", "OTHER"), cell_uid=np.array(rna_uid),
        eligible=eligible, n_cells=cells))
    project = SimpleNamespace(release_dir=release, release="pilot", path=tmp_path / "p.yaml",
                              settings={"same_slide_source_manifest": str(manifest),
                                        "rna_acquisitions": [acquisition]},
                              source=lambda name: manifest if name == "same_slide_sources" else None,
                              previous=lambda stage: "05_normalization", automatic_gates=())
    result = tasks.run_stage(pilot, project)
    assert result.gate_passed
    assert (result.path / f"protein/{entry}.parquet").is_file()
    assert (result.path / f"rna/{sample}.parquet").is_file()
    assert not any("paired_cell" in path.name or "final_class" in path.name
                   for path in result.path.rglob("*"))
    protein_calls = pq.read_table(result.path / f"protein/{entry}.parquet")
    rna_calls = pq.read_table(result.path / f"rna/{sample}.parquet")
    assert protein_calls["marker_call:CD8"][-1].as_py() == -1
    assert protein_calls["rule:Immune"][-1].as_py() == -1
    assert protein_calls["candidate_state"][-1].as_py() == "unavailable"
    assert protein_calls["candidate_phenotype"][-1].as_py() == ""
    assert rna_calls["gene:CD3E"][-1].as_py() == -1
    assert protein_calls["cell_uid"][0].as_py() != rna_calls["cell_uid"][0].as_py()
    assert (result.path / "gate/assay_evidence.png").is_file()
    assert (result.path / "gate/spleen_candidate_status.png").is_file()


def test_stage06_refuses_a_different_xenium_output_for_the_same_donor(tmp_path, monkeypatch):
    manifest = tmp_path / "pairs.yaml"
    manifest.write_text("schema: test\n")
    expected = tmp_path / "correct-output"
    pair = SimpleNamespace(protein_entry="HDL090-PC67", donor_id="HDL090",
                           xenium_acquisition_id="xenium:0110618:090",
                           xenium_output=expected)
    monkeypatch.setattr(pilot, "load_pairs", lambda _: (pair,))
    stage04 = tmp_path / "release/04_rna_ingest"
    stage04.mkdir(parents=True)
    with (stage04 / "regions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["acquisition_id", "donor_id",
                                                   "tissue", "source_path"])
        writer.writeheader()
        writer.writerow({"acquisition_id": pair.xenium_acquisition_id,
                         "donor_id": pair.donor_id, "tissue": "spleen",
                         "source_path": str(tmp_path / "other-output")})
    project = SimpleNamespace(release_dir=tmp_path / "release",
                              settings={"same_slide_source_manifest": str(manifest),
                                        "rna_acquisitions": [pair.xenium_acquisition_id]},
                              source=lambda name: manifest if name == "same_slide_sources" else None)
    with pytest.raises(pilot.SpleenTypingPilotError, match="Stage 04 mismatch"):
        pilot._selected(project)
