"""The pilot's full-panel fitter uses one pinned UniFORM path and preserves gaps."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from multiplex5000 import gates, normalize, roots, tasks
from multiplex5000.stages import normalization


def _accepted(stage: Path, name: str) -> None:
    stage.mkdir(parents=True)
    (stage / "receipt.json").write_text(json.dumps({"stage": name, "status": "complete"}))
    gate_dir = stage / "gate"
    gate_dir.mkdir()
    (gate_dir / "gate.json").write_text(json.dumps({"stage": name, "release": "test",
                                                   "checks": [], "passed": True}))
    gates.accept(stage, by="test")


def _project(tmp_path: Path, monkeypatch):
    release = tmp_path / "release"
    protein = release / "03_protein_measurement"
    _accepted(protein, "03_protein_measurement")
    _accepted(release / "04_rna_ingest", "04_rna_ingest")
    entries = ("HDL090_PC67_0110618-S1", "HDL082_PC67_0110618-S2")
    with (protein / "sections.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["protein_entry"])
        writer.writeheader()
        for entry in entries:
            writer.writerow({"protein_entry": entry})
    (protein / "cells").mkdir()
    for index, entry in enumerate(entries):
        pq.write_table(pa.table({"cell_uid": [f"{entry}:one", f"{entry}:two"],
                                 "object_id": ["one", "two"],
                                 "protein_eligible": [True, False],
                                 "marker:CD8": [20.0 * (index + 1), 3.0],
                                 "marker:CD20": [100.0 * (index + 1), None]}),
                       protein / "cells" / f"{entry}.parquet")
    checkout = tmp_path / "uniform"
    (checkout / "UniFORM").mkdir(parents=True)
    (checkout / "UniFORM/__init__.py").write_text("")
    pairs = tmp_path / "pairs.yaml"
    pairs.write_text("schema: fixture\n")
    from multiplex5000.same_slide_sources import SameSlidePair

    recorded = [SameSlidePair("spleen", entry.split("_")[0], entry,
                              tmp_path / f"{entry}.csv.gz", tmp_path / "image.qptiff",
                              f"xenium:0110618:0{index}", tmp_path / f"output-{index}",
                              "same_section_after_phenocycler", "crop_local_um",
                              "region_local_um") for index, entry in enumerate(entries)]
    monkeypatch.setattr(normalization, "load_pairs", lambda _: tuple(recorded))
    project = SimpleNamespace(
        release_dir=release, release="test", path=tmp_path / "pilot.yaml",
        settings={"same_slide_source_manifest": "fixture"},
        source=lambda name: {"same_slide_sources": pairs,
                             "uniform_checkout": checkout}.get(name),
        previous=lambda stage: "04_rna_ingest", automatic_gates=())
    return project


def test_stage05_fits_each_measured_marker_and_all_cd8_sections(tmp_path, monkeypatch):
    monkeypatch.setattr(roots, "tmp_dir", lambda: tmp_path / "stage_tmp")
    project = _project(tmp_path, monkeypatch)

    class Registration:
        @staticmethod
        def automatic_registration(hist, markers, _names, samples, _refs, _landmarks, _bins):
            return {markers[0]: [0] * len(samples)}, [samples[0]], None

    class Normalization:
        @staticmethod
        def calculate_shift_in_log_pixels(ranges, markers, bins, shifts):
            return {markers[0]: [0.0] * len(shifts[markers[0]])}

    monkeypatch.setattr(normalize, "load_uniform", lambda _: (Registration, Normalization))
    result = tasks.run_stage(normalization, project)
    assert result.gate_passed
    receipt = json.loads((result.path / "receipt.json").read_text())
    assert receipt["summary"]["markers_fitted"] == 2
    assert receipt["summary"]["CD8_shifted_sections"] == 2
    with (result.path / "scales.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {(row["protein_entry"], row["marker"]) for row in rows} == {
        (entry, marker) for entry in ("HDL090_PC67_0110618-S1", "HDL082_PC67_0110618-S2")
        for marker in ("CD8", "CD20")}
    assert (result.path / "gate/mode_spread.png").is_file()


def test_histogram_refuses_an_absent_marker_instead_of_inserting_zero():
    with pytest.raises(normalize.NormalizeError, match="no finite marker values"):
        normalize.marker_histograms({"A": np.array([None, np.nan], dtype=float),
                                     "B": np.array([5.0, 8.0])})
