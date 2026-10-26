"""A bounded pLN Stage 03 keeps native identity and missing measurements."""

import csv

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from multiplex5000 import gates, project, stages, tasks
from multiplex5000 import normalize


def _csv(path, columns, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _project(tmp_path):
    root = tmp_path / "exports"
    root.mkdir()
    roster, index = tmp_path / "roster.csv", tmp_path / "index.csv"
    roster_rows, index_rows = [], []
    images = ["9001pLN_Scan1.er.qptiff", "9002pLN_Scan1.er.qptiff"]
    for image, donor in zip(images, ("9001", "9002"), strict=True):
        source = root / f"{image}.parquet"
        pq.write_table(pa.table({
            "Image": [f"{image} - resolution #1"] * 3,
            "Object ID": ["native-0", "native-1", "native-2"],
            "Centroid X µm": [1.0, 2.0, 3.0],
            "Centroid Y µm": [4.0, 5.0, 6.0],
            "Nucleus: Area µm^2": [2.0, None, 0.0],
            "Cell: DAPI: Mean": [20.0, 21.0, 22.0],
            "Cell: CD8: Mean": [8.0 if donor == "9001" else 16.0, None, 0.0],
            "Cell: CD20: Mean": [None, 3.0, 0.0],
        }), source)
        roster_rows.append({"image_file": image, "donor_id": donor, "roi": "pln_1",
                            "n_cells": "2", "excluded": "False", "excluded_reason": ""})
        index_rows.append({"image_file": image, "source_path": source.name,
                           "segmentation_id": "v3-test"})
    _csv(roster, ["image_file", "donor_id", "roi", "n_cells", "excluded",
                  "excluded_reason"], roster_rows)
    _csv(index, ["image_file", "source_path", "segmentation_id"], index_rows)
    uniform = tmp_path / "uniform"
    (uniform / "UniFORM").mkdir(parents=True)
    (uniform / "UniFORM/__init__.py").write_text("")
    spec = {"schema": "multiplex5000-project/1", "id": "pln_test",
            "release": "pln-test-1", "cohort": "io60_v3",
            "stages": ["03_protein_measurement", "05_normalization", "06_typing"],
            "automatic_gates": ["03_protein_measurement", "05_normalization"],
            "settings": {"pln_export_images": images},
            "sources": {"pln_v3_export_cells": str(root),
                        "pln_export_roster": str(roster), "pln_export_index": str(index),
                        "uniform_checkout": str(uniform)},
            "releases_root": str(tmp_path / "releases")}
    path = tmp_path / "pln_test.yaml"
    path.write_text(yaml.safe_dump(spec, sort_keys=False))
    return project.load(path), images


def test_pilot_stage_03_preserves_object_ids_nucleus_and_missing_values(tmp_path):
    pilot, images = _project(tmp_path)
    result = tasks.run_stage(stages.get("03_protein_measurement"), pilot)
    assert result.gate_passed
    rows = list(csv.DictReader((result.path / "sections.csv").open(newline="")))
    assert len(rows) == 2
    assert {row["protein_entry"] for row in rows} == {"9001__pln_1", "9002__pln_1"}
    assert {(int(row["cells"]), int(row["eligible"]), int(row["no_nucleus"]),
             int(row["missing_nucleus"])) for row in rows} == {(3, 1, 1, 1)}
    table = pq.read_table(result.path / "cells/9001__pln_1.parquet")
    assert table["object_id"].to_pylist() == ["native-0", "native-1", "native-2"]
    assert table["cell_uid"].to_pylist() == [
        f"phenocycler:{images[0]}:v3-test:native-{i}" for i in range(3)]
    assert table["protein_eligible"].to_pylist() == [True, False, False]
    values = table["marker:CD20"].to_pylist()
    assert np.isnan(values[0]) and values[1:] == [3.0, 0.0]
    assert "marker:DAPI" not in table.column_names
    gate = gates.load(result.path)
    assert gate["passed"]
    assert "No REDSEA" in gate["summary"]
    assert (result.path / "gate/protein_cells.png").stat().st_size > 1000


def test_pilot_routes_stage03_stage05_stage06_through_one_project(tmp_path, monkeypatch):
    pilot, _ = _project(tmp_path)

    class Registration:
        @staticmethod
        def automatic_registration(hist, markers, _names, samples, _refs, _landmarks, _bins):
            return {markers[0]: [0] * len(samples)}, [samples[0]], None

    class Normalization:
        @staticmethod
        def calculate_shift_in_log_pixels(ranges, markers, bins, shifts):
            return {markers[0]: [0.0] * len(shifts[markers[0]])}

    monkeypatch.setattr(normalize, "load_uniform", lambda _: (Registration, Normalization))
    protein = tasks.run_stage(stages.get("03_protein_measurement"), pilot)
    assert protein.gate_passed and protein.accepted_automatically
    shifted = tasks.run_stage(stages.get("05_normalization"), pilot)
    assert shifted.gate_passed and shifted.accepted_automatically
    typed = tasks.run_stage(stages.get("06_typing"), pilot)
    assert typed.gate_passed and not typed.accepted_automatically
    receipt = tasks.read_receipt(typed.path)
    assert receipt["summary"]["protein_rules"] == 2 * 69
    assert receipt["summary"]["final_cell_classes"] == 0
    table = pq.read_table(typed.path / "protein/9001__pln_1.parquet")
    assert table["object_id"].to_pylist() == ["native-0", "native-1", "native-2"]
    assert all(table[name].to_pylist()[1:] == [-1, -1]
               for name in table.column_names if name.startswith("rule:"))
    assert table["rule:Naive_B_cell"].to_pylist()[0] == -1  # CD38 unmeasured
    assert (typed.path / "gate/pln_rule_evidence.png").stat().st_size > 1000
    assert (typed.path / "gate/pln_marker_support.png").stat().st_size > 1000
