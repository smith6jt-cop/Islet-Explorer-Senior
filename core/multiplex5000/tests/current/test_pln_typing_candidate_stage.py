"""A pLN stage publishes candidate status with every measured cell retained."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from multiplex5000 import gates, tasks
from multiplex5000.stages import pln_typing_pilot


def _accepted(folder: Path, stage: str) -> None:
    folder.mkdir(parents=True)
    (folder / "receipt.json").write_text(json.dumps({"stage": stage, "status": "complete"}))
    gate = folder / "gate"
    gate.mkdir()
    (gate / "gate.json").write_text(json.dumps({"stage": stage, "release": "pilot",
                                               "checks": [], "passed": True}))
    gates.accept(folder, by="test")


def test_pln_stage_retains_unavailable_and_writes_candidate_figure(tmp_path):
    release = tmp_path / "release"
    _accepted(release / "03_protein_measurement", "03_protein_measurement")
    _accepted(release / "05_normalization", "05_normalization")
    entry = "section-1"
    count = 1200
    section_path = release / "03_protein_measurement/sections.csv"
    with section_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["protein_entry", "tissue", "cells",
                                                   "eligible"])
        writer.writeheader()
        writer.writerow({"protein_entry": entry, "tissue": "pancreatic_lymph_node",
                         "cells": count, "eligible": count - 1})
    cells = release / "03_protein_measurement/cells"
    cells.mkdir()
    eligible = np.ones(count, dtype=bool)
    eligible[-1] = False
    pq.write_table(pa.table({
        "cell_uid": [f"protein:section-1:{n}" for n in range(count)],
        "object_id": [str(n) for n in range(count)],
        "protein_eligible": eligible,
        "marker:CD20": np.r_[np.linspace(5, 6, count - 20), np.full(20, 100.0)],
    }), cells / f"{entry}.parquet")
    with (release / "05_normalization/scales.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["protein_entry", "marker", "scale"])
        writer.writeheader()
        writer.writerow({"protein_entry": entry, "marker": "CD20", "scale": 1.0})
    project = SimpleNamespace(release_dir=release, release="pilot", path=tmp_path / "p.yaml",
                              settings={"pln_export_images": [entry]},
                              previous=lambda stage: "05_normalization", automatic_gates=())
    result = tasks.run_stage(pln_typing_pilot, project)
    assert result.gate_passed
    table = pq.read_table(result.path / f"protein/{entry}.parquet")
    assert len(table) == count
    assert table["candidate_state"][-1].as_py() == "unavailable"
    assert table["candidate_phenotype"][-1].as_py() == ""
    with (result.path / "protein_candidate_status.csv").open(newline="") as handle:
        assert sum(int(row["cells"]) for row in csv.DictReader(handle)) == count
    assert (result.path / "gate/pln_candidate_status.png").is_file()
