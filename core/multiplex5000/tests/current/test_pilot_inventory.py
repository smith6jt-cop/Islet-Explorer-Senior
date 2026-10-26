"""A same-slide source manifest enters the existing Stage 00/04 contract."""

import csv
import gzip
import hashlib
import json

import pyarrow.parquet as pq
import pytest
import yaml

from multiplex5000 import project, roots, tasks
from multiplex5000.stages import inventory, measurement, pilot_protein
from multiplex5000.xenium_sources import CELL_FILES, INSTRUMENT_FILES, build_roster


def test_spleen_pilot_project_uses_both_workflow_routes_stage_functions():
    pilot = project.load("spleen_pilot")
    assert pilot.stages == ("00_inventory", "03_protein_measurement", "04_rna_ingest")
    assert "protein_entries" not in pilot.settings
    assert pilot.settings["rna_acquisitions"] == ["xenium:0110618:090"]


def _source_project(tmp_path):
    section = tmp_path / "sections.csv"
    with section.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["slide", "section", "bbox_x", "bbox_y", "donor_id", "entry_name"])
        writer.writerow(["0110618", "S1", "720", "416", "HDL090", "HDL090_PC67_0110618-S1"])
    header = "object_id\tx_um\ty_um\tNucleus: Area µm^2\tCell: CD20: Mean\n".encode()
    export = tmp_path / "HDL090_PC67_0110618-S1.csv.gz"
    with gzip.open(export, "wb") as handle:
        handle.write(header)
        handle.write(b"shared-1\t10\t20\t7\t5\n")
        handle.write(b"zero-1\t11\t21\t0\t2\n")
        handle.write(b"missing-1\t12\t22\t\t\n")
    image = tmp_path / "0110618_Scan1.er.qptiff"
    image.write_bytes(b"image")
    output = tmp_path / "output-XETG00298__0110618__090__20260826__224731"
    output.mkdir()
    for name in (*INSTRUMENT_FILES, *CELL_FILES):
        (output / name).write_bytes(b"source")
    manifest = tmp_path / "pairs.yaml"
    manifest.write_text(yaml.safe_dump({
        "schema": "multiplex5000-same-slide-sources/1", "tissue": "spleen",
        "relationship": "same_section_after_phenocycler",
        "protein_frame": "crop_local_um", "xenium_frame": "region_local_um",
        "expected_pairs": 1,
        "section_table": {"path": str(section),
                          "sha256": hashlib.sha256(section.read_bytes()).hexdigest()},
        "protein_header_sha256": hashlib.sha256(header).hexdigest(),
        "pairs": [{"donor_id": "HDL090", "protein_entry": "HDL090_PC67_0110618-S1",
                   "protein_export": str(export), "protein_image": str(image),
                   "xenium_acquisition_id": "xenium:0110618:090",
                   "xenium_output": str(output)}],
    }))
    return project.Project(tmp_path / "spleen_pilot.yaml", "spleen_pilot", "pilot-test",
                           "spleen_same_slide_pilot",
                           ("00_inventory", "03_protein_measurement", "04_rna_ingest"),
                           ("00_inventory", "03_protein_measurement"), tmp_path / "releases",
                           {"same_slide_source_manifest": "unused-tracked-default",
                            "protein_entries": ["HDL090_PC67_0110618-S1"]},
                           {"same_slide_sources": manifest})


def test_stage00_publishes_exact_pairs_and_stage04_resolves_them(tmp_path, monkeypatch):
    monkeypatch.setattr(roots, "tmp_dir", lambda: tmp_path / "stage_tmp")
    built = _source_project(tmp_path)
    result = tasks.run_stage(inventory, built)
    stage = built.release_dir / "00_inventory"
    assert result.status == "published" and result.gate_passed
    assert (stage / "gate/accepted.json").is_file()
    with (stage / "acquisitions.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["platform"] for row in rows} == {"phenocycler", "xenium"}
    assert {row["donor_id"] for row in rows} == {"HDL090"}
    assert rows[1]["acquisition_id"] == "xenium:0110618:090"
    with (stage / "relationships.csv").open(newline="") as handle:
        relationship, = csv.DictReader(handle)
    assert relationship["relationship"] == "same_section_after_phenocycler"
    assert relationship["crop_origin_x_px"] == "720"
    receipt = json.loads((stage / "receipt.json").read_text())
    assert len([item for item in receipt["inputs"] if item["name"].startswith("xenium:")
                and item.get("sha256")]) == len(INSTRUMENT_FILES)
    source, = build_roster(stage)
    assert source.status == "ready"
    assert source.source_path == rows[1]["path"]


def test_stage03_preserves_objects_and_excludes_cells_without_a_nucleus(tmp_path, monkeypatch):
    monkeypatch.setattr(roots, "tmp_dir", lambda: tmp_path / "stage_tmp")
    built = _source_project(tmp_path)
    tasks.run_stage(inventory, built)
    result = tasks.run_stage(measurement, built)
    assert result.gate_passed and result.accepted_automatically
    table = pq.read_table(result.path / "cells/HDL090_PC67_0110618-S1.parquet").to_pydict()
    assert table["object_id"] == ["shared-1", "zero-1", "missing-1"]
    assert len(set(table["cell_uid"])) == 3
    assert table["protein_eligible"] == [True, False, False]
    assert table["marker:CD20"][:2] == [5.0, 2.0]
    assert table["nucleus_area_um2"][:2] == [7.0, 0.0]
    summary = json.loads((result.path / "receipt.json").read_text())["summary"]
    assert (summary["cells"], summary["eligible"], summary["no_nucleus"],
            summary["missing_nucleus"]) == (3, 1, 1, 1)
    assert summary["scope"] == "selected"


def test_stage03_refuses_an_export_changed_after_stage00(tmp_path, monkeypatch):
    monkeypatch.setattr(roots, "tmp_dir", lambda: tmp_path / "stage_tmp")
    built = _source_project(tmp_path)
    tasks.run_stage(inventory, built)
    export = tmp_path / "HDL090_PC67_0110618-S1.csv.gz"
    with gzip.open(export, "ab") as handle:
        handle.write(b"extra-1\t15\t25\t8\t4\n")
    with pytest.raises(pilot_protein.PilotProteinError, match="differs from accepted Stage 00"):
        tasks.run_stage(measurement, built)
