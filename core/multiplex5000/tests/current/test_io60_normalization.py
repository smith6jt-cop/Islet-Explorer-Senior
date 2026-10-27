"""The IO60 stage reads recorded block relations and never aligns another block's cells."""

from __future__ import annotations

import json
from types import SimpleNamespace

import pandas as pd
from shapely.geometry import box, mapping

from multiplex5000 import gates, normalization_block as block, normalize, tasks
from multiplex5000.stages import io60_normalization as stage
from multiplex5000.stages import normalization
from multiplex5000.tasks import Input, Task


def test_different_block_cells_stay_present_but_have_no_aligned_marker(tmp_path):
    image = "mixed-node-section"
    same, different = "node-same", "node-different"
    path = tmp_path / "expression.parquet"
    px = block.EXPORT_PIXEL_UM
    pd.DataFrame({
        "donor_id": ["112"] * 3, "object_id": ["tissue", "same", "different"],
        "qc_analysis_eligible": [True] * 3,
        "X_centroid": [50 * px, 150 * px, 350 * px],
        "Y_centroid": [50 * px, 150 * px, 350 * px],
        "cell_region": ["tissue", "pLN", "pLN"],
        "islet_num": [None] * 3, "CD8": [10.0, 20.0, 30.0],
    }).to_parquet(path, index=False)
    geometry = [
        {"type": "Feature", "id": name, "geometry": mapping(outline),
         "properties": {"classification": {"name": "pLN"}}}
        for name, outline in ((same, box(100, 100, 200, 200)),
                              (different, box(300, 300, 400, 400)))
    ]
    annotations = tmp_path / "annotations.geojson"
    annotations.write_text(json.dumps({"type": "FeatureCollection", "features": geometry}))
    task = Task(stage.STAGE, "candidate", tmp_path / "stage",
                {"expression:112": Input("expression:112", path),
                 "annotations:112": Input("annotations:112", annotations)})
    rows = [{"donor_id": "112", "image_id": image}]
    counts = stage._population(task, rows,
                               {image: {same: block.SAME, different: block.DIFFERENT}})
    assert counts == {"same-block": 1, "different-block": 1, "unassigned": 0}
    sidecars = task.staging / "alignment_block"
    assert stage._mask(path, "112", sidecars).tolist() == [True, True, False]
    assert stage._values({"112": path}, {"112": stage._mask(path, "112", sidecars)},
                         "CD8")["112"].tolist() == [10.0, 20.0]
    stage._aligned(task, {"112": path}, ["CD8"], {"CD8": {"112": 2.0}}, sidecars)
    aligned = pd.read_parquet(task.staging / "data/aligned/donor_id=112/data.parquet")
    assert aligned["object_id"].tolist() == ["tissue", "same", "different"]
    assert aligned["CD8"].tolist()[:2] == [20.0, 40.0]
    assert pd.isna(aligned.loc[2, "CD8"])
    assert aligned["in_alignment_block"].tolist() == [True, True, False]


def test_complete_io60_stage05_route_uses_declared_release_inputs(tmp_path, monkeypatch):
    release = tmp_path / "release"
    stage03 = release / "03_protein_measurement"
    stage03.mkdir(parents=True)
    (stage03 / "receipt.json").write_text(json.dumps({
        "stage": "03_protein_measurement", "status": "complete",
        "summary": {"donors": ["A", "B"]}}))
    gates.write(stage03, gates.Gate("03_protein_measurement", "fixture", [], []),
                release="candidate")
    gates.accept(stage03, by="test")
    for donor, values in (("A", [10.0, 20.0]), ("B", [12.0, 22.0])):
        path = stage03 / f"data/expression/donor_id={donor}/data.parquet"
        path.parent.mkdir(parents=True)
        pd.DataFrame({
            "donor_id": [donor] * 2, "object_id": [f"{donor}1", f"{donor}2"],
            "qc_analysis_eligible": [True, True], "X_centroid": [1.0, 2.0],
            "Y_centroid": [1.0, 2.0], "cell_region": ["tissue", "tissue"],
            "islet_num": [None, None], "CD8": values,
        }).to_parquet(path, index=False)
    spec = tmp_path / "cohort.json"
    spec.write_text(json.dumps({"images": [
        {"donor_id": donor, "image_id": donor} for donor in ("A", "B")]}))
    blocks = tmp_path / "blocks.csv"
    blocks.write_text("image,object_id,block_from_notes\n")
    annotations = tmp_path / "annotations"
    annotations.mkdir()
    references = tmp_path / "references.json"
    references.write_text(json.dumps({
        "donors": ["A", "B"], "markers": ["CD8"],
        "shifts": {"CD8": {"reference_donor": "A"}}}))
    usability = tmp_path / "usability.parquet"
    pd.DataFrame({"donor_id": ["A", "B"], "marker": ["CD8", "CD8"],
                  "acq_ratio": [0.5, 0.5]}).to_parquet(usability, index=False)
    checkout = tmp_path / "uniform"
    (checkout / "UniFORM").mkdir(parents=True)
    (checkout / "UniFORM/__init__.py").write_text("")
    sources = {"cohort_spec": spec, "pln_blocks": blocks,
               "pancreas_annotations": annotations,
               "normalization_references": references,
               "normalization_usability": usability, "uniform_checkout": checkout}
    project = SimpleNamespace(
        cohort="io60_v3", release_dir=release, release="candidate", settings={},
        path=tmp_path / "project.yaml", previous=lambda _stage: "03_protein_measurement",
        automatic_gates=(stage.STAGE,), source=lambda name: sources.get(name))

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
    assert result.gate_passed and result.accepted_automatically
    summary = json.loads((result.path / "receipt.json").read_text())["summary"]
    assert summary["CD8_shifted_donors"] == 2
    assert summary["markers_assessed"] == 1
    assert len(pd.read_parquet(result.path / "data/aligned/donor_id=A/data.parquet")) == 2
    assert (result.path / "gate/markers/CD8.png").is_file()
    assert (result.path / "gate/accepted.json").is_file()
