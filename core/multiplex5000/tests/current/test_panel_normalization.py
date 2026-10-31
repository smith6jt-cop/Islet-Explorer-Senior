"""Per-panel Normalization: UniFORM never sees two panels at once, and missing stays missing."""

from __future__ import annotations

import csv
import json
import math
from types import SimpleNamespace

import numpy as np
import pandas as pd

from multiplex5000 import gates, normalize, tasks
from multiplex5000.stages import normalization
from multiplex5000.stages import panel_normalization as stage

PANELS = {"pc67": ["A", "B"], "pc33": ["C", "D"]}


def _project(tmp_path, *, pinned=None):
    release = tmp_path / "release"
    measured = release / "03_protein_measurement"
    measured.mkdir(parents=True)
    donors = [d for ds in PANELS.values() for d in ds]
    (measured / "receipt.json").write_text(json.dumps({
        "stage": "03_protein_measurement", "status": "complete",
        "summary": {"donors": sorted(donors)}}))
    gates.write(measured, gates.Gate("03_protein_measurement", "fixture", [], []),
                release="candidate")
    gates.accept(measured, by="test")
    rng = np.random.default_rng(0)
    for panel, ds in PANELS.items():
        for d in ds:
            n = 50
            frame = {"donor_id": [d] * n, "object_id": [f"{d}{i}" for i in range(n)],
                     "qc_analysis_eligible": [True] * (n - 1) + [False],
                     "qc_estimation_eligible": [True] * n, "image": [d] * n,
                     "X_centroid": rng.uniform(0, 1e4, n), "Y_centroid": rng.uniform(0, 1e4, n),
                     "cell_region": ["tissue"] * n, "islet_num": [None] * n,
                     "CD8": rng.uniform(5, 50, n)}
            # only the 67-plex panel acquires CD20; only section A has a measurable SOX2
            frame["CD20"] = rng.uniform(5, 50, n) if panel == "pc67" else [np.nan] * n
            frame["SOX2"] = rng.uniform(5, 50, n) if d == "A" else [np.nan] * n
            path = measured / f"data/expression/donor_id={d}/data.parquet"
            path.parent.mkdir(parents=True)
            pd.DataFrame(frame).to_parquet(path, index=False)
    spec = tmp_path / "cohort.json"
    spec.write_text(json.dumps({"spec_version": 1, "images": [
        {"donor_id": d, "image_id": d, "panel_id": p} for p, ds in PANELS.items() for d in ds]}))
    checkout = tmp_path / "uniform"
    (checkout / "UniFORM").mkdir(parents=True)
    (checkout / "UniFORM/__init__.py").write_text("")
    sources = {"cohort_spec": spec, "uniform_checkout": checkout}
    if pinned:
        sources["normalization_references"] = pinned
    return SimpleNamespace(
        cohort="spleen_fixture", release_dir=release, release="candidate",
        settings={"normalization": "per_panel"}, path=tmp_path / "project.yaml",
        previous=lambda _stage: "03_protein_measurement",
        automatic_gates=(stage.STAGE,), source=lambda name: sources.get(name))


class _Uniform:
    """Records every sample set a fit sees; the second sample is shifted by log 2 (scale 0.5)."""

    def __init__(self):
        self.seen = []
        self.references = []

    def automatic_registration(self, hist, markers, _names, samples, refs, _landmarks, _bins):
        self.seen.append((markers[0], tuple(samples)))
        self.references.append(refs)
        return {markers[0]: [0] * len(samples)}, [refs[0] if refs else samples[0]], None

    @staticmethod
    def calculate_shift_in_log_pixels(ranges, markers, bins, shifts):
        n = len(shifts[markers[0]])
        return {markers[0]: [0.0] + [math.log(2.0)] * (n - 1)}


def _run(tmp_path, monkeypatch, **kw):
    uniform = _Uniform()
    monkeypatch.setattr(normalize, "load_uniform", lambda _: (uniform, uniform))
    result = tasks.run_stage(normalization, _project(tmp_path, **kw))
    return result, uniform


def test_the_route_is_chosen_by_the_project_setting(tmp_path):
    assert normalization._route(_project(tmp_path)) is stage
    other = SimpleNamespace(settings={}, cohort="io60_v3")
    assert normalization._route(other) is not stage


def test_each_panel_is_fitted_alone_and_the_gate_accepts_itself(tmp_path, monkeypatch):
    result, uniform = _run(tmp_path, monkeypatch)
    assert result.gate_passed and result.accepted_automatically
    for _marker, samples in uniform.seen:
        assert set(samples) <= set(PANELS["pc67"]) or set(samples) <= set(PANELS["pc33"])
    assert ("CD8", ("A", "B")) in uniform.seen and ("CD8", ("C", "D")) in uniform.seen
    summary = json.loads((result.path / "receipt.json").read_text())["summary"]
    assert summary["panels"]["pc67"]["CD8_aligned_sections"] == 2
    assert summary["panels"]["pc33"]["markers_unavailable"] == ["CD20", "SOX2"]
    assert summary["panels"]["pc67"]["markers_unavailable"] == ["SOX2"]


def test_scales_apply_and_unavailable_markers_are_nan_not_raw(tmp_path, monkeypatch):
    result, _ = _run(tmp_path, monkeypatch)
    raw = pd.read_parquet(tmp_path / "release/03_protein_measurement/data/expression/"
                          "donor_id=B/data.parquet")
    aligned = pd.read_parquet(result.path / "data/aligned/donor_id=B/data.parquet")
    assert aligned["object_id"].tolist() == raw["object_id"].tolist()   # every cell keeps its row
    assert aligned["X_centroid"].equals(raw["X_centroid"])            # metadata is never scaled
    assert np.allclose(aligned["CD8"], raw["CD8"] * 0.5)              # B is shifted by log 2
    assert aligned["SOX2"].isna().all()                                  # one measurable section
    c = pd.read_parquet(result.path / "data/aligned/donor_id=C/data.parquet")
    assert c["CD20"].isna().all()                                        # not acquired on the panel
    with (result.path / "scales.csv").open() as handle:
        rows = {(r["donor_id"], r["marker"]): r for r in csv.DictReader(handle)}
    assert rows[("A", "SOX2")]["status"] == "unavailable" and rows[("A", "SOX2")]["scale"] == ""
    assert rows[("B", "CD8")]["status"] == "aligned"
    assert rows[("C", "CD8")]["panel_id"] == "pc33"


def test_a_previous_scales_file_pins_the_references(tmp_path, monkeypatch):
    pinned = tmp_path / "pinned.csv"
    pinned.write_text("donor_id,panel_id,marker,reference_donor,scale,shift_log,n_fit_cells,"
                      "status,mode_sd_before,mode_sd_after\n"
                      "B,pc67,CD8,B,1,0,49,aligned,0.1,0.0\n")
    _, uniform = _run(tmp_path, monkeypatch, pinned=pinned)
    refs = dict(zip([m for m, _ in uniform.seen], uniform.references))
    by_fit = list(zip(uniform.seen, uniform.references))
    assert [r for (m, s), r in by_fit if m == "CD8" and s == ("A", "B")] == [["B"]]
    assert refs["CD20"] is None


def test_a_column_that_is_neither_metadata_nor_a_marker_is_refused(tmp_path):
    path = tmp_path / "expression.parquet"
    pd.DataFrame({"donor_id": ["A"], "object_id": ["A1"], "qc_analysis_eligible": [True],
                  "Z_offset": [1.0], "CD8": [10.0]}).to_parquet(path, index=False)
    import pytest
    with pytest.raises(stage.PanelNormalizationError, match="Z_offset"):
        stage._markers({"A": path}, stage.REGISTRY)
