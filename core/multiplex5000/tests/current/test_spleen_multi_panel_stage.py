"""Panel selection must not share marker availability or fitted levels across PC33/PC67."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from multiplex5000.positivity import lineage_levels, registry, sources
from multiplex5000.stages import spleen_typing as stage
from multiplex5000.tissues import spleen


class ReviewTask:
    def __init__(self, root, panels):
        self.staging = root
        self.settings = {"panels": panels, "typed": [t for tags in panels.values() for t in tags], "workers": 2}

    def input(self, name):
        if name == "protein_receipt":
            return self.staging / "protein/receipt.json"
        if name.startswith("context:"):
            return self.staging / "context" / (name.split(":")[1] + ".parquet")
        return self.staging / name

    def output(self, name):
        path = self.staging / name
        path.parent.mkdir(parents=True, exist_ok=True)
        return path


def _run_fixture(tmp_path, monkeypatch, required_stop=False):
    panels = {"spleen_pc67": ["a67", "b67"], "spleen_pc33": ["c33"]}
    task = ReviewTask(tmp_path, panels)
    chosen, named, covered, api_panels = [], [], [], []
    image_list = [SimpleNamespace(donor_id=tag, panel_id=panel) for panel, tags in panels.items() for tag in tags]
    monkeypatch.setattr(stage, "_check_prerequisites", lambda task: None)
    monkeypatch.setattr(stage, "_review_code_identity", lambda: {"commit": "fixture", "files": {}})
    monkeypatch.setattr(stage, "_write_section_review_ready", lambda *args: None)
    monkeypatch.setattr(stage, "_manifest_from", lambda task: SimpleNamespace(images=image_list))
    monkeypatch.setattr(stage, "_level_figures", lambda *args: None)
    monkeypatch.setattr(stage, "_class_figure", lambda *args: None)

    def load(panel, tissue):
        return SimpleNamespace(spec=(panel, tissue), table=pd.DataFrame(
            {"required": ["yes", "no"]}, index=["CD3e", "FOXP3"]))

    monkeypatch.setattr(registry, "load", load)
    monkeypatch.setattr(registry, "require_reviewed", lambda reg: None)

    def choose(source, spec, tags, workers):
        chosen.append((spec, tags.copy()))
        is_pc33 = spec[0] == "phenocycler-pc33"
        selection = pd.DataFrame({"marker": ["CD3e", "FOXP3"],
                                  "decision": ["stop" if is_pc33 and required_stop else "use",
                                               "stop" if is_pc33 else "use"]})
        levels = pd.DataFrame({"image": tags, "marker": ["CD3e"] * len(tags),
                               "level_log": [3.3 if is_pc33 else 6.7] * len(tags)})
        return {"selection": selection, "options": levels.assign(b_cells=10, t_cells=12),
                "differentiation_selection": selection.copy(), "levels": levels,
                "levels_differentiation": levels.copy(),
                "stops": selection.loc[selection.decision == "stop", "marker"].tolist()}

    monkeypatch.setattr(lineage_levels, "choose", choose)
    monkeypatch.setattr(lineage_levels, "canary_levels", lambda tables, spec: tables["levels"].copy())

    class Source:
        def __init__(self, *, stage_dir, context_dir, coverage_dir=None):
            self.context_dir = context_dir

        def redsea_ids(self, tag):
            return [tag]

    monkeypatch.setattr(sources, "ProteinMeasurementSource", Source)

    def rules(path, *, panel_id):
        api_panels.append(panel_id)
        return SimpleNamespace(missing_by_phenotype=lambda observed: {"missing": {"never_acquired"}})

    monkeypatch.setattr(spleen, "load_rule_source", rules)
    monkeypatch.setattr(spleen, "observed_markers", lambda path, *, panel_id: {"CD3e", panel_id})
    monkeypatch.setattr(spleen, "structure_only_branches", lambda path, *, panel_id: {"Endothelial"})
    monkeypatch.setattr(spleen, "panel_adaptation_audit", lambda path, *, panel_id: {"panel_id": panel_id})
    monkeypatch.setattr(stage, "_coverage", lambda task, image, tag, levels, reg, ids: covered.append((tag, reg.spec)))

    def name(task, source, tag, image, reg, tables, canary, rules, branches, *, panel_id):
        named.append((tag, panel_id, reg.spec, tables["selection"].copy()))
        return pd.DataFrame({"cell_class": ["T_cell"], "refined_class": ["Treg"], "status": ["unique"],
                             "start_class": ["Tissue"], "violations": [0], "unsupported_t_cell": [False]})

    monkeypatch.setattr(stage, "_named", name)
    return stage.run(task), task, chosen, named, covered, api_panels


def test_all_donors_fit_within_their_panel_and_optional_stop_does_not_cross_panels(tmp_path, monkeypatch):
    summary, task, chosen, named, covered, api_panels = _run_fixture(tmp_path, monkeypatch)
    assert chosen == [(("phenocycler-pc33", "spleen"), ["c33"]),
                      (("phenocycler-pc67", "spleen"), ["a67", "b67"])]
    assert api_panels == ["spleen_pc33", "spleen_pc67"]
    assert summary["typed"] and set(summary["per_section"]) == {"a67", "b67", "c33"}
    assert len(covered) == len(named) == 3
    assert summary["panels"]["spleen_pc67"]["markers"]["stop"] == []
    assert summary["panels"]["spleen_pc33"]["markers"]["stop"] == ["FOXP3"]
    assert [row[3].set_index("marker").loc["FOXP3", "decision"] for row in named] == ["stop", "use", "use"]
    aggregate = pd.read_csv(tmp_path / "tables/selection.csv")
    pc67 = pd.read_csv(tmp_path / "tables/by_panel/spleen_pc67/selection.csv")
    extracted = stage._panel_tables({"selection": aggregate}, "spleen_pc67")["selection"]
    pd.testing.assert_frame_equal(extracted.reset_index(drop=True), pc67)
    levels = pd.read_csv(tmp_path / "tables/levels.csv")
    assert levels.set_index("image")["level_log"].to_dict() == {"a67": 6.7, "b67": 6.7, "c33": 3.3}
    checks = stage.gate(task, summary).checks
    blocking = [check for check in checks if check.blocking]
    assert len(blocking) == 5 and all(check.passed for check in blocking)


def test_required_level_failure_in_either_panel_stops_naming_and_is_reported_by_panel(tmp_path, monkeypatch):
    summary, task, chosen, named, covered, _ = _run_fixture(tmp_path, monkeypatch, required_stop=True)
    assert len(chosen) == 2 and not named and not covered and not summary["typed"]
    check = stage.gate(task, summary).checks[0]
    assert check.blocking and not check.passed
    assert check.value == {"spleen_pc67": [], "spleen_pc33": ["CD3e"]}


def test_single_panel_table_keeps_the_original_schema(tmp_path):
    task = ReviewTask(tmp_path, {"spleen_pc67": ["a"]})
    table = pd.DataFrame({"image": ["a"], "marker": ["CD3e"], "level_log": [1.2345]})
    stage._write_panel_table(task, "levels", {"spleen_pc67": table})
    assert (tmp_path / "tables/levels.csv").read_bytes() == (tmp_path / "tables/by_panel/spleen_pc67/levels.csv").read_bytes()


@pytest.mark.parametrize("shape,axes,indices", [((2, 8, 8, 8), "TCYX", [0, 1]),
                                               ((8, 8, 2), "YXC", [0, 1]),
                                               ((2, 8, 8), "CYX", [0, 2]),
                                               ((2, 8, 8), "CYX", [-1, 1]),
                                               ((2, 8, 8), "CYX", [1, 1])])
def test_coverage_rejects_ambiguous_shapes_or_channel_maps(shape, axes, indices):
    image = SimpleNamespace(donor_id="pc33", channel_map=[
        SimpleNamespace(channel_name=name, channel_index=i) for name, i in zip(["DAPI", "CD3e"], indices)])
    with pytest.raises(stage.SpleenTypingError):
        stage._coverage_channel_index(image, shape, axes)


def test_coverage_reads_the_declared_channel_from_an_ome_tiff(tmp_path):
    import tifffile

    tag = "pc33"
    task = ReviewTask(tmp_path, {"spleen_pc33": [tag]})
    path = tmp_path / "small.ome.tiff"
    pixels = np.zeros((2, 24, 24), dtype=np.uint16)
    pixels[1] = 100
    tifffile.imwrite(path, pixels, metadata={"axes": "CYX"}, photometric="minisblack")
    geometry = {"type": "Polygon", "coordinates": [[[3, 3], [20, 3], [20, 20], [3, 20], [3, 3]]]}
    geo = tmp_path / "cells.geojson"
    geo.write_text(json.dumps({"type": "FeatureCollection", "features": [{"type": "Feature", "id": "cell1", "geometry": geometry}]}))
    task.input = lambda name: geo if name.startswith("cell_geojson:") else path
    image = SimpleNamespace(donor_id=tag, qptiff=SimpleNamespace(path=str(path)), channel_map=[
        SimpleNamespace(channel_name="DAPI", channel_index=0), SimpleNamespace(channel_name="CD3e", channel_index=1)])
    reg = SimpleNamespace(typing_markers=lambda: ["CD3e"], table=pd.DataFrame({"channel": ["CD3e"]}, index=["CD3e"]))
    levels = pd.DataFrame({"image": [tag], "marker": ["CD3e"], "level_raw": [50.0], "measurement": ["whole_cell"]})
    frame = stage._coverage(task, image, tag, levels, reg, ["cell1"])
    assert frame.object_id.tolist() == ["cell1"]
    assert frame["covint__CD3e"].tolist() == [1.0]
    assert (tmp_path / "data/coverage/pc33.parquet").is_file()


@pytest.mark.parametrize("include_excluded,dirty,tamper", [(False, False, False), (True, False, False),
                                                           (False, True, False), (False, False, True)])
def test_provisional_ready_marker_requires_exact_eligible_uuids_and_binds_finished_files(tmp_path, include_excluded, dirty, tamper):
    from multiplex5000.tasks import Input, Task

    tag, panel = "c33", "spleen_pc33"
    inputs = {}
    for name in set(stage.PANEL_FILES) | {"protein_receipt", "protein_acceptance", "normalization_receipt",
                                         "normalization_acceptance", "cohort_spec", "qupath_manifest",
                                         f"context:{tag}", f"registry:{panel}"}:
        path = tmp_path / "inputs" / name.replace(":", "_")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
        inputs[name] = Input(name, path)
    expression = tmp_path / "expression.parquet"
    pd.DataFrame({"object_id": ["admitted", "excluded"], "qc_analysis_eligible": [True, False]}).to_parquet(expression)
    inputs[f"expression:{tag}"] = Input(f"expression:{tag}", expression)
    redsea = tmp_path / "redsea.parquet"
    redsea.write_bytes(expression.read_bytes())
    inputs[f"redsea:{tag}"] = Input(f"redsea:{tag}", redsea)
    inputs["protein_receipt"].path.write_text(json.dumps({"outputs": [stage._file_record(expression), stage._file_record(redsea)]}))
    if tamper:
        redsea.write_bytes(b"changed after acceptance")
    task = Task(stage.STAGE, "fixture", tmp_path / ".06_typing.incomplete-fixture", inputs)
    cells = pd.DataFrame({"object_id": ["excluded" if include_excluded else "admitted"],
                          "cell_class": ["T_cell"], "refined_class": ["Treg"], "status": ["unique"],
                          "start_class": ["Tissue"], "violations": [0], "unsupported_t_cell": [False]})
    cells.to_parquet(task.output(f"data/cells/donor_id={tag}/data.parquet"), index=False)
    task.output(f"data/coverage/{tag}.parquet").write_bytes(b"coverage")
    for name in (*stage.TABLES, "canary_levels"):
        task.output(f"tables/by_panel/{panel}/{name}.csv").write_text("fixture")
    task.output(f"tables/by_panel/{panel}/rule_adaptations.json").write_text("{}")
    args = (task, tag, panel, cells, {"markers": {"required_stopped": []}},
            {"commit": "f" * 40, "dirty": dirty, "files": {}})
    if tamper:
        with pytest.raises(stage.SpleenTypingError, match="redsea input differs"):
            stage._write_section_review_ready(*args)
        return
    result = stage._write_section_review_ready(*args)
    marker = task.staging / f"review/ready/{tag}.json"
    if include_excluded or dirty:
        assert result is None and not marker.exists()
    else:
        assert result == marker
        record = json.loads(marker.read_text())
        assert record["status"] == "provisional" and record["code_commit"] == "f" * 40
        assert len(record["checks"]) == 5 and all(c["passed"] for c in record["checks"])
        assert record["eligibility"]["expected_cells"] == 1 and record["eligibility"]["uuid_set_matches"]
        assert record["eligibility"]["expression_output"] == stage._file_record(expression)
        assert record["eligibility"]["redsea_output"] == stage._file_record(redsea)
        assert record["eligibility"]["source_hashes_verified"]
        assert record["outputs"]["cells"] == stage._file_record(task.staging / f"data/cells/donor_id={tag}/data.parquet")
        assert not marker.with_suffix(".json.tmp").exists()
