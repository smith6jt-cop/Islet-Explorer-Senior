"""The spleen primary Cell Typing route (typing_route: spleen_protein): sections, windows and its own checks."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest

from multiplex5000.stages import spleen_typing as stage
from multiplex5000.stages import typing
from multiplex5000.tissues.lymph_node import Rule, RuleSource


def _project(tmp_path, sources=None):
    spec = tmp_path / "cohort.json"
    spec.write_text(json.dumps({"spec_version": 1, "images": [
        {"donor_id": d, "image_id": d, "panel_id": p}
        for p, ds in {"spleen_pc33": ["c33"], "spleen_pc67": ["a67", "b67"]}.items() for d in ds]}))
    named = {"cohort_spec": spec, **(sources or {})}
    return SimpleNamespace(cohort="spleen_fixture", release_dir=tmp_path / "release", settings={"typing_route": "spleen_protein"},
                           source=lambda name: named.get(name))


def test_the_route_is_chosen_by_name_and_only_pc67_is_typed(tmp_path):
    project = _project(tmp_path)
    assert typing._route(project) is stage
    s = stage.settings(project)
    assert s["typed"] == ["a67", "b67"]
    assert s["not_typed"] == {"spleen_pc33": ["c33"]} and "lineage roots" in s["not_typed_reason"]
    assert s["registry"] == ["phenocycler-pc67", "spleen"]


def test_inputs_need_the_typing_context_and_outlines(tmp_path):
    with pytest.raises(stage.SpleenTypingError, match="qupath_manifest"):
        stage.inputs(_project(tmp_path))


def test_the_coverage_window_holds_every_cell_and_the_shift_moves_them_into_it():
    feats = [{"id": "a", "geometry": {"type": "Polygon", "coordinates": [[[100, 200], [110, 200], [110, 215], [100, 200]]]}},
             {"id": "b", "geometry": {"type": "MultiPolygon", "coordinates": [
                 [[[5, 7], [9, 7], [9, 12], [5, 7]]], [[[300, 50], [305, 50], [305, 58], [300, 50]]]]}}]
    x0, y0, x1, y1 = stage._window(feats, (1000, 320))
    assert (x0, y0) == (0, 0)                     # margin clipped at the image edge
    assert x1 == 320 and y1 == 215 + stage.WINDOW_MARGIN_PX + 1
    moved = stage._shifted(feats, 3, 4)
    assert moved[0]["geometry"]["coordinates"][0][0] == [97, 196]
    assert moved[1]["geometry"]["coordinates"][1][0][0] == [297, 46]
    assert [f["id"] for f in moved] == ["a", "b"]


def test_a_class_whose_positive_clauses_are_not_met_is_counted():
    rules = RuleSource("0" * 64, (
        Rule("all", "Immune", (("CD45", "pos"),)),
        Rule("Immune", "T_cell", (("CD3e", "pos"),)),
        Rule("T_cell", "CD8_T", (("CD8", "pos"), ("CD4", "allneg"))),
        Rule("Immune", "Myeloid", (("CD68", "anypos"), ("CD163", "anypos"))),
    ), (), {}, "all", "any", ("Immune",))
    names = {"CD45": "CD45", "CD3e": "CD3e", "CD8": "CD8", "CD4": "CD4", "CD68": "CD68", "CD163": "CD163"}
    calls = {m: np.array(v, np.int8) for m, v in {
        "CD45": [1, 1, 1, 1], "CD3e": [1, 1, 0, 0], "CD8": [1, 0, 0, 0], "CD4": [0, 0, 0, 0],
        "CD68": [0, 0, 0, 0], "CD163": [0, 0, 1, 0]}.items()}
    refined = np.array(["CD8_T", "CD8_T", "Myeloid", "Myeloid"], dtype=object)
    assert stage._below_level(rules, refined, calls, names).tolist() == [0, 1, 0, 1]
