"""Identity through the Protein Measurement source, the tissue hooks, and what a source admits."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from multiplex5000.positivity import identity as ID
from multiplex5000.positivity import registry as RG
from multiplex5000.positivity import sources
from multiplex5000.positivity import tissue as TS

PANCREAS = ("codex-uniform-32", "pancreas")
HIGH, LOW = 1000.0, 10.0
#: object_id -> markers above both levels; c3 has no nucleus (not analysis-eligible)
CELLS = {"c0": ["CD45", "CD3e", "CD8a"], "c1": ["CD45"], "c2": [], "c3": ["CD45", "CD3e"],
         "c4": ["CD45", "CD3e", "CD68"], "c5": ["panCK", "ECAD"]}


def _stage(tmp_path, reg):
    """A Protein Measurement stage with one section, `s`, holding CELLS."""
    columns = sorted({reg.column(m, meas) for m in reg.measured_markers()
                      for meas in reg.candidates(m) + ([reg.differentiation(m)] if reg.differentiation(m) else [])})
    redsea = {"object_id": list(CELLS), "donor_id": ["s"] * len(CELLS)}
    for c in columns:
        marker = c.split("__", 1)[1]
        redsea[c] = [HIGH if marker in on else LOW for on in CELLS.values()]
    expression = {"donor_id": ["s"] * len(CELLS), "object_id": list(CELLS),
                  "qc_analysis_eligible": [oid != "c3" for oid in CELLS],
                  "X_centroid": [i * 100.0 for i in range(len(CELLS))], "Y_centroid": [0.0] * len(CELLS)}
    for table, frame in (("redsea", redsea), ("expression", expression)):
        path = tmp_path / "stage" / "data" / table / "donor_id=s" / "data.parquet"
        path.parent.mkdir(parents=True)
        pd.DataFrame(frame).to_parquet(path, index=False)
    return tmp_path / "stage"


def _tables(reg):
    markers = reg.typing_markers()
    level = np.log1p(100.0)
    rows = [{"image": "s", "marker": m, "measurement": reg.measurement(m), "level_log": level,
             "level_raw": 100.0} for m in markers + reg.markers("state")]
    diff = [{"image": "s", "marker": m, "measurement": reg.differentiation(m), "level_log": level,
             "level_raw": 100.0} for m in markers if reg.differentiation(m)]
    canary = [{"image": "s", "marker": m, "canary_log": np.log1p(200.0)} for m in markers]
    selection = [{"marker": m, "decision": "use"} for m in markers]
    return {"levels": pd.DataFrame(rows), "levels_differentiation": pd.DataFrame(diff),
            "canary_levels": pd.DataFrame(canary), "selection": pd.DataFrame(selection)}


def _coverage(tmp_path, reg, ids):
    frame = {"object_id": ids, "n_interior_pixels": [10] * len(ids)}
    for m in reg.typing_markers():
        frame[f"cov__{m}"] = [0.5] * len(ids)
        frame[f"covint__{m}"] = [0.5] * len(ids)
    path = tmp_path / "coverage" / "s.parquet"
    path.parent.mkdir()
    pd.DataFrame(frame).to_parquet(path, index=False)
    return path.parent


def test_identity_runs_end_to_end_on_protein_measurement_cells(tmp_path):
    reg = RG.load(*PANCREAS)
    stage = _stage(tmp_path, reg)
    source = sources.ProteinMeasurementSource(stage, coverage_dir=_coverage(tmp_path, reg, list(CELLS)))
    out, audit = ID.classify_image(source, "s", reg.spec, _tables(reg),
                                   start_class=TS.hooks("pancreas").start_class)
    out = out.set_index("object_id")
    assert "c3" not in out.index                          # no nucleus: never typed
    assert audit["cells"] == 5
    assert out.loc["c0", "label"] == "CD8a Tcell" and out.loc["c0", "compartment_status"] == "single"
    assert out.loc["c1", "label"] == "Hematopoietic CD45+"
    assert out.loc["c2", "label"] == "Tissue" and out.loc["c2", "status"] == "no_marker_above_background"
    assert out.loc["c4", "label"] == "Conflict T/Myeloid"  # equal depth, equal brightness and coverage
    assert out.loc["c5", "label"] == "Exocrine"


def test_a_section_level_missing_stops_identity(tmp_path):
    reg = RG.load(*PANCREAS)
    stage = _stage(tmp_path, reg)
    tables = _tables(reg)
    tables["levels"].loc[tables["levels"]["marker"] == "CD3e", "level_log"] = np.nan
    source = sources.ProteinMeasurementSource(stage, coverage_dir=_coverage(tmp_path, reg, list(CELLS)))
    with pytest.raises(ID.IdentityError, match="CD3e"):
        ID.classify_image(source, "s", reg.spec, tables)


def test_coverage_must_cover_every_redsea_cell(tmp_path):
    reg = RG.load(*PANCREAS)
    stage = _stage(tmp_path, reg)
    source = sources.ProteinMeasurementSource(stage, coverage_dir=_coverage(tmp_path, reg, ["c0", "c1", "c2"]))
    with pytest.raises(ValueError, match="differ from REDSEA"):
        ID.classify_image(source, "s", reg.spec, _tables(reg))


def test_the_source_admits_only_analysis_eligible_cells_and_fills_context(tmp_path):
    reg = RG.load(*PANCREAS)
    source = sources.ProteinMeasurementSource(_stage(tmp_path, reg))
    d = source.levels_frame("s", ["Cell__CD45"], list(sources.OVERLAP_COLUMNS))
    assert list(d["object_id"]) == ["c0", "c1", "c2", "c4", "c5"]
    assert (d["context"] == "tissue").all() and d["endocrine_object_class"].isna().all()
    assert (d[list(sources.OVERLAP_COLUMNS)] == 0).all().all()
    assert list(d["x_um"]) == [0.0, 100.0, 200.0, 400.0, 500.0]
    with pytest.raises(sources.SourceError, match="REDSEA has no"):
        source.levels_frame("s", ["Cell__NOT_A_MARKER"], [])


def test_a_context_file_supplies_objects_and_the_start_object(tmp_path):
    reg = RG.load(*PANCREAS)
    ctx = tmp_path / "context"
    ctx.mkdir()
    pd.DataFrame({"object_id": ["c0", "c1"], "context": ["Blood vessel", "tissue"],
                  "start_object": ["Follicle", "PALS"]}).to_parquet(ctx / "s.parquet", index=False)
    source = sources.ProteinMeasurementSource(_stage(tmp_path, reg), context_dir=ctx)
    d = source.identity_frame("s", ["Cell__CD45"], ["in_Blood vessel"])
    assert list(d["context"][:2]) == ["Blood vessel", "tissue"]
    assert list(TS.hooks("spleen").start_class(d)) == ["Follicle", "PALS", "Tissue", "Tissue", "Tissue"]
    assert list(TS.hooks("pancreas").start_class(d)) == ["Tissue"] * 5


def test_tissue_hooks_choose_their_selection_sections_and_references():
    pancreas, spleen = TS.hooks("pancreas"), TS.hooks("spleen")
    roster = {"ND_1": {"donor_status": "ND"}, "T1D_2": {"donor_status": "T1D"}, "ND_3": {"donor_status": "ND"}}
    assert pancreas.selection_sections(roster) == ["ND_1", "ND_3"]
    assert spleen.selection_sections(roster) == ["ND_1", "T1D_2", "ND_3"]      # all sections, blind to group
    assert pancreas.structure_objects() == ("Blood vessel", "Nerve", "Lymphatic")
    assert spleen.structure_objects() == ()
    reg = RG.load(*PANCREAS)
    TS.check_registry(reg, pancreas)
    with pytest.raises(RG.RegistryError, match="exocrine"):
        TS.check_registry(reg, spleen)                    # endocrine and exocrine references are pancreas-only


def test_a_tissue_without_hooks_is_refused():
    with pytest.raises(RG.RegistryError, match="typing hook"):
        TS.hooks("lymph_node")


def test_level_selection_runs_every_stage_and_stops_methods_it_cannot_show_working(tmp_path):
    """Two tiny sections: no background peak can be read (fewer than 1000 cells), so every option fails on them and the
    selection reports a stop for every marker instead of inventing a level; the canary step then refuses to run."""
    from multiplex5000.positivity import canary as CN
    from multiplex5000.positivity import levels as LV

    reg = RG.load(*PANCREAS)
    stage = _stage(tmp_path, reg)
    second = stage / "data"
    for table in ("redsea", "expression"):
        frame = pd.read_parquet(second / table / "donor_id=s" / "data.parquet")
        frame["donor_id"] = "t"
        (second / table / "donor_id=t").mkdir()
        frame.to_parquet(second / table / "donor_id=t" / "data.parquet", index=False)
    source = sources.ProteinMeasurementSource(stage)
    out = LV.choose(source, reg.spec, ["s", "t"], ["s", "t"], workers=1)
    assert set(out["stops"]) >= set(reg.typing_markers()) - {"PDPN"}
    assert (out["selection"]["decision"] == "stop").any()
    assert out["levels"]["level_log"].isna().all()
    assert set(out["levels"]["image"]) == {"s", "t"}
    with pytest.raises(LV.PositivityError, match="selection sections outside"):
        LV.choose(source, reg.spec, ["s"], ["t"], workers=1)
    with pytest.raises(LV.PositivityError, match="too few for the 99th percentile"):
        CN.levels(source, reg.spec, ["s"], {**out, "selection": out["selection"].assign(decision="use")})
