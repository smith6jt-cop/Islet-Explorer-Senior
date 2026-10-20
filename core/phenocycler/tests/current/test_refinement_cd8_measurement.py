"""Verify the refinement-only CD8 measurement ruling on keyed input artifacts."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
SCRIPT = REPO / "core/phenocycler/scripts/typing/refine.py"
sys.path.insert(0, str(SCRIPT.parent))
sys.path.insert(0, str(SCRIPT.parent.parent / "celltype"))
spec = importlib.util.spec_from_file_location("refine_cd8_test", SCRIPT)
refine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(refine)


def test_only_cd8_is_restored_without_mutating_primary_values():
    primary = pd.DataFrame({"object_id": ["c", "a", "b"],
                            "CD8": [24.0, np.nan, 0.0], "CD3e": [3.0, 9.0, 2.0]})
    original = primary.copy(deep=True)
    got = refine.restore_unshifted_cd8(primary, "0112",
                                     {"shifts": {"CD8": {"scale": {"0112": 6.0}}}})
    pd.testing.assert_frame_equal(primary, original)
    assert got.object_id.tolist() == ["c", "a", "b"]
    assert got.CD3e.tolist() == original.CD3e.tolist()
    assert got.CD8.iloc[0] == 4.0
    assert np.isnan(got.CD8.iloc[1])
    assert got.CD8.iloc[2] == 0.0


@pytest.mark.parametrize("multiplier", [0.0, -2.0, float("nan"), float("inf"), "invalid"])
def test_invalid_multiplier_is_refused(multiplier):
    with pytest.raises(ValueError, match="CD8 alignment multiplier"):
        refine.restore_unshifted_cd8(pd.DataFrame({"CD8": [10.0]}), "112",
                                    {"shifts": {"CD8": {"scale": {"112": multiplier}}}})


def test_donor_alias_is_not_guessed():
    with pytest.raises(ValueError, match="unavailable"):
        refine.restore_unshifted_cd8(pd.DataFrame({"CD8": [10.0]}), "0112",
                                    {"shifts": {"CD8": {"scale": {"112": 2.0}}}})


def test_frozen_reference_uses_same_units_as_unshifted_measurement():
    scale = pd.DataFrame({"marker": ["CD8", "CD3e"], "floor": [70.0, 3.0],
                          "top": [770.0, 13.0], "span": [700.0, 10.0]}).set_index("marker")
    original = scale.copy(deep=True)
    shifts = {"shifts": {"CD8": {"scale": {"0112": 7.0}}}}
    primary = pd.DataFrame({"CD8": [0.0, 140.0, 840.0, np.nan],
                            "CD3e": [2.0, 8.0, 15.0, 4.0]})
    unshifted = refine.restore_unshifted_cd8(primary, "0112", shifts)
    converted = refine.scale_for_unshifted_cd8(scale, "0112", shifts)
    before = refine.S.apply_scale(primary.to_numpy(), ["CD8", "CD3e"], scale)
    after = refine.S.apply_scale(unshifted.to_numpy(), ["CD8", "CD3e"], converted)
    np.testing.assert_allclose(after, before, equal_nan=True)
    pd.testing.assert_frame_equal(scale, original)
    pd.testing.assert_series_equal(converted.loc["CD3e"], original.loc["CD3e"])
    assert converted.loc["CD8", "floor"] == 10.0
    assert converted.loc["CD8", "span"] == 100.0


def test_donor_reader_uses_matching_cd8_record_and_keyed_geometry(tmp_path, monkeypatch):
    aligned = tmp_path / "aligned/donor_id=0112/data.parquet"
    aligned.parent.mkdir(parents=True)
    pd.DataFrame({"object_id": ["b", "a"], "CD8": [28.0, 7.0],
                  "CD3e": [8.0, 12.0]}).to_parquet(aligned, index=False)
    geometry = tmp_path / "geometry.parquet"
    pd.DataFrame({"object_id": ["a", "b"],
                  "qc_analysis_eligible": [False, True]}).to_parquet(geometry, index=False)
    shifts = tmp_path / "shifts.json"
    shifts.write_text(json.dumps({"shifts": {"CD8": {"scale": {"0112": 7.0}}}}))
    monkeypatch.setattr(refine.S, "ALIGNED", tmp_path / "aligned")
    monkeypatch.setattr(refine.S, "cell_source", lambda donor: (geometry, "revised"))
    monkeypatch.setattr(refine.S, "mask_fragment_hormones", lambda values, donor: values)
    monkeypatch.setattr(refine, "SHIFTS", shifts)
    monkeypatch.setattr(refine, "CD8_ALIGNMENT", "unshifted")
    got, have = refine.donor_values("0112", ["CD8", "CD3e"])
    assert have == ["CD8", "CD3e"]
    assert got.object_id.tolist() == ["b", "a"]
    assert got.CD8.tolist() == [4.0, 1.0]
    assert got.CD3e.tolist() == [8.0, 12.0]
    assert got.qc_analysis_eligible.tolist() == [True, False]
    monkeypatch.setattr(refine, "CD8_ALIGNMENT", "aligned")
    unchanged, _ = refine.donor_values("0112", ["CD8", "CD3e"])
    assert unchanged.CD8.tolist() == [28.0, 7.0]


def test_cached_representation_cannot_cross_cd8_measurement_modes(tmp_path):
    with pytest.raises(ValueError, match="fresh representation"):
        refine.main(["--out", str(tmp_path / "out"), "--cd8-alignment", "unshifted",
                     "--params", str(tmp_path / "old.parquet")])


def test_refinement_audit_distinguishes_broad_endocrine_from_validated_subtype(tmp_path, monkeypatch):
    uncertain = ["Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)"]
    labels = np.array(["Beta", *uncertain, "Exocrine"])
    strata = np.array(["islet"] * 4)
    frame = pd.DataFrame({"size_stratum": strata})
    monkeypatch.setattr(refine, "CD8_ALIGNMENT", "aligned")
    monkeypatch.setattr(refine, "SHIFTS", tmp_path / "absent-shifts.json")
    log = refine._log("fixture", frame, np.ones(4, bool), [], {}, labels, labels,
                      np.zeros(4, bool), np.zeros(4, bool), np.ones(4), np.ones(4),
                      strata, 0.0)
    evidence = log["endocrine_subtype_uncertainty"]
    assert evidence["classes"] == uncertain
    assert evidence["n_before"] == evidence["n_after"] == 2
    assert evidence["subtype_defining_marker_test_applied"] is False
    assert not (set(uncertain) & set(log["own_markers"]))
