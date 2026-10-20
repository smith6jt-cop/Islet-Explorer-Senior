"""Ordinary endocrine admission and hormone subtype evidence are distinct decisions."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/typing"))
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/celltype"))
import subpanel as SUB  # noqa: E402

NAMES = {"INS": "Beta", "GCG": "Alpha", "SST": "Delta"}
UNRESOLVED = "Endocrine (unresolved subtype)"
AMBIGUOUS = "Endocrine (ambiguous subtype)"


def test_zero_negative_and_missing_comparisons_cannot_name_an_ordinary_subtype():
    values = pd.DataFrame({"INS": [0., -.1, .8, .8, np.inf, np.nan],
                           "GCG": [0., -.2, np.nan, .1, .1, np.nan],
                           "SST": [0., -.3, .1, .1, .1, np.nan]})
    before = values.copy(deep=True)
    names, status = SUB.ordinary_endocrine_identity(values, NAMES)
    assert names.tolist() == [UNRESOLVED, UNRESOLVED, UNRESOLVED, "Beta",
                              UNRESOLVED, UNRESOLVED]
    assert status.tolist() == ["no_usable_hormone", "no_usable_hormone",
                               "incomplete_hormone_comparison", "unique_highest_hormone",
                               "incomplete_hormone_comparison", "incomplete_hormone_comparison"]
    pd.testing.assert_frame_equal(values, before)


def test_positive_exact_ties_are_ambiguous_but_unique_winners_keep_their_names():
    values = pd.DataFrame({"INS": [1., .3, .8, .1, .1],
                           "GCG": [1., .3, .2, .8, .2],
                           "SST": [.1, .3, .1, .2, .8]})
    names, status = SUB.ordinary_endocrine_identity(values, NAMES)
    assert names.tolist() == [AMBIGUOUS, AMBIGUOUS, "Beta", "Alpha", "Delta"]
    assert status.tolist() == ["tied_highest_hormones", "tied_highest_hormones",
                               "unique_highest_hormone", "unique_highest_hormone",
                               "unique_highest_hormone"]


def test_a_wholly_unmeasured_hormone_channel_cannot_supply_complete_comparison():
    values = pd.DataFrame({"INS": [.9, .8], "GCG": [.01, .02], "SST": [np.nan, np.nan]})
    assert SUB.hormone_evidence(values)["complete"].tolist() == [False, False]
    assert SUB.hormone_evidence(values)["confirmed"].tolist() == [False, False]
    names, status = SUB.ordinary_endocrine_identity(values.drop(columns="SST"), NAMES)
    assert names.tolist() == [UNRESOLVED, UNRESOLVED]
    assert status.tolist() == ["incomplete_hormone_comparison"] * 2


def test_subtype_uncertainty_keeps_broad_endocrine_admission_and_each_route():
    values = pd.DataFrame({"INS": [0., 1., .9], "GCG": [0., 1., np.nan],
                           "SST": [0., .1, .1], "E_cadherin": [.95] * 3})
    names, _ = SUB.ordinary_endocrine_identity(values, NAMES)
    classes, definitions, resolved, reason = SUB._run_chain(
        values, {"bottom_tier_immune": []}, {}, np.array([None] * 3, object),
        np.array([""] * 3, object), np.ones(3, bool), np.ones(3, bool), names,
        np.array(["cd99", "elongation", "hormone"], object), np.ones(3, bool))
    assert classes.tolist() == [UNRESOLVED, AMBIGUOUS, UNRESOLVED]
    assert resolved.tolist() == ["cd99", "elongation", "hormone"]
    assert definitions.tolist() == [""] * 3
    assert reason.tolist() == [""] * 3


def test_the_separately_ruled_direct_single_policy_is_preserved():
    values = pd.DataFrame({"INS": [0., 1., .8], "GCG": [0., 1., np.nan],
                           "SST": [0., .1, .1]})
    names, status = SUB.direct_single_identity(values, NAMES)
    assert names.tolist() == ["Excluded", "Excluded", "Beta"]
    assert status.tolist() == ["no_usable_hormone", "tied_highest_hormones",
                               "unique_highest_hormone"]


def test_classify_preserves_keys_routes_and_broad_counts_while_reporting_subtype_evidence(monkeypatch):
    n = 13
    values = pd.DataFrame({
        "INS": [0., 1., .8, .8, .1, 0., .9] + [.8] * 6,
        "GCG": [0., 1., np.nan, .1, .8, 0., np.nan] + [.1] * 6,
        "SST": [0., .1, .1, .1, .1, 0., .1] + [.1] * 6,
        "CD99": [.8, .8, .8, .1, .1, .8, .8] + [.8] * 6,
        "CD3e": [.95, 0., 0., 0., .95, 0., 0.] + [0.] * 6,
        "CD38": [np.nan] * n, "E_cadherin": [0.] * n})
    ids = [f"literal-id-{i}" for i in reversed(range(n))]
    ctx = pd.DataFrame({"object_id": ids, "in_islet": [True] * n,
                        "islet_num": ["7"] * n, "context": ["islet"] * n,
                        "object_context": ["tissue"] * n, "cell_region": ["core"] * n,
                        "endocrine_object_class": [""] * 5 + ["EndocrineSingle"] * 2 + [""] * 6})
    levels = {"hormone": .5, "CD99": .5}
    log = {"subpanels": {}, "bottom_tier_immune": ["CD3e"],
           "hormones": {name: marker for marker, name in NAMES.items()}}
    monkeypatch.setattr(SUB, "immune_rules", lambda: {"Imm-fixture": (["CD3e"], [], ["CD38"])})
    monkeypatch.setattr(SUB, "scaled_markers", lambda *args: values.copy())
    monkeypatch.setattr(SUB, "identity_argmax", lambda values: {
        SUB.CORRECTION_C_COLUMN: np.array(["INS"] * n, object)})
    monkeypatch.setattr(SUB, "identity_argmax_audit", lambda *args: {})
    monkeypatch.setattr(SUB, "crossed_area", lambda *args: {})
    monkeypatch.setattr(SUB, "_structure_calls", lambda *args: (
        np.array([None] * n, object), np.array([""] * n, object), {}, {}))
    monkeypatch.setattr(SUB, "nucleus_shape", lambda *args: pd.DataFrame({
        SUB.SHAPE_FEATURE: [1., 1., 1., 2., 1., 1., 1.] + [1.] * 6}))
    monkeypatch.setattr(SUB, "fit_islet_levels", lambda *args: {
        s: {"n_confirmed": 6, "CD99": .5,
            SUB.SHAPE_FEATURE: 1.5 if s == SUB.S.LARGE_STRATUM else None}
        for s in SUB.S.SIZE_STRATA})
    audit = {}
    frame = SUB.classify("fixture-donor", ctx, levels, pd.DataFrame(), log, True, audit=audit)
    assert frame["object_id"].tolist() == ids
    assert frame["cell_class"].tolist()[:7] == [
        UNRESOLVED, AMBIGUOUS, UNRESOLVED, "Beta", "Immune", "Excluded", "Beta"]
    assert frame["resolved_by"].tolist()[:7] == [
        "cd99", "cd99", "cd99", "elongation", "supporting",
        "excluded_single_hormone", "single_hormone"]
    assert frame["endocrine_subtype_status"].tolist()[:7] == [
        "no_usable_hormone", "tied_highest_hormones", "incomplete_hormone_comparison",
        "unique_highest_hormone", "", "no_usable_hormone", "unique_highest_hormone"]
    assert audit["n_endocrine_new"] == 11
    assert audit["n_immune_to_endocrine"] == 1
    assert audit["ordinary_endocrine_subtype_status"]["no_usable_hormone"] == 1
    assert audit["ordinary_endocrine_subtype_status"]["tied_highest_hormones"] == 1
    assert audit["ordinary_endocrine_subtype_status"]["incomplete_hormone_comparison"] == 1
    assert frame["immune_rules_unavailable_count"].tolist()[:5] == [1, 0, 0, 0, 1]


def test_builder_retains_uncertain_endocrine_admission_and_evidence_columns(tmp_path, monkeypatch):
    # Reuse the established synthetic geometry fixture to exercise the actual cell builder.
    from test_build_cells_carry import producer_frame, build
    import build_cells as BC

    calls = producer_frame([UNRESOLVED, AMBIGUOUS, "Beta"],
                           ["cd99", "elongation", "hormone"])
    calls["endocrine_subtype_status"] = ["no_usable_hormone", "tied_highest_hormones",
                                         "unique_highest_hormone"]
    calls["immune_rules_unavailable_count"] = np.array([2, 1, 0], dtype=np.int16)
    frame = build(tmp_path, monkeypatch, calls)
    assert frame["is_endocrine"].tolist() == [True, True, True]
    assert frame["cell_class"].tolist() == [UNRESOLVED, AMBIGUOUS, "Beta"]
    assert frame["endocrine_route"].tolist() == ["cd99", "elongation", "hormone"]
    assert frame["endocrine_subtype_status"].tolist() == calls["endocrine_subtype_status"].tolist()
    assert frame["immune_rules_unavailable_count"].tolist() == [2, 1, 0]
    assert {UNRESOLVED, AMBIGUOUS} <= set(BC.ENDOCRINE)
