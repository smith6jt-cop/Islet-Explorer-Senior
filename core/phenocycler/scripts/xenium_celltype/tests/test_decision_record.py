"""One shape for three lanes is only worth having if it cannot launder an absence.

Each adapter is checked against its lane's REAL output rather than a hand-typed row -- the
whole failure mode being adapted away is that each lane's private convention for "nothing was
called" (``lineage == ""``, ``call == "Unscored"``, ``broad_n_identity_genes == 0``) is
invisible to a reader who did not write it, so a test that retypes those conventions by hand
tests my reading of them rather than the lanes.

The load-bearing test is the last group: an unavailable record must survive a DataFrame and a
CSV without any of its NaNs becoming 0, because that conversion is the specific defect that
reported donor 6539's spatial arm as ``Immune_T`` 0.0% when Immune_T had cleared its null and
been outranked in six groups.
"""
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1] / "celltype"))

from decision_record import (  # noqa: E402
    AMBIGUOUS, AWARDED, CLEARED_BUT_OUTRANKED, NO_GROUP_CLEARED, UNAVAILABLE, UNSCORED,
    DecisionRecord, DecisionRecordError, assert_no_absence_as_zero,
    frame_to_records, from_broad_call, from_broad_obs, from_cell_call, from_cell_calls,
    from_composition, from_composition_row, from_group_call, from_group_calls,
    records_to_frame,
)
from nullgate import annotate_groups, MIN_GROUP_CELLS  # noqa: E402

GENES = [f"g{i}" for i in range(200)]


def _counts(n_cells=400, seed=0, hot=None, hot_rows=None):
    """Sparse counts like Xenium: ~96% zeros, a few genes hot in a few cells."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.04, size=(n_cells, len(GENES))).astype(float)
    if hot is not None:
        X[np.ix_(hot_rows, hot)] += rng.poisson(6.0, size=(len(hot_rows), len(hot)))
    return X


# ---------------------------------------------------------------------------
# adapter 1 -- nullgate
# ---------------------------------------------------------------------------
def test_an_awarded_group_keeps_its_null_and_is_given_no_confidence():
    """The gate's z ranks; it does not measure a probability, and the lane says so.

    Copying it into ``confidence`` would sit it beside a scANVI posterior as though a reader
    could compare the two. It is carried in ``null.z``, where its name says what it is.
    """
    X = _counts(hot=list(range(10)), hot_rows=np.arange(200))
    groups = np.where(np.arange(400) < 200, "A", "B")
    ann = annotate_groups(X, GENES, groups,
                          {"Real": GENES[:10], "Other": GENES[100:110]},
                          n_draws=40, seed=1).set_index("group")
    rec = from_group_call(ann.loc["A"], unit=("k20", "A"))

    assert rec.status == AWARDED
    assert rec.label == "Real"
    assert rec.n == 200
    assert rec.statistic.name == "contrast"          # not a cosine, not a posterior
    assert rec.statistic.value > rec.null.p95        # the award is "beat its own null"
    assert rec.null.kind == "size_and_detection_matched"
    assert np.isfinite(rec.null.z)
    assert np.isnan(rec.confidence), "a ranking device must not be rendered as a confidence"


def test_a_refused_group_carries_its_reason_and_no_label():
    X = _counts(n_cells=600, seed=3)                 # pure background, no lineage present
    groups = np.where(np.arange(600) < 300, "A", "B")
    ann = annotate_groups(X, GENES, groups,
                          {"A": GENES[:10], "B": GENES[20:30], "C": GENES[40:48]},
                          n_draws=40, seed=2)
    recs = list(from_group_calls(ann))

    assert all(r.status == UNAVAILABLE for r in recs)
    assert all(r.label is None for r in recs)        # never "" and never the nearest core
    assert all(r.reason == "no_core_beats_its_matched_null" for r in recs)
    assert all(np.isnan(r.confidence) for r in recs)
    # the losing contrast IS measured and stays measured -- refusing a call is not refusing
    # to report the number the refusal rests on
    assert all(np.isfinite(r.statistic.value) for r in recs)


def test_a_group_too_small_to_pool_has_a_null_that_was_never_computed():
    """``null=None`` (the lane computes no null) and a null of NaNs (one belonged and was not
    computed) are different states, and a reader needs to tell them apart."""
    X = _counts(n_cells=300)
    groups = np.where(np.arange(300) < MIN_GROUP_CELLS - 1, "tiny", "big")
    ann = annotate_groups(X, GENES, groups, {"A": GENES[:10]}, n_draws=5).set_index("group")
    rec = from_group_call(ann.loc["tiny"])

    assert rec.status == UNAVAILABLE and rec.reason == "too_few_cells"
    assert rec.null is not None                      # a null belonged here
    assert np.isnan(rec.null.mean) and np.isnan(rec.statistic.value)   # and was not computed


def test_cores_that_cleared_but_were_not_separable_keep_every_contender():
    X = _counts(n_cells=600, seed=4, hot=range(10), hot_rows=range(300))
    groups = np.array(["A"] * 300 + ["B"] * 300)
    ann = annotate_groups(X, GENES, groups,
                          {"CoreOne": GENES[:10], "CoreTwo": list(GENES[:10])},
                          n_draws=40, seed=0).set_index("group")
    rec = from_group_call(ann.loc["A"])

    assert rec.status == UNAVAILABLE
    assert rec.reason == "cores_cleared_but_are_not_separable"
    assert set(rec.runners_up) == {"CoreOne", "CoreTwo"}   # the top one included: none won


# ---------------------------------------------------------------------------
# adapter 2 -- celltype/score.call_from_scores
# ---------------------------------------------------------------------------
def _score_lane_calls():
    """Real ``call_from_scores`` output holding one of each outcome.

    Cell 0 matches reference P0 and is allowed it; cell 1 met no requirement (``Ambiguous``);
    cell 2 had no finite score at all (``Unscored``).
    """
    from score import call_from_scores
    S = np.array([[0.90, 0.40], [0.50, 0.10], [np.nan, np.nan]], dtype=np.float32)
    allowed = np.array([[True, True], [False, False], [True, True]])
    return call_from_scores(S, ["P0", "P1"], allowed=allowed)


def test_a_scored_cell_puts_its_cosine_in_the_statistic_not_the_confidence():
    calls = _score_lane_calls()
    rec = from_cell_call(calls.iloc[0], unit="obj-1")

    assert rec.status == AWARDED and rec.label == "P0" and rec.n == 1
    assert rec.statistic.name == "cosine"
    assert rec.statistic.value == pytest.approx(0.90, abs=1e-6)
    assert rec.statistic.margin == pytest.approx(0.50, abs=1e-6)
    assert rec.runners_up == ("P1",)
    assert rec.null is None, "this lane compares against no null; NaN would imply one belonged"
    assert np.isnan(rec.confidence), "a ceiling-normalised cosine is a similarity, not a probability"


def test_a_cell_meeting_no_requirement_is_unavailable_with_a_measured_zero_beside_it():
    """``support == 0`` is a measurement -- zero references met their rules -- and it stands
    in the same record as a NaN score. Those two columns disagreeing is the point."""
    calls = _score_lane_calls()
    rec = from_cell_call(calls.iloc[1], unit="obj-2")

    assert calls.iloc[1]["call"] == AMBIGUOUS
    assert rec.status == UNAVAILABLE and rec.label is None
    assert rec.reason == "no_reference_requirements_met"
    assert rec.statistic.support == 0                # measured zero, written as 0
    assert np.isnan(rec.statistic.value)             # unmeasured, written as NaN
    assert np.isnan(rec.confidence)


def test_a_cell_with_nothing_to_compare_says_so_rather_than_scoring_zero():
    calls = _score_lane_calls()
    recs = list(from_cell_calls(calls, units=["a", "b", "c"]))

    assert calls.iloc[2]["call"] == UNSCORED
    assert recs[2].status == UNAVAILABLE and recs[2].label is None
    assert recs[2].reason == "nothing_to_compare"
    assert np.isnan(recs[2].statistic.value)
    assert [r.unit for r in recs] == ["a", "b", "c"]


def test_a_positional_identifier_is_refused():
    calls = _score_lane_calls()
    with pytest.raises(DecisionRecordError, match="must not be positional"):
        list(from_cell_calls(calls, units=["a", "b"]))


# ---------------------------------------------------------------------------
# adapter 3 -- hierarchical_phenotype tier 1
# ---------------------------------------------------------------------------
def _broad_obs():
    """Tier-1 obs as ``hierarchical_phenotype`` writes it.

    Cell 2 is the case that matters: ``broad_lowsignal`` marks a cell expressing no
    identity-core gene of ANY lineage, which still receives the absolute-score argmax and a
    posterior. The lane has no unavailable state; the adapter is where that becomes visible.
    """
    return pd.DataFrame(
        {"celltype_broad": pd.Categorical(["Immune_T", "Exocrine", "Immune_B"]),
         "broad_confidence": np.array([0.94, 0.61, 0.88], dtype=np.float32),
         "broad_n_identity_genes": np.array([5, 2, 0], dtype=np.int16),
         "broad_lowsignal": np.array([False, False, True]),
         "broad_marker_margin": np.array([0.70, 0.12, 0.30], dtype=np.float32)},
        index=["c0", "c1", "c2"])


def test_a_broad_call_is_the_one_lane_with_a_real_probability():
    rec = from_broad_call(_broad_obs().iloc[0], unit="c0")

    assert rec.status == AWARDED and rec.label == "Immune_T"
    assert rec.confidence == pytest.approx(0.94, abs=1e-6)     # a posterior, genuinely [0, 1]
    assert rec.statistic.name == "scanvi_posterior"
    assert rec.statistic.support == 5                          # identity genes behind it
    assert rec.null is None, "and it is the lane with no null -- not comparable to a cleared one"


def test_a_cell_expressing_no_identity_gene_is_not_an_award():
    """The lane's label is Indeterminate-free by contract, so a cell with no identity
    evidence still gets one. Here it becomes ``unavailable`` and the lane's pick moves to
    ``runners_up``, where it reads as a contender rather than a call."""
    recs = list(from_broad_obs(_broad_obs()))
    low = recs[2]

    assert low.unit == "c2"
    assert low.status == UNAVAILABLE and low.label is None
    assert low.reason == "no_identity_gene_detected"
    assert low.runners_up == ("Immune_B",)          # what the lane picked, not laundered away
    assert low.statistic.support == 0               # measured: zero genes detected
    assert np.isnan(low.statistic.value)            # the 0.88 posterior is dropped, not zeroed
    assert np.isnan(low.confidence)


def test_obs_without_the_tier_one_column_is_refused():
    with pytest.raises(DecisionRecordError, match="not tier-1 output"):
        list(from_broad_obs(pd.DataFrame({"celltype": ["Acinar"]})))


# ---------------------------------------------------------------------------
# adapter 4 -- run_donor.composition, where the vocabulary comes from
# ---------------------------------------------------------------------------
def _composition():
    from run_donor import composition
    ann = pd.DataFrame({
        "group": [0, 1, 2, 3],
        "lineage": ["Exocrine", "", "Endocrine", ""],
        "reason": ["", "no_core_beats_its_matched_null", "", "too_few_cells"],
        "runners_up": ["Immune_B", "", "", ""],
    })
    groups = np.array([0] * 500 + [1] * 300 + [2] * 150 + [3] * 50)
    return composition(ann, groups,
                       lineages=["Exocrine", "Endocrine", "Immune_B", "Immune_T"])


def test_all_four_statuses_survive_the_adapter():
    recs = {(r.label, r.reason): r for r in from_composition(_composition(), unit="6539")}
    by_label = {r.label: r for r in recs.values() if r.label}

    assert by_label["Exocrine"].status == AWARDED
    assert by_label["Exocrine"].statistic.name == "frac"
    assert np.isfinite(by_label["Exocrine"].statistic.value)

    # outranked is not absent, and neither is a zero
    assert by_label["Immune_B"].status == CLEARED_BUT_OUTRANKED
    assert np.isnan(by_label["Immune_B"].statistic.value)
    assert by_label["Immune_T"].status == NO_GROUP_CLEARED
    assert np.isnan(by_label["Immune_T"].statistic.value)

    # the two refusals stay apart: one is about pooling, the other about evidence
    reasons = {r.reason for r in recs.values() if r.status == UNAVAILABLE}
    assert reasons == {"too_few_cells", "no_core_beats_its_matched_null"}
    # an unavailable row's share IS measured -- those cells were counted
    assert all(np.isfinite(r.statistic.value) for r in recs.values()
               if r.status == UNAVAILABLE)


def test_a_label_that_won_nothing_cannot_be_given_a_measured_share():
    """The Immune_T 0.0% regression, refused at construction."""
    row = {"lineage": "Immune_T", "status": NO_GROUP_CLEARED, "n_cells": 0,
           "frac": 0.0, "ci_lo": 0.0, "ci_hi": 0.0}
    with pytest.raises(DecisionRecordError, match="won nothing here"):
        from_composition_row(row, unit="6539")


# ---------------------------------------------------------------------------
# the schema's own rule
# ---------------------------------------------------------------------------
def test_a_non_awarded_record_cannot_carry_a_confidence_of_zero():
    with pytest.raises(DecisionRecordError, match="absence rendered as"):
        DecisionRecord(unit="g1", n=100, label=None, reason="too_few_cells",
                       confidence=0.0, status=UNAVAILABLE)


def test_an_unavailable_record_cannot_keep_a_label():
    with pytest.raises(DecisionRecordError, match="runners_up"):
        DecisionRecord(unit="g1", n=100, label="Exocrine", reason="too_few_cells",
                       status=UNAVAILABLE)


def test_a_refusal_must_keep_its_reason():
    with pytest.raises(DecisionRecordError, match="keeps its reason"):
        DecisionRecord(unit="g1", n=100, label=None, reason="", status=UNAVAILABLE)


def test_only_a_probability_may_be_called_a_confidence():
    with pytest.raises(DecisionRecordError, match=r"outside \[0, 1\]"):
        DecisionRecord(unit="c0", n=1, label="Immune_T", confidence=3.6, status=AWARDED)


def test_a_lane_reason_is_not_a_status():
    with pytest.raises(DecisionRecordError, match="not one of"):
        DecisionRecord(unit="g1", n=100, label=None, status="too_few_cells",
                       reason="too_few_cells")


# ---------------------------------------------------------------------------
# the point of the module: an absence stays an absence
# ---------------------------------------------------------------------------
def test_an_unavailable_record_round_trips_without_becoming_zero(tmp_path):
    """Through a DataFrame and through CSV, which is how these lanes actually persist.

    ``run_donor`` writes ``composition.csv`` and ``group_calls.csv``, so the round trip is
    not hypothetical: a NaN that survives the constructor and dies in ``read_csv`` is the same
    bug one step later.
    """
    src = from_broad_call(_broad_obs().iloc[2], unit="c2")    # lowsignal: nothing measured
    assert src.status == UNAVAILABLE

    frame = records_to_frame([src])
    assert frame.loc[0, "confidence"] != 0
    assert pd.isna(frame.loc[0, "confidence"])
    assert pd.isna(frame.loc[0, "statistic_value"])
    assert frame.loc[0, "statistic_support"] == 0             # the measured zero survives as 0

    path = tmp_path / "records.csv"
    frame.to_csv(path, index=False)
    back = frame_to_records(pd.read_csv(path))[0]

    assert back.status == UNAVAILABLE
    assert back.label is None                                 # not "", not "nan", not a call
    assert back.reason == "no_identity_gene_detected"
    assert np.isnan(back.confidence), "an absence came back as a number"
    assert np.isnan(back.statistic.value)
    assert back.statistic.support == 0
    assert back.runners_up == ("Immune_B",)
    assert back.null is None                                  # the lane still has no null
    assert_no_absence_as_zero(pd.read_csv(path))


def test_every_lane_round_trips_and_the_guard_catches_a_fillna():
    """One frame, four lanes -- and ``fillna(0)`` on it is caught rather than shipped."""
    recs = [from_group_call(annotate_groups(
                _counts(n_cells=300), GENES,
                np.where(np.arange(300) < MIN_GROUP_CELLS - 1, "tiny", "big"),
                {"A": GENES[:10]}, n_draws=5).set_index("group").loc["tiny"]),
            from_cell_call(_score_lane_calls().iloc[1], unit="obj-2"),
            from_broad_call(_broad_obs().iloc[2], unit="c2"),
            *from_composition(_composition(), unit="6539")]

    frame = records_to_frame(recs)
    assert_no_absence_as_zero(frame)
    # NaN-aware: a round trip that turned any absence into a number fails here
    pd.testing.assert_frame_equal(records_to_frame(frame_to_records(frame)), frame)

    with pytest.raises(DecisionRecordError, match="not zero, not negative"):
        assert_no_absence_as_zero(frame.fillna(0))


def test_one_shape_does_not_make_the_evidence_uniform():
    """The frame is uniform; the numbers in it are not, and the record keeps that legible.

    A reader comparing 0.9 to 0.9 across these rows has to pass ``statistic`` and ``null_kind``
    first: a contrast that beat a matched null, a ceiling-normalised cosine and a scANVI
    posterior are not the same quantity however similar the columns look.
    """
    X = _counts(hot=list(range(10)), hot_rows=np.arange(200))
    groups = np.where(np.arange(400) < 200, "A", "B")
    gated = annotate_groups(X, GENES, groups, {"Real": GENES[:10], "Other": GENES[100:110]},
                            n_draws=40, seed=1).set_index("group")
    frame = records_to_frame([
        from_group_call(gated.loc["A"]),
        from_cell_call(_score_lane_calls().iloc[0], unit="obj-1"),
        from_broad_call(_broad_obs().iloc[0], unit="c0"),
    ])

    assert list(frame["statistic"]) == ["contrast", "cosine", "scanvi_posterior"]
    # only the gated lane compared against anything
    assert list(frame["null_kind"]) == ["size_and_detection_matched", "", ""]
    # and only the lane with no null produced a probability
    assert list(frame["confidence"].notna()) == [False, False, True]
    assert (frame["status"] == AWARDED).all()


def test_records_to_frame_is_empty_not_absent():
    frame = records_to_frame([])
    assert len(frame) == 0
    assert list(frame.columns)[:4] == ["unit", "n", "label", "reason"]
