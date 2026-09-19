"""The islet chain: what may decide a cell in an islet, and what may not.

The risk this chain carries is not an arithmetic slip.  It is a criterion firing where it was
measured not to separate, a level invented for a stratum that could not supply one, or an
unmeasured channel read as a zero and so as evidence.  Each of those produces a table that looks
right and a number that is about the method rather than the tissue, so what is asserted here is
mostly about REFUSALS:

  * elongation must not be able to fire outside the large stratum -- barred twice, by the mask
    and by there being no level fitted below it;
  * a stratum that cannot supply a reference must report its level UNAVAILABLE, never zero,
    never the cohort's, never a neighbouring stratum's;
  * a hormone that was not measured cannot make a cell CONFIRMED endocrine, because the
    runner-up it would have to beat is unknown;
  * a cell with no nucleus measurement cannot clear a shape level;
  * exocrine inside an islet outline is an annotation error and keeps its reason in the row.

Nothing here asserts against a frozen copy of a result, and no test passes by reproducing a
biological finding.  The measured 20-donor numbers appear only as the REASON a rule exists, in
the modules' own comments.

Run with
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
      core/phenocycler/tests/current/test_islet_chain.py -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/celltype"))

import context as CX            # noqa: E402
import subpanel as SUB          # noqa: E402
import subpanel_io as S         # noqa: E402


# --------------------------------------------------------------------------------------------
# The size stratum.  One definition, and the large boundary is STRICTLY greater.
# --------------------------------------------------------------------------------------------

def test_the_large_stratum_is_strictly_more_than_ten_cells():
    out = S.size_stratum([9, 10, 11])
    assert list(out) == ["4-10", "4-10", S.LARGE_STRATUM], (
        "ten cells is not the large stratum -- islet_identity_mixture's --min-islet-cells is "
        "strictly greater and this must be the same boundary, not a second one")


def test_every_size_lands_in_exactly_one_named_stratum():
    sizes = np.arange(0, 60)
    out = S.size_stratum(sizes)
    assert set(out) <= set(S.SIZE_STRATA)
    assert list(S.size_stratum([0, 1])) == ["1", "1"], (
        "an islet that contributes no admitted cell has a size of zero, which is a size -- it "
        "belongs in the smallest stratum, not in a fifth one invented for it")
    assert list(S.size_stratum([2, 3])) == ["2-3", "2-3"]
    assert list(S.size_stratum([4, 7, 10])) == ["4-10"] * 3


def test_a_cell_in_no_islet_has_no_stratum():
    assert list(S.size_stratum([np.nan])) == [""], (
        "outside an islet there is no islet size, and an empty label is the honest answer -- "
        "not the smallest stratum, which would make it eligible for an islet criterion")


# --------------------------------------------------------------------------------------------
# The endocrine level: the mirror of the reference level, fitted per donor AND per stratum.
# --------------------------------------------------------------------------------------------

def test_the_endocrine_quantile_is_the_mirror_of_the_reference_quantile():
    assert S.ENDOCRINE_QUANTILE == 100.0 - S.REFERENCE_QUANTILE, (
        "the endocrine level must not be a second threshold chosen on its own -- it is the "
        "reference convention read in the other direction")


def test_the_level_is_the_low_tail_of_the_confirmed_cells():
    values = np.linspace(0.0, 1.0, 10_001)
    level = S.endocrine_level(values)
    assert level == pytest.approx(np.percentile(values, S.ENDOCRINE_QUANTILE))
    assert 0.0 < level < 0.5, "a level fitted from endocrine cells sits in their LOW tail"


def test_a_stratum_that_cannot_supply_a_reference_reports_unavailable():
    too_few = np.linspace(0.0, 1.0, S.MIN_REFERENCE_CELLS - 1)
    assert S.endocrine_level(too_few) is None, (
        "fewer confirmed cells than the guard means the level is unavailable -- returning a "
        "number anyway is how a stratum with no evidence comes to decide cells")


def test_the_level_ignores_cells_whose_measurement_is_missing():
    values = np.concatenate([np.linspace(0.2, 1.0, S.MIN_REFERENCE_CELLS),
                             np.full(5_000, np.nan)])
    level = S.endocrine_level(values)
    assert level is not None and level >= 0.2, (
        "NaN is an unavailable measurement; counting it as a low value would drag the level "
        "down and make the criterion pass cells it never measured")


def test_a_reference_of_only_missing_measurements_is_unavailable():
    assert S.endocrine_level(np.full(10_000, np.nan)) is None


def test_a_level_that_lands_on_the_scale_floor_is_refused():
    # `apply_scale` clips at the measured background, so a marker at or below background reads
    # exactly 0.0.  On donor 6534, 20.3 % of islet cells sit there and the 1st percentile of its
    # confirmed endocrine cells is 0.0 -- a "level" every measured cell in the islet clears.
    at_floor = np.concatenate([np.zeros(2_000), np.linspace(0.1, 0.9, 8_000)])
    assert S.endocrine_level(at_floor) is None, (
        "a level of 0.0 admits the entire islet while looking like a fitted number; it is the "
        "absence of a level and has to be reported as one")
    assert S.endocrine_level(at_floor, floor=float("-inf")) == pytest.approx(0.0), (
        "the diagnostic path must still be able to SHOW the degenerate value")


def test_the_elongation_floor_is_one_not_zero():
    # Nucleus elongation is max diameter over min diameter, so it cannot be below 1.0.
    barely = np.concatenate([np.ones(2_000), np.linspace(1.2, 3.0, 8_000)])
    assert S.endocrine_level(barely, floor=SUB.SHAPE_FLOOR) is None
    assert S.endocrine_level(barely) == pytest.approx(1.0), (
        "the marker floor of 0.0 would have let a vacuous elongation level through -- each "
        "feature is refused against its OWN floor")


# --------------------------------------------------------------------------------------------
# Fitting the levels: per stratum, and the shape level ONLY where shape was measured to separate.
# --------------------------------------------------------------------------------------------

def _pool(n_per_stratum: dict, cd99_by_stratum: dict, elong: float = 2.0):
    """A synthetic donor: n confirmed cells per stratum with a known CD99 distribution."""
    stratum, cd99, confirmed = [], [], []
    for name, n in n_per_stratum.items():
        stratum += [name] * n
        cd99 += list(np.linspace(*cd99_by_stratum[name], n))
        confirmed += [True] * n
    return (np.array(confirmed), np.array(stratum, dtype=object), np.array(cd99, dtype=float),
            np.full(len(cd99), elong, dtype=float))


def test_each_stratum_gets_its_own_cd99_level():
    """CORRECTED 2026-09-18. CD99 does not get a per-stratum level, and the reason is not a
    change of mind about the per-donor-per-stratum ruling -- it is what the ruling means once the
    ESTABLISHED estimator is used. A marker's reference population is that donor's ordinary
    tissue, which lies outside every islet and therefore carries no islet size, so there is one
    such population per donor and one level from it. `fit_islet_levels` is handed that level; it
    does not fit it. The predecessor fitted CD99 from a low quantile of the confirmed cells
    instead, which on 6534 returned exactly 0.0000 in all four strata because 20.3 % of its islet
    cells sit at the scale floor -- and `CD99 >= 0.0` admits every measured cell."""
    confirmed, stratum, cd99, elong = _pool(
        {"1": 2_000, "2-3": 2_000, "4-10": 2_000, ">10": 2_000},
        {s: (0.1, 0.9) for s in S.SIZE_STRATA})
    levels = SUB.fit_islet_levels(0.0518, confirmed, stratum, cd99, elong)
    assert levels["donor_reference_level"] == 0.0518, (
        "the donor's ordinary-tissue reference level is carried for comparison, not substituted")
    for name in S.SIZE_STRATA:
        assert levels[name]["CD99"] is not None, (
            f"{name}: per donor AND per size stratum -- a large islet's endocrine cells are the "
            f"only fair yardstick for a large islet's uncertain ones")


def test_the_shape_level_exists_only_in_the_large_stratum():
    confirmed, stratum, cd99, elong = _pool(
        {"1": 2_000, "2-3": 2_000, "4-10": 2_000, ">10": 2_000},
        {s: (0.1, 0.9) for s in S.SIZE_STRATA})
    levels = SUB.fit_islet_levels(0.05, confirmed, stratum, cd99, elong)
    for name in S.SIZE_STRATA:
        if name == S.LARGE_STRATUM:
            assert levels[name][SUB.SHAPE_FEATURE] is not None
        else:
            assert levels[name][SUB.SHAPE_FEATURE] is None, (
                f"{name}: no elongation level may be fitted where elongation was measured not "
                f"to separate -- the surest bar against it firing is that it does not exist")


def test_a_stratum_with_too_few_confirmed_cells_has_no_shape_level():
    confirmed, stratum, cd99, elong = _pool(
        {"1": 10, ">10": 10}, {"1": (0.1, 0.9), ">10": (0.1, 0.9)})
    levels = SUB.fit_islet_levels(0.05, confirmed, stratum, cd99, elong)
    assert levels[">10"][SUB.SHAPE_FEATURE] is None, (
        "ten confirmed cells cannot set a level; it is unavailable, never invented")
    assert levels[">10"]["n_confirmed"] == 10


def test_an_unavailable_donor_reference_level_is_carried_as_none():
    """A donor whose reference population cannot supply a level falls back to the hormone
    question. It must never read `None` as a level of zero."""
    confirmed, stratum, cd99, elong = _pool(
        {">10": 2_000}, {">10": (0.1, 0.9)})
    levels = SUB.fit_islet_levels(None, confirmed, stratum, cd99, elong)
    assert levels["donor_reference_level"] is None


def test_the_shape_feature_is_a_nucleus_feature():
    assert SUB.SHAPE_FEATURE.startswith("nucleus_"), (
        "cell-boundary shape is a tessellation artefact -- every class sits at cell circularity "
        "0.514-0.523 -- so the criterion has to read the nucleus, which is segmented from DAPI")


# --------------------------------------------------------------------------------------------
# The decision itself.
# --------------------------------------------------------------------------------------------

# ONE CD99 level for the donor -- the p99 of its ordinary tissue -- and a shape level only in
# the large stratum. Corrected 2026-09-18: CD99 used to be keyed per stratum here, which is the
# shape the broken estimator had.
LEVELS = {"donor_reference_level": 0.0518,
          "1": {"CD99": 0.20, SUB.SHAPE_FEATURE: None, "n_confirmed": 1_000},
          "2-3": {"CD99": 0.20, SUB.SHAPE_FEATURE: None, "n_confirmed": 1_000},
          "4-10": {"CD99": 0.20, SUB.SHAPE_FEATURE: None, "n_confirmed": 1_000},
          ">10": {"CD99": 0.20, SUB.SHAPE_FEATURE: 1.50, "n_confirmed": 1_000}}


def test_cd99_at_or_above_the_level_is_endocrine_and_below_it_is_not():
    stratum = np.array(["4-10", "4-10"], dtype=object)
    d = SUB.endocrine_decision(np.array([True, True]), stratum,
                               np.array([0.20, 0.19]), np.array([9.0, 9.0]), LEVELS,
                               np.array([True, True]))
    assert list(d["endocrine"]) == [True, False]
    assert list(d["cd99_ok"]) == [True, False]


def test_elongation_cannot_fire_below_the_large_stratum():
    """CORRECTED 2026-09-18, on the first run of islet_identity_mixture.py over 20 donors.

    Elongation was a second route that could call a cell endocrine where CD99 did not. The
    mixture test measured it at 3 of 16 fitted donors firing, median separation p 0.119 against
    a Bonferroni bar of 0.0083, and a median target separation of 2.194 against its own null of
    2.085 -- the second-weakest of the six features, where CD99 is the strongest (p_LR 0.0010,
    p_sep 0.0100). It is now computed and reported and decides nothing.
    """
    stratum = np.array(["4-10", ">10"], dtype=object)
    d = SUB.endocrine_decision(np.array([True, True]), stratum,
                               np.array([0.01, 0.01]), np.array([99.0, 99.0]), LEVELS,
                               np.array([True, True]))
    assert list(d["shape_ok"]) == [False, True]
    assert list(d["endocrine"]) == [False, True], (
        "below eleven cells per islet the shape was measured not to separate -- 1 of 20 donors "
        "in the one-cell stratum and 5 of 20 in 2-3 -- so it must not decide a cell there")


def test_elongation_only_adds_and_never_takes_endocrine_away():
    stratum = np.array([">10", ">10"], dtype=object)
    d = SUB.endocrine_decision(np.array([True, True]), stratum,
                               np.array([0.9, 0.9]),         # both clear CD99
                               np.array([9.0, 1.0]), LEVELS,  # one elongated, one round
                               np.array([True, True]))
    assert list(d["endocrine"]) == [True, True], (
        "a round cell that clears the key differentiator is still endocrine -- elongation is "
        "neither a requirement nor a route, so it can neither demote nor rescue")


def test_an_unfittable_stratum_never_reads_as_a_level_of_zero():
    """The level is one per donor now, so an unfittable one silences the islet question for the
    whole donor rather than for one stratum -- and it must be silence, not zero. `CD99 >= 0.0`
    would call every measured cell in every islet endocrine, which is exactly what the previous
    estimator did on 6534 before the floor guard refused it."""
    levels = dict(LEVELS, **{"1": {"CD99": None, SUB.SHAPE_FEATURE: None, "n_confirmed": 3}})
    stratum = np.array(["1", "1"], dtype=object)
    d = SUB.endocrine_decision(np.array([True, True]), stratum,
                               np.array([0.9, 0.0]), np.array([9.0, 9.0]), levels,
                               np.array([True, True]))
    assert not d["have_level"].any()
    assert not d["endocrine"].any(), (
        "no level means the question cannot be asked; a level read as zero would call every "
        "cell in the stratum endocrine, including one with no CD99 at all")


def test_a_cell_with_no_cd99_measurement_does_not_clear_the_level():
    stratum = np.array([">10"], dtype=object)
    d = SUB.endocrine_decision(np.array([True]), stratum,
                               np.array([np.nan]), np.array([1.0]), LEVELS,
                               np.array([True]))
    assert not d["cd99_ok"][0] and not d["endocrine"][0], (
        "an unavailable marker is not a low value and not a high one -- it is no evidence")


def test_a_cell_with_no_nucleus_cannot_clear_the_shape_level():
    stratum = np.array([">10"], dtype=object)
    d = SUB.endocrine_decision(np.array([True]), stratum,
                               np.array([0.01]), np.array([np.nan]), LEVELS,
                               np.array([True]))
    assert not d["shape_ok"][0] and not d["endocrine"][0], (
        "a cell with no nucleus is out of every calculation that needs one")


def test_a_cell_with_no_marker_measured_is_never_called_endocrine():
    # Donor 6534: 25,472 islet cells have no row in the aligned table at all, so all three
    # hormones are NaN.  The argmax that would NAME such a cell reads only those three channels,
    # so an unguarded chain names it by tie-break order -- Beta, every time.
    stratum = np.array([">10", ">10"], dtype=object)
    d = SUB.endocrine_decision(np.array([True, True]), stratum,
                               np.array([0.9, 0.9]), np.array([9.0, 9.0]), LEVELS,
                               np.array([True, False]))
    assert list(d["endocrine"]) == [True, False], (
        "a shape is not a name: whatever calls the cell endocrine, the hormone argmax has to "
        "have something to read before Beta, Alpha or Delta can be written on the row")
    assert list(d["have_level"]) == [True, False]


def test_elongation_alone_cannot_carry_an_unmeasured_cell():
    stratum = np.array([">10"], dtype=object)
    d = SUB.endocrine_decision(np.array([True]), stratum,
                               np.array([np.nan]), np.array([99.0]), LEVELS,
                               np.array([False]))
    assert not d["shape_ok"][0] and not d["endocrine"][0]


def test_the_exocrine_reason_is_only_written_where_exocrine_markers_were_read():
    cell_class = np.array(["Exocrine", "Exocrine", "Exocrine"], dtype=object)
    reason = SUB._finish(np.array([True, True, False]), cell_class,
                         np.array([True, False, False]))
    assert list(cell_class) == ["Excluded", "Excluded", "Exocrine"]
    assert list(reason) == [SUB.ISLET_EXOCRINE_REASON, SUB.NO_MEASUREMENT_REASON, ""], (
        "both reach Excluded by the same fall-through but for opposite causes -- exocrine "
        "evidence, and no evidence. Folding them together corrupts the count that is read as "
        "'how often the islet annotation caught exocrine tissue'")


def test_an_unmeasured_cell_that_the_chain_typed_keeps_its_class():
    cell_class = np.array(["Immune", "Beta"], dtype=object)
    reason = SUB._finish(np.array([True, True]), cell_class, np.array([False, False]))
    assert list(cell_class) == ["Immune", "Beta"] and list(reason) == ["", ""], (
        "only the REASON moves here. A cell the ordinary chain called Immune on markers that "
        "WERE measured keeps that call -- the non-endocrine cells inside the islet are the "
        "analysis, not a bucket to empty")


def test_a_cell_with_no_hormone_but_real_immune_markers_keeps_its_immune_call():
    # 1,163 cells of donor 6534 are in the aligned table with all three hormones missing and
    # real supporting-tier markers.  Ruling those out with the genuinely unmeasured ones would
    # delete immune cells from inside islets, which is the analysis.
    cell_class = np.array(["Immune"], dtype=object)
    reason = SUB._finish(np.array([True]), cell_class, np.array([True]))
    assert list(cell_class) == ["Immune"] and list(reason) == [""]
    values = pd.DataFrame({"INS": [np.nan, 0.8], "GCG": [np.nan, 0.0], "SST": [np.nan, 0.0],
                           "CD163": [0.95, 0.0]})
    assert list(SUB.hormone_evidence(values)["measured"]) == [False, True], (
        "no hormone was read on the first cell, so no hormone can name it -- but a marker WAS "
        "read, so it is not an unmeasured cell")
    assert list(values.notna().any(axis=1)) == [True, True], (
        "the exclusion reason keys on whether ANY marker was read, not on the hormones")


def test_a_cell_outside_an_islet_is_never_decided_by_this_chain():
    stratum = np.array([">10"], dtype=object)
    d = SUB.endocrine_decision(np.array([False]), stratum,
                               np.array([0.99]), np.array([9.0]), LEVELS,
                               np.array([True]))
    assert not d["endocrine"][0] and not d["have_level"][0]


# --------------------------------------------------------------------------------------------
# The hormone evidence, and what "confirmed" refuses to be built from.
# --------------------------------------------------------------------------------------------

def test_a_cell_missing_one_hormone_is_never_confirmed_endocrine():
    values = pd.DataFrame({"INS": [0.90, 0.90], "GCG": [0.01, np.nan], "SST": [0.01, 0.01]})
    out = SUB.hormone_evidence(values)
    assert list(out["confirmed"]) == [True, False], (
        "the runner-up an unmeasured channel would have had is unknown, so the margin cannot "
        "be shown -- treating it as a zero runner-up manufactures confidence from an absence")
    assert list(out["complete"]) == [True, False]


def test_the_confirmed_definition_is_islet_identity_s_and_not_a_second_one():
    import islet_identity as ID
    values = pd.DataFrame({"INS": [ID.REFERENCE_MIN, ID.REFERENCE_MIN - 1e-6],
                           "GCG": [0.0, 0.0], "SST": [0.0, 0.0]})
    assert list(SUB.hormone_evidence(values)["confirmed"]) == [True, False]
    margin = pd.DataFrame({"INS": [0.90, 0.90], "GCG": [0.90 / ID.REFERENCE_MARGIN, 0.60],
                           "SST": [0.0, 0.0]})
    assert list(SUB.hormone_evidence(margin)["confirmed"]) == [True, False]


def test_the_argmax_names_the_hormone_and_a_donor_without_one_is_refused():
    values = pd.DataFrame({"INS": [0.1, 0.9], "GCG": [0.9, 0.1], "SST": [0.0, 0.0]})
    assert list(SUB.hormone_evidence(values)["which"]) == ["GCG", "INS"]
    with pytest.raises(ValueError, match="no hormone channel"):
        SUB.hormone_evidence(pd.DataFrame({"INS": [np.nan], "GCG": [np.nan], "SST": [np.nan]}))


# --------------------------------------------------------------------------------------------
# Exocrine inside an islet outline, and the reason staying in the row.
# --------------------------------------------------------------------------------------------

def test_exocrine_inside_an_islet_is_excluded_and_says_why():
    cell_class = np.array(["Exocrine", "Exocrine", "Immune"], dtype=object)
    reason = SUB._finish(np.array([True, False, True]), cell_class)
    assert list(cell_class) == ["Excluded", "Exocrine", "Immune"], (
        "exocrine markers get into an islet annotation only by mistake, so the cell is ruled "
        "out of analysis rather than counted as an exocrine cell inside an islet")
    assert list(reason) == [SUB.ISLET_EXOCRINE_REASON, "", ""]


def test_the_non_endocrine_islet_population_is_typed_rather_than_discarded():
    cell_class = np.array(["Immune", "Blood vessel", "Nerve"], dtype=object)
    reason = SUB._finish(np.array([True, True, True]), cell_class)
    assert list(cell_class) == ["Immune", "Blood vessel", "Nerve"]
    assert not any(reason), (
        "the non-endocrine cells INSIDE the islet are the analysis -- only the ones the chain "
        "can say nothing about beyond Exocrine are ruled out")


# --------------------------------------------------------------------------------------------
# Object membership recorded beside the context, not instead of it.
# --------------------------------------------------------------------------------------------

def test_object_label_does_not_consult_islets_and_assign_overwrites_it():
    calls = {}

    def fake(x, y, stem, pixel):
        calls["n"] = calls.get("n", 0) + 1
        return {"Blood vessel": np.array([True, False]),
                "Nerve": np.array([False, False]),
                "Lymphatic": np.array([False, True])}

    original = CX.memberships
    CX.memberships = fake
    try:
        label = CX.object_label(np.zeros(2), np.zeros(2), "stem", 1.0)
        context = CX.assign(np.zeros(2), np.zeros(2), np.array([True, False]), "stem", 1.0)
    finally:
        CX.memberships = original

    assert list(label) == ["Blood vessel", "Lymphatic"], (
        "the object answer must survive on its own: a cell in a vessel INSIDE an islet is the "
        "case an earlier analysis of this lane could not see")
    assert list(context) == ["islet", "Lymphatic"]


def test_the_emitted_table_carries_the_object_answer_without_a_flag():
    ctx = pd.DataFrame({"object_id": ["a", "b"], "context": ["islet", "tissue"],
                        "object_context": ["Blood vessel", "tissue"],
                        "cell_region": ["core", "none"], "islet_num": ["7", None]})
    frame = SUB._emit("115", ctx,
                      np.array(["Beta", "Exocrine"], dtype=object),
                      np.array(["", ""], dtype=object),
                      np.array(["cd99", "fallthrough"], dtype=object),
                      np.array(["", ""], dtype=object),
                      np.array([42.0, np.nan]), np.array([">10", ""], dtype=object),
                      np.array(["cd99", ""], dtype=object),
                      np.array(["INS", "Ker8_18"], dtype=object), {}, detail=False)
    for column in ("object_context", "excluded_reason", "islet_num", "islet_cells",
                   "islet_stratum", "islet_route"):
        assert column in frame.columns, (
            f"{column} has to be on every row -- a column that exists only under --emit-detail "
            f"leaves the gap that misled the previous analysis open")
    assert frame["object_context"].tolist() == ["Blood vessel", "tissue"]


# --------------------------------------------------------------------------------------------
# The chain's own wiring: the baseline pass and the islet pass differ only where they should.
# --------------------------------------------------------------------------------------------

def test_the_islet_answer_never_overrules_a_structure():
    values = pd.DataFrame({"CD4": [0.0], "CD68": [0.0]})
    structure_class = np.array(["Blood vessel"], dtype=object)
    cell_class, _, resolved, reason = SUB._run_chain(
        values, {"bottom_tier_immune": [], "subpanels": {}}, {}, structure_class,
        np.array(["structure"], dtype=object), np.array([True]),
        np.array([True]), np.array(["Beta"], dtype=object),
        np.array(["cd99"], dtype=object))
    assert cell_class[0] == "Blood vessel" and resolved[0] == "structure", (
        "the classified objects separate vessels and nerves, and they do so inside an islet "
        "too -- the structure question is asked first and its answer stands")
    assert reason[0] == ""


def test_a_cell_the_islet_chain_declines_goes_down_the_ordinary_chain():
    values = pd.DataFrame({"CD4": [0.0], "CD68": [0.0]})
    cell_class, _, resolved, reason = SUB._run_chain(
        values, {"bottom_tier_immune": [], "subpanels": {}}, {},
        np.array([None], dtype=object), np.array([""], dtype=object), np.array([True]),
        np.array([False]), np.array(["Beta"], dtype=object),
        np.array([""], dtype=object))
    assert cell_class[0] == "Excluded" and reason[0] == SUB.ISLET_EXOCRINE_REASON
    assert resolved[0] == "fallthrough"
