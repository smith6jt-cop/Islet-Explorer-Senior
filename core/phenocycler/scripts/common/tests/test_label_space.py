"""The label space must be auditable and must fail closed.

The danger this file guards is the one the crosswalk's own policy names: a restriction
that quietly improves a result.  So the tests assert the *algebra* (renormalizing equals
masking the logits, which is what makes this a re-read of the model's own output rather
than a new decision), the *provenance* (every excluded class carries a written reason),
and the *failure* behaviour (bad input raises; it never returns a confident answer).
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from label_space import AMBIGUOUS, load_label_space  # noqa: E402

CROSSWALK = json.loads((HERE.parent / "dct" / "broad_class_crosswalk.json").read_text())  # DCT artifact
CLASSES = tuple(CROSSWALK["entries"])


def _space():
    return load_label_space(CLASSES)


def _softmax(z):
    e = np.exp(z - z.max(axis=1, keepdims=True))
    return e / e.sum(axis=1, keepdims=True)


# ---------------------------------------------------------------- the algebra
def test_restriction_equals_masking_the_logits():
    """Renormalizing the softmax over kept columns == setting excluded logits to -inf.

    This is the whole justification for doing it as post-processing.  If it ever stops
    holding, the restriction has become a new decision and needs re-arguing, not a fix.
    """
    ls = _space()
    rng = np.random.default_rng(0)
    z = rng.normal(size=(500, len(CLASSES))) * 3.0
    P = _softmax(z)

    masked = z.copy()
    drop = [CLASSES.index(c) for c in ls.excluded]
    masked[:, drop] = -np.inf
    expect = _softmax(masked)[:, [CLASSES.index(c) for c in ls.kept]]

    got = ls.restrict(P, CLASSES).probabilities
    np.testing.assert_allclose(got, expect, rtol=0, atol=1e-12)


def test_restricted_rows_are_probability_vectors():
    ls = _space()
    rng = np.random.default_rng(1)
    P = _softmax(rng.normal(size=(200, len(CLASSES))))
    R = ls.restrict(P, CLASSES)
    np.testing.assert_allclose(R.probabilities.sum(axis=1), 1.0, atol=1e-12)
    assert R.probabilities.shape[1] == len(ls.kept)


def test_impossible_mass_is_the_removed_probability():
    ls = _space()
    rng = np.random.default_rng(2)
    P = _softmax(rng.normal(size=(200, len(CLASSES))))
    R = ls.restrict(P, CLASSES)
    drop = [CLASSES.index(c) for c in ls.excluded]
    np.testing.assert_allclose(R.impossible_mass, P[:, drop].sum(axis=1), atol=1e-12)


def test_unrestricted_call_is_untouched():
    """The restriction reports what it changed; it does not overwrite the original."""
    ls = _space()
    rng = np.random.default_rng(3)
    P = _softmax(rng.normal(size=(300, len(CLASSES))) * 4)
    R = ls.restrict(P, CLASSES)
    np.testing.assert_array_equal(R.unrestricted_call, np.asarray(CLASSES)[P.argmax(axis=1)])


# ---------------------------------------------------------------- fail closed
def test_non_finite_raises_rather_than_renormalizing():
    ls = _space()
    P = np.full((3, len(CLASSES)), 1.0 / len(CLASSES))
    P[1, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        ls.restrict(P, CLASSES)


def test_rows_that_are_not_probabilities_raise():
    ls = _space()
    P = np.full((3, len(CLASSES)), 1.0 / len(CLASSES))
    P[2] *= 0.5
    with pytest.raises(ValueError, match="do not sum to 1"):
        ls.restrict(P, CLASSES)


def test_wrong_vocabulary_raises():
    ls = _space()
    wrong = ("NotAClass",) + CLASSES[1:]
    P = np.full((2, len(CLASSES)), 1.0 / len(CLASSES))
    with pytest.raises(ValueError, match="not the vocabulary"):
        ls.restrict(P, wrong)


def test_cell_with_no_surviving_mass_is_ambiguous_not_called():
    """All the probability on impossible classes must not become a confident call."""
    ls = _space()
    P = np.zeros((1, len(CLASSES)))
    P[0, CLASSES.index(ls.excluded[0])] = 1.0
    R = ls.restrict(P, CLASSES)
    assert R.restricted_call[0] == AMBIGUOUS
    assert R.ambiguous_reason[0] == "no_mass_in_label_space"
    assert np.isnan(R.probabilities[0]).all()


def test_removed_mass_outweighing_the_survivor_is_ambiguous():
    ls = _space()
    P = np.zeros((1, len(CLASSES)))
    P[0, CLASSES.index(ls.excluded[0])] = 0.7   # impossible
    P[0, CLASSES.index(ls.kept[0])] = 0.3       # best survivor
    R = ls.restrict(P, CLASSES)
    assert R.restricted_call[0] == AMBIGUOUS
    assert R.ambiguous_reason[0] == "mass_on_impossible_classes"
    # ...and the survivor is still reported, so nothing is lost
    assert R.probabilities[0, 0] == pytest.approx(1.0)


def test_a_class_with_no_decision_is_not_a_pass(tmp_path):
    d = json.loads(json.dumps(CROSSWALK))
    del d["entries"]["Tcell"]["in_label_space"]
    f = tmp_path / "cw.json"
    f.write_text(json.dumps(d))
    with pytest.raises(ValueError, match="no in_label_space decision"):
        load_label_space(CLASSES, f)


def test_an_exclusion_without_a_reason_is_rejected(tmp_path):
    d = json.loads(json.dumps(CROSSWALK))
    ex = next(c for c, v in d["entries"].items() if not v["in_label_space"])
    d["entries"][ex].pop("reason", None)
    d["entries"][ex].pop("panel_reason", None)
    f = tmp_path / "cw.json"
    f.write_text(json.dumps(d))
    with pytest.raises(ValueError, match="excluded with no written reason"):
        load_label_space(CLASSES, f)


# ------------------------------------------------------- the decisions themselves
def test_every_class_carries_a_decision_and_every_exclusion_a_reason():
    for c, v in CROSSWALK["entries"].items():
        assert isinstance(v["in_label_space"], bool), c
        if not v["in_label_space"]:
            assert v["exclusion_reason"] in {
                "class_not_present_in_tissue", "class_unsupported_by_panel"}, c
            assert (v.get("reason") or v.get("panel_reason")), c


def test_the_pancreas_specific_calls_are_the_right_way_round():
    """Two traps this project has already tripped over, pinned as tests.

    Stellate IS the pancreatic stellate cell and must survive; Langerhans is the
    EPIDERMAL dendritic cell, not the islet of Langerhans, and must not.
    """
    E = CROSSWALK["entries"]
    assert E["Stellate"]["in_label_space"] is True
    assert E["Langerhans"]["in_label_space"] is False
    for endocrine in ("AlphaCell", "BetaCell", "Endocrine"):
        assert E[endocrine]["in_label_space"] is True, endocrine


def test_label_space_size_is_pinned():
    """A silent change to the label space would move every gate number in the lane."""
    ls = _space()
    assert len(ls.kept) == 26
    assert len(ls.excluded) == 25
    assert CROSSWALK["method_version"] == "dct-broad-class-2026-08-24"


def test_broad_groups_cover_the_kept_classes_exactly():
    ls = _space()
    flat = [c for v in ls.broad_groups().values() for c in v]
    assert sorted(flat) == sorted(ls.kept)
    assert len(flat) == len(set(flat)), "a class landed in two compartments"
