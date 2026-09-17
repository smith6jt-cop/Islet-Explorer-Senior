"""The two fixes of 2026-09-05 to the rule reader and the gate, on tables whose answer is known.

1. ``rules.requirements`` returns the negative constraints (``none_of``) when asked, and
   ``score.meets_requirements`` enforces them -- before, a negative entered only the signed
   cosine as a direction, so a cell bright in a forbidden marker could still win a reference.
2. ``rules.reference_matrix`` keys an OR-group by ``(declaring row, token)`` -- before, a
   leaf inheriting ``Immune``'s ``anypos`` group and declaring its own merged the two into one
   group of n1 + n2 members, each weighed 1/(n1 + n2).

Nothing here touches the cohort.  Tables are built as ``parent, phenotype, <registry marker>``
frames, as the older tests do.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import rules as R                              # noqa: E402
import score as SC                             # noqa: E402


def _epithelial_table() -> pd.DataFrame:
    """Epithelial declares two epithelial OR-members and two immune negatives; Endo adds a
    hormone; Odd OVERRIDES an inherited negative, which must then leave its none_of."""
    return pd.DataFrame({
        "parent": ["all", "Epithelial", "Epithelial"],
        "phenotype": ["Epithelial", "Endo", "Odd"],
        "E_cadherin": ["anypos", None, None],
        "Ker8_18": ["anypos", None, None],
        "CD68": ["allneg", None, "allpos"],
        "CD3e": ["allneg", None, None],
        "INS": [None, "allpos", None],
    })


# ------------------------------------------------------------------ fix 1: negatives
def test_requirements_default_shape_is_unchanged():
    """Every existing caller unpacks two values; the default must keep giving two."""
    rq = R.requirements(_epithelial_table())
    assert all(len(v) == 2 for v in rq.values())
    assert rq["Endo"] == (["INS"], [["E_cadherin", "Ker8_18"]])


def test_none_of_is_inherited_and_the_child_override_wins():
    rq = R.requirements_with_negatives(_epithelial_table())
    assert rq["Epithelial"] == ([], [["E_cadherin", "Ker8_18"]], ["CD3e", "CD68"])
    assert rq["Endo"] == (["INS"], [["E_cadherin", "Ker8_18"]], ["CD3e", "CD68"]), \
        "an ancestor's allneg is a constraint on every leaf"
    assert rq["Odd"] == (["CD68"], [["E_cadherin", "Ker8_18"]], ["CD3e"]), \
        "a child committing the marker to allpos removes it from none_of"
    assert R.requirements(_epithelial_table(), with_negatives=True) == rq


def test_anyneg_is_folded_into_none_of_and_recorded_as_such():
    t = pd.DataFrame({"parent": ["all"], "phenotype": ["A"],
                      "CD3e": ["allpos"], "CD68": ["anyneg"], "CD20": ["anyneg"], "CD4": ["neg"]})
    assert R.requirements_with_negatives(t)["A"] == (["CD3e"], [], ["CD20", "CD4", "CD68"])


def test_the_gate_refuses_a_cell_bright_in_a_forbidden_marker():
    X = np.array([[1.0, 0.0], [1.0, 0.5], [1.0, 0.09], [1.0, 0.10]], dtype=np.float32)
    reqs = {"A": (["m0"], [], ["m1"])}
    allowed = SC.meets_requirements(X, ["m0", "m1"], reqs, ["A"], 0.10)
    assert allowed[:, 0].tolist() == [True, False, True, False], \
        "at or above require_min on a none_of marker fails; below it passes"


def test_a_two_tuple_still_enforces_no_negative():
    X = np.array([[1.0, 0.5]], dtype=np.float32)
    assert SC.meets_requirements(X, ["m0", "m1"], {"A": (["m0"], [])}, ["A"], 0.10)[0, 0]


def test_an_unmeasured_or_nan_negative_fails_closed():
    """Absence cannot be verified on a marker nobody measured; a missing measurement is
    unavailable, not negative."""
    X = np.array([[1.0, 0.0], [1.0, np.nan]], dtype=np.float32)
    unmeasured = SC.meets_requirements(X, ["m0", "m1"], {"A": (["m0"], [], ["m9"])}, ["A"], 0.10)
    assert not unmeasured.any()
    nan_cell = SC.meets_requirements(X, ["m0", "m1"], {"A": (["m0"], [], ["m1"])}, ["A"], 0.10)
    assert nan_cell[:, 0].tolist() == [True, False]


def test_the_negative_gates_the_call_not_just_the_angle():
    """A is 'm0 and not m1'; B is 'm1'.  A cell bright in both would have won A on angle with
    the old gate (m0 is shown); the negative now refuses A and the cell goes to B.  Where the
    two references forbid each other, the cell meets nothing and is Ambiguous with a reason."""
    markers = ["m0", "m1"]
    X = np.array([[1.0, 1.0]], dtype=np.float32)
    Rm = np.array([[1.0, 0.0], [-1.0, 1.0]], dtype=np.float32)          # signed encoding
    reqs = {"A": (["m0"], [], ["m1"]), "B": (["m1"], [], [])}
    allowed = SC.meets_requirements(X, markers, reqs, ["A", "B"], 0.10)
    out = SC.call_from_scores(SC.cosine_scores(X, Rm), ["A", "B"], allowed=allowed)
    assert out.loc[0, "call"] == "B"
    reqs_mutual = {"A": (["m0"], [], ["m1"]), "B": (["m1"], [], ["m0"])}
    allowed = SC.meets_requirements(X, markers, reqs_mutual, ["A", "B"], 0.10)
    out = SC.call_from_scores(SC.cosine_scores(X, Rm), ["A", "B"], allowed=allowed)
    assert out.loc[0, "call"] == "Ambiguous" and out.loc[0, "reason"] == "no_reference_requirements_met"


def test_score_and_gate_carries_the_negatives_through():
    X = np.array([[1.0, 1.0], [1.0, 0.0]], dtype=np.float32)
    Rm = np.array([[1.0], [-1.0]], dtype=np.float32)
    S, allowed = SC.score_and_gate(X, None, Rm, ["m0", "m1"], {"A": (["m0"], [], ["m1"])},
                                   ["A"], 0.10, np.zeros(1, dtype=bool))
    assert allowed[:, 0].tolist() == [False, True]


# ------------------------------------------------------------------ fix 2: group keys
def _two_level_table() -> pd.DataFrame:
    return pd.DataFrame({
        "parent": ["all", "Immune"], "phenotype": ["Immune", "Tsub"],
        "CD68": ["anypos", None], "CD20": ["anypos", None], "CD3e": ["anypos", None],
        "CD4": [None, "anypos"], "CD8": [None, "anypos"],
    })


MARKERS = ["CD68", "CD20", "CD3e", "CD4", "CD8"]


def test_flatten_carries_the_declaring_row():
    flat = R.flatten(_two_level_table())
    assert isinstance(flat, R.FlatRules) and isinstance(flat, dict)
    assert flat["Tsub"] == {"CD68": "anypos", "CD20": "anypos", "CD3e": "anypos",
                            "CD4": "anypos", "CD8": "anypos"}, "the plain mapping is unchanged"
    assert flat.source["Tsub"] == {"CD68": "Immune", "CD20": "Immune", "CD3e": "Immune",
                                   "CD4": "Tsub", "CD8": "Tsub"}


def test_or_groups_declared_at_different_levels_keep_their_own_weight():
    """The leaf inherits a 3-member group from Immune and declares a 2-member group of its own
    under the same token.  Keyed by token they merged: every member 1/5.  Keyed by declaring
    row: 1/3 and 1/2, so the leaf's own commitment is no longer diluted by its ancestor's."""
    flat = R.flatten(_two_level_table())
    M, dropped = R.reference_matrix(flat, MARKERS, "zero")
    assert not dropped
    assert M.loc[["CD68", "CD20", "CD3e"], "Tsub"].tolist() == pytest.approx([1 / 3] * 3)
    assert M.loc[["CD4", "CD8"], "Tsub"].tolist() == pytest.approx([1 / 2] * 2)
    assert M.loc[["CD68", "CD20", "CD3e"], "Immune"].tolist() == pytest.approx([1 / 3] * 3)
    assert (M.loc[["CD4", "CD8"], "Immune"] == 0).all()

    # the with_source shape gives the same matrix
    Ms, _ = R.reference_matrix(R.flatten(_two_level_table(), with_source=True), MARKERS, "zero")
    assert np.allclose(Ms.to_numpy(), M.to_numpy())

    # what it was, and what a dict without provenance still is: one level, one group of five
    plain = {ph: dict(ref) for ph, ref in flat.items()}
    Mp, _ = R.reference_matrix(plain, MARKERS, "zero")
    assert Mp.loc[MARKERS, "Tsub"].tolist() == pytest.approx([1 / 5] * 5)


def test_a_one_level_table_is_unchanged_by_the_new_key():
    t = pd.DataFrame({"parent": ["all"], "phenotype": ["A"],
                      "CD68": ["anypos"], "CD20": ["anypos"], "CD3e": ["allneg"]})
    M, _ = R.reference_matrix(R.flatten(t), ["CD68", "CD20", "CD3e"], "signed")
    assert M["A"].tolist() == pytest.approx([0.5, 0.5, -1.0])


def test_the_count_is_over_the_members_that_survive_in_the_leaf():
    """Recorded consequence, not a new decision: the counting basis is unchanged.  Epi-endocrine
    declares 'any of B3TUBB|GCG|INS|SST'; Alpha commits GCG to allpos and INS/SST to allneg, so
    B3TUBB is the only member of that group left in the flattened leaf and carries its unit
    weight alone.  Under the old token key it shared the token with Epithelial's two members and
    weighed 1/3."""
    t = pd.DataFrame({
        "parent": ["all", "Epithelial", "Epi-endocrine"],
        "phenotype": ["Epithelial", "Epi-endocrine", "Alpha"],
        "E_cadherin": ["anypos", None, None], "Ker8_18": ["anypos", None, None],
        "B3TUBB": [None, "anypos", None], "GCG": [None, "anypos", "allpos"],
        "INS": [None, "anypos", "allneg"], "SST": [None, "anypos", "allneg"],
    })
    mk = ["E_cadherin", "Ker8_18", "B3TUBB", "GCG", "INS", "SST"]
    M, _ = R.reference_matrix(R.flatten(t), mk, "signed")
    assert M.loc[["E_cadherin", "Ker8_18"], "Alpha"].tolist() == pytest.approx([0.5, 0.5])
    assert M.loc["B3TUBB", "Alpha"] == pytest.approx(1.0)
    assert M.loc["GCG", "Alpha"] == 1.0 and M.loc[["INS", "SST"], "Alpha"].tolist() == [-1.0, -1.0]
    assert M.loc[["B3TUBB", "GCG", "INS", "SST"], "Epi-endocrine"].tolist() == pytest.approx([0.25] * 4)
