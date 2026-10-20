"""Missing measurements cannot supply a required positive or a required negative.

These assertions distinguish an observed OR witness from an unmeasured alternative,
and keep legitimate broad supporting-tier calls separate from specific definitions.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/celltype"))
import subpanel as SUB  # noqa: E402


def test_a_required_negative_must_be_measured_and_below_the_existing_level():
    values = pd.DataFrame({"positive": [.2] * 7,
                           "negative": [np.nan, -np.inf, np.inf, .099, .1, .2, 0.]})
    assert SUB.fires(values, (["positive"], [], ["negative"])).tolist() == [
        False, False, False, True, False, False, True]


@pytest.mark.parametrize("rule", [(["missing"], [], []),
                                  ([], [["missing", "also_missing"]], []),
                                  ([], [], ["missing"])])
def test_a_wholly_unmeasured_required_term_does_not_fire(rule):
    values = pd.DataFrame({"unrelated": [.9, .0]}, index=["key-b", "key-a"])
    assert SUB.fires(values, rule).tolist() == [False, False]


def test_an_observed_or_witness_is_sufficient_with_an_unmeasured_alternative():
    values = pd.DataFrame({"witness": [.1, .099, np.nan, np.inf],
                           "alternative": [np.nan, np.nan, .3, np.nan]})
    before = values.copy(deep=True)
    assert SUB.fires(values, ([], [["witness", "alternative", "absent"]], [])).tolist() == [
        True, False, True, False]
    pd.testing.assert_frame_equal(values, before)


def test_positive_terms_also_require_finite_measurements():
    values = pd.DataFrame({"positive": [np.inf, np.nan, -np.inf, .1]})
    assert SUB.fires(values, (["positive"], [], [])).tolist() == [False, False, False, True]


def test_rule_unknown_is_separate_from_a_measured_contradiction():
    values = pd.DataFrame({"positive": [.2, .01, .2, .2],
                           "negative": [np.nan, np.nan, .2, 0.]})
    fired, unavailable = SUB.rule_evidence(values, (["positive"], [], ["negative"]))
    assert fired.tolist() == [False, False, False, True]
    assert unavailable.tolist() == [True, False, False, False]


def test_or_unknown_requires_no_positive_witness_and_an_unmeasured_alternative():
    values = pd.DataFrame({"a": [.2, .01, np.nan], "b": [np.nan, .01, np.nan]})
    fired, unavailable = SUB.rule_evidence(values, ([], [["a", "b"]], []))
    assert fired.tolist() == [True, False, False]
    assert unavailable.tolist() == [False, False, True]


def test_missing_negative_withholds_definition_but_keeps_measured_supporting_call():
    values = pd.DataFrame({"CD3e": [.95], "CD38": [np.nan], "E_cadherin": [0.]},
                          index=["original-cell-key"])
    classes = np.array([None], object)
    definitions = np.array([""], object)
    routes = np.array([""], object)
    SUB._tissue_chain(values, {"bottom_tier_immune": ["CD3e"]},
                      {"Imm-fixture": (["CD3e"], [], ["CD38"])},
                      classes, definitions, routes)
    assert classes.tolist() == ["Immune"]
    assert definitions.tolist() == [""]
    assert routes.tolist() == ["supporting"]
    assert values.index.tolist() == ["original-cell-key"]


def test_unavailable_rule_evidence_is_audited_without_inventing_a_type():
    values = pd.DataFrame({"CD3e": [.95, .01], "CD38": [np.nan, np.nan],
                           "E_cadherin": [0., .95]})
    classes = np.array([None, None], object)
    definitions = np.array(["", ""], object)
    routes = np.array(["", ""], object)
    audit = {}
    unavailable_count = np.zeros(2, dtype=np.int16)
    SUB._tissue_chain(values, {"bottom_tier_immune": ["CD3e"]},
                      {"Imm-fixture": (["CD3e"], [], ["CD38"])},
                      classes, definitions, routes, audit, unavailable_count)
    assert classes.tolist() == ["Immune", "Exocrine"]
    assert routes.tolist() == ["supporting", "fallthrough"]
    assert unavailable_count.tolist() == [1, 0]
    assert audit["Imm-fixture"]["n_unavailable"] == 1
    assert audit["Imm-fixture"]["n_firing"] == 0
