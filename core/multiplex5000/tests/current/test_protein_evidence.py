"""The review-only protein cutoff uses one physical section at a time."""

import numpy as np
import pytest

from multiplex5000.protein_evidence import ProteinEvidenceError, candidate_calls


def _section():
    background = np.linspace(5.0, 6.0, 1200)
    values = np.r_[background, np.full(20, 100.0), 0.0, np.nan, -1.0, 100.0]
    eligible = np.ones(len(values), dtype=bool)
    eligible[-1] = False
    return values, eligible


def test_candidate_cutoff_is_section_local_and_preserves_missing():
    values, eligible = _section()
    first = candidate_calls("HDL090-S1", "CD8", values, eligible, shift_scale=1.0)
    second = candidate_calls("HDL082-S2", "CD8", values * 10, eligible, shift_scale=1.0)
    assert first.unavailable_reason is None
    assert second.unavailable_reason is None
    assert 5 < first.floor < 7
    assert first.top == pytest.approx(100.0)
    assert first.cutoff == pytest.approx(first.floor + 0.10 * (first.top - first.floor))
    assert second.cutoff == pytest.approx(first.cutoff * 10, rel=0.002)
    assert first.states.tolist() == second.states.tolist()
    assert first.states[0] == 0
    assert first.states[1200] == 1
    assert first.states[-4:].tolist() == [0, -1, -1, -1]


def test_missing_marker_shift_or_background_are_unavailable():
    values, eligible = _section()
    absent = candidate_calls("HDL090-S1", "CD8", None, eligible, shift_scale=1.0)
    no_shift = candidate_calls("HDL090-S1", "CD8", values, eligible, shift_scale=None)
    few = candidate_calls("HDL090-S1", "CD8", values[:999], eligible[:999],
                          shift_scale=1.0)
    assert absent.unavailable_reason == "marker_not_measured"
    assert no_shift.unavailable_reason == "uniform_shift_unavailable"
    assert few.unavailable_reason == "insufficient_background_support"
    assert np.all(absent.states == -1)
    assert np.all(no_shift.states == -1)
    assert np.all(few.states == -1)


def test_uniform_shift_changes_units_without_a_fixed_positive_fraction():
    values, eligible = _section()
    first = candidate_calls("HDL090-S1", "CD8", values, eligible, shift_scale=1.0)
    shifted = candidate_calls("HDL090-S1", "CD8", values, eligible, shift_scale=2.0)
    assert shifted.cutoff == pytest.approx(first.cutoff * 2, rel=0.005)
    assert shifted.states.tolist() == first.states.tolist()
    assert int(np.count_nonzero(first.states == 1)) == 20


def test_input_identity_and_scale_must_be_valid():
    values, eligible = _section()
    with pytest.raises(ProteinEvidenceError, match="cell order"):
        candidate_calls("HDL090-S1", "CD8", values[:-1], eligible, shift_scale=1.0)
    with pytest.raises(ProteinEvidenceError, match="shift"):
        candidate_calls("HDL090-S1", "CD8", values, eligible, shift_scale=0.0)
