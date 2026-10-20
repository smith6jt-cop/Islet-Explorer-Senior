"""The structure call: centroid, then MAJORITY OF CELL AREA, then the object's own marker.

The maintainer's rule, stated in their own words and restated five times before it was
implemented:

    FOR ANY CELL WITH ITS CENTROID EXACTLY WITHIN THE BOUNDS OF THE OBJECT, THE ENTIRE CELL GETS
    CHECKED AND IF THE MAJORITY OF THE CELL BOUNDARIES ARE WITHIN THE OBJECT IT IS THEN A
    CANDIDATE TO BE CLASSIFIED ACCORDING TO WHAT THE OBJECT IS. IF THE SAME MARKER THAT CREATED
    THE OBJECT IS SIMPLY HIGHER THAN ALL THE REST, IT GETS CLASSIFIED. IF NOT, IT GETS A
    "PASSTHROUGH" FEATURE INDICATING THE FRACTION OF THE OBJECT THAT OVERLAPS THE CELL AREA.

What is asserted here is mostly REFUSALS, because the ways this goes wrong all look like success:

  * a centroid inside an object must NOT be candidacy on its own -- that was the behaviour until
    2026-09-19 and it called cells a vessel merely ran past;
  * a cell whose overlap was never measured is UNAVAILABLE: not a candidate, and not counted as a
    measured minority either;
  * the passthrough carries the FRACTION, not a flag, because being crossed is a measurement;
  * a SMALL number of calls is the expected result and no test may assert a floor on it -- most
    vessels, lymphatics and nerves are thin processes running past cells rather than cell bodies.

No test here reproduces a biological result, and none asserts a count from a build.

Run with
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
      core/multiplex5000/tests/current/test_structure_rule.py -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/typing"))
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/celltype"))

import subpanel as SUB  # noqa: E402

VESSEL = "Blood vessel"
#: One marker per object, as the rule has it: the marker the pixel classifier was BUILT ON.
CLF = "CD34"


def _log():
    """The minimum `_structure_calls` reads: which markers answer for each structure."""
    return {"subpanels": {VESSEL: [CLF], "Nerve": ["B3TUBB"], "Lymphatic": ["Podoplanin"]}}


def _frame(n, frac, clf_value=1.0, other=0.0):
    """`n` cells, all with a centroid inside the vessel, with the given per-cell area fraction."""
    ctx = pd.DataFrame({
        "object_id": [f"c{i}" for i in range(n)],
        "in_islet": np.zeros(n, dtype=bool),
        f"in_{VESSEL}": np.ones(n, dtype=bool),
        "in_Nerve": np.zeros(n, dtype=bool),
        "in_Lymphatic": np.zeros(n, dtype=bool),
    })
    values = pd.DataFrame({CLF: np.full(n, clf_value), "Ker8_18": np.full(n, other)})
    area = {VESSEL: np.asarray(frac, dtype=float),
            "Nerve": np.full(n, np.nan), "Lymphatic": np.full(n, np.nan)}
    return ctx, values, area


def _call(frac, top, clf_value=1.0):
    n = len(frac)
    ctx, values, area = _frame(n, frac, clf_value=clf_value)
    levels = {VESSEL: 0.5, "Nerve": 0.5, "Lymphatic": 0.5}
    return SUB._structure_calls(ctx, values, levels, _log(),
                                np.array(top, dtype=object), True, area)


# ---------------------------------------------------------------------------------------------
# Step (b) -- the majority of the cell's area, which is the step that was missing
# ---------------------------------------------------------------------------------------------
def test_a_centroid_inside_the_object_is_not_by_itself_a_candidate():
    """The defect this rule exists to stop: a cell the vessel merely runs past."""
    cls, _, passed, audit = _call([0.10, 0.49], [CLF, CLF])
    assert cls.tolist() == [None, None], (
        "a cell with a centroid inside the object but a MINORITY of its area inside it is not a "
        "candidate; calling it treats 'the vessel runs through me' as 'I am the vessel'")
    assert audit[VESSEL]["inside"] == 2
    assert audit[VESSEL]["inside_and_majority"] == 0
    assert audit[VESSEL]["majority_withholds"] == 2


def test_the_majority_is_a_half_and_the_boundary_is_strict():
    cls, _, _, _ = _call([0.5, 0.5000001], [CLF, CLF])
    assert cls[0] is None, "exactly half is not a majority"
    assert cls[1] == VESSEL
    assert SUB.MAJORITY_AREA == 0.5, "the majority is what the word means; it is not tuned"


def test_a_majority_cell_whose_own_marker_leads_is_called():
    cls, resolved, _, audit = _call([0.90], [CLF])
    assert cls.tolist() == [VESSEL]
    assert resolved.tolist() == ["structure"]
    assert audit[VESSEL]["inside_and_majority"] == 1
    assert audit[VESSEL]["called"] == 1


# ---------------------------------------------------------------------------------------------
# Step (c) -- the marker that CREATED the object
# ---------------------------------------------------------------------------------------------
def test_a_majority_cell_whose_own_marker_does_not_lead_is_not_called():
    cls, _, passed, audit = _call([0.90], ["Ker8_18"])
    assert cls.tolist() == [None], (
        "an exocrine cell that happens to sit inside a vessel keeps its own identity")
    assert audit[VESSEL]["top_marker_withholds"] == 1


def test_only_the_markers_that_drew_the_object_can_carry_it():
    """Three objects, four markers. Nothing else may stand in for them."""
    assert SUB.CLASSIFIER_MARKER == {
        "Nerve": ("B3TUBB",),
        "Lymphatic": ("Podoplanin",),
        "Blood vessel": ("Caveolin", "CD34"),
    }
    assert set(SUB.STRUCTURES) == set(SUB.CLASSIFIER_MARKER)


# ---------------------------------------------------------------------------------------------
# Step (d) -- the passthrough is a MEASUREMENT
# ---------------------------------------------------------------------------------------------
def test_the_passthrough_carries_the_fraction_not_a_flag():
    _, _, passed, _ = _call([0.10, 0.30], [CLF, CLF])
    assert passed[VESSEL].dtype.kind == "f", "being crossed is a measurement, not a flag"
    np.testing.assert_allclose(passed[VESSEL], [0.10, 0.30])


def test_a_called_cell_is_not_also_a_passthrough_of_its_own_object():
    _, _, passed, _ = _call([0.90], [CLF])
    assert np.isnan(passed[VESSEL][0]), (
        "a cell the object CLASSIFIES is not a cell the object passes through")


# ---------------------------------------------------------------------------------------------
# Unavailable is a third outcome, and is never a zero
# ---------------------------------------------------------------------------------------------
def test_an_unmeasured_overlap_is_neither_a_candidate_nor_a_measured_minority():
    cls, _, passed, audit = _call([np.nan], [CLF])
    assert cls.tolist() == [None], "a cell we cannot measure cannot be shown to be majority-inside"
    assert audit[VESSEL]["inside_area_unavailable"] == 1
    assert audit[VESSEL]["majority_withholds"] == 0, (
        "an unmeasured cell must not be tallied as 'measured and not majority' -- that would "
        "report a missing measurement as a finding")
    assert np.isnan(passed[VESSEL][0]), "unavailable overlap stays NaN; it is not a zero overlap"


def test_a_measured_zero_and_an_unmeasured_cell_are_different_answers():
    _, _, passed, audit = _call([0.0, np.nan], [CLF, CLF])
    assert passed[VESSEL][0] == 0.0, "a measured zero means the object was tested and did not touch"
    assert np.isnan(passed[VESSEL][1])
    assert audit[VESSEL]["majority_withholds"] == 1
    assert audit[VESSEL]["inside_area_unavailable"] == 1


def test_the_overlap_is_not_restricted_to_islet_associated_cells():
    """The rule is asked of every cell; the producer's scope must not narrow it."""
    source = (REPO / "core/multiplex5000/scripts/typing/build_crossed_fraction.py").read_text()
    assert 'cells["cell_region_20um"].isin(["core", "expansion"])' not in source, (
        "correction D restricted itself to islet core and expansion, which left most of the "
        "tissue with no overlap measurement and made the majority test unanswerable there")
    assert "measured_cells" in source


# ---------------------------------------------------------------------------------------------
# The expectation itself, so nobody 'fixes' a small class
# ---------------------------------------------------------------------------------------------
def test_no_floor_is_asserted_on_how_many_cells_a_structure_calls():
    """Most of these objects are thin processes, so few cells are majority-inside. That is the
    rule working, and a test that demanded a minimum would be a test against the tissue."""
    cls, _, _, audit = _call([0.02] * 100 + [0.95], [CLF] * 101)
    assert audit[VESSEL]["called"] == 1 and audit[VESSEL]["inside"] == 101
    assert (cls == VESSEL).sum() == 1
