"""The pLN v3 ``pct_immune`` denominator is the cells WITH a nucleus.

``per_donor_summary.csv`` records ``no_nucleus_excluded`` on the same row and then, until
this test, divided ``immune`` by every cell of the frame including those.  The
cross-modality stage carried that ratio into ``xen_pct_immune`` and ``delta_pct_immune``.
Both sites now go through one pure function, asserted here on two-row synthetic frames.
No number below is biological; the arithmetic is the claim.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multiplex5000.partial_run_guard import RunScope
from multiplex5000.tissues import pln_clustering

XENIUM = Path(__file__).resolve().parents[1]
SCRIPT = XENIUM / "scripts" / "pln_immune_typing_v3.py"


@pytest.fixture(scope="module")
def v3():
    spec = importlib.util.spec_from_file_location("pln_immune_typing_v3", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_worker_uses_package_scope_and_partition(v3):
    assert v3.RunScope is RunScope
    assert v3.cluster_module() is pln_clustering


def _counts(frame: pd.DataFrame) -> tuple[int, int, int]:
    """Exactly what ``run_calls`` reads off a donor frame before it writes the summary row."""
    s1 = frame["stage1"].value_counts()
    return int(s1.get("Immune", 0)), len(frame), int(s1.get("Unavailable_no_nucleus", 0))


def test_a_no_nucleus_cell_leaves_the_pct_immune_denominator(v3):
    # Two rows: one immune cell with a nucleus, one cell with no nucleus.
    frame = pd.DataFrame({"stage1": ["Immune", "Unavailable_no_nucleus"]})
    immune, cells, no_nucleus = _counts(frame)
    assert (immune, cells, no_nucleus) == (1, 2, 1)
    assert v3.pct_immune_of_nucleated(immune, cells, no_nucleus) == 100.0, (
        "the one cell with a nucleus is immune; the cell with no nucleus is not a cell we "
        "can be sure of and is excluded from every calculation")
    # The line this replaces divided by every row of the frame.
    assert round(100.0 * immune / max(cells, 1), 2) == 50.0


def test_a_frame_with_every_nucleus_present_is_numerically_unchanged(v3):
    frame = pd.DataFrame({"stage1": ["Immune", "Non_immune"]})
    immune, cells, no_nucleus = _counts(frame)
    assert no_nucleus == 0
    assert v3.pct_immune_of_nucleated(immune, cells, no_nucleus) == 50.0
    assert v3.pct_immune_of_nucleated(immune, cells, no_nucleus) == round(
        100.0 * immune / max(cells, 1), 2), "a row that had every measurement does not move"


def test_no_nucleated_cell_is_unavailable_not_a_measured_zero(v3):
    frame = pd.DataFrame({"stage1": ["Unavailable_no_nucleus", "Unavailable_no_nucleus"]})
    immune, cells, no_nucleus = _counts(frame)
    value = v3.pct_immune_of_nucleated(immune, cells, no_nucleus)
    assert math.isnan(value), "0 immune over 0 nucleated cells is unavailable, never 0.0"


def test_the_cross_modality_delta_uses_the_nucleated_denominator(v3):
    # Two summary rows as run_cross reads them back: one donor with a no-nucleus cell,
    # one without. The PhenoCycler side is a synthetic percentage.
    xen = pd.DataFrame([
        {"donor_id": "A", "cells": 4, "immune": 2, "no_nucleus_excluded": 2},
        {"donor_id": "B", "cells": 4, "immune": 2, "no_nucleus_excluded": 0},
    ])
    pheno_pct = {"A": 40.0, "B": 40.0}
    got = {}
    for r in xen.itertuples():
        xen_pct = v3.pct_immune_of_nucleated(r.immune, r.cells, r.no_nucleus_excluded)
        got[r.donor_id] = (xen_pct, xen_pct - pheno_pct[r.donor_id])
    assert got["A"] == (100.0, 60.0), "2 immune of 2 nucleated cells, not of 4"
    assert got["B"] == (50.0, 10.0), "no no-nucleus cell: unchanged from the old arithmetic"


# ---------------------------------------------------------------------------
# The nucleus rule is MEASURED on a donor only when every cell of its frame carries
# ``nucleus_measured`` True.  ``call_section`` forces ``has_nucleus`` all-True when no
# vendor bundle covers a section, so on such a frame no cell can carry
# ``Unavailable_no_nucleus`` and a 0 excluded is the rule never having run -- not a 0.
# ---------------------------------------------------------------------------
def _frame(nucleus_measured, stage1=("Immune", "Immune", "Non_immune", "Non_immune")):
    return pd.DataFrame({"stage1": list(stage1),
                         "nucleus_measured": [nucleus_measured] * len(stage1)})


def test_an_unmeasured_nucleus_rule_makes_the_exclusion_and_the_share_unavailable(v3):
    row = v3.immune_share_of_donor(_frame(False))
    assert row["cells"] == 4 and row["immune"] == 2, "what WAS called is still a count"
    assert math.isnan(row["no_nucleus_excluded"]), (
        "0 excluded would be a measurement; on this frame the rule never ran")
    assert math.isnan(row["pct_immune"]), "a share of nucleated cells needs the nucleus measured"
    assert row["nucleus_not_measured"] == 4
    assert row["nucleus_rule_state"] == v3.NUCLEUS_RULE_UNAVAILABLE
    reason = row["nucleus_rule_reason"]
    assert "all 4 cells" in reason and "never measured" in reason
    assert "not a measured absence" in reason


def test_the_same_frame_with_the_nucleus_measured_gives_the_nucleated_denominator(v3):
    row = v3.immune_share_of_donor(_frame(True))
    assert row["cells"] == 4 and row["immune"] == 2
    assert row["no_nucleus_excluded"] == 0 and row["pct_immune"] == 50.0
    assert row["nucleus_not_measured"] == 0
    assert row["nucleus_rule_state"] == v3.NUCLEUS_RULE_MEASURED
    assert row["nucleus_rule_reason"] == ""
    # One no-nucleus cell among the four leaves both numerator and denominator.
    row = v3.immune_share_of_donor(
        _frame(True, ("Immune", "Immune", "Non_immune", "Unavailable_no_nucleus")))
    assert row["no_nucleus_excluded"] == 1
    assert row["pct_immune"] == round(100.0 * 2 / 3, 2)
    assert row["nucleus_rule_state"] == v3.NUCLEUS_RULE_MEASURED


def test_a_partly_measured_donor_is_unavailable_with_both_counts_in_the_reason(v3):
    frame = pd.DataFrame({
        "stage1": ["Immune", "Immune", "Non_immune", "Unavailable_no_nucleus"],
        "nucleus_measured": [True, True, False, True],
    })
    row = v3.immune_share_of_donor(frame)
    assert math.isnan(row["no_nucleus_excluded"]) and math.isnan(row["pct_immune"])
    assert row["nucleus_not_measured"] == 1
    assert row["nucleus_rule_state"] == v3.NUCLEUS_RULE_UNAVAILABLE
    reason = row["nucleus_rule_reason"]
    assert "1 of 4 cells" in reason, "how many cells the rule never saw"
    assert "excluded 1" in reason, "how many the rule did exclude where it ran"


def test_a_cell_found_without_a_nucleus_that_was_never_measured_is_a_contradiction(v3):
    frame = pd.DataFrame({"stage1": ["Immune", "Unavailable_no_nucleus"],
                          "nucleus_measured": [True, False]})
    with pytest.raises(ValueError, match="Unavailable_no_nucleus"):
        v3.immune_share_of_donor(frame)


def test_run_cross_recompute_carries_an_unavailable_exclusion_as_nan(v3, tmp_path):
    # The summary round trip run_cross makes: two donors, one measured, one not.
    rows = [
        {"donor_id": "A", **v3.immune_share_of_donor(_frame(True))},
        {"donor_id": "B", **v3.immune_share_of_donor(_frame(False))},
    ]
    summary = v3.summary_frame(rows)
    path = tmp_path / "per_donor_summary.csv"
    summary.to_csv(path, index=False)
    text = path.read_text()
    assert "A,4,2,0,50.0" in text, "a measured donor keeps integer counts on disk"
    back = pd.read_csv(path, dtype={"donor_id": str})
    a, b = (back.set_index("donor_id").loc[d] for d in ("A", "B"))
    assert v3.pct_immune_of_nucleated(a.immune, a.cells, a.no_nucleus_excluded) == 50.0
    assert math.isnan(b.no_nucleus_excluded)
    assert math.isnan(v3.pct_immune_of_nucleated(b.immune, b.cells, b.no_nucleus_excluded)), (
        "an unavailable exclusion count makes the recomputed share unavailable")
    assert b.nucleus_rule_state == v3.NUCLEUS_RULE_UNAVAILABLE


def test_a_summary_without_the_state_column_derives_it_from_the_calls(v3, tmp_path):
    # A per_donor_summary.csv written before the column existed: no_nucleus_excluded 0 on
    # a donor whose nucleus was never measured. The calls partition says which it is.
    calls_root = tmp_path / "calls"
    for donor, measured in (("A", True), ("B", False)):
        dest = calls_root / f"donor_id={donor}/roi={v3.ROI}"
        dest.mkdir(parents=True)
        _frame(measured).to_parquet(dest / "data.parquet", index=False)
    old = pd.DataFrame([
        {"donor_id": "A", "cells": 4, "immune": 2, "no_nucleus_excluded": 0, "pct_immune": 50.0},
        {"donor_id": "B", "cells": 4, "immune": 2, "no_nucleus_excluded": 0, "pct_immune": 50.0},
    ])
    fixed = v3.with_nucleus_rule_state(old, calls_root).set_index("donor_id")
    assert fixed.at["A", "nucleus_rule_state"] == v3.NUCLEUS_RULE_MEASURED
    assert fixed.at["A", "no_nucleus_excluded"] == 0 and fixed.at["A", "pct_immune"] == 50.0
    assert fixed.at["B", "nucleus_rule_state"] == v3.NUCLEUS_RULE_UNAVAILABLE
    assert math.isnan(fixed.at["B", "no_nucleus_excluded"])
    assert math.isnan(fixed.at["B", "pct_immune"])
    assert "never measured" in fixed.at["B", "nucleus_rule_reason"]
    # A summary that already carries the column is read as written.
    same = v3.with_nucleus_rule_state(fixed.reset_index(), tmp_path / "nowhere")
    assert same.set_index("donor_id").at["B", "nucleus_rule_state"] == v3.NUCLEUS_RULE_UNAVAILABLE
    # Without the column AND without the calls partition, nothing can say which it is.
    with pytest.raises(SystemExit, match="re-run"):
        v3.with_nucleus_rule_state(old, tmp_path / "nowhere")


# ---------------------------------------------------------------------------
# ``nucleus_measured`` is read at ONE boundary and refused there unless boolean.
# ``astype(bool)`` reads the string 'False' as True and an int8 -1 as True, so until
# 2026-09-22 a CSV round trip or a wrong dtype reported the donor MEASURED with 0
# excluded -- the exact silent 0 the state column exists to prevent.
# ---------------------------------------------------------------------------
def test_nucleus_measured_as_strings_is_refused_at_the_boundary(v3):
    flags = pd.Series(["False"] * 4)
    # object under pandas 2, the str (StringDtype) default under pandas 3; both have kind "O"
    # and both reach the same refusal in nucleus_flags_of.
    assert flags.dtype.kind == "O"
    with pytest.raises(ValueError, match="object") as err:
        v3.nucleus_rule_of(flags)
    assert "nucleus_measured" in str(err.value) and "'False'" in str(err.value), (
        "the refusal names the column and the value that is not a bool")
    frame = pd.DataFrame({"stage1": ["Immune"] * 4, "nucleus_measured": ["False"] * 4})
    with pytest.raises(ValueError, match="object"):
        v3.immune_share_of_donor(frame)


def test_nucleus_measured_as_int8_is_refused_at_the_boundary(v3):
    flags = pd.Series([-1] * 4, dtype="int8")
    with pytest.raises(ValueError, match="int8"):
        v3.nucleus_rule_of(flags)
    frame = pd.DataFrame({"stage1": ["Immune"] * 4,
                          "nucleus_measured": pd.Series([-1] * 4, dtype="int8")})
    with pytest.raises(ValueError, match="int8"):
        v3.immune_share_of_donor(frame)


@pytest.mark.parametrize("dtype", ["bool", "boolean", "object"])
def test_a_boolean_nucleus_measured_is_read_as_before(v3, dtype):
    # numpy bool, pandas nullable boolean, and object holding only bools are all
    # boolean columns; every one reads as it did before the boundary went in.
    assert v3.nucleus_rule_of(pd.Series([True] * 4, dtype=dtype)) == (
        v3.NUCLEUS_RULE_MEASURED, "", 0)
    state, reason, n_not = v3.nucleus_rule_of(pd.Series([True, False, True, True], dtype=dtype))
    assert (state, n_not) == (v3.NUCLEUS_RULE_UNAVAILABLE, 1) and "1 of 4 cells" in reason
    row = v3.immune_share_of_donor(pd.DataFrame({
        "stage1": ["Immune", "Immune", "Non_immune", "Unavailable_no_nucleus"],
        "nucleus_measured": pd.Series([True] * 4, dtype=dtype)}))
    assert row["nucleus_rule_state"] == v3.NUCLEUS_RULE_MEASURED
    assert row["no_nucleus_excluded"] == 1 and row["pct_immune"] == round(100.0 * 2 / 3, 2)


def test_a_null_nucleus_measured_is_admitted_and_counted_as_not_measured(v3):
    # A null in a boolean column is a rule that did not run on that cell -- the
    # donor is unavailable, and the reason says the value was null, not False.
    for flags in (pd.Series([True, None, True, True], dtype="boolean"),
                  pd.Series([True, None, True, True], dtype=object)):
        state, reason, n_not = v3.nucleus_rule_of(flags)
        assert (state, n_not) == (v3.NUCLEUS_RULE_UNAVAILABLE, 1)
        assert "1 of 4 cells" in reason and "null" in reason


# ---------------------------------------------------------------------------
# Each count argument of ``pct_immune_of_nucleated`` is read at the boundary: a
# non-negative integer, or null.  Until 2026-09-22 ``int()`` truncated 2.5 to 2 and
# the arithmetic ran on a negative or an impossible count: (-1, 4, 0) gave -25.0,
# (5, 4, 0) gave 125.0 and (2, 4, -1) gave 40.0, every one of them a number on a row.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("args, named", [
    ((-1, 4, 0), "immune"),            # a negative count
    ((5, 4, 0), "immune"),             # more immune cells than nucleated cells
    ((2, 4, -1), "no_nucleus_excluded"),  # a negative exclusion count
    ((2, 4, 5), "no_nucleus_excluded"),   # more excluded than there are cells
    ((2.5, 4, 0), "immune"),           # not an integer: 2.5 is not a count of cells
    ((2, 4.5, 0), "cells"),            # not an integer
    ((2, float("inf"), 0), "cells"),   # not finite
    (("two", 4, 0), "immune"),         # not a number
    ((True, 4, 0), "immune"),          # a bool is not a count
])
def test_an_input_that_cannot_be_a_count_gives_nan_with_a_reason(v3, args, named):
    assert math.isnan(v3.pct_immune_of_nucleated(*args))
    value, reason = v3.pct_immune_of_nucleated_with_reason(*args)
    assert math.isnan(value)
    assert named in reason, f"the reason names the refused argument: {reason!r}"


def test_integral_floats_and_numeric_strings_are_counts(v3):
    # A CSV column with one blank reads as float, so 4.0 is the count 4; a numeric
    # string is what a dtype=str read returns. Both are admitted and read as before.
    assert v3.pct_immune_of_nucleated("2", "4", "0") == 50.0
    assert v3.pct_immune_of_nucleated(2.0, 4.0, 0.0) == 50.0
    assert v3.pct_immune_of_nucleated(np.int64(2), np.float64(4.0), np.int32(0)) == 50.0
    assert v3.pct_immune_of_nucleated_with_reason(2, 4, 0) == (50.0, "")
    # A null argument is unavailable, not refused: NaN with no reason, as before.
    for args in ((None, 4, 0), (2, float("nan"), 0), (2, 4, pd.NA)):
        assert v3.pct_immune_of_nucleated_with_reason(*args) == (
            pytest.approx(float("nan"), nan_ok=True), "")
    # The share of nucleated cells is bounded: immune <= nucleated is the last check.
    assert v3.pct_immune_of_nucleated(4, 4, 0) == 100.0
    assert v3.pct_immune_of_nucleated(0, 4, 0) == 0.0


def test_each_refused_argument_carries_its_own_reason(v3):
    value, reason = v3.pct_immune_of_nucleated_with_reason(-1, 4.5, 0)
    assert math.isnan(value)
    assert "immune" in reason and "cells" in reason, "one refusal per column"


# ---------------------------------------------------------------------------
# ``nucleus_rule_state`` on a summary row carries one of its two declared values or
# is derived.  Until 2026-09-22 a summary that HAD the column was read as written,
# so a null, an empty or a 'bogus' state passed through with pct_immune 50.0.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("written", [float("nan"), "", "bogus"])
def test_a_null_or_out_of_vocabulary_state_is_derived_or_unavailable(v3, tmp_path, written):
    # Three donors with the column present and unreadable: A's partition says
    # measured, B's says never measured, C has no partition. D is read as written.
    calls_root = tmp_path / "calls"
    for donor, measured in (("A", True), ("B", False)):
        dest = calls_root / f"donor_id={donor}/roi={v3.ROI}"
        dest.mkdir(parents=True)
        _frame(measured).to_parquet(dest / "data.parquet", index=False)
    xen = pd.DataFrame([
        {"donor_id": d, "cells": 4, "immune": 2, "no_nucleus_excluded": 0,
         "pct_immune": 50.0, "nucleus_rule_state": written, "nucleus_rule_reason": ""}
        for d in ("A", "B", "C")
    ] + [{"donor_id": "D", "cells": 4, "immune": 2, "no_nucleus_excluded": 0,
          "pct_immune": 50.0, "nucleus_rule_state": v3.NUCLEUS_RULE_MEASURED,
          "nucleus_rule_reason": ""}])
    fixed = v3.with_nucleus_rule_state(xen, calls_root).set_index("donor_id")
    assert set(fixed["nucleus_rule_state"]) <= {v3.NUCLEUS_RULE_MEASURED,
                                                v3.NUCLEUS_RULE_UNAVAILABLE}
    # A: the partition exists and every cell was measured -> measured, counts kept.
    assert fixed.at["A", "nucleus_rule_state"] == v3.NUCLEUS_RULE_MEASURED
    assert fixed.at["A", "no_nucleus_excluded"] == 0 and fixed.at["A", "pct_immune"] == 50.0
    assert repr(written) in fixed.at["A", "nucleus_rule_reason"], (
        "the row says what state it carried before it was derived")
    # B: the partition exists and says the rule never ran -> unavailable, both NaN.
    assert fixed.at["B", "nucleus_rule_state"] == v3.NUCLEUS_RULE_UNAVAILABLE
    assert math.isnan(fixed.at["B", "no_nucleus_excluded"])
    assert math.isnan(fixed.at["B", "pct_immune"])
    assert "never measured" in fixed.at["B", "nucleus_rule_reason"]
    # C: no partition to derive from -> the ROW is unavailable (no SystemExit: the
    # summary has the column; it is this one row that cannot be read).
    assert fixed.at["C", "nucleus_rule_state"] == v3.NUCLEUS_RULE_UNAVAILABLE
    assert math.isnan(fixed.at["C", "no_nucleus_excluded"])
    assert math.isnan(fixed.at["C", "pct_immune"])
    reason = fixed.at["C", "nucleus_rule_reason"]
    assert repr(written) in reason and "donor_id=C" in reason
    # D: a readable state is read as written.
    assert fixed.at["D", "nucleus_rule_state"] == v3.NUCLEUS_RULE_MEASURED
    assert fixed.at["D", "pct_immune"] == 50.0 and fixed.at["D", "nucleus_rule_reason"] == ""


def test_a_summary_without_the_column_still_refuses_when_no_partition_exists(v3, tmp_path):
    # The contrast with the row-level case above: WITHOUT the column at all and
    # without the partition, the whole summary is refused (SystemExit, as before).
    old = pd.DataFrame([
        {"donor_id": "C", "cells": 4, "immune": 2, "no_nucleus_excluded": 0, "pct_immune": 50.0},
    ])
    with pytest.raises(SystemExit, match="re-run"):
        v3.with_nucleus_rule_state(old, tmp_path / "nowhere")


def test_run_cross_masks_the_xenium_share_wherever_the_state_is_not_measured(v3):
    # The Xenium side of one cross-modality row, from one per_donor_summary row as
    # run_cross reads it after with_nucleus_rule_state. Whatever the counts say, a
    # state other than 'measured' has no share of nucleated cells.
    base = {"cells": 4, "immune": 2, "no_nucleus_excluded": 0, "T_NK_compartment": 1,
            "T_CD4": 1, "T_CD8": 0, "T_Treg": 0, "T_CD4_CD8_both": 0,
            "T_subset_unresolved": 0, "NK_ILC": 0, "T_NK_unresolved": 0}
    rows = pd.DataFrame([
        {"donor_id": "A", **base, "nucleus_rule_state": v3.NUCLEUS_RULE_MEASURED,
         "nucleus_rule_reason": ""},
        {"donor_id": "B", **base, "nucleus_rule_state": v3.NUCLEUS_RULE_UNAVAILABLE,
         "nucleus_rule_reason": "the rule never ran"},
    ])
    got = {d["donor_id"]: d for d in (v3.xen_cross_row(r) for r in rows.itertuples())}
    assert got["A"]["xen_pct_immune"] == 50.0
    assert got["A"]["xen_nucleus_rule_state"] == v3.NUCLEUS_RULE_MEASURED
    assert math.isnan(got["B"]["xen_pct_immune"]), (
        "counts that would divide to 50.0 do not make a share under an unavailable state")
    assert got["B"]["xen_nucleus_rule_state"] == v3.NUCLEUS_RULE_UNAVAILABLE
    assert got["B"]["xen_nucleus_rule_reason"] == "the rule never ran"
    assert got["B"]["xen_immune"] == 2 and got["B"]["xen_cells"] == 4, "the counts stay counts"


def test_run_cross_writes_a_refused_count_reason_on_the_row(v3):
    # A measured state with immune 5 of 4 cells: the share is NaN and the reason
    # column of the row says which argument was refused.
    row = pd.DataFrame([{
        "donor_id": "A", "cells": 4, "immune": 5, "no_nucleus_excluded": 0,
        "T_NK_compartment": 1, "T_CD4": 1, "T_CD8": 0, "T_Treg": 0, "T_CD4_CD8_both": 0,
        "T_subset_unresolved": 0, "NK_ILC": 0, "T_NK_unresolved": 0,
        "nucleus_rule_state": v3.NUCLEUS_RULE_MEASURED, "nucleus_rule_reason": ""}])
    d = v3.xen_cross_row(next(row.itertuples()))
    assert math.isnan(d["xen_pct_immune"])
    assert "immune" in d["xen_nucleus_rule_reason"] and "5" in d["xen_nucleus_rule_reason"]
