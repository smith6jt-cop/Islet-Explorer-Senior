"""Immune subtyping: what the lane refuses to say, and the arithmetic that holds it together.

The risk this build carries is not an arithmetic slip.  It is a number appearing where a
measurement does not exist -- a rule scored as a negative on a donor that never carried its
marker, a NaN quietly reading as "below the level", a count of a channel reported as a count of a
cell type, or a cohort fact written as a literal so nothing can see it drift.  So most of what is
asserted here is about DEFINITIONS and REFUSALS, against the modules that are their authority
(``cohort``, ``lineages``, ``rules``), never against a frozen copy of a result.

The tests that touch the 2.9M-row artifact skip when it is not on disk.

Run with
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
      core/multiplex5000/tests/current/test_immune_subtype.py -q
"""

from __future__ import annotations

import glob
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
SCRIPT = REPO / "core/multiplex5000/scripts/typing/immune_subtype.py"
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/typing"))
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/celltype"))
SPEC = importlib.util.spec_from_file_location("immune_subtype", SCRIPT)
sub = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sub)

import cohort as COHORT      # noqa: E402
import lineages as LIN       # noqa: E402

OUT = REPO / "data/development/typing/immune"


def frame(**columns) -> pd.DataFrame:
    return pd.DataFrame(columns)


# --------------------------------------------------------------------------------------------
# A missing measurement is unavailable.  Never zero, never negative, never False.
# --------------------------------------------------------------------------------------------

def test_a_rule_naming_a_marker_the_donor_lacks_is_unavailable_not_false():
    scaled = frame(CD3e=[0.9, 0.0])
    fires, undecided, missing = sub.evaluate(scaled, (["CD3e", "CD38"], [], []), {"CD3e"})
    assert fires is None and undecided is None
    assert missing == ["CD38"], "the missing marker has to be named, not merely refused"


def test_a_nan_on_one_cell_makes_that_cell_undecided_and_never_negative():
    scaled = frame(CD3e=[0.9, np.nan, 0.01], CD8=[0.9, 0.9, 0.9])
    fires, undecided, missing = sub.evaluate(scaled, (["CD3e", "CD8"], [], []), {"CD3e", "CD8"})
    assert missing == []
    assert list(fires) == [True, False, False]
    assert list(undecided) == [False, True, False], \
        "a NaN must be undecided; `value >= level` is False for NaN and would read as a negative"


def test_a_negative_constraint_on_a_nan_is_also_undecided():
    # The failure mode this is about: `none_of` on a marker the cell has no value for would score
    # True -- the cell is "not positive" -- and the rule would fire on a cell nothing was read on.
    scaled = frame(CD3e=[0.9, 0.9], CD8=[0.01, np.nan])
    fires, undecided, _ = sub.evaluate(scaled, (["CD3e"], [], ["CD8"]), {"CD3e", "CD8"})
    assert list(fires) == [True, False]
    assert list(undecided) == [False, True]


def test_a_cell_a_rule_named_is_not_reported_undecided_by_another_rule_it_could_not_read():
    reqs = {"named": (["CD3e"], [], []), "unreadable": (["CD4"], [], [])}
    scaled = frame(CD3e=[0.9], CD4=[np.nan])
    call, n, _, undecided, skipped = sub.resolve(
        scaled, ["named", "unreadable"], reqs, scaled, {"CD3e", "CD4"},
        {"named": 1, "unreadable": 1})
    assert list(call) == ["named"] and list(n) == [1] and skipped == {}
    assert list(undecided) == [False], "a cell that WAS named is decided"


def test_an_unavailable_rule_is_counted_EMPTY_and_never_zero():
    counts = pd.Series(dtype=int)
    value, status, why = sub._count_or_empty(counts, "Imm-plasma-cell", {"Imm-plasma-cell": ["CD38"]},
                                             True, "Immune", {})
    assert value == "" and status == "EMPTY_rule_unavailable"
    assert "CD38" in why and "not a negative answer" in why


def test_a_subtype_whose_parent_was_never_named_is_EMPTY_not_zero():
    value, status, why = sub._count_or_empty(pd.Series(dtype=int), "T-CD8", {}, False, "Imm-T-cell", {})
    assert value == "" and status == "EMPTY_parent_absent" and "Imm-T-cell" in why


def test_a_rule_that_could_be_read_and_did_not_fire_is_a_measured_zero():
    value, status, why = sub._count_or_empty(pd.Series({"Imm-B-cell": 4}), "Imm-NK-cell",
                                             {}, True, "Immune", {})
    assert value == 0 and status == "measured" and why == ""


# --------------------------------------------------------------------------------------------
# Acquisition is measured from the donor's own partition, not assumed from a schema.
# --------------------------------------------------------------------------------------------

def test_a_column_that_is_entirely_null_is_not_an_acquired_marker(tmp_path, monkeypatch):
    part = tmp_path / "donor_id=T1"
    part.mkdir()
    pq.write_table(pa.table({"CD3e": pa.array([1.0, 2.0]),
                             "CD38": pa.array([None, None], type=pa.float64())}),
                   part / "data.parquet")
    monkeypatch.setattr(sub.S, "ALIGNED", tmp_path)
    got, absent, empty = sub.acquired_markers("T1", {"CD3e", "CD38", "CD57"})
    assert got == ["CD3e"]
    assert absent == ["CD57"], "a column that is not there at all is absent"
    assert empty == ["CD38"], "a column that is there and entirely null was never acquired"


def test_a_donor_with_no_partition_is_credited_with_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(sub.S, "ALIGNED", tmp_path)
    got, absent, empty = sub.acquired_markers("nobody", {"CD3e", "CD38"})
    assert got == [] and absent == ["CD38", "CD3e"] and empty == []


# --------------------------------------------------------------------------------------------
# A SUPPORTING_ONLY marker may not carry a call alone, so the label has to say what it names.
# --------------------------------------------------------------------------------------------

def test_foxp3_and_cd56_are_the_supporting_only_markers_this_lane_meets():
    assert {"FOXP3", "CD56"} <= LIN.SUPPORTING_ONLY


def test_a_rule_requiring_foxp3_says_what_it_actually_names():
    note = sub.supporting_only_note("T-regulatory",
                                    (["CD3e", "CD4", "FOXP3"], [["CD39", "ICOS"]], []))
    assert "FOXP3" in note and "SUPPORTING_ONLY" in note
    assert "reading highest on FOXP3" in note
    assert "not as a count of T-regulatory by an independent marker" in note
    assert "CD3e" in note and "CD4" in note and "CD39" in note


def test_a_supporting_only_marker_used_negatively_gets_its_own_sentence():
    note = sub.supporting_only_note("T-CD4-helper", (["CD3e", "CD4"], [], ["CD8", "FOXP3"]))
    assert "excludes cells above the level on FOXP3" in note


def test_a_supporting_only_marker_with_alternatives_in_its_or_group_is_not_load_bearing():
    # CD11b sits in a nine-member OR group; any other member satisfies it, so the rule does not
    # rest on it and claiming otherwise would put a caveat on almost every label.
    assert sub.supporting_only_note("Imm-macrophage",
                                    ([], [["CD68", "Iba1"], ["CD11b", "CD68", "MPO"]], [])) == ""


# --------------------------------------------------------------------------------------------
# Terms: what separates the cells, and what only looks like it does.
# --------------------------------------------------------------------------------------------

def test_a_term_the_parent_already_required_is_not_reported_as_a_majority_term():
    fractions = {"CD68": 0.98, "HLA_DR": 0.9}
    parent = ([], [["CD68", "Iba1"]], [])
    got = sub.majority_terms((["HLA_DR"], [["CD68", "Iba1"]], []), fractions, parent)
    assert [t["marker"] for t in got] == ["HLA_DR"], \
        "CD68 is high in macrophages because that is what named them -- reporting it reports the selection"


def test_restrict_drops_only_or_group_members_and_never_empties_a_group():
    kept = sub.restrict((["CD56"], [["CD57", "Granzyme_B"], ["CD11b"]], ["CD3e"]),
                        {"CD57", "CD11b", "CD56", "CD3e"})
    all_of, groups, none_of = kept
    assert all_of == ["CD56"], "a required term is the rule's meaning and is never dropped"
    assert none_of == ["CD3e"], "a negative term is never dropped either"
    assert groups[0] == ["Granzyme_B"]
    assert groups[1] == ["CD11b"], "a group cannot be emptied -- it keeps its only member"


def test_eta_squared_is_one_for_a_perfect_split_and_the_null_is_reported():
    values = np.array([0.0, 0.0, 1.0, 1.0])
    groups = np.array(["a", "a", "b", "b"])
    e2, p, null = sub.eta_squared_with_null(values, groups, draws=200, seed=1)
    assert e2 == pytest.approx(1.0)
    assert 0.0 <= p <= 1.0 and 0.0 <= null <= 1.0
    assert null < e2, "a shuffled label must not reach a perfect split as often as the truth does"


def test_the_permutation_null_is_what_makes_the_effect_size_mean_anything():
    # Ten groups of one: every arrangement separates perfectly, so eta-squared is 1 and carries
    # no information at all.  Only the null says so.
    values = np.arange(10.0)
    groups = np.array([f"run{i}" for i in range(10)])
    e2, p, null = sub.eta_squared_with_null(values, groups, draws=200, seed=1)
    assert e2 == pytest.approx(1.0) and null == pytest.approx(1.0) and p == pytest.approx(1.0)


# --------------------------------------------------------------------------------------------
# No cohort fact is written as a literal.
# --------------------------------------------------------------------------------------------

def test_no_donor_id_is_written_into_this_module_as_a_string():
    import ast
    tree = ast.parse(SCRIPT.read_text())
    docstrings = {id(n.body[0].value) for n in ast.walk(tree)
                  if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef))
                  and getattr(n, "body", None) and isinstance(n.body[0], ast.Expr)
                  and isinstance(n.body[0].value, ast.Constant)
                  and isinstance(n.body[0].value.value, str)}
    known = set(COHORT.CLEAN) | set(COHORT.MARGINAL)
    literals = {n.value for n in ast.walk(tree)
                if isinstance(n, ast.Constant) and isinstance(n.value, str)
                and id(n) not in docstrings}
    assert not (literals & known), \
        "a donor id in the code is a cohort fact the test suite cannot see change"


def test_the_donor_list_comes_from_the_cohort_module():
    assert COHORT.resolve([], include_marginal=True) == list(COHORT.CLEAN) + sorted(COHORT.MARGINAL)


# --------------------------------------------------------------------------------------------
# The artifact, when it is on disk.
# --------------------------------------------------------------------------------------------

def _coverage():
    path = OUT / "coverage.json"
    if not path.exists():
        pytest.skip("immune coverage.json not built")
    return json.loads(path.read_text())


@pytest.mark.data
def test_every_requested_donor_was_built_or_carries_a_reason():
    cov = _coverage()
    for donor in cov["donors_requested"]:
        record = cov["per_donor"].get(donor)
        assert record is not None, f"{donor} was requested and is not in the coverage report"
        if record["state"] != "built":
            assert record.get("reason"), f"{donor} was not built and no reason was written"


@pytest.mark.data
def test_a_donor_left_out_of_the_run_is_named_with_the_cohort_module_s_reason():
    cov = _coverage()
    for donor, why in cov["donors_not_requested"].items():
        assert why == COHORT.MARGINAL[donor]


@pytest.mark.data
def test_no_count_is_written_as_zero_where_the_rule_could_not_be_evaluated():
    path = OUT / "per_donor_counts.csv"
    if not path.exists():
        pytest.skip("per_donor_counts.csv not built")
    table = pd.read_csv(path, dtype=str).fillna("")
    empty = table[table["status"].str.startswith("EMPTY")]
    assert (empty["count"] == "").all(), "an unevaluable rule must be EMPTY, never 0"
    assert (empty["reason"].str.len() > 0).all(), "every EMPTY carries its reason"
    measured = table[table["status"].str.startswith("measured")]
    assert measured["count"].str.match(r"^\d+$").all()


@pytest.mark.data
def test_a_count_that_is_not_cross_donor_comparable_says_so_with_its_numbers():
    path = OUT / "per_donor_counts.csv"
    if not path.exists():
        pytest.skip("per_donor_counts.csv not built")
    table = pd.read_csv(path, dtype=str).fillna("")
    flagged = table[table["status"] == "measured_not_cross_donor_comparable"]
    for _, row in flagged.iterrows():
        assert "eta2 by acq_run" in row["reason"]
        assert "permutation p" in row["reason"]
        assert "NOT comparable" in row["reason"]


@pytest.mark.data
def test_the_written_cells_are_exactly_the_Immune_cells_and_the_subsets_sum_to_their_parent():
    parts = sorted(glob.glob(str(OUT / "donor_id=*/data.parquet")))
    if not parts:
        pytest.skip("no immune partitions built")
    cov = _coverage()
    for part in parts:
        donor = Path(part).parent.name.split("=", 1)[1]
        got = pd.read_parquet(part, columns=["immune_class", "immune_class_undecided",
                                             "immune_subtype"])
        record = cov["per_donor"][donor]
        assert len(got) == record["n_immune"]
        named = int(got["immune_class"].notna().sum())
        undecided = int(got["immune_class_undecided"].sum())
        assert named + undecided + record["uncalled_no_rule_fired"] == len(got)
        for parent in sub.LEVEL3_PARENTS:
            inside = got["immune_class"] == parent
            named_here = got.loc[inside, "immune_subtype"].notna().sum()
            assert named_here <= int(inside.sum()), "a subset cannot exceed its parent"
        outside = got["immune_subtype"].notna() & ~got["immune_class"].isin(sub.LEVEL3_PARENTS)
        assert not outside.any(), "a subtype may only be written under its own parent"
