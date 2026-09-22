"""Silent-failure contracts for the grid builder's cross-modality T_NK table.

A T_NK lineage a lane did not resolve is a composition row with ``n_cells`` NULL, a
state and a reason.  ``crossmodal_tnk`` must carry that NULL through to
``crossmodal_tnk.csv`` -- never a measured zero that a rank can order among the
donors both modalities measured.  Everything asserted here is arithmetic on a
synthetic frame; no number is biological and no result is a target.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
BUILDER = REPO / "scripts" / "senior" / "build_grid_composition.py"


def _load():
    spec = importlib.util.spec_from_file_location("build_grid_composition", BUILDER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def grid():
    return _load()


def _lineage(grid, tissue, modality, donor, lineage, n, *, state=None, reason="",
             is_immune="yes"):
    """One composition row as ``build`` writes it: a measured count, or NULL with a state."""
    measured = state is None or state == grid.STATE_MEASURED
    return {
        "tissue": tissue, "modality": modality, "donor_id": donor,
        "lineage_common": lineage, "bucket_kind": grid.BK_LINEAGE,
        "n_cells": int(n) if measured else pd.NA,
        "state": grid.STATE_MEASURED if measured else state,
        "reason": reason, "is_immune": is_immune,
    }


def _comp(rows):
    comp = pd.DataFrame(rows)
    comp["n_cells"] = comp["n_cells"].astype("Int64")
    return comp


def _by_donor(wide: pd.DataFrame, tissue: str) -> pd.DataFrame:
    return wide[wide["tissue"].eq(tissue)].set_index("donor_id")


# ---------------------------------------------------------------------------
# Two rows: one donor with T_NK measured, one with T_NK unavailable.
# ---------------------------------------------------------------------------
def test_an_unavailable_tnk_stays_null_with_its_state_and_reason(grid):
    reason = "this lane emitted no T_NK call although the crosswalk records pheno_resolution=yes"
    comp = _comp([
        _lineage(grid, "pancreas", "protein", "A", "T_NK", 10),
        _lineage(grid, "pancreas", "protein", "A", "B_Plasma", 30),
        _lineage(grid, "pancreas", "protein", "B", "T_NK", 0,
                 state=grid.STATE_UNAVAILABLE, reason=reason),
        _lineage(grid, "pancreas", "protein", "B", "B_Plasma", 50),
        # a non-immune lineage never enters the immune denominator
        _lineage(grid, "pancreas", "protein", "B", "Endocrine", 999, is_immune="no"),
    ])
    wide = _by_donor(grid.crossmodal_tnk(comp), "pancreas")

    # The donor that HAD a measurement is unchanged arithmetic.
    assert wide.at["A", "n_resolved_immune_protein"] == 40
    assert wide.at["A", "frac_T_NK_protein"] == pytest.approx(10 / 40)
    assert wide.at["A", "n_T_NK_protein"] == 10
    assert wide.at["A", "state_T_NK_protein"] == grid.STATE_MEASURED

    # The donor whose T_NK the lane did not resolve is NULL, not a measured zero.
    assert wide.at["B", "n_resolved_immune_protein"] == 50, (
        "the immune cells the lane DID resolve are still a count")
    assert pd.isna(wide.at["B", "n_T_NK_protein"]), "an unavailable T_NK count is NULL, not 0"
    assert pd.isna(wide.at["B", "frac_T_NK_protein"]), "NULL / 50 is NULL, not 0.0"
    assert wide.at["B", "state_T_NK_protein"] == grid.STATE_UNAVAILABLE
    assert wide.at["B", "reason_T_NK_protein"] == reason
    # RNA was never in this frame: absent on both donors, and said so.
    assert pd.isna(wide.at["A", "frac_T_NK_RNA"]) and pd.isna(wide.at["B", "frac_T_NK_RNA"])
    assert not wide["both_modalities"].any()


# ---------------------------------------------------------------------------
# A NULL must not enter the rank, and the ranked donors must not move because of it.
# ---------------------------------------------------------------------------
def _both_modalities_frame(grid, with_c: bool):
    rows = [
        _lineage(grid, "pancreas", "protein", "A", "T_NK", 10),
        _lineage(grid, "pancreas", "protein", "A", "B_Plasma", 30),
        _lineage(grid, "pancreas", "RNA", "A", "T_NK", 20),
        _lineage(grid, "pancreas", "RNA", "A", "B_Plasma", 20),
        _lineage(grid, "pancreas", "protein", "B", "T_NK", 5),
        _lineage(grid, "pancreas", "protein", "B", "B_Plasma", 45),
        _lineage(grid, "pancreas", "RNA", "B", "T_NK", 30),
        _lineage(grid, "pancreas", "RNA", "B", "B_Plasma", 10),
    ]
    if with_c:
        rows += [
            _lineage(grid, "pancreas", "protein", "C", "T_NK", 0,
                     state=grid.STATE_UNAVAILABLE, reason="no T_NK class in this lane"),
            _lineage(grid, "pancreas", "protein", "C", "B_Plasma", 50),
            _lineage(grid, "pancreas", "RNA", "C", "T_NK", 25),
            _lineage(grid, "pancreas", "RNA", "C", "B_Plasma", 25),
        ]
    return _comp(rows)


def test_a_null_tnk_never_enters_the_rank_and_the_ranked_donors_do_not_move(grid):
    with_c = _by_donor(grid.crossmodal_tnk(_both_modalities_frame(grid, True)), "pancreas")
    without_c = _by_donor(grid.crossmodal_tnk(_both_modalities_frame(grid, False)), "pancreas")

    # C: T_NK unavailable on the protein side, measured on the RNA side.
    c = with_c.loc["C"]
    assert pd.isna(c["frac_T_NK_protein"]) and c["frac_T_NK_RNA"] == pytest.approx(0.5)
    assert c["state_T_NK_protein"] == grid.STATE_UNAVAILABLE
    assert c["state_T_NK_RNA"] == grid.STATE_MEASURED
    assert not c["both_modalities"], "a donor with a NULL on one side is not measured in both"
    assert pd.isna(c["rank_protein"]) and pd.isna(c["rank_RNA"]), "a NULL has no rank"
    assert pd.isna(c["rank_disagrees"])
    assert pd.isna(c["offset_pp_RNA_minus_protein"])

    # A and B: exactly the arithmetic, and identical whether or not C is in the frame.
    for donor in ("A", "B"):
        pd.testing.assert_series_equal(
            with_c.loc[donor].drop("tissue"), without_c.loc[donor].drop("tissue"),
            check_names=False, obj=f"donor {donor}")
    assert with_c.at["A", "frac_T_NK_protein"] == pytest.approx(0.25)
    assert with_c.at["B", "frac_T_NK_protein"] == pytest.approx(0.10)
    assert with_c.at["A", "frac_T_NK_RNA"] == pytest.approx(0.50)
    assert with_c.at["B", "frac_T_NK_RNA"] == pytest.approx(0.75)
    assert with_c.at["A", "rank_protein"] == 2 and with_c.at["B", "rank_protein"] == 1
    assert with_c.at["A", "rank_RNA"] == 1 and with_c.at["B", "rank_RNA"] == 2
    assert with_c.at["A", "rank_disagrees"] is np.True_ or with_c.at["A", "rank_disagrees"] == True  # noqa: E712
    assert with_c.at["B", "rank_disagrees"] == True  # noqa: E712
    # The rank was taken over the two donors measured in both, not three.
    assert with_c["rank_protein"].max() == 2


def test_a_donor_unavailable_on_both_sides_keeps_its_row_and_is_never_compared(grid):
    comp = _comp([
        _lineage(grid, "pln", "protein", "D", "T_NK", 0,
                 state=grid.STATE_UNAVAILABLE, reason="protein: no T_NK class"),
        _lineage(grid, "pln", "protein", "D", "B_Plasma", 40),
        _lineage(grid, "pln", "RNA", "D", "T_NK", 0,
                 state=grid.STATE_BLANK, reason="RNA: never typed"),
        _lineage(grid, "pln", "RNA", "D", "B_Plasma", 60),
    ])
    wide = _by_donor(grid.crossmodal_tnk(comp), "pln")
    assert list(wide.index) == ["D"], "a row with a NULL on both sides is still a row"
    assert pd.isna(wide.at["D", "frac_T_NK_protein"]) and pd.isna(wide.at["D", "frac_T_NK_RNA"])
    assert wide.at["D", "state_T_NK_protein"] == grid.STATE_UNAVAILABLE
    assert wide.at["D", "state_T_NK_RNA"] == grid.STATE_BLANK
    assert wide.at["D", "n_resolved_immune_protein"] == 40
    assert wide.at["D", "n_resolved_immune_RNA"] == 60
    assert not wide.at["D", "both_modalities"]
    assert pd.isna(wide.at["D", "rank_disagrees"])


# ---------------------------------------------------------------------------
# Two rows for one grid cell AND one lineage is a contradiction in the composition
# table.  It is refused by key -- never resolved by taking the first row (what
# ``drop_duplicates`` did until 2026-09-21) and never by summing.  The check runs
# over every immune lineage BEFORE the groupby-sum that builds ``n_immune``: until
# 2026-09-22 only the T_NK rows were checked, so a duplicated B_Plasma row was
# added into the denominator silently.
# ---------------------------------------------------------------------------
def test_a_duplicated_tnk_key_is_refused_by_name_not_resolved_to_the_first_row(grid):
    comp = _comp([
        _lineage(grid, "pancreas", "protein", "A", "T_NK", 10),
        _lineage(grid, "pancreas", "protein", "A", "T_NK", 30),
        _lineage(grid, "pancreas", "protein", "A", "B_Plasma", 60),
        _lineage(grid, "pancreas", "protein", "B", "T_NK", 5),
        _lineage(grid, "pancreas", "protein", "B", "B_Plasma", 45),
    ])
    with pytest.raises(grid.GridError) as err:
        grid.crossmodal_tnk(comp)
    msg = str(err.value)
    assert "('pancreas', 'protein', 'A', 'T_NK')" in msg, (
        "the refusal names the duplicated key down to the lineage")
    assert "'B'" not in msg, "a donor with one T_NK row is not named"
    assert "B_Plasma" not in msg, "a lineage with one row is not named"


def test_a_duplicated_non_tnk_immune_key_is_refused_before_it_is_summed(grid):
    # Two B_Plasma rows for one grid cell. Until 2026-09-22 the groupby-sum added
    # them into n_immune = 70 and no refusal fired; only T_NK was checked.
    comp = _comp([
        _lineage(grid, "pancreas", "protein", "A", "T_NK", 10),
        _lineage(grid, "pancreas", "protein", "A", "B_Plasma", 30),
        _lineage(grid, "pancreas", "protein", "A", "B_Plasma", 30),
        _lineage(grid, "pancreas", "protein", "B", "T_NK", 5),
        _lineage(grid, "pancreas", "protein", "B", "B_Plasma", 45),
    ])
    with pytest.raises(grid.GridError) as err:
        grid.crossmodal_tnk(comp)
    msg = str(err.value)
    assert "('pancreas', 'protein', 'A', 'B_Plasma')" in msg, (
        "the refusal names the duplicated key down to the lineage")
    assert "'B'" not in msg and "T_NK" not in msg, "only the duplicated key is named"


def test_a_duplicated_non_immune_key_does_not_reach_the_immune_check(grid):
    # The check is on the immune lineage rows the sum reads; a duplicated
    # non-immune row never enters n_immune and is not this table's contradiction.
    comp = _comp([
        _lineage(grid, "pancreas", "protein", "A", "T_NK", 10),
        _lineage(grid, "pancreas", "protein", "A", "B_Plasma", 30),
        _lineage(grid, "pancreas", "protein", "A", "Endocrine", 100, is_immune="no"),
        _lineage(grid, "pancreas", "protein", "A", "Endocrine", 100, is_immune="no"),
    ])
    wide = _by_donor(grid.crossmodal_tnk(comp), "pancreas")
    assert wide.at["A", "n_resolved_immune_protein"] == 40


# ---------------------------------------------------------------------------
# A row in state ``measured`` carries a count.  A NULL count under that state is a
# contradiction, and the groupby-sum that builds ``n_immune`` would skip it silently.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("null_lineage", ["B_Plasma", "T_NK"])
def test_a_measured_row_with_a_null_count_is_refused(grid, null_lineage):
    comp = _comp([
        _lineage(grid, "pln", "RNA", "A", "T_NK", 10),
        _lineage(grid, "pln", "RNA", "A", "B_Plasma", 30),
        _lineage(grid, "pln", "RNA", "B", "T_NK", 5),
        _lineage(grid, "pln", "RNA", "B", "B_Plasma", 45),
    ])
    hit = comp["donor_id"].eq("A") & comp["lineage_common"].eq(null_lineage)
    comp.loc[hit, "n_cells"] = pd.NA          # the state stays "measured": the contradiction
    assert comp.loc[hit, "state"].eq(grid.STATE_MEASURED).all()
    with pytest.raises(grid.GridError) as err:
        grid.crossmodal_tnk(comp)
    msg = str(err.value)
    assert grid.STATE_MEASURED in msg and "n_cells" in msg
    assert "('pln', 'RNA', 'A', " in msg and null_lineage in msg, "the refusal names the row"
    assert "'B'" not in msg
