"""What the producer decided must reach the table the app reads, unchanged and undropped.

`build_cells.py` is a join, not a decision.  Everything that can go wrong here is a LOSS: a
column the producer writes and this file never carries, a decision re-derived from something
that looks like it, or two answers to one question merged under one name.  A loss of that kind
is invisible downstream -- a re-derived flag goes on reading the same whatever the chain
decided, and a dropped column is indistinguishable from a decision that was never made.  So the
assertions here are mostly about CARRYING and about REFUSING:

  * every column the producer emits survives the build, with its values;
  * ``is_endocrine`` is the islet chain's decision read off ``resolved_by``, and a cell NAMED
    endocrine with no decision behind it stops the build rather than being counted;
  * a route this file has never been told about stops the build instead of defaulting to
    not-endocrine;
  * a level a stratum could not supply arrives as UNAVAILABLE, never as 0.0;
  * the refinement's verdict is carried BESIDE the class, never over it, and a sidecar from
    another generation of the typing is refused.

The producer's column set is not written down here.  It is obtained by CALLING ``subpanel._emit``,
so a column added upstream is covered the moment it exists.

Run with
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
      core/phenocycler/tests/current/test_build_cells_carry.py -q
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/celltype"))

import build_cells as BC        # noqa: E402
import islet_membership as MEMBERSHIP   # noqa: E402
import subpanel as SUB          # noqa: E402
import subpanel_io as S         # noqa: E402

DONOR = "TESTDONOR"


# --------------------------------------------------------------------------------------------
# A donor, staged on disk, small enough to read by eye.
# --------------------------------------------------------------------------------------------

def producer_frame(cell_class, resolved_by, *, islet_num=None, stratum=None, route=None,
                   cell_region=None, detail=True) -> pd.DataFrame:
    """The producer's OWN table for these cells, built by calling ``subpanel._emit``.

    Never a hand-written column list: the point of the carry test is that a column added to the
    producer is covered without anyone remembering to add it here too.
    """
    n = len(cell_class)
    islet_num = ["1"] * n if islet_num is None else islet_num
    in_islet = np.array([bool(v) for v in islet_num])
    ctx = pd.DataFrame({
        "object_id": [f"uuid-{i:04d}" for i in range(n)],
        "context": np.where(in_islet, "islet", "tissue"),
        "object_context": ["tissue"] * n,
        "cell_region": cell_region if cell_region is not None else (
            np.where(in_islet, "core", "tissue")),
        "islet_num": np.array(islet_num, dtype=object)})
    frame = SUB._emit(
        DONOR, ctx,
        np.array(cell_class, dtype=object), np.array([""] * n, dtype=object),
        np.array(resolved_by, dtype=object), np.array([""] * n, dtype=object),
        np.where(in_islet, 40.0, np.nan),
        np.array(stratum if stratum is not None else np.where(in_islet, ">10", ""), dtype=object),
        np.array(route if route is not None else [""] * n, dtype=object),
        np.array(["FOXP3"] * n, dtype=object),
        {name: np.zeros(n, bool) for name in SUB.STRUCTURES},
        detail)
    # `subpanel.main` drops donor_id before writing: it is the partition key.
    return frame.drop(columns=["donor_id"])


def stage(tmp_path, monkeypatch, calls: pd.DataFrame, *, member=None, cell_region=None,
          distance=None):
    """Write the producer's table and a matching export, and point the module at them."""
    n = len(calls)
    typing = tmp_path / "typing"
    (typing / f"donor_id={DONOR}").mkdir(parents=True)
    calls.to_parquet(typing / f"donor_id={DONOR}/data.parquet", index=False)

    revised = tmp_path / "revised"
    (revised / f"donor_id={DONOR}").mkdir(parents=True)
    pd.DataFrame({
        "object_id": calls["object_id"].to_numpy(),
        "X_centroid": np.arange(n, dtype=float),
        "Y_centroid": np.arange(n, dtype=float),
        "cell_region": (cell_region if cell_region is not None
                        else calls["cell_region"].to_numpy()),
        "qc_analysis_eligible": np.ones(n, bool),
        "analysis_reason": ["ok"] * n,
        "signed_dist_islet_um": (distance if distance is not None else
                                 np.where(calls["islet_num"].notna(), -5.0, 50.0)),
    }).to_parquet(revised / f"donor_id={DONOR}/data.parquet", index=False)

    monkeypatch.setattr(S, "TYPING", typing)
    monkeypatch.setattr(S, "REVISED", revised)
    islets = (member if member is not None
              else calls["islet_num"].fillna("").astype(str).to_numpy(dtype=object))
    monkeypatch.setattr(MEMBERSHIP, "membership",
                        lambda *a, **k: pd.DataFrame({"islet_num": islets,
                                                      "in_islet": islets != ""}))
    return {"pixel": {DONOR: 0.5}, "stem": {DONOR: "stem"},
            "runs": pd.DataFrame({"donor_id": [DONOR], "acq_run": ["r1"],
                                  "stain_group": ["g1"]})}


def build(tmp_path, monkeypatch, calls, **kwargs):
    env = stage(tmp_path, monkeypatch, calls,
                **{k: kwargs.pop(k) for k in ("member", "cell_region", "distance")
                   if k in kwargs})
    return BC.build_donor(DONOR, env["pixel"], env["stem"], env["runs"], **kwargs)


# --------------------------------------------------------------------------------------------
# 1. The carry contract.  Nothing the producer writes may be lost here.
# --------------------------------------------------------------------------------------------

def test_every_column_the_producer_writes_reaches_the_table(tmp_path, monkeypatch):
    calls = producer_frame(["Beta", "Exocrine"], ["cd99", "fallthrough"])
    out = build(tmp_path, monkeypatch, calls)
    renamed = {BC.CARRY_RENAME.get(c, c) for c in calls.columns}
    missing = sorted(renamed - set(out.columns))
    assert not missing, (
        f"{missing} are written by subpanel._emit and do not reach typing/cells. A column that "
        f"exists upstream and is dropped here is indistinguishable downstream from a decision "
        f"that was never made. Produced: {sorted(calls.columns)}")


def test_the_producers_values_survive_verbatim(tmp_path, monkeypatch):
    calls = producer_frame(["Beta", "Immune", "Exocrine"],
                           ["cd99", "definition", "fallthrough"],
                           stratum=[">10", "4-10", ""], route=["cd99", "", ""],
                           islet_num=["7", "7", None])
    out = build(tmp_path, monkeypatch, calls)
    for column in calls.columns:
        if column == "cell_region":
            continue                      # reconciled against the export, checked on its own
        carried = BC.CARRY_RENAME.get(column, column)
        before = calls[column].fillna("").astype(str).tolist()
        after = out[carried].fillna("").astype(str).tolist()
        assert before == after, f"{column} changed value on the way through: {before} -> {after}"


def test_object_id_is_preserved_verbatim_and_every_cell_comes_back(tmp_path, monkeypatch):
    calls = producer_frame(["Beta", "Exocrine", "Immune"],
                           ["cd99", "fallthrough", "definition"])
    out = build(tmp_path, monkeypatch, calls)
    assert out["object_id"].tolist() == calls["object_id"].tolist()
    assert len(out) == len(calls)


def test_a_producer_column_that_collides_with_a_derived_one_stops_the_build(tmp_path,
                                                                            monkeypatch):
    calls = producer_frame(["Beta"], ["cd99"])
    calls["phenotype"] = "Beta_cell"
    with pytest.raises(AssertionError, match="also writes"):
        build(tmp_path, monkeypatch, calls)


# --------------------------------------------------------------------------------------------
# 2. is_endocrine is the DECISION, not the name.
# --------------------------------------------------------------------------------------------

def test_is_endocrine_is_read_off_the_route_the_chain_took(tmp_path, monkeypatch):
    calls = producer_frame(["Beta", "Alpha", "Delta", "Exocrine", "Immune"],
                           ["cd99", "elongation", "hormone", "fallthrough", "definition"])
    out = build(tmp_path, monkeypatch, calls)
    assert out["is_endocrine"].tolist() == [True, True, True, False, False]
    assert out["endocrine_route"].tolist() == ["cd99", "elongation", "hormone", "", ""]


def test_the_three_routes_stay_distinguishable_in_the_table(tmp_path, monkeypatch):
    """A cell decided by CD99 and a cell that kept the hormone question are both endocrine.

    They are not the same fact: the second is in a stratum whose CD99 level could not be fitted,
    so it was never asked the maintainer's key question at all. Collapsing them to one boolean
    is what made the old flag unable to say which cells the chain actually decided.
    """
    calls = producer_frame(["Beta", "Beta"], ["cd99", "hormone"])
    out = build(tmp_path, monkeypatch, calls)
    assert out["is_endocrine"].all()
    assert set(out["endocrine_route"]) == {"cd99", "hormone"}


def test_a_cell_named_endocrine_with_no_decision_behind_it_stops_the_build(tmp_path,
                                                                           monkeypatch):
    calls = producer_frame(["Beta"], ["fallthrough"])
    with pytest.raises(AssertionError, match="no endocrine route"):
        build(tmp_path, monkeypatch, calls)


def test_a_route_written_on_a_cell_that_was_not_named_endocrine_stops_the_build(tmp_path,
                                                                                monkeypatch):
    calls = producer_frame(["Exocrine"], ["cd99"])
    with pytest.raises(AssertionError, match="not named endocrine"):
        build(tmp_path, monkeypatch, calls)


def test_a_route_this_file_was_never_told_about_stops_the_build(tmp_path, monkeypatch):
    """A new route must be ruled endocrine or not, never defaulted to not-endocrine."""
    calls = producer_frame(["Beta"], ["mixture"])
    with pytest.raises(AssertionError, match="never been told"):
        build(tmp_path, monkeypatch, calls)


def test_calls_without_the_decision_column_are_refused(tmp_path, monkeypatch):
    calls = producer_frame(["Beta"], ["cd99"], detail=False)
    assert "resolved_by" not in calls.columns
    with pytest.raises(AssertionError, match="emit-detail"):
        build(tmp_path, monkeypatch, calls)


def test_the_endocrine_routes_are_the_ones_the_producer_writes():
    """The three strings are the producer's, and every one of them is in its source."""
    source = (REPO / "core/phenocycler/scripts/typing/subpanel.py").read_text()
    for route in BC.ENDOCRINE_ROUTES:
        assert f'"{route}"' in source, (
            f"{route!r} is not written anywhere in subpanel.py, so build_cells is counting a "
            "route the producer does not produce")
    assert not set(BC.ENDOCRINE_ROUTES) & set(BC.NON_ENDOCRINE_RESOLVED)
    for value in BC.NON_ENDOCRINE_RESOLVED:
        assert f'"{value}"' in source


# --------------------------------------------------------------------------------------------
# 3. Two answers to one question are reconciled, never merged.
# --------------------------------------------------------------------------------------------

def test_the_producers_islet_number_is_kept_beside_this_files(tmp_path, monkeypatch):
    """The stratum was counted on the producer's islet number; the app reports this file's.

    They are the same call under ``--membership geometry`` and a different call under
    ``--membership parent``. Keeping only one makes it unanswerable whether a cell's stratum was
    counted on the islet the cell is reported in.
    """
    calls = producer_frame(["Beta", "Beta"], ["cd99", "cd99"], islet_num=["7", "7"])
    out = build(tmp_path, monkeypatch, calls, member=np.array(["7", "9"], dtype=object))
    assert out["typing_islet_num"].tolist() == ["7", "7"]
    assert out["islet_num"].tolist() == ["7", "9"]
    assert out["combined_islet_id"].tolist() == [f"{DONOR}_Islet_7", f"{DONOR}_Islet_9"]
    assert out.attrs["islet_num_conflicts"] == 1, (
        "one cell is reported in islet 9 and had its stratum counted in islet 7; that is a "
        "conflict and has to be countable")
    assert out.attrs["islet_num_only_here"] == 0
    assert out.attrs["islet_num_only_in_typing"] == 0


def test_a_peri_islet_cell_the_typing_does_not_number_is_not_a_conflict(tmp_path, monkeypatch):
    """`membership` numbers a cell in the peri-islet expansion; the producer numbers only the
    core.  Measured on donor 115: 72,824 cells, every one of them in the expansion, and 0 cells
    where the two carry different non-empty numbers.  Counting those together would bury the
    one case that IS a defect under a difference that is the definition."""
    calls = producer_frame(["Beta", "Exocrine"], ["cd99", "fallthrough"],
                           islet_num=["7", None])
    out = build(tmp_path, monkeypatch, calls, member=np.array(["7", "7"], dtype=object))
    assert out.attrs["islet_num_only_here"] == 1
    assert out.attrs["islet_num_conflicts"] == 0
    assert out.attrs["islet_num_only_in_typing"] == 0


def test_a_cell_region_the_two_tables_disagree_about_stops_the_build(tmp_path, monkeypatch):
    calls = producer_frame(["Beta"], ["cd99"])
    with pytest.raises(AssertionError, match="disagree about `cell_region`"):
        build(tmp_path, monkeypatch, calls, cell_region=np.array(["tissue"], dtype=object))


# --------------------------------------------------------------------------------------------
# 4. The per-stratum level a cell was measured against.
# --------------------------------------------------------------------------------------------

def _audit(tmp_path, levels) -> Path:
    path = tmp_path / "audit.json"
    path.write_text(json.dumps([{"donor_id": DONOR, "islet_levels": levels}]))
    return path


def test_the_level_each_cell_was_measured_against_reaches_the_cell(tmp_path, monkeypatch):
    calls = producer_frame(["Beta", "Beta"], ["cd99", "cd99"], stratum=[">10", "4-10"])
    audit = _audit(tmp_path, {
        ">10": {"n_confirmed": 8927, "CD99": 0.0745, BC.SHAPE_FEATURE: 1.0890},
        "4-10": {"n_confirmed": 569, "CD99": 0.0324, BC.SHAPE_FEATURE: None}})
    out = build(tmp_path, monkeypatch, calls, audit=audit)
    assert out["islet_cd99_level"].tolist() == [0.0745, 0.0324]
    assert out["islet_elongation_level"].iloc[0] == 1.0890


def test_a_stratum_that_could_supply_no_level_is_unavailable_and_not_zero(tmp_path,
                                                                          monkeypatch):
    """0.0 is not "no level".  A scaled marker clips at its measured background, so a level of
    0.0 is one every measured cell clears -- which is how a whole islet was once called
    endocrine on a number that looked fitted."""
    calls = producer_frame(["Beta", "Beta"], ["hormone", "cd99"], stratum=["2-3", ">10"])
    audit = _audit(tmp_path, {
        "2-3": {"n_confirmed": 131, "CD99": None, BC.SHAPE_FEATURE: None},
        ">10": {"n_confirmed": 12306, "CD99": 0.0508, BC.SHAPE_FEATURE: 1.0896}})
    out = build(tmp_path, monkeypatch, calls, audit=audit)
    assert np.isnan(out["islet_cd99_level"].iloc[0]), (
        "an unfittable stratum must arrive as unavailable; 0.0 would read as a level every "
        "measured cell clears")
    assert out["islet_cd99_level"].iloc[1] == 0.0508
    assert np.isnan(out["islet_elongation_level"].iloc[0])


def test_an_audit_from_another_run_is_refused(tmp_path, monkeypatch):
    calls = producer_frame(["Beta"], ["cd99"])
    path = tmp_path / "audit.json"
    path.write_text(json.dumps([{"donor_id": "6442", "islet_levels": {}}]))
    with pytest.raises(KeyError, match="different run"):
        build(tmp_path, monkeypatch, calls, audit=path)


# --------------------------------------------------------------------------------------------
# 5. The refinement pass reaches the app, beside the class and never over it.
# --------------------------------------------------------------------------------------------

def _sidecar(tmp_path, calls, refined=None, demoted=None, note=None, cell_class=None,
             drop=None) -> Path:
    out = tmp_path / "refine"
    (out / f"donor_id={DONOR}").mkdir(parents=True)
    n = len(calls)
    frame = pd.DataFrame({
        "object_id": calls["object_id"].to_numpy(),
        "donor_id": DONOR,
        "cell_class": (cell_class if cell_class is not None
                       else calls["cell_class"].to_numpy()),
        "refined_class": (refined if refined is not None else calls["cell_class"].to_numpy()),
        "size_stratum": calls["islet_stratum"].to_numpy(),
        "islet_cells": calls["islet_cells"].to_numpy(),
        "immune_leaf": np.zeros(n, bool),
        "concordance": np.full(n, 0.9),
        "concordance_before": np.full(n, 0.8),
        "demoted": demoted if demoted is not None else np.zeros(n, bool),
        "refine_note": note if note is not None else [""] * n,
        "refine_pass": np.ones(n, int)})
    if drop:
        frame = frame.iloc[drop:]
    frame.to_parquet(out / f"donor_id={DONOR}/data.parquet", index=False)
    return out


def test_the_refinements_verdict_reaches_the_table(tmp_path, monkeypatch):
    calls = producer_frame(["Exocrine", "Beta"], ["fallthrough", "cd99"])
    refine = _sidecar(tmp_path, calls, refined=[BC.UNRESOLVED, "Beta"],
                      demoted=np.array([True, False]),
                      note=["below Exocrine bound -1.300 SD in stratum none", ""])
    out = build(tmp_path, monkeypatch, calls, refine=refine)
    for column in BC.REFINE_CARRY:
        assert column in out.columns, f"{column} does not reach typing/cells"
    assert out["refined_class"].tolist() == [BC.UNRESOLVED, "Beta"]
    assert out["refined_phenotype"].tolist() == ["Unresolved", "Beta_cell"]
    assert out["refine_note"].iloc[0].startswith("below Exocrine bound")
    assert out.attrs["refine_demoted"] == 1


def test_the_refinement_never_overwrites_the_class_it_was_run_on(tmp_path, monkeypatch):
    """The pass is a sidecar on purpose: a cell it could not measure keeps its class."""
    calls = producer_frame(["Exocrine", "Beta"], ["fallthrough", "cd99"])
    refine = _sidecar(tmp_path, calls, refined=[BC.UNRESOLVED, "Beta"],
                      demoted=np.array([True, False]))
    out = build(tmp_path, monkeypatch, calls, refine=refine)
    assert out["cell_class"].tolist() == ["Exocrine", "Beta"]
    assert out["phenotype"].tolist() == ["Exocrine", "Beta_cell"]
    assert out["is_endocrine"].tolist() == [False, True], (
        "is_endocrine is the islet chain's decision; the refinement's verdict on the same cell "
        "is `demoted` and `refined_class`, and merging the two would answer a question the "
        "maintainer has not been asked")


def test_a_sidecar_from_another_generation_of_the_typing_is_refused(tmp_path, monkeypatch):
    calls = producer_frame(["Exocrine", "Beta"], ["fallthrough", "cd99"])
    refine = _sidecar(tmp_path, calls, cell_class=["Immune", "Beta"])
    with pytest.raises(AssertionError, match="different generation"):
        build(tmp_path, monkeypatch, calls, refine=refine)


def test_a_typed_cell_with_no_row_in_the_sidecar_is_refused(tmp_path, monkeypatch):
    calls = producer_frame(["Exocrine", "Beta"], ["fallthrough", "cd99"])
    refine = _sidecar(tmp_path, calls, drop=1)
    with pytest.raises(AssertionError, match="no row in the refinement sidecar"):
        build(tmp_path, monkeypatch, calls, refine=refine)


def test_the_stratum_is_not_carried_twice_under_two_names(tmp_path, monkeypatch):
    calls = producer_frame(["Beta"], ["cd99"])
    refine = _sidecar(tmp_path, calls)
    out = build(tmp_path, monkeypatch, calls, refine=refine)
    assert "size_stratum" not in out.columns
    assert "islet_stratum" in out.columns
    assert out.attrs["refine_size_stratum_differs"] == 0


def test_a_sidecar_column_with_no_counterpart_is_not_dropped_for_being_redundant(tmp_path,
                                                                                  monkeypatch):
    """The pre-islet-chain calls carry no stratum, so the sidecar's is the only one there is.

    Dropping a column whose only justification for going is that something else says the same
    thing -- when nothing else says it -- is the loss this whole file exists to prevent. Measured
    on donor 6539 against the deployed calls: `size_stratum` and `islet_cells` come only from the
    sidecar there.
    """
    calls = producer_frame(["Beta"], ["cd99"]).drop(columns=["islet_stratum", "islet_cells"])
    refine = _sidecar(tmp_path, producer_frame(["Beta"], ["cd99"]))
    out = build(tmp_path, monkeypatch, calls, refine=refine)
    assert "size_stratum" in out.columns and out["size_stratum"].iloc[0] == ">10"
    assert "islet_cells" in out.columns and out["islet_cells"].iloc[0] == 40.0


def test_a_stratum_the_two_passes_disagree_about_is_kept_not_dropped(tmp_path, monkeypatch):
    calls = producer_frame(["Beta"], ["cd99"], stratum=[">10"])
    side_view = producer_frame(["Beta"], ["cd99"], stratum=["4-10"])
    refine = _sidecar(tmp_path, side_view)
    out = build(tmp_path, monkeypatch, calls, refine=refine)
    assert out.attrs["refine_size_stratum_differs"] == 1
    assert out["islet_stratum"].iloc[0] == ">10"
    assert out["refine_size_stratum"].iloc[0] == "4-10", (
        "the bound was fitted in one stratum and the endocrine level in another; dropping the "
        "sidecar's would hide the only thing a reader needs")


def test_a_build_without_the_refinement_says_so_on_every_row(tmp_path, monkeypatch):
    """A table with no refinement columns reads exactly like one where the pass demoted nothing.

    The flag cannot be made required -- ``refine.py`` READS ``typing/cells``, so the first build
    of a generation necessarily precedes it -- so the absent case has to be written down instead.
    """
    calls = producer_frame(["Exocrine", "Beta"], ["fallthrough", "cd99"])
    out = build(tmp_path, monkeypatch, calls)
    for column in BC.REFINE_CARRY:
        assert column in out.columns, f"{column} is missing, so the skip is invisible"
    assert (out["refine_note"] == BC.NOT_REFINED).all()
    assert out["refine_pass"].tolist() == [0, 0]
    assert out["refined_class"].isna().all(), (
        "refined_class = cell_class would assert the pass looked at the cell and kept it")
    assert out["demoted"].isna().all(), "demoted = False would be a measured zero"
    assert out["concordance"].isna().all()
    assert out.attrs["refine_applied"] is False
    assert "refine_demoted" not in out.attrs


def test_the_refinement_columns_cannot_arrive_from_the_producer(tmp_path, monkeypatch):
    calls = producer_frame(["Beta"], ["cd99"])
    calls["refined_class"] = "Beta"
    with pytest.raises(AssertionError, match="also writes"):
        build(tmp_path, monkeypatch, calls)


def test_the_unresolved_class_has_a_display_name_and_it_is_refines_own_label():
    refine = pytest.importorskip("refine")
    assert BC.UNRESOLVED == refine.UNRESOLVED
    assert BC.DISPLAY[BC.UNRESOLVED] == "Unresolved"


def test_the_shape_feature_is_the_producers():
    assert BC.SHAPE_FEATURE == SUB.SHAPE_FEATURE


# --------------------------------------------------------------------------------------------
# 6. The corrections this file already made still hold.
# --------------------------------------------------------------------------------------------

def test_the_expansion_is_derived_from_the_signed_distance(tmp_path, monkeypatch):
    calls = producer_frame(["Beta", "Exocrine", "Exocrine"],
                           ["cd99", "fallthrough", "fallthrough"],
                           islet_num=["1", None, None])
    out = build(tmp_path, monkeypatch, calls,
                distance=np.array([-3.0, 19.9, 23.1]))
    assert out["cell_region_20um"].tolist() == ["core", "expansion", "tissue"], (
        "every expansion is 20 um by construction: one donor's annotation is drawn at 23.1 um "
        "and carries the same _exp20um name")


# --------------------------------------------------------------------------------------------
# 7. The served table carries the refinement join (a data-state guard).
# --------------------------------------------------------------------------------------------

@pytest.mark.data
def test_the_served_cells_table_carries_the_refinement_join():
    """Every donor partition of ``typing/cells`` was built with ``--refine``.

    The refinement pass is part of producer 1 since the maintainer's ruling of 2026-09-23, by
    one route only: ``refine.py --out data/development/typing/refine``, then
    ``build_cells.py --include-marginal --refine data/development/typing/refine``, which joins
    the pass's verdict BESIDE ``cell_class`` and leaves the islet chain's columns alone. A build
    without the flag writes the refinement columns as unavailable and records
    ``refine_applied`` False; this guard is what says such a table is not the served one.

    Read off the build ledger the build itself writes, not off 25M rows. RED until the DATA
    step: on 2026-09-23 the ledger lists 3 of the 20 partitions (the last of two invocations on
    2026-09-19 replaced the file whole), each with ``refine_applied`` False. The DATA step runs
    the join once over all twenty, which is also what makes the ledger whole again.
    """
    cells = REPO / "data/development/typing/cells"
    partitions = sorted(p.name.split("=", 1)[1] for p in cells.glob("donor_id=*")
                        if (p / "data.parquet").exists())
    if not partitions:
        pytest.skip(f"no typing/cells partitions under {cells}")
    ledger_path = cells / "build_ledger.csv"
    assert ledger_path.exists(), f"{ledger_path} is missing: the build left no ledger"
    ledger = pd.read_csv(ledger_path, dtype={"donor_id": str}).set_index("donor_id")
    unlisted = [d for d in partitions if d not in ledger.index]
    applied = (ledger["refine_applied"].astype(str).str.lower() == "true"
               if "refine_applied" in ledger.columns
               else pd.Series(False, index=ledger.index))
    not_joined = [d for d in partitions if d in ledger.index and not applied.loc[d]]
    assert not unlisted and not not_joined, (
        f"typing/cells is not the refined build: {len(unlisted)} of {len(partitions)} "
        f"partitions have no ledger row {unlisted} and {len(not_joined)} were built without "
        f"--refine {not_joined}. Run refine.py --out data/development/typing/refine, then "
        "build_cells.py --include-marginal --refine data/development/typing/refine, once over "
        "every donor (the DATA step, on a go).")
