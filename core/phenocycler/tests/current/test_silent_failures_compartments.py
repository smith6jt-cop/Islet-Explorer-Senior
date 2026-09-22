"""Places where a missing measurement used to read as a number.

Each test is a two-row synthetic frame: one row carrying the value, one row with
it null or the column absent.  The row that HAD a measurement must come out
numerically identical to what the builder produced before the fix; the row that
did not must come out unavailable, counted under its own name, and never as a
zero, a False or an unreduced denominator.  No number here is biological.

Run from core/phenocycler with
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider -q <this file>.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
BUILDER = REPO / "scripts/senior/build_compartment_gradient.py"
spec = importlib.util.spec_from_file_location("build_compartment_gradient", BUILDER)
builder = importlib.util.module_from_spec(spec)
sys.modules["build_compartment_gradient"] = builder
spec.loader.exec_module(builder)

sys.path.insert(0, str(REPO))
from scripts.senior import build_availability as ledger  # noqa: E402


def _rows_by_lineage(rows):
    return {r["lineage"]: r for r in rows}


# ---------------------------------------------------------------------------
# site 1 -- a null immune lineage is unavailable, not non-immune


def test_a_null_immune_lineage_leaves_the_denominator_and_is_counted_by_name():
    df = pd.DataFrame({
        "is_immune": pd.array([True, None], dtype="boolean"),
        "common": ["T_NK", None],
    })
    lineage_unavailable = df["is_immune"].isna()
    rows = builder._composition_rows(
        "D", "pancreatic_lymph_node", "phenocycler", "pln",
        df.loc[~lineage_unavailable].assign(is_immune=lambda d: d["is_immune"].astype(bool)),
        "not_applicable_no_islets",
        n_lineage_unavailable=int(lineage_unavailable.sum()),
        lineage_unavailable_reason="is_immune is null",
    )
    total = _rows_by_lineage(rows)["Immune_total"]
    # The measured row is unchanged: one immune cell, one immune cell in the count.
    assert total["n_immune"] == 1
    assert total["n_lineage"] == 1
    # The null row is not in the denominator ...
    assert total["n_cells"] == 1
    assert total["fraction_of_compartment"] == 1.0
    # ... and it is reported as its own quantity, with a reason, on every row.
    assert all(r["n_lineage_unavailable"] == 1 for r in rows)
    assert all(r["lineage_unavailable_reason"] == "is_immune is null" for r in rows)
    assert total["state"] == "measured"


def test_rows_that_carried_a_lineage_are_numerically_identical_to_before():
    # Pre-fix arithmetic on a frame with no null: n_cells 2, n_immune 1, 0.5.
    df = pd.DataFrame({"is_immune": [True, False], "common": ["Myeloid", None]})
    rows = builder._composition_rows("D", "pancreas", "phenocycler", "core", df, "all_islets")
    by = _rows_by_lineage(rows)
    assert by["Immune_total"]["n_cells"] == 2
    assert by["Immune_total"]["n_immune"] == 1
    assert by["Immune_total"]["fraction_of_compartment"] == 0.5
    assert by["Myeloid"]["n_lineage"] == 1
    assert by["Myeloid"]["fraction_of_immune"] == 1.0
    assert by["T_NK"]["n_lineage"] == 0
    # The new column is a measured zero here: every row carried a lineage.
    assert all(r["n_lineage_unavailable"] == 0 for r in rows)
    assert all(r["lineage_unavailable_reason"] == "" for r in rows)


def _pln_tree(tmp_path, cells: dict, *, status_row=None, counts_row=None) -> Path:
    root = tmp_path / "pheno"
    part = root / "donor_id=D" / "roi=pln_1"
    part.mkdir(parents=True)
    pd.DataFrame(cells).to_parquet(part / "data.parquet", index=False)
    status_row = status_row or {"donor_id": "D", "state": "typed", "n_eligible": len(cells["is_immune"]),
                                "n_clusters": 1, "n_immune": "", "reason": ""}
    pd.DataFrame([status_row]).to_csv(root / "per_donor_status.csv", index=False)
    if counts_row is not None:
        pd.DataFrame([counts_row]).to_csv(root / "per_donor_counts.csv", index=False)
    return root


def test_pln_partition_with_a_null_is_immune_reports_it_and_keeps_it_out(tmp_path, monkeypatch):
    root = _pln_tree(tmp_path, {
        "is_immune": pd.array([True, None], dtype="boolean"),
        "immune_class": ["T_cell", ""],
        "qc_analysis_eligible": [True, True],
        "excluded_reason": ["", ""],
    })
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    total = _rows_by_lineage(rows)["Immune_total"]
    assert total["n_immune"] == 1, "the measured immune cell is unchanged"
    assert total["n_cells"] == 1, "the null-lineage cell is not in the denominator"
    assert total["n_lineage_unavailable"] == 1
    assert total["lineage_unavailable_reason"]
    assert note["n_lineage_unavailable"] == 1
    assert note["n_cells"] == 1


# ---------------------------------------------------------------------------
# site 2 -- the pLN protein path applies the pancreas eligibility filter


def test_pln_partition_carrying_eligibility_columns_is_filtered_like_the_pancreas(tmp_path, monkeypatch):
    root = _pln_tree(tmp_path, {
        "is_immune": pd.array([True, True], dtype="boolean"),
        "immune_class": ["T_cell", "T_cell"],
        "qc_analysis_eligible": pd.array([True, None], dtype="boolean"),
        "excluded_reason": ["", "missing_nucleus_geometry"],
    })
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    total = _rows_by_lineage(rows)["Immune_total"]
    assert total["n_cells"] == 1, "the ineligible cell leaves every calculation"
    assert total["n_immune"] == 1
    assert note["n_rows"] == 2
    assert note["n_eligible"] == 1
    # Counted the same way the pancreas path counts it: by its recorded reason.
    assert note["dropped"] == {"missing_nucleus_geometry": 1}
    assert note["eligibility_filter"].startswith("row-level")


def test_pln_partition_without_eligibility_columns_names_the_upstream_filter(tmp_path, monkeypatch):
    # Today's partitions carry neither column: the producer filtered on
    # qc_analysis_eligible before writing and recorded the counts beside the tree.
    root = _pln_tree(
        tmp_path,
        {"is_immune": pd.array([True, False], dtype="boolean"), "immune_class": ["T_cell", ""]},
        status_row={"donor_id": "D", "state": "typed", "n_eligible": 2, "n_clusters": 1,
                    "n_immune": "", "reason": ""},
        counts_row={"donor_id": "D", "n_rows": 5, "n_eligible": 2, "n_clustered": 2},
    )
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    total = _rows_by_lineage(rows)["Immune_total"]
    assert total["state"] == "measured"
    assert total["n_cells"] == 2 and total["n_immune"] == 1
    assert note["n_rows"] == 5
    assert note["n_eligible"] == 2
    assert note["dropped"] == {"ineligible_recorded_by_producer": 3}
    assert note["eligibility_filter"].startswith("upstream")


def test_pln_partition_holding_more_rows_than_the_recorded_eligible_count_is_unavailable(tmp_path, monkeypatch):
    root = _pln_tree(
        tmp_path,
        {"is_immune": pd.array([True, False], dtype="boolean"), "immune_class": ["T_cell", ""]},
        status_row={"donor_id": "D", "state": "typed", "n_eligible": 1, "n_clusters": 1,
                    "n_immune": "", "reason": ""},
    )
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    assert rows and all(r["state"] == "unavailable" for r in rows)
    assert all(pd.isna(r["n_cells"]) for r in rows)
    assert all("n_eligible" in r["reason"] for r in rows)
    assert note["state"] == "unavailable"


# ---------------------------------------------------------------------------
# site 3 -- a Xenium summary without no_nucleus_excluded is unavailable, never unreduced


class _FakeRelease:
    available = True
    release_id = "fake-release"
    roots = ()
    rejected = ()

    def __init__(self, root: Path):
        self.root = root

    def sidecar(self, name):
        p = self.root / name
        return p if p.exists() else None

    def describe(self):
        return {"release_id": self.release_id}


def _xenium_rows(tmp_path, monkeypatch, summary: pd.DataFrame):
    root = tmp_path / "release"
    root.mkdir()
    summary.to_csv(root / "per_donor_summary.csv", index=False)
    monkeypatch.setattr(builder, "xenium_source", lambda tissue: _FakeRelease(root))
    return builder.pln_xenium()


def test_xenium_donor_with_a_null_no_nucleus_excluded_is_unavailable(tmp_path, monkeypatch):
    summary = pd.DataFrame({
        "donor_id": ["A", "B"], "cells": [10, 10], "immune": [4, 4],
        "no_nucleus_excluded": [2, None],
        "T_NK_compartment": [3, 3], "Immune_B_Plasma": [1, 1], "Myeloid": [0, 0],
        "Granulocyte": [0, 0],
    })
    rows, note = _xenium_rows(tmp_path, monkeypatch, summary)
    a = _rows_by_lineage([r for r in rows if r["donor_id"] == "A"])
    b = _rows_by_lineage([r for r in rows if r["donor_id"] == "B"])
    # The donor that reported the exclusion is unchanged: 10 - 2 = 8.
    assert a["Immune_total"]["n_cells"] == 8
    assert a["Immune_total"]["n_immune"] == 4
    assert a["Immune_total"]["fraction_of_compartment"] == 0.5
    assert a["T_NK"]["n_lineage"] == 3
    assert a["Immune_total"]["state"] == "measured"
    # The donor that did not is unavailable on every row -- never 10.
    assert set(b) == {"Immune_total", *builder.COMMON_LINEAGES}
    assert all(r["state"] == "unavailable" for r in b.values())
    assert all(pd.isna(r["n_cells"]) for r in b.values())
    assert all(r["reason"] == "no_nucleus_excluded not reported" for r in b.values())
    assert note["donors_without_no_nucleus_excluded"] == ["B"]


def test_xenium_summary_lacking_the_column_altogether_is_unavailable(tmp_path, monkeypatch):
    summary = pd.DataFrame({
        "donor_id": ["A", "B"], "cells": [10, 10], "immune": [4, 4],
        "no_nucleus_excluded_": [2, 2],  # a misspelling is an absent column
        "T_NK_compartment": [3, 3], "Immune_B_Plasma": [1, 1], "Myeloid": [0, 0],
        "Granulocyte": [0, 0],
    })
    rows, note = _xenium_rows(tmp_path, monkeypatch, summary)
    assert len(rows) == 2 * (1 + len(builder.COMMON_LINEAGES))
    assert all(r["state"] == "unavailable" for r in rows)
    assert all(pd.isna(r["n_cells"]) for r in rows)
    assert all(r["reason"] == "no_nucleus_excluded not reported" for r in rows)
    assert note["donors_without_no_nucleus_excluded"] == ["A", "B"]


# ---------------------------------------------------------------------------
# site 4 -- an all-null qc_analysis_eligible column is unmeasured, not zero


def _typing_root(tmp_path, values_by_donor: dict) -> Path:
    root = tmp_path / "cells"
    for donor, values in values_by_donor.items():
        part = root / f"donor_id={donor}"
        part.mkdir(parents=True)
        pq.write_table(pa.table({"qc_analysis_eligible": pa.array(values, type=pa.bool_())}),
                       part / "data.parquet")
    return root


def test_all_null_qc_eligible_is_unmeasured_with_its_own_reason(tmp_path):
    root = _typing_root(tmp_path, {"112": [True, None], "115": [None, None]})
    rows = ledger.typing_availability(root)
    # One value present: the count is the measured value, unchanged.
    assert rows["112"]["state"] == "measured"
    assert rows["112"]["n"] == 1
    assert rows["112"]["reason"] == "qc_analysis_eligible_cells; n_null=1"
    # No value present anywhere: the column was never measured for this donor.
    assert rows["115"]["state"] == "unavailable"
    assert rows["115"]["n"] is None
    assert rows["115"]["reason"] == "qc_analysis_eligible_unmeasured; n_null=2"


def test_a_measured_zero_keeps_its_current_reason(tmp_path):
    root = _typing_root(tmp_path, {"112": [False, False]})
    rows = ledger.typing_availability(root)
    assert rows["112"]["state"] == "unavailable"
    assert rows["112"]["n"] is None
    assert rows["112"]["reason"] == "no_qc_eligible_typing_cells; n_null=0"


# ---------------------------------------------------------------------------
# site 1, pancreas half -- _composition_rows is shared, so the pancreas path
# counts a null cell_class the same way (0 of 25,230,109 on disk are null, so
# no pancreas number changes today; the column must still mean one thing)


def _pancreas_tree(tmp_path, cells: dict) -> Path:
    root = tmp_path / "cells"
    part = root / "donor_id=D"
    part.mkdir(parents=True)
    pd.DataFrame(cells).to_parquet(part / "data.parquet", index=False)
    return root


_PANC_COLUMNS = dict(
    object_id=["a", "b"], immune_definition=["Imm-T-cell", ""], excluded_reason=["", ""],
    qc_analysis_eligible=[True, True], signed_dist_islet_um=[-1.0, -1.0], islet_num=[1, 1],
    cell_region_20um=["core", "core"], acq_run=["r1", "r1"], stain_group=["s1", "s1"],
)


def test_pancreas_null_cell_class_leaves_the_denominator_of_its_compartment(tmp_path, monkeypatch):
    root = _pancreas_tree(tmp_path, {**_PANC_COLUMNS, "cell_class": ["Immune", None]})
    monkeypatch.setattr(builder, "CELLS", root)
    monkeypatch.setattr(builder, "MEMBERSHIP", tmp_path / "no_membership")
    res = builder.pancreas_donor("D", min_islet_cells=1)
    core = _rows_by_lineage([r for r in res["rows"] if r["compartment"] == "core"])
    assert core["Immune_total"]["n_immune"] == 1, "the classed cell is unchanged"
    assert core["Immune_total"]["n_cells"] == 1, "the unclassed cell is not in the denominator"
    assert core["Immune_total"]["n_lineage_unavailable"] == 1
    assert core["Immune_total"]["lineage_unavailable_reason"]
    exo = _rows_by_lineage([r for r in res["rows"] if r["compartment"] == "exocrine"])
    assert exo["Immune_total"]["n_lineage_unavailable"] == 0, "counted in its own compartment only"
    assert res["diagnostics"]["n_lineage_unavailable"] == 1
    # The size stratum carries it too, restricted to the islets in the stratum.
    strat = _rows_by_lineage([r for r in res["strat_rows"] if r["compartment"] == "core"])
    assert strat["Immune_total"]["n_lineage_unavailable"] == 1
    # Every row is rectangular once framed: the card reads columns by name.
    framed = pd.DataFrame(res["rows"])
    assert {"n_lineage_unavailable", "lineage_unavailable_reason"} <= set(framed.columns)


def test_pancreas_rows_with_every_cell_classed_are_numerically_identical_to_before(tmp_path, monkeypatch):
    root = _pancreas_tree(tmp_path, {**_PANC_COLUMNS, "cell_class": ["Immune", "Exocrine"]})
    monkeypatch.setattr(builder, "CELLS", root)
    monkeypatch.setattr(builder, "MEMBERSHIP", tmp_path / "no_membership")
    res = builder.pancreas_donor("D", min_islet_cells=1)
    core = _rows_by_lineage([r for r in res["rows"] if r["compartment"] == "core"])
    # Pre-fix arithmetic: 2 eligible core cells, 1 immune, 0.5.
    assert core["Immune_total"]["n_cells"] == 2
    assert core["Immune_total"]["n_immune"] == 1
    assert core["Immune_total"]["fraction_of_compartment"] == 0.5
    assert core["T_NK"]["n_lineage"] == 1
    assert all(r["n_lineage_unavailable"] == 0 for r in res["rows"] + res["strat_rows"])
    assert res["diagnostics"]["n_eligible"] == 2


# ---------------------------------------------------------------------------
# site 5 -- a pLN is_immune column that is not boolean-valued is refused by dtype
#
# ``.astype(bool)`` on a string column reads every non-empty value as True and on
# an int column every non-zero value as True, so ['True', 'NA'] and [1, -1] each
# came out as two immune cells, measured.  The pancreas path raises for a
# partition it cannot read (read_parquet_retrying); this does the same, by name.


def test_pln_bool_is_immune_column_is_read_as_before(tmp_path, monkeypatch):
    root = _pln_tree(tmp_path, {
        "is_immune": [True, False],  # plain numpy bool, the dtype on disk for donor 115
        "immune_class": ["T_cell", ""],
    })
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    total = _rows_by_lineage(rows)["Immune_total"]
    assert total["state"] == "measured"
    assert total["n_cells"] == 2 and total["n_immune"] == 1
    assert total["fraction_of_compartment"] == 0.5
    assert note["n_cells"] == 2


def test_pln_arrow_bool_with_nulls_and_no_pandas_metadata_is_still_boolean(tmp_path, monkeypatch):
    # The producer writes a nullable boolean; a table written by pyarrow alone
    # carries no pandas dtype metadata and reads back as object.  Boolean-valued
    # is boolean-valued whichever way it was written: the null leaves, the True
    # is counted, and nothing is refused.
    root = tmp_path / "pheno"
    part = root / "donor_id=D" / "roi=pln_1"
    part.mkdir(parents=True)
    pq.write_table(pa.table({"is_immune": pa.array([True, None], type=pa.bool_()),
                             "immune_class": pa.array(["T_cell", ""])}), part / "data.parquet")
    pd.DataFrame([{"donor_id": "D", "state": "typed", "n_eligible": 2, "n_clusters": 1,
                   "n_immune": "", "reason": ""}]).to_csv(root / "per_donor_status.csv", index=False)
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    total = _rows_by_lineage(rows)["Immune_total"]
    assert total["n_cells"] == 1 and total["n_immune"] == 1
    assert total["n_lineage_unavailable"] == 1


@pytest.mark.parametrize("values", [
    pytest.param(["True", "NA"], id="string"),
    pytest.param(pd.array([1, -1], dtype="int8"), id="int8"),
])
def test_pln_non_boolean_is_immune_column_is_refused_naming_its_dtype(tmp_path, monkeypatch, values):
    root = _pln_tree(tmp_path, {"is_immune": values, "immune_class": ["T_cell", "T_cell"]})
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    on_disk = str(pd.read_parquet(root / "donor_id=D" / "roi=pln_1" / "data.parquet")["is_immune"].dtype)
    with pytest.raises(ValueError, match="is_immune") as excinfo:
        builder.pln_phenocycler("D")
    assert on_disk in str(excinfo.value), "the refusal names the dtype it found"
    assert "not boolean" in str(excinfo.value)


# ---------------------------------------------------------------------------
# site 6 -- a pancreas cell_class that is string-shaped missing, or outside the
# declared vocabulary, is lineage-unavailable, not non-immune
#
# ``cell_class`` on disk carries the KEYS of build_cells.DISPLAY ('Alpha',
# 'Blood vessel', 'Nerve', ...): build_cells.py:394 maps those keys to the
# display labels, and build_cells.py:396 reports a value the table lacks.  That
# table is the one vocabulary; the builder derives its own from it.


def _typing_build_cells():
    sys.path.insert(0, str(REPO / "core" / "phenocycler" / "scripts" / "typing"))
    import build_cells  # noqa: E402
    return build_cells


def test_pancreas_class_vocabulary_is_the_declared_table_not_a_restated_list():
    cells = _typing_build_cells()
    assert builder.PANC_CLASS_VOCABULARY == frozenset(cells.DISPLAY)
    assert {"Immune", "Exocrine", "Blood vessel", "Nerve", "Lymphatic",
            "Alpha", "Beta", "Delta"} <= builder.PANC_CLASS_VOCABULARY
    # The declared absence is in the table so it can be displayed, and it is not a lineage.
    assert builder.PANC_LINEAGE_ABSENT == frozenset({cells.UNRESOLVED})


@pytest.mark.parametrize("bad", ["", "nan", "None", "NA"], ids=repr)
def test_pancreas_string_shaped_missing_cell_class_is_lineage_unavailable(tmp_path, monkeypatch, bad):
    root = _pancreas_tree(tmp_path, {**_PANC_COLUMNS, "cell_class": ["Immune", bad]})
    monkeypatch.setattr(builder, "CELLS", root)
    monkeypatch.setattr(builder, "MEMBERSHIP", tmp_path / "no_membership")
    res = builder.pancreas_donor("D", min_islet_cells=1)
    core = _rows_by_lineage([r for r in res["rows"] if r["compartment"] == "core"])
    assert core["Immune_total"]["n_immune"] == 1, "the classed cell is unchanged"
    assert core["Immune_total"]["n_cells"] == 1, "the missing-class cell is not in the denominator"
    assert core["Immune_total"]["fraction_of_compartment"] == 1.0
    assert core["Immune_total"]["n_lineage_unavailable"] == 1
    assert "string" in core["Immune_total"]["lineage_unavailable_reason"]
    assert res["diagnostics"]["n_lineage_unavailable"] == 1
    assert res["diagnostics"]["lineage_unavailable_by_kind"]["missing_string"] == 1


def test_pancreas_cell_class_outside_the_declared_vocabulary_is_lineage_unavailable(tmp_path, monkeypatch):
    root = _pancreas_tree(tmp_path, {**_PANC_COLUMNS, "cell_class": ["Immune", "Bogus"]})
    monkeypatch.setattr(builder, "CELLS", root)
    monkeypatch.setattr(builder, "MEMBERSHIP", tmp_path / "no_membership")
    res = builder.pancreas_donor("D", min_islet_cells=1)
    core = _rows_by_lineage([r for r in res["rows"] if r["compartment"] == "core"])
    assert core["Immune_total"]["n_immune"] == 1
    assert core["Immune_total"]["n_cells"] == 1, "'Bogus' is not a non-immune cell"
    assert core["Immune_total"]["n_lineage_unavailable"] == 1
    assert "Bogus" in core["Immune_total"]["lineage_unavailable_reason"]
    assert "vocabulary" in core["Immune_total"]["lineage_unavailable_reason"]
    assert res["diagnostics"]["lineage_unavailable_by_kind"]["outside_vocabulary"] == 1
    assert res["diagnostics"]["cell_class_values_outside_vocabulary"] == ["Bogus"]


def test_pancreas_every_declared_non_immune_class_still_counts_as_non_immune(tmp_path, monkeypatch):
    # The eight classes are unchanged: a declared non-immune class is in the
    # denominator exactly as before, and only 'Immune' is immune.
    root = _pancreas_tree(tmp_path, {**_PANC_COLUMNS, "cell_class": ["Immune", "Blood vessel"]})
    monkeypatch.setattr(builder, "CELLS", root)
    monkeypatch.setattr(builder, "MEMBERSHIP", tmp_path / "no_membership")
    res = builder.pancreas_donor("D", min_islet_cells=1)
    core = _rows_by_lineage([r for r in res["rows"] if r["compartment"] == "core"])
    assert core["Immune_total"]["n_cells"] == 2 and core["Immune_total"]["n_immune"] == 1
    assert core["Immune_total"]["n_lineage_unavailable"] == 0
    assert core["Immune_total"]["lineage_unavailable_reason"] == ""


# ---------------------------------------------------------------------------
# site 7 -- a partially null qc_analysis_eligible column carries its null count
#
# availability() (build_availability.py:78-88) forbids n on a non-measured row,
# and the layer's six columns are pinned by test_availability_layer.py:29, so
# the count travels in the reason text, with its number, on EVERY typing row.


def test_partially_null_qc_eligible_reports_the_null_count_beside_n(tmp_path):
    root = _typing_root(tmp_path, {
        "111": [False, None, None, None],   # a measured zero with three unmeasured rows
        "112": [True, None, None],          # a measured one with two unmeasured rows
        "113": [None, None],                # nothing measured
        "114": [True, False],               # everything measured
    })
    rows = ledger.typing_availability(root)
    assert rows["111"]["state"] == "unavailable" and rows["111"]["n"] is None
    assert rows["111"]["reason"] == "no_qc_eligible_typing_cells; n_null=3"
    assert rows["112"]["state"] == "measured" and rows["112"]["n"] == 1
    assert rows["112"]["reason"] == "qc_analysis_eligible_cells; n_null=2"
    assert rows["113"]["state"] == "unavailable" and rows["113"]["n"] is None
    assert rows["113"]["reason"] == "qc_analysis_eligible_unmeasured; n_null=2"
    assert rows["114"]["state"] == "measured" and rows["114"]["n"] == 1
    assert rows["114"]["reason"] == "qc_analysis_eligible_cells; n_null=0"


# ---------------------------------------------------------------------------
# site 8 -- a pLN partition with FEWER rows than the producer's recorded
# n_eligible is flagged in the note, never carried silently beside a smaller n_cells


def test_pln_partition_holding_fewer_rows_than_the_recorded_eligible_count_is_flagged(tmp_path, monkeypatch):
    root = _pln_tree(
        tmp_path,
        {"is_immune": pd.array([True, False], dtype="boolean"), "immune_class": ["T_cell", ""]},
        status_row={"donor_id": "D", "state": "typed", "n_eligible": 5, "n_clusters": 1,
                    "n_immune": "", "reason": ""},
        counts_row={"donor_id": "D", "n_rows": 7, "n_eligible": 5, "n_clustered": 5},
    )
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    total = _rows_by_lineage(rows)["Immune_total"]
    # What is on disk is measured as before ...
    assert total["state"] == "measured"
    assert total["n_cells"] == 2 and total["n_immune"] == 1
    # ... and the note cannot carry n_eligible=5 beside n_cells=2 without saying so.
    assert note["n_eligible"] == 5
    assert note["n_rows_on_disk"] == 2
    assert note["eligible_count_disagreement"]
    assert "2" in note["eligible_count_disagreement"] and "5" in note["eligible_count_disagreement"]


def test_pln_partition_matching_the_recorded_eligible_count_carries_no_flag(tmp_path, monkeypatch):
    root = _pln_tree(
        tmp_path,
        {"is_immune": pd.array([True, False], dtype="boolean"), "immune_class": ["T_cell", ""]},
        status_row={"donor_id": "D", "state": "typed", "n_eligible": 2, "n_clusters": 1,
                    "n_immune": "", "reason": ""},
        counts_row={"donor_id": "D", "n_rows": 5, "n_eligible": 2, "n_clustered": 2},
    )
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    assert note["n_eligible"] == 2 and note["n_rows_on_disk"] == 2
    assert note["eligible_count_disagreement"] == ""


# ---------------------------------------------------------------------------
# site 9 -- a Xenium summary row without an immune (or cells) count is unavailable


def test_xenium_donor_with_a_null_immune_count_is_unavailable(tmp_path, monkeypatch):
    summary = pd.DataFrame({
        "donor_id": ["A", "B"], "cells": [10, 10], "immune": [4, None],
        "no_nucleus_excluded": [2, 2],
        "T_NK_compartment": [3, 3], "Immune_B_Plasma": [1, 1], "Myeloid": [0, 0],
        "Granulocyte": [0, 0],
    })
    rows, note = _xenium_rows(tmp_path, monkeypatch, summary)
    a = _rows_by_lineage([r for r in rows if r["donor_id"] == "A"])
    b = _rows_by_lineage([r for r in rows if r["donor_id"] == "B"])
    assert a["Immune_total"]["n_cells"] == 8 and a["Immune_total"]["n_immune"] == 4
    assert a["Immune_total"]["state"] == "measured"
    assert set(b) == {"Immune_total", *builder.COMMON_LINEAGES}
    assert all(r["state"] == "unavailable" for r in b.values())
    assert all(pd.isna(r["n_cells"]) and pd.isna(r["n_immune"]) for r in b.values())
    assert all(r["reason"] == "immune not reported" for r in b.values())
    assert note["donors_unavailable"] == {"B": "immune not reported"}


def test_xenium_summary_lacking_the_immune_column_altogether_is_unavailable(tmp_path, monkeypatch):
    summary = pd.DataFrame({
        "donor_id": ["A", "B"], "cells": [10, 10], "immune_": [4, 4],  # a misspelling is absent
        "no_nucleus_excluded": [2, 2],
        "T_NK_compartment": [3, 3], "Immune_B_Plasma": [1, 1], "Myeloid": [0, 0],
        "Granulocyte": [0, 0],
    })
    rows, note = _xenium_rows(tmp_path, monkeypatch, summary)
    assert len(rows) == 2 * (1 + len(builder.COMMON_LINEAGES))
    assert all(r["state"] == "unavailable" for r in rows)
    assert all(r["reason"] == "immune not reported" for r in rows)
    assert note["donors_unavailable"] == {"A": "immune not reported", "B": "immune not reported"}


def test_xenium_donor_with_a_null_cells_count_is_unavailable(tmp_path, monkeypatch):
    summary = pd.DataFrame({
        "donor_id": ["A", "B"], "cells": [10, None], "immune": [4, 4],
        "no_nucleus_excluded": [2, 2],
        "T_NK_compartment": [3, 3], "Immune_B_Plasma": [1, 1], "Myeloid": [0, 0],
        "Granulocyte": [0, 0],
    })
    rows, note = _xenium_rows(tmp_path, monkeypatch, summary)
    b = _rows_by_lineage([r for r in rows if r["donor_id"] == "B"])
    assert all(r["state"] == "unavailable" for r in b.values())
    assert all(r["reason"] == "cells not reported" for r in b.values())
    assert note["donors_unavailable"] == {"B": "cells not reported"}


# ---------------------------------------------------------------------------
# site 10 -- a pLN status row with no recorded n_eligible, and no
# per_donor_counts.csv, records NO eligible count: n_rows is what is on disk,
# not what the producer said was eligible, and it is never written as n_eligible


@pytest.mark.parametrize("recorded", ["", "abc"], ids=["blank", "unparseable"])
def test_pln_status_without_a_recorded_eligible_count_reports_none_not_n_rows(tmp_path, monkeypatch, recorded):
    root = _pln_tree(
        tmp_path,
        {"is_immune": pd.array([True, False], dtype="boolean"), "immune_class": ["T_cell", ""]},
        status_row={"donor_id": "D", "state": "typed", "n_eligible": recorded, "n_clusters": 1,
                    "n_immune": "", "reason": ""},
    )
    monkeypatch.setattr(builder, "PLN_PHENO", root)
    rows, note = builder.pln_phenocycler("D")
    total = _rows_by_lineage(rows)["Immune_total"]
    # What is on disk is measured as before ...
    assert total["state"] == "measured"
    assert total["n_cells"] == 2 and total["n_immune"] == 1
    assert note["n_rows_on_disk"] == 2
    # ... and the eligible count nobody recorded is absent, by name, not 2.
    assert note["n_eligible"] is None
    assert "n_eligible" in note["dropped_reason"]
    assert "absent" in note["dropped_reason"] or "not recorded" in note["dropped_reason"]
    # The disagreement check reads nothing it was not given.
    assert note["eligible_count_disagreement"] == ""


# ---------------------------------------------------------------------------
# site 11 -- a Xenium per_donor_summary.csv count is refused at the read when
# it cannot be a count: negative, fractional, non-numeric, or a decomposition
# that does not hold (no_nucleus_excluded <= cells; immune and each lineage
# <= nucleated cells).  Donor A carries the same consistent row in every case
# and must come out numerically identical to before.


_XENIUM_CONSISTENT = {
    "cells": 10, "immune": 4, "no_nucleus_excluded": 2,
    "T_NK_compartment": 3, "Immune_B_Plasma": 1, "Myeloid": 0, "Granulocyte": 0,
}


@pytest.mark.parametrize("bad, named", [
    ({"immune": -1}, ["immune", "-1"]),
    ({"immune": 4.7}, ["immune", "4.7"]),
    ({"cells": "ten"}, ["cells", "ten"]),
    ({"no_nucleus_excluded": 12}, ["no_nucleus_excluded", "12", "cells", "10"]),
    ({"immune": 9}, ["immune", "9", "8"]),
    ({"T_NK_compartment": 9}, ["T_NK_compartment", "9", "8"]),
], ids=["negative", "fractional", "non_numeric", "nne_exceeds_cells",
        "immune_exceeds_nucleated", "lineage_exceeds_nucleated"])
def test_xenium_count_that_cannot_be_a_count_is_refused_at_the_read(tmp_path, monkeypatch, bad, named):
    b_row = {**_XENIUM_CONSISTENT, **bad}
    summary = pd.DataFrame({
        "donor_id": ["A", "B"],
        **{c: [_XENIUM_CONSISTENT[c], b_row[c]] for c in _XENIUM_CONSISTENT},
    })
    rows, note = _xenium_rows(tmp_path, monkeypatch, summary)
    a = _rows_by_lineage([r for r in rows if r["donor_id"] == "A"])
    b = _rows_by_lineage([r for r in rows if r["donor_id"] == "B"])
    # The consistent row is unchanged: 10 - 2 = 8 nucleated, 4 immune, 0.5, T_NK 3 of 4.
    assert a["Immune_total"]["state"] == "measured"
    assert a["Immune_total"]["n_cells"] == 8 and a["Immune_total"]["n_immune"] == 4
    assert a["Immune_total"]["fraction_of_compartment"] == 0.5
    assert a["T_NK"]["n_lineage"] == 3 and a["T_NK"]["fraction_of_immune"] == 0.75
    assert a["B_Plasma"]["n_lineage"] == 1
    # The row that cannot be a measurement is unavailable on every lineage ...
    assert set(b) == {"Immune_total", *builder.COMMON_LINEAGES}
    assert all(r["state"] == "unavailable" for r in b.values())
    assert all(pd.isna(r["n_cells"]) and pd.isna(r["n_immune"]) and pd.isna(r["n_lineage"])
               for r in b.values())
    assert all(pd.isna(r["fraction_of_compartment"]) for r in b.values())
    # ... with one reason naming the column and the values.
    reasons = {r["reason"] for r in b.values()}
    assert len(reasons) == 1
    reason = reasons.pop()
    for token in named:
        assert token in reason, (token, reason)
    assert note["donors_unavailable"] == {"B": reason}


def test_xenium_consistent_rows_are_numerically_identical_to_before(tmp_path, monkeypatch):
    # No refusal fires on a frame that is all integers and decomposes.
    summary = pd.DataFrame({"donor_id": ["A", "B"],
                            **{c: [v, v] for c, v in _XENIUM_CONSISTENT.items()}})
    rows, note = _xenium_rows(tmp_path, monkeypatch, summary)
    assert len(rows) == 2 * (1 + len(builder.COMMON_LINEAGES))
    assert all(r["state"] == "measured" for r in rows)
    assert note["donors_unavailable"] == {}


# ---------------------------------------------------------------------------
# site 12 -- typing_availability refuses a qc_analysis_eligible column that is
# not boolean: pc.sum over int8 [1, -1] is 0 and read as a measured zero, and
# over a string column it raises out of the build


@pytest.mark.parametrize("values, dtype", [
    ([1, -1], pa.int8()),
    (["True", "False"], pa.string()),
], ids=["int8", "string"])
def test_non_boolean_qc_eligible_column_is_refused_naming_its_dtype(tmp_path, values, dtype):
    root = tmp_path / "cells"
    for donor, column in (("112", pa.array(values, type=dtype)),
                          ("115", pa.array([True, False], type=pa.bool_()))):
        part = root / f"donor_id={donor}"
        part.mkdir(parents=True)
        pq.write_table(pa.table({"qc_analysis_eligible": column}), part / "data.parquet")
    rows = ledger.typing_availability(root)
    # The column that is not boolean is refused, with its dtype in the reason.
    assert rows["112"]["state"] == "unavailable"
    assert rows["112"]["n"] is None
    assert str(dtype) in rows["112"]["reason"]
    assert "qc_analysis_eligible" in rows["112"]["reason"]
    # The boolean column is read exactly as before.
    assert rows["115"]["state"] == "measured" and rows["115"]["n"] == 1
    assert rows["115"]["reason"] == "qc_analysis_eligible_cells; n_null=0"
