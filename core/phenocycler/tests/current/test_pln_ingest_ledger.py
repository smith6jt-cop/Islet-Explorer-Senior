"""What the pLN ingest ingested, and for every donor it did not, why not.

Two failures this lane has actually had, and one it could have had:

* **A donor the annotation manifest admits went missing and nothing said so.** The 2026-09-17
  run withheld 6579 because the tracked pLN manifest still carried the exclusion `cohort.py` had
  lifted 45 minutes earlier. The withholding was correct on the day; the silence was not. The
  result was 18 partitions where 19 were expected, visible only as a directory that is not there.
  A per-image listing cannot report a donor with no rows, so the ledger is built from the
  manifest's donor list and prints one line per donor either way.

* **"No intensity measurement" is a statement about columns, and it was being carried forward.**
  Four donors held that status from a run that predates nothing -- it is re-measured here, per
  image, per column, every run, and the column census is reported beside the row-level rule so a
  reader can see whether a re-export filled what was empty.

* **There are two copies of the 43.76 GB export on this machine**, same size and same mtime,
  different inodes. Which one a number came from has to be recorded, or it is unanswerable.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("duckdb")

PLN = Path(__file__).resolve().parents[2] / "scripts/pln"
sys.path.insert(0, str(PLN))

import ingest_pln_cells as I  # noqa: E402

MARKERS = ("DAPI", "CD3e")
HEADER = (["Image", "Object ID", "Object type", "Centroid X µm", "Centroid Y µm", "Parent",
           "Classification", "Nucleus: Area µm^2", "Cell: Area µm^2", "Cell: Solidity",
           "Nucleus/Cell area ratio"]
          + [f"Cell: {m}: Mean" for m in MARKERS]
          + [f"Nucleus: {m}: Mean" for m in MARKERS])

MEASURED = ["30.0", "110.0", "0.94", "0.27", "500.0", "900.0", "480.0", "870.0"]
EMPTY = ["", "", "", "", "", "", "", ""]


def row(image: str, index: int, values: list[str]) -> str:
    return "\t".join([image, f"uuid-{image[:4]}-{index:03d}", "Cell",
                      str(10.0 + index), str(20.0 + index), "Annotation (Tissue)", ""] + values)


@pytest.fixture()
def lane(tmp_path, monkeypatch):
    """Six donors in the manifest; four reach the export, one of them measuring nothing."""
    nodes = tmp_path / "pln_node_annotations.csv"
    nodes.write_text(
        "image_file,donor_id,roi,qupath_entry,node_classes,n_nodes,n_cells,excluded,"
        "recorded_by,recorded_on,excluded_reason\n"
        "6476pLN_Scan1.er.qptiff,6476,pln_1,1,pLN_top,1,3,False,t,2026-09-16,\n"
        "6579pLN_Scan1.er.qptiff,6579,pln_1,2,pLN_top,1,2,False,t,2026-09-16,\n"
        "6551pLN_Scan1.er.qptiff,6551,pln_1,3,pLN_top,1,2,False,t,2026-09-16,\n"
        "6450pLN_Scan1.er.qptiff,6450,pln_1,4,pLN_top,1,9,False,t,2026-09-16,\n"
        "6623pLN_Scan1.er.qptiff,6623,pln_1,5,pLN_top,1,1,True,t,2026-09-16,folded section\n"
        "6457pLN_Scan1.er.qptiff,6457,pln_1,6,pLN_top,1,1,False,t,2026-09-16,\n"
    )
    export = tmp_path / "IO60pLN_cells.tsv"
    lines = ["\t".join(HEADER)]
    for index in range(3):
        lines.append(row("6476pLN_Scan1.er.qptiff - resolution #1", index, MEASURED))
    for index in range(2):
        lines.append(row("6579pLN_Scan1.er.qptiff - resolution #1", index, MEASURED))
    for index in range(2):
        lines.append(row("6551pLN_Scan1.er.qptiff - resolution #1", index, EMPTY))
    lines.append(row("6623pLN_Scan1.er.qptiff - resolution #1", 0, MEASURED))
    lines.append(row("6457pLN_Scan1.er.qptiff - resolution #1", 0, MEASURED))
    export.write_text("\n".join(lines) + "\n", encoding="utf-8")

    monkeypatch.setattr(I, "NODES", nodes)
    return export, tmp_path


def run(lane, *extra: str) -> dict:
    export, tmp_path = lane
    out = tmp_path / "cells_pln"
    assert I.main(["--source", str(export), "--out", str(out),
                   "--temp-dir", str(tmp_path / "duck"), "--threads", "2",
                   "--memory-limit", "2GB", *extra]) == 0
    return json.loads((tmp_path / "cells_pln_ingest_report.json").read_text())


def ledger_of(report: dict) -> dict[str, dict]:
    return {entry["donor_id"]: entry for entry in report["per_donor"]}


def test_every_manifest_donor_gets_a_line(lane):
    """Six donors recorded, six lines -- including the two that were never written."""
    report = run(lane)
    assert set(ledger_of(report)) == {"6450", "6457", "6476", "6551", "6579", "6623"}
    assert {entry["outcome"] for entry in report["per_donor"]} == {
        "ingested", "refused", "absent"}


def test_a_donor_both_authorities_admit_is_ingested(lane):
    """6579's exclusion was lifted; a re-run writes it rather than skipping it in silence."""
    entry = ledger_of(run(lane))["6579"]
    assert entry["outcome"] == "ingested"
    assert entry["n_rows_written"] == 2


def test_a_donor_the_pln_manifest_withholds_is_named(lane):
    entry = ledger_of(run(lane))["6623"]
    assert entry["outcome"] == "refused"
    assert "pln_node_annotations" in entry["reason"]
    assert entry["n_rows_written"] is None


def test_a_centrally_excluded_donor_is_named_with_the_central_reason(lane):
    """The exclusion list is asked, never restated."""
    entry = ledger_of(run(lane))["6457"]
    assert entry["outcome"] == "refused"
    assert "DAPI" in entry["reason"]


def test_a_donor_the_export_does_not_contain_is_absent_not_silent(lane):
    """This is the shape of the 6579 failure: admitted, expected, and simply not there."""
    entry = ledger_of(run(lane))["6450"]
    assert entry["outcome"] == "absent"
    assert entry["n_rows_in_export"] is None
    assert "carries no row in this export" in entry["reason"]


def test_all_null_columns_are_measured_per_donor_on_this_export(lane):
    """Not inherited from an earlier verdict: counted here, column by column, and named."""
    report = run(lane)
    per_image = {entry["donor_id"]: entry for entry in report["per_image"]}
    empty = per_image["6551"]
    assert empty["columns_state"].startswith("unavailable")
    assert empty["n_all_null_measurement_columns"] == empty["n_measurement_columns"]
    assert set(empty["all_null_measurement_columns"]) == {
        "nucleus_area_um2", "cell_area", "cell_solidity", "nucleus_cell_area_ratio",
        "DAPI", "CD3e", "Nucleus__DAPI", "Nucleus__CD3e"}

    measured = per_image["6476"]
    assert measured["columns_state"] == "measured"
    assert measured["all_null_measurement_columns"] == []

    census = report["columns_measured_on_this_export"]
    assert census["donors_with_every_column_all_null"] == ["6551"]
    assert census["donors_with_some_columns_all_null"] == {}


def test_a_column_that_fills_on_a_re_export_stops_being_reported_as_empty(lane):
    """The point of re-measuring: the verdict follows the file, not the previous run."""
    export, tmp_path = lane
    assert ledger_of(run(lane))["6551"]["outcome"] == "ingested"
    before = {e["donor_id"]: e for e in run(lane)["per_image"]}["6551"]
    assert before["columns_state"].startswith("unavailable")

    refreshed = export.read_text().replace(
        "\t".join(["Annotation (Tissue)", ""] + EMPTY),
        "\t".join(["Annotation (Tissue)", ""] + MEASURED))
    export.write_text(refreshed, encoding="utf-8")
    after = {e["donor_id"]: e for e in run(lane)["per_image"]}["6551"]
    assert after["columns_state"] == "measured"
    assert after["all_null_measurement_columns"] == []


def test_the_export_the_run_read_is_recorded(lane):
    export, _ = lane
    report = run(lane)
    assert report["source"]["path"] == str(export)
    assert [entry["path"] for entry in report["source"]["candidates_seen"]] == [str(export)]


def test_resolve_source_takes_the_newest_and_reports_both(tmp_path, monkeypatch):
    import os
    older, newer = tmp_path / "a.tsv", tmp_path / "b.tsv"
    older.write_text("x")
    newer.write_text("x")
    os.utime(older, (1_000_000, 1_000_000))
    os.utime(newer, (2_000_000, 2_000_000))
    monkeypatch.setattr(I, "SOURCE_CANDIDATES", (older, newer))
    chosen, seen = I.resolve_source()
    assert chosen == newer
    assert [entry["path"] for entry in seen] == [str(older), str(newer)]


def test_resolve_source_honours_an_explicit_choice(tmp_path, monkeypatch):
    import os
    older, newer = tmp_path / "a.tsv", tmp_path / "b.tsv"
    older.write_text("x")
    newer.write_text("x")
    os.utime(older, (1_000_000, 1_000_000))
    os.utime(newer, (2_000_000, 2_000_000))
    monkeypatch.setattr(I, "SOURCE_CANDIDATES", (older, newer))
    chosen, seen = I.resolve_source(str(older))
    assert chosen == older
    assert len(seen) == 1


def test_resolve_source_refuses_rather_than_reading_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(I, "SOURCE_CANDIDATES", (tmp_path / "absent.tsv",))
    with pytest.raises(SystemExit, match="not found"):
        I.resolve_source()
