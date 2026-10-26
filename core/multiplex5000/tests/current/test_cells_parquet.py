"""Unit tests for the canonical cells-parquet build."""

from __future__ import annotations

import csv
from pathlib import Path

import duckdb
import pandas as pd
import pytest

from phenocycler import cells_parquet as cp

HEADER = [
    "Image", "Object ID", "Centroid X µm", "Centroid Y µm", "Parent", "Classification",
    "Cell: Area µm^2", "Nucleus: Area µm^2", "Cell: Length µm", "Cell: Circularity",
    "Cell: Solidity", "Cell: Max diameter µm", "Cell: Min diameter µm",
    "Nucleus: Length µm", "Nucleus: Circularity", "Nucleus: Solidity",
    "Nucleus: Max diameter µm", "Nucleus: Min diameter µm", "Nucleus/Cell area ratio",
    "Signed distance to annotation Islet µm", "Signed distance to annotation Tissue µm",
    "Cell: DAPI: Mean", "Cell: E-cadherin: Mean", "Nucleus: DAPI: Mean",
    *(
        f"{compartment}: {marker}: Mean (norm)"
        for marker in cp.QUPATH_NORMALIZED_MARKERS
        for compartment in cp.QUPATH_NORMALIZED_COMPARTMENTS
    ),
]


def _row(oid, area, parent="Annotation (Tissue)", image="6374_Scan1.er.qptiff", nuc=30.0):
    try:
        ratio = nuc / float(area)
    except (TypeError, ValueError, ZeroDivisionError):
        ratio = ""            # unparseable area -> unparseable ratio, as QuPath would emit
    return {
        **{c: "" for c in HEADER},
        "Image": image, "Object ID": oid, "Centroid X µm": "1.0", "Centroid Y µm": "2.0",
        "Parent": parent, "Classification": "",
        "Cell: Area µm^2": str(area), "Nucleus: Area µm^2": str(nuc),
        "Cell: Length µm": "10", "Cell: Circularity": "0.5", "Cell: Solidity": "0.9",
        "Cell: Max diameter µm": "10", "Cell: Min diameter µm": "7",
        "Nucleus: Length µm": "6", "Nucleus: Circularity": "0.8", "Nucleus: Solidity": "0.95",
        "Nucleus: Max diameter µm": "6", "Nucleus: Min diameter µm": "5",
        "Nucleus/Cell area ratio": str(ratio),
        "Signed distance to annotation Islet µm": "-5",
        "Signed distance to annotation Tissue µm": "3",
        "Cell: DAPI: Mean": "500", "Cell: E-cadherin: Mean": "700", "Nucleus: DAPI: Mean": "600",
        **{
            f"{compartment}: {marker}: Mean (norm)": str(
                -0.25 if compartment == "Cell" and marker == "Iba1" else 0.5
            )
            for marker in cp.QUPATH_NORMALIZED_MARKERS
            for compartment in cp.QUPATH_NORMALIZED_COMPARTMENTS
        },
    }


def _csv(tmp_path, rows) -> Path:
    p = tmp_path / "Cellmeasurements.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return p


def _run(tmp_path, rows, min_cell_area, image_to_donor=None):
    csv_path = _csv(tmp_path, rows)
    out = tmp_path / "out"
    (out / "cells").mkdir(parents=True, exist_ok=True)
    mapping = (
        {"6374_Scan1.er.qptiff": "6374"}
        if image_to_donor is None
        else image_to_donor
    )
    sql, mm = cp.build_sql(
        csv_path,
        out,
        None,
        min_cell_area,
        image_to_donor=mapping,
    )
    duckdb.connect().execute(sql)
    # donor_id is a hive PARTITION key (a `donor_id=<id>` directory), not a stored column, so it has
    # to come back off the path — exactly as `config.discover_donors` reads it in production.
    files = sorted((out / "cells").rglob("*.parquet"))
    parts = []
    for f in files:
        part = pd.read_parquet(f)
        part["donor_id"] = f.parent.name.split("=", 1)[1]
        parts.append(part)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return df, mm


def test_marker_means_maps_only_whole_cell_means():
    mm = cp.marker_means(HEADER)
    assert mm == {"Cell: DAPI: Mean": "DAPI", "Cell: E-cadherin: Mean": "E-cadherin"}
    assert "Nucleus: DAPI: Mean" not in mm
    assert cp.compartment_means(HEADER, "Nucleus") == {
        "Nucleus: DAPI: Mean": "DAPI"
    }


def test_no_marker_requests_a_qupath_normalized_column(tmp_path):
    """The signed ``Mean (norm)`` columns were a RESTORE-era contract for Iba1, CD68 and
    M2Gal3.  With that retired no marker selects the source, so nothing is requested -- even
    from an export that still happens to carry the columns.
    """
    assert cp.QUPATH_NORMALIZED_MARKERS == ()
    stale_header = [
        *HEADER,
        "Cell: Iba1: Mean (norm)",
        "Cytoplasm: CD68: Mean (norm)",
        "Nucleus: M2Gal3: Mean (norm)",
    ]
    assert cp.qupath_normalized_means(stale_header) == {}


def test_an_export_without_norm_columns_ingests(tmp_path):
    """This is the unblock: while Iba1/CD68/M2Gal3 were contracted to ``qupath_norm``, an
    export lacking ``Mean (norm)`` columns aborted here -- and no export we hold has them, so
    no donor could be re-ingested and the upstream was unreproducible.  It must now succeed.
    """
    path = _csv_named(tmp_path, "no_norm.csv", [_row("a", 60.0)], list(HEADER))
    assert cp.reconcile_qupath_normalized_headers([path]) == {}


def test_duckdb_overflow_goes_with_python_s_temporary_files(tmp_path, monkeypatch):
    """Not DuckDB's default, ``.tmp`` in the working directory: a stage run points Python's
    temporary files at the site's tmp_dir, and the ingest's overflow follows them."""
    import tempfile

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    con = cp._connect()
    setting = con.execute("SELECT current_setting('temp_directory')").fetchall()[0][0]
    assert Path(setting) == tmp_path / "duckdb"


def test_sanitize_matches_the_redsea_column_convention():
    assert cp.sanitize("E-cadherin") == "E_cadherin"
    assert cp.sanitize("Pan-Cytokeratin") == "Pan_Cytokeratin"
    assert cp.sanitize("Ker8/18") == "Ker8_18"


def test_min_cell_area_predicate_is_inclusive_at_the_threshold(tmp_path):
    """The filter is `>=`, so a cell exactly at the floor is retained."""
    rows = [_row("a", 19.999), _row("b", 20.0), _row("c", 20.001), _row("d", 100.0)]
    df, _ = _run(tmp_path, rows, 20.0)
    assert sorted(df.object_id) == ["b", "c", "d"]


def test_ingest_writes_no_qupath_normalized_columns(tmp_path):
    """Iba1, CD68 and M2Gal3 now take the REDSEA whole-cell mean like the other 57 markers,
    so no second scale reaches the cell table under any name.
    """
    df, _ = _run(tmp_path, [_row("a", 60.0)], 20.0)
    assert [c for c in df.columns if c.startswith("QuPathNorm__")] == []
    assert "Iba1_norm" not in df


def test_min_cell_area_zero_disables_the_filter(tmp_path):
    rows = [_row("a", 1.0), _row("b", 100.0)]
    df, _ = _run(tmp_path, rows, 0.0)
    assert sorted(df.object_id) == ["a", "b"]


def test_non_finite_area_is_always_dropped(tmp_path):
    """`isfinite` rides inside the filter, so unparseable or infinite areas never survive."""
    rows = [_row("ok", 60.0), _row("blank", ""), _row("junk", "not-a-number"), _row("inf", "inf")]
    df, _ = _run(tmp_path, rows, 20.0)
    assert sorted(df.object_id) == ["ok"]


def test_donor_id_comes_from_the_exact_manifest_mapping(tmp_path):
    rows = [
        _row("a", 60.0, image="6374_Scan1.er.qptiff"),
        _row("b", 60.0, image="HDL115pancLN_Scan1.er.qptiff"),
        _row("c", 60.0, image="6450panc_Scan1.er.qptiff"),
    ]
    df, _ = _run(
        tmp_path,
        rows,
        20.0,
        image_to_donor={
            "6374_Scan1.er.qptiff": "6374",
            "HDL115pancLN_Scan1.er.qptiff": "115",
            "6450panc_Scan1.er.qptiff": "6450",
        },
    )
    assert dict(zip(df.object_id, df.donor_id, strict=True)) == {
        "a": "6374",
        "b": "115",
        "c": "6450",
    }


def test_out_of_cohort_image_in_a_shared_export_is_ignored(tmp_path):
    rows = [
        _row("included", 60.0, image="6374_Scan1.er.qptiff"),
        _row("excluded", 60.0, image="6457_Scan1.er.qptiff"),
    ]
    df, _ = _run(tmp_path, rows, 20.0)
    assert df.object_id.tolist() == ["included"]


def test_empty_manifest_image_mapping_still_fails_closed(tmp_path):
    csv_path = _csv(tmp_path, [_row("a", 60.0)])
    with pytest.raises(Exception, match="image_to_donor cannot be empty"):
        cp.build_sql(
            csv_path,
            tmp_path / "out",
            None,
            20.0,
            image_to_donor={},
        )


def test_cell_region_scheme(tmp_path):
    rows = [
        _row("core", 60.0, parent="Islet_7"),
        _row("peri", 60.0, parent="Islet_7_exp20um"),
        _row("tissue", 60.0, parent="Annotation (Tissue)"),
        _row("root", 60.0, parent="Root object (Image)"),
        _row("other", 60.0, parent="Something Else"),
    ]
    df, _ = _run(tmp_path, rows, 20.0)
    got = dict(zip(df.object_id, df.cell_region, strict=True))
    assert got == {"core": "core", "peri": "peri", "tissue": "tissue",
                   "root": "tissue", "other": "other"}
    assert dict(zip(df.object_id, df.islet_num, strict=True))["core"] == "7"


def test_morphology_columns_needed_by_geometry_qc_are_all_present(tmp_path):
    """Geometry QC reads these from canonical cells; none may go missing."""
    df, _ = _run(tmp_path, [_row("a", 60.0)], 20.0)
    from phenocycler.geometry_stage import GEOMETRY_SOURCE_COLUMNS

    missing = set(GEOMETRY_SOURCE_COLUMNS) - set(df.columns)
    assert not missing, f"cells parquet is missing QC inputs: {sorted(missing)}"


def test_build_sql_rejects_a_csv_with_no_marker_columns(tmp_path):
    p = tmp_path / "bad.csv"
    with open(p, "w", newline="") as f:
        csv.writer(f).writerow(cp.REQUIRED_CONTEXT_COLUMNS)
    with pytest.raises(SystemExit, match="No 'Cell: <marker>: Mean' columns"):
        cp.build_sql(
            p,
            tmp_path / "out",
            None,
            20.0,
            image_to_donor={"image": "D"},
        )


# --------------------------------------------------------------- multi-export ingest (2026-07-29)


def _csv_named(tmp_path, name, rows, header):
    """One QuPath export with an explicit header (so batches can differ)."""
    p = tmp_path / name
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return p


def _two_batch_exports(tmp_path):
    """Batch 1 carries CD38, batch 2 carries b-Catenin1 — the real panel difference in this cohort."""
    h1 = HEADER + ["Cell: CD38: Mean"]
    h2 = HEADER + ["Cell: b-Catenin1: Mean"]
    r1 = {**_row("a", 50.0, image="6374_Scan1.er.qptiff"), "Cell: CD38: Mean": "111"}
    r2 = {**_row("b", 50.0, image="6523_Scan1.er.qptiff"), "Cell: b-Catenin1: Mean": "222"}
    return (_csv_named(tmp_path, "batch1.csv", [r1], h1),
            _csv_named(tmp_path, "batch2.csv", [r2], h2))


def test_reconcile_headers_keeps_union_and_records_availability(tmp_path):
    """Batch-specific markers remain usable for donors that acquired them."""
    p1, p2 = _two_batch_exports(tmp_path)
    union, report = cp.reconcile_headers([p1, p2])

    assert set(union.values()) == {"DAPI", "E_cadherin", "CD38", "b_Catenin1"}
    assert report["union"] == 4
    assert report["intersection"] == 2
    assert report["files"] == {str(p1): 3, str(p2): 3}
    assert report["availability"][str(p1)]["CD38"]
    assert not report["availability"][str(p1)]["b_Catenin1"]
    assert report["availability"][str(p2)]["b_Catenin1"]
    assert not report["availability"][str(p2)]["CD38"]


def test_reconcile_headers_fails_loudly_on_a_missing_export(tmp_path):
    p1, _ = _two_batch_exports(tmp_path)
    with pytest.raises(SystemExit, match="not found"):
        cp.reconcile_headers([p1, tmp_path / "absent.csv"])
    with pytest.raises(SystemExit, match="cohort manifest.*no QuPath measurement CSV"):
        cp.reconcile_headers([])


def test_build_sql_reads_every_export_in_one_pass(tmp_path):
    """Both batches must reach the output, aligned by NAME (their column counts differ)."""
    p1, p2 = _two_batch_exports(tmp_path)
    out = tmp_path / "out"
    (out / "cells").mkdir(parents=True, exist_ok=True)
    mapping = {
        "6374_Scan1.er.qptiff": "6374",
        "6523_Scan1.er.qptiff": "6523",
    }
    sql, mm = cp.build_sql(
        [p1, p2], out, None, 20.0, image_to_donor=mapping
    )
    assert "UNION ALL BY NAME" in sql
    duckdb.connect().execute(sql)

    donors = sorted(d.name.split("=", 1)[1] for d in (out / "cells").glob("donor_id=*"))
    assert donors == ["6374", "6523"]                 # a donor from each batch
    frames = [pd.read_parquet(f) for f in sorted((out / "cells").rglob("*.parquet"))]
    for frame in frames:
        assert {"DAPI", "E_cadherin", "CD38", "b_Catenin1"} <= set(frame.columns)
        assert "Nucleus__DAPI" in frame
    by_object = pd.concat(frames, ignore_index=True).set_index("object_id")
    assert by_object.loc["a", "CD38"] == 111
    assert pd.isna(by_object.loc["a", "b_Catenin1"])
    assert by_object.loc["b", "b_Catenin1"] == 222
    assert pd.isna(by_object.loc["b", "CD38"])


def test_build_sql_still_accepts_a_single_path(tmp_path):
    """Back-compat: one export is the degenerate case of the list, not a separate code path."""
    p1, _ = _two_batch_exports(tmp_path)
    mapping = {"6374_Scan1.er.qptiff": "6374"}
    sql_single, mm_single = cp.build_sql(
        p1, tmp_path / "out", None, 20.0, image_to_donor=mapping
    )
    sql_list, mm_list = cp.build_sql(
        [p1], tmp_path / "out", None, 20.0, image_to_donor=mapping
    )
    assert sql_single == sql_list and mm_single == mm_list
    assert "CD38" in mm_single.values()               # nothing dropped when there is nothing to reconcile


def test_annotation_coverage_names_a_donor_whose_ring_was_never_drawn(tmp_path):
    """The check that would have caught donor 6539.

    `cell_region` comes from QuPath's `Parent`, so a cell lands in `peri` only where an
    `Islet_N_exp20um` annotation exists. 6539 reached the cohort with 541 islets and 124 peri
    cells -- the expansion had run on 3 of 539 -- and nothing flagged it, so it read
    downstream as a low-immune ND donor in the one compartment the study is about.
    """
    import duckdb
    import pandas as pd
    from phenocycler.cells_parquet import annotation_coverage, report_annotation_coverage

    def _write(donor, rows):
        d = tmp_path / f"donor_id={donor}"
        d.mkdir(parents=True)
        pd.DataFrame(rows).to_parquet(d / "data.parquet", index=False)

    # good donor: every islet with core cells also has ring cells, and a tissue boundary
    _write("good", {
        "parent_raw": ["Islet_1", "Islet_1_exp20um", "Islet_2", "Islet_2_exp20um",
                       "Annotation (Tissue)"],
        "cell_region": ["core", "peri", "core", "peri", "tissue"],
        "islet_num": ["1", "1", "2", "2", None],
    })
    # broken donor: two islets, one stray ring cell, and no tissue annotation at all
    _write("broken", {
        "parent_raw": ["Islet_1", "Islet_2", "Islet_1_exp20um", "Root object (Image)"],
        "cell_region": ["core", "core", "peri", "tissue"],
        "islet_num": ["1", "2", "1", None],
    })

    con = duckdb.connect()
    cov = annotation_coverage(con, tmp_path)

    assert cov["good"]["ring_coverage"] == 1.0
    assert cov["good"]["has_tissue_annotation"] is True
    assert cov["broken"]["ring_coverage"] == 0.5
    assert cov["broken"]["has_tissue_annotation"] is False

    warnings = report_annotation_coverage(cov)
    assert len(warnings) == 1
    assert "donor broken" in warnings[0]
    assert "1/2 islets" in warnings[0]
    assert "no Annotation (Tissue)" in warnings[0]


def test_an_image_in_two_exports_is_read_from_only_its_contracted_one():
    """The guard against silently doubling a donor.

    Donor 6539 was re-exported after its QuPath annotation was corrected, but the shared
    batch-1 `Cellmeasurements.csv` still contains the old rows for that image. A filter
    applied globally over the union admits both copies -- the corrected and the broken export
    interleaved, one donor at double size, with nothing raised. The manifest already records
    which export each image belongs to; the filter now honours it per file.
    """
    shared = "/exports/Cellmeasurements.csv"
    fixed = "/exports/Cellmeasurements_6539parentFixed.csv"
    image_to_donor = {"6374_Scan1": "6374", "6539_Scan1": "6539"}
    image_to_csv = {"6374_Scan1": shared, "6539_Scan1": fixed}

    src = cp.contracted_image_source([Path(shared), Path(fixed)], image_to_donor, image_to_csv)

    # Each read_csv carries its own image filter. Inspect the per-file clauses only -- the
    # trailing global filter still names every contracted image, deliberately, as a second
    # line of defence.
    clauses = src.split("read_csv(")[1:]
    assert len(clauses) == 2
    shared_clause = next(c for c in clauses if "Cellmeasurements.csv" in c[:60])
    fixed_clause = next(c for c in clauses if "6539parentFixed" in c[:60])
    # each clause is truncated at the start of the next one, so the filter it carries is its own
    assert "'6374_Scan1'" in shared_clause.split("AS all_qupath_exports")[0]
    assert "'6539_Scan1'" not in shared_clause.split("UNION ALL")[0]
    assert "'6539_Scan1'" in fixed_clause.split("AS all_qupath_exports")[0]

    # Without the assignment no per-file filter is emitted at all -- exactly the behaviour
    # that would have admitted 6539 from both exports.
    legacy = cp.contracted_image_source([Path(shared), Path(fixed)], image_to_donor)
    first = legacy.split("read_csv(")[1].split("UNION ALL")[0]
    assert "WHERE" not in first


def test_an_incomplete_image_to_csv_assignment_is_refused():
    """Fail closed: a partial assignment would read some export in full."""
    image_to_donor = {"a": "1", "b": "2"}
    with pytest.raises(Exception, match="does not assign an export to every image"):
        cp.contracted_image_source([Path("/x.csv")], image_to_donor, {"a": "/x.csv"})

    with pytest.raises(Exception, match="absent from image_to_donor"):
        cp.contracted_image_source(
            [Path("/x.csv")], image_to_donor,
            {"a": "/x.csv", "b": "/x.csv", "ghost": "/x.csv"})

    with pytest.raises(Exception, match="no contracted image is assigned to export"):
        cp.contracted_image_source(
            [Path("/x.csv"), Path("/unused.csv")], image_to_donor,
            {"a": "/x.csv", "b": "/x.csv"})
