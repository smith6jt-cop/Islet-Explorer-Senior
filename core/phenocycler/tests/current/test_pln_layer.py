"""pLN node-membership layer: the totals reconcile, absence is named, the UUID survives.

The fixture is a six-cell synthetic slide built in tmp_path -- roster, annotation geojson, pixel
calibration, the derivation's own output and QuPath's annotation audit log.  Nothing here reads
the 8.9M-row artifact; the one test that does is skipped when it is not on disk.

Run with
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
      core/phenocycler/tests/current/test_pln_layer.py -q
"""

from __future__ import annotations

import glob
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
SCRIPT = REPO / "core/phenocycler/scripts/pln/build_pln_layer.py"
SPEC = importlib.util.spec_from_file_location("build_pln_layer", SCRIPT)
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)

ROSTER = REPO / "core/phenocycler/data/integration/pln_node_annotations.csv"
#: The published root, resolved the way the builder and every R consumer resolve it.  It is NOT
#: the in-repo data/islet_layers: an artifact written there is invisible to the cards.
LIVE = builder.DEFAULT_OUT

#: QuPath UUIDs, verbatim.  The point of the artifact is that these survive the whole chain, so
#: the fixture uses real-shaped ones and the test compares them character for character.
CELL_IDS = [
    "3f0a1c2e-1111-4a00-8000-000000000001",
    "3f0a1c2e-2222-4a00-8000-000000000002",
    "3f0a1c2e-3333-4a00-8000-000000000003",
    "3f0a1c2e-4444-4a00-8000-000000000004",
    "3f0a1c2e-5555-4a00-8000-000000000005",
    "3f0a1c2e-6666-4a00-8000-000000000006",
]
TOP, BOTTOM, DEAD = "a-top", "a-bottom", "a-dead"
TISSUE_TOP, TISSUE_BOTTOM = "a-tissue-top", "a-tissue-bottom"

#: donor 7001 is measured; 7002 is eligible but its derivation never ran; 7003 the tracked
#: manifest withholds.  None of the three is in phenocycler.cohort's exclusions, so the cohort
#: authority admits all three and the manifest authority is what withholds 7003.
ROSTER_HEADER = ("image_file,donor_id,roi,qupath_entry,node_classes,n_nodes,n_cells,excluded,"
                 "recorded_by,recorded_on,excluded_reason")
ROSTER_CSV = f"""# synthetic
{ROSTER_HEADER}
A.qptiff,7001,pln_1,101,pLN_top;pLN_bottom;pLN_dead,3,6,False,test,2026-09-17,
B.qptiff,7002,pln_1,102,pLN_only,1,5,False,test,2026-09-17,
C.qptiff,7003,pln_1,103,pLN_top,1,3,True,test,2026-09-17,synthetic: withheld by this fixture
"""

#: One image whose single node holds every cell on it and outlines more than the slide's Tissue.
#: That is donor 6534's situation in the real cohort, where `pLN_bottom` is the slide wearing a
#: node's name -- and 115's and 6380's too.
SLIDE_ROSTER_CSV = f"""# synthetic
{ROSTER_HEADER}
D.qptiff,7004,pln_1,104,pLN_only,1,3,False,test,2026-09-17,
"""
SLIDE_CELL_IDS = [
    "7f0a1c2e-1111-4a00-8000-00000000000a",
    "7f0a1c2e-2222-4a00-8000-00000000000b",
    "7f0a1c2e-3333-4a00-8000-00000000000c",
]


def _polygon(x0, y0, x1, y1):
    return {"type": "Polygon", "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]]}


def _feature(fid, name, geometry):
    return {"type": "Feature", "id": fid,
            "properties": {"classification": {"name": name}}, "geometry": geometry}


@pytest.fixture
def world(tmp_path):
    """A synthetic slide: two area nodes, one zero-area node, six cells, one of them outside."""
    roster = tmp_path / "pln_node_annotations.csv"
    roster.write_text(ROSTER_CSV)

    geojson_dir = tmp_path / "node_geojson"
    geojson_dir.mkdir()
    (geojson_dir / "entry101.geojson").write_text(json.dumps({
        "type": "FeatureCollection",
        "features": [
            _feature(TOP, "pLN_top", _polygon(0, 0, 10, 10)),
            _feature(BOTTOM, "pLN_bottom", _polygon(0, 20, 10, 30)),
            # a line bounds no region; this is donor 6534's situation in the real cohort
            _feature(DEAD, "pLN_dead",
                     {"type": "LineString", "coordinates": [[0, 15], [10, 15]]}),
            _feature(TISSUE_TOP, "Tissue", _polygon(1, 1, 9, 9)),
            _feature(TISSUE_BOTTOM, "Tissue", _polygon(1, 21, 9, 29)),
        ],
    }))

    project = tmp_path / "project"
    (project / "data/101").mkdir(parents=True)
    (project / "data/101/server.json").write_text(json.dumps(
        {"metadata": {"pixelCalibration": {"pixelWidth": {"value": 1.0},
                                           "pixelHeight": {"value": 1.0}}}}))

    derived = tmp_path / "derived"
    (derived / "node_membership").mkdir(parents=True)
    pq.write_table(
        pa.table({
            "object_id": pa.array(CELL_IDS, type=pa.string()),
            "donor_id": pa.array(["7001"] * 6, type=pa.string()),
            "image_file": pa.array(["A.qptiff"] * 6, type=pa.string()),
            "node_id": pa.array(
                ["pLN_top", "pLN_top", "pLN_top", "pLN_bottom", "pLN_bottom", "unassigned"],
                type=pa.string()),
            "x_um": pa.array([2.0, 3.0, 4.0, 2.0, 3.0, 50.0], type=pa.float64()),
            "y_um": pa.array([2.0, 3.0, 4.0, 22.0, 23.0, 50.0], type=pa.float64()),
            "n_node_regions_containing": pa.array([1, 1, 1, 1, 1, 0], type=pa.int64()),
            "tissue_annotation_index": pa.array([0, 0, 0, 1, 1, -1], type=pa.int64()),
        }),
        derived / "node_membership/pln_node_membership_7001.parquet",
    )
    pd.DataFrame([
        {"donor_id": "7001", "qupath_entry": 101, "image_file": "A.qptiff",
         "annotation_id": TOP, "classification": "pLN_top", "geom_type": "Polygon",
         "area_px": 100.0, "area_um2": 100.0, "is_axis_aligned_rectangle": True,
         "n_vertices": 5.0, "role": "node"},
        {"donor_id": "7001", "qupath_entry": 101, "image_file": "A.qptiff",
         "annotation_id": BOTTOM, "classification": "pLN_bottom", "geom_type": "Polygon",
         "area_px": 100.0, "area_um2": 100.0, "is_axis_aligned_rectangle": True,
         "n_vertices": 5.0, "role": "node"},
        {"donor_id": "7001", "qupath_entry": 101, "image_file": "A.qptiff",
         "annotation_id": DEAD, "classification": "pLN_dead", "geom_type": "LineString",
         "area_px": 0.0, "area_um2": 0.0, "is_axis_aligned_rectangle": False,
         "n_vertices": None, "role": "refused_no_area"},
        {"donor_id": "7001", "qupath_entry": 101, "image_file": "A.qptiff",
         "annotation_id": TISSUE_TOP, "classification": "Tissue", "geom_type": "Polygon",
         "area_px": 64.0, "area_um2": 64.0, "is_axis_aligned_rectangle": True,
         "n_vertices": 5.0, "role": "not_a_node"},
        {"donor_id": "7001", "qupath_entry": 101, "image_file": "A.qptiff",
         "annotation_id": TISSUE_BOTTOM, "classification": "Tissue", "geom_type": "Polygon",
         "area_px": 64.0, "area_um2": 64.0, "is_axis_aligned_rectangle": True,
         "n_vertices": 5.0, "role": "not_a_node"},
    ]).to_csv(derived / "node_geometry_audit.csv", index=False)

    # QuPath's own child counts: the top Tissue agrees, the bottom one is one cell short of
    # agreeing.  A cross-check that disagrees is reported, never adopted.
    audit_log = tmp_path / "qupath_annotation_audit.log"
    audit_log.write_text(
        "ENTRY=101\n"
        "READ_MS=12 N_OBJECTS=6 N_ANNOTATIONS=3\n"
        f"ANNOT\tid={TOP}\tclass=pLN_top\troi=AREA/RectangleROI\tnPoints=4\tareaPx=100.0\t"
        "nChildren=0\n"
        f"ANNOT\tid={TISSUE_TOP}\tclass=Tissue\troi=AREA/PolygonROI\tnPoints=5\tareaPx=64.0\t"
        "nChildren=3\n"
        f"ANNOT\tid={TISSUE_BOTTOM}\tclass=Tissue\troi=AREA/PolygonROI\tnPoints=5\tareaPx=64.0\t"
        "nChildren=3\n"
    )

    out = tmp_path / "layer"
    assert builder.main([
        "--roster", str(roster), "--derived", str(derived), "--out", str(out),
        "--geojson-dir", str(geojson_dir), "--qupath-project", str(project),
        "--audit-log", str(audit_log),
        "--ingest-report", str(tmp_path / "no_such_ingest_report.json"),
    ]) == 0
    manifest = json.loads((out / "MANIFEST.json").read_text())
    return {"out": out, "manifest": manifest, "roster": roster}


def _donor(manifest, donor_id):
    rows = [d for d in manifest["donors"] if d["donor_id"] == donor_id]
    assert len(rows) == 1, f"{donor_id}: {len(rows)} manifest rows"
    return rows[0]


def test_per_donor_total_equals_the_roster_total(world):
    """assigned + unassigned == the hand-recorded n_cells, per donor and over the artifact."""
    manifest = world["manifest"]
    roster = pd.read_csv(world["roster"], comment="#", dtype={"donor_id": str})

    measured = [d for d in manifest["donors"] if d["state"] == "measured"]
    assert [d["donor_id"] for d in measured] == ["7001"]
    for donor in measured:
        recorded = int(roster.loc[roster["donor_id"] == donor["donor_id"], "n_cells"].iloc[0])
        assert donor["cells_assigned"] + donor["cells_unassigned"] == recorded
        assert donor["cells_total"] == recorded
        assert donor["reconciles_with_roster"] is True

    written = sum(
        pq.ParquetFile(f).metadata.num_rows
        for f in glob.glob(str(world["out"] / "donor_id=*/*.parquet"))
    )
    assert written == manifest["totals"]["cells_written"] == sum(
        d["cells_total"] for d in measured)
    assert manifest["totals"]["reconciles_with_roster"] is True


def test_a_cell_in_no_node_is_kept_with_a_reason(world):
    """Not dropped, not labelled with a node it is not in, not called zero."""
    table = pq.read_table(world["out"] / "donor_id=7001/pln_node_membership.parquet")
    rows = table.to_pylist()
    outside = [r for r in rows if r["node_state"] == "unassigned"]
    assert len(outside) == 1
    assert outside[0]["node_name"] is None
    assert outside[0]["node_index"] is None
    assert outside[0]["reason"] and "outside" in outside[0]["reason"]
    assert outside[0]["object_id"] == CELL_IDS[5]
    assert all(r["reason"] == "" for r in rows if r["node_state"] == "measured")


def test_a_donor_with_no_eligible_export_is_unavailable_with_a_reason(world):
    """7002 is eligible and has no derivation; 7003 the tracked manifest withholds.

    Neither may appear as a donor with 0 cells: a 0 would read as a measurement.
    """
    manifest = world["manifest"]
    nodes = pq.read_table(world["out"] / "nodes.parquet").to_pylist()

    for donor_id, state in (("7002", "unavailable"), ("7003", "withheld")):
        donor = _donor(manifest, donor_id)
        assert donor["state"] == state
        assert donor["reason"], f"{donor_id}: state {state} with no reason"
        for field in ("cells_total", "cells_assigned", "cells_unassigned"):
            assert donor[field] is None, f"{donor_id}.{field} is {donor[field]}, not None"
        assert not glob.glob(str(world["out"] / f"donor_id={donor_id}/*"))
        own = [n for n in nodes if n["donor_id"] == donor_id]
        assert own and all(n["state"] == state and n["n_cells"] is None for n in own)

    not_written = {d["donor_id"]: d["state"] for d in manifest["totals"]["donors_not_written"]}
    assert not_written == {"7002": "unavailable", "7003": "withheld"}


def test_a_zero_area_node_is_unavailable_not_a_zero_cell_node(world):
    """donor 6534's situation: a PolylineROI bounds no region, so nothing can be tested."""
    nodes = pq.read_table(world["out"] / "nodes.parquet").to_pylist()
    dead = [n for n in nodes if n["donor_id"] == "7001" and n["node_name"] == "pLN_dead"]
    assert len(dead) == 1
    assert dead[0]["state"] == "unavailable"
    assert dead[0]["n_cells"] is None
    assert "zero_area" in dead[0]["reason"]
    # it keeps its roster position rather than being dropped from the node list
    assert dead[0]["node_index"] == 2
    assert [n["node_index"] for n in nodes if n["donor_id"] == "7001"] == [0, 1, 2]


def test_object_id_is_a_string_verbatim_and_unique_within_the_donor(world):
    table = pq.read_table(world["out"] / "donor_id=7001/pln_node_membership.parquet")
    assert table.schema.field("object_id").type == pa.string()
    ids = table.column("object_id").to_pylist()
    assert ids == CELL_IDS                      # verbatim, in the derivation's own order
    assert len(set(ids)) == len(ids) == table.num_rows
    assert table.schema.field("donor_id").type == pa.string()
    assert set(table.column("donor_id").to_pylist()) == {"7001"}


def test_the_qupath_cross_check_reports_a_disagreement_rather_than_adopting_it(world):
    check = world["manifest"]["qupath_cross_check"]
    assert check["nodes_checked"] == 2
    assert check["nodes_reproducing_exactly"] == 1
    assert [(d["node_name"], d["derived_n_cells"], d["qupath_tissue_child_count"], d["delta"])
            for d in check["nodes_differing"]] == [("pLN_bottom", 2, 3, -1)]


def test_an_absent_audit_log_leaves_the_check_unavailable_never_zero(world, tmp_path):
    """No log is no check.  It is never a count of 0 and never a silent 'reproduced'."""
    out = tmp_path / "layer_no_log"
    assert builder.main([
        "--roster", str(world["roster"]),
        "--derived", str(tmp_path / "derived"), "--out", str(out),
        "--geojson-dir", str(tmp_path / "node_geojson"),
        "--qupath-project", str(tmp_path / "project"),
        "--audit-log", str(tmp_path / "absent.log"),
        "--ingest-report", str(tmp_path / "no_such_ingest_report.json"),
    ]) == 0
    manifest = json.loads((out / "MANIFEST.json").read_text())
    assert manifest["qupath_cross_check"]["nodes_checked"] == 0
    nodes = pq.read_table(out / "nodes.parquet").to_pylist()
    for node in nodes:
        assert node["qupath_tissue_child_count"] is None
        assert node["delta_vs_qupath"] is None
        if node["state"] == "measured":
            assert node["qupath_check"].startswith("unavailable(")
    # the membership itself is unaffected: it never depended on the cross-check
    assert manifest["totals"]["cells_written"] == 6


@pytest.mark.skipif(not (LIVE / "MANIFEST.json").exists(),
                    reason="the cohort artifact is not built in this checkout")
def test_the_built_artifact_reconciles_with_the_tracked_roster():
    """Arithmetic only: the artifact's own total against the hand-recorded roster it came from."""
    manifest = json.loads((LIVE / "MANIFEST.json").read_text())
    roster = pd.read_csv(ROSTER, comment="#", dtype={"donor_id": str})
    roster["donor_id"] = roster["donor_id"].map(builder.norm_donor)
    written = {d["donor_id"] for d in manifest["donors"] if d["state"] == "measured"}
    expected = int(roster.loc[roster["donor_id"].isin(written), "n_cells"].sum())

    assert manifest["totals"]["cells_written"] == expected
    assert (manifest["totals"]["cells_assigned_to_a_node"]
            + manifest["totals"]["cells_unassigned"]) == expected
    on_disk = sum(pq.ParquetFile(f).metadata.num_rows
                  for f in glob.glob(str(LIVE / "donor_id=*/*.parquet")))
    assert on_disk == expected


# --------------------------------------------------------------------------------------------
# The six defects an adversarial pass proved on 2026-09-17.  Each test fails on the old builder.
# --------------------------------------------------------------------------------------------


def _rebuild(world, roster_csv: str, tmp_path, name: str) -> dict:
    """Build again into the SAME root, from a roster this call writes."""
    roster = tmp_path / name
    roster.write_text(roster_csv)
    assert builder.main([
        "--roster", str(roster),
        "--derived", str(tmp_path / "derived"), "--out", str(world["out"]),
        "--geojson-dir", str(tmp_path / "node_geojson"),
        "--qupath-project", str(tmp_path / "project"),
        "--audit-log", str(tmp_path / "qupath_annotation_audit.log"),
        "--ingest-report", str(tmp_path / "no_such_ingest_report.json"),
    ]) == 0
    return json.loads((world["out"] / "MANIFEST.json").read_text())


def test_a_rebuild_into_a_populated_root_removes_the_partition_it_no_longer_writes(
        world, tmp_path):
    """Defect 1.  7001 was written; withhold it and rebuild into the same root.

    The old builder left 7001's partition behind and reported the smaller total anyway, because
    the total came from its own counter rather than from the files.
    """
    assert (world["out"] / "donor_id=7001" / builder.PARTITION_FILE).exists()

    withheld = ROSTER_CSV.replace(
        "A.qptiff,7001,pln_1,101,pLN_top;pLN_bottom;pLN_dead,3,6,False,test,2026-09-17,",
        "A.qptiff,7001,pln_1,101,pLN_top;pLN_bottom;pLN_dead,3,6,True,test,2026-09-17,"
        "synthetic: withheld to prove the root reconciles")
    manifest = _rebuild(world, withheld, tmp_path, "roster_no7001.csv")

    assert not (world["out"] / "donor_id=7001").exists()
    recon = manifest["output_root_reconciliation"]
    assert recon["donors_written_this_run"] == []
    assert recon["donors_on_disk_after"] == []
    assert recon["root_matches_this_run"] is True
    assert [(r["donor_id"], r["rows_removed"]) for r in recon["stale_partitions_removed"]] \
        == [("7001", 6)]
    assert manifest["totals"]["cells_written"] == 0


def test_a_donor_gained_or_lost_is_named_rather_than_appearing_silently(world, tmp_path):
    """Defect 5's other half.  6579's 1,231,818 cells joined a built layer; that must be stated.

    Before and after are both measured off the files, so a rebuild cannot change the cohort
    without saying which donors moved and by how many cells.
    """
    first = world["manifest"]["output_root_reconciliation"]
    assert first["before_this_run"]["donors"] == []      # the fixture root started empty
    assert first["donors_gained"] == ["7001"]
    assert first["cells_delta"] == 6

    withheld = ROSTER_CSV.replace(
        "A.qptiff,7001,pln_1,101,pLN_top;pLN_bottom;pLN_dead,3,6,False,test,2026-09-17,",
        "A.qptiff,7001,pln_1,101,pLN_top;pLN_bottom;pLN_dead,3,6,True,test,2026-09-17,"
        "synthetic: withheld to prove the loss is named")
    again = _rebuild(world, withheld, tmp_path, "roster_lost7001.csv")
    recon = again["output_root_reconciliation"]
    assert recon["before_this_run"]["donors"] == ["7001"]
    assert recon["before_this_run"]["cells"] == 6
    assert recon["before_this_run"]["previous_manifest_totals"]["cells_written"] == 6
    assert recon["donors_lost"] == ["7001"]
    assert recon["donors_gained"] == []
    assert recon["cells_delta"] == -6


def test_the_manifest_total_is_measured_off_the_files_not_off_a_counter(world):
    """Defect 1.  The total has to be a measurement OF THE FILES, and say so."""
    manifest, totals = world["manifest"], world["manifest"]["totals"]
    on_disk = sum(pq.ParquetFile(f).metadata.num_rows
                  for f in glob.glob(str(world["out"] / "donor_id=*/*.parquet")))
    assert totals["cells_written"] == on_disk
    assert "metadata.num_rows" in totals["cells_written_measured_by"]
    assert totals["cells_counted_in_memory_by_this_run"] == on_disk
    assert totals["files_match_the_in_memory_counter"] is True
    assert [p["rows"] for p in manifest["output_root_reconciliation"]["partitions"]] == [on_disk]


def test_how_to_read_is_probed_on_this_build_not_transcribed(world):
    """Defect 2.  The route strings carried a cohort row count that a rebuild re-emitted."""
    routes = world["manifest"]["how_to_read"]
    ok = routes["ds.dataset(glob.glob(root + '/donor_id=*/*.parquet'))"]
    assert ok.startswith("OK, 6 rows, donor_id string")
    assert "8,961,168" not in json.dumps(routes)
    assert routes["measured_files"] == 1
    assert routes["ds.dataset(root)"].startswith("FAILS ")


def test_the_unassigned_note_is_derived_from_this_builds_own_cells(world):
    """Defect 2.  '~1.3 um' was baked in and survived a rebuild over other cells."""
    check = world["manifest"]["qupath_cross_check"]
    note = check["unassigned_cells_note"]
    assert "1.3" not in note
    assert note.startswith("1 cell(s) fell in no node on this build")
    measured = check["unassigned_cells"][0]["distance_outside_um"]
    assert f"{measured:g}" in note


def test_the_annotation_count_is_counted_not_transcribed(world):
    """Defect 3.  84 ANNOT records plus 18 READ_MS lines were reported as 102 annotations."""
    check = world["manifest"]["qupath_cross_check"]
    log = check["audit_log_measured"]
    assert log["state"] == "measured"
    assert log["annotation_records"] == 3            # the fixture log holds three ANNOT lines
    assert log["read_ms_lines"] == 1
    assert log["annot_plus_read_ms_lines"] == 4      # a LINE count, and named as one
    assert "a line count" in log["not_the_annotation_count"]
    assert log["exporter_declared_total"] == 3
    assert log["agrees_with_exporter"] is True
    # the geojson is the second count, reported beside it rather than overriding it
    assert check["geojson_annotations_measured"]["features_total"] == 5
    assert check["audit_log_agrees_with_geojson"] is False


def test_an_excluded_row_must_say_why_and_may_not_restate_a_lapsed_cohort_ruling(tmp_path):
    """Defect 5.  6579 kept excluded=True after phenocycler.cohort reinstated it."""
    from phenocycler.cohort import REINSTATED

    reinstated = sorted(REINSTATED)
    assert reinstated, "this test needs at least one reinstated donor to be meaningful"
    donor = reinstated[0]

    no_reason = pd.DataFrame([{
        "donor_id": "7003", "excluded": True, "excluded_reason": "",
        "image_file": "C.qptiff", "roi": "pln_1", "qupath_entry": 103,
        "node_classes": ["pLN_top"], "n_nodes": 1, "n_cells": 3,
    }])
    with pytest.raises(SystemExit) as err:
        builder.check_roster_exclusions(no_reason, tmp_path / "roster.csv")
    assert "empty excluded_reason" in str(err.value)

    lapsed = no_reason.assign(donor_id=donor, excluded_reason="centrally excluded")
    with pytest.raises(SystemExit) as err:
        builder.check_roster_exclusions(lapsed, tmp_path / "roster.csv")
    assert "REINSTATED" in str(err.value) and "lapsed" in str(err.value)


def test_the_tracked_roster_does_not_withhold_a_reinstated_donor():
    """Defect 5, on the tracked file itself rather than on a fixture."""
    builder.check_roster_exclusions(builder.read_roster(ROSTER), ROSTER)


@pytest.fixture
def slide_world(tmp_path):
    """One image, one node, and that node holds the slide.  Donor 6534's situation."""
    roster = tmp_path / "slide_roster.csv"
    roster.write_text(SLIDE_ROSTER_CSV)

    geojson_dir = tmp_path / "slide_geojson"
    geojson_dir.mkdir()
    (geojson_dir / "entry104.geojson").write_text(json.dumps({
        "type": "FeatureCollection",
        "features": [
            _feature("d-node", "pLN_only", _polygon(0, 0, 10, 10)),      # area 100
            _feature("d-tissue", "Tissue", _polygon(1, 1, 9, 9)),        # area 64
        ],
    }))

    project = tmp_path / "slide_project"
    (project / "data/104").mkdir(parents=True)
    (project / "data/104/server.json").write_text(json.dumps(
        {"metadata": {"pixelCalibration": {"pixelWidth": {"value": 1.0},
                                           "pixelHeight": {"value": 1.0}}}}))

    derived = tmp_path / "slide_derived"
    (derived / "node_membership").mkdir(parents=True)
    pq.write_table(
        pa.table({
            "object_id": pa.array(SLIDE_CELL_IDS, type=pa.string()),
            "donor_id": pa.array(["7004"] * 3, type=pa.string()),
            "image_file": pa.array(["D.qptiff"] * 3, type=pa.string()),
            "node_id": pa.array(["pLN_only"] * 3, type=pa.string()),
            "x_um": pa.array([2.0, 3.0, 4.0], type=pa.float64()),
            "y_um": pa.array([2.0, 3.0, 4.0], type=pa.float64()),
            "n_node_regions_containing": pa.array([1, 1, 1], type=pa.int64()),
            "tissue_annotation_index": pa.array([0, 0, 0], type=pa.int64()),
        }),
        derived / "node_membership/pln_node_membership_7004.parquet",
    )
    pd.DataFrame([
        {"donor_id": "7004", "qupath_entry": 104, "image_file": "D.qptiff",
         "annotation_id": "d-node", "classification": "pLN_only", "geom_type": "Polygon",
         "area_px": 100.0, "area_um2": 100.0, "is_axis_aligned_rectangle": True,
         "n_vertices": 5.0, "role": "node"},
        {"donor_id": "7004", "qupath_entry": 104, "image_file": "D.qptiff",
         "annotation_id": "d-tissue", "classification": "Tissue", "geom_type": "Polygon",
         "area_px": 64.0, "area_um2": 64.0, "is_axis_aligned_rectangle": True,
         "n_vertices": 5.0, "role": "not_a_node"},
    ]).to_csv(derived / "node_geometry_audit.csv", index=False)

    out = tmp_path / "slide_layer"
    assert builder.main([
        "--roster", str(roster), "--derived", str(derived), "--out", str(out),
        "--geojson-dir", str(geojson_dir), "--qupath-project", str(project),
        "--audit-log", str(tmp_path / "absent.log"),
        "--ingest-report", str(tmp_path / "no_such_ingest_report.json"),
    ]) == 0
    return {"out": out, "manifest": json.loads((out / "MANIFEST.json").read_text())}


def test_a_node_that_is_the_slide_is_named_as_one_in_the_node_table(slide_world):
    """Defect 6.  A per-node number for such a node is a per-slide number."""
    nodes = pq.read_table(slide_world["out"] / "nodes.parquet").to_pylist()
    assert len(nodes) == 1
    node = nodes[0]
    assert node["node_scope"] == "slide"
    assert node["n_usable_nodes_on_image"] == 1
    assert node["holds_every_assigned_cell_on_image"] is True
    assert node["image_tissue_area_um2"] == 64.0
    assert node["area_ratio_node_over_image_tissue"] == pytest.approx(100.0 / 64.0)
    assert "per-slide" in node["node_scope_reason"]

    named = slide_world["manifest"]["node_scope"]["slide_scope_nodes"]
    assert [(n["donor_id"], n["node_name"], n["n_cells"]) for n in named] == [
        ("7004", "pLN_only", 3)]
    assert slide_world["manifest"]["totals"]["nodes_slide_scope"] == 1


def test_the_slide_scope_caveat_is_carried_per_row_not_only_in_prose(slide_world, world):
    """Defect 6.  A consumer pooling nodes across donors must meet it in the data."""
    slide = pq.read_table(slide_world["out"] / "donor_id=7004/pln_node_membership.parquet")
    assert "node_scope" in slide.schema.names
    assert slide.column("node_scope").to_pylist() == ["slide"] * 3

    ordinary = pq.read_table(world["out"] / "donor_id=7001/pln_node_membership.parquet")
    scopes = ordinary.column("node_scope").to_pylist()
    # two usable nodes on that image, so neither is the slide; the unassigned cell has no node
    assert scopes == ["node"] * 5 + [None]


@pytest.mark.skipif(not (LIVE / "MANIFEST.json").exists(),
                    reason="the cohort artifact is not built in this checkout")
def test_the_built_artifact_reports_a_root_it_actually_has():
    """The published totals are the files', the root matches the run, and 6579 is in it."""
    manifest = json.loads((LIVE / "MANIFEST.json").read_text())
    on_disk = sorted(p.parent.name.split("=", 1)[1]
                     for p in LIVE.glob(f"donor_id=*/{builder.PARTITION_FILE}"))
    recon = manifest["output_root_reconciliation"]
    assert recon["donors_on_disk_after"] == on_disk
    assert recon["root_matches_this_run"] is True
    assert manifest["totals"]["cells_written"] == sum(
        pq.ParquetFile(LIVE / f"donor_id={d}" / builder.PARTITION_FILE).metadata.num_rows
        for d in on_disk)
    assert "6579" in on_disk, "6579 was reinstated on 2026-09-17; the layer must carry it"


# --------------------------------------------------------------------------------------------
# D4, 2026-09-17: a partial run must never overwrite a cohort-wide artifact.
#
# It bit twice in one day.  `build_pln_immune_pheno.py --donors 6539 --stage standout` rewrote
# cluster_calls.csv (293 rows -> 19), cluster_heterogeneity.csv (293 -> 19) and
# per_donor_separation.csv (14 -> 1) with the one donor it had processed, and earlier the same
# day `protein_rna --donors 6539` dropped 6436 and 6442 from its summary index the same way.
# The rule and the reasoning are in core/phenocycler/scripts/pln/partial_run_guard.py.
# --------------------------------------------------------------------------------------------

CLUSTER_SCRIPT = REPO / "core/phenocycler/scripts/pln/cluster_pln.py"
GUARD_SCRIPT = REPO / "core/phenocycler/scripts/pln/partial_run_guard.py"


def _load(path: Path, name: str):
    """Import a lane script by path, registered in sys.modules so dataclasses resolve."""
    import sys

    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


guard = _load(GUARD_SCRIPT, "partial_run_guard")


def test_a_named_subset_is_partial_and_a_run_naming_every_unit_is_not():
    every = guard.RunScope.build("donor", ["115", "6380"], ["115", "6380"])
    some = guard.RunScope.build("donor", ["115"], ["115", "6380"])
    whole = guard.RunScope.build("donor", None, ["115", "6380"])
    assert not every.is_partial and not whole.is_partial
    assert some.is_partial
    # a stage that cannot enumerate its units treats a named selection as partial
    assert guard.RunScope.undiscovered("donor", ["115"]).is_partial
    assert not guard.RunScope.undiscovered("donor", None).is_partial


def test_a_partial_run_is_routed_off_the_cohort_name(tmp_path):
    some = guard.RunScope.build("donor", ["6539"], ["6539", "6436", "6442"])
    dest = some.cohort_path(tmp_path, "cluster_calls.csv")
    assert dest != tmp_path / "cluster_calls.csv"
    assert dest == tmp_path / "_partial/donors=6539/cluster_calls.csv"
    # order of the selection cannot make a second directory holding the same thing
    assert (guard.RunScope.build("donor", ["b", "a"], ["a", "b", "c"]).slug
            == guard.RunScope.build("donor", ["a", "b"], ["a", "b", "c"]).slug)


def test_the_refusal_names_every_unit_it_would_have_lost():
    with pytest.raises(SystemExit) as err:
        guard.refuse_shrink(artifact="the layer", unit_kind="donor",
                            losing={"6436": "unavailable: no derivation",
                                    "6442": "unavailable: no derivation"},
                            hint="run the derivation")
    assert "6436" in str(err.value) and "6442" in str(err.value)
    guard.refuse_shrink(artifact="the layer", unit_kind="donor", losing={}, hint="")


def test_a_one_donor_cluster_run_leaves_the_cohort_profile_intact(tmp_path):
    """The `--stage profile --donors X` shape, end to end, against a populated root.

    Two donors are clustered and profiled; then the profile stage is run again for ONE of them.
    Before this rule that second run concatenated its single donor onto the cohort-wide
    `naming_profiles.parquet` and the other donor left the file without a word.
    """
    cluster = _load(CLUSTER_SCRIPT, "cluster_pln")

    cells = tmp_path / "cells_pln"
    rng = __import__("numpy").random.default_rng(7)
    ids = {}
    for donor in ("7001", "7002"):
        part = cells / f"donor_id={donor}" / "roi=pln_1"
        part.mkdir(parents=True)
        n = 40
        oid = [f"{donor}-{i:04d}" for i in range(n)]
        ids[donor] = oid
        pq.write_table(
            pa.table({
                "object_id": pa.array(oid, type=pa.string()),
                "qc_analysis_eligible": pa.array([True] * n),
                "CD3e": pa.array(rng.normal(10, 2, n)),
                "Nucleus__CD3e": pa.array(rng.normal(10, 2, n)),
                "CD20": pa.array(rng.normal(8, 2, n)),
                "Nucleus__CD20": pa.array(rng.normal(8, 2, n)),
            }),
            part / "cells.parquet",
        )

    out_root = tmp_path / "clusters"
    out_dir = out_root / "phenocycler"
    for donor in ("7001", "7002"):
        dest = out_dir / "assignments" / f"donor_id={donor}"
        dest.mkdir(parents=True)
        pq.write_table(
            pa.table({
                "object_id": pa.array(ids[donor], type=pa.string()),
                "donor_id": pa.array([donor] * len(ids[donor]), type=pa.string()),
                "cluster_id": pa.array([f"{donor}_{i % 2}" for i in range(len(ids[donor]))],
                                       type=pa.string()),
            }),
            dest / "clusters.parquet",
        )

    cluster.PHENO_CELLS = cells
    cluster.OUT_ROOT = out_root

    class Args:
        donors = None

    cluster.set_scope(guard.RunScope.undiscovered("donor", None))
    cluster.naming_profile_pheno(Args())
    cohort = out_dir / "naming_profiles.parquet"
    whole = pd.read_parquet(cohort)
    assert sorted(set(whole["donor_id"])) == ["7001", "7002"]
    before = cohort.read_bytes()

    class OneDonor:
        donors = ["7001"]

    cluster.set_scope(guard.RunScope.undiscovered("donor", OneDonor.donors))
    cluster.naming_profile_pheno(OneDonor())

    # the cohort-wide file is byte-for-byte what it was
    assert cohort.read_bytes() == before
    after = pd.read_parquet(cohort)
    assert sorted(set(after["donor_id"])) == ["7001", "7002"]
    assert len(after) == len(whole)

    # and the partial run's own output is written, in full, where its scope is visible
    scoped = out_dir / "_partial/donors=7001/naming_profiles.parquet"
    assert scoped.exists()
    assert sorted(set(pd.read_parquet(scoped)["donor_id"])) == ["7001"]
    note = json.loads((out_dir / "_partial/donors=7001/SCOPE.json").read_text())
    assert note["is_partial"] is True
    assert note["not_in_this_run"] == ["7002"]


def test_a_missing_derivation_may_not_shrink_the_published_layer(world, tmp_path):
    """D4 in this lane: 7001 is on disk, its derivation is gone, and the build must refuse.

    A donor removed by a DECISION is still removed -- that is
    `test_a_rebuild_into_a_populated_root_removes_the_partition_it_no_longer_writes`, which
    still passes.  This is the other case: the input is simply absent, and republishing a
    smaller layer under the same name is the defect.
    """
    partition = world["out"] / "donor_id=7001" / builder.PARTITION_FILE
    assert partition.exists()
    rows_before = pq.ParquetFile(partition).metadata.num_rows

    src = tmp_path / "derived/node_membership/pln_node_membership_7001.parquet"
    moved = tmp_path / "pln_node_membership_7001.parquet.moved"
    src.rename(moved)                      # the derivation was not run for this build
    try:
        with pytest.raises(SystemExit) as err:
            _rebuild(world, ROSTER_CSV, tmp_path, "roster_same.csv")
        message = str(err.value)
        assert "7001" in message and "refusing" in message.lower()
        assert partition.exists()
        assert pq.ParquetFile(partition).metadata.num_rows == rows_before
    finally:
        moved.rename(src)
