"""The v3 QuPath re-export: a per-image (cells, context) Parquet pair.

Every fixture here is synthetic and built in ``tmp_path``. The v2 CSV/TSV path is covered by
``test_cells_parquet.py``; the tests below that touch it assert only that it is unchanged
beside the new source.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import duckdb
import pandas as pd
import pytest

from phenocycler import cells_parquet as cp
from phenocycler.artifacts import (
    ChannelMapping,
    ContractError,
    QuPathCohortManifest,
    QuPathImageManifest,
    measurement_header,
)

IMAGE = "6374_v3.ome.tiff"

#: The v3 cells export: QuPath's own names, every compartment and statistic, no "(norm)".
#: Its Parent is written before islets are parented, so it says Root object for every cell.
CELLS_COLUMNS = [
    "Image", "Object ID", "Object type", "Name", "Classification", "Parent", "ROI",
    "Centroid X µm", "Centroid Y µm",
    "Cell: Area µm^2", "Cell: Solidity", "Nucleus: Area µm^2", "Nucleus/Cell area ratio",
    *(
        f"{compartment}: {marker}: {stat}"
        for compartment in ("Nucleus", "Cytoplasm", "Membrane", "Cell")
        for marker in ("DAPI", "E-cadherin")
        for stat in ("Mean", "Min", "Max", "Std.Dev.")
    ),
]

#: v3 Parent names, one per region the ingest CASE distinguishes.
PARENTS = {
    "c-core": ("Islet_7", "core"),
    "c-peri": ("Islet_7_exp20um", "peri"),
    "c-cluster": ("Cluster_12", "cluster"),
    "c-single": ("Single_3", "single"),
    "c-pln": ("Annotation (pLN)", "pLN"),
    "c-artifact": ("Annotation (Artifact)", "artifact"),
    "c-tissue": ("Annotation (Tissue)", "tissue"),
    "c-root": ("Root object (Image)", "tissue"),
    "c-other": ("Annotation (Fat)", "other"),
}


def _cells_frame(object_ids) -> pd.DataFrame:
    rows = []
    for index, oid in enumerate(object_ids):
        row = {column: 1.0 for column in CELLS_COLUMNS}
        row.update({
            "Image": IMAGE, "Object ID": oid, "Object type": "Cell", "Name": "",
            "Classification": "", "Parent": "Root object (Image)", "ROI": "Polygon",
            # Deliberately NOT the context centroids: the context file is authoritative.
            "Centroid X µm": -1.0, "Centroid Y µm": -1.0,
            "Cell: Area µm^2": 60.0 + index, "Cell: Solidity": 0.9,
            "Cell: DAPI: Mean": 500.0 + index, "Cell: E-cadherin: Mean": 700.0,
            "Nucleus: DAPI: Mean": 600.0,
        })
        rows.append(row)
    return pd.DataFrame(rows, columns=CELLS_COLUMNS)


def _context_frame(object_ids, parents=None, *, with_word=True) -> pd.DataFrame:
    word = "with " if with_word else ""
    parents = parents or {}
    return pd.DataFrame({
        "Object ID": list(object_ids),
        "Parent": [parents.get(oid, "Annotation (Tissue)") for oid in object_ids],
        "Centroid X µm": [10.0 + i for i in range(len(object_ids))],
        "Centroid Y µm": [20.0 + i for i in range(len(object_ids))],
        f"Signed distance to annotation {word}Islet µm": [-2.0] * len(object_ids),
        f"Signed distance to annotation {word}Tissue µm": [5.0] * len(object_ids),
    })


def _pair(tmp_path: Path, cell_ids, context_ids=None, parents=None, *, with_word=True):
    (tmp_path / "cells").mkdir(exist_ok=True)
    (tmp_path / "context").mkdir(exist_ok=True)
    cells = tmp_path / "cells" / "6374.parquet"
    context = tmp_path / "context" / "6374.parquet"
    _cells_frame(cell_ids).to_parquet(cells, index=False)
    _context_frame(
        cell_ids if context_ids is None else context_ids, parents, with_word=with_word
    ).to_parquet(context, index=False)
    return cells, context


def _read_partitions(cells_dir: Path) -> pd.DataFrame:
    parts = []
    for f in sorted(cells_dir.rglob("*.parquet")):
        part = pd.read_parquet(f)
        part["donor_id"] = f.parent.name.split("=", 1)[1]
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def _ingest_pair(tmp_path, cells, context, *, extra_csv=None, extra_map=None):
    out = tmp_path / "out"
    (out / "cells").mkdir(parents=True, exist_ok=True)
    paths = [cells] + ([extra_csv] if extra_csv else [])
    image_to_donor = {IMAGE: "6374", **(extra_map or {})}
    image_to_csv = {IMAGE: str(cells)}
    if extra_csv:
        image_to_csv.update({image: str(extra_csv) for image in extra_map})
    sql, mm = cp.build_sql(
        paths, out, None, 5.0,
        image_to_donor=image_to_donor, image_to_csv=image_to_csv,
        contexts={str(cells): str(context)},
    )
    duckdb.connect().execute(sql)
    return _read_partitions(out / "cells"), mm


# -- header reading -----------------------------------------------------------------------


def test_parquet_header_is_the_schema_verbatim(tmp_path):
    cells, context = _pair(tmp_path, ["a", "b"])
    assert cp.header_cols(cells) == CELLS_COLUMNS
    assert measurement_header(context)[:2] == ("Object ID", "Parent")
    assert "Signed distance to annotation with Islet µm" in cp.header_cols(context)


def test_csv_and_tsv_headers_are_still_read_by_suffix(tmp_path):
    for name, sep in (("x.csv", ","), ("x.tsv", "\t")):
        path = tmp_path / name
        path.write_text(sep.join(["Image", "Object ID", "Cell: DAPI: Mean"]) + "\n")
        assert cp.header_cols(path) == ["Image", "Object ID", "Cell: DAPI: Mean"]
        assert measurement_header(path) == ("Image", "Object ID", "Cell: DAPI: Mean")


# -- the pair join ------------------------------------------------------------------------


def test_the_pair_takes_parent_centroids_and_distances_from_context(tmp_path):
    ids = list(PARENTS)
    cells, context = _pair(tmp_path, ids, parents={k: v[0] for k, v in PARENTS.items()})
    source = cp.paired_parquet_source(cells, context)
    joined = duckdb.connect().execute(f"SELECT * FROM {source} AS s").fetchdf()
    assert sorted(joined["Object ID"]) == sorted(ids)
    by_id = joined.set_index("Object ID")
    assert by_id.loc["c-core", "Parent"] == "Islet_7"                # not "Root object"
    assert by_id.loc["c-core", "Centroid X µm"] != -1.0              # context, not cells
    assert "Signed distance to annotation with Islet µm" in joined.columns
    assert list(joined.columns).count("Parent") == 1


def test_ingest_of_a_pair_writes_regions_from_the_context_parent(tmp_path):
    ids = list(PARENTS)
    cells, context = _pair(tmp_path, ids, parents={k: v[0] for k, v in PARENTS.items()})
    df, mm = _ingest_pair(tmp_path, cells, context)
    assert set(df.donor_id) == {"6374"}
    assert sorted(df.object_id) == sorted(ids)
    assert dict(zip(df.object_id, df.cell_region, strict=True)) == {
        k: v[1] for k, v in PARENTS.items()
    }
    by_id = df.set_index("object_id")
    assert by_id.loc["c-core", "X_centroid"] >= 10.0
    assert by_id.loc["c-core", "parent_raw"] == "Islet_7"
    assert set(mm.values()) == {"DAPI", "E_cadherin"}


@pytest.mark.parametrize(
    "context_ids, fragment",
    [
        (["a", "b", "c", "d"], "row count differs: cells 3 != context 4"),
        (["a", "b"], "row count differs: cells 3 != context 2"),
        (["a", "b", "z"], "1 only in cells (e.g. ['c']), 1 only in context (e.g. ['z'])"),
        (["a", "b", "b"], "duplicate Object ID"),
    ],
)
def test_a_pair_whose_object_ids_differ_is_refused(tmp_path, context_ids, fragment):
    cells, context = _pair(tmp_path, ["a", "b", "c"], context_ids=context_ids)
    with pytest.raises(ContractError) as err:
        cp.verify_context_pair(cells, context)
    assert fragment in str(err.value)
    # ...and build_sql refuses it before any partition could be written.
    with pytest.raises(ContractError, match="does not describe one object set"):
        _ingest_pair(tmp_path, cells, context)
    assert not any((tmp_path / "out" / "cells").rglob("*.parquet"))


def test_a_matching_pair_reports_its_counts(tmp_path):
    cells, context = _pair(tmp_path, ["a", "b", "c"], context_ids=["c", "a", "b"])
    assert cp.verify_context_pair(cells, context) == {"cells_rows": 3, "context_rows": 3}


def test_a_context_without_parent_is_refused(tmp_path):
    cells, context = _pair(tmp_path, ["a"])
    pd.read_parquet(context).drop(columns=["Parent"]).to_parquet(context, index=False)
    with pytest.raises(ContractError, match="Parent"):
        cp.paired_parquet_source(cells, context)


def test_a_context_for_an_unread_export_is_refused(tmp_path):
    cells, context = _pair(tmp_path, ["a"])
    with pytest.raises(ContractError, match="not read"):
        cp.build_sql([cells], tmp_path, None, 5.0, image_to_donor={IMAGE: "6374"},
                     contexts={str(tmp_path / "elsewhere.parquet"): str(context)})


# -- region mapping -----------------------------------------------------------------------


def test_region_and_id_columns_for_every_v3_parent():
    parents = [
        "Islet_7", "Islet_7_exp20um", "Islet_ab12_exp20um", "Cluster_12", "Single_3",
        "Annotation (pLN)", "Annotation (Artifact)", "Annotation (Tissue)",
        "Root object (Image)", "Annotation (Fat)", "Cluster_x",
    ]
    frame = pd.DataFrame({"Parent": parents})  # noqa: F841  (read by duckdb)
    got = duckdb.connect().execute(
        f"SELECT \"Parent\", {cp.CELL_REGION_SQL} AS r, {cp.ISLET_NUM_SQL} AS i, "
        f"{cp.CLUSTER_NUM_SQL} AS c, {cp.SINGLE_NUM_SQL} AS s FROM frame"
    ).fetchall()
    assert {p: (r, i, c, s) for p, r, i, c, s in got} == {
        "Islet_7": ("core", "7", "", ""),
        "Islet_7_exp20um": ("peri", "7", "", ""),
        "Islet_ab12_exp20um": ("peri", "ab12", "", ""),     # a v3 uid need not be numeric
        "Cluster_12": ("cluster", "", "12", ""),
        "Single_3": ("single", "", "", "3"),
        "Annotation (pLN)": ("pLN", "", "", ""),
        "Annotation (Artifact)": ("artifact", "", "", ""),
        "Annotation (Tissue)": ("tissue", "", "", ""),
        "Root object (Image)": ("tissue", "", "", ""),
        "Annotation (Fat)": ("other", "", "", ""),
        "Cluster_x": ("other", "", "", ""),
    }


# -- (norm) is not required ---------------------------------------------------------------


def test_no_marker_is_contracted_to_qupath_norm_so_none_is_required(tmp_path):
    from phenocycler.marker_registry import load_registry

    assert not any(m.measurement_source == "qupath_norm" for m in load_registry().markers)
    assert cp.QUPATH_NORMALIZED_MARKERS == ()
    cells, _ = _pair(tmp_path, ["a"])
    assert not any("(norm)" in c for c in cp.header_cols(cells))
    assert cp.reconcile_qupath_normalized_headers([cells]) == {}


def test_a_v3_pair_with_no_norm_columns_ingests(tmp_path):
    cells, context = _pair(tmp_path, ["a", "b"])
    df, _ = _ingest_pair(tmp_path, cells, context)
    assert len(df) == 2
    assert not [c for c in df.columns if c.startswith("QuPathNorm__")]


# -- signed distance spellings ------------------------------------------------------------


@pytest.mark.parametrize("with_word", [True, False])
def test_signed_distance_column_accepts_both_spellings(tmp_path, with_word):
    _, context = _pair(tmp_path, ["a"], with_word=with_word)
    columns = cp.header_cols(context)
    word = "with " if with_word else ""
    assert cp.signed_distance_column(columns, "Islet") == (
        f"Signed distance to annotation {word}Islet µm"
    )
    assert cp.signed_distance_column(columns, "Tissue") == (
        f"Signed distance to annotation {word}Tissue µm"
    )
    assert cp.signed_distance_column(columns, "Nerve") is None


def test_the_with_spelling_is_preferred_when_both_are_present():
    both = ["Signed distance to annotation Islet µm",
            "Signed distance to annotation with Islet µm"]
    assert cp.signed_distance_column(both, "Islet") == both[1]


def test_a_v2_spelled_context_still_supplies_its_distances(tmp_path):
    cells, context = _pair(tmp_path, ["a"], with_word=False)
    joined = duckdb.connect().execute(
        f"SELECT * FROM {cp.paired_parquet_source(cells, context)} AS s"
    ).fetchdf()
    assert "Signed distance to annotation Islet µm" in joined.columns


# -- v2 beside v3 -------------------------------------------------------------------------


def _v2_csv(tmp_path: Path) -> Path:
    path = tmp_path / "Cellmeasurements.csv"
    header = ["Image", "Object ID", "Centroid X µm", "Centroid Y µm", "Parent",
              "Cell: Area µm^2", "Cell: Solidity", "Cell: DAPI: Mean", "Nucleus: DAPI: Mean",
              "Signed distance to annotation Islet µm"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerow(["0112_Scan1", "v2-a", "1", "2", "Islet_4", "60", "0.9", "5", "6", "-1"])
        w.writerow(["0112_Scan1", "v2-b", "1", "2", "Islet_4_exp20um", "60", "0.9", "5", "6",
                    "3"])
    return path


def test_a_v2_csv_and_a_v3_pair_ingest_together(tmp_path):
    cells, context = _pair(tmp_path, ["a", "b"], parents={"a": "Cluster_2"})
    v2 = _v2_csv(tmp_path)
    df, mm = _ingest_pair(tmp_path, cells, context, extra_csv=v2,
                          extra_map={"0112_Scan1": "0112"})
    got = df.set_index("object_id")
    assert got.loc["v2-a", ["donor_id", "cell_region", "islet_num"]].tolist() == [
        "0112", "core", "4"]
    assert got.loc["v2-b", "cell_region"] == "peri"
    assert got.loc["a", ["donor_id", "cell_region", "cluster_num"]].tolist() == [
        "6374", "cluster", "2"]
    # The union keeps E-cadherin, which only the v3 export carries, as null for the v2 donor.
    assert pd.isna(got.loc["v2-a", "E_cadherin"])
    assert got.loc["v2-a", "DAPI"] == 5.0


def test_a_v2_only_build_is_the_same_sql_with_or_without_the_new_argument(tmp_path):
    v2 = _v2_csv(tmp_path)
    kwargs = dict(image_to_donor={"0112_Scan1": "0112"})
    assert cp.build_sql([v2], tmp_path, None, 5.0, **kwargs) == cp.build_sql(
        [v2], tmp_path, None, 5.0, contexts={}, **kwargs)


# -- the manifest entry and the whole ingest stage ----------------------------------------


def _geojson(path: Path, ids) -> Path:
    ring = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]]
    path.write_text(json.dumps({"type": "FeatureCollection", "features": [
        {"type": "Feature", "id": oid, "geometry": {"type": "Polygon", "coordinates": [ring]},
         "nucleusGeometry": {"type": "Polygon", "coordinates": [ring]},
         "properties": {"objectType": "cell"}}
        for oid in ids
    ]}))
    return path


def _v3_image(tmp_path: Path, ids) -> QuPathImageManifest:
    cells, context = _pair(tmp_path, ids, parents={ids[0]: "Islet_1"})
    return _v3_create(tmp_path, ids, cells, context)


def _v3_create(tmp_path: Path, ids, cells: Path, context: Path) -> QuPathImageManifest:
    qptiff = tmp_path / "image.qptiff"
    qptiff.write_bytes(b"synthetic")
    return QuPathImageManifest.create(
        donor_id="6374", image_id=IMAGE, measurement_csv=cells, qptiff=qptiff,
        cell_geojson=_geojson(tmp_path / "cells.geojson", ids), nucleus_geojson=None,
        pixel_size_um_x=0.5, pixel_size_um_y=0.5, panel_id="panel-v3",
        channel_map=(ChannelMapping("DAPI", "DAPI", 0),
                     ChannelMapping("E_cadherin", "E-cadherin", 1)),
        segmentation_version="instanseg-v3", measurement_context=context,
    )


def test_a_manifest_entry_declares_and_fingerprints_both_files(tmp_path):
    image = _v3_image(tmp_path, ["a", "b"])
    context = tmp_path / "context" / "6374.parquet"
    assert image.measurement_context is not None
    assert image.measurement_context.path == str(context.resolve())
    assert image.measurement_context.size_bytes == context.stat().st_size
    assert len(image.measurement_context.sha256) == 64
    assert "Signed distance to annotation with Islet µm" in image.measurement_columns
    reloaded = QuPathImageManifest.from_dict(json.loads(json.dumps(image.to_dict())))
    assert reloaded == image
    # The context file is hashed: re-parenting one cell, cells file untouched, moves the id.
    frame = pd.read_parquet(context)
    frame.loc[1, "Parent"] = "Single_9"
    frame.to_parquet(context, index=False)
    cells = tmp_path / "cells" / "6374.parquet"
    other = _v3_create(tmp_path, ["a", "b"], cells, context)
    assert other.measurement_csv == image.measurement_csv
    assert other.measurement_context.sha256 != image.measurement_context.sha256
    assert other.content_id != image.content_id


def test_a_v2_entry_is_serialised_without_the_new_key(tmp_path):
    path = tmp_path / "cells.csv"
    path.write_text("Image,Object ID,Cell: DAPI: Mean\nimg,a,1\n")
    qptiff = tmp_path / "image.qptiff"
    qptiff.write_bytes(b"synthetic")
    image = QuPathImageManifest.create(
        donor_id="6374", image_id="img", measurement_csv=path, qptiff=qptiff,
        cell_geojson=_geojson(tmp_path / "c.geojson", ["a"]), nucleus_geojson=None,
        pixel_size_um_x=0.5, pixel_size_um_y=0.5, panel_id="p",
        channel_map=(ChannelMapping("DAPI", "DAPI", 0),), segmentation_version="v2",
    )
    assert "measurement_context" not in image.to_dict()
    cohort = QuPathCohortManifest.create(
        cohort_name="c", expected_donors=("6374",), images=(image,))
    assert "measurement_context" not in cohort.to_dict()["images"][0]


def test_a_context_may_not_pair_with_a_csv(tmp_path):
    path = tmp_path / "cells.csv"
    path.write_text("Image,Object ID\nimg,a\n")
    _, context = _pair(tmp_path, ["a"])
    with pytest.raises(ContractError, match="pairs only with a cells Parquet"):
        QuPathImageManifest.create(
            donor_id="6374", image_id="img", measurement_csv=path, qptiff=path,
            cell_geojson=_geojson(tmp_path / "c.geojson", ["a"]), nucleus_geojson=None,
            pixel_size_um_x=0.5, pixel_size_um_y=0.5, panel_id="p",
            channel_map=(ChannelMapping("DAPI", "DAPI", 0),), segmentation_version="v3",
            measurement_context=context,
        )


def test_the_ingest_stage_reads_a_v3_cohort_manifest_end_to_end(tmp_path):
    image = _v3_image(tmp_path, ["a", "b", "c"])
    cohort = QuPathCohortManifest.create(
        cohort_name="v3", expected_donors=("6374",), images=(image,))
    manifest = cohort.write_json(tmp_path / "qupath_manifest.json")
    data_dir = tmp_path / "data"
    cfg = SimpleNamespace(
        qupath_manifest=manifest, data_dir=data_dir, cells_dir=data_dir / "cells",
        cells_min_cell_area=5.0, duckdb_threads=1,
    )
    out = cp.build_cells_parquet(cfg)
    df = _read_partitions(out)
    assert sorted(df.object_id) == ["a", "b", "c"]
    assert df.set_index("object_id").loc["a", "cell_region"] == "core"
    assert (data_dir / "ingest_rejects.parquet").exists()
    assert len(pd.read_parquet(data_dir / "ingest_rejects.parquet")) == 0


def test_the_cohort_spec_accepts_a_context_beside_the_cells_parquet(tmp_path, monkeypatch):
    from phenocycler import qupath_manifest

    ids = ["a", "b"]
    _pair(tmp_path, ids)
    (tmp_path / "image.qptiff").write_bytes(b"synthetic")
    _geojson(tmp_path / "cells.geojson", ids)
    monkeypatch.setattr(
        qupath_manifest, "_channels_from_qptiff",
        lambda path: (ChannelMapping("DAPI", "DAPI", 0),
                      ChannelMapping("E_cadherin", "E-cadherin", 1)),
    )
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps({
        "spec_version": qupath_manifest.SPEC_VERSION, "cohort_name": "v3",
        "images": [{
            "donor_id": "6374", "image_id": IMAGE,
            "measurement_csv": "cells/6374.parquet",
            "measurement_context": "context/6374.parquet",
            "qptiff": "image.qptiff", "cell_geojson": "cells.geojson",
            "nucleus_geojson": None, "pixel_size_um_x": 0.5, "pixel_size_um_y": 0.5,
            "panel_id": "panel-v3", "channel_map": None, "segmentation_version": "v3",
        }],
    }))
    cohort = qupath_manifest.create_cohort_manifest(spec, tmp_path / "manifest.json")
    image = cohort.images[0]
    assert image.measurement_context.path == str(tmp_path / "context" / "6374.parquet")
    assert QuPathCohortManifest.read_json(tmp_path / "manifest.json") == cohort
