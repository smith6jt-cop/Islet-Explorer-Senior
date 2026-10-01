"""Where the typing lane reads its inputs from, once the cells come from the v3 QuPath exports.

Four readers changed on 2026-10-01 and each could put a v2 object beside a v3 cell without any
error if it regressed:

* the structure and islet objects must come from the v3 export root, never ``results/Geojson``
  (the v2 objects, where "the newest matching file wins");
* a donor's cell outlines are the file the QuPath manifest declares, not a folder globbed by name;
* a donor's ``cell_area`` comes from the same source as its cells -- the revised rebuild when
  there is one -- not from the historical run tree unconditionally;
* the fundamentals extractor reads a v3 cells Parquet into the same table the CSV path writes.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
TYPING = Path(__file__).resolve().parents[2] / "scripts/typing"
sys.path.insert(0, str(TYPING))
sys.path.insert(0, str(REPO / "scripts"))


@pytest.fixture()
def S(tmp_path, monkeypatch):
    import subpanel_io as module
    monkeypatch.setattr(module, "REPO", tmp_path)
    monkeypatch.setattr(module, "REVISED", tmp_path / "revised_cells")
    monkeypatch.setattr(module, "RUN", tmp_path / "run")
    return module


def test_objects_come_from_the_v3_export_and_never_from_the_v2_folder():
    import subpanel_io as S
    for path in (S.GEOJSON, S.ISLETS):
        assert S.EXPORT in path.parents, path
        assert "results/v3/export" in path.as_posix()
        assert "results/Geojson" not in path.as_posix()
    assert S.GEOJSON != S.ISLETS


def test_the_context_reader_takes_exactly_one_structure_file_per_image(tmp_path, monkeypatch):
    import context as CX
    monkeypatch.setattr(CX, "GEOJSON", tmp_path)
    square = {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {},
              "geometry": {"type": "Polygon", "coordinates": [[[0, 0], [0, 4], [4, 4], [4, 0], [0, 0]]]}}]}
    folder = tmp_path / "Nerve"
    folder.mkdir()
    (folder / "Nerve__6539_Scan1.er.qptiff_-_resolution_1.geojson").write_text(json.dumps(square))
    (tmp_path / "scratch").mkdir()                    # the driver's resume files sit beside, not inside
    (tmp_path / "scratch" / "Nerve__6539_Scan1.er.qptiff_-_resolution_1.scratch.jsonl").write_text("{}")
    assert len(CX.load_objects("6539_Scan1.er", "Nerve")) == 1
    (folder / "Nerve__6539_Scan1.er.qptiff_-_resolution__1.geojson").write_text(json.dumps(square))
    with pytest.raises(FileNotFoundError, match="expected 1 file, found 2"):
        CX.load_objects("6539_Scan1.er", "Nerve")


def _manifest(root: Path, rows: list[dict]) -> None:
    target = root / "core/phenocycler/data/qupath_manifest.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps({"images": rows}))


def test_cell_outlines_are_the_file_the_manifest_declares(S, tmp_path):
    declared = tmp_path / "anywhere/6539.geojson"
    _manifest(tmp_path, [{"donor_id": "6539", "cell_geojson": {"path": str(declared)}},
                         {"donor_id": "0112", "cell_geojson": {"path": str(tmp_path / "b.geojson")}}])
    assert S.cell_outlines("6539") == declared
    assert S.cell_outlines("0112").name == "b.geojson"
    with pytest.raises(FileNotFoundError, match="0 images"):
        S.cell_outlines("6457")


def _areas(path: Path, rows: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(path, index=False)


def test_cell_area_follows_the_cell_source(S, tmp_path):
    _areas(tmp_path / "run/expression/donor_id=6380/data.parquet", {"object_id": ["a"]})
    _areas(tmp_path / "run/cells/donor_id=6380/data_0.parquet", {"object_id": ["a"], "cell_area": [40.0]})
    assert S.cell_area_source("6380") == tmp_path / "run/cells/donor_id=6380/data_0.parquet"
    assert S.cell_areas("6380")["cell_area"].tolist() == [40.0]

    # a donor with a revised rebuild reads ITS areas, even though the run still has a partition
    _areas(tmp_path / "revised_cells/donor_id=6505/data.parquet",
           {"object_id": ["a", "b"], "cell_area": [55.0, 61.0], "X_centroid": [1.0, 2.0]})
    _areas(tmp_path / "run/cells/donor_id=6505/data_0.parquet", {"object_id": ["a"], "cell_area": [999.0]})
    assert S.cell_area_source("6505") == tmp_path / "revised_cells/donor_id=6505/data.parquet"
    assert S.cell_areas("6505").set_index("object_id")["cell_area"].to_dict() == {"a": 55.0, "b": 61.0}


def test_fundamentals_from_a_v3_parquet_match_the_csv_table(tmp_path, monkeypatch):
    import extract_qupath_fundamentals as E
    monkeypatch.setattr(E, "OUT", tmp_path / "out")
    columns = {c: [] for c in E.IDENTITY + E.MORPHOMETRY}
    for i, image in enumerate(["img-a", "img-a", "img-b"]):
        columns["Object ID"].append(f"id-{i}")
        columns["Image"].append(image)
        columns["Classification"].append(None)
        columns["Parent"].append("Annotation (Tissue)")
        for c in E.MORPHOMETRY:
            columns[c].append(None if (c == "Nucleus: Area µm^2" and i == 1) else float(i) + 0.5)
    columns["Cell: DAPI: Mean"] = [1.0, 2.0, 3.0]          # a column the extractor must leave behind
    source = tmp_path / "cells.parquet"
    pq.write_table(pa.table(columns), source)

    E.extract_parquet(str(source), ["A", "B"], {"img-a": "A", "img-b": "B"}, wanted={"A"})
    assert not (tmp_path / "out/donor_id=B").exists()       # not asked for
    got = pd.read_parquet(tmp_path / "out/donor_id=A/data.parquet")
    assert list(got.columns) == [E.RENAME[c] for c in E.IDENTITY + E.MORPHOMETRY]
    assert got["object_id"].tolist() == ["id-0", "id-1"]
    assert got["parent_raw"].tolist() == ["Annotation (Tissue)"] * 2
    assert got["qupath_classification"].isna().all()
    assert got["nucleus_area_um2"].isna().tolist() == [False, True]   # no nucleus stays unavailable, not 0
    assert str(got["cell_area_um2"].dtype) == "float32"

    with pytest.raises(KeyError, match="missing"):
        pq.write_table(pa.table({"Object ID": ["x"], "Image": ["img-a"]}), tmp_path / "short.parquet")
        E.extract_parquet(str(tmp_path / "short.parquet"), ["A"], {"img-a": "A"}, wanted=None)


def test_a_declared_slide_takes_its_donor_from_the_specification(tmp_path):
    import build_acquisition_runs as B
    spec = tmp_path / "cohort-spec.json"
    spec.write_text(json.dumps({"images": [
        {"donor_id": "6533", "qptiff": "/x/Images/0080870_6533_pLNpanc_Scan3.er.qptiff"},
        {"donor_id": "0112", "qptiff": "/x/Images/0081366_HDL112_pLNpanc_Scan3.er.qptiff"},
        {"donor_id": "115", "qptiff": "/x/Images/HDL115pancLN_Scan1.er.qptiff"}]}))
    declared = B.declared_donors(spec)
    assert declared == {"0080870_6533_pLNpanc_Scan3.er.qptiff": "6533",
                        "0081366_HDL112_pLNpanc_Scan3.er.qptiff": "0112",
                        "HDL115pancLN_Scan1.er.qptiff": "115"}
    # the name pattern alone cannot read a same-slide file name, which is why they had no row
    assert B.NAME.match("0080870_6533_pLNpanc_Scan3.er.qptiff") is None
    assert B.declared_donors(tmp_path / "absent.json") == {}


def test_node_membership_reads_the_v3_centroids_like_the_tsv_cut(tmp_path):
    """The IO60pLN project moved and its v3 export is one Parquet table per slide; the centroid
    columns must come back exactly as the TSV cut gave them -- as text, micro sign and all."""
    sys.path.insert(0, str(REPO / "scripts/senior"))
    import derive_pln_node_membership as N
    assert N.QP_PROJECT == Path.home() / "IO60pLN"
    assert N.GEOJSON_DIR == N.QP_PROJECT / "results/v3/export/node_geojson"

    rows = {"Image": ["a - resolution #1", "a - resolution #1", "b - resolution #1"],
            "Object ID": ["u1", "u2", "u3"], "Parent": ["Annotation (Tissue)"] * 3,
            "Centroid X µm": [10.5, 20.25, 30.0], "Centroid Y µm": [1.0, 2.0, 3.5],
            "Cell: DAPI: Mean": [5.0, 6.0, 7.0]}                    # a column the reader must leave behind
    frame = pd.DataFrame(rows)
    tables = tmp_path / "cells"
    tables.mkdir()
    for image, part in frame.groupby("Image"):
        part.to_parquet(tables / f"{image.split(' ')[0]}.parquet", index=False)
    cut = tmp_path / "cut.tsv"
    frame.drop(columns=["Cell: DAPI: Mean"]).to_csv(cut, sep="\t", index=False)

    from_parquet, from_tsv = N.read_centroids(tables), N.read_centroids(cut)
    assert list(from_parquet.columns) == list(from_tsv.columns)
    assert from_parquet["Object ID"].tolist() == ["u1", "u2", "u3"]
    assert from_parquet["Centroid X µm"].astype(float).tolist() == from_tsv["Centroid X µm"].astype(float).tolist()
    assert all(isinstance(v, str) for v in from_parquet["Centroid Y µm"])
    with pytest.raises(SystemExit, match="no \\*.parquet table"):
        (tmp_path / "empty").mkdir()
        N.read_centroids(tmp_path / "empty")


def test_the_pln_lane_reads_the_layer_and_the_project_its_own_builders_use(monkeypatch):
    """cluster_pln.py read the node-membership layer from a fixed /data/islet_layers while
    build_pln_layer.py wrote it under $ISLET_LAYERS_DIR, and build_pln_layer.py looked for the
    IO60pLN project where it no longer is."""
    import importlib
    sys.path.insert(0, str(REPO / "core/phenocycler/scripts/pln"))
    sys.path.insert(0, str(REPO / "scripts/senior"))
    monkeypatch.setenv("ISLET_LAYERS_DIR", "/somewhere/else")
    import build_pln_layer as L
    import cluster_pln as C
    import derive_pln_node_membership as N
    L, C = importlib.reload(L), importlib.reload(C)
    assert C.PLN_LAYER == Path("/somewhere/else/pln") == L.DEFAULT_OUT
    assert L.QP_PROJECT == N.QP_PROJECT and L.GEOJSON_DIR == N.GEOJSON_DIR
    monkeypatch.delenv("ISLET_LAYERS_DIR")
    assert importlib.reload(C).PLN_LAYER == Path("/data/islet_layers/pln") == importlib.reload(L).DEFAULT_OUT

