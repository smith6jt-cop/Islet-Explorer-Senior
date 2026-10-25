"""Stage 03 -- protein measurement -- on a synthetic project: one tiny image, four cells.

    /data/envs/multiplex5000/bin/python -m pytest core/multiplex5000/tests/current/test_stage_measurement.py -q

The image is a real (tiny) TIFF with three named channels; its 16 cells are QuPath's v3 export, a
(cells, context) Parquet pair with every compartment's mean, and one GeoJSON with each cell's
outline and nucleus; the manifest is built from those files by ``create_cohort_manifest``. Each
cell is filled with its own intensity and the export's means are those intensities, so REDSEA's
alignment contract (more than ten cells, r >= 0.8) holds. Cell c15 has no nucleus; the recorded
exclusions take out c05.
"""

from __future__ import annotations

import csv
import hashlib
import json

import numpy as np
import pandas as pd
import pytest
import yaml

from multiplex5000 import gates, stages, tasks
from multiplex5000 import project as P

STAGE = stages.get("03_protein_measurement")

DONOR = "6374"
IMAGE = "6374_Scan1.er.qptiff - resolution #1"
CHANNELS = ("DAPI", "CD3e", "E-cadherin")
PIXEL_UM = 0.5
#: Sixteen 10 x 10 px cells in a 4 x 4 grid (25 um^2 each), each with a 4 x 4 px nucleus but c15.
CELLS = {f"c{k:02d}": (10 * (k % 4), 10 * (k // 4), k != 15) for k in range(16)}
NO_NUCLEUS, EXCLUDED = "c15", "c05"
OFFSET = 10          # the grid's corner in the image, so every cell has slide around it
SIZE = 60


def _intensity(k: int, channel_index: int) -> float:
    """What cell k holds in a channel: its own level, so every cell differs."""
    return 100.0 + 15.0 * k + 40.0 * channel_index


def _square(x0: int, y0: int, size: int) -> list:
    x1, y1 = x0 + size, y0 + size
    return [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]]


def _qptiff(path) -> None:
    import tifffile

    for index, name in enumerate(CHANNELS):
        description = ("<PerkinElmer-QPI-ImageDescription><ImageType>FullResolution</ImageType>"
                       f"<Name>{name}</Name></PerkinElmer-QPI-ImageDescription>")
        pixels = np.full((SIZE, SIZE), 5, np.uint16)
        for k, (x, y, _) in enumerate(CELLS.values()):
            pixels[OFFSET + y:OFFSET + y + 10, OFFSET + x:OFFSET + x + 10] = _intensity(k, index)
        tifffile.imwrite(path, pixels, append=True, description=description,
                         resolution=(1e4 / PIXEL_UM, 1e4 / PIXEL_UM), resolutionunit="CENTIMETER",
                         metadata=None)


def _geojson(path) -> None:
    features = []
    for oid, (x, y, has_nucleus) in CELLS.items():
        feature = {"type": "Feature", "id": oid,
                   "geometry": {"type": "Polygon",
                                "coordinates": _square(OFFSET + x, OFFSET + y, 10)},
                   "properties": {"objectType": "cell"}}
        if has_nucleus:
            feature["nucleusGeometry"] = {"type": "Polygon",
                                          "coordinates": _square(OFFSET + x + 3, OFFSET + y + 3, 4)}
        features.append(feature)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


def _exports(cells_path, context_path) -> None:
    rows, context = [], []
    for index, (oid, (x, y, has_nucleus)) in enumerate(CELLS.items()):
        row = {"Image": IMAGE, "Object ID": oid, "Object type": "Cell", "Name": "",
               "Classification": "", "Parent": "Root object (Image)", "ROI": "Polygon",
               "Centroid X µm": -1.0, "Centroid Y µm": -1.0, "Cell: Area µm^2": 25.0,
               "Cell: Solidity": 1.0, "Nucleus: Area µm^2": 4.0 if has_nucleus else np.nan,
               "Nucleus/Cell area ratio": 0.16 if has_nucleus else np.nan}
        for compartment in ("Nucleus", "Cytoplasm", "Membrane", "Cell"):
            for channel_index, channel in enumerate(CHANNELS):
                for stat in ("Mean", "Min", "Max", "Std.Dev."):
                    row[f"{compartment}: {channel}: {stat}"] = (
                        0.0 if stat == "Std.Dev." else _intensity(index, channel_index))
        rows.append(row)
        context.append({"Object ID": oid, "Parent": "Annotation (Tissue)",
                        "Centroid X µm": (OFFSET + x + 5) * PIXEL_UM,
                        "Centroid Y µm": (OFFSET + y + 5) * PIXEL_UM,
                        "Signed distance to annotation with Islet µm": 50.0,
                        "Signed distance to annotation with Tissue µm": -5.0})
    pd.DataFrame(rows).to_parquet(cells_path, index=False)
    pd.DataFrame(context).to_parquet(context_path, index=False)


def _sha(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def measurement_project(tmp_path):
    """A project whose only stage is 03, with every source it reads."""
    from multiplex5000.qupath_manifest import create_cohort_manifest

    source = tmp_path / "sources"
    source.mkdir()
    qptiff, geojson = source / "6374.qptiff", source / "6374.geojson"
    cells, context = source / "6374_cells.parquet", source / "6374_context.parquet"
    _qptiff(qptiff)
    _geojson(geojson)
    _exports(cells, context)
    spec = {"spec_version": 1, "cohort_name": "synthetic", "expected_donors": [DONOR],
            "images": [{"donor_id": DONOR, "image_id": IMAGE, "measurement_csv": str(cells),
                        "measurement_context": str(context), "qptiff": str(qptiff),
                        "cell_geojson": str(geojson), "nucleus_geojson": str(geojson),
                        "pixel_size_um_x": PIXEL_UM, "pixel_size_um_y": PIXEL_UM,
                        "panel_id": "synthetic-3ch", "channel_map": None,
                        "segmentation_version": "synthetic"}]}
    (source / "manifest_spec.json").write_text(json.dumps(spec))
    create_cohort_manifest(source / "manifest_spec.json", source / "qupath_manifest.json")
    (source / "cohort_spec.json").write_text(json.dumps({
        "spec_version": 1, "cohort_name": "synthetic", "expected_donors": [DONOR],
        "defaults": {}, "images": [{"donor_id": DONOR, "image_id": IMAGE}]}))
    # the recorded exclusions: B is a reviewed cell; D has no nucleus (geometry QC's own reason)
    exclusions = tmp_path / "exclusions"
    (exclusions / f"donor_id={DONOR}").mkdir(parents=True)
    table = exclusions / f"donor_id={DONOR}" / "data.parquet"
    reasons = {EXCLUDED: "exocrine_inside_islet_border", NO_NUCLEUS: "missing_nucleus_geometry"}
    pd.DataFrame({"object_id": list(CELLS), "analysis_eligible": [c not in reasons for c in CELLS],
                  "analysis_exclusion_reason": [reasons.get(c) for c in CELLS]}).to_parquet(
        table, index=False)
    (exclusions / "manifest.json").write_text(json.dumps({"donors": [
        {"donor_id": DONOR, "donor_id_native": DONOR, "path": str(table), "sha256": _sha(table)}]}))
    project = {"schema": "multiplex5000-project/1", "id": "synthetic", "release": "r1",
               "cohort": "io60_v3", "stages": ["03_protein_measurement"], "automatic_gates": [],
               "settings": {}, "releases_root": str(tmp_path / "releases"),
               "sources": {"cohort_spec": str(source / "cohort_spec.json"),
                           "qupath_manifest": str(source / "qupath_manifest.json"),
                           "analysis_exclusions": str(exclusions)}}
    path = tmp_path / "synthetic.yaml"
    path.write_text(yaml.safe_dump(project, sort_keys=False))
    return path


def _run(path):
    project = P.load(path)
    result = tasks.run_stage(STAGE, project)
    checks = {c["name"]: c for c in gates.load(result.path)["checks"]}
    return result.path, checks


def _partition(folder, table):
    return pd.read_parquet(next((folder / "data" / table / f"donor_id={DONOR}").glob("*.parquet")))


def test_a_measured_image_has_every_table_and_its_gate_passes(measurement_project):
    folder, checks = _run(measurement_project)
    assert all(c["passed"] for c in checks.values() if c["blocking"]), {
        n: c["value"] for n, c in checks.items() if c["blocking"] and not c["passed"]}
    for table in ("cells", "geometry_qc", "redsea", "expression", "availability"):
        assert (folder / "data" / table / f"donor_id={DONOR}").is_dir(), table
    assert not (folder / "data" / "scratch").exists(), "REDSEA's scratch is not an output"
    expression = _partition(folder, "expression").set_index("object_id")
    # c15 has no nucleus: never admitted; c05 is a recorded exclusion: taken out after REDSEA
    eligible = expression["qc_analysis_eligible"].to_dict()
    assert {c for c, ok in eligible.items() if not ok} == {NO_NUCLEUS, EXCLUDED}
    admission = pd.read_csv(folder / "nucleus_admission.csv").iloc[0]
    assert (admission["cells"], admission["geometry_admitted"], admission["excluded_by_ruling"],
            admission["admitted"]) == (16, 15, 1, 14)
    assert (folder / "gate" / "nucleus_admission.png").stat().st_size > 1000
    assert (folder / "gate" / "markers_before_after.png").stat().st_size > 1000


def test_each_marker_is_read_from_the_compartment_the_ruling_names(measurement_project):
    folder, checks = _run(measurement_project)
    with STAGE.RULING.open() as handle:
        ruling = {r["marker"]: r["measurement_compartment"] for r in csv.DictReader(handle)}
    availability = _partition(folder, "availability").set_index("marker")
    for marker in ("CD3e", "E_cadherin", "DAPI"):
        expected = ruling[marker] if ruling.get(marker) in STAGE.RULED else "Nucleus"
        assert availability.loc[marker, "compartment"] == expected, marker
    expression = _partition(folder, "expression").set_index("object_id")
    redsea = _partition(folder, "redsea").set_index("object_id")
    compartment = availability.loc["CD3e", "compartment"]
    # the value is REDSEA's, in the ruled compartment
    admitted = ["c00", "c06", "c10"]
    assert np.allclose(expression.loc[admitted, "CD3e"], redsea.loc[admitted, f"{compartment}__CD3e"],
                       equal_nan=True)
    assert sum(checks["compartments measured"]["value"].values()) == 72


def test_the_ruling_moves_registry_compartments_as_the_october_chain_did():
    payload = json.loads(STAGE.REGISTRY.read_text())
    ruling = STAGE._ruling(STAGE.RULING)
    own = STAGE.own_compartment_registry(payload, ruling)
    from collections import Counter

    assert Counter(m["compartment"] for m in own["markers"]) == {
        "Membrane": 33, "Cytoplasm": 18, "Cell": 12, "Nucleus": 9}
    assert STAGE.ruling_problems(payload, ruling) == []
    assert STAGE.ruling_problems(payload, {**ruling, "NotAMarker": "Membrane"}) == [
        "NotAMarker: not in the registry"]


def test_an_exclusion_naming_an_unknown_cell_blocks(measurement_project, tmp_path):
    table = tmp_path / "exclusions" / f"donor_id={DONOR}" / "data.parquet"
    pd.DataFrame({"object_id": ["c00", "Z"], "analysis_eligible": [True, False],
                  "analysis_exclusion_reason": [None, "exocrine_inside_islet_border"]}).to_parquet(
        table, index=False)
    manifest = tmp_path / "exclusions" / "manifest.json"
    record = json.loads(manifest.read_text())
    record["donors"][0]["sha256"] = _sha(table)
    manifest.write_text(json.dumps(record))
    _, checks = _run(measurement_project)
    check = checks["every recorded exclusion names a cell the export has"]
    assert not check["passed"] and check["blocking"]
    assert check["value"] == [f"{DONOR}: 1 excluded cell(s) the export does not have, e.g. Z"]


def test_the_stage_refuses_to_run_without_its_exclusions(measurement_project, monkeypatch):
    from multiplex5000 import roots

    monkeypatch.setattr(roots, "SITE", None)          # neither the project nor a site names them
    data = yaml.safe_load(measurement_project.read_text())
    del data["sources"]["analysis_exclusions"]
    measurement_project.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(STAGE.MeasurementError, match="analysis_exclusions"):
        _run(measurement_project)


def test_an_exclusion_of_a_cell_the_ingest_dropped_is_consistent_not_unknown():
    # the export has c99 but the ingest dropped it (under its minimum area): not an error
    stray, dropped = STAGE.exclusion_matches({"c01", "c99", "zz"}, {"c01", "c99"}, {"c01"})
    assert (stray, dropped) == (["zz"], ["c99"])
