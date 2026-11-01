"""Empty-slide nulls on a whole slide that holds two sections: the null is read on glass only."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("cv2")
pytest.importorskip("shapely")

from multiplex5000.positivity import (
    nulls,  # noqa: E402
    sources,  # noqa: E402
)

PIXEL_UM = 0.5                       # the 20 um margin is 40 px
SIZE = 600
CHANNELS = ("DAPI", "CD3e")
GLASS = {"DAPI": 50, "CD3e": 100}
#: two sections on one slide, each a tissue square: (x0, y0, size) and what its tissue reads
SECTIONS = {"1001": ((80, 80, 150), {"DAPI": 3000, "CD3e": 5000}),
            "1002": ((330, 330, 190), {"DAPI": 3000, "CD3e": 9000})}
UNSCANNED = 30                       # columns x < 30 were never scanned: exact zeros


def _square(x0, y0, size):
    return [[[x0, y0], [x0 + size, y0], [x0 + size, y0 + size], [x0, y0 + size], [x0, y0]]]


def _outline(donor):
    (x0, y0, size), _ = SECTIONS[donor]
    return [{"type": "Feature", "geometry": {"type": "Polygon", "coordinates": _square(x0, y0, size)},
             "properties": {"classification": {"name": "Tissue"}}}]


def _channel(name):
    pixels = np.full((SIZE, SIZE), GLASS[name], np.uint16)
    for (x0, y0, size), tissue in SECTIONS.values():
        pixels[y0:y0 + size, x0:x0 + size] = tissue[name]
    pixels[:, :UNSCANNED] = 0
    return pixels


def test_the_null_is_glass_when_every_section_on_the_slide_is_left_out():
    outlines = _outline("1001") + _outline("1002")
    per = nulls.slide_nulls(_channel, {c: c for c in CHANNELS}, (SIZE, SIZE), outlines, PIXEL_UM)
    for c in CHANNELS:
        assert per[c]["windows"] > 0
        assert per[c]["null_p50"] == per[c]["null_p99"] == per[c]["null_p999"] == GLASS[c]


def test_another_sections_tissue_is_not_glass():
    """Typing section 1001 with only its own outline would read 1002's tissue as the empty slide."""
    own_only = nulls.slide_nulls(_channel, {c: c for c in CHANNELS}, (SIZE, SIZE), _outline("1001"), PIXEL_UM)
    assert own_only["CD3e"]["null_p999"] > GLASS["CD3e"]


def test_unscanned_slide_is_not_glass():
    outlines = _outline("1001") + _outline("1002")
    loose = nulls.empty_windows(outlines, SIZE, SIZE, PIXEL_UM)                       # the CODEX assumption
    strict = nulls.empty_windows(outlines, SIZE, SIZE, PIXEL_UM, scanned=_channel("DAPI") > 0)
    assert nulls.window_means(_channel("CD3e"), *loose).min() < GLASS["CD3e"]          # zeros entered as glass
    assert nulls.window_means(_channel("CD3e"), *strict).min() == GLASS["CD3e"]
    assert (strict[1] >= UNSCANNED + 40).all()             # and 20 um clear of the unscanned edge


def test_no_glass_left_is_refused():
    whole = [{"type": "Feature", "geometry": {"type": "Polygon", "coordinates": _square(0, 0, SIZE)}}]
    with pytest.raises(nulls.NullsError, match="no empty window"):
        nulls.slide_nulls(_channel, {c: c for c in CHANNELS}, (SIZE, SIZE), whole, PIXEL_UM)


# --------------------------------------------------------------------------- through a cohort manifest
def _slide(path):
    import tifffile

    for name in CHANNELS:
        description = ("<PerkinElmer-QPI-ImageDescription><ImageType>FullResolution</ImageType>"
                       f"<Name>{name}</Name></PerkinElmer-QPI-ImageDescription>")
        tifffile.imwrite(path, _channel(name), append=True, description=description,
                         resolution=(1e4 / PIXEL_UM, 1e4 / PIXEL_UM), resolutionunit="CENTIMETER", metadata=None)


def _manifest(tmp_path):
    """Two sections, 1001 and 1002, on one qptiff, each with one cell in its tissue."""
    from multiplex5000.qupath_manifest import create_cohort_manifest

    src = tmp_path / "sources"
    src.mkdir()
    slide = src / "slide.qptiff"
    _slide(slide)
    images = []
    for donor, ((x0, y0, _), _tissue) in SECTIONS.items():
        image = f"{donor}_Scan1.er.qptiff - resolution #1"
        geojson = src / f"{donor}.geojson"
        geojson.write_text(json.dumps({"type": "FeatureCollection", "features": [
            {"type": "Feature", "id": f"{donor}c", "properties": {"objectType": "cell"},
             "geometry": {"type": "Polygon", "coordinates": _square(x0 + 10, y0 + 10, 10)},
             "nucleusGeometry": {"type": "Polygon", "coordinates": _square(x0 + 13, y0 + 13, 4)}}]}))
        row = {"Image": image, "Object ID": f"{donor}c", "Parent": "Root object (Image)",
               "Centroid X µm": (x0 + 15) * PIXEL_UM, "Centroid Y µm": (y0 + 15) * PIXEL_UM,
               "Cell: Area µm^2": 25.0, "Cell: Solidity": 1.0, "Nucleus: Area µm^2": 4.0}
        for compartment in ("Nucleus", "Cytoplasm", "Membrane", "Cell"):
            for c in CHANNELS:
                row[f"{compartment}: {c}: Mean"] = 1.0
        cells = src / f"{donor}.parquet"
        pd.DataFrame([row]).to_parquet(cells, index=False)
        images.append({"donor_id": donor, "image_id": image, "measurement_csv": str(cells), "qptiff": str(slide),
                       "cell_geojson": str(geojson), "nucleus_geojson": str(geojson),
                       "pixel_size_um_x": PIXEL_UM, "pixel_size_um_y": PIXEL_UM, "panel_id": "synthetic-2ch",
                       "channel_map": None, "segmentation_version": "synthetic"})
    spec = src / "spec.json"
    spec.write_text(json.dumps({"spec_version": 1, "cohort_name": "two-sections", "expected_donors": list(SECTIONS),
                                "images": images}))
    create_cohort_manifest(spec, src / "manifest.json")
    tissue = tmp_path / "tissue"
    tissue.mkdir()
    for donor in SECTIONS:
        (tissue / f"{donor}.geojson").write_text(json.dumps({"type": "FeatureCollection",
                                                             "features": _outline(donor)}))
    return src / "manifest.json", tissue


def test_every_section_on_a_slide_takes_the_slides_glass_null(tmp_path):
    manifest, tissue = _manifest(tmp_path)
    table = nulls.cohort_nulls(manifest, tissue)
    assert set(table["image"]) == set(SECTIONS)
    by = table.set_index(["image", "marker"])
    for c in CHANNELS:
        assert by.loc[("1001", c), "null_p99"] == by.loc[("1002", c), "null_p99"] == GLASS[c]
    assert table["slide"].nunique() == 1
    only = nulls.cohort_nulls(manifest, tissue, sections=["1001"])
    assert set(only["image"]) == {"1001"} and (only.set_index("marker").loc["CD3e", "null_p99"] == GLASS["CD3e"])


def test_a_section_on_the_slide_without_an_outline_stops_the_nulls(tmp_path):
    manifest, tissue = _manifest(tmp_path)
    (tissue / "1002.geojson").unlink()
    with pytest.raises(nulls.NullsError, match=r"1002"):
        nulls.cohort_nulls(manifest, tissue, sections=["1001"])     # 1002's tissue would be read as glass


def test_tissue_of_no_cohort_section_can_be_named_for_its_slide(tmp_path):
    manifest, tissue = _manifest(tmp_path)
    (tissue / "1002.geojson").rename(tissue / "_slide__slide.geojson")
    for donor in ("1002",):
        (tissue / f"{donor}.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": []}))
    table = nulls.cohort_nulls(manifest, tissue)
    assert (table.loc[table["marker"] == "CD3e", "null_p99"] == GLASS["CD3e"]).all()


def test_the_protein_measurement_source_measures_caches_and_refuses_without_inputs(tmp_path):
    manifest, tissue = _manifest(tmp_path)
    cache = tmp_path / "nulls" / "empty_area_nulls.csv"
    source = sources.ProteinMeasurementSource(tmp_path / "stage", empty_nulls_csv=cache, manifest=manifest,
                                              tissue_dir=tissue)
    got = source.empty_nulls(99)
    assert got[("1001", "CD3e")] == got[("1002", "CD3e")] == GLASS["CD3e"]
    assert cache.is_file()
    cached = sources.ProteinMeasurementSource(tmp_path / "stage", empty_nulls_csv=cache)
    assert cached.empty_nulls(99.9) == {k: float(v) for k, v in source.empty_nulls(99.9).items()}
    with pytest.raises(sources.SourceError, match="no empty-slide nulls"):
        sources.ProteinMeasurementSource(tmp_path / "stage").empty_nulls(99)
