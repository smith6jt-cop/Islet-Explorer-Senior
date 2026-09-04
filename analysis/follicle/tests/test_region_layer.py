"""Tests for the analysis region layer.

The region layer decides which region every cell and vessel belongs to, so a
silent change here moves every downstream number. Nothing covered these
functions before; `tests/` held only the U-Net utilities, neither file of which
imports data_utils.

Everything runs on synthetic polygons written inline. The tests must not need
the 11.8 GB Cells.csv, the QuPath exports, or a built layer on disk.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Polygon, box, mapping

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "analysis"))
sys.path.insert(0, str(PROJECT / "scripts"))

import data_utils as du  # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures — a two-image layer with disjoint regions
# ---------------------------------------------------------------------------
@pytest.fixture
def layer():
    """Two images. Follicle inside PALS inside RedPulp, all disjoint rings."""
    def rings():
        return {
            "Follicle": [box(0, 0, 100, 100)],
            "PALS": [box(100, 0, 200, 100)],
            "Perifollicular": [box(0, 100, 100, 200)],
            "RedPulp": [box(100, 100, 200, 200)],
        }
    return {"HDL011_PC33": rings(), "HDL043_PC33": rings()}


@pytest.fixture
def objects():
    """One object per region per image, plus one outside every polygon."""
    rows = []
    for image in ("HDL011_PC33.ome.tiff", "HDL043_PC33.ome.tiff"):
        for x, y in ((50, 50), (150, 50), (50, 150), (150, 150), (500, 500)):
            rows.append({"Image": image, "Centroid X µm": float(x),
                         "Centroid Y µm": float(y)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# assign_objects_to_polygons / assign_regions_by_image
# ---------------------------------------------------------------------------
def test_assign_objects_to_polygons_labels_each_region(layer):
    coords = np.array([[50, 50], [150, 50], [50, 150], [150, 150]], float)
    got = du.assign_objects_to_polygons(coords, layer["HDL011_PC33"])
    assert list(got) == ["Follicle", "PALS", "Perifollicular", "RedPulp"]


def test_points_outside_every_polygon_are_unassigned(layer):
    """Not NaN — "Unassigned" survives a groupby, and 4-9% of real vessels
    land here, just outside the KDE tissue contour."""
    got = du.assign_objects_to_polygons(np.array([[500.0, 500.0]]),
                                        layer["HDL011_PC33"])
    assert list(got) == ["Unassigned"]


def test_assign_regions_by_image_is_per_image(layer, objects):
    labels = du.assign_regions_by_image(objects, region_layer=layer)
    assert len(labels) == len(objects)
    expected = ["Follicle", "PALS", "Perifollicular", "RedPulp", "Unassigned"]
    assert list(labels[:5]) == expected
    assert list(labels[5:]) == expected


def test_assign_regions_by_image_honours_a_shuffled_index(layer, objects):
    """Row order must follow the frame, not the groupby. A positional write
    against a non-monotonic index would scramble labels across images."""
    shuffled = objects.sample(frac=1.0, random_state=0)
    labels = pd.Series(du.assign_regions_by_image(shuffled, region_layer=layer),
                       index=shuffled.index)
    reference = pd.Series(du.assign_regions_by_image(objects, region_layer=layer),
                          index=objects.index)
    pd.testing.assert_series_equal(labels.sort_index(), reference.sort_index())


def test_image_with_no_polygons_is_all_unassigned(layer, objects):
    extra = pd.DataFrame({"Image": ["HDL018_PC25.ome.tiff"] * 3,
                          "Centroid X µm": [50.0, 150.0, 50.0],
                          "Centroid Y µm": [50.0, 50.0, 150.0]})
    labels = du.assign_regions_by_image(pd.concat([objects, extra],
                                                  ignore_index=True),
                                        region_layer=layer)
    assert set(labels[-3:]) == {"Unassigned"}


def test_alternate_coordinate_column_names(layer, objects):
    """load_single_cells spells the centroids x/y, not "Centroid X µm"."""
    renamed = objects.rename(columns={"Centroid X µm": "x", "Centroid Y µm": "y"})
    got = du.assign_regions_by_image(renamed, region_layer=layer,
                                     x_col="x", y_col="y")
    assert list(got) == list(du.assign_regions_by_image(objects, region_layer=layer))


# ---------------------------------------------------------------------------
# signed_distance_to_regions / signed_distance_by_image
# ---------------------------------------------------------------------------
def test_signed_distance_is_negative_inside():
    sq = [box(0, 0, 100, 100)]
    d = du.signed_distance_to_regions(np.array([[50.0, 50.0]]), sq)
    assert d[0] == pytest.approx(-50.0, abs=4.0)


def test_signed_distance_is_positive_outside():
    sq = [box(0, 0, 100, 100)]
    d = du.signed_distance_to_regions(np.array([[150.0, 50.0]]), sq)
    assert d[0] == pytest.approx(50.0, abs=4.0)


def test_signed_distance_without_polygons_is_nan():
    """A region a donor does not have is unmeasured, not infinitely far away."""
    d = du.signed_distance_to_regions(np.array([[0.0, 0.0], [1.0, 1.0]]), [])
    assert np.isnan(d).all()


def test_signed_distance_accuracy_tracks_boundary_step():
    """Boundaries are sampled, so error is bounded by half a step."""
    sq = [box(0, 0, 1000, 1000)]
    pt = np.array([[1200.0, 500.0]])
    fine = du.signed_distance_to_regions(pt, sq, boundary_step=1.0)
    assert fine[0] == pytest.approx(200.0, abs=1.0)


def test_signed_distance_agrees_with_the_labels(layer, objects):
    """The reason these exist: a cell labelled Follicle must be inside the
    Follicle polygon. QuPath's own columns measure to unfiltered annotations
    and so can disagree with the labels."""
    labels = du.assign_regions_by_image(objects, region_layer=layer)
    dist = du.signed_distance_by_image(objects, "Follicle", region_layer=layer)
    assert (dist[labels == "Follicle"] < 0).all()
    assert (dist[labels == "PALS"] > 0).all()


# ---------------------------------------------------------------------------
# get_vessels — polygon by default, Parent only on request
# ---------------------------------------------------------------------------
def _vessel_frame():
    return pd.DataFrame({
        "Image": ["HDL011_PC33.ome.tiff"] * 3,
        "Classification": ["SmallVessel", "SmallVessel", "Follicle"],
        "Parent": ["Root object (Image)", "Annotation (PALS)", "Root object (Image)"],
        "Centroid X µm": [50.0, 150.0, 50.0],
        "Centroid Y µm": [50.0, 150.0, 50.0],
    })


def test_get_vessels_uses_polygons_not_parent(layer):
    v = du.get_vessels(_vessel_frame(), region_layer=layer)
    assert list(v["Classification"]) == ["SmallVessel", "SmallVessel"]
    # First vessel is orphaned at the image root but sits inside the Follicle;
    # this is exactly the 1901HBMP004 failure the polygon route fixes.
    assert list(v["Region"]) == ["Follicle", "RedPulp"]


def test_get_vessels_from_parent_reproduces_the_old_labels(layer):
    v = du.get_vessels(_vessel_frame(), region_layer=layer, from_parent=True)
    assert list(v["Region"])[1] == "PALS"
    assert pd.isna(list(v["Region"])[0])   # stranded at root, region lost


# ---------------------------------------------------------------------------
# compute_density_from_regions
# ---------------------------------------------------------------------------
def _density_frame():
    rows = []
    for x, y in ((10, 10), (20, 20), (30, 30), (150, 150)):
        rows.append({"Image": "HDL011_PC33.ome.tiff", "Sample": "HDL011",
                     "Genotype": "C/C", "Classification": "SmallVessel",
                     "Centroid X µm": float(x), "Centroid Y µm": float(y),
                     "Area µm^2": 12.0})
    return pd.DataFrame(rows)


def test_density_counts_by_containment_and_areas_by_polygon(layer):
    d = du.compute_density_from_regions(_density_frame(), region_layer=layer,
                                        regions=["Follicle", "RedPulp"])
    fol = d[d["Region"] == "Follicle"].iloc[0]
    assert fol["Vessel_Count"] == 3
    assert fol["Region_Area_mm2"] == pytest.approx(100 * 100 / 1e6)
    assert fol["Density_per_mm2"] == pytest.approx(3 / (100 * 100 / 1e6))


def test_density_schema_matches_the_legacy_one(layer):
    """compute_density delegates here; downstream code reads these names."""
    d = du.compute_density_from_regions(_density_frame(), region_layer=layer)
    for col in ("Image", "Sample", "Genotype", "Region", "Vessel_Count",
                "Region_Area_um2", "Region_Area_mm2", "Density_per_mm2",
                "RedPulp_Density", "Density_Normalized"):
        assert col in d.columns, col


def test_density_reports_regions_with_no_vessels(layer):
    """A region with zero vessels must appear as 0, not vanish — otherwise a
    genuine absence reads as a missing sample."""
    d = du.compute_density_from_regions(_density_frame(), region_layer=layer,
                                        regions=["Follicle", "PALS", "RedPulp"])
    assert set(d["Region"]) == {"Follicle", "PALS", "RedPulp"}
    assert d[d["Region"] == "PALS"].iloc[0]["Vessel_Count"] == 0


def test_compute_density_delegates(layer):
    df = _density_frame()
    legacy = du.compute_density(df, region_layer=layer)
    direct = du.compute_density_from_regions(df, region_layer=layer,
                                             regions=du.MAIN_REGIONS)
    pd.testing.assert_frame_equal(legacy, direct)


# ---------------------------------------------------------------------------
# Pixel size — the values are not readable from the OME-XML
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("sample,expected", [
    ("HDL011", 0.5077663810243286),
    ("HDL043", 0.5077663810243286),
    # Their headers declare 0.5078, written by convert_to_ome_tiff.py. Wrong.
    ("1901HBMP004", 0.374),
    ("HDL073", 0.374),
    # The ForSH2B3 CODEX sections.
    ("HBMP006", 0.37740007578193524),
    ("HBMP007", 0.37740007578193524),
    ("HDL098", 0.37740007578193524),
])
def test_get_pixel_size(sample, expected):
    assert du.get_pixel_size(sample) == pytest.approx(expected)


def test_unknown_sample_gets_the_phenocycler_default():
    assert du.get_pixel_size("HDL999") == pytest.approx(du.DEFAULT_PIXEL_SIZE)


# ---------------------------------------------------------------------------
# CLASS_TO_REGION — the export/builder contract
# ---------------------------------------------------------------------------
def test_both_trabecula_spellings_map_to_one_region():
    """The CODEX sections use the singular. A mapping that handles only the
    plural silently drops their trabecular geometry into RedPulp."""
    from build_region_layer_qupath import CLASS_TO_REGION
    assert CLASS_TO_REGION["Trabecula"] == "Trabeculae"
    assert CLASS_TO_REGION["Trabeculae"] == "Trabeculae"


def test_peripheral_white_pulp_is_perifollicular_not_pals():
    """_CODEX_REGION_MAP folds it into PALS, which is wrong here: it is the
    zone around the follicle, while PALS sheaths the central arteriole.
    Folding it in doubled the CODEX donors' PALS area and ate their band."""
    from build_region_layer_qupath import CLASS_TO_REGION
    assert CLASS_TO_REGION["Peripheral_White_Pulp"] == "Perifollicular"
    assert du._CODEX_REGION_MAP.get("Peripheral_White_Pulp") == "PALS"


def test_exporter_matches_every_class_the_builder_consumes():
    """The Groovy filter and CLASS_TO_REGION have to agree, or a re-export
    silently drops a class the layer depends on."""
    from build_region_layer_qupath import CLASS_TO_REGION
    groovy = (PROJECT / "scripts" / "export_all_regions_geojson.groovy").read_text()
    line = next(l for l in groovy.splitlines() if "REGION_CLASSES" in l and "=" in l)
    block = groovy[groovy.index(line):groovy.index("] as Set") + 8]
    for cls in CLASS_TO_REGION:
        assert f"'{cls}'" in block, f"{cls} missing from the exporter's filter"


# ---------------------------------------------------------------------------
# load_region_layer round-trip
# ---------------------------------------------------------------------------
def test_load_region_layer_round_trip(tmp_path):
    feats = [{"type": "Feature",
              "properties": {"region": "Follicle", "sample": "HDL011"},
              "geometry": mapping(box(0, 0, 100, 100))},
             {"type": "Feature",
              "properties": {"region": "PALS", "sample": "HDL011"},
              "geometry": mapping(box(100, 0, 200, 100))}]
    (tmp_path / "HDL011_PC33.geojson").write_text(
        json.dumps({"type": "FeatureCollection", "features": feats}))

    layer = du.load_region_layer(regions_dir=tmp_path)
    assert set(layer) == {"HDL011_PC33"}
    assert set(layer["HDL011_PC33"]) == {"Follicle", "PALS"}
    assert layer["HDL011_PC33"]["Follicle"][0].area == pytest.approx(10_000)


def test_load_region_layer_names_the_builder_when_absent(tmp_path):
    with pytest.raises(FileNotFoundError, match="build_region_layer_qupath.py"):
        du.load_region_layer(regions_dir=tmp_path / "nope")


# ---------------------------------------------------------------------------
# Area floors
# ---------------------------------------------------------------------------
def test_follicle_floor_is_near_the_fitted_crossover():
    """scripts/fit_area_floors.py reports ~8,490 µm² against the committed
    7,945. Within 15% the choice moves neither the object count nor the area
    dropped materially; beyond that it wants looking at."""
    assert 6_500 < du.MIN_QUPATH_FOLLICLE_AREA < 9_500


def test_pals_floor_stays_far_below_the_gmm_crossover():
    """The PALS floor is geometric, not a GMM crossover. The crossover lands at
    4,141 µm² and would discard 67% of the objects and 9.5% of the PALS area,
    breaking the validation that PALS reproduce the CSV at 0.95-1.00x."""
    assert du.MIN_QUPATH_PALS_AREA < 1_500


# ---------------------------------------------------------------------------
# Sample id resolution — the layer is keyed on it
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("image,sample", [
    ("HDL011_PC33.ome.tiff", "HDL011"),
    ("HDL052SPLN_2025Aug6_Scan1.er.qptiff - resolution #1", "HDL052"),
    ("1901HBMP004_PC29.ome.tiff", "1901HBMP004"),
    ("HDL073_PC29.ome.tiff", "HDL073"),
    # The same section, spelled differently by the two CSVs. Both must resolve,
    # or the region layer and the annotations key on different names.
    ("2006_CC2B", "HBMP006"),
    ("20_06_CC2B.tif", "HBMP006"),
])
def test_extract_sample_id(image, sample):
    assert du.extract_sample_id(image) == sample


def test_extract_sample_id_rejects_unknown_alias():
    with pytest.raises(du.SampleManifestError, match="unlisted sample alias"):
        du.extract_sample_id("HDL011_typo.ome.tiff")


def test_manifest_covers_tracked_identity_sources_exactly():
    aliases = set(du.load_sample_manifest()["alias"])
    he_images = set(pd.read_csv(
        PROJECT / "27_donor_spleen_measurements1.csv", usecols=["Image"]
    )["Image"].astype(str))
    project = json.loads((PROJECT / "project.qpproj").read_text())
    qupath_images = {
        image["serverBuilder"]["metadata"]["name"]
        for image in project["images"]
    }
    assert he_images <= aliases
    assert qupath_images <= aliases


def test_manifest_rejects_duplicate_aliases(tmp_path):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "alias,sample_id,alias_type,evidence\n"
        "same,A,test,one\n"
        "same,A,test,two\n"
    )
    with pytest.raises(du.SampleManifestError, match="repeats aliases"):
        du.load_sample_manifest(manifest)


def test_validation_compares_against_the_parent_route():
    """validate_region_layer.py must ask get_vessels for the *Parent* labels.

    get_vessels defaults to the region layer, so a bare get_vessels(df) there
    makes the script compare the region layer against itself: it reports 0.0%
    disagreement for every donor and every region, which looks like a
    spectacular pass and checks nothing. This happened.
    """
    src = (PROJECT / "scripts" / "validate_region_layer.py").read_text()
    body = src[src.index("def qupath_counts"):src.index("def main")]
    assert "from_parent=True" in body, (
        "qupath_counts must use get_vessels(..., from_parent=True) or the "
        "validation is vacuous")
