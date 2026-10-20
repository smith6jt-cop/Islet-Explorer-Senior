"""Focused contracts for the fail-closed image-seed fallback."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import run_image_seed as IS  # noqa: E402

from phenocycler.integration.transform import Transform  # noqa: E402


def test_resized_keypoints_keep_the_exported_pixel_centre_origin():
    points = [SimpleNamespace(pt=(0.0, 0.0)), SimpleNamespace(pt=(3.0, 5.0))]
    got = IS.keypoints_um(points, native_um=1.0, origin_um=0.25, resize_scale=0.5)
    # OpenCV output p samples native (p+.5)/scale-.5, then export.json supplies the origin.
    assert np.allclose(got, [[0.75, 0.75], [6.75, 10.75]])


def test_similarity_qc_refuses_reflection_and_nonphysical_scale():
    valid = np.array([[0.999, 0.01, 200.0], [-0.01, 0.999, -300.0]])
    metrics = IS._validate_similarity(valid)
    assert not metrics["mirrored"]
    with pytest.raises(IS.ImageSeedError, match="mirrored"):
        IS._validate_similarity(np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]))
    with pytest.raises(IS.ImageSeedError, match="scale"):
        IS._validate_similarity(np.array([[0.8, 0.0, 0.0], [0.0, 0.8, 0.0]]))


def test_xenium_patch_sampling_uses_forward_xenium_to_if_direction():
    source = np.zeros((20, 20), np.float32)
    source[5, 6] = 100.0
    # Xenium (x=6.5,y=5.5) -> IF (x=9.5,y=9.5).
    affine = np.array([[1.0, 0.0, 3.0], [0.0, 1.0, 4.0]])
    meta = {"xenium_origin_um": 0.5, "if_origin_um": 0.5,
            "xenium_um_per_px": 1.0, "if_um_per_px": 1.0}
    sampled = IS._sample_xenium_patch(source, affine, meta, 0, 0, 20, 20)
    assert np.unravel_index(int(np.argmax(sampled)), sampled.shape) == (9, 9)


def test_camera_field_shape_comes_from_xml_and_exact_export_factor(tmp_path, monkeypatch):
    import tifffile

    qptiff = tmp_path / "slide.qptiff"
    qptiff.touch()
    xml = ("<PerkinElmer-QPI-ImageDescription><Name>DAPI</Name>"
           "<PixelSizeMicrons>0.25</PixelSizeMicrons><CameraSettings>"
           "<Width>1920</Width><Height>1440</Height>"
           "</CameraSettings></PerkinElmer-QPI-ImageDescription>")
    page = SimpleNamespace(description=xml)
    fake = SimpleNamespace(series=[SimpleNamespace(
        levels=[SimpleNamespace(shape=(80, 120), pages=[page])])])

    class FakeTiffFile:
        def __init__(self, _path):
            pass

        def __enter__(self):
            return fake

        def __exit__(self, *_exc):
            return False

    monkeypatch.setattr(tifffile, "TiffFile", FakeTiffFile)
    got = IS.camera_field_shape({"qptiff": str(qptiff), "if_pixel_um": 0.25,
                                 "if_um_per_px": 1.0, "if_dapi_page": 0,
                                 "if_channel_names": ["DAPI"],
                                 "if_pyramid_shapes": [[80, 120]]})
    assert (got["export_width_px"], got["export_height_px"]) == (480, 360)
    assert (got["export_width_um"], got["export_height_um"]) == (480.0, 360.0)


def _gaussian_image(shape, points):
    rows, cols = np.mgrid[:shape[0], :shape[1]]
    image = np.full(shape, 10.0)
    for x, y in points:
        image += 180.0 * np.exp(-((cols - x) ** 2 + (rows - y) ** 2) / 4.0)
    return np.clip(image, 0, 255).astype(np.uint8)


def test_tile_detector_uses_local_fields_and_excludes_camera_seams(monkeypatch):
    monkeypatch.setattr(IS, "MIN_XENIUM_NUCLEI_PER_FIELD", 5)
    monkeypatch.setattr(IS, "DETECT_MIN_Z", 1.0)
    width, height = 96, 72
    planted = [(x0 + x, y) for x0 in (0, width)
               for x in (25, 45, 65) for y in (25, 45)]
    image = _gaussian_image((height, 2 * width), planted)
    # A bright acquisition seam must not become evidence.
    image[:, width - 1:width + 1] = 255
    nuclei = np.asarray([[x + 0.5, y + 0.5] for x, y in planted])
    meta = {"if_um_per_px": 1.0, "if_origin_um": 0.5}
    camera = {"export_width_px": width, "export_height_px": height}
    peaks, _z, _field, report = IS.tile_aware_peaks(
        image, Transform(), nuclei, meta, camera)
    assert report["n_fields_measured"] == 2
    # Every planted nucleus has a recovered peak, while none can sit in the 14-um seam margin.
    from scipy.spatial import cKDTree
    assert np.max(cKDTree(peaks).query(nuclei)[0]) < 1.0
    local_x = (peaks[:, 0] - 0.5) % width
    assert np.all((local_x >= IS.DETECT_FIELD_MARGIN_UM)
                  & (local_x < width - IS.DETECT_FIELD_MARGIN_UM))


def test_enrichment_requires_replicated_fields_not_just_small_point_distances(monkeypatch):
    monkeypatch.setattr(IS, "MIN_XENIUM_NUCLEI_PER_FIELD", 5)
    monkeypatch.setattr(IS, "MIN_EVIDENCE_NUCLEI", 40)
    monkeypatch.setattr(IS, "MIN_DETECTED_PEAKS", 40)
    monkeypatch.setattr(IS, "MIN_EXCESS_MATCHES_2UM", 40)
    monkeypatch.setattr(IS, "MIN_EXCESS_MUTUAL_PAIRS", 40)
    monkeypatch.setattr(IS, "MIN_BLOCKS", 8)
    width, height = 200, 160
    points = []
    for iy in range(4):
        for ix in range(2):
            for x in (50, 90, 130):
                for y in (50, 80, 110):
                    points.append((ix * width + x + 0.5, iy * height + y + 0.5))
    nuclei = np.asarray(points)
    peaks = nuclei + np.array([0.2, -0.1])
    image = np.ones((4 * height, 2 * width), np.uint8)
    meta = {"if_um_per_px": 1.0, "if_origin_um": 0.5}
    camera = {"export_width_px": width, "export_height_px": height}
    evidence = IS.enrichment_evidence(
        Transform(), nuclei, peaks, image, meta, camera)
    assert evidence["qc"] == {"verdict": "PASS", "reasons": []}
    assert evidence["excess_mutual_pairs"] >= 40
    assert evidence["positive_blocks"] == 8


def test_local_enrichment_cannot_hide_missing_whole_slide_if_support():
    image = np.ones((200, 200), np.uint8)
    local = np.asarray([(x + 0.5, y + 0.5)
                        for x in range(50, 151, 10) for y in range(50, 151, 10)])
    outside = local + np.array([1_000.0, 0.0])
    nuclei = np.concatenate([local, outside, outside + np.array([0.0, 1_000.0])])
    evidence = IS.enrichment_evidence(
        Transform(), nuclei, local + np.array([0.1, -0.1]), image,
        {"if_um_per_px": 1.0, "if_origin_um": 0.5},
        {"export_width_px": 200, "export_height_px": 200})
    support = evidence["whole_slide_support"]
    assert support["in_if_bounds_fraction"] == pytest.approx(1 / 3)
    assert evidence["qc"]["verdict"] == "FAIL"
    assert any("nuclei in IF bounds" in reason for reason in evidence["qc"]["reasons"])


@pytest.mark.parametrize("verdict,creates", [("PASS", True), ("FAIL", False)])
def test_only_a_passing_image_seed_is_persisted_and_refined_is_untouched(
        tmp_path, monkeypatch, verdict, creates):
    donor = tmp_path / "donor_id=112"
    donor.mkdir()
    refined = donor / "refined.json"
    refined.write_text("sentinel")
    transform = Transform(affine=np.array([[1.0, 0.0, 10.0], [0.0, 1.0, 20.0]]))
    report = {"stage": "image_seed", "donor_id": "112", "serial": "0041396",
              "roi": "panc", "qc": {"verdict": verdict,
                                       "reasons": [] if creates else ["no enrichment"]}}
    monkeypatch.setattr(IS, "fit_image_seed", lambda _path: (transform, report))
    got = IS.image_seed_one(donor)
    assert (donor / "image_seed.json").exists() is creates
    assert json.loads((donor / "image_seed_report.json").read_text())["qc"]["verdict"] == verdict
    assert got["output"]["releasable"] is False
    assert got["output"]["image_seed_qc_pass"] is creates
    assert refined.read_text() == "sentinel"


def test_exact_grid_metadata_is_required_and_must_be_self_consistent():
    x = np.zeros((3, 4), np.uint16)
    y = np.zeros((5, 6), np.uint8)
    meta = {
        "donor_id": "112", "serial": "0041396", "roi": "panc",
        "xenium_um_per_px": 0.85, "if_um_per_px": 1.0,
        "xenium_origin_um": 0.10625, "if_origin_um": 0.5,
        "qptiff": "/tmp/test.qptiff", "if_pixel_um": 0.25, "if_dapi_page": 0,
        "if_channel_names": ["DAPI"],
        "pixel_size_source": "qptiff XML PixelSizeMicrons",
        "xenium_origin_check": {"brackets_cell_table": True, "base_um_per_px": 0.2125},
        "if_pyramid_shapes": [[20, 24], [10, 12], [5, 6]],
        "xenium": {"shape": [3, 4], "um_per_px": 0.85, "origin_um": 0.10625,
                    "factor": 4, "level": 2, "pyramid": "decimated",
                    "levels": [[12, 16], [6, 8], [3, 4]]},
        "if": {"shape": [5, 6], "um_per_px": 1.0, "origin_um": 0.5,
               "factor": 4, "level": 2, "pyramid": "box_mean",
               "qptiff": "/tmp/test.qptiff"},
    }
    assert IS._required_metadata(meta, x, y)["if"]["shape"] == [5, 6]
    stale = json.loads(json.dumps(meta))
    stale["if"]["origin_um"] = 0.0
    with pytest.raises(IS.ImageSeedError, match="contradicts itself"):
        IS._required_metadata(stale, x, y)


def test_failed_overwrite_invalidates_prior_candidate_after_durable_failure_report(
        tmp_path, monkeypatch):
    donor = tmp_path / "donor_id=112"
    donor.mkdir()
    prior_transform = donor / "image_seed.json"
    prior_report = donor / "image_seed_report.json"
    prior_transform.write_text("prior transform")
    prior_report.write_text("prior report")
    failed = {"stage": "image_seed", "donor_id": "112", "serial": "0041396",
              "roi": "panc", "qc": {"verdict": "FAIL", "reasons": ["weak evidence"]}}
    monkeypatch.setattr(
        IS, "fit_image_seed", lambda _path: (Transform(), failed))

    got = IS.image_seed_one(donor, overwrite=True)

    assert not prior_transform.exists()
    durable = json.loads(prior_report.read_text())
    assert durable["qc"]["verdict"] == "FAIL"
    assert durable["output"]["invalidated_prior_transform"] is True
    assert got["report_file"] == str(prior_report)
