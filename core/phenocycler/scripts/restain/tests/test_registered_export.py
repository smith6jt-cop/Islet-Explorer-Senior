"""Registration-only export contracts on small synthetic rasters."""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import registration_inputs as RI  # noqa: E402
import run_registered_export as RX  # noqa: E402

from phenocycler.integration.rasterize import GridSpec  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402


def _record(
    tmp_path,
    *,
    donor="6539",
    serial="0059865",
    pixel_um=RX.OUTPUT_UM_PER_PX,
    geojson=None,
):
    bundle = tmp_path / "xenium_bundle"
    (bundle / "morphology_focus").mkdir(parents=True)
    for relative in (
        "cells.parquet",
        "nucleus_boundaries.parquet",
        "morphology_focus/morphology_focus_0000.ome.tif",
    ):
        (bundle / relative).touch()
    return RI.RegistrationInput(
        donor_id=donor,
        serial=serial,
        roi="panc",
        qptiff=tmp_path / f"{serial}_Scan1.qptiff",
        bundle=bundle,
        geojson=geojson,
        pixel_um=float(pixel_um),
        dapi_page=0,
        channel_names=("DAPI", "Opal 570"),
        pyramid_shapes=((32, 32),),
        provenance_source="synthetic test",
    )


def _passed_gate(record, root, transform, *, stem="refined", window_fold="all"):
    donor_dir = root / f"donor_id={record.donor_id}"
    transform.save(donor_dir, stem=stem)
    transform_path = donor_dir / f"{stem}.json"
    displacement_path = donor_dir / f"{stem}_displacement.npz"
    gate_path = donor_dir / f"alignment_{stem}.json"
    gate = {
        "donor_id": record.donor_id,
        "serial": record.serial,
        "roi": record.roi,
        "registration_input": {
            "qptiff": str(record.qptiff),
            "bundle": str(record.bundle),
        },
        "alignment": {
            "transform_source": stem,
            "transform_file": str(transform_path),
            "transform_sha256": RX._sha256(transform_path),
            "displacement_file": (
                str(displacement_path) if transform.displacement is not None else None
            ),
            "displacement_sha256": (
                RX._sha256(displacement_path) if transform.displacement is not None else None
            ),
            "has_displacement": transform.displacement is not None,
            "affine": transform.affine.tolist(),
            "window_fold": window_fold,
            "verdict": "PASS",
            "verdict_why": "synthetic held-out pass",
            "gate_file": str(gate_path),
            "certifies_mapping": {
                "pixel_um_full_res": record.pixel_um,
                "source": "qptiff XML PixelSizeMicrons",
                "dapi_page": record.dapi_page,
                "dapi_channel_name": record.channel_names[record.dapi_page],
            },
        },
    }
    gate_path.write_text(json.dumps(gate))
    return gate_path


def test_pass_gate_is_identity_bound_and_transform_is_not_reconstructed(tmp_path):
    record = _record(tmp_path)
    alignment_root = tmp_path / "valis"
    transform = Transform(affine=[[1.0, 0.01, 5.0], [-0.01, 1.0, 7.0]])
    gate = _passed_gate(record, alignment_root, transform)

    passed = RX.load_passed_alignment(record, alignment_root)
    assert passed.gate_path == gate
    assert np.array_equal(passed.transform.affine, transform.affine)
    assert passed.displacement_path is None

    payload = json.loads(gate.read_text())
    payload["alignment"]["verdict"] = "INCONCLUSIVE"
    gate.write_text(json.dumps(payload))
    with pytest.raises(RX.RegisteredExportError, match="not 'PASS'"):
        RX.load_passed_alignment(record, alignment_root)


def test_transform_stem_selects_one_exact_gate_without_fallback(tmp_path):
    record = _record(tmp_path)
    alignment_root = tmp_path / "valis"
    _passed_gate(record, alignment_root, Transform(), stem="refined")
    with pytest.raises(RX.RegisteredExportError, match="alignment_intensity_refined.json"):
        RX.load_passed_alignment(record, alignment_root, transform_stem="intensity_refined")

    expected = Transform(affine=[[1.0, 0.0, 3.0], [0.0, 1.0, -2.0]])
    _passed_gate(record, alignment_root, expected, stem="intensity_refined")
    selected = RX.load_passed_alignment(record, alignment_root, transform_stem="intensity_refined")
    assert selected.stem == "intensity_refined"
    assert np.array_equal(selected.transform.affine, expected.affine)


def test_spatial_cv_gate_must_use_the_distinct_withheld_fold(tmp_path):
    record = _record(tmp_path)
    alignment_root = tmp_path / "valis"
    transform = Transform(meta={"training_fold": "even", "withheld_gate_fold": "odd"})
    gate_path = _passed_gate(
        record,
        alignment_root,
        transform,
        stem="intensity_refined",
        window_fold="odd",
    )
    passed = RX.load_passed_alignment(record, alignment_root, transform_stem="intensity_refined")
    assert passed.transform.meta["training_fold"] == "even"

    payload = json.loads(gate_path.read_text())
    payload["alignment"]["window_fold"] = "even"
    gate_path.write_text(json.dumps(payload))
    with pytest.raises(RX.RegisteredExportError, match="expected independently withheld fold"):
        RX.load_passed_alignment(record, alignment_root, transform_stem="intensity_refined")


def test_gate_hash_refuses_transform_bytes_changed_after_validation(tmp_path):
    record = _record(tmp_path)
    alignment_root = tmp_path / "valis"
    _passed_gate(record, alignment_root, Transform())
    transform_path = alignment_root / "donor_id=6539/refined.json"
    payload = json.loads(transform_path.read_text())
    payload["meta"]["changed_after_gate"] = True
    transform_path.write_text(json.dumps(payload))
    with pytest.raises(RX.RegisteredExportError, match="transform_sha256 does not match"):
        RX.load_passed_alignment(record, alignment_root)


def test_gate_displacement_path_must_name_the_exact_selected_sibling(tmp_path):
    record = _record(tmp_path)
    alignment_root = tmp_path / "valis"
    displacement = np.zeros((2, 4, 5), dtype=np.float32)
    transform = Transform(
        displacement=displacement,
        grid=GridSpec(0.0, 0.0, 1.0, 4, 5),
    )
    gate_path = _passed_gate(record, alignment_root, transform)
    payload = json.loads(gate_path.read_text())
    payload["alignment"]["displacement_file"] = str(gate_path.parent / "other.npz")
    gate_path.write_text(json.dumps(payload))

    with pytest.raises(RX.RegisteredExportError, match="declared displacement"):
        RX.load_passed_alignment(record, alignment_root)


def test_identity_transform_reproduces_raw_channels_in_native_order():
    source = np.arange(2 * 21 * 19, dtype=np.uint8).reshape(2, 21, 19)
    warped, rec = RX._warp_store_to_xenium_grid(
        source, Transform(), input_um_per_px=RX.OUTPUT_UM_PER_PX, height=21, width=19, tile=16
    )
    assert np.array_equal(warped, source)
    assert rec["includes_displacement"] is False
    assert rec["outside_if_fraction"] == 0.0


def test_pull_warp_direction_and_half_pixel_convention_match_the_affine_reference():
    """Xenium output centres must be mapped forward through Xenium->IF Transform.apply."""
    rng = np.random.default_rng(4)
    source = rng.integers(0, 256, size=(1, 47, 41), dtype=np.uint8)
    affine = np.array([[0.997, 0.015, 0.4], [-0.012, 1.002, 0.7]])
    transform = Transform(affine)
    warped, _ = RX._warp_store_to_xenium_grid(
        source, transform, input_um_per_px=0.25, height=45, width=43, tile=16
    )

    # explorer.alignment_matrix expresses the same Xenium-um -> IF-um affine as the inverse
    # IF-pixel -> Xenium-pixel matrix used by the previously validated affine-only exporter.
    matrix = RX.EX.alignment_matrix(affine, if_um_per_px=0.25, if_origin_um=0.125)
    reference = RX.EX.warp_plane_to_xenium_grid(source[0], matrix, 45, 43, tile=16)
    assert np.array_equal(warped[0], reference)


def test_displacement_field_is_applied_instead_of_silently_dropped():
    source = np.arange(18 * 20, dtype=np.uint8).reshape(1, 18, 20)
    displacement = np.zeros((2, 18, 20), dtype=np.float32)
    displacement[0] = RX.OUTPUT_UM_PER_PX
    transform = Transform(
        displacement=displacement,
        grid=GridSpec(0.0, 0.0, RX.OUTPUT_UM_PER_PX, 18, 20),
    )
    warped, rec = RX._warp_store_to_xenium_grid(
        source, transform, input_um_per_px=RX.OUTPUT_UM_PER_PX, height=18, width=20, tile=16
    )
    expected = np.zeros_like(source)
    expected[:, :, :-1] = source[:, :, 1:]
    assert np.array_equal(warped, expected)
    assert rec["includes_displacement"] is True


def test_displacement_is_sampled_on_its_micron_grid_at_output_pixel_centres():
    """Pin the field grid convention separately from the constant-field direction test."""
    displacement = np.zeros((2, 5, 7), dtype=np.float32)
    displacement[0] = np.arange(7, dtype=np.float32)[None, :]
    displacement[1] = 0.25 * np.arange(5, dtype=np.float32)[:, None]
    transform = Transform(
        displacement=displacement,
        grid=GridSpec(0.0, 0.0, 1.0, 5, 7),
    )
    map_x, map_y = RX._source_map_for_tile(
        transform,
        x0=0,
        y0=0,
        width=3,
        height=2,
        input_um_per_px=1.0,
        output_um_per_px=1.0,
    )
    # Output centres x=.5,1.5,2.5 sample dx=.5,1.5,2.5 on the displacement grid.
    np.testing.assert_allclose(map_x, [[0.5, 2.5, 4.5], [0.5, 2.5, 4.5]])
    # Output centres y=.5,1.5 sample dy=.125,.375 before the IF half-pixel subtraction.
    np.testing.assert_allclose(map_y, [[0.125, 0.125, 0.125], [1.375, 1.375, 1.375]])


def test_export_is_atomic_raw_only_and_reopens_every_level(tmp_path, monkeypatch):
    import tifffile

    record = _record(tmp_path)
    source = np.stack(
        [
            np.arange(32 * 32, dtype=np.uint8).reshape(32, 32),
            np.flipud(np.arange(32 * 32, dtype=np.uint8).reshape(32, 32)),
        ]
    )
    tifffile.imwrite(
        record.qptiff,
        source,
        ome=True,
        tile=(16, 16),
        photometric="minisblack",
        metadata={
            "axes": "CYX",
            "Channel": {"Name": list(record.channel_names)},
            "PhysicalSizeX": record.pixel_um,
            "PhysicalSizeY": record.pixel_um,
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeYUnit": "µm",
        },
    )

    alignment_root = tmp_path / "valis"
    _passed_gate(record, alignment_root, Transform())
    morphology = {
        "file": str(record.bundle / "morphology_focus/morphology_focus_0000.ome.tif"),
        "height": 32,
        "width": 32,
        "n_levels": 2,
        "level_shapes": [[32, 32], [16, 16]],
        "um_per_px": RX.OUTPUT_UM_PER_PX,
    }
    monkeypatch.setattr(RX.EX, "morphology_grid", lambda _bundle: morphology)
    bundle_root = tmp_path / "bundles"
    plan = RX.prepare_export(
        record, alignment_root=alignment_root, bundle_root=bundle_root, compression="zlib", tile=16
    )
    result = RX.export_one(plan, compression="zlib", workers=1, tile=16)

    assert result["status"] == "EXPORTED"
    destination = bundle_root / f"donor_{record.donor_id}_{record.serial}"
    assert destination.is_dir()
    assert not (destination / "cell_groups").exists()
    assert not (bundle_root / f".{destination.name}.partial").exists()
    manifest = json.loads((destination / "manifest.json").read_text())
    qc = json.loads((destination / "qc.json").read_text())
    assert manifest["mode"] == "registration_only"
    assert [entry["name"] for entry in manifest["input"]["raw_channels"]] == list(
        record.channel_names
    )
    assert all(entry["role"] is None for entry in manifest["input"]["raw_channels"])
    assert manifest["channel_roles"]["status"] == "unassigned"
    assert manifest["positivity"]["status"] == "not_generated"
    assert manifest["cell_groups"] == {"status": "not_generated", "path": None}
    assert qc["status"] == "PASS"
    assert len(qc["output"]["decoded_tiles"]) == 2 * len(record.channel_names)

    image = destination / "images" / f"{record.serial}_IF_on_xenium_grid.ome.tif"
    with tifffile.TiffFile(image) as tif:
        assert tif.series[0].shape == (2, 32, 32)
        assert [level.shape for level in tif.series[0].levels] == [(2, 32, 32), (2, 16, 16)]
        assert np.array_equal(tif.series[0].levels[0].asarray(), source)

    with pytest.raises(RX.RegisteredExportError, match="already exists"):
        RX.prepare_export(
            record,
            alignment_root=alignment_root,
            bundle_root=bundle_root,
            compression="zlib",
            tile=16,
        )
    skip = RX.prepare_export(
        record,
        alignment_root=alignment_root,
        bundle_root=bundle_root,
        skip_existing=True,
        compression="zlib",
        tile=16,
    )
    assert (
        RX.export_one(skip, compression="zlib", workers=1, tile=16)["status"] == "SKIPPED_EXISTING"
    )

    manifest_path = destination / "manifest.json"
    stale_manifest = json.loads(manifest_path.read_text())
    stale_manifest["input"]["qptiff_bytes"] += 1
    manifest_path.write_text(json.dumps(stale_manifest))
    with pytest.raises(RX.RegisteredExportError, match="input.qptiff_bytes"):
        RX.export_one(skip, compression="zlib", workers=1, tile=16)


def test_completed_reference_donors_are_protected_even_with_overwrite(tmp_path):
    record = _record(tmp_path, donor="6436", serial="0059907")
    with pytest.raises(RX.RegisteredExportError, match="protected completed"):
        RX.prepare_export(record, alignment_root=tmp_path, bundle_root=tmp_path, overwrite=True)


def test_export_carries_same_section_islets_into_an_explorer_sidecar(tmp_path, monkeypatch):
    import tifffile

    source_geojson = tmp_path / "0059865_islets.geojson"
    source_geojson.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "id": "source-islet",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [
                                [[2, 2], [8, 2], [8, 8], [2, 8], [2, 2]]
                            ],
                        },
                        "properties": {
                            "objectType": "annotation",
                            "classification": {"name": "Islet", "color": [255, 0, 0]},
                        },
                    }
                ],
            }
        )
    )
    # A deliberately large source calibration and scale keep a >50-um source polygon inside
    # the tiny synthetic 32-pixel Xenium output while exercising both unit conversions.
    record = _record(tmp_path, pixel_um=10.0, geojson=source_geojson)
    source_image = np.arange(2 * 32 * 32, dtype=np.uint8).reshape(2, 32, 32)
    tifffile.imwrite(
        record.qptiff,
        source_image,
        ome=True,
        tile=(16, 16),
        photometric="minisblack",
        metadata={
            "axes": "CYX",
            "Channel": {"Name": list(record.channel_names)},
            "PhysicalSizeX": record.pixel_um,
            "PhysicalSizeY": record.pixel_um,
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeYUnit": "µm",
        },
    )

    alignment_root = tmp_path / "valis"
    transform = Transform(affine=[[20, 0, 0], [0, 20, 0]])
    _passed_gate(record, alignment_root, transform)
    morphology = {
        "file": str(record.bundle / "morphology_focus/morphology_focus_0000.ome.tif"),
        "height": 32,
        "width": 32,
        "n_levels": 2,
        "level_shapes": [[32, 32], [16, 16]],
        "um_per_px": RX.OUTPUT_UM_PER_PX,
    }
    monkeypatch.setattr(RX.EX, "morphology_grid", lambda _bundle: morphology)
    bundle_root = tmp_path / "bundles"
    plan = RX.prepare_export(record, alignment_root=alignment_root, bundle_root=bundle_root, tile=16)
    RX.export_one(plan, workers=1, tile=16)

    destination = bundle_root / "donor_6539_0059865"
    manifest = json.loads((destination / "manifest.json").read_text())
    qc = json.loads((destination / "qc.json").read_text())
    annotation = manifest["islet_annotations"]
    annotation_path = destination / annotation["file"]
    assert annotation["status"] == "generated"
    assert annotation["source_um_per_px"] == 10.0
    assert annotation["output_um_per_px"] == RX.OUTPUT_UM_PER_PX
    assert annotation["counts"]["exported_islet_polygons"] == 1
    assert annotation["transform_sha256"] == manifest["alignment"]["transform_sha256"]
    assert manifest["uses_qupath_annotations"] is True
    assert annotation_path.is_file()
    assert qc["islet_annotations"]["validation"]["status"] == "PASS"
    readme = (destination / "README.md").read_text()
    assert "**Annotations**" in readme and "**pixels**" in readme

    skipped = RX.prepare_export(
        record,
        alignment_root=alignment_root,
        bundle_root=bundle_root,
        skip_existing=True,
        tile=16,
    )
    assert RX.export_one(skipped, workers=1, tile=16)["annotation_validation"]["status"] == "PASS"

    manifest["islet_annotations"]["minimum_equivalent_diameter_um"] = 40.0
    (destination / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(RX.RegisteredExportError, match="minimum_equivalent_diameter_um.*stale"):
        RX.export_one(skipped, workers=1, tile=16)


def test_malformed_annotation_aborts_before_image_render_and_leaves_no_partial(
    tmp_path, monkeypatch
):
    source_geojson = tmp_path / "0059865_bad.geojson"
    source_geojson.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [
                                [[0, 0], [300, 300], [0, 300], [300, 0], [0, 0]]
                            ],
                        },
                        "properties": {"classification": {"name": "Islet"}},
                    }
                ],
            }
        )
    )
    record = _record(tmp_path, pixel_um=0.25, geojson=source_geojson)
    record.qptiff.touch()
    alignment_root = tmp_path / "valis"
    _passed_gate(record, alignment_root, Transform())
    morphology = {
        "file": str(record.bundle / "morphology_focus/morphology_focus_0000.ome.tif"),
        "height": 1_000,
        "width": 1_000,
        "n_levels": 1,
        "level_shapes": [[1_000, 1_000]],
        "um_per_px": RX.OUTPUT_UM_PER_PX,
    }
    monkeypatch.setattr(RX.EX, "morphology_grid", lambda _bundle: morphology)

    def image_render_must_not_start(*_args, **_kwargs):
        raise AssertionError("image render started before annotation validation")

    monkeypatch.setattr(RX, "warp_qptiff_to_xenium_level0", image_render_must_not_start)
    bundle_root = tmp_path / "bundles"
    plan = RX.prepare_export(record, alignment_root=alignment_root, bundle_root=bundle_root)
    with pytest.raises(
        RX.RegisteredExportError,
        match="annotation export failed.*(?:simple valid|empty or has zero area)",
    ):
        RX.export_one(plan, workers=1)
    destination = bundle_root / "donor_6539_0059865"
    assert not destination.exists()
    assert not (bundle_root / f".{destination.name}.partial").exists()
