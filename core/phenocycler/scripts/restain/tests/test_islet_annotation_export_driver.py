"""Annotation-only augmentation of existing Xenium Explorer bundles."""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import islet_annotations as IA  # noqa: E402
import registration_inputs as RI  # noqa: E402
import run_islet_annotation_export as AX  # noqa: E402

from phenocycler.integration.transform import Transform  # noqa: E402


def _source_geojson(path: pathlib.Path) -> pathlib.Path:
    feature = {
        "type": "Feature",
        "id": "source-islet",
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[10, 10], [70, 10], [70, 70], [10, 70], [10, 10]]],
        },
        "properties": {
            "objectType": "annotation",
            "classification": {"name": "Islet", "color": [255, 0, 0]},
        },
    }
    path.write_text(json.dumps({"type": "FeatureCollection", "features": [feature]}))
    return path


def _record(tmp_path: pathlib.Path, *, donor: str = "6539", serial: str = "0059865"):
    qptiff = tmp_path / f"{serial}_Scan1.qptiff"
    qptiff.touch()
    geojson = _source_geojson(tmp_path / f"{serial}_islets.geojson")
    xenium = tmp_path / "xenium" / serial
    (xenium / "morphology_focus").mkdir(parents=True)
    morphology_file = xenium / "morphology_focus" / "morphology_focus_0000.ome.tif"
    morphology_file.touch()
    record = RI.RegistrationInput(
        donor_id=donor,
        serial=serial,
        roi="panc",
        qptiff=qptiff,
        bundle=xenium,
        geojson=geojson,
        pixel_um=1.0,
        dapi_page=0,
        channel_names=("DAPI", "Opal 570"),
        pyramid_shapes=((100, 100),),
        provenance_source="synthetic inventory",
    )
    morphology = {
        "file": str(morphology_file),
        "height": 100,
        "width": 100,
        "n_levels": 1,
        "level_shapes": [[100, 100]],
        "um_per_px": 1.0,
    }
    return record, morphology


def _registered_image(path: pathlib.Path, *, shape=(100, 100), um_per_px=1.0) -> None:
    import tifffile

    path.parent.mkdir(parents=True, exist_ok=True)
    data = np.zeros((2, *shape), dtype=np.uint8)
    tifffile.imwrite(
        path,
        data,
        ome=True,
        bigtiff=True,
        photometric="minisblack",
        metadata={
            "axes": "CYX",
            "Channel": {"Name": ["DAPI", "Opal 570"]},
            "PhysicalSizeX": float(um_per_px),
            "PhysicalSizeY": float(um_per_px),
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeYUnit": "µm",
        },
    )


def _save_transform(record, alignment_root, transform=None, *, stem="refined"):
    transform = transform or Transform()
    donor_dir = alignment_root / f"donor_id={record.donor_id}"
    transform.save(donor_dir, stem=stem)
    transform_path = donor_dir / f"{stem}.json"
    displacement_path = donor_dir / f"{stem}_displacement.npz"
    return transform, transform_path, displacement_path


def _new_bundle(tmp_path, monkeypatch, *, donor="6539", serial="0059865"):
    record, morphology = _record(tmp_path, donor=donor, serial=serial)
    monkeypatch.setattr(AX.EX, "morphology_grid", lambda _bundle: dict(morphology))
    bundle_root = tmp_path / "explorer_bundle"
    destination = bundle_root / f"donor_{donor}_{serial}"
    image = destination / "images" / f"{serial}_IF_on_xenium_grid.ome.tif"
    _registered_image(image)
    (destination / "cell_groups").mkdir()
    (destination / "cell_groups" / "biological.csv").write_text(
        "cell_id,group,color\na,Positive,#ff0000\n"
    )

    alignment_root = tmp_path / "valis"
    transform, transform_path, displacement_path = _save_transform(record, alignment_root)
    gate_path = transform_path.parent / "alignment_refined.json"
    gate = {
        "donor_id": donor,
        "serial": serial,
        "roi": "panc",
        "registration_input": {
            "qptiff": str(record.qptiff),
            "bundle": str(record.bundle),
        },
        "alignment": {
            "transform_source": "refined",
            "transform_file": str(transform_path),
            "transform_sha256": AX._sha256(transform_path),
            "displacement_file": None,
            "displacement_sha256": None,
            "has_displacement": False,
            "affine": transform.affine.tolist(),
            "verdict": "PASS",
        },
    }
    gate_path.write_text(json.dumps(gate))
    manifest = {
        "schema_version": 1,
        "mode": "registration_only",
        "donor_id": donor,
        "xenium_serial": serial,
        "roi": "panc",
        "input": {"qptiff": str(record.qptiff), "pixel_size_um": record.pixel_um},
        "xenium_bundle": str(record.bundle),
        "morphology": dict(morphology),
        "alignment": {
            "gate_file": str(gate_path),
            "gate_sha256": AX._sha256(gate_path),
            "verdict": "PASS",
            "transform_source": "refined",
            "transform_file": str(transform_path),
            "transform_sha256": AX._sha256(transform_path),
            "displacement_file": None,
            "displacement_sha256": None,
            "has_displacement": False,
            "affine": transform.affine.tolist(),
        },
        "image_on_xenium_grid": {
            "file": f"images/{image.name}",
            "width": 100,
            "height": 100,
            "n_levels": 1,
            "level_shapes": [[100, 100]],
            "um_per_px": 1.0,
        },
        "positivity": {"status": "not_generated"},
        "cell_groups": {"status": "not_generated", "path": None},
        "uses_qupath_annotations": False,
    }
    (destination / "manifest.json").write_text(json.dumps(manifest))
    (destination / "README.md").write_text("# Existing bundle\n\nOriginal biological text.\n")
    (destination / "qc.json").write_text(
        json.dumps({"status": "PASS", "existing_image_check": "preserve me"})
    )
    return record, morphology, alignment_root, bundle_root, destination


def _legacy_bundle(tmp_path, monkeypatch, *, donor="6436", serial="0059907"):
    record, morphology = _record(tmp_path, donor=donor, serial=serial)
    monkeypatch.setattr(AX.EX, "morphology_grid", lambda _bundle: dict(morphology))
    bundle_root = tmp_path / "explorer_bundle"
    destination = bundle_root / f"donor_{donor}_{serial}"
    image = destination / "images" / f"{serial}_IF_on_xenium_grid.ome.tif"
    _registered_image(image)
    (destination / "cell_groups").mkdir()
    (destination / "cell_groups" / "hormone_positivity.csv").write_text("legacy biology\n")
    alignment_root = tmp_path / "valis"
    transform, transform_path, _ = _save_transform(record, alignment_root)
    manifest = {
        "donor_id": donor,
        "xenium_serial": serial,
        "roi": "panc",
        "qptiff": str(record.qptiff),
        "xenium_bundle": str(record.bundle),
        "morphology": dict(morphology),
        "transform": {"source": "refined", **transform.to_dict()},
        "image": {"levels": [{"level": 0, "um_per_px": record.pixel_um}]},
        "alignment": {"if_um_per_px": record.pixel_um},
        "image_on_xenium_grid": {
            "file": f"images/{image.name}",
            "width": 100,
            "height": 100,
            "n_levels": 1,
            "level_shapes": [[100, 100]],
            "um_per_px": 1.0,
        },
        "positivity": {"file": "cell_groups/hormone_positivity.csv"},
        "uses_qupath_annotations": False,
    }
    (destination / "manifest.json").write_text(json.dumps(manifest))
    (destination / "README.md").write_text(
        "# Protected reviewed bundle\n\n"
        f"{AX.LEGACY_NO_ANNOTATION_TEXT}\n"
    )
    return record, morphology, alignment_root, bundle_root, destination, transform_path


def test_registration_only_bundle_gets_hashed_validated_sidecar_without_payload_changes(
    tmp_path, monkeypatch
):
    record, _, alignment_root, bundle_root, destination = _new_bundle(tmp_path, monkeypatch)
    image = next((destination / "images").iterdir())
    group = destination / "cell_groups" / "biological.csv"
    image_before = (image.read_bytes(), image.stat().st_mtime_ns)
    group_before = group.read_bytes()
    original_manifest = json.loads((destination / "manifest.json").read_text())

    plan = AX.prepare_annotation_export(
        record, alignment_root=alignment_root, bundle_root=bundle_root
    )
    result = AX.export_one(plan)

    sidecar = destination / "annotations" / f"{record.serial}_islets_on_xenium_grid.geojson"
    assert result["status"] == "EXPORTED" and sidecar.is_file()
    validation = IA.validate_islet_geojson(sidecar, width=100, height=100, expected_features=1)
    manifest = json.loads((destination / "manifest.json").read_text())
    annotation = manifest["islet_annotations"]
    assert annotation["source_sha256"] == AX._sha256(record.geojson)
    assert annotation["output_sha256"] == AX._sha256(sidecar) == validation["sha256"]
    assert annotation["transform"]["sha256"] == AX._sha256(
        alignment_root / f"donor_id={record.donor_id}/refined.json"
    )
    assert annotation["alignment_gate_sha256"] == manifest["alignment"]["gate_sha256"]
    assert annotation["transform_sha256"] == manifest["alignment"]["transform_sha256"]
    assert annotation["displacement_sha256"] is None
    assert annotation["evidence_policy"]["pass_gate_revalidated"] is True
    assert annotation["output_coordinate_space"] == "Xenium morphology level-0 pixel edges"
    assert manifest["positivity"] == original_manifest["positivity"]
    assert manifest["cell_groups"] == original_manifest["cell_groups"]
    assert image.read_bytes() == image_before[0]
    assert image.stat().st_mtime_ns == image_before[1]
    assert group.read_bytes() == group_before
    qc = json.loads((destination / "qc.json").read_text())
    assert qc["existing_image_check"] == "preserve me"
    assert qc["islet_annotations"]["status"] == "PASS"
    readme = (destination / "README.md").read_text()
    assert "Original biological text." in readme
    assert "Associate it with the **Xenium In Situ morphology image**" in readme
    assert "choose **pixels** if a unit prompt appears" in readme
    assert not (destination / ".islet_annotations.partial").exists()
    assert not (destination / ".islet_annotations.previous").exists()


def test_existing_annotations_refused_by_default_and_narrow_overwrite_preserves_other_files(
    tmp_path, monkeypatch
):
    record, _, alignment_root, bundle_root, destination = _new_bundle(tmp_path, monkeypatch)
    AX.export_one(
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)
    )
    unrelated = destination / "annotations" / "reviewer_notes.geojson"
    unrelated.write_bytes(b"reviewer-owned annotation bytes")
    image = next((destination / "images").iterdir())
    image_before = (image.stat().st_size, image.stat().st_mtime_ns)

    with pytest.raises(AX.AnnotationAugmentationError, match="already exists"):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)

    plan = AX.prepare_annotation_export(
        record,
        alignment_root=alignment_root,
        bundle_root=bundle_root,
        overwrite_annotations=True,
    )
    AX.export_one(plan)
    assert unrelated.read_bytes() == b"reviewer-owned annotation bytes"
    assert (image.stat().st_size, image.stat().st_mtime_ns) == image_before
    assert (destination / "README.md").read_text().count(AX.README_BEGIN) == 1


def test_protected_legacy_bundle_uses_exact_serialized_transform_without_gate(
    tmp_path, monkeypatch
):
    record, _, alignment_root, bundle_root, destination, transform_path = _legacy_bundle(
        tmp_path, monkeypatch
    )
    group = destination / "cell_groups" / "hormone_positivity.csv"
    group_before = group.read_bytes()

    plan = AX.prepare_annotation_export(
        record, alignment_root=alignment_root, bundle_root=bundle_root
    )
    assert plan.transform_evidence["manifest_schema"] == "protected_legacy"
    assert plan.transform_evidence["sha256"] == AX._sha256(transform_path)
    assert plan.transform_evidence["evidence_policy"]["protected_legacy_exception"] is True
    AX.export_one(plan)

    manifest = json.loads((destination / "manifest.json").read_text())
    assert manifest["positivity"]["file"] == "cell_groups/hormone_positivity.csv"
    assert manifest["islet_annotations"]["transform"]["gate"] is None
    assert manifest["islet_annotations"]["alignment_gate_sha256"] is None
    assert group.read_bytes() == group_before
    readme = (destination / "README.md").read_text()
    assert AX.LEGACY_IMAGE_ANNOTATION_TEXT in readme
    assert AX.LEGACY_NO_ANNOTATION_TEXT not in readme
    dedicated_qc = destination / AX.ANNOTATION_QC_FILENAME
    assert json.loads(dedicated_qc.read_text())["scope"] == "islet annotation augmentation only"


def test_protected_legacy_ignores_but_records_unreferenced_displacement_file(tmp_path, monkeypatch):
    record, _, alignment_root, bundle_root, _, transform_path = _legacy_bundle(
        tmp_path, monkeypatch
    )
    stale_displacement = transform_path.with_name(f"{transform_path.stem}_displacement.npz")
    stale_displacement.write_bytes(b"unreferenced legacy artifact")

    plan = AX.prepare_annotation_export(
        record, alignment_root=alignment_root, bundle_root=bundle_root
    )

    assert plan.transform.displacement is None
    assert plan.transform_evidence["has_displacement"] is False
    assert plan.transform_evidence["ignored_unreferenced_displacement_file"] == str(
        stale_displacement
    )


def test_registration_only_refuses_unreferenced_displacement_file(tmp_path, monkeypatch):
    record, _, alignment_root, bundle_root, _ = _new_bundle(tmp_path, monkeypatch)
    stale_displacement = alignment_root / f"donor_id={record.donor_id}" / "refined_displacement.npz"
    stale_displacement.write_bytes(b"unreferenced artifact")

    with pytest.raises(
        AX.AnnotationAugmentationError,
        match="displacement declaration/file presence disagree",
    ):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)


def test_tampered_legacy_transform_is_refused_even_for_protected_donor(tmp_path, monkeypatch):
    record, _, alignment_root, bundle_root, destination, _ = _legacy_bundle(tmp_path, monkeypatch)
    manifest_path = destination / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["transform"]["affine"][0][2] += 1.0
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(AX.AnnotationAugmentationError, match="does not exactly match"):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)


def test_gate_free_legacy_exception_does_not_extend_to_other_donors(tmp_path, monkeypatch):
    record, _, alignment_root, bundle_root, _, _ = _legacy_bundle(
        tmp_path, monkeypatch, donor="6539", serial="0059865"
    )
    with pytest.raises(AX.AnnotationAugmentationError, match="only for protected donors"):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("donor_id", "wrong", "donor_id"),
        ("xenium_serial", "wrong", "xenium_serial"),
        ("roi", "wrong", "roi"),
        ("xenium_bundle", "/wrong/bundle", "xenium_bundle"),
    ],
)
def test_bundle_manifest_identity_mismatch_is_refused(tmp_path, monkeypatch, field, value, message):
    record, _, alignment_root, bundle_root, destination = _new_bundle(tmp_path, monkeypatch)
    manifest_path = destination / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest[field] = value
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(AX.AnnotationAugmentationError, match=message):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)


def test_registration_gate_hash_change_is_refused(tmp_path, monkeypatch):
    record, _, alignment_root, bundle_root, destination = _new_bundle(tmp_path, monkeypatch)
    manifest = json.loads((destination / "manifest.json").read_text())
    pathlib.Path(manifest["alignment"]["gate_file"]).write_text("{}")
    with pytest.raises(AX.AnnotationAugmentationError, match="gate_sha256"):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)


def test_registration_only_bundle_without_exact_pass_gate_is_refused(tmp_path, monkeypatch):
    record, _, alignment_root, bundle_root, destination = _new_bundle(tmp_path, monkeypatch)
    manifest_path = destination / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for key in ("gate_file", "gate_sha256", "verdict"):
        manifest["alignment"].pop(key)
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(AX.AnnotationAugmentationError, match="must declare its exact PASS gate"):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)


@pytest.mark.parametrize("legacy", [False, True])
def test_manifest_if_pixel_calibration_must_match_qptiff_xml(tmp_path, monkeypatch, legacy):
    if legacy:
        record, _, alignment_root, bundle_root, destination, _ = _legacy_bundle(
            tmp_path, monkeypatch
        )
        manifest_path = destination / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["image"]["levels"][0]["um_per_px"] = 0.5
    else:
        record, _, alignment_root, bundle_root, destination = _new_bundle(tmp_path, monkeypatch)
        manifest_path = destination / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["input"]["pixel_size_um"] = 0.5
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(AX.AnnotationAugmentationError, match="disagrees with qptiff XML"):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)


def test_registered_image_must_still_equal_morphology_grid(tmp_path, monkeypatch):
    record, _, alignment_root, bundle_root, destination = _new_bundle(tmp_path, monkeypatch)
    image = next((destination / "images").iterdir())
    _registered_image(image, shape=(99, 100))
    with pytest.raises(AX.AnnotationAugmentationError, match="do not equal morphology"):
        AX.prepare_annotation_export(record, alignment_root=alignment_root, bundle_root=bundle_root)


def test_publication_failure_rolls_back_manifest_and_removes_new_annotations(tmp_path, monkeypatch):
    record, _, alignment_root, bundle_root, destination = _new_bundle(tmp_path, monkeypatch)
    originals = {
        name: (destination / name).read_bytes()
        for name in ("manifest.json", "README.md", "qc.json")
    }
    image = next((destination / "images").iterdir())
    image_before = (image.read_bytes(), image.stat().st_mtime_ns)
    plan = AX.prepare_annotation_export(
        record, alignment_root=alignment_root, bundle_root=bundle_root
    )
    real_move = AX._move
    failed = False

    def fail_installing_manifest(source, target):
        nonlocal failed
        if not failed and source.name == "manifest.json" and source.parent.name.endswith("partial"):
            failed = True
            raise OSError("injected publication failure")
        real_move(source, target)

    monkeypatch.setattr(AX, "_move", fail_installing_manifest)
    with pytest.raises(AX.AnnotationAugmentationError, match="rolled back"):
        AX.export_one(plan)

    assert not (destination / "annotations").exists()
    for name, content in originals.items():
        assert (destination / name).read_bytes() == content
    assert image.read_bytes() == image_before[0]
    assert image.stat().st_mtime_ns == image_before[1]
    assert not (destination / ".islet_annotations.partial").exists()
    assert not (destination / ".islet_annotations.previous").exists()


def test_explicit_missing_geojson_and_missing_existing_bundle_are_actionable(tmp_path):
    record, _ = _record(tmp_path)
    no_geojson = RI.RegistrationInput(**{**record.__dict__, "geojson": None})
    with pytest.raises(AX.AnnotationAugmentationError, match="no source islet GeoJSON"):
        AX.prepare_annotation_export(no_geojson, bundle_root=tmp_path / "bundles")
    with pytest.raises(AX.AnnotationAugmentationError, match="bundle is absent"):
        AX.prepare_annotation_export(record, bundle_root=tmp_path / "bundles")


def test_main_normalizes_inventory_and_annotation_errors_to_exit_two(monkeypatch, capsys):
    monkeypatch.setattr(
        AX.RI,
        "load_registration_inputs",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            RI.RegistrationInputError("synthetic inventory refusal")
        ),
    )
    assert AX.main([]) == 2
    assert "[islet-annotations] ERROR: synthetic inventory refusal" in capsys.readouterr().err
