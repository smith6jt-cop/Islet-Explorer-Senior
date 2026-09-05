#!/usr/bin/env python
"""Add transformed QuPath islet annotations to an existing Explorer bundle.

This is an annotation-only augmentation lane.  It never rewrites the registered OME-TIFF,
alignment files, keypoints, or cell groups.  The existing bundle establishes the exact
registration transform and Xenium morphology grid; a same-section QuPath GeoJSON is converted
through that transform and published as a Xenium Explorer 4.1 annotation sidecar.

Two evidence paths are accepted:

* a registration-only manifest whose declared transform is still certified by its exact PASS
  gate and content hashes; or
* the protected reviewed legacy bundles for donors 6436 and 6442, which predate that manifest
  schema.  They are accepted only when the complete transform embedded in the manifest equals
  the authoritative serialized transform and the registered image still exactly matches the
  declared Xenium morphology grid.

Publication stages a complete replacement ``annotations`` directory and metadata files inside
the donor bundle, then swaps only those paths with rollback.  Existing annotation files other
than the generated islet sidecar are copied forward byte-for-byte when
``--overwrite-annotations`` is used.

Run from ``core/phenocycler``::

    ../../.venv/bin/python scripts/restain/run_islet_annotation_export.py --donors 6436 6442
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import explorer as EX  # noqa: E402
import islet_annotations as IA  # noqa: E402
import registration_inputs as RI  # noqa: E402

from phenocycler.integration.transform import Transform  # noqa: E402

OUT = RI.REPO / "data" / "development" / "restain"
ALIGNMENT_ROOT = OUT / "valis"
BUNDLE_ROOT = OUT / "explorer_bundle"
PROTECTED_LEGACY_DONORS = frozenset({"6436", "6442"})

ANNOTATION_FILENAME = "{serial}_islets_on_xenium_grid.geojson"
ANNOTATION_QC_FILENAME = "islet_annotations_qc.json"
README_BEGIN = "<!-- BEGIN ISLET ANNOTATION AUGMENTATION -->"
README_END = "<!-- END ISLET ANNOTATION AUGMENTATION -->"
LEGACY_NO_ANNOTATION_TEXT = "* **No QuPath annotation of any class**, on either modality."
LEGACY_IMAGE_ANNOTATION_TEXT = (
    "* **No QuPath annotation was used to generate the original images or cell groups.**"
)


class AnnotationAugmentationError(RuntimeError):
    """An existing Explorer bundle cannot safely receive the annotation sidecar."""


@dataclass(frozen=True)
class AnnotationPlan:
    record: RI.RegistrationInput
    destination: Path
    manifest: dict[str, Any]
    morphology: dict[str, Any]
    transform: Transform
    transform_evidence: dict[str, Any]
    image_identity: dict[str, Any]
    input_calibration: dict[str, Any]
    overwrite_annotations: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json_object(path: Path, what: str) -> dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise AnnotationAugmentationError(f"cannot read {what} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise AnnotationAugmentationError(f"{what} {path} must contain one JSON object")
    return value


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _same_path(left: Path, right: Path) -> bool:
    return Path(left).expanduser().resolve(strict=False) == Path(right).expanduser().resolve(
        strict=False
    )


def _declared_path(value: Any, *, relative_to: Path, what: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise AnnotationAugmentationError(f"{what} is absent")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = Path(relative_to) / path
    return path.resolve(strict=False)


def _safe_bundle_file(destination: Path, value: Any, *, what: str) -> tuple[Path, Path]:
    if not isinstance(value, str) or not value.strip():
        raise AnnotationAugmentationError(f"{what} is absent")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise AnnotationAugmentationError(f"{what} must be a safe bundle-relative path: {value!r}")
    absolute = destination / relative
    if not absolute.is_file():
        raise AnnotationAugmentationError(f"{what} does not exist: {absolute}")
    return relative, absolute


def _finite_affine(value: Any, *, what: str) -> np.ndarray:
    try:
        affine = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise AnnotationAugmentationError(f"{what} must be one finite 2x3 affine") from exc
    if affine.shape != (2, 3) or not np.isfinite(affine).all():
        raise AnnotationAugmentationError(f"{what} must be one finite 2x3 affine")
    return affine


def _compare_manifest_transform(
    section: Mapping[str, Any], transform: Transform, *, what: str
) -> None:
    affine = _finite_affine(section.get("affine"), what=f"{what}.affine")
    if not np.array_equal(affine, transform.affine):
        raise AnnotationAugmentationError(
            f"{what}.affine does not exactly match the authoritative serialized transform"
        )
    has_displacement = section.get("has_displacement")
    if not isinstance(has_displacement, bool):
        raise AnnotationAugmentationError(f"{what}.has_displacement must be a JSON boolean")
    if has_displacement != (transform.displacement is not None):
        raise AnnotationAugmentationError(
            f"{what}.has_displacement disagrees with the authoritative serialized transform"
        )


def _serialized_transform(
    record: RI.RegistrationInput,
    stem: str,
    alignment_root: Path,
    *,
    allow_unreferenced_displacement: bool = False,
) -> tuple[Transform, Path, str, Path | None, str | None, dict[str, Any]]:
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", stem) is None:
        raise AnnotationAugmentationError(f"invalid transform stem declared by manifest: {stem!r}")
    donor_dir = Path(alignment_root) / f"donor_id={record.donor_id}"
    transform_path = donor_dir / f"{stem}.json"
    displacement_path = donor_dir / f"{stem}_displacement.npz"
    if not transform_path.is_file():
        raise AnnotationAugmentationError(
            f"authoritative serialized transform is absent: {transform_path}"
        )
    serialized = _read_json_object(transform_path, "serialized transform")
    has_displacement = serialized.get("has_displacement")
    if not isinstance(has_displacement, bool):
        raise AnnotationAugmentationError(
            f"{transform_path}: has_displacement must be a JSON boolean"
        )
    if has_displacement and not displacement_path.is_file():
        raise AnnotationAugmentationError(
            f"{transform_path}: displacement declaration/file presence disagree"
        )
    if not has_displacement and displacement_path.is_file() and not allow_unreferenced_displacement:
        raise AnnotationAugmentationError(
            f"{transform_path}: displacement declaration/file presence disagree"
        )
    try:
        transform = Transform.load(donor_dir, stem=stem)
    except Exception as exc:  # noqa: BLE001 - normalize malformed serialization
        raise AnnotationAugmentationError(f"cannot load {transform_path}: {exc}") from exc
    return (
        transform,
        transform_path,
        _sha256(transform_path),
        displacement_path if has_displacement else None,
        _sha256(displacement_path) if has_displacement else None,
        serialized,
    )


def _verify_gate(
    record: RI.RegistrationInput,
    stem: str,
    section: Mapping[str, Any],
    *,
    transform_path: Path,
    transform_sha256: str,
    displacement_path: Path | None,
    displacement_sha256: str | None,
) -> dict[str, Any] | None:
    missing_gate_fields = [
        key for key in ("gate_file", "gate_sha256", "verdict") if key not in section
    ]
    if missing_gate_fields:
        raise AnnotationAugmentationError(
            "registration-only manifest must declare its exact PASS gate; missing "
            f"{missing_gate_fields}"
        )
    if section.get("verdict") != "PASS":
        raise AnnotationAugmentationError(
            f"registration-only manifest verdict is {section.get('verdict')!r}, not 'PASS'"
        )
    expected_gate = transform_path.parent / f"alignment_{stem}.json"
    gate_path = _declared_path(
        section.get("gate_file"), relative_to=expected_gate.parent, what="alignment.gate_file"
    )
    if not _same_path(gate_path, expected_gate) or not gate_path.is_file():
        raise AnnotationAugmentationError(
            f"alignment.gate_file is not the authoritative sibling gate {expected_gate}"
        )
    gate_sha256 = _sha256(gate_path)
    if section.get("gate_sha256") != gate_sha256:
        raise AnnotationAugmentationError("alignment.gate_sha256 does not match the gate bytes")
    gate = _read_json_object(gate_path, "alignment gate")
    for key, expected in (
        ("donor_id", record.donor_id),
        ("serial", record.serial),
        ("roi", record.roi),
    ):
        if str(gate.get(key, "")) != expected:
            raise AnnotationAugmentationError(
                f"{gate_path}: {key}={gate.get(key)!r}, expected {expected!r}"
            )
    registration_input = gate.get("registration_input")
    if not isinstance(registration_input, Mapping):
        raise AnnotationAugmentationError(f"{gate_path}: registration_input is absent")
    for key, expected in (("qptiff", record.qptiff), ("bundle", record.bundle)):
        declared = registration_input.get(key)
        if not declared or not _same_path(Path(str(declared)), expected):
            raise AnnotationAugmentationError(
                f"{gate_path}: registration_input.{key} does not match inventory"
            )
    audited = gate.get("alignment")
    if not isinstance(audited, Mapping):
        raise AnnotationAugmentationError(f"{gate_path}: alignment object is absent")
    if audited.get("transform_source") != stem or audited.get("verdict") != "PASS":
        raise AnnotationAugmentationError(
            f"{gate_path}: gate does not certify PASS for transform {stem!r}"
        )
    audited_path = _declared_path(
        audited.get("transform_file"), relative_to=gate_path.parent, what="gate transform_file"
    )
    if not _same_path(audited_path, transform_path):
        raise AnnotationAugmentationError(f"{gate_path}: transform_file is not {transform_path}")
    if audited.get("transform_sha256") != transform_sha256:
        raise AnnotationAugmentationError(f"{gate_path}: transform hash is stale")
    if audited.get("displacement_sha256") != displacement_sha256:
        raise AnnotationAugmentationError(f"{gate_path}: displacement hash is stale")
    expected_has_displacement = displacement_path is not None
    if audited.get("has_displacement") is not expected_has_displacement:
        raise AnnotationAugmentationError(f"{gate_path}: has_displacement is stale")
    if not np.array_equal(
        _finite_affine(audited.get("affine"), what=f"{gate_path}: alignment.affine"),
        _finite_affine(section.get("affine"), what="manifest alignment.affine"),
    ):
        raise AnnotationAugmentationError(f"{gate_path}: affine disagrees with bundle manifest")
    return {"file": str(gate_path), "sha256": gate_sha256, "verdict": "PASS"}


def _registration_transform_evidence(
    record: RI.RegistrationInput,
    manifest: Mapping[str, Any],
    alignment_root: Path,
) -> tuple[Transform, dict[str, Any]]:
    section = manifest.get("alignment")
    if not isinstance(section, Mapping):
        raise AnnotationAugmentationError("registration-only manifest has no alignment object")
    stem = section.get("transform_source")
    if not isinstance(stem, str):
        raise AnnotationAugmentationError("alignment.transform_source is absent")
    transform, path, digest, displacement, displacement_digest, _ = _serialized_transform(
        record, stem, alignment_root
    )
    _compare_manifest_transform(section, transform, what="alignment")

    if "transform_file" in section:
        declared = _declared_path(
            section.get("transform_file"), relative_to=path.parent, what="alignment.transform_file"
        )
        if not _same_path(declared, path):
            raise AnnotationAugmentationError(
                f"alignment.transform_file is not the authoritative sibling {path}"
            )
    if "transform_sha256" in section and section.get("transform_sha256") != digest:
        raise AnnotationAugmentationError(
            "alignment.transform_sha256 does not match transform bytes"
        )
    if "displacement_file" in section:
        declared_displacement = section.get("displacement_file")
        if displacement is None:
            if declared_displacement is not None:
                raise AnnotationAugmentationError(
                    "alignment.displacement_file must be null without a displacement"
                )
        else:
            declared = _declared_path(
                declared_displacement,
                relative_to=path.parent,
                what="alignment.displacement_file",
            )
            if not _same_path(declared, displacement):
                raise AnnotationAugmentationError(
                    f"alignment.displacement_file is not the authoritative sibling {displacement}"
                )
    if (
        "displacement_sha256" in section
        and section.get("displacement_sha256") != displacement_digest
    ):
        raise AnnotationAugmentationError(
            "alignment.displacement_sha256 does not match displacement bytes"
        )
    gate = _verify_gate(
        record,
        stem,
        section,
        transform_path=path,
        transform_sha256=digest,
        displacement_path=displacement,
        displacement_sha256=displacement_digest,
    )
    policy = {
        "name": "registration_only_exact_manifest_transform",
        "scope": "registration geometry only; no cell identity or marker positivity",
        "protected_legacy_exception": False,
        "pass_gate_revalidated": True,
        "requirements": [
            "bundle donor/serial/roi and source paths equal registration inventory",
            "registered image equals the declared Xenium morphology grid",
            "manifest affine and displacement presence equal serialized Transform",
            "PASS gate, transform, and displacement content hashes revalidated",
        ],
    }
    evidence = {
        "manifest_schema": "registration_only",
        "stem": stem,
        "file": str(path),
        "sha256": digest,
        "displacement_file": str(displacement) if displacement is not None else None,
        "displacement_sha256": displacement_digest,
        "has_displacement": displacement is not None,
        "gate": gate,
        "evidence_policy": policy,
    }
    return transform, evidence


def _legacy_transform_evidence(
    record: RI.RegistrationInput,
    manifest: Mapping[str, Any],
    alignment_root: Path,
) -> tuple[Transform, dict[str, Any]]:
    if record.donor_id not in PROTECTED_LEGACY_DONORS:
        raise AnnotationAugmentationError(
            "a gate-free legacy manifest is allowed only for protected donors 6436 and 6442"
        )
    section = manifest.get("transform")
    if not isinstance(section, Mapping):
        raise AnnotationAugmentationError("protected legacy manifest has no transform object")
    stem = section.get("source")
    if not isinstance(stem, str):
        raise AnnotationAugmentationError("legacy transform.source is absent")
    transform, path, digest, displacement, displacement_digest, _ = _serialized_transform(
        record, stem, alignment_root, allow_unreferenced_displacement=True
    )
    _compare_manifest_transform(section, transform, what="transform")
    serialized_view = transform.to_dict()
    mismatched = [key for key, value in serialized_view.items() if section.get(key) != value]
    if mismatched:
        raise AnnotationAugmentationError(
            "protected legacy manifest transform does not exactly equal serialized Transform "
            f"for fields {mismatched}"
        )
    policy = {
        "name": "protected_legacy_exact_serialized_transform",
        "scope": "registration geometry only; no cell identity or marker positivity",
        "protected_legacy_exception": True,
        "pass_gate_revalidated": False,
        "requirements": [
            "donor is protected legacy donor 6436 or 6442",
            "complete manifest transform equals authoritative serialized Transform",
            (
                "an unreferenced displacement file is ignored when both serialized and manifest "
                "transforms declare has_displacement=false"
            ),
            "registered image equals the declared Xenium morphology grid",
            "bundle donor/serial/roi and source paths equal registration inventory",
        ],
    }
    evidence = {
        "manifest_schema": "protected_legacy",
        "stem": stem,
        "file": str(path),
        "sha256": digest,
        "displacement_file": str(displacement) if displacement is not None else None,
        "displacement_sha256": displacement_digest,
        "ignored_unreferenced_displacement_file": (
            str(path.parent / f"{stem}_displacement.npz")
            if displacement is None and (path.parent / f"{stem}_displacement.npz").is_file()
            else None
        ),
        "has_displacement": displacement is not None,
        "gate": None,
        "evidence_policy": policy,
    }
    return transform, evidence


def _morphology_identity(
    record: RI.RegistrationInput,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    declared = manifest.get("morphology")
    if not isinstance(declared, Mapping):
        raise AnnotationAugmentationError("bundle manifest has no morphology object")
    try:
        actual = EX.morphology_grid(record.bundle)
    except Exception as exc:  # noqa: BLE001 - add donor/bundle context
        raise AnnotationAugmentationError(
            f"cannot inspect Xenium morphology for donor {record.donor_id}: {exc}"
        ) from exc
    for key in ("height", "width", "n_levels"):
        if int(declared.get(key, -1)) != int(actual[key]):
            raise AnnotationAugmentationError(
                f"manifest morphology.{key}={declared.get(key)!r}, actual is {actual[key]!r}"
            )
    try:
        declared_um = float(declared.get("um_per_px"))
    except (TypeError, ValueError) as exc:
        raise AnnotationAugmentationError("manifest morphology.um_per_px is invalid") from exc
    if declared_um != float(actual["um_per_px"]):
        raise AnnotationAugmentationError(
            f"manifest morphology.um_per_px={declared_um}, actual is {actual['um_per_px']}"
        )
    declared_shapes = [list(map(int, shape)) for shape in declared.get("level_shapes", [])]
    actual_shapes = [list(map(int, shape)) for shape in actual["level_shapes"]]
    if declared_shapes != actual_shapes:
        raise AnnotationAugmentationError("manifest morphology.level_shapes are stale")
    declared_file = declared.get("file")
    if declared_file and not _same_path(Path(str(declared_file)), Path(actual["file"])):
        raise AnnotationAugmentationError("manifest morphology.file does not match Xenium bundle")
    return {
        "file": str(actual["file"]),
        "height": int(actual["height"]),
        "width": int(actual["width"]),
        "n_levels": int(actual["n_levels"]),
        "level_shapes": actual_shapes,
        "um_per_px": float(actual["um_per_px"]),
    }


def _input_pixel_calibration(
    record: RI.RegistrationInput,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    """Bind QuPath pixel-to-micron scaling to calibration already recorded by the bundle."""
    recorded: list[tuple[str, Any]] = []
    if manifest.get("mode") == "registration_only":
        input_section = manifest.get("input")
        if isinstance(input_section, Mapping):
            recorded.append(("input.pixel_size_um", input_section.get("pixel_size_um")))
    else:
        image = manifest.get("image")
        levels = image.get("levels") if isinstance(image, Mapping) else None
        if isinstance(levels, list) and levels and isinstance(levels[0], Mapping):
            recorded.append(("image.levels[0].um_per_px", levels[0].get("um_per_px")))
        alignment = manifest.get("alignment")
        if isinstance(alignment, Mapping) and "if_um_per_px" in alignment:
            recorded.append(("alignment.if_um_per_px", alignment.get("if_um_per_px")))
    if not recorded:
        raise AnnotationAugmentationError(
            "bundle manifest has no IF level-0 pixel calibration to bind GeoJSON coordinates"
        )
    values: dict[str, float] = {}
    for name, raw in recorded:
        try:
            value = float(raw)
        except (TypeError, ValueError) as exc:
            raise AnnotationAugmentationError(f"manifest {name} is not a valid pixel size") from exc
        if not np.isfinite(value) or value <= 0:
            raise AnnotationAugmentationError(f"manifest {name} is not a positive finite value")
        if value != float(record.pixel_um):
            raise AnnotationAugmentationError(
                f"manifest {name}={value} disagrees with qptiff XML pixel size {record.pixel_um}"
            )
        values[name] = value
    return {
        "qptiff_xml_pixel_size_um": float(record.pixel_um),
        "matching_manifest_fields": values,
        "status": "MATCH",
    }


def _ome_physical_size(ome_xml: str, path: Path) -> tuple[float, float]:
    try:
        root = ET.fromstring(ome_xml)
    except ET.ParseError as exc:
        raise AnnotationAugmentationError(f"registered image {path} has malformed OME-XML") from exc
    pixels = next((node for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "Pixels"), None)
    if pixels is None:
        raise AnnotationAugmentationError(f"registered image {path} has no OME Pixels element")
    try:
        x = float(pixels.attrib["PhysicalSizeX"])
        y = float(pixels.attrib["PhysicalSizeY"])
    except (KeyError, TypeError, ValueError) as exc:
        raise AnnotationAugmentationError(
            f"registered image {path} has no valid physical pixel size"
        ) from exc
    return x, y


def _registered_image_identity(
    destination: Path,
    manifest: Mapping[str, Any],
    morphology: Mapping[str, Any],
) -> dict[str, Any]:
    import tifffile

    image_section = manifest.get("image_on_xenium_grid")
    if not isinstance(image_section, Mapping):
        raise AnnotationAugmentationError("bundle manifest has no image_on_xenium_grid object")
    relative, image_path = _safe_bundle_file(
        destination, image_section.get("file"), what="image_on_xenium_grid.file"
    )
    expected_shapes = [tuple(shape) for shape in morphology["level_shapes"]]
    try:
        with tifffile.TiffFile(str(image_path)) as tif:
            if not tif.is_ome or len(tif.series) != 1:
                raise AnnotationAugmentationError(
                    f"registered image {image_path} must be one OME-TIFF series"
                )
            series = tif.series[0]
            actual_shapes = [tuple(map(int, level.shape[-2:])) for level in series.levels]
            if actual_shapes != expected_shapes:
                raise AnnotationAugmentationError(
                    f"registered image levels {actual_shapes} do not equal morphology "
                    f"levels {expected_shapes}"
                )
            physical_x, physical_y = _ome_physical_size(tif.ome_metadata or "", image_path)
    except AnnotationAugmentationError:
        raise
    except Exception as exc:  # noqa: BLE001 - normalize TIFF/codec errors
        raise AnnotationAugmentationError(
            f"cannot inspect registered image {image_path}: {exc}"
        ) from exc
    expected_um = float(morphology["um_per_px"])
    if physical_x != expected_um or physical_y != expected_um:
        raise AnnotationAugmentationError(
            f"registered image physical size {(physical_x, physical_y)} does not equal "
            f"morphology {expected_um}"
        )
    for key, expected in (
        ("width", morphology["width"]),
        ("height", morphology["height"]),
        ("n_levels", morphology["n_levels"]),
    ):
        if key in image_section and int(image_section[key]) != int(expected):
            raise AnnotationAugmentationError(
                f"image_on_xenium_grid.{key} disagrees with morphology"
            )
    if "level_shapes" in image_section:
        recorded_shapes = [list(map(int, shape)) for shape in image_section["level_shapes"]]
        if recorded_shapes != morphology["level_shapes"]:
            raise AnnotationAugmentationError(
                "image_on_xenium_grid.level_shapes disagree with morphology"
            )
    if "um_per_px" in image_section and float(image_section["um_per_px"]) != expected_um:
        raise AnnotationAugmentationError(
            "image_on_xenium_grid.um_per_px disagrees with morphology"
        )
    stat = image_path.stat()
    return {
        "file": relative.as_posix(),
        "bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "level_shapes": [list(shape) for shape in actual_shapes],
        "um_per_px": expected_um,
        "identity_check": "OME dimensions, pyramid levels, and physical scale equal morphology",
    }


def prepare_annotation_export(
    record: RI.RegistrationInput,
    *,
    alignment_root: Path = ALIGNMENT_ROOT,
    bundle_root: Path = BUNDLE_ROOT,
    overwrite_annotations: bool = False,
) -> AnnotationPlan:
    """Validate all identity and registration evidence before staging any output."""
    if record.geojson is None or not Path(record.geojson).is_file():
        raise AnnotationAugmentationError(
            f"donor {record.donor_id} ({record.serial}) has no source islet GeoJSON"
        )
    if not Path(record.qptiff).is_file():
        raise AnnotationAugmentationError(f"qptiff is absent: {record.qptiff}")
    destination = Path(bundle_root) / f"donor_{record.donor_id}_{record.serial}"
    if not destination.is_dir():
        raise AnnotationAugmentationError(f"existing Explorer bundle is absent: {destination}")
    partial = destination / ".islet_annotations.partial"
    backup = destination / ".islet_annotations.previous"
    if partial.exists() or backup.exists():
        stale = partial if partial.exists() else backup
        raise AnnotationAugmentationError(
            f"stale annotation transaction path exists: {stale}; inspect it before retrying"
        )
    annotations = destination / "annotations"
    if annotations.exists() and not overwrite_annotations:
        raise AnnotationAugmentationError(
            f"{annotations} already exists; use --overwrite-annotations to replace only the "
            "generated islet sidecar while preserving its other files"
        )

    manifest_path = destination / "manifest.json"
    manifest = _read_json_object(manifest_path, "Explorer bundle manifest")
    for key, expected in (
        ("donor_id", record.donor_id),
        ("xenium_serial", record.serial),
        ("roi", record.roi),
    ):
        if str(manifest.get(key, "")) != expected:
            raise AnnotationAugmentationError(
                f"{manifest_path}: {key}={manifest.get(key)!r}, expected {expected!r}"
            )
    declared_bundle = manifest.get("xenium_bundle")
    if not declared_bundle or not _same_path(Path(str(declared_bundle)), record.bundle):
        raise AnnotationAugmentationError(
            f"{manifest_path}: xenium_bundle does not match registration inventory"
        )
    input_section = manifest.get("input")
    declared_qptiff = (
        input_section.get("qptiff")
        if isinstance(input_section, Mapping)
        else manifest.get("qptiff")
    )
    if not declared_qptiff or not _same_path(Path(str(declared_qptiff)), record.qptiff):
        raise AnnotationAugmentationError(
            f"{manifest_path}: qptiff does not match registration inventory"
        )

    morphology = _morphology_identity(record, manifest)
    image_identity = _registered_image_identity(destination, manifest, morphology)
    input_calibration = _input_pixel_calibration(record, manifest)
    if manifest.get("mode") == "registration_only":
        transform, evidence = _registration_transform_evidence(record, manifest, alignment_root)
    else:
        transform, evidence = _legacy_transform_evidence(record, manifest, alignment_root)
    return AnnotationPlan(
        record=record,
        destination=destination,
        manifest=manifest,
        morphology=morphology,
        transform=transform,
        transform_evidence=evidence,
        image_identity=image_identity,
        input_calibration=input_calibration,
        overwrite_annotations=bool(overwrite_annotations),
    )


def _readme_with_annotation(existing: str, annotation: Mapping[str, Any]) -> str:
    counts = annotation["counts"]
    section = f"""{README_BEGIN}
## Import islet annotations

In Xenium Explorer 4.1+, open the **Annotations** layer, choose **Import**, and select
`{annotation["file"]}`. Associate it with the **Xenium In Situ morphology image**. Coordinates
are already transformed into morphology level-0 pixels; choose **pixels** if a unit prompt appears,
and do not apply an image alignment file.

The source contains {counts["islet_source_features"]} Islet features. This sidecar retains
{counts["exported_islet_polygons"]} polygons at least
{annotation["minimum_equivalent_diameter_um"]:g} µm in equivalent diameter. It is registration
geometry only: it does not change or support any cell identity, marker role, or positivity call.
The existing images and cell-group files were not rewritten.

Validation and exact source/output/transform hashes are recorded in `{annotation["qc_file"]}`
and `manifest.json`.
{README_END}"""
    # The protected legacy README accurately said there was no QuPath product in the original
    # bundle.  Once a sidecar is added, retain the historical claim but narrow it to the payloads
    # whose provenance it described so the current README does not contradict itself.
    text = existing.replace(LEGACY_NO_ANNOTATION_TEXT, LEGACY_IMAGE_ANNOTATION_TEXT).rstrip()
    if README_BEGIN in text or README_END in text:
        if text.count(README_BEGIN) != 1 or text.count(README_END) != 1:
            raise AnnotationAugmentationError(
                "README has malformed annotation augmentation markers"
            )
        start = text.index(README_BEGIN)
        end = text.index(README_END, start) + len(README_END)
        text = text[:start].rstrip() + "\n\n" + section + text[end:]
    else:
        text = (text + "\n\n" if text else "") + section
    return text.rstrip() + "\n"


def _move(source: Path, destination: Path) -> None:
    """Same-filesystem atomic rename, split out so rollback is fault-injectable in tests."""
    os.replace(source, destination)


def _publish_transaction(
    destination: Path,
    partial: Path,
    names: Sequence[str],
) -> None:
    backup = destination / ".islet_annotations.previous"
    backup.mkdir()
    backed_up: list[str] = []
    installed: list[str] = []
    try:
        for name in names:
            source = partial / name
            target = destination / name
            if not source.exists():
                raise AnnotationAugmentationError(f"staged transaction is missing {source}")
            if target.exists():
                _move(target, backup / name)
                backed_up.append(name)
            _move(source, target)
            installed.append(name)
    except BaseException as exc:
        rollback_errors: list[str] = []
        for name in reversed(installed):
            target = destination / name
            if target.exists():
                try:
                    _move(target, partial / name)
                except BaseException as rollback_exc:  # noqa: BLE001
                    rollback_errors.append(f"remove new {name}: {rollback_exc}")
        for name in reversed(backed_up):
            old = backup / name
            target = destination / name
            if old.exists():
                try:
                    _move(old, target)
                except BaseException as rollback_exc:  # noqa: BLE001
                    rollback_errors.append(f"restore old {name}: {rollback_exc}")
        if rollback_errors:
            raise AnnotationAugmentationError(
                f"annotation publication failed ({exc}) and rollback was incomplete; preserve "
                f"{partial} and {backup}: {rollback_errors}"
            ) from exc
        if backup.exists():
            shutil.rmtree(backup)
        raise AnnotationAugmentationError(
            f"annotation publication failed and was rolled back: {exc}"
        ) from exc
    if partial.exists():
        shutil.rmtree(partial)
    shutil.rmtree(backup)


def export_one(plan: AnnotationPlan) -> dict[str, Any]:
    """Stage, validate, and transactionally add one annotation sidecar and its metadata."""
    destination = plan.destination
    partial = destination / ".islet_annotations.partial"
    annotations_source = destination / "annotations"
    partial.mkdir()
    try:
        staged_annotations = partial / "annotations"
        if annotations_source.exists():
            shutil.copytree(annotations_source, staged_annotations, copy_function=shutil.copy2)
        else:
            staged_annotations.mkdir()
        relative_annotation = Path("annotations") / ANNOTATION_FILENAME.format(
            serial=plan.record.serial
        )
        staged_annotation = partial / relative_annotation
        annotation = IA.export_islet_geojson(
            plan.record.geojson,
            staged_annotation,
            plan.transform,
            input_um_per_px=float(plan.record.pixel_um),
            output_um_per_px=float(plan.morphology["um_per_px"]),
            width=int(plan.morphology["width"]),
            height=int(plan.morphology["height"]),
        )
        # Reopen once more at the transaction boundary rather than trusting the in-memory record.
        validation = IA.validate_islet_geojson(
            staged_annotation,
            width=int(plan.morphology["width"]),
            height=int(plan.morphology["height"]),
            expected_features=annotation["counts"]["exported_islet_polygons"],
        )
        if validation["sha256"] != annotation["output_sha256"]:
            raise AnnotationAugmentationError("staged annotation changed after export validation")

        generated_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
        annotation = copy.deepcopy(annotation)
        annotation["file"] = relative_annotation.as_posix()
        annotation["generated_utc"] = generated_utc
        annotation["validation"] = validation
        annotation["transform"] = copy.deepcopy(plan.transform_evidence)
        annotation["alignment_gate_sha256"] = (
            plan.transform_evidence["gate"]["sha256"]
            if plan.transform_evidence["gate"] is not None
            else None
        )
        annotation["transform_sha256"] = plan.transform_evidence["sha256"]
        annotation["displacement_sha256"] = plan.transform_evidence["displacement_sha256"]
        annotation["evidence_policy"] = copy.deepcopy(plan.transform_evidence["evidence_policy"])
        annotation["registered_image_identity"] = copy.deepcopy(plan.image_identity)
        annotation["input_calibration"] = copy.deepcopy(plan.input_calibration)

        existing_qc = destination / "qc.json"
        qc_name = "qc.json" if existing_qc.is_file() else ANNOTATION_QC_FILENAME
        annotation["qc_file"] = qc_name
        annotation_qc = {
            "status": "PASS",
            "scope": "islet annotation augmentation only",
            "generated_utc": generated_utc,
            "donor_id": plan.record.donor_id,
            "xenium_serial": plan.record.serial,
            "roi": plan.record.roi,
            "source_file": str(plan.record.geojson),
            "source_sha256": annotation["source_sha256"],
            "output_file": relative_annotation.as_posix(),
            "output_sha256": annotation["output_sha256"],
            "transform": copy.deepcopy(plan.transform_evidence),
            "evidence_policy": copy.deepcopy(annotation["evidence_policy"]),
            "registered_image_identity": copy.deepcopy(plan.image_identity),
            "input_calibration": copy.deepcopy(plan.input_calibration),
            "output_validation": validation,
            "preservation": {
                "images_rewritten": False,
                "cell_groups_rewritten": False,
                "existing_non_target_annotations_copied_forward": annotations_source.exists(),
            },
        }

        manifest = copy.deepcopy(plan.manifest)
        manifest["islet_annotations"] = annotation
        manifest["uses_qupath_annotations"] = True
        _write_json(partial / "manifest.json", manifest)

        readme_path = destination / "README.md"
        try:
            readme = readme_path.read_text() if readme_path.is_file() else ""
        except (OSError, UnicodeError) as exc:
            raise AnnotationAugmentationError(f"cannot read {readme_path}: {exc}") from exc
        (partial / "README.md").write_text(_readme_with_annotation(readme, annotation))

        if existing_qc.is_file():
            qc = _read_json_object(existing_qc, "existing bundle QC")
            qc["islet_annotations"] = annotation_qc
        else:
            qc = annotation_qc
        _write_json(partial / qc_name, qc)

        names = ["annotations", "manifest.json", "README.md", qc_name]
        _publish_transaction(destination, partial, names)
    except BaseException:
        # A failed rollback deliberately leaves its transaction directories for manual recovery.
        if partial.exists() and not (destination / ".islet_annotations.previous").exists():
            shutil.rmtree(partial)
        raise

    print(
        f"[islet-annotations] donor {plan.record.donor_id}: "
        f"{annotation['counts']['exported_islet_polygons']} polygons -> "
        f"{destination / relative_annotation}",
        flush=True,
    )
    return {
        "status": "EXPORTED",
        "destination": str(destination),
        "annotation": annotation,
    }


def select_annotation_records(
    records: Mapping[str, RI.RegistrationInput],
    donors: Sequence[str] | None,
) -> list[RI.RegistrationInput]:
    selected = RI.select_registration_inputs(
        records, donors, context="post-Xenium islet annotation export"
    )
    missing = [record.donor_id for record in selected if record.geojson is None]
    if donors and missing:
        raise AnnotationAugmentationError(
            f"requested donors have no source islet GeoJSON: {sorted(missing)}"
        )
    return [record for record in selected if record.geojson is not None]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--roi", default="panc")
    parser.add_argument("--inventory", type=Path, default=RI.REGISTRATION_INPUTS)
    parser.add_argument("--images-dir", type=Path, default=RI.IMAGES)
    parser.add_argument("--geojson-dir", type=Path, default=RI.GEOJSON_DIR)
    parser.add_argument("--xenium-paths", type=Path, default=RI.XENIUM_PATHS)
    parser.add_argument("--xenium-root", type=Path, default=RI.XENIUM_ROOT)
    parser.add_argument("--alignment-root", type=Path, default=ALIGNMENT_ROOT)
    parser.add_argument("--bundle-root", type=Path, default=BUNDLE_ROOT)
    parser.add_argument(
        "--overwrite-annotations",
        action="store_true",
        help=(
            "replace only the generated islet sidecar/annotation metadata; preserve images, "
            "cell groups, and other annotation files"
        ),
    )
    args = parser.parse_args(argv)

    try:
        records = RI.load_registration_inputs(
            args.inventory,
            args.xenium_paths,
            args.images_dir,
            args.geojson_dir,
            args.xenium_root,
            roi=args.roi,
        )
        selected = select_annotation_records(records, args.donors)
        if not selected:
            raise AnnotationAugmentationError("no eligible inventory records have a source GeoJSON")

        # Preflight every donor before publishing the first one.
        plans = [
            prepare_annotation_export(
                record,
                alignment_root=args.alignment_root,
                bundle_root=args.bundle_root,
                overwrite_annotations=args.overwrite_annotations,
            )
            for record in selected
        ]
        for plan in plans:
            export_one(plan)
    except (RI.RegistrationInputError, IA.IsletAnnotationError, AnnotationAugmentationError) as exc:
        print(f"[islet-annotations] ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"[islet-annotations] complete: {len(plans)} donor(s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
