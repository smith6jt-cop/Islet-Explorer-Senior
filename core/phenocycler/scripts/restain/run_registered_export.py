#!/usr/bin/env python
"""Export raw post-Xenium IF channels already registered to the Xenium morphology grid.

This is deliberately a registration-only lane.  It reads donor/serial/ROI provenance from
``restain_registration_inputs.csv`` and raw calibration/channel names from each qptiff.  When a
same-section post-Xenium QuPath GeoJSON is present, its >=50 um ``Islet`` polygons are carried
through the same transform into a Xenium Explorer annotation sidecar.  It does not import the
antibody-role registry, make positivity calls, or create cell groups.

Only a donor whose per-donor ``alignment_<stem>.json`` gate says ``PASS`` is eligible.  The
authoritative sibling ``<stem>.json`` (and optional ``<stem>_displacement.npz``) supplies the
micron-space Transform.  ``--transform-stem`` selects one exact gate/transform pair; there is no
fallback to another transform.  All native qptiff channels are pulled through that transform
onto the exact level-0 Xenium morphology grid, then a floor-halved box-mean pyramid is written
at 0.2125 um/px for Xenium Explorer's "Yes, skip alignment" import route.

Publication is transactional: the complete donor bundle is built and validated in a sibling
``.partial`` directory before one rename makes it visible.  Existing bundles are refused by
default, may be validated and skipped with ``--skip-existing``, and are replaced only with an
explicit ``--overwrite``.  The two already reviewed bundles (6436 and 6442) are always protected.

Run from ``core/phenocycler``::

    ../../.venv/bin/python scripts/restain/run_registered_export.py --donors 115 6591
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import time
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import explorer as EX  # noqa: E402
import islet_annotations as IA  # noqa: E402
import registration_inputs as RI  # noqa: E402

from phenocycler.integration.transform import Transform  # noqa: E402
# The one file checksum (phenocycler/artifacts.py): every byte, hex digest. Bound to the
# local name because other modules and tests call it through this one.
from phenocycler.artifacts import sha256_file as _sha256  # noqa: E402

OUT = RI.REPO / "data" / "development" / "restain"
ALIGNMENT_ROOT = OUT / "valis"
BUNDLE_ROOT = OUT / "explorer_bundle"

# These are the two bundles already reviewed in Xenium Explorer.  This new lane must not
# replace them even if --overwrite is supplied: their biological exports have a different
# contract, and preserving the user's verified reference is more important than CLI symmetry.
PROTECTED_COMPLETED_DONORS = frozenset({"6436", "6442"})

MODE = "registration_only"
OUTPUT_UM_PER_PX = EX.XENIUM_UM_PER_PX
DEFAULT_TILE = EX.EXPLORER_TILE


class RegisteredExportError(RuntimeError):
    """A registration-only bundle cannot be truthfully or atomically exported."""


@dataclass(frozen=True)
class PassedAlignment:
    """A PASS gate paired with the authoritative serialized Transform."""

    stem: str
    gate_path: Path
    transform_path: Path
    displacement_path: Path | None
    gate: dict[str, Any]
    transform: Transform
    gate_sha256: str
    transform_sha256: str
    displacement_sha256: str | None


@dataclass(frozen=True)
class ExportPlan:
    """Everything checked before the first potentially expensive donor export starts."""

    record: RI.RegistrationInput
    alignment: PassedAlignment
    morphology: dict[str, Any]
    destination: Path
    action: str


def _read_json_object(path: Path, what: str) -> dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RegisteredExportError(f"cannot read {what} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RegisteredExportError(f"{what} {path} must contain one JSON object")
    return value


def _same_path(left: Path, right: Path) -> bool:
    """Compare paths without requiring either target to exist."""
    return Path(left).expanduser().resolve(strict=False) == Path(right).expanduser().resolve(
        strict=False
    )


def load_passed_alignment(
    record: RI.RegistrationInput,
    alignment_root: Path = ALIGNMENT_ROOT,
    *,
    transform_stem: str = "refined",
) -> PassedAlignment:
    """Require a donor-specific PASS gate and load its sibling authoritative Transform.

    ``alignment_<stem>.json`` is evidence that the mapping passed the held-out registration
    checks; it is not itself transform serialization.  Reconstructing an affine from the gate
    would silently drop a displacement field, so the transform always comes from the exact
    selected stem.  The gate's affine and content hashes are compared with that object as
    independent direction/unit/identity checks.
    """
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", str(transform_stem)) is None:
        raise RegisteredExportError(
            f"transform stem must be one filename token, got {transform_stem!r}"
        )
    transform_stem = str(transform_stem)
    donor_dir = Path(alignment_root) / f"donor_id={record.donor_id}"
    gate_path = donor_dir / f"alignment_{transform_stem}.json"
    transform_path = donor_dir / f"{transform_stem}.json"
    displacement_path = donor_dir / f"{transform_stem}_displacement.npz"
    if not gate_path.is_file():
        raise RegisteredExportError(
            f"donor {record.donor_id} ({record.serial}) has no PASS gate: {gate_path}"
        )
    if not transform_path.is_file():
        raise RegisteredExportError(
            f"donor {record.donor_id} ({record.serial}) has no authoritative transform: "
            f"{transform_path}"
        )

    gate = _read_json_object(gate_path, "alignment gate")
    expected = {
        "donor_id": record.donor_id,
        "serial": record.serial,
        "roi": record.roi,
    }
    for key, value in expected.items():
        if str(gate.get(key, "")) != value:
            raise RegisteredExportError(f"{gate_path}: {key}={gate.get(key)!r}, expected {value!r}")

    alignment = gate.get("alignment")
    if not isinstance(alignment, dict):
        raise RegisteredExportError(f"{gate_path}: missing object-valued 'alignment'")
    if alignment.get("transform_source") != transform_stem:
        raise RegisteredExportError(
            f"{gate_path}: transform_source={alignment.get('transform_source')!r}, expected "
            f"{transform_stem!r}"
        )
    if alignment.get("verdict") != "PASS":
        raise RegisteredExportError(
            f"{gate_path}: alignment verdict is {alignment.get('verdict')!r}, not 'PASS'"
        )

    declared_gate = alignment.get("gate_file")
    if not declared_gate:
        raise RegisteredExportError(f"{gate_path}: alignment.gate_file is absent")
    declared_gate_path = Path(str(declared_gate)).expanduser()
    if not declared_gate_path.is_absolute():
        declared_gate_path = gate_path.parent / declared_gate_path
    if not _same_path(declared_gate_path, gate_path):
        raise RegisteredExportError(
            f"{gate_path}: alignment.gate_file points to {declared_gate_path}"
        )

    registration_input = gate.get("registration_input")
    if not isinstance(registration_input, dict):
        raise RegisteredExportError(f"{gate_path}: missing object-valued registration_input")
    for key, expected_path in (("qptiff", record.qptiff), ("bundle", record.bundle)):
        declared = registration_input.get(key)
        if not declared or not _same_path(Path(str(declared)), expected_path):
            raise RegisteredExportError(
                f"{gate_path}: registration_input.{key}={declared!r}, expected {expected_path}"
            )
    certified = alignment.get("certifies_mapping")
    if not isinstance(certified, dict):
        raise RegisteredExportError(f"{gate_path}: missing alignment.certifies_mapping")
    certified_values = {
        "pixel_um_full_res": float(record.pixel_um),
        "source": "qptiff XML PixelSizeMicrons",
        "dapi_page": int(record.dapi_page),
        "dapi_channel_name": record.channel_names[record.dapi_page],
    }
    for key, expected_value in certified_values.items():
        if certified.get(key) != expected_value:
            raise RegisteredExportError(
                f"{gate_path}: certifies_mapping.{key}={certified.get(key)!r}, "
                f"expected {expected_value!r}"
            )

    declared_transform = alignment.get("transform_file")
    if not declared_transform:
        raise RegisteredExportError(f"{gate_path}: alignment.transform_file is absent")
    declared_path = Path(str(declared_transform)).expanduser()
    if not declared_path.is_absolute():
        declared_path = gate_path.parent / declared_path
    if not _same_path(declared_path, transform_path):
        raise RegisteredExportError(
            f"{gate_path}: declared transform {declared_path} is not sibling {transform_path}"
        )

    serialized = _read_json_object(transform_path, "serialized transform")
    if not isinstance(serialized.get("has_displacement"), bool):
        raise RegisteredExportError(f"{transform_path}: has_displacement must be a JSON boolean")
    serialized_has_displacement = serialized["has_displacement"]
    if serialized_has_displacement and not displacement_path.is_file():
        raise RegisteredExportError(
            f"{transform_path} declares a displacement field but {displacement_path} is absent"
        )
    if not serialized_has_displacement and displacement_path.exists():
        raise RegisteredExportError(
            f"{displacement_path} exists but {transform_path} declares no displacement field"
        )

    declared_displacement = alignment.get("displacement_file")
    if serialized_has_displacement:
        if not declared_displacement:
            raise RegisteredExportError(f"{gate_path}: alignment.displacement_file is absent")
        declared_displacement_path = Path(str(declared_displacement)).expanduser()
        if not declared_displacement_path.is_absolute():
            declared_displacement_path = gate_path.parent / declared_displacement_path
        if not _same_path(declared_displacement_path, displacement_path):
            raise RegisteredExportError(
                f"{gate_path}: declared displacement {declared_displacement_path} is not sibling "
                f"{displacement_path}"
            )
    elif declared_displacement is not None:
        raise RegisteredExportError(
            f"{gate_path}: alignment.displacement_file must be null without a displacement"
        )

    transform_sha256 = _sha256(transform_path)
    displacement_sha256 = _sha256(displacement_path) if serialized_has_displacement else None
    for hash_key in ("transform_sha256", "displacement_sha256"):
        if hash_key not in alignment:
            raise RegisteredExportError(f"{gate_path}: alignment.{hash_key} is absent")
    if alignment.get("transform_sha256") != transform_sha256:
        raise RegisteredExportError(
            f"{gate_path}: transform_sha256 does not match {transform_path}"
        )
    if alignment.get("displacement_sha256") != displacement_sha256:
        raise RegisteredExportError(
            f"{gate_path}: displacement_sha256 does not match {displacement_path}"
        )

    try:
        transform = Transform.load(donor_dir, stem=transform_stem)
    except Exception as exc:  # noqa: BLE001 - normalize all malformed serialization failures
        raise RegisteredExportError(
            f"cannot load authoritative transform {transform_path}: {exc}"
        ) from exc

    training_fold = transform.meta.get("training_fold")
    withheld_gate_fold = transform.meta.get("withheld_gate_fold")
    gate_fold = alignment.get("window_fold")
    if training_fold is None and withheld_gate_fold is None:
        if gate_fold != "all":
            raise RegisteredExportError(
                f"{gate_path}: non-CV transform requires alignment.window_fold='all', "
                f"found {gate_fold!r}"
            )
    else:
        valid_folds = {"even", "odd"}
        if training_fold not in valid_folds or withheld_gate_fold not in valid_folds:
            raise RegisteredExportError(
                f"{transform_path}: training_fold and withheld_gate_fold must both be "
                f"'even' or 'odd', found {training_fold!r}/{withheld_gate_fold!r}"
            )
        if training_fold == withheld_gate_fold:
            raise RegisteredExportError(
                f"{transform_path}: training_fold and withheld_gate_fold must differ"
            )
        if gate_fold != withheld_gate_fold:
            raise RegisteredExportError(
                f"{gate_path}: alignment.window_fold={gate_fold!r}, expected independently "
                f"withheld fold {withheld_gate_fold!r} (training fold {training_fold!r})"
            )

    try:
        affine = np.asarray(alignment.get("affine"), dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise RegisteredExportError(
            f"{gate_path}: alignment.affine must be one finite 2x3 array"
        ) from exc
    if affine.shape != (2, 3) or not np.isfinite(affine).all():
        raise RegisteredExportError(f"{gate_path}: alignment.affine must be one finite 2x3 array")
    if not np.allclose(affine, transform.affine, rtol=1e-12, atol=1e-12):
        raise RegisteredExportError(
            f"{gate_path}: audited affine differs from authoritative {transform_path}"
        )
    if (
        not np.isfinite(transform.affine).all()
        or abs(np.linalg.det(transform.affine[:, :2])) < 1e-12
    ):
        raise RegisteredExportError(f"{transform_path}: affine is non-finite or singular")

    declared_has_displacement = alignment.get("has_displacement")
    if not isinstance(declared_has_displacement, bool):
        raise RegisteredExportError(
            f"{gate_path}: alignment.has_displacement must be a JSON boolean"
        )
    if declared_has_displacement != serialized_has_displacement:
        raise RegisteredExportError(
            f"{gate_path}: has_displacement={declared_has_displacement!r} disagrees with "
            f"{transform_path}"
        )
    if transform.displacement is not None:
        if transform.grid is None:
            raise RegisteredExportError(f"{transform_path}: displacement has no grid")
        expected_shape = (2, transform.grid.height, transform.grid.width)
        if transform.displacement.shape != expected_shape:
            raise RegisteredExportError(
                f"{displacement_path}: shape {transform.displacement.shape}, expected {expected_shape}"
            )
        if not np.isfinite(transform.displacement).all():
            raise RegisteredExportError(
                f"{displacement_path}: displacement contains non-finite values"
            )

    return PassedAlignment(
        stem=transform_stem,
        gate_path=gate_path,
        transform_path=transform_path,
        displacement_path=displacement_path if serialized_has_displacement else None,
        gate=gate,
        transform=transform,
        gate_sha256=_sha256(gate_path),
        transform_sha256=transform_sha256,
        displacement_sha256=displacement_sha256,
    )


def _source_map_for_tile(
    transform: Transform,
    *,
    x0: int,
    y0: int,
    width: int,
    height: int,
    input_um_per_px: float,
    output_um_per_px: float = OUTPUT_UM_PER_PX,
) -> tuple[np.ndarray, np.ndarray]:
    """Map one output tile's pixel centres to level-0 qptiff pixel-centre coordinates.

    The Transform maps Xenium microns to IF microns.  Xenium output pixel ``(x, y)`` is centred
    at ``((x + .5) * 0.2125, (y + .5) * 0.2125)`` microns, while qptiff level-0 pixel zero is
    centred at ``0.5 * input_um_per_px``.  Subtracting that half-pixel before division makes an
    identity physical transform reproduce an input raster pixel-for-pixel when scales match.
    """
    yy, xx = np.mgrid[y0 : y0 + height, x0 : x0 + width]
    fixed_um = np.column_stack(
        [
            (xx.ravel() + 0.5) * float(output_um_per_px),
            (yy.ravel() + 0.5) * float(output_um_per_px),
        ]
    )
    source_um = transform.apply(fixed_um)
    source_px = source_um / float(input_um_per_px) - 0.5
    if not np.isfinite(source_px).all():
        raise RegisteredExportError("the transform produced non-finite qptiff sample coordinates")
    return (
        source_px[:, 0].reshape(height, width).astype(np.float32),
        source_px[:, 1].reshape(height, width).astype(np.float32),
    )


def _warp_store_to_xenium_grid(
    store,
    transform: Transform,
    *,
    input_um_per_px: float,
    height: int,
    width: int,
    tile: int = DEFAULT_TILE,
    progress_label: str | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Pull a ``CYX`` array onto the Xenium grid, reading only needed source crops.

    Mapping is calculated once per output tile and reused for every channel.  This keeps source
    I/O to qptiff-sized crops and, unlike an affine-only 3x3 shortcut, applies an optional
    displacement field through :meth:`Transform.apply`.
    """
    import cv2

    if len(store.shape) != 3:
        raise RegisteredExportError(f"qptiff level 0 must be CYX, found shape {store.shape}")
    if np.dtype(store.dtype) != np.dtype("uint8"):
        raise RegisteredExportError(f"qptiff level 0 must be uint8, found {store.dtype}")
    channels, source_h, source_w = (int(value) for value in store.shape)
    if channels < 1 or source_h < 1 or source_w < 1:
        raise RegisteredExportError(f"qptiff level 0 has invalid shape {store.shape}")
    if not float(input_um_per_px) > 0 or not float(OUTPUT_UM_PER_PX) > 0:
        raise RegisteredExportError("input and output pixel sizes must be positive")
    if int(tile) < 16:
        raise RegisteredExportError(f"tile must be at least 16, got {tile}")

    output = np.zeros((channels, int(height), int(width)), dtype=np.uint8)
    n_y = (int(height) + int(tile) - 1) // int(tile)
    n_x = (int(width) + int(tile) - 1) // int(tile)
    n_tiles = n_y * n_x
    announce_every = max(1, n_tiles // 20)
    outside = 0
    largest_crop = 0
    started = time.monotonic()

    for tile_index, (y0, x0) in enumerate(
        ((y, x) for y in range(0, int(height), int(tile)) for x in range(0, int(width), int(tile))),
        start=1,
    ):
        th = min(int(tile), int(height) - y0)
        tw = min(int(tile), int(width) - x0)
        map_x, map_y = _source_map_for_tile(
            transform, x0=x0, y0=y0, width=tw, height=th, input_um_per_px=input_um_per_px
        )
        in_field = (
            (map_x >= -0.5)
            & (map_x <= source_w - 0.5)
            & (map_y >= -0.5)
            & (map_y <= source_h - 0.5)
        )
        outside += int(in_field.size - np.count_nonzero(in_field))
        if not np.any(in_field):
            continue

        # Two pixels on each side cover bilinear neighbours and OpenCV's fixed-point rounding.
        cx0 = max(int(np.floor(float(map_x.min()))) - 2, 0)
        cy0 = max(int(np.floor(float(map_y.min()))) - 2, 0)
        cx1 = min(int(np.ceil(float(map_x.max()))) + 3, source_w)
        cy1 = min(int(np.ceil(float(map_y.max()))) + 3, source_h)
        if cx1 <= cx0 or cy1 <= cy0:
            continue
        crop_pixels = (cx1 - cx0) * (cy1 - cy0)
        largest_crop = max(largest_crop, crop_pixels)
        # A PASS same-section transform should map one tile to a similarly sized crop.  A huge
        # crop is evidence of a direction/unit error and can otherwise trigger an enormous read.
        if crop_pixels > 64 * int(tile) * int(tile):
            raise RegisteredExportError(
                f"output tile {(y0, x0)} maps to implausible source crop "
                f"{cy1 - cy0}x{cx1 - cx0}; transform direction or calibration is wrong"
            )
        if cx1 - cx0 >= 32767 or cy1 - cy0 >= 32767:
            raise RegisteredExportError("an OpenCV source crop exceeds its 32767-pixel limit")

        crop = np.asarray(store[:, cy0:cy1, cx0:cx1])
        if crop.shape != (channels, cy1 - cy0, cx1 - cx0) or crop.dtype != np.uint8:
            raise RegisteredExportError(
                f"qptiff crop has shape/dtype {crop.shape}/{crop.dtype}, expected "
                f"{(channels, cy1 - cy0, cx1 - cx0)}/uint8"
            )
        local_x = map_x - float(cx0)
        local_y = map_y - float(cy0)
        for channel in range(channels):
            output[channel, y0 : y0 + th, x0 : x0 + tw] = cv2.remap(
                np.ascontiguousarray(crop[channel]),
                local_x,
                local_y,
                interpolation=cv2.INTER_LINEAR,
                borderMode=cv2.BORDER_CONSTANT,
                borderValue=0,
            )

        if progress_label and (tile_index == n_tiles or tile_index % announce_every == 0):
            elapsed = time.monotonic() - started
            print(
                f"[registered] {progress_label}: warped {tile_index}/{n_tiles} tiles "
                f"({100.0 * tile_index / n_tiles:.0f}%, {elapsed:.0f} s)",
                flush=True,
            )

    elapsed = time.monotonic() - started
    return output, {
        "seconds": round(elapsed, 1),
        "outside_if_fraction": float(outside / (int(height) * int(width))),
        "largest_source_crop_pixels": int(largest_crop),
        "output_tile": int(tile),
        "interpolation": "OpenCV INTER_LINEAR",
        "border": "constant zero outside qptiff",
        "transform_application": "Xenium output pixel centre -> Transform.apply -> IF level-0 pixel centre",
        "includes_displacement": transform.displacement is not None,
    }


def warp_qptiff_to_xenium_level0(
    record: RI.RegistrationInput,
    transform: Transform,
    morphology: Mapping[str, Any],
    *,
    tile: int = DEFAULT_TILE,
    workers: int = 8,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Decode qptiff level 0 and warp every native channel in native order.

    ``tifffile``'s page decoder is used directly instead of its optional Zarr adapter.  The
    adapter's synchronous bridge can deadlock with Zarr 3 event loops; page decoding is the
    simpler contract here and the machine running this batch has room for one raw source stack
    plus one registered stack.  Channels are still decoded one at a time into a preallocated
    array, so there is no transient second copy of the five-channel qptiff.
    """
    import tifffile

    with tifffile.TiffFile(str(record.qptiff)) as tif:
        series = tif.series[0]
        if series.axes != "CYX":
            raise RegisteredExportError(
                f"{record.qptiff.name}: series 0 axes are {series.axes!r}, expected 'CYX'"
            )
        if np.dtype(series.dtype) != np.dtype("uint8"):
            raise RegisteredExportError(
                f"{record.qptiff.name}: series 0 is {series.dtype}, expected uint8"
            )
        level0 = series.levels[0]
        if len(level0.pages) != len(record.channel_names):
            raise RegisteredExportError(
                f"{record.qptiff.name}: {len(level0.pages)} pages but "
                f"{len(record.channel_names)} recorded names"
            )
        source_shape = (len(record.channel_names), int(level0.shape[-2]), int(level0.shape[-1]))
        source = np.empty(source_shape, dtype=np.uint8)
        for channel_index, page in enumerate(level0.pages):
            plane = page.asarray(maxworkers=int(workers))
            if plane.shape != source_shape[-2:] or plane.dtype != np.uint8:
                raise RegisteredExportError(
                    f"{record.qptiff.name}: channel {channel_index} decoded as "
                    f"{plane.shape}/{plane.dtype}, expected {source_shape[-2:]}/uint8"
                )
            source[channel_index] = plane
            del plane
            print(
                f"[registered] donor {record.donor_id}: decoded source channel "
                f"{channel_index + 1}/{len(record.channel_names)} "
                f"({record.channel_names[channel_index]})",
                flush=True,
            )
        result = _warp_store_to_xenium_grid(
            source,
            transform,
            input_um_per_px=record.pixel_um,
            height=int(morphology["height"]),
            width=int(morphology["width"]),
            tile=tile,
            progress_label=f"donor {record.donor_id}",
        )
        del source
    return result


def _ome_pixels(ome_xml: str, path: Path) -> tuple[ET.Element, list[str]]:
    try:
        root = ET.fromstring(ome_xml)
    except ET.ParseError as exc:
        raise RegisteredExportError(f"{path}: invalid OME-XML: {exc}") from exc
    pixels = root.find(".//{*}Pixels")
    if pixels is None:
        raise RegisteredExportError(f"{path}: OME-XML contains no Pixels element")
    names = [str(channel.get("Name", "")) for channel in pixels.findall("./{*}Channel")]
    return pixels, names


def _compression_names(compression: str) -> frozenset[str]:
    return {
        "zlib": frozenset({"ADOBE_DEFLATE", "DEFLATE"}),
        "lzw": frozenset({"LZW"}),
        "jpeg2000": frozenset({"JPEG2000", "JPEG2000_LOSSY"}),
        "none": frozenset({"NONE"}),
    }[compression]


def validate_registered_ome(
    path: Path,
    *,
    channel_names: Sequence[str],
    morphology: Mapping[str, Any],
    compression: str,
    tile: int,
) -> dict[str, Any]:
    """Reopen the written image and prove its Explorer-facing contract, including tile decode."""
    import tifffile

    path = Path(path)
    expected_channels = list(channel_names)
    expected_shapes = [tuple(int(v) for v in shape) for shape in morphology["level_shapes"]]
    expected_compression = _compression_names(compression)
    decoded: list[dict[str, Any]] = []

    try:
        with tifffile.TiffFile(str(path)) as tif:
            if not tif.is_ome or not tif.is_bigtiff:
                raise RegisteredExportError(f"{path}: output must be pyramidal OME BigTIFF")
            if len(tif.series) != 1:
                raise RegisteredExportError(
                    f"{path}: expected one TIFF series, found {len(tif.series)}"
                )
            series = tif.series[0]
            if series.axes != "CYX":
                raise RegisteredExportError(f"{path}: axes are {series.axes!r}, expected 'CYX'")
            expected_shape = (len(expected_channels), *expected_shapes[0])
            if tuple(series.shape) != expected_shape:
                raise RegisteredExportError(
                    f"{path}: shape {series.shape}, expected {expected_shape}"
                )
            if np.dtype(series.dtype) != np.dtype("uint8"):
                raise RegisteredExportError(f"{path}: dtype is {series.dtype}, expected uint8")
            if len(series.levels) != len(expected_shapes):
                raise RegisteredExportError(
                    f"{path}: {len(series.levels)} pyramid levels, expected {len(expected_shapes)}"
                )
            actual_shapes = [tuple(level.shape[-2:]) for level in series.levels]
            if actual_shapes != expected_shapes:
                raise RegisteredExportError(
                    f"{path}: pyramid shapes {actual_shapes}, expected {expected_shapes}"
                )

            ome_xml = tif.ome_metadata or ""
            pixels, actual_channels = _ome_pixels(ome_xml, path)
            if actual_channels != expected_channels:
                raise RegisteredExportError(
                    f"{path}: channels {actual_channels}, expected native order {expected_channels}"
                )
            for axis in ("X", "Y"):
                raw = pixels.get(f"PhysicalSize{axis}")
                unit = pixels.get(f"PhysicalSize{axis}Unit")
                try:
                    value = float(raw) if raw is not None else float("nan")
                except ValueError:
                    value = float("nan")
                if not np.isfinite(value) or abs(value - OUTPUT_UM_PER_PX) > 1e-12:
                    raise RegisteredExportError(
                        f"{path}: PhysicalSize{axis}={raw!r}, expected {OUTPUT_UM_PER_PX}"
                    )
                if unit != "µm":
                    raise RegisteredExportError(
                        f"{path}: PhysicalSize{axis}Unit={unit!r}, expected 'µm'"
                    )

            for level_index, level in enumerate(series.levels):
                if len(level.pages) != len(expected_channels):
                    raise RegisteredExportError(
                        f"{path}: level {level_index} has {len(level.pages)} pages, expected "
                        f"{len(expected_channels)}"
                    )
                for page in level.pages:
                    # tifffile represents later same-layout channels as TiffFrame objects;
                    # their keyframe owns the structural tile tags while offsets remain local.
                    layout = page.keyframe
                    if (
                        not layout.is_tiled
                        or layout.tilewidth != int(tile)
                        or layout.tilelength != int(tile)
                    ):
                        raise RegisteredExportError(
                            f"{path}: level {level_index} is not tiled {tile}x{tile}"
                        )
                    if page.compression.name not in expected_compression:
                        raise RegisteredExportError(
                            f"{path}: level {level_index} compression {page.compression.name}, "
                            f"expected {sorted(expected_compression)}"
                        )

                # Decode one actual compressed TIFF tile from every channel and level.  This
                # intentionally bypasses the optional Zarr adapter: validation is of the TIFF
                # bytes and installed codec, and a header-only check cannot establish that.
                height, width = expected_shapes[level_index]
                tiles_across = (width + int(tile) - 1) // int(tile)
                tile_row = (height // 2) // int(tile)
                tile_col = (width // 2) // int(tile)
                tile_index = tile_row * tiles_across + tile_col
                y0, x0 = tile_row * int(tile), tile_col * int(tile)
                y1, x1 = min(height, y0 + int(tile)), min(width, x0 + int(tile))
                for channel_index, (channel_name, page) in enumerate(
                    zip(expected_channels, level.pages, strict=True)
                ):
                    if tile_index >= len(page.dataoffsets):
                        raise RegisteredExportError(
                            f"{path}: level {level_index} channel {channel_index} lacks "
                            f"tile index {tile_index}"
                        )
                    byte_count = int(page.databytecounts[tile_index])
                    if byte_count <= 0:
                        raise RegisteredExportError(
                            f"{path}: level {level_index} channel {channel_index} tile "
                            f"{tile_index} has no encoded bytes"
                        )
                    tif.filehandle.seek(int(page.dataoffsets[tile_index]))
                    encoded = tif.filehandle.read(byte_count)
                    tile_array, _, _ = page.decode(
                        encoded, tile_index, jpegtables=page.keyframe.jpegtables
                    )
                    if tile_array is None:
                        raise RegisteredExportError(
                            f"{path}: codec returned no pixels for level {level_index} "
                            f"channel {channel_index} tile {tile_index}"
                        )
                    patch = np.squeeze(np.asarray(tile_array))
                    if patch.shape != (int(tile), int(tile)) or patch.dtype != np.uint8:
                        raise RegisteredExportError(
                            f"{path}: decoded level {level_index} channel {channel_index} "
                            f"tile as {patch.shape}/{patch.dtype}"
                        )
                    decoded.append(
                        {
                            "level": int(level_index),
                            "channel_index": int(channel_index),
                            "channel_name": channel_name,
                            "tile_index": int(tile_index),
                            "window_yx": [int(y0), int(y1), int(x0), int(x1)],
                            "encoded_bytes": byte_count,
                            "min": int(patch.min()),
                            "max": int(patch.max()),
                            "blake2b_128": hashlib.blake2b(
                                patch.tobytes(order="C"), digest_size=16
                            ).hexdigest(),
                        }
                    )
    except RegisteredExportError:
        raise
    except Exception as exc:  # noqa: BLE001 - codec/XML/TIFF errors need one donor-context type
        raise RegisteredExportError(f"cannot validate registered OME-TIFF {path}: {exc}") from exc

    return {
        "status": "PASS",
        "checks": {
            "ome": True,
            "bigtiff": True,
            "axes": "CYX",
            "dtype": "uint8",
            "shape": list(expected_shape),
            "channel_names_native_order": expected_channels,
            "level_shapes": [list(shape) for shape in expected_shapes],
            "n_levels": len(expected_shapes),
            "physical_size_um": OUTPUT_UM_PER_PX,
            "tile": int(tile),
            "compression": compression,
            "decoded_every_channel_every_level": True,
        },
        "decoded_tiles": decoded,
        "bytes": int(path.stat().st_size),
    }


def _alignment_manifest(alignment: PassedAlignment) -> dict[str, Any]:
    transform = alignment.transform
    return {
        "gate_file": str(alignment.gate_path),
        "gate_sha256": alignment.gate_sha256,
        "verdict": "PASS",
        "verdict_why": alignment.gate["alignment"].get("verdict_why"),
        "transform_source": alignment.stem,
        "transform_file": str(alignment.transform_path),
        "transform_sha256": alignment.transform_sha256,
        "displacement_file": (
            str(alignment.displacement_path) if alignment.displacement_path is not None else None
        ),
        "displacement_sha256": alignment.displacement_sha256,
        "direction": "Xenium microns -> IF microns",
        **transform.to_dict(),
    }


def _manifest(
    record: RI.RegistrationInput,
    alignment: PassedAlignment,
    morphology: Mapping[str, Any],
    image_record: Mapping[str, Any],
    warp_record: Mapping[str, Any],
    validation: Mapping[str, Any],
    image_name: str,
    annotation_record: Mapping[str, Any],
) -> dict[str, Any]:
    channels = [
        {"index": int(index), "name": name, "role": None, "role_status": "unassigned"}
        for index, name in enumerate(record.channel_names)
    ]
    return {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "mode": MODE,
        "donor_id": record.donor_id,
        "xenium_serial": record.serial,
        "roi": record.roi,
        "input": {
            "qptiff": str(record.qptiff),
            "qptiff_bytes": int(record.qptiff.stat().st_size),
            "qptiff_mtime_ns": int(record.qptiff.stat().st_mtime_ns),
            "pixel_size_um": float(record.pixel_um),
            "pixel_size_source": "qptiff XML PixelSizeMicrons",
            "dtype": "uint8",
            "pyramid_shapes": [list(shape) for shape in record.pyramid_shapes],
            "raw_channels": channels,
            "channel_order": "native qptiff series-0 page order",
            "dapi_page": int(record.dapi_page),
            "provenance_source": record.provenance_source,
        },
        "xenium_bundle": str(record.bundle),
        "morphology": dict(morphology),
        "alignment": _alignment_manifest(alignment),
        "image_on_xenium_grid": {
            "file": f"images/{image_name}",
            "import_route": "Images -> Add image -> Yes, skip alignment",
            "maps": "identity in Xenium morphology pixel space; never apply an alignment CSV",
            "um_per_px": OUTPUT_UM_PER_PX,
            "raw_channels_native_order": list(record.channel_names),
            "resampling": (
                "Each Xenium level-0 pixel centre is mapped Xenium-um -> IF-um by the full "
                "selected Transform and bilinearly sampled from qptiff level 0; outside IF is 0. "
                "Coarser levels are recursive 2x2 box means, floor-halved, round-half-up."
            ),
            "warp": dict(warp_record),
            "tiff": dict(image_record),
        },
        "channel_roles": {
            "status": "unassigned",
            "reason": (
                "The qptiff records technical fluorophore names but no Biomarker metadata; "
                "registration does not infer antibody identities."
            ),
        },
        "positivity": {
            "status": "not_generated",
            "reason": "registration-only export has no verified biological channel roles",
        },
        "cell_groups": {
            "status": "not_generated",
            "path": None,
        },
        "islet_annotations": dict(annotation_record),
        "uses_qupath_annotations": annotation_record.get("status") == "generated",
        "qc": {
            "status": "PASS",
            "file": "qc.json",
            "output_validation": dict(validation["checks"]),
        },
    }


def _readme(manifest: Mapping[str, Any]) -> str:
    channels = "\n".join(
        f"| {channel['index']} | `{channel['name']}` | unassigned |"
        for channel in manifest["input"]["raw_channels"]
    )
    image = manifest["image_on_xenium_grid"]["file"]
    morphology = manifest["morphology"]
    transform = manifest["alignment"]
    annotations = manifest["islet_annotations"]
    if annotations["status"] == "generated":
        annotation_section = f"""
## Import islet annotations

In Xenium Explorer 4.1 or newer, open the **Annotations** layer, choose **Import**, and select
`{annotations["file"]}`. Associate it with the **Xenium In Situ morphology image** and select
**pixels** if Explorer asks for coordinate units. The polygons have already been transformed to
that image's level-0 pixel grid. Do not apply an image-alignment CSV.

The source contains {annotations["counts"]["islet_source_features"]} Islet features. This export
keeps {annotations["counts"]["exported_islet_polygons"]} polygons at least
{annotations["minimum_equivalent_diameter_um"]:g} µm in equivalent diameter; excluded source
objects and geometry types are itemized in `manifest.json` and `qc.json`.
"""
    else:
        annotation_section = """
## Islet annotations

No same-section post-Xenium IF GeoJSON was available for this serial when the bundle was built,
so this bundle has no annotation sidecar.
"""
    return f"""# Donor {manifest["donor_id"]} — registration-only post-Xenium IF

Generated {manifest["generated_utc"]} for Xenium serial `{manifest["xenium_serial"]}`,
ROI `{manifest["roi"]}`.

## Import in Xenium Explorer

Open the Xenium bundle `{Path(manifest["xenium_bundle"]).name}`. Add `{image}` and choose
**Yes, skip alignment**. The image is already resampled onto the exact Xenium morphology grid:
{morphology["width"]} × {morphology["height"]} pixels at {OUTPUT_UM_PER_PX} µm/px with
{morphology["n_levels"]} pyramid levels. Do not apply an alignment CSV to this image.
{annotation_section}

## Channels

These are the qptiff's raw technical names, preserved in native page order.

| page | OME channel name | biological role |
|---:|---|---|
{channels}

No biomarker identity was inferred: the source qptiff contains no Biomarker metadata. This
bundle therefore contains no positivity calls and no `cell_groups` directory.

## Provenance and checks

Input calibration is {manifest["input"]["pixel_size_um"]} µm/px from qptiff XML
`PixelSizeMicrons`. The `{transform["transform_source"]}` transform maps Xenium microns to IF
microns and was consumed
only after `{Path(transform["gate_file"]).name}` reported `PASS`. The full transform
(`has_displacement={str(transform["has_displacement"]).lower()}`) was used for resampling.

`qc.json` records reopening checks for OME metadata, dimensions, uint8 dtype, channel order,
pyramid levels, physical scale, tiling, compression, and one decoded tile from every channel at
every level. `manifest.json` records the exact input, gate, transform, and output contracts.
"""


def _export_annotations(
    record: RI.RegistrationInput,
    alignment: PassedAlignment,
    morphology: Mapping[str, Any],
    root: Path,
) -> dict[str, Any]:
    """Export the optional annotation before the expensive image render starts."""
    if record.geojson is None:
        return {
            "status": "not_available",
            "file": None,
            "source_file": None,
            "reason": "no same-serial GeoJSON found in the configured post-Xenium IF directory",
        }

    relative = Path("annotations") / f"{record.serial}_islets_on_xenium_grid.geojson"
    try:
        annotation = IA.export_islet_geojson(
            record.geojson,
            Path(root) / relative,
            alignment.transform,
            input_um_per_px=float(record.pixel_um),
            output_um_per_px=OUTPUT_UM_PER_PX,
            width=int(morphology["width"]),
            height=int(morphology["height"]),
        )
    except IA.IsletAnnotationError as exc:
        raise RegisteredExportError(
            f"donor {record.donor_id} ({record.serial}) annotation export failed: {exc}"
        ) from exc
    annotation["file"] = str(relative)
    annotation["alignment_gate_sha256"] = alignment.gate_sha256
    annotation["transform_sha256"] = alignment.transform_sha256
    annotation["displacement_sha256"] = alignment.displacement_sha256
    return annotation


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _destination(record: RI.RegistrationInput, bundle_root: Path) -> Path:
    return Path(bundle_root) / f"donor_{record.donor_id}_{record.serial}"


def _partial_and_backup(destination: Path) -> tuple[Path, Path]:
    return (
        destination.parent / f".{destination.name}.partial",
        destination.parent / f".{destination.name}.previous",
    )


def _validate_existing_bundle(
    plan: ExportPlan,
    *,
    compression: str,
    tile: int,
) -> dict[str, Any]:
    destination = plan.destination
    manifest_path = destination / "manifest.json"
    qc_path = destination / "qc.json"
    manifest = _read_json_object(manifest_path, "existing registration manifest")
    for key, expected in (
        ("mode", MODE),
        ("donor_id", plan.record.donor_id),
        ("xenium_serial", plan.record.serial),
        ("roi", plan.record.roi),
    ):
        if manifest.get(key) != expected:
            raise RegisteredExportError(
                f"{manifest_path}: {key}={manifest.get(key)!r}, expected {expected!r}"
            )
    if (destination / "cell_groups").exists():
        raise RegisteredExportError(f"{destination}: registration-only bundle has cell_groups")
    if manifest.get("alignment", {}).get("gate_sha256") != plan.alignment.gate_sha256:
        raise RegisteredExportError(
            f"{destination}: existing export was not made from the current PASS gate"
        )
    if manifest.get("alignment", {}).get("transform_sha256") != plan.alignment.transform_sha256:
        raise RegisteredExportError(
            f"{destination}: existing export was not made from the current selected transform"
        )
    if (
        manifest.get("alignment", {}).get("displacement_sha256")
        != plan.alignment.displacement_sha256
    ):
        raise RegisteredExportError(
            f"{destination}: existing export was not made from the current displacement field"
        )
    input_record = manifest.get("input")
    if not isinstance(input_record, Mapping):
        raise RegisteredExportError(f"{manifest_path}: input object is absent")
    declared_qptiff = input_record.get("qptiff")
    if not isinstance(declared_qptiff, str) or not _same_path(
        Path(declared_qptiff), plan.record.qptiff
    ):
        raise RegisteredExportError(
            f"{destination}: existing export names a different source qptiff"
        )
    qptiff_stat = plan.record.qptiff.stat()
    for key, expected in (
        ("qptiff_bytes", int(qptiff_stat.st_size)),
        ("pyramid_shapes", [list(shape) for shape in plan.record.pyramid_shapes]),
        ("dapi_page", int(plan.record.dapi_page)),
    ):
        if input_record.get(key) != expected:
            raise RegisteredExportError(
                f"{destination}: existing input.{key}={input_record.get(key)!r}, "
                f"current qptiff metadata is {expected!r}"
            )
    if (
        "qptiff_mtime_ns" in input_record
        and input_record["qptiff_mtime_ns"] != int(qptiff_stat.st_mtime_ns)
    ):
        raise RegisteredExportError(
            f"{destination}: existing input.qptiff_mtime_ns is stale"
        )
    raw = input_record.get("raw_channels", [])
    names = [entry.get("name") for entry in raw if isinstance(entry, dict)]
    if names != list(plan.record.channel_names):
        raise RegisteredExportError(
            f"{destination}: existing raw channel order {names}, expected {plan.record.channel_names}"
        )
    if float(input_record.get("pixel_size_um", float("nan"))) != plan.record.pixel_um:
        raise RegisteredExportError(f"{destination}: existing qptiff calibration is stale")
    annotation = manifest.get("islet_annotations")
    if plan.record.geojson is not None:
        if not isinstance(annotation, dict) or annotation.get("status") != "generated":
            raise RegisteredExportError(
                f"{destination}: same-section islet GeoJSON is now available but the existing "
                "bundle has no generated annotation; run run_islet_annotation_export.py"
            )
        annotation_rel = annotation.get("file")
        if (
            not isinstance(annotation_rel, str)
            or Path(annotation_rel).is_absolute()
            or ".." in Path(annotation_rel).parts
        ):
            raise RegisteredExportError(
                f"{manifest_path}: invalid islet annotation path {annotation_rel!r}"
            )
        if annotation.get("source_sha256") != _sha256(plan.record.geojson):
            raise RegisteredExportError(
                f"{destination}: existing annotation was not made from the current GeoJSON"
            )
        for key, expected in (
            ("source_um_per_px", float(plan.record.pixel_um)),
            ("output_um_per_px", OUTPUT_UM_PER_PX),
            ("morphology_width_px", int(plan.morphology["width"])),
            ("morphology_height_px", int(plan.morphology["height"])),
        ):
            if annotation.get(key) != expected:
                raise RegisteredExportError(
                    f"{destination}: existing annotation {key}={annotation.get(key)!r}, "
                    f"expected {expected!r}"
                )
        if annotation.get("alignment_gate_sha256") != plan.alignment.gate_sha256:
            raise RegisteredExportError(
                f"{destination}: existing annotation was not made from the current PASS gate"
            )
        if annotation.get("transform_sha256") != plan.alignment.transform_sha256:
            raise RegisteredExportError(
                f"{destination}: existing annotation was not made from the current transform"
            )
        if annotation.get("displacement_sha256") != plan.alignment.displacement_sha256:
            raise RegisteredExportError(
                f"{destination}: existing annotation was not made from the current displacement"
            )
        try:
            expected_annotation_features = annotation.get("counts", {}).get(
                "exported_islet_polygons"
            )
            if type(expected_annotation_features) is not int or expected_annotation_features < 1:
                raise RegisteredExportError(
                    f"{destination}: existing annotation has no valid exported polygon count"
                )
            annotation_validation = IA.validate_islet_geojson(
                destination / annotation_rel,
                width=int(plan.morphology["width"]),
                height=int(plan.morphology["height"]),
                expected_features=expected_annotation_features,
            )
            canonical_document, canonical_provenance = IA.build_islet_geojson(
                plan.record.geojson,
                plan.alignment.transform,
                input_um_per_px=float(plan.record.pixel_um),
                output_um_per_px=OUTPUT_UM_PER_PX,
                width=int(plan.morphology["width"]),
                height=int(plan.morphology["height"]),
            )
        except RegisteredExportError:
            raise
        except IA.IsletAnnotationError as exc:
            raise RegisteredExportError(
                f"{destination}: existing islet annotation failed validation: {exc}"
            ) from exc
        if annotation_validation["sha256"] != annotation.get("output_sha256"):
            raise RegisteredExportError(
                f"{destination}: existing annotation bytes differ from its manifest hash"
            )
        canonical_bytes = (
            json.dumps(canonical_document, separators=(",", ":"), allow_nan=False) + "\n"
        ).encode()
        canonical_sha256 = hashlib.sha256(canonical_bytes).hexdigest()
        if annotation_validation["sha256"] != canonical_sha256:
            raise RegisteredExportError(
                f"{destination}: existing annotation is not the canonical current "
                f">={IA.MIN_ISLET_DIAMETER_UM:g} um full-transform export"
            )
        for key in (
            "source_coordinate_space",
            "source_um_per_px",
            "class",
            "minimum_equivalent_diameter_um",
            "output_coordinate_space",
            "output_um_per_px",
            "morphology_width_px",
            "morphology_height_px",
            "transform_direction_stored",
            "transform_operation_used",
            "includes_displacement",
            "counts",
        ):
            if annotation.get(key) != canonical_provenance[key]:
                raise RegisteredExportError(
                    f"{destination}: existing annotation provenance field {key!r} is stale"
                )
        if manifest.get("uses_qupath_annotations") is not True:
            raise RegisteredExportError(
                f"{manifest_path}: generated islet annotations are not declared at top level"
            )
    else:
        if isinstance(annotation, dict) and annotation.get("status") == "generated":
            raise RegisteredExportError(
                f"{destination}: annotation source is no longer available for validation"
            )
        if (destination / "annotations").exists():
            raise RegisteredExportError(
                f"{destination}: annotation files exist without a current same-section source"
            )
        annotation_validation = None
    image_rel = manifest.get("image_on_xenium_grid", {}).get("file")
    if (
        not isinstance(image_rel, str)
        or Path(image_rel).is_absolute()
        or ".." in Path(image_rel).parts
    ):
        raise RegisteredExportError(f"{manifest_path}: invalid registered image path {image_rel!r}")
    qc = _read_json_object(qc_path, "existing registration QC")
    if qc.get("status") != "PASS":
        raise RegisteredExportError(f"{qc_path}: status is {qc.get('status')!r}, expected PASS")
    validation = validate_registered_ome(
        destination / image_rel,
        channel_names=plan.record.channel_names,
        morphology=plan.morphology,
        compression=compression,
        tile=tile,
    )
    return {
        "status": "SKIPPED_EXISTING",
        "destination": str(destination),
        "validation": validation,
        "annotation_validation": annotation_validation,
        "manifest": manifest,
    }


def prepare_export(
    record: RI.RegistrationInput,
    *,
    alignment_root: Path = ALIGNMENT_ROOT,
    bundle_root: Path = BUNDLE_ROOT,
    transform_stem: str = "refined",
    overwrite: bool = False,
    skip_existing: bool = False,
    compression: str = "zlib",
    tile: int = DEFAULT_TILE,
) -> ExportPlan:
    """Perform all cheap gates before allocating or decoding a full registered image."""
    if record.donor_id in PROTECTED_COMPLETED_DONORS:
        raise RegisteredExportError(
            f"donor {record.donor_id} is a protected completed Explorer export and is never "
            "rewritten by this registration-only lane"
        )
    if overwrite and skip_existing:
        raise RegisteredExportError("overwrite and skip_existing are mutually exclusive")
    RI.require_registration_files(record)
    alignment = load_passed_alignment(record, alignment_root, transform_stem=transform_stem)
    morphology = EX.morphology_grid(record.bundle)
    if abs(float(morphology["um_per_px"]) - OUTPUT_UM_PER_PX) > 1e-12:
        raise RegisteredExportError(
            f"donor {record.donor_id}: morphology is {morphology['um_per_px']} um/px, "
            f"expected {OUTPUT_UM_PER_PX}"
        )

    destination = _destination(record, bundle_root)
    partial, backup = _partial_and_backup(destination)
    if partial.exists() or backup.exists():
        stale = partial if partial.exists() else backup
        raise RegisteredExportError(
            f"{stale} remains from an interrupted export; inspect/remove it before retrying"
        )
    if destination.exists():
        if skip_existing:
            action = "skip"
        elif overwrite:
            action = "overwrite"
        else:
            raise RegisteredExportError(
                f"{destination} already exists; use --skip-existing to validate/keep it or "
                "--overwrite to replace it"
            )
    else:
        action = "write"
    return ExportPlan(record, alignment, morphology, destination, action)


def _publish(partial: Path, destination: Path, *, overwrite: bool) -> None:
    """Publish one validated directory by rename, restoring the previous one on failure."""
    _, backup = _partial_and_backup(destination)
    if not destination.exists():
        partial.rename(destination)
        return
    if not overwrite:
        raise RegisteredExportError(f"refusing to replace existing {destination}")
    if backup.exists():
        raise RegisteredExportError(f"refusing overwrite because backup path exists: {backup}")
    destination.rename(backup)
    try:
        partial.rename(destination)
    except BaseException:
        backup.rename(destination)
        raise
    shutil.rmtree(backup)


def export_one(
    plan: ExportPlan,
    *,
    compression: str = "zlib",
    workers: int = 8,
    tile: int = DEFAULT_TILE,
) -> dict[str, Any]:
    """Build, validate, and atomically publish one registration-only donor bundle."""
    if plan.record.donor_id in PROTECTED_COMPLETED_DONORS:
        raise RegisteredExportError(
            f"donor {plan.record.donor_id} is protected and cannot be exported"
        )
    if plan.action == "skip":
        print(f"[registered] donor {plan.record.donor_id}: validating existing bundle", flush=True)
        return _validate_existing_bundle(plan, compression=compression, tile=tile)

    destination = plan.destination
    partial, _ = _partial_and_backup(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial.mkdir()
    image_name = f"{plan.record.serial}_IF_on_xenium_grid.ome.tif"
    image_path = partial / "images" / image_name
    image_path.parent.mkdir()

    try:
        # Geometry is small and strict to validate.  Do it before allocating the registered
        # raster so a malformed annotation cannot waste a multi-gigabyte image render.
        annotation_record = _export_annotations(
            plan.record,
            plan.alignment,
            plan.morphology,
            partial,
        )
        level0, warp_record = warp_qptiff_to_xenium_level0(
            plan.record,
            plan.alignment.transform,
            plan.morphology,
            tile=tile,
            workers=workers,
        )
        levels = EX.xenium_pyramid(level0, int(plan.morphology["n_levels"]))
        del level0
        expected_shapes = [tuple(shape) for shape in plan.morphology["level_shapes"]]
        if [tuple(level.shape[-2:]) for level in levels] != expected_shapes:
            raise RegisteredExportError(
                "generated pyramid does not exactly match morphology levels"
            )

        image_record = EX.write_warped_ome_tiff(
            levels,
            image_path,
            list(plan.record.channel_names),
            um_per_px=OUTPUT_UM_PER_PX,
            compression=compression,
            tile=tile,
            workers=int(workers),
        )
        del levels
        image_record = dict(image_record)
        image_record["path"] = f"images/{image_name}"
        validation = validate_registered_ome(
            image_path,
            channel_names=plan.record.channel_names,
            morphology=plan.morphology,
            compression=compression,
            tile=tile,
        )

        manifest = _manifest(
            plan.record,
            plan.alignment,
            plan.morphology,
            image_record,
            warp_record,
            validation,
            image_name,
            annotation_record,
        )
        qc = {
            "status": "PASS",
            "mode": MODE,
            "donor_id": plan.record.donor_id,
            "xenium_serial": plan.record.serial,
            "roi": plan.record.roi,
            "alignment_gate": {
                "file": str(plan.alignment.gate_path),
                "sha256": plan.alignment.gate_sha256,
                "verdict": "PASS",
            },
            "output": validation,
            "islet_annotations": annotation_record,
            "forbidden_products": {
                "cell_groups_absent": True,
                "positivity_not_generated": True,
                "biological_channel_roles_unassigned": True,
            },
        }
        _write_json(partial / "manifest.json", manifest)
        _write_json(partial / "qc.json", qc)
        (partial / "README.md").write_text(_readme(manifest))
        if (partial / "cell_groups").exists():
            raise RegisteredExportError(
                "registration-only partial unexpectedly contains cell_groups"
            )
        _publish(partial, destination, overwrite=plan.action == "overwrite")
    except BaseException:
        if partial.exists():
            shutil.rmtree(partial)
        raise

    print(f"[registered] donor {plan.record.donor_id}: PASS -> {destination}", flush=True)
    return {"status": "EXPORTED", "destination": str(destination), "manifest": manifest}


def _select_remaining(
    records: Mapping[str, RI.RegistrationInput],
    donors: Sequence[str] | None,
) -> tuple[list[RI.RegistrationInput], list[RI.RegistrationInput]]:
    selected = RI.select_registration_inputs(
        records, donors, context="registration-only prewarped OME export"
    )
    protected = [record for record in selected if record.donor_id in PROTECTED_COMPLETED_DONORS]
    remaining = [record for record in selected if record.donor_id not in PROTECTED_COMPLETED_DONORS]
    return remaining, protected


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
        "--transform-stem",
        default="refined",
        help=(
            "exact transform/gate stem: reads <stem>.json and alignment_<stem>.json; "
            "never falls back to another stem"
        ),
    )
    parser.add_argument(
        "--compression", default="zlib", choices=["zlib", "lzw", "jpeg2000", "none"]
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--tile", type=int, default=DEFAULT_TILE)
    existing = parser.add_mutually_exclusive_group()
    existing.add_argument(
        "--skip-existing",
        action="store_true",
        help="revalidate and keep an existing registration-only bundle",
    )
    existing.add_argument(
        "--overwrite",
        action="store_true",
        help="replace a non-protected bundle only after its replacement passes QC",
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
        selected, protected = _select_remaining(records, args.donors)
        for record in protected:
            print(
                f"[registered] donor {record.donor_id}: protected completed bundle untouched",
                flush=True,
            )
        if not selected:
            if protected:
                return 0
            raise RegisteredExportError("no eligible remaining donor matched")

        # Preflight the entire requested batch so a missing PASS gate or an existing destination
        # is reported before any earlier donor spends hours writing a multi-gigabyte image.
        plans = [
            prepare_export(
                record,
                alignment_root=args.alignment_root,
                bundle_root=args.bundle_root,
                transform_stem=args.transform_stem,
                overwrite=args.overwrite,
                skip_existing=args.skip_existing,
                compression=args.compression,
                tile=args.tile,
            )
            for record in selected
        ]
        results = [
            export_one(plan, compression=args.compression, workers=args.workers, tile=args.tile)
            for plan in plans
        ]
    except (RI.RegistrationInputError, RegisteredExportError) as exc:
        print(f"[registered] ERROR: {exc}", file=sys.stderr)
        return 2

    exported = sum(result["status"] == "EXPORTED" for result in results)
    skipped = sum(result["status"] == "SKIPPED_EXISTING" for result in results)
    print(
        f"[registered] complete: {exported} exported, {skipped} validated/skipped, "
        f"{len(protected)} protected",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
