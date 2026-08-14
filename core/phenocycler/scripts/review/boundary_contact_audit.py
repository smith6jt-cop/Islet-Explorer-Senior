#!/usr/bin/env python3
"""Audit segmentation-boundary overlap in a mutable development workspace.

This review-only tool reads a freshly projected subset of the upstream QuPath
cell export, current geometry-QC, raw-image, and explicitly supplied scratch
files.  It never reads RESTORE output or live REDSEA morphology.  It does not
validate production stage manifests, adoption receipts, file fingerprints, or
historical snapshots, and it never mutates ``data/runs``.

The scratch contact graph is authoritative for development diagnostics.  Only
object IDs and a structurally valid CSR graph are required; old scratch files
may omit every compartment array.  Coverage against the current QuPath object
table is reported, but a partial overlap does not suppress valid contact rows.
Likewise, the upstream QuPath projection is inner-joined to current geometry-QC
by exact object ID, with unmatched IDs on either side reported rather than
treated as a development-audit failure.
Raw target/reference/DAPI crops remain independent of scratch availability.
The default workspace is ``../data/development/boundary_contact_diagnostics``;
rerunning the command replaces only that tool-owned development directory.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import sys
import tempfile
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
import tifffile

REPOSITORY = Path(__file__).resolve().parents[2]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from phenocycler.artifacts import (  # noqa: E402
    ContractError,
    iter_qupath_geojson_features,
)
from phenocycler.redsea import validated_channel_markers  # noqa: E402
from scripts.review.fundamental_pair_survey import (  # noqa: E402
    ROSTER_PATH as DEFAULT_ROSTER,
)
from scripts.review.fundamental_pair_survey import (  # noqa: E402
    load_diagnostic_roster,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    DEFAULT_OUTPUT_ROOT as DEFAULT_QUPATH_PROJECTION_ROOT,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    DEFAULT_SPEC as DEFAULT_COHORT_SPEC,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    DIAGNOSTIC_SOURCE_KINDS,
    SOURCE_MEASUREMENTS,
    load_qupath_projection_donor,
    refresh_qupath_projection,
)

SCRATCH_AVAILABLE = "AVAILABLE"
SCRATCH_MISSING = "MISSING"
SCRATCH_INVALID = "INVALID_SCRATCH"
DEFAULT_DEVELOPMENT_ROOT = REPOSITORY.parent / "data" / "development"
DEFAULT_LABEL = "boundary_contact_diagnostics"

REQUIRED_SCRATCH_KEYS = frozenset({"oids", "m_data", "m_indices", "m_indptr", "m_shape"})
STRATA = ("dual_bright", "target_bright", "reference_bright")
DISPLAY_MARKER_OVERRIDES = {"Iba1_norm": "Iba1"}
WORKSPACE_MARKER = ".boundary_contact_development"


class ScratchDataError(ContractError):
    """The development scratch graph cannot be interpreted safely."""


class TiledPageReader:
    """Read arbitrary level-0 qptiff pixels without materializing a channel."""

    def __init__(self, filehandle, page, *, cache_tiles: int = 512) -> None:
        if not bool(page.is_tiled):
            raise ContractError("boundary qptiff measurement requires tiled pages")
        if cache_tiles <= 0:
            raise ContractError("cache_tiles must be positive")
        self.filehandle = filehandle
        self.page = page
        self.tile_height = int(page.tilelength)
        self.tile_width = int(page.tilewidth)
        self.image_height = int(page.imagelength)
        self.image_width = int(page.imagewidth)
        self.tiles_across = int(math.ceil(self.image_width / self.tile_width))
        self.cache_tiles = int(cache_tiles)
        self._cache: OrderedDict[int, np.ndarray] = OrderedDict()

    def _tile(self, index: int) -> np.ndarray:
        cached = self._cache.pop(index, None)
        if cached is not None:
            self._cache[index] = cached
            return cached
        offsets = self.page.dataoffsets
        counts = self.page.databytecounts
        if index < 0 or index >= len(offsets) or counts[index] <= 0:
            raise ContractError(f"qptiff tile {index} is absent")
        self.filehandle.seek(offsets[index])
        raw = self.filehandle.read(counts[index])
        decoded = self.page.decode(raw, index)
        tile = np.asarray(decoded[0]).squeeze()
        if tile.ndim != 2:
            raise ContractError(f"qptiff tile {index} decoded to shape {tile.shape}")
        self._cache[index] = tile
        while len(self._cache) > self.cache_tiles:
            self._cache.popitem(last=False)
        return tile

    def values(self, rows: Sequence[int], columns: Sequence[int]) -> np.ndarray:
        rows_array = np.asarray(rows, dtype=int)
        columns_array = np.asarray(columns, dtype=int)
        if rows_array.shape != columns_array.shape:
            raise ContractError("qptiff row and column vectors must align")
        if (
            np.any(rows_array < 0)
            or np.any(rows_array >= self.image_height)
            or np.any(columns_array < 0)
            or np.any(columns_array >= self.image_width)
        ):
            raise ContractError("qptiff pixel request falls outside the image")
        tile_row = rows_array // self.tile_height
        tile_column = columns_array // self.tile_width
        tile_index = tile_row * self.tiles_across + tile_column
        output = np.empty(len(rows_array), dtype=float)
        for index in np.unique(tile_index):
            selected = tile_index == index
            tile = self._tile(int(index))
            local_rows = rows_array[selected] % self.tile_height
            local_columns = columns_array[selected] % self.tile_width
            output[selected] = tile[local_rows, local_columns]
        return output


@dataclass(frozen=True, order=True)
class MarkerPair:
    target: str
    reference: str

    @property
    def pair_id(self) -> str:
        return f"{self.target}:{self.reference}"


@dataclass(frozen=True)
class LivePath:
    path: str


@dataclass(frozen=True)
class LiveChannel:
    channel_index: int
    channel_name: str
    marker: str


@dataclass(frozen=True)
class LiveImage:
    donor_id: str
    image_id: str
    qptiff: LivePath
    cell_geojson: LivePath
    nucleus_geojson: LivePath
    combined_geometry_file: bool
    pixel_size_um_x: float
    pixel_size_um_y: float
    channel_map: tuple[LiveChannel, ...]


def pairs_from_roster(path: Path = DEFAULT_ROSTER) -> tuple[MarkerPair, ...]:
    roster = load_diagnostic_roster(path)
    return tuple(MarkerPair(pair.target_column, pair.reference_column) for pair in roster.pairs)


DEFAULT_PAIRS = pairs_from_roster()


@dataclass
class ScratchBundle:
    oids: np.ndarray
    contact: sp.csr_matrix


@dataclass
class ScratchValidation:
    bundle: ScratchBundle
    current_object_ids: np.ndarray
    current_scratch_indices: np.ndarray
    current_only_n: int
    scratch_only_n: int
    current_only_examples: tuple[str, ...]
    scratch_only_examples: tuple[str, ...]


@dataclass
class ScratchAudit:
    state: str
    reason: str
    validation: ScratchValidation | None = None


@dataclass(frozen=True)
class ProjectionGeometryCoverage:
    projection_cell_n: int
    current_cell_n: int
    matched_cell_n: int
    projection_only_n: int
    current_only_n: int
    projection_only_examples: tuple[str, ...]
    current_only_examples: tuple[str, ...]


def display_marker(marker: str) -> str:
    """Map a diagnostic measurement name to an explicitly reviewed raw channel."""

    return DISPLAY_MARKER_OVERRIDES.get(str(marker), str(marker))


def selection_source_measurement(marker: str) -> str:
    """Return the exact upstream QuPath measurement behind a selection alias."""

    try:
        return SOURCE_MEASUREMENTS[str(marker)]
    except KeyError as exc:
        raise ContractError(
            f"unsupported boundary selection marker {marker!r}; the mutable QuPath "
            "projection deliberately contains only the fundamental inputs"
        ) from exc


def selection_source_kind(marker: str) -> str:
    """Label a selection value without conflating it with its raw display channel."""

    if str(marker) in {"DAPI", "X_centroid", "Y_centroid", "cell_area", "cell_solidity"}:
        return "qupath_raw_cell_measurement"
    try:
        return DIAGNOSTIC_SOURCE_KINDS[str(marker)]
    except KeyError as exc:
        raise ContractError(f"no declared QuPath source kind for {marker!r}") from exc


def load_live_images(path: Path) -> tuple[LiveImage, ...]:
    """Read current source paths and channel maps without content-ID validation."""

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        records = payload["images"]
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ContractError(f"cannot read current QuPath image registrations: {path}") from exc
    if not isinstance(records, list) or not records:
        raise ContractError("current QuPath image registrations must be a nonempty list")
    images: list[LiveImage] = []
    for index, record in enumerate(records):
        try:
            channels = tuple(
                LiveChannel(
                    channel_index=int(channel["channel_index"]),
                    channel_name=str(channel["channel_name"]),
                    marker=str(channel["marker"]),
                )
                for channel in record["channel_map"]
            )
            image = LiveImage(
                donor_id=str(record["donor_id"]),
                image_id=str(record["image_id"]),
                qptiff=LivePath(str(record["qptiff"]["path"])),
                cell_geojson=LivePath(str(record["cell_geojson"]["path"])),
                nucleus_geojson=LivePath(str(record["nucleus_geojson"]["path"])),
                combined_geometry_file=bool(record["combined_geometry_file"]),
                pixel_size_um_x=float(record["pixel_size_um_x"]),
                pixel_size_um_y=float(record["pixel_size_um_y"]),
                channel_map=channels,
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ContractError(f"invalid current QuPath image registration {index}") from exc
        if not image.donor_id or not image.image_id:
            raise ContractError(f"current QuPath image registration {index} has an empty ID")
        if image.pixel_size_um_x <= 0 or image.pixel_size_um_y <= 0:
            raise ContractError(f"{image.donor_id}: pixel calibration must be positive")
        channel_indices = tuple(channel.channel_index for channel in channels)
        if channel_indices != tuple(range(len(channels))):
            raise ContractError(
                f"{image.donor_id}: channel indices must exactly cover 0..N-1 in order"
            )
        if len({channel.marker for channel in channels}) != len(channels):
            raise ContractError(f"{image.donor_id}: channel markers must be unique")
        images.append(image)
    if len({image.donor_id for image in images}) != len(images):
        raise ContractError("current QuPath image registrations contain duplicate donors")
    return tuple(images)


def parse_pairs(
    values: Sequence[str] | None,
    *,
    roster_path: Path = DEFAULT_ROSTER,
) -> tuple[MarkerPair, ...]:
    """Parse explicit pairs or derive them from the editable shared roster."""

    if not values:
        return pairs_from_roster(roster_path)
    pairs: list[MarkerPair] = []
    for value in values:
        pieces = str(value).split(":")
        if len(pieces) != 2 or not all(piece.strip() for piece in pieces):
            raise ContractError(f"invalid pair {value!r}; expected TARGET:REFERENCE")
        pair = MarkerPair(*(piece.strip() for piece in pieces))
        if pair.target == pair.reference:
            raise ContractError(f"pair {value!r} uses the same marker twice")
        pairs.append(pair)
    if len(set(pairs)) != len(pairs):
        raise ContractError("diagnostic pairs must be unique")
    return tuple(pairs)


def _as_string_vector(values: Any, *, name: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1:
        raise ScratchDataError(f"{name} must be one-dimensional")
    if pd.isna(array).any():
        raise ScratchDataError(f"{name} contains a null value")
    strings = np.asarray([str(value) for value in array], dtype=object)
    if any(not value for value in strings):
        raise ScratchDataError(f"{name} contains an empty value")
    if len(set(strings)) != len(strings):
        raise ScratchDataError(f"{name} contains duplicates")
    return strings


def _scratch_bundle(values: Mapping[str, Any]) -> ScratchBundle:
    missing = sorted(REQUIRED_SCRATCH_KEYS.difference(values))
    if missing:
        raise ScratchDataError(f"missing contact arrays: {missing}")

    oids = _as_string_vector(values["oids"], name="scratch object IDs")
    n_cells = len(oids)
    shape = np.asarray(values["m_shape"])
    m_data = np.asarray(values["m_data"])
    m_indices = np.asarray(values["m_indices"])
    m_indptr = np.asarray(values["m_indptr"])
    if (
        shape.shape != (2,)
        or not np.issubdtype(shape.dtype, np.number)
        or not np.equal(shape, np.floor(shape)).all()
        or tuple(map(int, shape)) != (n_cells, n_cells)
    ):
        raise ScratchDataError(
            f"scratch contact shape is {tuple(shape)}, expected {(n_cells, n_cells)}"
        )
    if m_data.ndim != 1 or m_indices.ndim != 1 or m_indptr.shape != (n_cells + 1,):
        raise ScratchDataError("scratch contact CSR arrays have incompatible dimensions")
    if len(m_data) != len(m_indices):
        raise ScratchDataError("scratch contact CSR data/index lengths differ")
    if not np.issubdtype(m_indices.dtype, np.integer) or not np.issubdtype(
        m_indptr.dtype, np.integer
    ):
        raise ScratchDataError("scratch contact CSR indices must be integer")
    if not np.issubdtype(m_data.dtype, np.number):
        raise ScratchDataError("scratch contact CSR values must be numeric")
    if not np.isfinite(m_data).all() or np.any(m_data < 0):
        raise ScratchDataError("scratch contact CSR values must be finite and nonnegative")
    try:
        contact = sp.csr_matrix((m_data, m_indices, m_indptr), shape=(n_cells, n_cells), copy=False)
        contact.check_format(full_check=True)
    except (TypeError, ValueError) as exc:
        raise ScratchDataError(f"scratch contact CSR cannot be reconstructed: {exc}") from exc
    if not contact.has_canonical_format:
        raise ScratchDataError("scratch contact CSR is not in canonical form")
    if np.any(contact.diagonal() != 0):
        raise ScratchDataError("scratch contact CSR has a nonzero diagonal")
    difference = contact - contact.T
    if difference.nnz and np.any(difference.data != 0):
        raise ScratchDataError("scratch contact CSR is not symmetric")

    return ScratchBundle(oids=oids, contact=contact)


def _current_scratch_mapping(
    bundle: ScratchBundle,
    current_object_ids: Sequence[object],
) -> ScratchValidation:
    current = _as_string_vector(current_object_ids, name="current object IDs")
    lookup = {str(object_id): index for index, object_id in enumerate(bundle.oids)}
    matched = [str(object_id) for object_id in current if str(object_id) in lookup]
    current_set = set(map(str, current))
    scratch_set = set(map(str, bundle.oids))
    current_only = sorted(current_set - scratch_set)
    scratch_only = sorted(scratch_set - current_set)
    return ScratchValidation(
        bundle=bundle,
        current_object_ids=np.asarray(matched, dtype=object),
        current_scratch_indices=np.asarray([lookup[object_id] for object_id in matched]),
        current_only_n=len(current_only),
        scratch_only_n=len(scratch_only),
        current_only_examples=tuple(current_only[:5]),
        scratch_only_examples=tuple(scratch_only[:5]),
    )


def validate_scratch_mapping(
    values: Mapping[str, Any],
    current_object_ids: Sequence[object],
) -> ScratchValidation:
    """Validate a contact graph and measure overlap with the current cell table."""

    bundle = _scratch_bundle(values)
    return _current_scratch_mapping(bundle, current_object_ids)


def audit_scratch_mapping(
    values: Mapping[str, Any],
    current_object_ids: Sequence[object],
) -> ScratchAudit:
    """Return an availability status without comparing historical snapshots."""

    try:
        validation = validate_scratch_mapping(values, current_object_ids)
    except (ContractError, ValueError, IndexError, KeyError) as exc:
        return ScratchAudit(SCRATCH_INVALID, str(exc))
    return ScratchAudit(SCRATCH_AVAILABLE, "", validation=validation)


def audit_scratch_file(
    path: Path,
    current_object_ids: Sequence[object],
) -> ScratchAudit:
    """Load the current NPZ once and validate only its contact graph."""

    if not path.is_file():
        return ScratchAudit(SCRATCH_MISSING, f"scratch archive is missing: {path}")
    try:
        with np.load(path, allow_pickle=True) as archive:
            values = {name: archive[name] for name in REQUIRED_SCRATCH_KEYS if name in archive}
        return audit_scratch_mapping(values, current_object_ids)
    except (OSError, ValueError, ContractError) as exc:
        return ScratchAudit(SCRATCH_INVALID, f"cannot read scratch archive: {exc}")


def contact_metrics(validation: ScratchValidation) -> pd.DataFrame:
    """Return scratch contact geometry for IDs present in the current cell table."""

    bundle = validation.bundle
    selected = validation.current_scratch_indices
    degree = np.diff(bundle.contact.indptr)
    total = np.asarray(bundle.contact.sum(axis=1)).ravel()
    maximum = np.zeros(len(bundle.oids), dtype=float)
    for row in range(len(bundle.oids)):
        values = bundle.contact.data[bundle.contact.indptr[row] : bundle.contact.indptr[row + 1]]
        maximum[row] = float(values.max()) if len(values) else 0.0
    return pd.DataFrame(
        {
            "object_id": validation.current_object_ids,
            "contact_neighbor_n": degree[selected].astype(np.int64),
            "contact_boundary_px": total[selected],
            "max_shared_boundary_px": maximum[selected],
            "contact_isolated": degree[selected] == 0,
        }
    )


def selected_contact_edges(
    validation: ScratchValidation,
    selected_object_ids: Sequence[object],
) -> pd.DataFrame:
    """Return one-hop scratch contacts for crop centers, including context extras."""

    bundle = validation.bundle
    lookup = {str(object_id): index for index, object_id in enumerate(bundle.oids)}
    current = set(map(str, validation.current_object_ids))
    rows: list[dict[str, object]] = []
    for center_id in sorted(set(map(str, selected_object_ids))):
        if center_id not in lookup:
            continue
        row_index = lookup[center_id]
        start, end = bundle.contact.indptr[row_index : row_index + 2]
        order = sorted(
            zip(
                bundle.contact.indices[start:end],
                bundle.contact.data[start:end],
                strict=True,
            ),
            key=lambda item: (-float(item[1]), str(bundle.oids[item[0]])),
        )
        for neighbor_index, shared in order:
            neighbor_id = str(bundle.oids[neighbor_index])
            rows.append(
                {
                    "object_id": center_id,
                    "neighbor_object_id": neighbor_id,
                    "neighbor_is_current": neighbor_id in current,
                    "shared_boundary_px": float(shared),
                }
            )
    return pd.DataFrame.from_records(
        rows,
        columns=(
            "object_id",
            "neighbor_object_id",
            "neighbor_is_current",
            "shared_boundary_px",
        ),
    )


def _percentile_rank(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").rank(method="average", pct=True)


def select_crop_manifest(
    candidates: pd.DataFrame,
    *,
    donor_id: str,
    image_id: str,
    pairs: Sequence[MarkerPair],
    image_shape: tuple[int, int],
    pixel_size_um_x: float,
    pixel_size_um_y: float,
    crops_per_stratum: int = 4,
    crop_size_px: int = 256,
) -> pd.DataFrame:
    """Choose stable raw-intensity crops without identity or outcome columns."""

    if crops_per_stratum <= 0 or crop_size_px <= 0:
        raise ContractError("crop counts and dimensions must be positive")
    required = {
        "object_id",
        "X_centroid",
        "Y_centroid",
        "DAPI",
        "analysis_eligible",
        *(pair.target for pair in pairs),
        *(pair.reference for pair in pairs),
    }
    missing = sorted(required.difference(candidates.columns))
    if missing:
        raise ContractError(f"crop candidates lack required columns {missing}")
    if pixel_size_um_x <= 0 or pixel_size_um_y <= 0:
        raise ContractError("pixel calibration must be positive")
    height, width = map(int, image_shape)
    if height <= 0 or width <= 0:
        raise ContractError("image shape must be positive")

    base_columns = [
        "object_id",
        "X_centroid",
        "Y_centroid",
        "DAPI",
        "analysis_eligible",
        "cell_area",
        "cell_solidity",
        "raster_area_ratio",
        "nucleus_cell_area_ratio",
    ]
    marker_columns = tuple(
        dict.fromkeys(marker for pair in pairs for marker in (pair.target, pair.reference))
    )
    available = [column for column in base_columns if column in candidates]
    frame = candidates[[*available, *marker_columns]].copy()
    frame["object_id"] = frame["object_id"].astype("string")
    if frame["object_id"].isna().any() or frame["object_id"].duplicated().any():
        raise ContractError("crop candidate object IDs must be non-null and unique")
    frame = frame[frame["analysis_eligible"].fillna(False).astype(bool)].copy()

    rows: list[dict[str, object]] = []
    for pair_index, pair in enumerate(pairs):
        usable = frame.copy()
        numeric_columns = [pair.target, pair.reference, "DAPI", "X_centroid", "Y_centroid"]
        for column in numeric_columns:
            try:
                usable[column] = pd.to_numeric(usable[column], errors="raise")
            except (TypeError, ValueError) as exc:
                raise ContractError(f"{pair.pair_id}: crop input {column} is not numeric") from exc
        usable = usable[np.isfinite(usable[numeric_columns]).all(axis=1)].copy()
        if usable.empty:
            raise ContractError(f"{pair.pair_id}: no finite crop candidates")
        usable["_target_rank"] = _percentile_rank(usable[pair.target])
        usable["_reference_rank"] = _percentile_rank(usable[pair.reference])
        scores = {
            "dual_bright": np.minimum(usable["_target_rank"], usable["_reference_rank"]),
            "target_bright": usable["_target_rank"] - usable["_reference_rank"],
            "reference_bright": usable["_reference_rank"] - usable["_target_rank"],
        }
        already_selected: set[str] = set()
        for stratum_index, stratum in enumerate(STRATA):
            scored = usable.assign(_score=np.asarray(scores[stratum], dtype=float))
            scored = scored.sort_values(
                ["_score", "object_id"], ascending=[False, True], kind="mergesort"
            )
            scored = scored[~scored["object_id"].astype(str).isin(already_selected)]
            chosen = scored.head(crops_per_stratum)
            if len(chosen) != crops_per_stratum:
                raise ContractError(
                    f"{pair.pair_id} {stratum}: requested {crops_per_stratum} crops but "
                    f"only {len(chosen)} non-overlapping candidates remain"
                )
            for crop_rank, (_, cell) in enumerate(chosen.iterrows(), start=1):
                object_id = str(cell["object_id"])
                already_selected.add(object_id)
                center_x = int(np.rint(float(cell["X_centroid"]) / pixel_size_um_x))
                center_y = int(np.rint(float(cell["Y_centroid"]) / pixel_size_um_y))
                crop_width = min(crop_size_px, width)
                crop_height = min(crop_size_px, height)
                x0 = int(np.clip(center_x - crop_width // 2, 0, width - crop_width))
                y0 = int(np.clip(center_y - crop_height // 2, 0, height - crop_height))
                record: dict[str, object] = {
                    "donor_id": str(donor_id),
                    "image_id": str(image_id),
                    "pair_id": pair.pair_id,
                    "target": pair.target,
                    "reference": pair.reference,
                    "target_selection_source_measurement": selection_source_measurement(
                        pair.target
                    ),
                    "reference_selection_source_measurement": selection_source_measurement(
                        pair.reference
                    ),
                    "target_selection_source_kind": selection_source_kind(pair.target),
                    "reference_selection_source_kind": selection_source_kind(pair.reference),
                    "target_display_marker": display_marker(pair.target),
                    "reference_display_marker": display_marker(pair.reference),
                    "target_display_channel": display_marker(pair.target),
                    "reference_display_channel": display_marker(pair.reference),
                    "stratum": stratum,
                    "pair_finite_n": int(len(usable)),
                    "stratum_requested_n": int(crops_per_stratum),
                    "stratum_selected_n": int(len(chosen)),
                    "crop_rank": crop_rank,
                    "object_id": object_id,
                    "center_x_px": center_x,
                    "center_y_px": center_y,
                    "x0_px": x0,
                    "y0_px": y0,
                    "width_px": crop_width,
                    "height_px": crop_height,
                    "selection_target_cell_value": float(cell[pair.target]),
                    "selection_reference_cell_value": float(cell[pair.reference]),
                    "raw_dapi_cell_mean": float(cell["DAPI"]),
                    "target_percentile_rank": float(cell["_target_rank"]),
                    "reference_percentile_rank": float(cell["_reference_rank"]),
                    "selection_score": float(cell["_score"]),
                    "_pair_index": pair_index,
                    "_stratum_index": stratum_index,
                }
                for source, output in (
                    ("cell_area", "qupath_cell_area_um2"),
                    ("cell_solidity", "qupath_cell_solidity"),
                    ("raster_area_ratio", "raster_area_ratio"),
                    ("nucleus_cell_area_ratio", "nucleus_cell_area_ratio"),
                ):
                    value = cell.get(source, np.nan)
                    record[output] = float(value) if pd.notna(value) else np.nan
                rows.append(record)
    output = pd.DataFrame.from_records(rows)
    if output.empty:
        return output
    output = output.sort_values(
        ["_pair_index", "_stratum_index", "crop_rank", "object_id"],
        kind="mergesort",
        ignore_index=True,
    )
    return output.drop(columns=["_pair_index", "_stratum_index"])


def extract_tiled_window(reader: TiledPageReader, *, x0: int, y0: int, width: int, height: int):
    """Read exactly one window through ``TiledPageReader.values``."""

    rows, columns = np.indices((height, width), dtype=np.int64)
    values = reader.values((rows + y0).ravel(), (columns + x0).ravel())
    return np.asarray(values).reshape(height, width)


def _geometry_bounds(
    geometry: Mapping[str, Any] | None,
) -> tuple[float, float, float, float] | None:
    if not isinstance(geometry, Mapping) or geometry.get("type") not in {
        "Polygon",
        "MultiPolygon",
    }:
        return None
    polygons = geometry.get("coordinates", ())
    if geometry.get("type") == "Polygon":
        polygons = [polygons]
    coordinates: list[np.ndarray] = []
    for polygon in polygons:
        for ring in polygon:
            array = np.asarray(ring, dtype=float)
            if array.ndim == 2 and array.shape[1] >= 2 and len(array):
                coordinates.append(array[:, :2])
    if not coordinates:
        return None
    joined = np.concatenate(coordinates)
    return (
        float(joined[:, 0].min()),
        float(joined[:, 1].min()),
        float(joined[:, 0].max()),
        float(joined[:, 1].max()),
    )


def _intersects(bounds: tuple[float, float, float, float], crop: pd.Series) -> bool:
    left, top, right, bottom = bounds
    x0, y0 = float(crop["x0_px"]), float(crop["y0_px"])
    x1 = x0 + float(crop["width_px"])
    y1 = y0 + float(crop["height_px"])
    return right >= x0 and left < x1 and bottom >= y0 and top < y1


def collect_crop_geometries(
    cell_geojson: Path,
    nucleus_geojson: Path,
    *,
    combined_geometry_file: bool,
    crops: pd.DataFrame,
) -> dict[str, list[dict[str, object]]]:
    """Stream GeoJSON and retain only geometries intersecting selected crops."""

    if crops.empty:
        return {}
    keys = [f"{row.pair_id}|{row.stratum}|{row.crop_rank}" for row in crops.itertuples()]
    records: dict[str, list[dict[str, object]]] = {key: [] for key in keys}
    crop_rows = list(crops.iterrows())
    selected_ids = set(crops["object_id"].astype(str))
    seen_selected: set[str] = set()

    def add_features(path: Path, *, role: str, use_embedded_nucleus: bool = False) -> None:
        for feature in iter_qupath_geojson_features(path):
            object_id = str(feature.get("id", ""))
            geometries = [(role, feature.get("geometry"))]
            if use_embedded_nucleus:
                geometries.append(("nucleus", feature.get("nucleusGeometry")))
            if role == "cell" and object_id in selected_ids:
                seen_selected.add(object_id)
            for geometry_role, geometry in geometries:
                bounds = _geometry_bounds(geometry)
                if bounds is None:
                    continue
                for position, (_, crop) in enumerate(crop_rows):
                    if _intersects(bounds, crop):
                        records[keys[position]].append(
                            {
                                "object_id": object_id,
                                "role": geometry_role,
                                "geometry": geometry,
                                "selected": object_id == str(crop["object_id"]),
                            }
                        )

    add_features(
        cell_geojson,
        role="cell",
        use_embedded_nucleus=combined_geometry_file,
    )
    if not combined_geometry_file:
        add_features(nucleus_geojson, role="nucleus")
    missing = sorted(selected_ids - seen_selected)
    if missing:
        raise ContractError(
            f"selected crop cells lack streamed cell geometry; examples={missing[:5]}"
        )
    return records


def _rings(geometry: Mapping[str, Any]):
    polygons = geometry["coordinates"]
    if geometry["type"] == "Polygon":
        polygons = [polygons]
    for polygon in polygons:
        for ring in polygon:
            array = np.asarray(ring, dtype=float)
            if array.ndim == 2 and array.shape[1] >= 2 and len(array) >= 2:
                yield array[:, :2]


def _display_channel(values: np.ndarray) -> np.ndarray:
    finite = values[np.isfinite(values)]
    if not len(finite):
        return np.zeros(values.shape, dtype=float)
    lower, upper = np.percentile(finite, [1.0, 99.5])
    if not np.isfinite(upper) or upper <= lower:
        return np.zeros(values.shape, dtype=float)
    return np.clip((values - lower) / (upper - lower), 0, 1)


def render_raw_crops(
    *,
    image,
    crops: pd.DataFrame,
    geometries: Mapping[str, Sequence[Mapping[str, object]]],
    output_root: Path,
) -> pd.DataFrame:
    """Render raw target/reference/DAPI panels using tiled level-0 reads."""

    if crops.empty:
        return crops.assign(render_path=pd.Series(dtype="string"))
    try:
        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise ContractError("crop rendering requires the notebook matplotlib extra") from exc

    output_root.mkdir(parents=True, exist_ok=True)
    rendered_paths: list[str] = []
    with tifffile.TiffFile(image.qptiff.path) as image_file:
        series = image_file.series[0]
        pages = series.pages
        channel_names = []
        for page in pages:
            description = page.description or ""
            match = re.search(r"<Biomarker>([^<]+)</Biomarker>", description) or re.search(
                r"<Name>([^<]+)</Name>", description
            )
            channel_names.append(match.group(1).strip() if match else f"ch{len(channel_names)}")
        markers = validated_channel_markers(image, channel_names)
        channel_lookup = {marker: index for index, marker in enumerate(markers)}
        required_markers = set(crops["target_display_marker"]).union(
            crops["reference_display_marker"], {"DAPI"}
        )
        missing = sorted(required_markers.difference(channel_lookup))
        if missing:
            raise ContractError(f"raw qptiff lacks requested crop markers {missing}")
        readers = {
            marker: TiledPageReader(image_file.filehandle, pages[channel_lookup[marker]])
            for marker in sorted(required_markers)
        }
        for crop in crops.itertuples(index=False):
            target_display = str(crop.target_display_marker)
            reference_display = str(crop.reference_display_marker)
            arrays = {
                marker: extract_tiled_window(
                    readers[marker],
                    x0=int(crop.x0_px),
                    y0=int(crop.y0_px),
                    width=int(crop.width_px),
                    height=int(crop.height_px),
                )
                for marker in (target_display, reference_display, "DAPI")
            }
            scaled = {marker: _display_channel(values) for marker, values in arrays.items()}
            composite = np.stack(
                [scaled[target_display], scaled[reference_display], scaled["DAPI"]], axis=-1
            )
            fig, axes = plt.subplots(1, 4, figsize=(12, 3.2), constrained_layout=True)
            axes[0].imshow(composite)
            axes[0].set_title(f"RGB: {target_display}/{reference_display}/DAPI")
            display_labels = (
                f"{crop.target} (raw {target_display})",
                f"{crop.reference} (raw {reference_display})",
                "DAPI",
            )
            for axis, marker, label in zip(
                axes[1:],
                (target_display, reference_display, "DAPI"),
                display_labels,
                strict=True,
            ):
                axis.imshow(scaled[marker], cmap="gray", vmin=0, vmax=1)
                axis.set_title(label)
            key = f"{crop.pair_id}|{crop.stratum}|{crop.crop_rank}"
            for axis in axes:
                for item in geometries.get(key, ()):
                    color = (
                        "#FFD700"
                        if item["selected"]
                        else ("#00FFFF" if item["role"] == "nucleus" else "#FFFFFF")
                    )
                    linewidth = 1.4 if item["selected"] else 0.45
                    for ring in _rings(item["geometry"]):
                        axis.plot(
                            ring[:, 0] - crop.x0_px,
                            ring[:, 1] - crop.y0_px,
                            color=color,
                            linewidth=linewidth,
                        )
                axis.set_xlim(0, crop.width_px)
                axis.set_ylim(crop.height_px, 0)
                axis.set_axis_off()
            fig.suptitle(
                f"{image.donor_id} · {crop.pair_id} · {crop.stratum} · {crop.object_id}",
                fontsize=9,
            )
            pair_slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", crop.pair_id)
            filename = f"{pair_slug}__{crop.stratum}__{int(crop.crop_rank):02d}.png"
            path = output_root / filename
            fig.savefig(path, dpi=150)
            plt.close(fig)
            rendered_paths.append(str(path.resolve()))
    result = crops.copy()
    result["render_path"] = rendered_paths
    return result


def _single_partition(root: Path, donor_id: str) -> Path:
    files = sorted((root / f"donor_id={donor_id}").glob("*.parquet"))
    if len(files) != 1:
        raise ContractError(
            f"expected one {root.name} parquet for donor {donor_id}, found {len(files)}"
        )
    return files[0]


def _candidate_frame(
    run_root: Path,
    projection_root: Path,
    donor_id: str,
    pairs: Sequence[MarkerPair],
) -> tuple[pd.DataFrame, ProjectionGeometryCoverage]:
    """Join exact upstream QuPath inputs to current geometry eligibility.

    The projection is produced from the full QuPath export by the shared
    development projection tool.  In particular, ``Iba1_norm`` is exactly
    ``Cell: Iba1: Mean (norm)``; neither raw Iba1 nor a RESTORE-derived value is
    an accepted substitute.  Projection/geometry coverage differences are
    reported while matched object IDs remain valid candidates.
    """

    markers = tuple(
        dict.fromkeys(
            ["DAPI", *(marker for pair in pairs for marker in (pair.target, pair.reference))]
        )
    )
    for marker in markers:
        selection_source_measurement(marker)
    try:
        qupath, _projection_paths, _projection_summary = load_qupath_projection_donor(
            projection_root,
            str(donor_id),
            [
                "object_id",
                "image",
                "X_centroid",
                "Y_centroid",
                "cell_area",
                "cell_solidity",
                *markers,
            ],
        )
    except (FileNotFoundError, ValueError) as exc:
        raise ContractError(str(exc)) from exc
    qupath = qupath[
        [
            "object_id",
            "image",
            "X_centroid",
            "Y_centroid",
            "cell_area",
            "cell_solidity",
            *markers,
        ]
    ].copy()
    geometry = pd.read_parquet(
        _single_partition(run_root / "geometry_qc", donor_id),
        columns=[
            "object_id",
            "analysis_eligible",
            "raster_area_ratio",
            "nucleus_cell_area_ratio",
        ],
    )
    for frame, label in (
        (qupath, "QuPath projection"),
        (geometry, "current geometry-QC"),
    ):
        object_ids = frame["object_id"]
        if object_ids.isna().any() or object_ids.astype("string").eq("").any():
            raise ContractError(f"{label} object IDs must be non-null and nonempty")
        if object_ids.astype("string").duplicated().any():
            raise ContractError(f"{label} object IDs must be unique")
        frame["object_id"] = object_ids.astype("string")

    joined = qupath.merge(
        geometry,
        on="object_id",
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    projection_only = sorted(joined.loc[joined["_merge"].eq("left_only"), "object_id"].astype(str))
    current_only = sorted(joined.loc[joined["_merge"].eq("right_only"), "object_id"].astype(str))
    result = joined.loc[joined["_merge"].eq("both")].drop(columns="_merge").copy()
    coverage = ProjectionGeometryCoverage(
        projection_cell_n=int(len(qupath)),
        current_cell_n=int(len(geometry)),
        matched_cell_n=int(len(result)),
        projection_only_n=len(projection_only),
        current_only_n=len(current_only),
        projection_only_examples=tuple(projection_only[:5]),
        current_only_examples=tuple(current_only[:5]),
    )
    return result, coverage


def _write_frame(frame: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    return path


def prepare_development_workspace(development_root: Path, label: str) -> Path:
    """Refresh only a boundary-audit-owned directory below a development root."""

    if Path(label).name != label or label in {"", ".", ".."}:
        raise ContractError("label must be one nonempty path component")
    root = development_root.resolve()
    output = root / label
    if output.is_symlink():
        raise ContractError(f"development workspace cannot be a symlink: {output}")
    if output.exists():
        marker = output / WORKSPACE_MARKER
        owned = marker.is_file()
        for summary_name in ("development_summary.json", "manifest.json"):
            summary_path = output / summary_name
            if summary_path.is_file():
                try:
                    payload = json.loads(summary_path.read_text(encoding="utf-8"))
                    owned = owned or str(payload.get("method", "")).startswith(
                        "boundary-contact-audit-"
                    )
                except (OSError, json.JSONDecodeError):
                    pass
        if not owned:
            raise ContractError(f"refusing to replace a workspace not owned by this tool: {output}")
        shutil.rmtree(output)
    output.mkdir(parents=True)
    (output / WORKSPACE_MARKER).write_text(
        json.dumps({"tool": "boundary_contact_audit", "status": "running"}) + "\n",
        encoding="utf-8",
    )
    return output


def complete_development_workspace(output: Path) -> None:
    (output / WORKSPACE_MARKER).write_text(
        json.dumps({"tool": "boundary_contact_audit", "status": "complete"}) + "\n",
        encoding="utf-8",
    )


def run_donor_audit(
    *,
    run_root: Path,
    projection_root: Path,
    scratch_root: Path,
    donor_id: str,
    image,
    pairs: Sequence[MarkerPair],
    output_root: Path,
    crops_per_stratum: int,
    crop_size_px: int,
) -> tuple[dict[str, object], list[Path]]:
    markers_for_crops = tuple(
        dict.fromkeys(
            ["DAPI", *(marker for pair in pairs for marker in (pair.target, pair.reference))]
        )
    )
    acquired = tuple(
        channel.marker for channel in sorted(image.channel_map, key=lambda x: x.channel_index)
    )
    display_markers = tuple(dict.fromkeys(display_marker(marker) for marker in markers_for_crops))
    missing = sorted(set(display_markers).difference(acquired))
    if missing:
        raise ContractError(f"donor {donor_id} has no acquired channels {missing}")

    # Crop evidence is intentionally prepared before and independently of scratch validation.
    with tifffile.TiffFile(image.qptiff.path) as image_file:
        shape = tuple(map(int, image_file.series[0].shape[-2:]))
    candidates, projection_geometry = _candidate_frame(run_root, projection_root, donor_id, pairs)
    if set(candidates["image"].astype(str)) != {image.image_id}:
        raise ContractError(f"donor {donor_id} cells declare a different image identity")
    crops = select_crop_manifest(
        candidates,
        donor_id=donor_id,
        image_id=image.image_id,
        pairs=pairs,
        image_shape=shape,
        pixel_size_um_x=image.pixel_size_um_x,
        pixel_size_um_y=image.pixel_size_um_y,
        crops_per_stratum=crops_per_stratum,
        crop_size_px=crop_size_px,
    )
    geometries = collect_crop_geometries(
        Path(image.cell_geojson.path),
        Path(image.nucleus_geojson.path),
        combined_geometry_file=image.combined_geometry_file,
        crops=crops,
    )
    crops = render_raw_crops(
        image=image,
        crops=crops,
        geometries=geometries,
        output_root=output_root / "crops" / f"donor_id={donor_id}",
    )
    crop_path = _write_frame(
        crops, output_root / "crop_manifest" / f"donor_id={donor_id}" / "data.parquet"
    )

    scratch_path = scratch_root / f"{donor_id}.npz"
    scratch = audit_scratch_file(scratch_path, candidates["object_id"])
    output_paths = [
        crop_path,
        *map(Path, crops.get("render_path", ())),
    ]
    projection_path = _single_partition(projection_root, donor_id)
    status = {
        "run_id": run_root.name,
        "donor_id": str(donor_id),
        "qupath_projection_path": str(projection_path.resolve()),
        "geometry_qc_path": str(_single_partition(run_root / "geometry_qc", donor_id).resolve()),
        "scratch_state": scratch.state,
        "scratch_reason": scratch.reason,
        "scratch_path": str(scratch_path.resolve()),
        "projection_cell_n": projection_geometry.projection_cell_n,
        "current_cell_n": projection_geometry.current_cell_n,
        "projection_geometry_matched_n": projection_geometry.matched_cell_n,
        "projection_only_n": projection_geometry.projection_only_n,
        "current_only_n": projection_geometry.current_only_n,
        "projection_only_examples": json.dumps(projection_geometry.projection_only_examples),
        "current_only_examples": json.dumps(projection_geometry.current_only_examples),
        "candidate_cell_n": int(len(candidates)),
        "scratch_cell_n": (
            int(len(scratch.validation.bundle.oids)) if scratch.validation is not None else 0
        ),
        "scratch_matched_cell_n": (
            int(len(scratch.validation.current_object_ids)) if scratch.validation is not None else 0
        ),
        "candidate_only_from_scratch_n": (
            scratch.validation.current_only_n
            if scratch.validation is not None
            else int(len(candidates))
        ),
        "scratch_only_n": (
            scratch.validation.scratch_only_n if scratch.validation is not None else 0
        ),
        "candidate_only_from_scratch_examples": (
            json.dumps(scratch.validation.current_only_examples)
            if scratch.validation is not None
            else "[]"
        ),
        "scratch_only_examples": (
            json.dumps(scratch.validation.scratch_only_examples)
            if scratch.validation is not None
            else "[]"
        ),
        "crop_n": int(len(crops)),
        "selected_crop_missing_from_scratch_n": (
            int((~crops["object_id"].astype(str).isin(scratch.validation.bundle.oids)).sum())
            if scratch.validation is not None and not crops.empty
            else int(len(crops))
            if scratch.validation is None
            else 0
        ),
    }
    if scratch.validation is not None:
        contacts_path = _write_frame(
            contact_metrics(scratch.validation),
            output_root / "contact_metrics" / f"donor_id={donor_id}" / "data.parquet",
        )
        selected_edges_path = _write_frame(
            selected_contact_edges(scratch.validation, crops["object_id"]),
            output_root / "selected_contacts" / f"donor_id={donor_id}" / "data.parquet",
        )
        output_paths.extend([contacts_path, selected_edges_path])
    return status, output_paths


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, help="current run directory name below data/runs")
    parser.add_argument("--donors", nargs="+", required=True)
    parser.add_argument(
        "--scratch-root",
        type=Path,
        required=True,
        help="explicit REDSEA intermediates directory (for this cohort: ../data/redsea_scratch/intermediates)",
    )
    parser.add_argument(
        "--label",
        default=DEFAULT_LABEL,
        help=f"mutable workspace name (default: {DEFAULT_LABEL})",
    )
    parser.add_argument(
        "--development-root",
        type=Path,
        default=DEFAULT_DEVELOPMENT_ROOT,
        help="development output root (default: outer ../data/development)",
    )
    parser.add_argument(
        "--pairs",
        nargs="+",
        help="explicit TARGET:REFERENCE pairs; omission uses the current roster",
    )
    parser.add_argument(
        "--roster",
        type=Path,
        default=DEFAULT_ROSTER,
        help="shared development roster used when --pairs is omitted",
    )
    parser.add_argument("--crops-per-stratum", type=int, default=4)
    parser.add_argument("--crop-size-px", type=int, default=256)
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_QUPATH_PROJECTION_ROOT,
        help="mutable upstream QuPath input projection",
    )
    parser.add_argument(
        "--refresh-inputs",
        action="store_true",
        help="refresh --input-root from the full QuPath exports before auditing",
    )
    parser.add_argument("--cohort-spec", type=Path, default=DEFAULT_COHORT_SPEC)
    parser.add_argument("--projection-block-size-mb", type=int, default=64)
    parser.add_argument(
        "--cohort-manifest", type=Path, default=REPOSITORY / "data/qupath_manifest.json"
    )
    args = parser.parse_args(argv)

    pairs = parse_pairs(args.pairs, roster_path=args.roster)
    run_root = REPOSITORY / "data" / "runs" / str(args.run)
    if not run_root.is_dir():
        parser.error(f"current run directory is missing: {run_root}")
    scratch_root = args.scratch_root.expanduser().resolve()
    if not scratch_root.is_dir():
        parser.error(f"--scratch-root is not a directory: {scratch_root}")

    try:
        images = load_live_images(args.cohort_manifest.expanduser().resolve())
    except ContractError as exc:
        parser.error(str(exc))
    image_by_donor = {image.donor_id: image for image in images}
    donors = tuple(dict.fromkeys(map(str, args.donors)))
    unknown = sorted(set(donors).difference(image_by_donor))
    if unknown:
        parser.error(f"donors absent from cohort manifest: {unknown}")
    projection_root = args.input_root.expanduser().resolve()
    if args.refresh_inputs:
        if args.projection_block_size_mb <= 0:
            parser.error("--projection-block-size-mb must be positive")
        try:
            refresh_qupath_projection(
                args.cohort_spec,
                projection_root,
                donors,
                block_size=args.projection_block_size_mb * 1024 * 1024,
            )
        except ValueError as exc:
            parser.error(str(exc))
    if not projection_root.is_dir():
        parser.error(
            f"upstream QuPath input projection is missing: {projection_root}; "
            "run with --refresh-inputs"
        )
    try:
        output_root = prepare_development_workspace(
            args.development_root.expanduser(), str(args.label)
        )
    except ContractError as exc:
        parser.error(str(exc))

    statuses: list[dict[str, object]] = []
    output_paths: list[Path] = []
    for donor_id in donors:
        status, donor_outputs = run_donor_audit(
            run_root=run_root,
            projection_root=projection_root,
            scratch_root=scratch_root,
            donor_id=donor_id,
            image=image_by_donor[donor_id],
            pairs=pairs,
            output_root=output_root,
            crops_per_stratum=args.crops_per_stratum,
            crop_size_px=args.crop_size_px,
        )
        statuses.append(status)
        output_paths.extend(donor_outputs)
        print(
            f"[{donor_id}] scratch={status['scratch_state']} crops={status['crop_n']}",
            flush=True,
        )
    status_path = _write_frame(
        pd.DataFrame.from_records(statuses), output_root / "audit_status.parquet"
    )
    output_paths.append(status_path)
    summary = {
        "method": "boundary-contact-audit-development-v4",
        "review_only": True,
        "workspace_semantics": "mutable; same label replaces this tool's prior outputs",
        "development_root": str(args.development_root.expanduser().resolve()),
        "label": str(args.label),
        "run_id": str(args.run),
        "run_root": str(run_root.resolve()),
        "cohort_manifest_path": str(args.cohort_manifest.expanduser().resolve()),
        "cohort_spec_path": str(args.cohort_spec.expanduser().resolve()),
        "qupath_input_projection_root": str(projection_root),
        "qupath_input_projection_refreshed": bool(args.refresh_inputs),
        "scratch_root": str(scratch_root),
        "pairs": [asdict(pair) for pair in pairs],
        "display_marker_overrides": DISPLAY_MARKER_OVERRIDES,
        "selection_source_measurements": {
            marker: SOURCE_MEASUREMENTS[marker]
            for marker in dict.fromkeys(
                marker for pair in pairs for marker in (pair.target, pair.reference)
            )
        },
        "selection_source_kinds": {
            marker: DIAGNOSTIC_SOURCE_KINDS[marker]
            for marker in dict.fromkeys(
                marker for pair in pairs for marker in (pair.target, pair.reference)
            )
        },
        "default_pair_roster": args.pairs is None,
        "selection": {
            "source": "mutable projection of exact upstream QuPath measurements plus current geometry QC",
            "restore_outputs_read": False,
            "projection_geometry_coverage": (
                "inner join exact object IDs; projection-only and current-only IDs are reported"
            ),
            "strata": list(STRATA),
            "crops_per_stratum": args.crops_per_stratum,
            "crop_size_px": args.crop_size_px,
            "tie_break": "lexicographic object_id",
        },
        "display": {
            "source": "raw qptiff channels",
            "selection_values_used_for_display": False,
        },
        "projection_geometry_contract": {
            "authority": "exact object ID inner join",
            "coverage_mismatch": "reported; matched candidate rows are retained",
            "invalid_ids": "null, empty, or duplicate IDs fail",
        },
        "scratch_contract": {
            "authority": "scratch object IDs and contact CSR",
            "required_arrays": sorted(REQUIRED_SCRATCH_KEYS),
            "coverage_mismatch": "reported; does not suppress matched contact rows",
            "compartment_arrays": "not read from scratch and not required",
        },
        "donors": statuses,
        "outputs": [str(path.resolve()) for path in sorted(output_paths)],
    }
    summary_path = output_root / "development_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    complete_development_workspace(output_root)
    print(f"wrote review audit: {output_root}", flush=True)
    return 0


if __name__ == "__main__":
    # Matplotlib cannot write its login-home cache in the managed environment.
    os.environ.setdefault(
        "MPLCONFIGDIR",
        str(Path(tempfile.gettempdir()) / "phenocycler-matplotlib-cache"),
    )
    raise SystemExit(main())
