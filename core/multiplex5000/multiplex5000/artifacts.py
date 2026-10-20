"""Input contracts for the Phenocycler workflow, and the hashing helpers behind them.

CORRECTED 2026-09-18.  This docstring used to describe a content-addressed regime that no longer
exists: the run framework was REMOVED on 2026-09-16, and with it the stage and run manifests, the
dataset snapshots and the source-code fingerprints.  Presence of output IS the completion check
now -- see ``pipeline.run_stage``, which states that nothing validates a fingerprint or a manifest.

What actually remains here:

* an explicit QuPath image and cohort contract (measurements, image, cell/nucleus geometry,
  pixel calibration, channel map and segmentation version) -- :class:`QuPathImageManifest` and
  :class:`QuPathCohortManifest`;
* streaming GeoJSON inspection and the geometry/schema contract errors those raise; and
* the hashing helpers (:func:`sha256_file`, :func:`sha256_json`, :func:`hash_identifiers`,
  :func:`dataframe_content_sha256`) the cohort manifest is written with.

The cohort manifest's ``content_id`` and per-file digests are DESCRIPTIVE.  Nothing validates
them -- ``PanelChannelMap.verify_identity`` says so in one line and returns -- which is
deliberate: the manifest is a description of where the raw inputs are, so it can be hand-edited
when a source file moves.

The module intentionally has no dependency on ``pipeline.py``.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
from dataclasses import asdict, dataclass, fields, is_dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


MANIFEST_VERSION = 1
_HASH_CHUNK_BYTES = 8 * 1024 * 1024
_SAMPLE_BYTES = 1024 * 1024
_SAMPLE_COUNT = 3


class ContractError(ValueError):
    """Base class for a violated input or artifact contract."""


class GeometryContractError(ContractError):
    """A QuPath geometry export is absent, incomplete or malformed."""


class DuplicateObjectError(GeometryContractError):
    """A supposedly canonical object-id universe contains duplicates."""


class PartialArtifactError(ContractError):
    """The exact expected donor set is not complete."""


class SchemaContractError(ContractError):
    """Parquet files do not carry one exact, stable schema."""


def _normalise_json(value: Any) -> Any:
    """Convert common Python/scientific values to deterministic JSON values."""
    if is_dataclass(value):
        return _normalise_json(asdict(value))
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {
            str(key): _normalise_json(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_normalise_json(item) for item in value]
    if isinstance(value, (set, frozenset)):
        normalised = [_normalise_json(item) for item in value]
        return sorted(normalised, key=canonical_json)
    if hasattr(value, "item") and callable(value.item):
        try:
            return _normalise_json(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            return {"__nonfinite_float__": "nan"}
        return {"__nonfinite_float__": "inf" if value > 0 else "-inf"}
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def canonical_json(value: Any) -> str:
    """Return a deterministic, whitespace-free JSON representation."""
    return json.dumps(
        _normalise_json(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def sha256_json(value: Any) -> str:
    """Hash a value after deterministic JSON normalisation."""
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def _content_identity(value: Any) -> Any:
    """Remove validation-cache/timestamp fields from a content identifier."""
    normalised = _normalise_json(value)
    if isinstance(normalised, dict):
        return {
            key: _content_identity(item)
            for key, item in normalised.items()
            if key not in {"created_at_utc", "mtime_ns"}
        }
    if isinstance(normalised, list):
        return [_content_identity(item) for item in normalised]
    return normalised


def sha256_file(path: Path | str) -> str:
    """Hash all bytes in ``path`` (no mtime/size shortcut)."""
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_HASH_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_sampled_file(
    path: Path | str,
    *,
    sample_size_bytes: int = _SAMPLE_BYTES,
    sample_count: int = _SAMPLE_COUNT,
) -> str:
    """Hash deterministic windows across a large file.

    This is an explicitly weaker, opt-in fingerprint for environments where an
    initial full read of a multi-hundred-GB qptiff is infeasible. The file size,
    sample layout and offsets are included in the digest so the result cannot be
    confused with a full SHA-256.
    """
    path = Path(path)
    if sample_size_bytes <= 0 or sample_count <= 0:
        raise ContractError("sample size and count must be positive")
    size = path.stat().st_size
    width = min(sample_size_bytes, size)
    if size <= width or sample_count == 1:
        offsets = (0,)
    else:
        last = size - width
        offsets = tuple(
            sorted({round(index * last / (sample_count - 1)) for index in range(sample_count)})
        )
    digest = hashlib.sha256()
    digest.update(b"phenocycler-sampled-sha256-v1\0")
    digest.update(size.to_bytes(16, "big"))
    digest.update(sample_size_bytes.to_bytes(8, "big"))
    digest.update(sample_count.to_bytes(8, "big"))
    with path.open("rb") as handle:
        for offset in offsets:
            handle.seek(offset)
            chunk = handle.read(width)
            digest.update(offset.to_bytes(16, "big"))
            digest.update(len(chunk).to_bytes(8, "big"))
            digest.update(chunk)
    return digest.hexdigest()


def hash_identifiers(values: Iterable[Any]) -> str:
    """Hash an exact identifier set independent of row/file order.

    Length-prefixing prevents ambiguous concatenations (``["ab", "c"]`` versus
    ``["a", "bc"]``). Nulls and duplicates are contract violations.
    """
    ids: list[str] = []
    for value in values:
        if value is None:
            raise ContractError("object_id contains null")
        text = str(value)
        if not text:
            raise ContractError("object_id contains an empty value")
        ids.append(text)
    ids.sort()
    for previous, current in zip(ids, ids[1:]):
        if previous == current:
            raise DuplicateObjectError(f"duplicate object_id: {current}")
    digest = hashlib.sha256()
    for value in ids:
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def dataframe_content_sha256(frame, *, sort_by: Sequence[str] = ()) -> str:
    """Hash tabular values deterministically without relying on pickle internals."""
    missing = set(sort_by).difference(frame.columns)
    if missing:
        raise ContractError(f"cannot hash frame: missing sort columns {sorted(missing)}")
    ordered = frame.copy()
    if sort_by:
        ordered = ordered.sort_values(list(sort_by), kind="mergesort")
    rows = []
    for row in ordered.itertuples(index=False, name=None):
        values = []
        for value in row:
            try:
                missing_value = bool(value is None or value != value)
            except (TypeError, ValueError):
                missing_value = False
            values.append(None if missing_value else _normalise_json(value))
        rows.append(values)
    return sha256_json({"columns": list(map(str, ordered.columns)), "rows": rows})


@dataclass(frozen=True)
class FileFingerprint:
    """An exact path + byte-content fingerprint."""

    path: str
    size_bytes: int
    mtime_ns: int
    sha256: str
    hash_mode: str = "full"
    sample_size_bytes: int = 0
    sample_count: int = 0

    @classmethod
    def capture(
        cls,
        path: Path | str,
        *,
        hash_mode: str = "full",
        sample_size_bytes: int = _SAMPLE_BYTES,
        sample_count: int = _SAMPLE_COUNT,
    ) -> "FileFingerprint":
        resolved = Path(path).expanduser().resolve()
        if not resolved.is_file():
            raise ContractError(f"required file does not exist: {resolved}")
        if hash_mode not in {"full", "sampled"}:
            raise ContractError("hash_mode must be 'full' or 'sampled'")
        stat = resolved.stat()
        digest = (
            sha256_file(resolved)
            if hash_mode == "full"
            else sha256_sampled_file(
                resolved,
                sample_size_bytes=sample_size_bytes,
                sample_count=sample_count,
            )
        )
        return cls(
            path=str(resolved),
            size_bytes=stat.st_size,
            mtime_ns=stat.st_mtime_ns,
            sha256=digest,
            hash_mode=hash_mode,
            sample_size_bytes=sample_size_bytes if hash_mode == "sampled" else 0,
            sample_count=sample_count if hash_mode == "sampled" else 0,
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "FileFingerprint":
        return cls(
            path=str(value["path"]),
            size_bytes=int(value["size_bytes"]),
            mtime_ns=int(value["mtime_ns"]),
            sha256=str(value["sha256"]),
            hash_mode=str(value.get("hash_mode", "full")),
            sample_size_bytes=int(value.get("sample_size_bytes", 0)),
            sample_count=int(value.get("sample_count", 0)),
        )


@dataclass(frozen=True)
class ChannelMapping:
    """One exact marker-to-qptiff-channel mapping."""

    marker: str
    channel_name: str
    channel_index: int

    def __post_init__(self) -> None:
        if not self.marker.strip():
            raise ContractError("channel-map marker cannot be empty")
        if not self.channel_name.strip():
            raise ContractError(f"{self.marker}: channel name cannot be empty")
        if self.channel_index < 0:
            raise ContractError(f"{self.marker}: channel index must be >= 0")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ChannelMapping":
        return cls(
            marker=str(value["marker"]),
            channel_name=str(value["channel_name"]),
            channel_index=int(value["channel_index"]),
        )


@dataclass(frozen=True)
class GeoJSONSummary:
    feature_count: int
    primary_geometry_count: int
    nucleus_geometry_count: int
    object_id_sha256: str

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "GeoJSONSummary":
        return cls(
            feature_count=int(value["feature_count"]),
            primary_geometry_count=int(value["primary_geometry_count"]),
            nucleus_geometry_count=int(value["nucleus_geometry_count"]),
            object_id_sha256=str(value["object_id_sha256"]),
        )


def geometry_present(geometry: Any) -> bool:
    if not isinstance(geometry, Mapping) or not geometry.get("type"):
        return False
    if geometry.get("type") == "GeometryCollection":
        return bool(geometry.get("geometries"))
    return bool(geometry.get("coordinates"))


class _JSONStream:
    """Incrementally decode JSON values without retaining the complete file."""

    def __init__(self, path: Path, *, chunk_size: int = 1024 * 1024):
        self.path = path
        self.handle = path.open(encoding="utf-8")
        self.chunk_size = chunk_size
        self.buffer = ""
        self.position = 0
        self.eof = False
        self.decoder = json.JSONDecoder()

    def close(self) -> None:
        self.handle.close()

    def _read_more(self) -> bool:
        if self.eof:
            return False
        if self.position:
            self.buffer = self.buffer[self.position :]
            self.position = 0
        chunk = self.handle.read(self.chunk_size)
        if not chunk:
            self.eof = True
            return False
        self.buffer += chunk
        return True

    def _skip_whitespace(self) -> None:
        while True:
            while (
                self.position < len(self.buffer)
                and self.buffer[self.position].isspace()
            ):
                self.position += 1
            if self.position < len(self.buffer) or not self._read_more():
                return

    def peek(self) -> str:
        self._skip_whitespace()
        if self.position >= len(self.buffer):
            raise json.JSONDecodeError(
                "unexpected end of JSON", self.buffer, self.position
            )
        return self.buffer[self.position]

    def expect(self, token: str) -> None:
        if self.peek() != token:
            raise json.JSONDecodeError(
                f"expected {token!r}", self.buffer, self.position
            )
        self.position += 1

    def decode(self) -> Any:
        self._skip_whitespace()
        while True:
            try:
                value, end = self.decoder.raw_decode(self.buffer, self.position)
            except json.JSONDecodeError:
                if not self._read_more():
                    raise
            else:
                self.position = end
                return value

    def require_end(self) -> None:
        self._skip_whitespace()
        if self.position < len(self.buffer) or self._read_more():
            self._skip_whitespace()
            if self.position < len(self.buffer):
                raise json.JSONDecodeError(
                    "extra data", self.buffer, self.position
                )


def iter_qupath_geojson_features(path: Path | str) -> Iterator[Mapping[str, Any]]:
    """Stream one QuPath FeatureCollection feature at a time.

    Cell exports in this cohort are several gigabytes each.  Using
    :func:`json.load` expands the full coordinate payload in memory even though
    manifest inspection and rasterization need only one feature at a time.
    """

    path = Path(path)
    if not path.is_file():
        raise GeometryContractError(f"GeoJSON is missing: {path}")
    stream = _JSONStream(path)
    top_level_type: Any = None
    found_features = False
    try:
        stream.expect("{")
        first_member = True
        while stream.peek() != "}":
            if not first_member:
                stream.expect(",")
            key = stream.decode()
            if not isinstance(key, str):
                raise GeometryContractError(
                    f"{path}: top-level GeoJSON key is not a string"
                )
            stream.expect(":")
            if key == "type":
                top_level_type = stream.decode()
            elif key == "features":
                if found_features:
                    raise GeometryContractError(
                        f"{path}: duplicate top-level features member"
                    )
                found_features = True
                stream.expect("[")
                first_feature = True
                while stream.peek() != "]":
                    if not first_feature:
                        stream.expect(",")
                    feature = stream.decode()
                    if not isinstance(feature, Mapping):
                        raise GeometryContractError(
                            f"{path}: feature is not an object"
                        )
                    yield feature
                    first_feature = False
                stream.expect("]")
            else:
                stream.decode()
            first_member = False
        stream.expect("}")
        stream.require_end()
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise GeometryContractError(f"cannot read GeoJSON {path}: {exc}") from exc
    finally:
        stream.close()
    if top_level_type != "FeatureCollection" or not found_features:
        raise GeometryContractError(f"{path}: expected a GeoJSON FeatureCollection")


def inspect_qupath_geojson(
    path: Path | str,
    *,
    require_primary_geometry: bool = True,
    require_nucleus_geometry: bool = False,
) -> tuple[GeoJSONSummary, tuple[str, ...]]:
    """Validate a QuPath FeatureCollection and return its exact object universe."""
    path = Path(path)
    object_ids: list[str] = []
    primary_count = 0
    nucleus_count = 0
    for index, feature in enumerate(iter_qupath_geojson_features(path)):
        object_id = feature.get("id")
        if object_id is None or not str(object_id):
            raise GeometryContractError(f"{path}: feature {index} has no object id")
        object_ids.append(str(object_id))
        primary = geometry_present(feature.get("geometry"))
        nucleus = geometry_present(feature.get("nucleusGeometry"))
        primary_count += int(primary)
        nucleus_count += int(nucleus)
        if require_primary_geometry and not primary:
            raise GeometryContractError(
                f"{path}: object {object_id} has no usable cell geometry"
            )
        if require_nucleus_geometry and not nucleus:
            raise GeometryContractError(
                f"{path}: object {object_id} has no usable nucleusGeometry"
            )
    object_hash = hash_identifiers(object_ids)
    return (
        GeoJSONSummary(
            feature_count=len(object_ids),
            primary_geometry_count=primary_count,
            nucleus_geometry_count=nucleus_count,
            object_id_sha256=object_hash,
        ),
        tuple(object_ids),
    )



def measurement_delimiter(path) -> str:
    """Field separator for a QuPath measurement export, from its suffix.

    QuPath writes tab-separated when the file is named .tsv and comma-separated for .csv, and
    the 2026-09 re-exports under data/RevisedQuPathExports are .tsv. `csv.reader` defaults to
    a comma, so a .tsv header parsed as one 994-character field and the manifest failed with
    "missing required column 'Image'" -- the column was there, the delimiter was not.

    Suffix, not a sniff: csv.Sniffer on a 994-column header is slow and can pick the wrong
    separator when field text contains commas, and both export formats are unambiguous by name.
    The bulk row read in cells_parquet.py goes through DuckDB's read_csv, which auto-detects and
    already handled both (verified: 994 columns, 'Image' present).
    """
    return "\t" if str(path).lower().endswith((".tsv", ".tab")) else ","


def is_parquet_measurement(path) -> bool:
    """A v3 QuPath export is Parquet; v2 exports are CSV/TSV. Decided by suffix, as above."""
    return str(path).lower().endswith(".parquet")


def measurement_header(path) -> tuple[str, ...]:
    """Column names of one QuPath measurement export: Parquet schema, else the CSV/TSV header.

    The Parquet schema is read from the footer only (``pyarrow.parquet.read_schema``); no row
    is read. QuPath's own column names are kept verbatim, micro sign included.
    """
    if is_parquet_measurement(path):
        import pyarrow.parquet as pq

        return tuple(pq.read_schema(str(path)).names)
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        try:
            return tuple(next(csv.reader(handle, delimiter=measurement_delimiter(path))))
        except StopIteration as exc:
            raise ContractError(f"QuPath measurement export {path} is empty") from exc


#: What the v3 per-image context Parquet must carry to stand in for the cells file's Parent,
#: centroids and signed distances. The cells file is written before islets are parented, so
#: its own Parent is not authoritative when a context file is declared.
CONTEXT_REQUIRED_COLUMNS = ("Object ID", "Parent", "Centroid X µm", "Centroid Y µm")

@dataclass(frozen=True)
class QuPathImageManifest:
    """The complete immutable contract for one QuPath image export."""

    manifest_version: int
    donor_id: str
    image_id: str
    measurement_csv: FileFingerprint
    qptiff: FileFingerprint
    cell_geojson: FileFingerprint
    nucleus_geojson: FileFingerprint
    pixel_size_um_x: float
    pixel_size_um_y: float
    panel_id: str
    channel_map: tuple[ChannelMapping, ...]
    segmentation_version: str
    measurement_columns: tuple[str, ...]
    cell_geometry: GeoJSONSummary
    nucleus_geometry: GeoJSONSummary
    combined_geometry_file: bool
    content_id: str
    #: v3 exports split one image into a cells Parquet (``measurement_csv``) and a context
    #: Parquet joined to it on ``Object ID``. ``None`` for a v2 CSV/TSV entry, and then absent
    #: from the serialised entry, so a v2 manifest's shape and ``content_id`` are unchanged.
    measurement_context: FileFingerprint | None = None

    @classmethod
    def create(
        cls,
        *,
        donor_id: str,
        image_id: str,
        measurement_csv: Path | str,
        qptiff: Path | str,
        cell_geojson: Path | str,
        nucleus_geojson: Path | str | None,
        pixel_size_um_x: float,
        pixel_size_um_y: float,
        panel_id: str,
        channel_map: Sequence[ChannelMapping],
        segmentation_version: str,
        qptiff_hash_mode: str = "full",
        qptiff_sample_size_bytes: int = _SAMPLE_BYTES,
        qptiff_sample_count: int = _SAMPLE_COUNT,
        measurement_fingerprint: FileFingerprint | None = None,
        qptiff_fingerprint: FileFingerprint | None = None,
        measurement_context: Path | str | None = None,
        context_fingerprint: FileFingerprint | None = None,
    ) -> "QuPathImageManifest":
        donor_id = str(donor_id).strip()
        image_id = str(image_id).strip()
        if not donor_id or not image_id:
            raise ContractError("donor_id and image_id are required")
        if pixel_size_um_x <= 0 or pixel_size_um_y <= 0:
            raise ContractError("pixel calibration must be finite and positive")
        if not math.isfinite(pixel_size_um_x) or not math.isfinite(pixel_size_um_y):
            raise ContractError("pixel calibration must be finite and positive")
        if not panel_id.strip() or not segmentation_version.strip():
            raise ContractError("panel_id and segmentation_version are required")
        channels = tuple(channel_map)
        if not channels:
            raise ContractError("channel_map cannot be empty")
        markers = [channel.marker for channel in channels]
        names = [channel.channel_name for channel in channels]
        indices = [channel.channel_index for channel in channels]
        if len(set(markers)) != len(markers):
            raise ContractError("channel_map contains duplicate markers")
        if len(set(names)) != len(names) or len(set(indices)) != len(indices):
            raise ContractError("channel_map maps more than one marker to a channel")

        measurement_path = Path(measurement_csv).expanduser().resolve()
        if measurement_fingerprint is None:
            measurement_fp = FileFingerprint.capture(measurement_path)
        else:
            measurement_fp = measurement_fingerprint
            if Path(measurement_fp.path) != measurement_path:
                raise ContractError(
                    "precomputed measurement fingerprint path differs from "
                    f"measurement_csv: {measurement_fp.path} != {measurement_path}"
                )
        measurement_columns = measurement_header(measurement_fp.path)
        if not measurement_columns:
            raise ContractError("QuPath measurement CSV is empty")
        context_fp: FileFingerprint | None = None
        if measurement_context is not None:
            if not is_parquet_measurement(measurement_fp.path):
                raise ContractError(
                    "measurement_context pairs only with a cells Parquet, not "
                    f"{measurement_fp.path}"
                )
            context_path = Path(measurement_context).expanduser().resolve()
            if not is_parquet_measurement(context_path):
                raise ContractError(f"measurement_context must be Parquet: {context_path}")
            if context_fingerprint is None:
                context_fp = FileFingerprint.capture(context_path)
            else:
                context_fp = context_fingerprint
                if Path(context_fp.path) != context_path:
                    raise ContractError(
                        "precomputed context fingerprint path differs from "
                        f"measurement_context: {context_fp.path} != {context_path}"
                    )
            context_columns = measurement_header(context_fp.path)
            absent = [c for c in CONTEXT_REQUIRED_COLUMNS if c not in context_columns]
            if absent:
                raise ContractError(
                    f"QuPath context Parquet {context_fp.path} is missing {absent}"
                )
            measurement_columns = measurement_columns + tuple(
                c for c in context_columns if c not in measurement_columns
            )
        elif context_fingerprint is not None:
            raise ContractError("context_fingerprint given without measurement_context")
        for required in ("Image", "Object ID"):
            if required not in measurement_columns:
                raise ContractError(
                    f"QuPath measurement CSV is missing required column {required!r}"
                )

        cell_path = Path(cell_geojson).expanduser().resolve()
        nucleus_path = (
            Path(nucleus_geojson).expanduser().resolve()
            if nucleus_geojson is not None
            else cell_path
        )
        combined = cell_path == nucleus_path
        if combined:
            cell_summary, cell_ids = inspect_qupath_geojson(
                cell_path,
                require_primary_geometry=True,
            )
            if cell_summary.nucleus_geometry_count == 0:
                raise GeometryContractError(
                    f"{cell_path}: no feature has a usable nucleusGeometry"
                )
            nucleus_summary = GeoJSONSummary(
                feature_count=cell_summary.feature_count,
                primary_geometry_count=cell_summary.nucleus_geometry_count,
                nucleus_geometry_count=cell_summary.nucleus_geometry_count,
                object_id_sha256=cell_summary.object_id_sha256,
            )
        else:
            cell_summary, cell_ids = inspect_qupath_geojson(
                cell_path, require_primary_geometry=True
            )
            nucleus_summary, nucleus_ids = inspect_qupath_geojson(
                nucleus_path, require_primary_geometry=True
            )
            if set(cell_ids) != set(nucleus_ids):
                raise GeometryContractError(
                    "cell and nucleus GeoJSON object-id universes differ"
                )
        draft = cls(
            manifest_version=MANIFEST_VERSION,
            donor_id=donor_id,
            image_id=image_id,
            measurement_csv=measurement_fp,
            qptiff=(
                FileFingerprint.capture(
                    qptiff,
                    hash_mode=qptiff_hash_mode,
                    sample_size_bytes=qptiff_sample_size_bytes,
                    sample_count=qptiff_sample_count,
                )
                if qptiff_fingerprint is None
                else qptiff_fingerprint
            ),
            cell_geojson=FileFingerprint.capture(cell_path),
            nucleus_geojson=FileFingerprint.capture(nucleus_path),
            pixel_size_um_x=float(pixel_size_um_x),
            pixel_size_um_y=float(pixel_size_um_y),
            panel_id=panel_id.strip(),
            channel_map=tuple(sorted(channels, key=lambda item: item.channel_index)),
            segmentation_version=segmentation_version.strip(),
            measurement_columns=measurement_columns,
            cell_geometry=cell_summary,
            nucleus_geometry=nucleus_summary,
            combined_geometry_file=combined,
            content_id="",
            measurement_context=context_fp,
        )
        if Path(draft.qptiff.path) != Path(qptiff).expanduser().resolve():
            raise ContractError(
                "precomputed qptiff fingerprint path differs from qptiff: "
                f"{draft.qptiff.path} != {Path(qptiff).expanduser().resolve()}"
            )
        return replace(draft, content_id=sha256_json(draft._identity_payload()))

    def _identity_payload(self) -> dict[str, Any]:
        value = self.to_dict()
        value.pop("content_id", None)
        return _content_identity(value)

    def to_dict(self) -> dict[str, Any]:
        value = _normalise_json(asdict(self))
        if value.get("measurement_context") is None:
            value.pop("measurement_context", None)
        return value

    def verify_identity(self) -> None:
        if self.manifest_version != MANIFEST_VERSION:
            raise ContractError(
                f"unsupported QuPath manifest version {self.manifest_version}"
            )

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "QuPathImageManifest":
        result = cls(
            manifest_version=int(value["manifest_version"]),
            donor_id=str(value["donor_id"]),
            image_id=str(value["image_id"]),
            measurement_csv=FileFingerprint.from_dict(value["measurement_csv"]),
            qptiff=FileFingerprint.from_dict(value["qptiff"]),
            cell_geojson=FileFingerprint.from_dict(value["cell_geojson"]),
            nucleus_geojson=FileFingerprint.from_dict(value["nucleus_geojson"]),
            pixel_size_um_x=float(value["pixel_size_um_x"]),
            pixel_size_um_y=float(value["pixel_size_um_y"]),
            panel_id=str(value["panel_id"]),
            channel_map=tuple(
                ChannelMapping.from_dict(item) for item in value["channel_map"]
            ),
            segmentation_version=str(value["segmentation_version"]),
            measurement_columns=tuple(map(str, value["measurement_columns"])),
            cell_geometry=GeoJSONSummary.from_dict(value["cell_geometry"]),
            nucleus_geometry=GeoJSONSummary.from_dict(value["nucleus_geometry"]),
            combined_geometry_file=bool(value["combined_geometry_file"]),
            content_id=str(value["content_id"]),
            measurement_context=(
                None
                if value.get("measurement_context") is None
                else FileFingerprint.from_dict(value["measurement_context"])
            ),
        )
        result.verify_identity()
        return result

    @classmethod
    def read_json(cls, path: Path | str) -> "QuPathImageManifest":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))


@dataclass(frozen=True)
class PanelChannelMap:
    """The one exact channel mapping permitted for a named acquisition panel."""

    panel_id: str
    channels: tuple[ChannelMapping, ...]
    content_id: str

    @classmethod
    def create(
        cls, panel_id: str, channels: Sequence[ChannelMapping]
    ) -> "PanelChannelMap":
        panel_id = str(panel_id).strip()
        if not panel_id:
            raise ContractError("panel_id cannot be empty")
        ordered = tuple(sorted(channels, key=lambda item: item.channel_index))
        if not ordered:
            raise ContractError(f"panel {panel_id!r} has no channels")
        markers = [channel.marker for channel in ordered]
        names = [channel.channel_name for channel in ordered]
        indices = [channel.channel_index for channel in ordered]
        if len(set(markers)) != len(markers):
            raise ContractError(f"panel {panel_id!r} contains duplicate markers")
        if len(set(names)) != len(names) or len(set(indices)) != len(indices):
            raise ContractError(
                f"panel {panel_id!r} maps multiple markers to one channel"
            )
        payload = {
            "panel_id": panel_id,
            "channels": [asdict(channel) for channel in ordered],
        }
        return cls(
            panel_id=panel_id,
            channels=ordered,
            content_id=sha256_json(payload),
        )

    def verify_identity(self) -> None:
        """Structural only. The recorded content_id is descriptive, not a gate."""
        return None

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PanelChannelMap":
        result = cls(
            panel_id=str(value["panel_id"]),
            channels=tuple(
                ChannelMapping.from_dict(item) for item in value["channels"]
            ),
            content_id=str(value["content_id"]),
        )
        result.verify_identity()
        return result


def _cohort_panel_maps(
    images: Sequence[QuPathImageManifest],
) -> tuple[PanelChannelMap, ...]:
    by_panel: dict[str, tuple[ChannelMapping, ...]] = {}
    for image in images:
        mapping = tuple(
            sorted(image.channel_map, key=lambda item: item.channel_index)
        )
        previous = by_panel.get(image.panel_id)
        if previous is not None and previous != mapping:
            raise ContractError(
                f"panel {image.panel_id!r} has contradictory channel mappings "
                f"between QuPath images"
            )
        by_panel[image.panel_id] = mapping
    return tuple(
        PanelChannelMap.create(panel_id, by_panel[panel_id])
        for panel_id in sorted(by_panel)
    )


def _validate_cohort_images(
    images: Sequence[QuPathImageManifest],
    expected_donors: Sequence[str],
) -> tuple[tuple[QuPathImageManifest, ...], tuple[str, ...]]:
    records = tuple(
        sorted(images, key=lambda item: (item.donor_id, item.image_id))
    )
    if not records:
        raise ContractError("a QuPath cohort manifest requires at least one image")
    for image in records:
        image.verify_identity()

    donor_ids = [image.donor_id for image in records]
    image_ids = [image.image_id for image in records]
    duplicate_donors = sorted(
        donor for donor in set(donor_ids) if donor_ids.count(donor) > 1
    )
    duplicate_images = sorted(
        image for image in set(image_ids) if image_ids.count(image) > 1
    )
    if duplicate_donors:
        raise ContractError(
            "QuPath cohort requires exactly one image per donor_id; duplicate donor_id="
            f"{duplicate_donors}. Note that donor_id is the SECTION key, not the donor: a "
            "donor contributes a pancreas and 1-4 lymph nodes, and each is its own image, "
            "its own partition and its own coordinate frame. Key them the way "
            "phenocycler.sections does -- '6539' (pancreas), '6539pln' (pln_1), '6539pln2' "
            "(pln_2) -- rather than repeating the bare donor id. Two images sharing a "
            "partition would be masked with one another's polygons and normalised as one "
            "image, with nothing raised downstream."
        )
    if duplicate_images:
        raise ContractError(
            f"QuPath cohort contains duplicate image_id values={duplicate_images}"
        )

    expected = tuple(sorted(map(str, expected_donors)))
    if not expected or len(set(expected)) != len(expected):
        raise ContractError("expected_donors must be a non-empty unique set")
    actual = set(donor_ids)
    expected_set = set(expected)
    missing = sorted(expected_set - actual)
    extra = sorted(actual - expected_set)
    if missing or extra:
        raise PartialArtifactError(
            f"QuPath cohort donor set differs: missing={missing}, extra={extra}"
        )
    return records, expected


@dataclass(frozen=True)
class QuPathCohortManifest:
    """Exact cohort-level collection of one QuPath image per donor.

    Image manifests are embedded rather than linked indirectly, so the cohort
    content identifier owns every source digest and geometry contract. Images
    sharing a ``panel_id`` must carry byte-for-byte equivalent marker, channel
    name and channel-index mappings. Different named panels may intentionally
    differ.
    """

    manifest_version: int
    cohort_name: str
    expected_donors: tuple[str, ...]
    images: tuple[QuPathImageManifest, ...]
    panels: tuple[PanelChannelMap, ...]
    created_at_utc: str
    content_id: str

    @classmethod
    def create(
        cls,
        *,
        cohort_name: str,
        expected_donors: Sequence[str],
        images: Sequence[QuPathImageManifest],
        validate_sources: bool = True,
        validation_mode: str = "fast",
    ) -> "QuPathCohortManifest":
        cohort_name = str(cohort_name).strip()
        if not cohort_name:
            raise ContractError("cohort_name is required")
        records, expected = _validate_cohort_images(images, expected_donors)
        panels = _cohort_panel_maps(records)
        draft = cls(
            manifest_version=MANIFEST_VERSION,
            cohort_name=cohort_name,
            expected_donors=expected,
            images=records,
            panels=panels,
            created_at_utc=datetime.now(timezone.utc).isoformat(),
            content_id="",
        )
        return replace(draft, content_id=sha256_json(draft._identity_payload()))

    @property
    def donor_ids(self) -> tuple[str, ...]:
        return tuple(image.donor_id for image in self.images)

    @property
    def image_ids(self) -> tuple[str, ...]:
        return tuple(image.image_id for image in self.images)

    @property
    def section_keys(self) -> tuple[str, ...]:
        """The partition keys, under their honest name.

        Identical to :attr:`donor_ids`; both return the ``donor_id`` field verbatim. The
        field is named for what it held when a cohort was one image per donor, but the
        value is a *section* key -- one image, one GeoJSON, one coordinate frame. A donor
        with a pancreas and two lymph nodes contributes ``6539``, ``6539pln``, ``6539pln2``.
        Renaming the field would change ``content_id`` and therefore every ``run_id``, so
        the name stays and this alias says what it means.
        """
        return self.donor_ids

    def donors(self) -> tuple[str, ...]:
        """The unique donors behind the sections, in sorted order.

        This is the thing to join metadata and the Xenium pairing manifest on. Sections are
        the work unit; donors are the cohort. Conflating them is how ``6539pln`` becomes a
        phantom donor with no disease status and no Xenium counterpart.
        """
        from .sections import SectionError, parse

        out = set()
        for key in self.donor_ids:
            try:
                out.add(parse(key).donor_id)
            except SectionError:
                # A key that names no section is its own donor. Only that failure is
                # tolerable here; anything else `parse` raises is a real defect and must
                # not be turned into a plausible-looking donor id.
                out.add(key)
        return tuple(sorted(out))

    def image_to_donor(self) -> dict[str, str]:
        """Return a fresh ingest mapping without reading any source file."""
        return image_id_to_donor_id(self)

    def _identity_payload(self) -> dict[str, Any]:
        value = self.to_dict()
        value.pop("content_id", None)
        return _content_identity(value)

    def to_dict(self) -> dict[str, Any]:
        value = _normalise_json(asdict(self))
        # Each image serialises itself, so a v2 entry keeps omitting measurement_context.
        value["images"] = [image.to_dict() for image in self.images]
        return value

    def verify_identity(self) -> None:
        if self.manifest_version != MANIFEST_VERSION:
            raise ContractError(
                f"unsupported QuPath cohort manifest version {self.manifest_version}"
            )
        records, expected = _validate_cohort_images(
            self.images, self.expected_donors
        )
        if records != self.images or expected != self.expected_donors:
            raise ContractError("QuPath cohort records are not canonically ordered")
        derived_panels = _cohort_panel_maps(self.images)
        for panel in self.panels:
            panel.verify_identity()
        if self.panels != derived_panels:
            raise ContractError(
                "QuPath cohort panel registry does not match its image manifests"
            )

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "QuPathCohortManifest":
        result = cls(
            manifest_version=int(value["manifest_version"]),
            cohort_name=str(value["cohort_name"]),
            expected_donors=tuple(map(str, value["expected_donors"])),
            images=tuple(
                QuPathImageManifest.from_dict(item) for item in value["images"]
            ),
            panels=tuple(
                PanelChannelMap.from_dict(item) for item in value["panels"]
            ),
            created_at_utc=str(value["created_at_utc"]),
            content_id=str(value["content_id"]),
        )
        result.verify_identity()
        return result

    @classmethod
    def read_json(cls, path: Path | str) -> "QuPathCohortManifest":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))


def image_id_to_donor_id(
    cohort_or_images: QuPathCohortManifest | Sequence[QuPathImageManifest],
) -> dict[str, str]:
    """Purely derive the exact QuPath ``Image`` value -> donor mapping.

    No path is opened and no global state is consulted. This is the intended
    ingest replacement for extracting the first digit run from an image name.
    """
    images = (
        cohort_or_images.images
        if isinstance(cohort_or_images, QuPathCohortManifest)
        else tuple(cohort_or_images)
    )
    result: dict[str, str] = {}
    for image in images:
        if image.image_id in result:
            raise ContractError(
                f"duplicate image_id cannot map to one donor: {image.image_id!r}"
            )
        result[image.image_id] = image.donor_id
    return {image_id: result[image_id] for image_id in sorted(result)}


def _atomic_write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


__all__ = [
    "ChannelMapping",
    "ContractError",
    "DuplicateObjectError",
    "FileFingerprint",
    "GeometryContractError",
    "PartialArtifactError",
    "PanelChannelMap",
    "QuPathCohortManifest",
    "QuPathImageManifest",
    "SchemaContractError",
    "canonical_json",
    "dataframe_content_sha256",
    "geometry_present",
    "hash_identifiers",
    "image_id_to_donor_id",
    "inspect_qupath_geojson",
    "iter_qupath_geojson_features",
    "sha256_file",
    "sha256_sampled_file",
    "sha256_json",
]
