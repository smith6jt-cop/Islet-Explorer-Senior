#!/usr/bin/env python3
"""Refresh a small mutable projection of upstream QuPath cell measurements.

The full QuPath CSV exports are the inputs to RESTORE.  They are very large, so
the fundamental diagnostics share one narrow donor-partitioned projection
instead of scanning the exports independently.  This development cache is
replaceable and carries no content ID, source hash, run receipt, or publication
semantics.

``Iba1_norm`` is read only from QuPath's exact
``Cell: Iba1: Mean (norm)`` measurement.  Raw Iba1 and RESTORE-derived values
are deliberately absent from the projection.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import tempfile
import time
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
WORKSPACE = REPO.parent
DEFAULT_SPEC = REPO / "cohort-spec.json"
DEFAULT_OUTPUT_ROOT = WORKSPACE / "data" / "development" / "qupath_fundamental_inputs"
METHOD_VERSION = "qupath-fundamental-input-projection-v3"
WORKSPACE_MARKER = ".qupath_fundamental_inputs_workspace"

# Alias -> exact column in the full QuPath cell export.  These are intentionally
# explicit: accepting every column ending in ``norm`` would recreate the source
# ambiguity this projection exists to eliminate.
SOURCE_MEASUREMENTS: dict[str, str] = {
    "image": "Image",
    "object_id": "Object ID",
    "X_centroid": "Centroid X µm",
    "Y_centroid": "Centroid Y µm",
    "cell_area": "Cell: Area µm^2",
    "cell_solidity": "Cell: Solidity",
    "DAPI": "Cell: DAPI: Mean",
    "CD3e": "Cell: CD3e: Mean",
    "Iba1_norm": "Cell: Iba1: Mean (norm)",
    "Ker8_18": "Cell: Ker8-18: Mean",
}
DIAGNOSTIC_SOURCE_KINDS: dict[str, str] = {
    "CD3e": "qupath_raw_cell_mean",
    "Iba1_norm": "qupath_acquisition_normalized_cell_mean",
    "Ker8_18": "qupath_raw_cell_mean",
}
NUMERIC_COLUMNS = tuple(
    column for column in SOURCE_MEASUREMENTS if column not in {"image", "object_id"}
)
PROJECTION_COLUMNS = ("donor_id", *SOURCE_MEASUREMENTS)
PROJECTION_SCHEMA = pa.schema(
    [
        ("donor_id", pa.string()),
        ("image", pa.string()),
        ("object_id", pa.string()),
        *((column, pa.float64()) for column in NUMERIC_COLUMNS),
    ]
)


@dataclass(frozen=True)
class QuPathSource:
    donor_id: str
    image_id: str
    measurement_csv: Path


def _path(base: Path, value: object) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def _measurement_path(
    base: Path, defaults: Mapping[str, object], image: Mapping[str, object]
) -> Path:
    value = image.get("measurement_csv", defaults.get("measurement_csv"))
    if isinstance(value, list):
        if not value:
            raise ValueError("measurement_csv list cannot be empty")
        index = image.get("measurement_csv_index")
        if not isinstance(index, int) or isinstance(index, bool):
            raise ValueError("list-valued measurement_csv requires measurement_csv_index")
        if not 0 <= index < len(value):
            raise ValueError(f"measurement_csv_index {index} is outside 0..{len(value) - 1}")
        value = value[index]
    elif image.get("measurement_csv_index") is not None:
        raise ValueError("measurement_csv_index is valid only for list-valued measurement_csv")
    if not value:
        raise ValueError("image has no measurement_csv")
    return _path(base, value)


def load_qupath_sources(
    spec_path: str | Path,
    donors: Sequence[str],
) -> tuple[QuPathSource, ...]:
    """Resolve selected donors from the editable cohort specification."""

    path = Path(spec_path).expanduser().resolve()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        defaults = payload.get("defaults", {})
        images = payload["images"]
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read QuPath cohort specification: {path}") from exc
    if not isinstance(defaults, Mapping) or not isinstance(images, list):
        raise ValueError("cohort specification defaults/images have invalid types")

    requested = tuple(map(str, donors))
    if not requested or len(set(requested)) != len(requested):
        raise ValueError("donors must be nonempty and unique")
    records: dict[str, QuPathSource] = {}
    for index, image in enumerate(images):
        if not isinstance(image, Mapping):
            raise ValueError(f"cohort image {index} must be an object")
        donor = str(image.get("donor_id", ""))
        if donor not in requested:
            continue
        if donor in records:
            raise ValueError(f"cohort specification contains duplicate donor {donor}")
        image_id = str(image.get("image_id", ""))
        if not image_id:
            raise ValueError(f"cohort donor {donor} has no image_id")
        duplicate_image = next(
            (existing.donor_id for existing in records.values() if existing.image_id == image_id),
            None,
        )
        if duplicate_image is not None:
            raise ValueError(
                f"cohort donors {duplicate_image} and {donor} share image_id {image_id!r}"
            )
        records[donor] = QuPathSource(
            donor_id=donor,
            image_id=image_id,
            measurement_csv=_measurement_path(path.parent, defaults, image),
        )
    missing = sorted(set(requested).difference(records))
    if missing:
        raise ValueError(f"donors absent from cohort specification: {missing}")
    return tuple(records[donor] for donor in requested)


def _header(path: Path) -> tuple[str, ...]:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            return tuple(next(csv.reader(handle)))
    except (OSError, StopIteration) as exc:
        raise ValueError(f"cannot read QuPath measurement header: {path}") from exc


def _validate_source_header(path: Path) -> None:
    columns = set(_header(path))
    missing = sorted(set(SOURCE_MEASUREMENTS.values()).difference(columns))
    if missing:
        raise ValueError(f"QuPath measurement export {path} is missing exact columns {missing}")


def _as_float(batch: pa.RecordBatch, source_name: str) -> pa.Array:
    try:
        return pc.cast(batch.column(source_name), pa.float64(), safe=False)
    except (pa.ArrowInvalid, pa.ArrowNotImplementedError) as exc:
        raise ValueError(f"QuPath measurement {source_name!r} contains a nonnumeric value") from exc


def _project_batch(
    batch: pa.RecordBatch,
    *,
    source: QuPathSource,
) -> pa.Table:
    image = batch.column(SOURCE_MEASUREMENTS["image"])
    selected = batch.filter(pc.equal(image, source.image_id))
    if not selected.num_rows:
        return pa.Table.from_batches([], schema=PROJECTION_SCHEMA)
    object_ids = selected.column(SOURCE_MEASUREMENTS["object_id"])
    if object_ids.null_count:
        raise ValueError(f"{source.donor_id}: QuPath export contains a null Object ID")
    arrays: list[pa.Array | pa.ChunkedArray] = [
        pa.array([source.donor_id] * selected.num_rows, type=pa.string()),
        selected.column(SOURCE_MEASUREMENTS["image"]),
        object_ids,
    ]
    arrays.extend(_as_float(selected, SOURCE_MEASUREMENTS[column]) for column in NUMERIC_COLUMNS)
    return pa.Table.from_arrays(arrays, schema=PROJECTION_SCHEMA)


def _validate_output_root(output_root: Path) -> None:
    root = output_root.expanduser().resolve()
    dangerous = {Path("/"), Path.home().resolve(), REPO.resolve(), WORKSPACE.resolve()}
    if root in dangerous or root.parent == root:
        raise ValueError("projection output must be a dedicated development subdirectory")
    if root.is_symlink():
        raise ValueError(f"projection output cannot be a symlink: {root}")
    if root.exists() and not (root / WORKSPACE_MARKER).is_file():
        raise ValueError(f"refusing to replace a projection not owned by this tool: {root}")


def validate_qupath_projection(
    input_root: str | Path,
    donors: Sequence[str] = (),
) -> dict[str, object]:
    """Validate the mutable projection's source semantics, without hashes.

    Column names alone are not provenance: the retired RESTORE flat tables also
    used an ``Iba1_norm`` alias.  Every diagnostic therefore requires the
    projection tool's live workspace marker and exact source declaration before
    reading a donor Parquet.
    """

    root = Path(input_root).expanduser().resolve()
    marker = root / WORKSPACE_MARKER
    summary_path = root / "projection_summary.json"
    try:
        marker_version = marker.read_text(encoding="utf-8").strip()
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"{root} is not a declared upstream QuPath input projection") from exc
    if marker_version != METHOD_VERSION or summary.get("method_version") != METHOD_VERSION:
        raise ValueError(f"{root}: unsupported QuPath projection method")
    expected = {
        "status": "complete",
        "source_kind": "upstream QuPath full-cell measurement export",
        "source_measurements": SOURCE_MEASUREMENTS,
        "restore_outputs_read": False,
    }
    drift = {key: summary.get(key) for key, value in expected.items() if summary.get(key) != value}
    if drift:
        raise ValueError(f"{root}: QuPath projection source contract differs: {sorted(drift)}")
    declared = summary.get("donors")
    if not isinstance(declared, Mapping):
        raise ValueError(f"{root}: QuPath projection has no donor inventory")
    donor_images = summary.get("donor_images")
    if not isinstance(donor_images, Mapping) or set(map(str, donor_images)) != set(
        map(str, declared)
    ):
        raise ValueError(f"{root}: QuPath projection donor/image inventory is invalid")
    for donor, count in declared.items():
        if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
            raise ValueError(f"{root}: donor {donor} has invalid declared row count {count!r}")
        if not isinstance(donor_images.get(str(donor)), str) or not donor_images[str(donor)]:
            raise ValueError(f"{root}: donor {donor} has no declared image identity")
    image_ids = list(map(str, donor_images.values()))
    if len(set(image_ids)) != len(image_ids):
        raise ValueError(f"{root}: QuPath projection assigns one image to multiple donors")
    requested = tuple(map(str, donors))
    missing = sorted(set(requested).difference(map(str, declared)))
    if missing:
        raise ValueError(f"{root}: donors absent from QuPath projection: {missing}")
    return dict(summary)


def load_qupath_projection_donor(
    input_root: str | Path,
    donor: str,
    columns: Sequence[str],
) -> tuple[pd.DataFrame, tuple[Path, ...], dict[str, object]]:
    """Load one declared donor partition and enforce its mutable source contract.

    This is deliberately identity-based rather than hash-based: development projections
    may be refreshed in place, but a consumer must still prove that the rows, donor, image,
    object IDs, and numeric schema agree with the projection summary it is claiming.
    """

    root = Path(input_root).expanduser().resolve()
    donor_id = str(donor)
    summary = validate_qupath_projection(root, (donor_id,))
    requested = tuple(dict.fromkeys(map(str, columns)))
    identity_columns = ("donor_id", "image", "object_id")
    read_columns = tuple(dict.fromkeys((*identity_columns, *requested)))
    paths = tuple(sorted((root / f"donor_id={donor_id}").glob("*.parquet")))
    if not paths:
        raise FileNotFoundError(f"{donor_id}: no QuPath input Parquets below {root}")

    numeric_columns = tuple(column for column in read_columns if column not in identity_columns)
    frames: list[pd.DataFrame] = []
    for path in paths:
        schema = pq.ParquetFile(path).schema_arrow
        missing = sorted(set(read_columns).difference(schema.names))
        if missing:
            raise ValueError(f"{path} is missing required QuPath input columns {missing}")
        for column in identity_columns:
            source_type = schema.field(column).type
            if not (pa.types.is_string(source_type) or pa.types.is_large_string(source_type)):
                raise ValueError(
                    f"{path}: QuPath identity {column} must use string storage, got {source_type}"
                )
        for column in numeric_columns:
            source_type = schema.field(column).type
            if not (pa.types.is_integer(source_type) or pa.types.is_floating(source_type)):
                raise ValueError(
                    f"{path}: QuPath input {column} must use integer or floating numeric "
                    f"storage, got {source_type}"
                )
        frames.append(pd.read_parquet(path, columns=list(read_columns)))

    frame = pd.concat(frames, ignore_index=True)
    expected_count = int(summary["donors"][donor_id])
    if len(frame) != expected_count:
        raise ValueError(
            f"{donor_id}: projection row count {len(frame)} differs from declared "
            f"count {expected_count}"
        )
    internal_donors = frame["donor_id"].astype("string")
    if internal_donors.isna().any() or set(internal_donors.astype(str)) != {donor_id}:
        raise ValueError(f"{donor_id}: projection rows carry the wrong internal donor_id")
    expected_image = str(summary["donor_images"][donor_id])
    internal_images = frame["image"].astype("string")
    if internal_images.isna().any() or set(internal_images.astype(str)) != {expected_image}:
        raise ValueError(
            f"{donor_id}: projection rows do not match declared image {expected_image!r}"
        )
    identifiers = frame["object_id"].astype("string")
    if identifiers.isna().any() or identifiers.eq("").any() or identifiers.duplicated().any():
        raise ValueError(f"{donor_id}: object IDs must be non-null, nonempty, and unique")
    frame["donor_id"] = internal_donors
    frame["image"] = internal_images
    frame["object_id"] = identifiers
    for column in numeric_columns:
        values = pd.to_numeric(frame[column], errors="raise").astype("float64")
        if np.isinf(values.to_numpy(dtype=float)).any():
            raise ValueError(f"{donor_id}: {column} contains infinite values")
        frame[column] = values
    return frame, paths, summary


def refresh_qupath_projection(
    spec_path: str | Path,
    output_root: str | Path,
    donors: Sequence[str],
    *,
    block_size: int = 64 * 1024 * 1024,
) -> Path:
    """Scan each required export once and replace the shared mutable projection."""

    if block_size <= 0:
        raise ValueError("block_size must be positive")
    sources = load_qupath_sources(spec_path, donors)
    for path in sorted({source.measurement_csv for source in sources}):
        _validate_source_header(path)

    output = Path(output_root).expanduser().resolve()
    _validate_output_root(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    writers: dict[str, pq.ParquetWriter] = {}
    counts = {source.donor_id: 0 for source in sources}
    started = time.time()
    try:
        by_path: dict[Path, list[QuPathSource]] = defaultdict(list)
        for source in sources:
            by_path[source.measurement_csv].append(source)
        include_columns = list(SOURCE_MEASUREMENTS.values())
        string_types = {column: pa.string() for column in include_columns}
        for measurement_path, path_sources in sorted(
            by_path.items(), key=lambda item: str(item[0])
        ):
            reader = pacsv.open_csv(
                measurement_path,
                read_options=pacsv.ReadOptions(block_size=block_size, use_threads=True),
                parse_options=pacsv.ParseOptions(newlines_in_values=True),
                convert_options=pacsv.ConvertOptions(
                    include_columns=include_columns,
                    column_types=string_types,
                    strings_can_be_null=True,
                ),
            )
            for batch in reader:
                for source in path_sources:
                    projected = _project_batch(batch, source=source)
                    if not projected.num_rows:
                        continue
                    writer = writers.get(source.donor_id)
                    if writer is None:
                        partition = staging / f"donor_id={source.donor_id}"
                        partition.mkdir(parents=True)
                        writer = pq.ParquetWriter(
                            partition / "data.parquet",
                            PROJECTION_SCHEMA,
                            compression="zstd",
                        )
                        writers[source.donor_id] = writer
                    writer.write_table(projected)
                    counts[source.donor_id] += projected.num_rows
        for writer in writers.values():
            writer.close()
        writers.clear()
        empty = sorted(donor for donor, count in counts.items() if not count)
        if empty:
            raise ValueError(f"selected QuPath images produced no rows for donors {empty}")

        (staging / WORKSPACE_MARKER).write_text(METHOD_VERSION + "\n", encoding="utf-8")
        summary = {
            "status": "complete",
            "method_version": METHOD_VERSION,
            "workspace_semantics": "mutable development projection; replaced on refresh",
            "source_kind": "upstream QuPath full-cell measurement export",
            "cohort_spec": str(Path(spec_path).expanduser().resolve()),
            "source_measurements": SOURCE_MEASUREMENTS,
            "source_files": sorted({str(source.measurement_csv) for source in sources}),
            "donors": counts,
            "donor_images": {source.donor_id: source.image_id for source in sources},
            "restore_outputs_read": False,
            "elapsed_seconds": time.time() - started,
        }
        (staging / "projection_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if output.exists():
            shutil.rmtree(output)
        staging.replace(output)
    except Exception:
        for writer in writers.values():
            writer.close()
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="+", required=True)
    parser.add_argument("--cohort-spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--block-size-mb", type=int, default=64)
    args = parser.parse_args(argv)
    try:
        output = refresh_qupath_projection(
            args.cohort_spec,
            args.output_root,
            tuple(map(str, args.donors)),
            block_size=args.block_size_mb * 1024 * 1024,
        )
    except ValueError as exc:
        parser.error(str(exc))
    print(f"refreshed upstream QuPath inputs: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
