"""Build the one authoritative per-cell expression matrix.

Earlier workflow versions exposed several competing versions of the same
marker: raw QuPath whole-cell means, multiple corrected compartment means,
normalised values, and binary gates.  This module makes the current selection
explicit and registry-driven.  Each marker appears exactly once:

* ``uncompensated`` markers (``measurement_source`` ``qupath_raw``) come from
  the raw QuPath mean in their declared compartment, read from ``cells`` under
  the name ingest wrote (:func:`raw_measurement_column`);
* ``redsea_corrected`` and ``redsea_passthrough`` markers (source ``redsea``)
  come from the marker's declared REDSEA compartment;
* ``qupath_norm`` would select QuPath's signed ``Mean (norm)``, and the
  registry refuses it; and
* markers not acquired in an image are present as null, never silently treated
  as zero.

Since 2026-09-29 every marker is ``uncompensated``, so the stage reads no REDSEA
partition at all (``redsea_compartments=None``); see ``DATA_README.md``
("Spillover policy").

The resulting wide table is intentionally economical.  One row is one cell and
one column is one selected marker intensity.  Calibration evidence is a
separate wide table so the pipeline does not materialise tens of millions of
duplicated long-form identity rows merely to pivot them back again.
"""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .artifacts import ContractError, DuplicateObjectError
from .cells_parquet import sanitize
from .marker_registry import (
    RAW_QUPATH_COMPARTMENTS,
    MarkerRegistry,
    MarkerSpec,
    load_registry,
)

# v3: every marker reads its raw QuPath mean (registry v4), the raw lookup uses the
# column name ingest wrote, and REDSEA is read only when some marker is routed to it.
METHOD_VERSION = "selected-expression-v3-raw-qupath-means"
KEY_COLUMNS = ("donor_id", "object_id")
CONTEXT_COLUMNS = (
    "image",
    "X_centroid",
    "Y_centroid",
    "cell_region",
    "islet_num",
)


@dataclass(frozen=True)
class ExpressionAvailability:
    """The declared source and availability of one donor-marker measurement."""

    donor_id: str
    marker: str
    available: bool
    source: str
    compartment: str
    measurement_source: str
    spillover_policy: str
    calibration_status: str
    reason: str


def _indexed(
    frame: pd.DataFrame,
    *,
    name: str,
    expected_ids: pd.Index | None = None,
) -> pd.DataFrame:
    if "object_id" not in frame:
        raise ContractError(f"{name} is missing object_id")
    if frame["object_id"].isna().any():
        raise ContractError(f"{name}.object_id contains null")
    result = frame.copy()
    result["object_id"] = result["object_id"].astype(str)
    if result["object_id"].duplicated().any():
        examples = result.loc[
            result["object_id"].duplicated(keep=False), "object_id"
        ].unique()[:5]
        raise DuplicateObjectError(
            f"{name} contains duplicate object_id values: {examples.tolist()}"
        )
    result = result.set_index("object_id", drop=False)
    if expected_ids is not None:
        missing = expected_ids.difference(result.index, sort=False)
        extra = result.index.difference(expected_ids, sort=False)
        if len(missing) or len(extra):
            raise ContractError(
                f"{name} object universe differs: "
                f"missing={len(missing):,}, extra={len(extra):,}"
            )
        result = result.reindex(expected_ids)
    return result


def _qc_columns(qc: pd.DataFrame) -> tuple[str, str]:
    analysis = next(
        (
            column
            for column in ("analysis_eligible", "qc_analysis_eligible", "qc_pass")
            if column in qc
        ),
        None,
    )
    estimation = next(
        (
            column
            for column in ("estimation_eligible", "qc_estimation_eligible")
            if column in qc
        ),
        None,
    )
    if analysis is None or estimation is None:
        raise ContractError(
            "geometry QC must provide analysis_eligible and estimation_eligible"
        )
    return analysis, estimation


def raw_measurement_column(spec: MarkerSpec) -> str:
    """The ``cells`` column holding this marker's raw QuPath mean in its compartment.

    Ingest (``cells_parquet.build_sql``) writes ``Cell: <m>: Mean`` under the bare
    marker name and ``Nucleus: <m>: Mean`` as ``Nucleus__<m>``, where ``<m>`` is
    :func:`cells_parquet.sanitize` of the QuPath name, and writes no other raw
    compartment. The registry spells a marker the way the qptiff channel was
    sanitised, which keeps punctuation that ingest drops, so the registry name is
    passed through the same ``sanitize``: ``Histone_H3_(p_Ser28)`` is stored as
    ``Nucleus__Histone_H3_p_Ser28``. For every other registered marker the two
    spellings already agree.

    The lookup is exact. It used to take the first present of ``<Compartment>__<m>``,
    ``Raw__<Compartment>__<m>`` and (for Cell) the bare name, without sanitising;
    ``Cell__<m>`` and ``Raw__...`` are names nothing writes, and a similarly named
    column is not a substitute for the declared one.
    """

    if spec.compartment not in RAW_QUPATH_COMPARTMENTS:
        raise ContractError(
            f"{spec.name}: ingest carries no raw QuPath {spec.compartment!r} mean; "
            f"raw sources exist only for {list(RAW_QUPATH_COMPARTMENTS)}"
        )
    name = sanitize(spec.name)
    return name if spec.compartment == "Cell" else f"{spec.compartment}__{name}"


def _resolve_measurement(
    *,
    spec: MarkerSpec,
    acquired: bool,
    cells: pd.DataFrame,
    redsea: pd.DataFrame | None,
) -> tuple[np.ndarray, str, str]:
    """Return values, source column and availability reason.

    ``redsea`` is ``None`` when the registry routes no marker to REDSEA, in which
    case the stage did not read the partition and a ``redsea`` source is an error.
    """

    n = len(cells)
    if not acquired:
        return (
            np.full(n, np.nan, dtype=float),
            "",
            "marker_not_acquired_in_panel",
        )

    if spec.measurement_source == "qupath_norm":
        source = spec.measurement_column
        if source not in cells:
            raise ContractError(
                f"acquired marker {spec.name!r} requires exact QuPath normalized "
                f"column {source!r}; raw and RESTORE-normalized fallbacks are forbidden"
            )
        values = pd.to_numeric(cells[source], errors="coerce").to_numpy(dtype=float)
        return values, f"qupath_norm:{source}", ""

    if spec.measurement_source == "redsea":
        source = spec.measurement_column
        if redsea is None:
            raise ContractError(
                f"acquired marker {spec.name!r} requires REDSEA column {source!r}, "
                "but no REDSEA partition was supplied"
            )
        if source not in redsea:
            raise ContractError(
                f"acquired marker {spec.name!r} requires REDSEA column {source!r}"
            )
        values = pd.to_numeric(redsea[source], errors="coerce").to_numpy(dtype=float)
        return values, f"redsea:{source}", ""

    if spec.measurement_source == "qupath_raw":
        source = raw_measurement_column(spec)
        if source not in cells:
            raise ContractError(
                f"acquired marker {spec.name!r} requires the raw QuPath "
                f"{spec.compartment} mean, cells column {source!r}"
            )
        values = pd.to_numeric(cells[source], errors="coerce").to_numpy(dtype=float)
        return values, f"qupath:{source}", ""

    raise ContractError(
        f"{spec.name}: unsupported measurement source {spec.measurement_source!r}"
    )


def build_selected_expression(
    cells: pd.DataFrame,
    redsea_compartments: pd.DataFrame | None,
    geometry_qc: pd.DataFrame,
    *,
    donor_id: str,
    acquired_markers: Iterable[str],
    registry: MarkerRegistry | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assemble one donor's selected expression and availability tables.

    ``acquired_markers`` must come from the versioned QuPath image contract,
    rather than being inferred from whichever columns happen to survive a
    partial stage.  If a declared acquired measurement is missing, the build is
    fatal.  A legitimately unacquired marker is represented by an all-null
    column and an explicit availability reason.

    ``redsea_compartments`` may be ``None`` when no registered marker reads
    REDSEA (``registry.redsea_markers`` is empty), which is the case for the
    shipped registry. A supplied frame must still cover exactly the cells'
    object universe.
    """

    donor = str(donor_id).strip()
    if not donor:
        raise ContractError("donor_id is required")
    active_registry = load_registry() if registry is None else registry
    acquired = {str(marker) for marker in acquired_markers}
    unknown = sorted(acquired.difference(active_registry.by_name))
    if unknown:
        raise ContractError(
            f"QuPath panel contains marker(s) absent from the registry: {unknown}"
        )

    cell_indexed = _indexed(cells, name="cells")
    ids = cell_indexed.index
    redsea_indexed = (
        None
        if redsea_compartments is None
        else _indexed(
            redsea_compartments, name="redsea_compartments", expected_ids=ids
        )
    )
    qc_indexed = _indexed(geometry_qc, name="geometry_qc", expected_ids=ids)
    analysis_column, estimation_column = _qc_columns(qc_indexed)

    out = pd.DataFrame(
        {
            "donor_id": donor,
            "object_id": ids.to_numpy(dtype=str),
            "qc_analysis_eligible": qc_indexed[analysis_column]
            .astype(bool)
            .to_numpy(),
            "qc_estimation_eligible": qc_indexed[estimation_column]
            .astype(bool)
            .to_numpy(),
        }
    )
    for column in CONTEXT_COLUMNS:
        if column in cell_indexed:
            out[column] = cell_indexed[column].to_numpy()

    availability: list[ExpressionAvailability] = []
    for spec in active_registry.markers:
        is_acquired = spec.name in acquired
        values, source, reason = _resolve_measurement(
            spec=spec,
            acquired=is_acquired,
            cells=cell_indexed,
            redsea=redsea_indexed,
        )
        out[spec.name] = values
        availability.append(
            ExpressionAvailability(
                donor_id=donor,
                marker=spec.name,
                available=is_acquired,
                source=source,
                compartment=spec.compartment,
                measurement_source=spec.measurement_source,
                spillover_policy=spec.spillover_policy,
                calibration_status=spec.calibration_status,
                reason=reason,
            )
        )

    availability_frame = pd.DataFrame(
        [record.__dict__ for record in availability]
    )
    return out, availability_frame


def write_donor_partition(
    frame: pd.DataFrame,
    root: Path,
    *,
    donor_id: str,
    filename: str = "data.parquet",
) -> Path:
    """Write one immutable donor partition.

    Production stages construct new content-addressed run directories.  An
    existing file therefore indicates an orchestration error and is never
    overwritten here.
    """

    donor = str(donor_id)
    partition = Path(root) / f"donor_id={donor}"
    path = partition / filename
    if path.exists():
        raise FileExistsError(f"immutable donor artifact already exists: {path}")
    partition.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{filename}.",
        suffix=".tmp",
        dir=partition,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        frame.to_parquet(temporary, index=False)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def read_redsea_for_expression(
    registry: MarkerRegistry, root: Path, donor_id: str
) -> pd.DataFrame | None:
    """The donor's REDSEA partition if a registered marker reads REDSEA, else ``None``.

    With no REDSEA-routed marker the partition is never opened, so ``expression``
    runs on a donor REDSEA never reached (or on a cohort with no ``redsea/`` at all).
    """

    if not registry.redsea_markers:
        return None
    return read_single_partition(root, donor_id)


def read_single_partition(root: Path, donor_id: str) -> pd.DataFrame:
    """Read the only Parquet file in one donor partition."""

    partition = Path(root) / f"donor_id={donor_id}"
    files = sorted(partition.glob("*.parquet"))
    if len(files) != 1:
        raise ContractError(
            f"expected exactly one parquet for donor {donor_id} in {root}; "
            f"found {len(files)}"
        )
    return pd.read_parquet(files[0])


__all__ = [
    "CONTEXT_COLUMNS",
    "ExpressionAvailability",
    "KEY_COLUMNS",
    "METHOD_VERSION",
    "build_selected_expression",
    "raw_measurement_column",
    "read_redsea_for_expression",
    "read_single_partition",
    "write_donor_partition",
]
