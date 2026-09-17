#!/usr/bin/env python3
"""Build the post-Xenium IF layer described in ISLET_LAYERS_SPEC.md §1.4.

Map native Xenium centroids INTO the IF frame with the exact ``refined``
transform named by the intensity producer, including its displacement field.
``warped_if_um.npy`` contains the older VALIS seed in that same direction;
it must not be paired with intensities measured through the refined transform.
The declared pixel size is read only from ``restain_sections.csv``; this
builder never opens a qptiff and never rescales the producer's coordinates.

The warp is in the recorded Xenium bundle's ``cells.parquet`` order. Join
that order to the geometry table by native ``cell_id``, as ``stamp_if_values``
does for IF intensities measured inside each Xenium cell's own segmentation.
These cells have no QuPath ``object_id``; its availability is recorded explicitly.

The output lives entirely under ``$ISLET_LAYERS_DIR/restain``. Its small
``availability.parquet`` records every donor, counts before/after the nucleus
filter, the CSV calibration, and the producer's hormone validation. The
``cells/<donor>.parquet`` files contain the IF-frame overlay. A missing value
always has an explicit availability state and reason beside it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))
from phenocycler.cohort import filter_eligible_donors  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402

MARKERS = ("INS", "GCG")
STATES = ("measured", "blank", "unavailable", "withheld", "failed")
# Half the final decimal place in the maintainer's four-decimal declaration.
# This allows its rounding precision, not the historical 0.25 -> 0.2505 change.
PIXEL_TOLERANCE_UM = 0.00005


def norm_donor(value: str) -> str:
    """The shared spelling prescribed by spec §0.3: '0112' -> '112'."""
    return str(value).strip().lstrip("0") or "0"


def read_sections(path: Path) -> pd.DataFrame:
    """Read the maintainer's declared calibration, without inspecting images."""
    sections = pd.read_csv(path, comment="#", dtype=str, keep_default_na=False)
    required = {"donor_id", "roi", "image_file", "xenium_serial", "pixel_size_um"}
    if not required.issubset(sections):
        raise ValueError(f"{path}: missing section columns {sorted(required - set(sections))}")
    if (sections[list(required)].apply(lambda col: col.str.strip()) == "").any().any():
        raise ValueError(f"{path}: blank required section value")
    sections["donor_id"] = sections["donor_id"].map(norm_donor)
    sections["pixel_size_um"] = pd.to_numeric(sections["pixel_size_um"], errors="raise")
    if not (np.isfinite(sections.pixel_size_um) & (sections.pixel_size_um > 0)).all():
        raise ValueError(f"{path}: pixel_size_um must be finite and positive")
    if sections.duplicated(["donor_id", "roi"]).any():
        raise ValueError(f"{path}: duplicate donor/roi calibration")
    return sections


def read_source(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    """Read file columns, never infer conflicting Hive partition columns."""
    return ds.dataset(path, format="parquet", partitioning=None).to_table(
        columns=columns
    ).to_pandas()


def source_inventory(root: Path) -> dict[str, list[Path]]:
    """Discover pancreas partitions from their recorded columns, not names."""
    inventory: dict[str, list[Path]] = {}
    paths = sorted(root.rglob("*.parquet"))
    if not paths:
        raise FileNotFoundError(f"No Xenium cell parquet files under {root}")
    for path in paths:
        identities = read_source(path, ["donor_id", "roi"]).drop_duplicates()
        identities["donor_id"] = identities.donor_id.map(norm_donor)
        identities = identities.drop_duplicates()
        if len(identities) != 1:
            raise ValueError(f"{path}: expected one recorded donor/roi per source partition")
        row = identities.iloc[0]
        if row.roi == "panc" and filter_eligible_donors([row.donor_id]):
            inventory.setdefault(row.donor_id, []).append(path)
    return inventory


def warp_inventory(root: Path) -> dict[str, Path]:
    """Use the export's explicit donor and ROI; filenames are not donor metadata."""
    found: dict[str, Path] = {}
    for path in sorted(root.glob("*/export.json")):
        record = json.loads(path.read_text())
        donor = norm_donor(record["donor_id"])
        if record["roi"] != "panc" or not filter_eligible_donors([donor]):
            continue
        if donor in found:
            raise ValueError(f"Duplicate pancreas warp exports for donor {donor}")
        found[donor] = path.parent
    return found


def validate_warp_count(warped: np.ndarray, n_xenium_cells: int, donor: str) -> None:
    """Assert the full source count BEFORE filtering nuclei or writing any file."""
    n_xenium_warped_to_if = len(warped)
    if n_xenium_warped_to_if != n_xenium_cells:
        raise AssertionError(
            f"{donor}: n_xenium_warped_to_if={n_xenium_warped_to_if} "
            f"!= n_xenium_cells={n_xenium_cells}"
        )
    if warped.shape != (n_xenium_cells, 2):
        raise ValueError(f"{donor}: expected an (n_xenium_cells, 2) warp, got {warped.shape}")
    if not np.isfinite(warped).all():
        raise ValueError(f"{donor}: nonfinite coordinates in the IF warp")


def hormone_records(path: Path) -> dict[str, dict]:
    if not path.is_file():
        return {}
    records: dict[str, dict] = {}
    for row in json.loads(path.read_text()):
        donor = norm_donor(row["donor_id"])
        if row["roi"] != "panc":
            continue
        if donor in records:
            raise ValueError(f"{path}: duplicate pancreas hormone validation for {donor}")
        records[donor] = row["hormone"]
    return records


def protein_records(path: Path) -> dict[str, dict]:
    """Read the measurement provenance; a seed is never a missing-data fallback."""
    if not path.is_file():
        return {}
    records = {}
    for record in json.loads(path.read_text()):
        if record["roi"] != "panc":
            continue
        donor = norm_donor(record["donor_id"])
        if donor in records:
            raise ValueError(f"{path}: duplicate pancreas protein measurement for {donor}")
        records[donor] = record
    return records


def calibration_record(export: dict, declared: float | None, protein: dict | None = None) -> dict:
    """Audit recorded producing scales without opening images or rescaling points.

    Historical exports lack the level-0 scale. In that case use the intensity
    producer's recorded sampling scale / pyramid factor and name that source.
    ``if_um_per_px`` describes the exported DAPI level, not the level-0 image.
    """
    def positive(value, name):
        value = float(value)
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"invalid producing calibration {name}={value}")
        return value

    result = dict(export_if_um_per_px=None, export_if_pixel_um=None,
                  export_pixel_size_source=export.get("pixel_size_source") or "not_recorded_in_export",
                  producing_pixel_size_um=None, producing_pixel_size_source=None,
                  pixel_size_delta_um=None, pixel_size_tolerance_um=PIXEL_TOLERANCE_UM,
                  calibration_state="unavailable", calibration_reason="no_producing_pixel_size")
    if export.get("if_um_per_px") is not None:
        result["export_if_um_per_px"] = positive(export["if_um_per_px"], "if_um_per_px")
    candidates = {}
    if export.get("if_pixel_um") is not None:
        result["export_if_pixel_um"] = positive(export["if_pixel_um"], "if_pixel_um")
        candidates["export.json:if_pixel_um"] = result["export_if_pixel_um"]
    grid = export.get("if", {})
    if grid.get("factor") is not None and result["export_if_um_per_px"] is not None:
        candidates["export.json:if_um_per_px/if.factor"] = (
            result["export_if_um_per_px"] / positive(grid["factor"], "if.factor"))
    if protein is not None:
        candidates["protein_rna_panc.json:protein_um_per_px/pyramid_factor"] = (
            positive(protein["protein_um_per_px"], "protein_um_per_px")
            / positive(protein["pyramid_factor"], "pyramid_factor"))
    if not candidates:
        return result
    source, producing = next(iter(candidates.items()))
    result.update(producing_pixel_size_um=producing, producing_pixel_size_source=source)
    if any(not np.isclose(value, producing, rtol=0, atol=1e-12) for value in candidates.values()):
        result.update(calibration_state="failed",
                      calibration_reason=f"inconsistent_producing_pixel_sizes:{candidates}")
        return result
    if declared is None:
        result["calibration_reason"] = "no_restain_section_calibration"
        return result
    delta = abs(float(declared) - producing)
    result["pixel_size_delta_um"] = delta
    if delta > PIXEL_TOLERANCE_UM + 1e-12:
        result.update(calibration_state="failed", calibration_reason=(
            f"pixel_size_mismatch: declared={declared:g} um/px; producing={producing:g} um/px; "
            f"difference={delta:g} > tolerance={PIXEL_TOLERANCE_UM:g} um/px; coordinates_not_rescaled"))
    else:
        result.update(calibration_state="measured", calibration_reason="")
    return result


def registration_gate(warp_dir: Path, donor: str) -> dict:
    """A donor-local measured failure takes precedence over a pending pairing."""
    path = warp_dir / "alignment_refined.json"
    if not path.is_file():
        return {}
    record = json.loads(path.read_text())
    export = json.loads((warp_dir / "export.json").read_text())
    if (norm_donor(record["donor_id"]), record["roi"], str(record["serial"])) != (
            donor, "panc", str(export["serial"])):
        raise ValueError(f"{path}: registration gate donor/roi/serial disagrees with export")
    alignment = record["alignment"]
    if alignment["transform_source"] != "refined":
        raise ValueError(f"{path}: registration gate does not describe refined transform")
    verdict = alignment["verdict"]
    if verdict not in {"PASS", "FAIL", "INCONCLUSIVE"}:
        raise ValueError(f"{path}: unknown registration verdict {verdict!r}")
    offset = alignment.get("median_offset_um")
    if offset is not None and (not np.isfinite(offset) or offset < 0):
        raise ValueError(f"{path}: invalid median registration offset")
    state = {"PASS": "measured", "FAIL": "failed", "INCONCLUSIVE": "unavailable"}[verdict]
    reason = ""
    if verdict != "PASS":
        measured = f"median registration offset {offset:.3f} um; " if offset is not None else ""
        reason = f"registration_gate_{verdict.lower()}: {measured}{alignment['verdict_why']}"
    return dict(registration_gate_source=str(path), registration_verdict=verdict,
                registration_state=state, registration_reason=reason,
                registration_median_offset_um=offset)


def hormone_availability(record: dict, marker: str) -> dict:
    """Preserve the recorded release decision; never fit or invent positivity."""
    evidence = record.get(marker, {})
    auroc = evidence.get("auroc")
    bar = record.get("preregistered_auroc_bar")
    if auroc is None or bar is None or not np.isfinite([auroc, bar]).all():
        return {"state": "unavailable", "reason": "no_hormone_validation",
                "auroc": None, "bar": None}
    if auroc < bar or evidence.get("clears_preregistered_bar") is not True:
        return {"state": "withheld", "reason": "auroc_below_preregistered_bar",
                "auroc": float(auroc), "bar": float(bar)}
    # A validation result does not itself supply any per-cell released calls.
    return {"state": "unavailable", "reason": "no_released_hormone_calls",
            "auroc": float(auroc), "bar": float(bar)}


def _schema() -> pa.Schema:
    fields = [(name, pa.string()) for name in
              ("donor_id", "state", "reason", "source", "frame", "image_file",
               "object_id_state", "object_id_reason", "transform_source",
               "transform_file", "transform_sha256", "transform_displacement_sha256",
               "transform_direction", "intensity_source", "producing_pixel_size_source",
               "export_pixel_size_source", "calibration_state", "calibration_reason",
               "registration_gate_source", "registration_verdict", "registration_state",
               "registration_reason", "nucleus_filter", "nucleus_filter_source")]
    fields += [(name, pa.float64()) for name in
               ("pixel_size_um", "producing_pixel_size_um", "export_if_um_per_px",
                "export_if_pixel_um", "pixel_size_delta_um", "pixel_size_tolerance_um",
                "registration_median_offset_um")]
    fields += [(name, pa.int64()) for name in
               ("n", "n_xenium_cells", "n_xenium_warped_to_if", "n_without_nucleus", "n_nucleated_cells")]
    for marker in MARKERS:
        fields += [(f"n_measured_if_{marker}{suffix}", pa.int64()) for suffix in
                   ("", "_before_nucleus_filter", "_excluded_without_nucleus")]
        prefix = f"hormone_positivity_{marker}_"
        fields += [(prefix + "state", pa.string()), (prefix + "reason", pa.string()),
                   (prefix + "auroc", pa.float64()), (prefix + "bar", pa.float64())]
    return pa.schema(fields)


def _status(donor: str, source: Path) -> dict:
    row = dict.fromkeys(_schema().names)
    row.update(donor_id=donor, state="unavailable", reason="no_restain_acquisition_record",
               source=str(source), frame=f"if_um:{donor}:panc",
               object_id_state="unavailable", object_id_reason="no_qupath_object_id",
               calibration_state="unavailable", calibration_reason="no_producing_pixel_size")
    for marker in MARKERS:
        row[f"hormone_positivity_{marker}_state"] = "unavailable"
        row[f"hormone_positivity_{marker}_reason"] = "no_if_intensity_measurement"
    return row


def inspect_donor(paths: list[Path], warp_dir: Path, donor: str, section: pd.Series,
                  protein: dict) -> tuple:
    """Validate source identity and full warp cardinality without aligning records."""
    frames = []
    for path in paths:
        names = set(pq.ParquetFile(path).schema_arrow.names)
        required = {"donor_id", "roi", "cell_id", *(f"feat__if_{m}" for m in MARKERS)}
        if not required.issubset(names):
            raise ValueError(f"{path}: missing IF source columns {sorted(required - names)}")
        frames.append(read_source(path, sorted(required)))
    cells = pd.concat(frames, ignore_index=True)
    cells["donor_id"] = cells.donor_id.map(norm_donor)
    if not cells.donor_id.eq(donor).all() or not cells.roi.eq("panc").all():
        raise ValueError(f"{donor}: source donor/roi disagrees with requested IF section")
    if cells.cell_id.isna().any() or cells.cell_id.duplicated().any() or cells.cell_id.eq("").any():
        raise ValueError(f"{donor}: missing or duplicate Xenium cell_id")
    export = json.loads((warp_dir / "export.json").read_text())
    if (norm_donor(export["donor_id"]), export["roi"], str(export["serial"])) != (
        donor, "panc", section.xenium_serial
    ):
        raise ValueError(f"{donor}: warp export donor/roi/serial disagrees with sections CSV")
    if (norm_donor(protein["donor_id"]), protein["roi"], str(protein["xenium_serial"]),
            Path(protein["xenium_bundle"]).resolve(), protein["n_cells"]) != (
            donor, "panc", section.xenium_serial,
            Path(export["xenium_bundle"]).resolve(), len(cells)):
        raise ValueError(f"{donor}: intensity producer disagrees with warp donor/bundle/count")
    if (Path(protein["qptiff"]).resolve() != Path(export["qptiff"]).resolve()
            or Path(export["qptiff"]).name != Path(section.image_file).name):
        raise ValueError(f"{donor}: intensity producer disagrees with IF image identity")
    if protein.get("transform_source") != "refined":
        raise ValueError(f"{donor}: unsupported intensity transform {protein.get('transform_source')!r}")
    path = warp_dir / "refined.json"
    payload = json.loads(path.read_text())
    displacement = warp_dir / "refined_displacement.npz"
    if payload.get("has_displacement") and not displacement.is_file():
        raise FileNotFoundError(2, "refined transform declares a missing displacement", str(displacement))
    native = np.load(warp_dir / "xenium_centroids_um.npy", allow_pickle=False)
    validate_warp_count(native, len(cells), donor)
    warped = Transform.load(warp_dir, stem="refined").apply(native)
    validate_warp_count(warped, len(cells), donor)
    provenance = dict(transform_source="refined", transform_file=str(path),
                      transform_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      transform_displacement_sha256=(hashlib.sha256(displacement.read_bytes()).hexdigest()
                                                     if payload.get("has_displacement") else None),
                      transform_direction="xenium_um_to_if_um")
    return cells, warped, provenance


def prepare_restain_cells(cells: pd.DataFrame, warped: np.ndarray,
                          warp_dir: Path, status: dict) -> tuple:
    """Join bundle-ordered IF coordinates to geometry by Xenium ``cell_id``.

    Neither the geometry table's row order nor a second segmentation is a join
    key. Refuse count/key mismatches before filtering out cells without nuclei.
    IF intensities are already stamped by ``postxen.stamp_if_values`` and stay
    with their recorded cell IDs, including NaN where no IF value was measured.
    """
    donor = status["donor_id"]
    validate_warp_count(warped, len(cells), donor)
    export = json.loads((warp_dir / "export.json").read_text())
    bundle = export.get("xenium_bundle")
    if not bundle:
        raise ValueError(f"{donor}: warp export has no recorded xenium_bundle")
    order = read_source(Path(bundle) / "cells.parquet", ["cell_id"]).cell_id
    if len(order) != len(warped):
        raise ValueError(f"{donor}: {len(warped)} warp rows for {len(order)} bundle cells; "
                         "refusing to pad or truncate")
    if order.isna().any() or order.duplicated().any() or order.eq("").any():
        raise ValueError(f"{donor}: missing or duplicate bundle cell_id")
    order = pd.Index(order.astype(str))
    # export_dapi writes this object-string array directly from the bundle.
    recorded = np.load(warp_dir / "xenium_cell_ids.npy", allow_pickle=True)
    if recorded.ndim != 1 or not np.array_equal(recorded, order.to_numpy()):
        raise ValueError(f"{donor}: warp cell_id order disagrees with bundle cells.parquet")
    keys = pd.Index(cells.cell_id.astype(str))
    positions = order.get_indexer(keys)
    if keys.has_duplicates or np.any(positions < 0):
        raise ValueError(f"{donor}: geometry cell_id keys disagree with bundle cells.parquet")

    nucleus_ids = np.load(warp_dir / "xenium_nuclei_cell_ids.npy", allow_pickle=False)
    if nucleus_ids.ndim != 1:
        raise ValueError(f"{donor}: expected a vector of Xenium nucleus cell_id values")
    nuclei = pd.Index(nucleus_ids)
    if (nuclei.hasnans or (nuclei == "").any()
            or not nuclei.isin(order).all()):
        raise ValueError(f"{donor}: invalid Xenium nucleus cell_id keys")
    # Multiple nuclei may belong to one cell; membership keeps that cell once.
    keep = keys.isin(nuclei)
    details = dict(n_without_nucleus=int((~keep).sum()), n_nucleated_cells=int(keep.sum()),
                   nucleus_filter="cell_id_in_xenium_nuclei_cell_ids",
                   nucleus_filter_source=str(warp_dir / "xenium_nuclei_cell_ids.npy"))
    for marker in MARKERS:
        measured = np.isfinite(cells[f"feat__if_{marker}"].to_numpy(dtype=float))
        details[f"n_measured_if_{marker}_before_nucleus_filter"] = int(measured.sum())
        details[f"n_measured_if_{marker}_excluded_without_nucleus"] = int((measured & ~keep).sum())
        details[f"n_measured_if_{marker}"] = int((measured & keep).sum())
    if not keep.any():
        return None, {**details, "state": "unavailable", "reason": "no_recorded_xenium_nuclei"}

    overlay = cells.loc[keep, ["donor_id", "cell_id", *(f"feat__if_{m}" for m in MARKERS)]].copy()
    overlay["x_um"] = warped[positions[keep], 0]
    overlay["y_um"] = warped[positions[keep], 1]
    overlay["frame"] = status["frame"]
    for field in ("transform_source", "transform_file", "transform_sha256",
                  "transform_displacement_sha256", "transform_direction", "intensity_source",
                  "pixel_size_um", "producing_pixel_size_um", "producing_pixel_size_source",
                  "pixel_size_tolerance_um"):
        overlay[field] = status[field]
    overlay["object_id"] = pd.array([None] * len(overlay), dtype="string")
    overlay["object_id_state"] = "unavailable"
    overlay["object_id_reason"] = "no_qupath_object_id"
    for marker in MARKERS:
        field = f"feat__if_{marker}"
        measured = np.isfinite(overlay[field].to_numpy(dtype=float))
        overlay[f"{field}_state"] = np.where(measured, "measured", "unavailable")
        overlay[f"{field}_reason"] = np.where(measured, "", "no_finite_if_measurement")
        details[f"n_measured_if_{marker}"] = int(measured.sum())
        # No per-cell hormone calls have been released by the source lane.
        prefix = f"hormone_positivity_{marker}"
        overlay[prefix] = pd.array([None] * len(overlay), dtype="boolean")
        for field in ("state", "reason"):
            overlay[f"{prefix}_{field}"] = status[f"{prefix}_{field}"]
    if not any(details[f"n_measured_if_{m}"] for m in MARKERS):
        return None, {**details, "state": "unavailable", "reason": "no_finite_if_measurements"}
    details.update(state="measured", reason="", n=len(overlay))
    return pa.Table.from_pandas(overlay, preserve_index=False), details


def build_restain_layer(*, cells_xen_root: Path, restain_root: Path,
                        sections_csv: Path, xenium_paths_csv: Path,
                        layers_dir: Path) -> pd.DataFrame:
    """Validate every paired donor before creating any output, then write one root."""
    sections = read_sections(sections_csv)
    sections = sections[sections.roi == "panc"].set_index("donor_id")
    inventory_csv = sections_csv.with_name("restain_registration_inputs.csv")
    acquisitions = pd.DataFrame()
    if inventory_csv.is_file():
        acquisitions = pd.read_csv(inventory_csv, dtype=str, comment="#", keep_default_na=False)
        acquisitions["donor_id"] = acquisitions.donor_id.map(norm_donor)
        acquisitions = acquisitions[acquisitions.roi == "panc"].set_index("donor_id")
        if acquisitions.index.has_duplicates:
            raise ValueError(f"{inventory_csv}: duplicate pancreas acquisition")
    sources = source_inventory(cells_xen_root)
    warps = warp_inventory(restain_root / "valis")
    validation = hormone_records(restain_root / "hormone_panc_refined.json")
    protein_path = restain_root / "protein_rna_panc.json"
    proteins = protein_records(protein_path)
    known = pd.read_csv(xenium_paths_csv, dtype=str, comment="#", keep_default_na=False)
    donors = sorted(set(filter_eligible_donors(
        set(known.donor_id.map(norm_donor)) | set(sections.index) | set(sources)
        | set(warps) | set(acquisitions.index)
    )))
    statuses, overlays = [], {}
    for donor in donors:
        row = _status(donor, inventory_csv if inventory_csv.is_file() else sections_csv)
        if donor in acquisitions.index or donor in warps or donor in sections.index:
            row.update(reason="pending_qupath_classifier_and_geojson_export")
            if donor in acquisitions.index:
                row["image_file"] = acquisitions.loc[donor, "image_file"]
        if donor in sections.index:
            section = sections.loc[donor]
            row.update(pixel_size_um=float(section.pixel_size_um), image_file=section.image_file)
        if donor in warps:
            export = json.loads((warps[donor] / "export.json").read_text())
            row.update(source=str(warps[donor] / "export.json"))
            row.update(calibration_record(export, row["pixel_size_um"], proteins.get(donor)))
            gate = registration_gate(warps[donor], donor)
            row.update(gate)
            if gate and gate["registration_state"] != "measured":
                row.update(state=gate["registration_state"], reason=gate["registration_reason"],
                           source=gate["registration_gate_source"])
                statuses.append(row)
                continue
        paths = sources.get(donor, [])
        paired = bool(paths) and all(
            all(f"feat__if_{m}" in pq.ParquetFile(path).schema_arrow.names for m in MARKERS)
            for path in paths
        )
        if not paired:
            statuses.append(row)
            continue
        if donor not in sections.index:
            row.update(reason="no_restain_section_calibration")
            statuses.append(row)
            continue
        if donor not in warps:
            row.update(reason="no_xenium_to_if_warp")
            statuses.append(row)
            continue
        if donor not in proteins:
            row.update(reason="no_if_intensity_provenance")
            statuses.append(row)
            continue
        section = sections.loc[donor]
        try:
            cells, warped, provenance = inspect_donor(paths, warps[donor], donor, section,
                                                       proteins[donor])
            row.update(provenance, intensity_source=str(protein_path))
            row.update(n_xenium_cells=len(cells), n_xenium_warped_to_if=len(warped),
                       source=";".join(str(p) for p in paths))
            for marker in MARKERS:
                evidence = hormone_availability(validation.get(donor, {}), marker)
                for field, value in evidence.items():
                    row[f"hormone_positivity_{marker}_{field}"] = value
            # All donor shape/identity assertions precede every artifact write.
            overlay, details = prepare_restain_cells(cells, warped, warps[donor], row)
        except FileNotFoundError as exc:
            row.update(reason=f"missing_restain_input:{Path(exc.filename or str(exc)).name}")
            statuses.append(row)
            continue
        row.update(details)
        if row["calibration_state"] != "measured":
            row.update(state=row["calibration_state"], reason=row["calibration_reason"], n=None)
            overlay = None
        if overlay is not None:
            overlays[donor] = overlay
        statuses.append(row)
    status_table = pa.Table.from_pylist(statuses, schema=_schema())
    destination = Path(layers_dir) / "restain"
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "cells").mkdir(exist_ok=True)
    for donor, table in overlays.items():
        pq.write_table(table, destination / "cells" / f"{donor}.parquet", compression="zstd")
    # A rejected rebuild must retire the former overlay, including direct readers.
    for row in statuses:
        if row["state"] != "measured":
            (destination / "cells" / f"{row['donor_id']}.parquet").unlink(missing_ok=True)
    # Publish availability last; every rejection retains its measured evidence.
    pq.write_table(status_table, destination / "availability.parquet", compression="zstd")
    return status_table.to_pandas()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    core_data = REPO / "core/phenocycler/data"
    parser.add_argument("--restain-root", type=Path,
                        default=REPO / "data/development/restain")
    parser.add_argument("--sections-csv", type=Path,
                        default=core_data / "integration/restain_sections.csv")
    parser.add_argument("--xenium-paths-csv", type=Path,
                        default=core_data / "integration/xenium_paths.csv")
    parser.add_argument("--layers-dir", type=Path,
                        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")))
    args = parser.parse_args(argv)
    result = build_restain_layer(cells_xen_root=core_data / "integration/cells_xen", **vars(args))
    print(result[["donor_id", "state", "reason", "n", "n_xenium_cells",
                  "n_xenium_warped_to_if", "pixel_size_um", "producing_pixel_size_um"]].to_string(index=False))
    failed = result[result.state == "failed"]
    for row in failed.itertuples():
        print(f"FAILED donor {row.donor_id}: {row.reason}", file=sys.stderr)
    return 1 if len(failed) else 0


if __name__ == "__main__":
    raise SystemExit(main())
