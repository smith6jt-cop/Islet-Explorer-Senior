#!/usr/bin/env python3
"""Quantify primary/aligned coordinate conventions; do not transform scientific outputs."""
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

BASE = Path(__file__).resolve().parents[2]
RUN = BASE / "outputs/dynamic_portability_20261102"
UPSTREAM = Path("/data/io60-current/releases/io60-candidate-1")
EXPORT = Path("/home/smith6jt/IO60panc2nd/results/v3/export")
COLUMNS = ["object_id", "X_centroid", "Y_centroid"]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def ref(path):
    return {"path": str(path), "sha256": sha(path)}


def coords(path, ids=None, raw=False):
    columns = ["Object ID", "Centroid X µm", "Centroid Y µm"] if raw else COLUMNS
    frame = pq.ParquetFile(path).read(columns=columns).to_pandas()
    frame.columns = COLUMNS
    assert not frame.object_id.duplicated().any() and not frame.object_id.isna().any()
    frame = frame.set_index("object_id")
    if ids is not None:
        frame = frame.loc[ids]
    assert np.isfinite(frame.to_numpy()).all()
    return frame


def producer_sources():
    """Archive the receipt-bound implementations, separately from today's checkout."""
    sources = []
    for stage, name in (("03_protein_measurement", "cells_parquet.py"),
                        ("05_normalization", "stages/io60_normalization.py")):
        commit = read(UPSTREAM / stage / "receipt.json")["code"]["commit"]
        relative = f"core/multiplex5000/multiplex5000/{name}"
        payload = subprocess.check_output(["git", "-C", str(BASE), "show", f"{commit}:{relative}"])
        path = RUN / "provenance/coordinate_sources" / f"{commit}__{name.replace('/', '_')}"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        sources.append({**ref(path), "producer_commit": commit, "source_path": relative})
    return sources


def main():
    source_audit = RUN / "protein_primary/per_acquisition.json"
    audits = read(source_audit)
    inner_path = RUN / "protein_primary/pancreas/run.json"
    inner = read(inner_path)
    receipt_path = UPSTREAM / "03_protein_measurement/receipt.json"
    receipt = read(receipt_path)
    normalization_receipt = UPSTREAM / "05_normalization/receipt.json"
    adapter_summary = EXPORT / "typing_islets/adapter_summary_codex_noexp.json"
    adapter = {r["donor_id"]: r for r in read(adapter_summary)["images"]}
    provenance = [source_audit, inner_path, receipt_path, normalization_receipt, adapter_summary,
                  BASE / "core/phenocycler/data/qupath_manifest_v3.json",
                  Path("/home/smith6jt/IO60panc2nd/scripts/v3/islet_export_adapter.py"),
                  BASE / "core/multiplex5000/multiplex5000/cells_parquet.py",
                  BASE / "core/multiplex5000/multiplex5000/stages/io60_normalization.py",
                  BASE / "core/multiplex5000/scripts/typing/build_cells.py"]
    rows = []
    for audit in audits:
        donor = audit["donor_id"]
        primary_path = Path(audit["primary_table"]["path"])
        aligned_path = Path(inner["aligned_root"]) / f"donor_id={donor}/data.parquet"
        geometry_path = Path(audit["geometry_source"]["path"])
        expression_path = UPSTREAM / f"03_protein_measurement/data/expression/donor_id={donor}/data.parquet"
        context_input = next(x for x in receipt["inputs"] if x["name"] == f"context:{donor}")
        context_path = Path(context_input["path"])
        primary = coords(primary_path)
        aligned = coords(aligned_path, primary.index)
        geometry = coords(geometry_path, primary.index)
        expression = coords(expression_path, primary.index)
        context = coords(context_path, primary.index, raw=True)
        p, a = primary.to_numpy(), aligned.to_numpy()
        delta = p - a
        norm = np.linalg.norm(delta, axis=1)
        declared = adapter[donor]
        scale = declared["pixel_um_manifest"] / declared["pixel_um_qupath"]
        scaled_residual = p - a * scale
        ratio = np.divide(p, a, out=np.full_like(p, np.nan), where=a != 0)
        current_context_sha = sha(context_path)
        context_expected = context_input.get("sha256") or context_input.get("recorded_sha256")
        row = {"donor_id": donor, "acquisition_id": audit["acquisition_id"], "primary_cells": len(primary),
               "primary_minus_aligned_units": "um in each table's declared convention",
               "coordinate_disagreements": int(np.any(p != a, axis=1).sum()),
               "dx_min": float(delta[:, 0].min()), "dx_max": float(delta[:, 0].max()),
               "dy_min": float(delta[:, 1].min()), "dy_max": float(delta[:, 1].max()),
               "dx_abs_max": float(np.abs(delta[:, 0]).max()), "dy_abs_max": float(np.abs(delta[:, 1]).max()),
               "displacement_max": float(norm.max()), "displacement_median": float(np.median(norm)),
               "displacement_p95": float(np.percentile(norm, 95)),
               "primary_over_aligned_ratio_min": float(np.nanmin(ratio)),
               "primary_over_aligned_ratio_max": float(np.nanmax(ratio)),
               "declared_qupath_um_per_pixel": declared["pixel_um_qupath"],
               "declared_manifest_um_per_pixel": declared["pixel_um_manifest"],
               "declared_primary_over_aligned_scale": scale,
               "declared_scale_residual_abs_max_um": float(np.abs(scaled_residual).max()),
               "primary_equals_frozen_geometry_exactly": bool(np.array_equal(p, geometry.to_numpy())),
               "aligned_equals_stage03_expression_exactly": bool(np.array_equal(a, expression.to_numpy())),
               "aligned_equals_current_export_context_exactly": bool(np.array_equal(a, context.to_numpy())),
               "context_hash_check": "pass" if current_context_sha == context_expected else "fail",
               "sources": {"primary": audit["primary_table"],
                           "aligned": {"path": str(aligned_path), "sha256": inner["aligned_files"][donor],
                                       "hash_basis": "verified full-file hash in protein_primary/receipt_checks.csv"},
                           "geometry": audit["geometry_source"],
                           "stage03_expression": {"path": str(expression_path),
                               "sha256": next(x["sha256"] for x in receipt["outputs"]
                                              if x["path"] == f"data/expression/donor_id={donor}/data.parquet"),
                               "hash_basis": "bound Stage03 output receipt"},
                           "context": {"path": str(context_path), "sha256": current_context_sha,
                                       "expected_sha256": context_expected}},
               "geometric_acceptance": "not established; differing coordinate calibrations must be explicit at each consumer; no new tolerance applied"}
        rows.append(row)
        print(donor, len(primary), row["displacement_max"], row["declared_scale_residual_abs_max_um"], flush=True)
    report = {"schema": "dynamic-portability-primary-coordinates/1", "acquisitions": len(rows),
              "primary_cells": sum(r["primary_cells"] for r in rows),
              "coordinate_disagreements": sum(r["coordinate_disagreements"] for r in rows),
              "maximum_displacement_um": max(r["displacement_max"] for r in rows),
              "maximum_declared_scale_residual_abs_um": max(r["declared_scale_residual_abs_max_um"] for r in rows),
              "coordinate_contract": {
                  "primary": "build_cells takes frozen revised geometry, previously rescaled by islet_export_adapter from QuPath calibration to image-manifest calibration",
                  "aligned": "Stage03 ingests current context Centroid X/Y µm verbatim; Stage05 carries coordinates unchanged",
                  "frame_relationship": "same native IDs and image origin; two explicitly documented um-per-pixel conventions",
                  "interpretation": "descriptive comparison only; matching IDs do not establish interchangeable frames or geometric acceptance; no coordinates were rewritten"},
              "provenance": [ref(p) for p in provenance], "producer_sources": producer_sources(), "rows": rows}
    (RUN / "protein_primary_coordinates.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
