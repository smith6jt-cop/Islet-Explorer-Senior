#!/usr/bin/env python3
"""Independently check nucleus admission using only ID/geometry columns for all25."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BASE = Path(__file__).resolve().parents[2]
RUN = BASE / "outputs/dynamic_portability_20261102"
UPSTREAM = Path("/data/io60-current/releases/io60-candidate-1")


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def ref(path):
    return {"path": str(path), "sha256": sha(path)}


def positive_finite(series):
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    return np.isfinite(values) & (values > 0)


def main():
    receipt_path = UPSTREAM / "03_protein_measurement/receipt.json"
    receipt = json.loads(receipt_path.read_text())
    inputs = {r["name"]: r for r in receipt["inputs"]}
    outputs = {r["path"]: r for r in receipt["outputs"]}
    min_cell_area = float(receipt["settings"]["config"]["cells"]["min_cell_area"])
    comparison_path = RUN / "protein_primary_expected_universe.json"
    comparison = json.loads(comparison_path.read_text())
    admission_path = UPSTREAM / "03_protein_measurement/nucleus_admission.csv"
    admission = pd.read_csv(admission_path, dtype={"donor_id": str}).set_index("donor_id")
    rows = []
    for item in comparison["rows"]:
        donor = item["donor_id"]
        source_record = inputs[f"cells:{donor}"]
        source_path = Path(source_record["path"])
        wanted = ["Object ID", "Nucleus: Area µm^2", "Nucleus/Cell area ratio", "Cell: Area µm^2"]
        names = pq.ParquetFile(source_path).schema_arrow.names
        missing = sorted(set(wanted) - set(names))
        row = {"donor_id": donor, "status": "unavailable", "unavailable_columns": missing,
               "source": {"path": str(source_path), "producer_recorded_sha256": source_record.get("recorded_sha256"),
                          "full_file_sha256_recomputed": False, "bytes": source_path.stat().st_size,
                          "recorded_bytes_match": source_record["bytes"] == source_path.stat().st_size,
                          "recorded_mtime_matches": source_record["mtime_ns"] == source_path.stat().st_mtime_ns}}
        if missing:
            rows.append(row)
            continue
        source = pq.ParquetFile(source_path).read(columns=wanted).to_pandas().rename(columns={"Object ID": "object_id"})
        assert source.object_id.notna().all() and source.object_id.is_unique
        raw_cells = len(source)
        source_area = pd.to_numeric(source["Cell: Area µm^2"], errors="coerce").to_numpy(float)
        source = source.loc[np.isfinite(source_area) & (source_area >= min_cell_area)]
        aligned = pq.ParquetFile(item["aligned_input"]["path"]).read(columns=["object_id", "qc_analysis_eligible"]).to_pandas()
        primary = pq.ParquetFile(item["primary_output"]["path"]).read(columns=["object_id"]).to_pandas()
        assert aligned.object_id.notna().all() and aligned.object_id.is_unique
        assert primary.object_id.notna().all() and primary.object_id.is_unique
        assert len(source) == len(aligned) and aligned.object_id.isin(source.object_id).all()
        ordered = source.set_index("object_id").reindex(aligned.object_id)
        area_present = positive_finite(ordered["Nucleus: Area µm^2"])
        ratio_present = positive_finite(ordered["Nucleus/Cell area ratio"])
        eligible = aligned.qc_analysis_eligible.to_numpy(bool)
        primary_mask = aligned.object_id.isin(primary.object_id).to_numpy()
        assert int(primary_mask.sum()) == len(primary)
        qc_rel = f"data/geometry_qc/donor_id={donor}/data.parquet"
        qc_path = UPSTREAM / "03_protein_measurement" / qc_rel
        qc = pq.ParquetFile(qc_path).read(columns=["object_id", "nucleus_cell_area_ratio"]).to_pandas()
        assert qc.object_id.notna().all() and qc.object_id.is_unique
        assert len(qc) == len(aligned) and aligned.object_id.isin(qc.object_id).all()
        raster_present = positive_finite(qc.set_index("object_id").reindex(aligned.object_id).nucleus_cell_area_ratio)
        qc_ref = ref(qc_path)
        qc_ref["producer_expected_sha256"] = outputs[qc_rel]["sha256"]
        qc_ref["receipt_match"] = qc_ref["sha256"] == qc_ref["producer_expected_sha256"]
        row.update(raw_export_cells=raw_cells, existing_ingest_filter_excluded=raw_cells-len(source),
                   cells=len(aligned), eligible=int(eligible.sum()), primary_cells=int(primary_mask.sum()),
                   missing_source_nucleus_area=int((~area_present).sum()),
                   admitted_without_source_nucleus_area=int((eligible & ~area_present).sum()),
                   primary_without_source_nucleus_area=int((primary_mask & ~area_present).sum()),
                   admitted_without_source_nucleus_ratio=int((eligible & ~ratio_present).sum()),
                   primary_without_source_nucleus_ratio=int((primary_mask & ~ratio_present).sum()),
                   admitted_without_raster_nucleus_ratio=int((eligible & ~raster_present).sum()),
                   primary_without_raster_nucleus_ratio=int((primary_mask & ~raster_present).sum()),
                   primary_ineligible=int((primary_mask & ~eligible).sum()), geometry_qc=qc_ref,
                   admission_summary_matches=(int(eligible.sum()) == int(admission.loc[donor, "admitted"]) and
                                              int((eligible & ~raster_present).sum()) == int(admission.loc[donor, "admitted_without_nucleus"])),
                   aligned_input=item["aligned_input"], primary_output=item["primary_output"])
        checks = [value == 0 for key, value in row.items() if "without_" in key or key == "primary_ineligible"]
        row["status"] = "pass" if all(checks) and row["admission_summary_matches"] and qc_ref["receipt_match"] else "fail"
        rows.append(row)
        print(donor, row["status"], row["admitted_without_source_nucleus_area"], flush=True)
    report = {"schema": "protein-nucleus-admission-audit/1", "acquisitions": len(rows),
              "criterion": "finite positive nucleus area and source nucleus/cell ratio; existing Stage03 gate separately requires finite positive raster nucleus/cell ratio",
              "criteria_source": ref(BASE / "core/multiplex5000/multiplex5000/stages/measurement.py"),
              "scope": "Read only native IDs, QC eligibility, source nucleus geometry and raster ratio; no intensity/images or new thresholds",
              "existing_ingest_min_cell_area_um2": min_cell_area,
              "source_hash_limitation": "Large raw measurement files were read by selected columns; original producer SHA is retained but not recomputed. Size and mtime comparison are not content-hash verification. Geometry QC SHA rechecked; aligned/primary SHA bound by independent preceding audit.",
              "measurement_receipt": ref(receipt_path), "published_admission": ref(admission_path),
              "input_binding_audit": ref(comparison_path), "rows": rows,
              "failures": [r["donor_id"] for r in rows if r["status"] == "fail"],
              "unavailable": [r["donor_id"] for r in rows if r["status"] == "unavailable"],
              "cells": sum(r.get("cells", 0) for r in rows),
              "eligible": sum(r.get("eligible", 0) for r in rows),
              "primary_cells": sum(r.get("primary_cells", 0) for r in rows)}
    (RUN / "protein_nucleus.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
