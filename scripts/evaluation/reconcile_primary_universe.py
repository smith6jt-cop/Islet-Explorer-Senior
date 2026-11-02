#!/usr/bin/env python3
"""Compare failed-run primary IDs to the exact mixed-route pancreas universe.

Read-only audit of primary tables, including panel/normalization availability.
This neither completes refinement nor promotes the incomplete stage to a release.
"""
import hashlib
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

BASE = Path(__file__).resolve().parents[2]
RUN = BASE / "outputs/dynamic_portability_20261102"


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def ref(path):
    return {"path": str(path), "sha256": sha(path)}


def main():
    stage = next((RUN / "releases/protein-evaluation-v2").glob(".06_typing.incomplete-*"))
    run_path = stage / "pancreas/run.json"
    run = json.loads(run_path.read_text())
    upstream = Path(run["aligned_root"]).parents[2]
    scales_path = upstream / "05_normalization/scales.csv"
    scales = pd.read_csv(scales_path, dtype={"donor_id": str})
    params_path = stage / "pancreas/refinement/representation_params.parquet"
    markers = pd.read_parquet(params_path).index.tolist()
    availability_receipt_path = upstream / "03_protein_measurement/receipt.json"
    availability_receipt = json.loads(availability_receipt_path.read_text())
    availability_hashes = {x["path"]: x["sha256"] for x in availability_receipt["outputs"]}
    rows = []
    for donor in run["donors"]:
        aligned_path = Path(run["aligned_root"]) / f"donor_id={donor}/data.parquet"
        primary_path = stage / f"pancreas/cells/donor_id={donor}/data.parquet"
        aligned = pq.ParquetFile(aligned_path).read(columns=["object_id", "qc_analysis_eligible", "cell_region"]).to_pandas()
        primary = pq.ParquetFile(primary_path).read(columns=["object_id"]).to_pandas()
        assert aligned.object_id.notna().all() and aligned.object_id.is_unique
        assert primary.object_id.notna().all() and primary.object_id.is_unique
        expected = set(aligned.loc[aligned.qc_analysis_eligible & aligned.cell_region.ne("pLN"), "object_id"])
        actual = set(primary.object_id)
        availability_path = upstream / f"03_protein_measurement/data/availability/donor_id={donor}/data.parquet"
        availability = pd.read_parquet(availability_path)
        availability_ref = ref(availability_path)
        availability_ref["producer_expected_sha256"] = availability_hashes[str(availability_path.relative_to(upstream / "03_protein_measurement"))]
        availability_ref["receipt_match"] = availability_ref["sha256"] == availability_ref["producer_expected_sha256"]
        donor_scales = scales.loc[scales.donor_id.eq(donor)]
        row = {"donor_id": donor, "aligned_cells": len(aligned), "expected_primary_cells": len(expected),
               "primary_cells": len(primary), "missing_expected_ids": len(expected - actual),
               "extra_primary_ids": len(actual - expected), "identity_check": "pass" if expected == actual else "fail",
               "panel_unavailable_markers": availability.loc[~availability.available, "marker"].tolist(),
               "refinement_markers_unavailable_in_panel": availability.loc[~availability.available & availability.marker.isin(markers), "marker"].tolist(),
               "normalization_status_counts": donor_scales.status.value_counts().to_dict(),
               "refinement_markers_without_alignment": donor_scales.loc[donor_scales.marker.isin(markers) & donor_scales.status.ne("aligned"), ["marker", "status", "n_fit_cells"]].to_dict("records"),
               "aligned_input": {**ref(aligned_path), "producer_expected_sha256": run["aligned_files"][donor]},
               "primary_output": ref(primary_path), "measurement_availability": availability_ref}
        row["aligned_input"]["receipt_match"] = row["aligned_input"]["sha256"] == row["aligned_input"]["producer_expected_sha256"]
        rows.append(row)
        print(donor, row["identity_check"], row["refinement_markers_unavailable_in_panel"], flush=True)
    report = {"schema": "protein-primary-universe-audit/1", "acquisitions": len(rows),
              "expected_universe": "aligned.qc_analysis_eligible AND aligned.cell_region != 'pLN'",
              "scope": "primary-only IDs and availability; no completed refinement or mixed-stage output",
              "primary_cells": sum(r["primary_cells"] for r in rows),
              "identity_failures": [r["donor_id"] for r in rows if r["identity_check"] != "pass"],
              "hash_failures": [r["donor_id"] for r in rows if not r["aligned_input"]["receipt_match"] or not r["measurement_availability"]["receipt_match"]],
              "refinement_markers": markers,
              "source_contract": ref(BASE / "core/multiplex5000/multiplex5000/tissues/pancreas_mixed.py"),
              "run_manifest": ref(run_path), "scales": ref(scales_path),
              "measurement_receipt": ref(availability_receipt_path), "representation_params": ref(params_path),
              "rows": rows}
    destination = RUN / "protein_primary_expected_universe.json"
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
