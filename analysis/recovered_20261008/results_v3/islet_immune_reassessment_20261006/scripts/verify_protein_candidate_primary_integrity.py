#!/usr/bin/env python3
"""Verify every prepared candidate primary UUID and mandatory exclusion exactly.

Cell counts and operational labels are accounting data, never accuracy targets.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

def read(path, columns=None):
    return pq.ParquetFile(path).read(columns=columns).to_pandas()

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate-root", type=Path, required=True)
    args = ap.parse_args(); root = args.candidate_root.resolve()
    manifest = json.loads((root / "source_manifest.json").read_text())
    rows = []
    for method in manifest["arms"]:
        if json.loads((root / method / "typing/run.json").read_text())["status"] not in ("primary_complete", "complete"):
            raise RuntimeError(f"Primary method is not complete: {method}")
        for donor in manifest["donors"]:
            calls_path = root / method / f"typing/cells/donor_id={donor}/data.parquet"
            geometry_path = root / f"geometry/revised_cells/donor_id={donor}/data.parquet"
            calls = read(calls_path, ["object_id", "cell_class", "qc_analysis_eligible", "X_centroid", "Y_centroid"])
            geometry = read(geometry_path, ["object_id", "qc_analysis_eligible", "nucleus_area_um2"])
            for name, frame in [("calls", calls), ("geometry", geometry)]:
                if frame.object_id.isna().any() or frame.object_id.duplicated().any():
                    raise ValueError(f"Invalid {name} UUID keys: {method} {donor}")
            expected = set(geometry.loc[geometry.qc_analysis_eligible.astype(bool), "object_id"])
            actual = set(calls.object_id)
            if actual != expected:
                raise ValueError(f"Primary UUID universe differs: {method} {donor}, missing={len(expected-actual)}, extra={len(actual-expected)}")
            if calls.cell_class.isna().any() or not calls.qc_analysis_eligible.astype(bool).all():
                raise ValueError(f"Null call or ineligible primary row: {method} {donor}")
            if not np.isfinite(calls[["X_centroid", "Y_centroid"]].to_numpy(float)).all():
                raise ValueError(f"Invalid coordinates: {method} {donor}")
            admitted = calls[["object_id"]].merge(geometry, on="object_id", validate="one_to_one")
            nucleus = admitted.nucleus_area_um2.to_numpy(float)
            if not (np.isfinite(nucleus) & (nucleus > 0)).all():
                raise ValueError(f"No-nucleus object survived: {method} {donor}")
            record = next(r for r in manifest["records"] if r["arm"] == method and r["donor_id_native"] == donor)
            assignment = read(record["assignment"]["path"], ["object_id", "analysis_eligible"])
            denied = set(assignment.loc[~assignment.analysis_eligible.astype(bool), "object_id"])
            violation = actual.intersection(denied)
            if violation:
                raise ValueError(f"Mandatory exclusion was returned as typed: {method} {donor}, n={len(violation)}")
            rows.append(dict(method=method, donor_id=donor, n_primary_calls=len(calls),
                n_source_eligible=len(expected), missing_source_UUIDs=0, extra_primary_UUIDs=0,
                duplicate_primary_UUIDs=0, null_classes=0, no_nucleus_survivors=0,
                mandatory_exclusion_survivors=0, calls_path=str(calls_path), calls_sha256=sha(calls_path),
                geometry_sha256=sha(geometry_path), assignment_sha256=sha(record["assignment"]["path"])))
            print("verified full primary universe", method, donor, len(calls), flush=True)
    table = root / "resume_primary_integrity_QC.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    (root / "resume_primary_integrity_QC.json").write_text(json.dumps(dict(
        status="passed", n_donor_method_partitions=len(rows), rows_checked=sum(r["n_primary_calls"] for r in rows),
        result_table=str(table), result_table_sha256=sha(table),
        exact_key_universe_matches_geometry=True, mandatory_exclusions_preserved=True,
        no_nucleus_exclusion_preserved=True, accuracy_validated=False,
        interpretation="Full scope key/QC/mandatory-exclusion integrity; no marker identity or biological accuracy validation"), indent=2) + "\n")

if __name__ == "__main__":
    main()
