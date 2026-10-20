#!/usr/bin/env python3
"""Verify prepared compartment values against exact immutable raw-source UUIDs.

The bounded, systematic UUID sample is independent of labels and values. This is
measurement-source integrity evidence, not cell-identity accuracy evidence.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def exact_rows(path, columns, ids):
    need = set(ids)
    found = []
    for batch in pq.ParquetFile(path).iter_batches(batch_size=32768, columns=["object_id", *columns]):
        frame = batch.to_pandas()
        hit = frame.object_id.astype(str).isin(need)
        if hit.any():
            part = frame.loc[hit].copy()
            found.append(part)
            need.difference_update(part.object_id.astype(str))
        if not need:
            break
    if need:
        raise ValueError(f"Missing {len(need)} exact source UUIDs: {path}")
    result = pd.concat(found, ignore_index=True)
    if result.object_id.duplicated().any():
        raise ValueError(f"Duplicate exact source UUIDs: {path}")
    return result.set_index("object_id").reindex(ids)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate-root", type=Path, required=True)
    ap.add_argument("--n-per-donor", type=int, default=2048)
    args = ap.parse_args()
    root = args.candidate_root.resolve()
    manifest_path = root / "source_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    rows = []
    source_entries = []
    for record in manifest["records"]:
        expression_path = Path(record["expression"])
        availability_path = Path(record["availability"])
        availability = pq.ParquetFile(availability_path).read().to_pandas()
        markers = availability.marker.tolist()
        selected = next(pq.ParquetFile(expression_path).iter_batches(
            batch_size=args.n_per_donor, columns=["object_id", *markers])).to_pandas().set_index("object_id")
        ids = selected.index
        originals = {}
        for prefix, key in [("qupath", "cells"), ("redsea", "redsea")]:
            columns = sorted(set(r.source.split(":", 1)[1] for r in availability.itertuples()
                                 if r.available and r.source.startswith(prefix + ":")))
            if columns:
                source_path = Path(record["source_files"][key]["path"])
                originals[prefix] = exact_rows(source_path, columns, ids)
                source_entries.append(dict(arm=record["arm"], donor_id=record["donor_id"],
                    path=str(source_path), recorded_sha256=record["source_files"][key]["sha256"],
                    columns_checked=columns))
        for marker in availability.itertuples():
            actual = selected[marker.marker].to_numpy(float)
            if not marker.available:
                expected = np.full(len(selected), np.nan)
            else:
                prefix, column = marker.source.split(":", 1)
                expected = originals[prefix][column].to_numpy(float)
            same = np.isclose(actual, expected, rtol=1e-7, atol=1e-8, equal_nan=True)
            row = dict(arm=record["arm"], donor_id=record["donor_id"], marker=marker.marker,
                declared_compartment=marker.compartment, acquired=bool(marker.available),
                source=marker.source, n_systematic_UUIDs=len(selected), n_equal=int(same.sum()),
                n_source_missing=int((~np.isfinite(expected)).sum()),
                n_prepared_missing=int((~np.isfinite(actual)).sum()),
                n_source_zero=int((np.isfinite(expected) & (expected == 0)).sum()),
                n_prepared_zero=int((np.isfinite(actual) & (actual == 0)).sum()),
                pass_source_and_missingness=bool(same.all()))
            rows.append(row)
            if not same.all():
                raise ValueError(f"Prepared measurement source mismatch: {row}")
        print("verified source values", record["arm"], record["donor_id"], len(selected), len(markers), flush=True)
    table = root / "resume_expression_source_parity_QC.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    (root / "resume_expression_source_parity_QC.json").write_text(json.dumps(dict(
        status="passed", source_manifest=str(manifest_path), source_manifest_sha256=sha(manifest_path),
        result_table=str(table), result_table_sha256=sha(table),
        sampling="First prepared UUIDs per donor in each comparison, exact UUID join to original recorded source; no label or intensity selection",
        n_per_donor=args.n_per_donor, n_marker_donor_comparisons=len(rows), sources=source_entries,
        all_marker_sources_and_missingness_reproduced=True, accuracy_validated=False,
        interpretation="Bounded measurement-source parity; no biological validation or promotion"), indent=2) + "\n")

if __name__ == "__main__":
    main()
