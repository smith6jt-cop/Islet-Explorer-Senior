#!/usr/bin/env python3
"""Recover the implicit pLN embedding target from bound Stage 04 counts/QC.

This records a reconstructed fitted value; it does not rerun clustering or alter calls.
The fixed 10000 target used for identity detection is a different calculation.
"""
import argparse
import hashlib
import inspect
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    calls = pd.read_csv(args.root / "06_typing/rna_acquisitions.csv", keep_default_na=False)
    for row in calls[calls.tissue == "pancreatic_lymph_node"].to_dict("records"):
        p = args.root / f"04_rna_ingest/cells/{row['sample_key']}.parquet"
        cells = pq.read_table(p, columns=["gene_counts", "detected_genes", "nucleus_count"]).to_pandas()
        keep = (cells.nucleus_count > 0) & (cells.gene_counts >= 10) & (cells.detected_genes >= 5)
        target = float(cells.loc[keep, "gene_counts"].median()) if row["status"] == "candidate" else None
        report = args.root / f"06_typing/rna/{row['sample_key']}/{row['sample_key']}.report.json"
        expected = json.loads(report.read_text())["partition"]["n_cells_partitioned"] if report.exists() else None
        if expected is not None and expected != int(keep.sum()):
            raise ValueError(f"{row['acquisition_id']}: reconstructed partition universe differs")
        rows.append({"acquisition_id": row["acquisition_id"], "status": row["status"],
                     "normalization_target": target, "unit": "gene counts/cell",
                     "partition_cells": int(keep.sum()) if expected is not None else None,
                     "source": str(p), "reconstruction": "median gene_counts after nucleus and 10-count/5-feature QC; default target_sum=None in normalize_total"})
    source = inspect.getsource(sc.pp.normalize_total)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"scanpy_version": sc.__version__,
                                   "normalizer_signature": str(inspect.signature(sc.pp.normalize_total)),
                                   "normalizer_source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                                   "acquisitions": rows}, indent=2) + "\n")
    args.out.with_suffix(".source.py").write_text(source)


if __name__ == "__main__":
    main()
