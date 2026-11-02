#!/usr/bin/env python3
"""Inspect every recorded spleen RNA bundle without inventing a typing specification."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def strings(data):
    return [s.decode() if isinstance(s, bytes) else str(s) for s in data]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--cores", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    cores = {p.stem: set(pd.read_csv(p).gene.dropna().astype(str)) for p in args.cores.glob("*.csv")}
    results = []
    for pair in json.loads(args.pairs.read_text())["pairs"]:
        root = Path(pair["xenium_output"])
        matrix = root / "cell_feature_matrix.h5"
        cells_path = root / "cells.parquet"
        with h5py.File(matrix) as data:
            features = strings(data["matrix/features/name"][:])
            types = strings(data["matrix/features/feature_type"][:])
            ids = strings(data["matrix/barcodes"][:])
            shape = data["matrix/shape"][:].tolist()
        names = pq.ParquetFile(cells_path).schema_arrow.names
        selected = [n for n in ("cell_id", "nucleus_count", "nucleus_area", "cell_area") if n in names]
        cells = pq.read_table(cells_path, columns=selected).to_pandas()
        gene_names = {name for name, kind in zip(features, types) if kind == "Gene Expression"}
        nucleus = cells.nucleus_count if "nucleus_count" in cells else None
        result = {"acquisition_id": pair["xenium_acquisition_id"], "physical_section_id": pair["protein_entry"],
                  "cells": len(cells), "matrix_shape": shape, "matrix_hash": sha(matrix), "cells_hash": sha(cells_path),
                  "feature_names": features, "feature_types": dict(pd.Series(types).value_counts()),
                  "gene_panel_sha256": hashlib.sha256(json.dumps(sorted(gene_names), separators=(",", ":")).encode()).hexdigest(),
                  "gene_count": len(gene_names), "cell_id_unique": bool(cells.cell_id.is_unique),
                  "matrix_cell_ids_match": set(ids) == set(cells.cell_id),
                  "nucleus_count_available": nucleus is not None,
                  "nucleus_eligible": int(np.count_nonzero(nucleus > 0)) if nucleus is not None else None,
                  "no_nucleus": int(np.count_nonzero(nucleus <= 0)) if nucleus is not None else None,
                  "core_coverage": {name: {"present": sorted(genes & gene_names), "absent": sorted(genes - gene_names)}
                                    for name, genes in cores.items()},
                  "typing_status": "blocked: no full spleen-specific independent RNA lineage specification/route"}
        # pandas integer scalars are normalized by JSON default below.
        results.append(result)
        print(pair["xenium_acquisition_id"], len(cells), "genes", len(gene_names), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"pairs_sha256": sha(args.pairs),
                                   "core_hashes": {p.name: sha(p) for p in args.cores.glob("*.csv")},
                                   "scope": "Inputs and generic immune-core coverage only; no transferred labels or biological acceptance",
                                   "acquisitions": results}, indent=2, default=lambda x: x.item()) + "\n")


if __name__ == "__main__":
    main()
