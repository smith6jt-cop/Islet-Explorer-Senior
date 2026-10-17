"""Exact native expression joins, retaining raw counts and missing measurements.

The H5AD decoding follows the existing regenerate_current_rna_umaps reader.
No embedding, identity model, or additional cell-admission rule is introduced.
"""
from __future__ import annotations

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse

from .common import require


GENES = ("INS", "GCG", "SST", "CHGA", "CD3E", "TRAC", "CD79A", "MS4A1",
         "LST1", "CD163", "COL1A1", "VWF", "CCL19", "CCL21", "LTB",
         "TCF7", "LEF1", "TOX", "GZMB", "PRF1", "NKG7", "FOXP3")


def read_vector(node, *, allow_null: bool = False) -> np.ndarray:
    if isinstance(node, h5py.Dataset):
        values = node[()]
        missing = np.zeros(len(values), dtype=bool)
    elif "values" in node:
        values = node["values"][()]
        missing = node["mask"][()] if "mask" in node else np.zeros(len(values), dtype=bool)
    elif "categories" in node and "codes" in node:
        categories = read_vector(node["categories"])
        codes = node["codes"][()]
        missing = codes < 0
        require(allow_null or not missing.any(), f"Null categorical identity at {node.name}")
        require(not ((codes >= len(categories)) & ~missing).any(), f"Invalid category code at {node.name}")
        values = np.empty(len(codes), dtype=categories.dtype)
        values[~missing] = categories[codes[~missing]]
    else:
        raise ValueError(f"Unsupported H5AD array encoding: {node.name}")
    require(values.ndim == 1, f"Expected a vector at {node.name}")
    require(allow_null or not missing.any(), f"Null identity/feature values at {node.name}")
    if values.dtype.kind in "OSU":
        values = np.asarray([value.decode() if isinstance(value, bytes) else str(value)
                             for value in values], dtype=object)
    if missing.any():
        values = values.astype(float if values.dtype.kind in "bifu" else object)
        values[missing] = np.nan if values.dtype.kind == "f" else None
    return values


def exact_index(source, requested, description: str) -> np.ndarray:
    index = pd.Index(source)
    require(not index.hasnans and index.is_unique, f"Missing/duplicate source keys: {description}")
    requested = pd.Index(requested)
    require(not requested.hasnans and requested.is_unique, f"Missing/duplicate requested keys: {description}")
    positions = index.get_indexer(requested)
    require((positions >= 0).all(), f"Unmatched exact keys: {description}")
    return positions


def project_expression(export: dict, annotation: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    matrix_binding = export["native_expression_matrix_binding"]
    path = matrix_binding["path"]
    with h5py.File(path, "r") as handle:
        obs, var = handle["obs"], handle["var"]
        cells = read_vector(obs[obs.attrs["_index"]])
        genes = read_vector(var[var.attrs["_index"]])
        require(pd.Index(genes).is_unique, f"Duplicate native feature names: {path}")
        positions = exact_index(cells, annotation.cell_id, export["scope_id"])
        node = handle["X"]
        encoding = node.attrs["encoding-type"]
        shape = tuple(int(value) for value in node.attrs["shape"])
        require(shape == (len(cells), len(genes)), f"Native matrix dimension mismatch: {path}")
        constructors = {"csr_matrix": sparse.csr_matrix, "csc_matrix": sparse.csc_matrix}
        require(encoding in constructors, f"Unsupported native sparse matrix: {encoding}")
        data = node["data"][:]
        require(np.isfinite(data).all() and (data >= 0).all() and np.equal(data, np.round(data)).all(),
                f"Native X is not finite, nonnegative raw integer counts: {path}")
        matrix = constructors[encoding]((data, node["indices"][:], node["indptr"][:]), shape=shape)
        matrix.check_format(full_check=True)
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        library = np.asarray(matrix.sum(axis=1)).ravel()[positions]
        require((library > 0).all(), f"An admitted cell has no raw gene library: {export['scope_id']}")
        result = annotation.copy()
        result["raw_library_counts"] = library
        result["n_genes_detected"] = np.asarray(matrix.getnnz(axis=1)).ravel()[positions]
        available = [gene for gene in GENES if gene in genes]
        missing = [gene for gene in GENES if gene not in genes]
        columns = pd.Index(genes).get_indexer(available)
        projected = matrix[:, columns].tocsr()[positions].toarray()
        for column, gene in enumerate(available):
            result["RNA_" + gene] = projected[:, column]
        for gene in missing:
            result["RNA_" + gene] = np.nan

        coordinate_source = None
        if "x_centroid" in obs and "y_centroid" in obs:
            result["x_um"] = read_vector(obs["x_centroid"], allow_null=True)[positions]
            result["y_um"] = read_vector(obs["y_centroid"], allow_null=True)[positions]
            coordinate_source = "bound native H5AD obs centroids"
        elif "spatial" in handle.get("obsm", {}):
            xy = handle["obsm/spatial"][:]
            require(xy.shape == (len(cells), 2), f"Invalid native spatial dimensions: {path}")
            result["x_um"], result["y_um"] = xy[positions].T
            coordinate_source = "bound native H5AD obsm/spatial"

        observed_fields = {}
        for source, target in (("nucleus_area", "nucleus_area_um2"), ("cell_area", "cell_area_um2"),
                               ("genomic_control_counts", "genomic_control_counts"),
                               ("control_probe_counts", "control_probe_counts"),
                               ("control_codeword_counts", "control_codeword_counts")):
            if source in obs:
                result[target] = read_vector(obs[source], allow_null=True)[positions]
                observed_fields[target] = f"native H5AD obs/{source}"

    vendor = export["native_geometry_bindings"].get("vendor_cell_table")
    if vendor is not None:
        pf = pq.ParquetFile(vendor["path"])
        fields = {"x_centroid": "x_um", "y_centroid": "y_um", "cell_area": "cell_area_um2",
                  "nucleus_area": "nucleus_area_um2", "genomic_control_counts": "genomic_control_counts",
                  "control_probe_counts": "control_probe_counts", "control_codeword_counts": "control_codeword_counts"}
        columns = [field for field in fields if field in pf.schema_arrow.names]
        frame = pf.read(columns=["cell_id", *columns]).to_pandas()
        indices = exact_index(frame.cell_id, result.cell_id, export["scope_id"] + " vendor geometry")
        for source in columns:
            target = fields[source]
            values = frame.iloc[indices][source].to_numpy()
            if target in result:
                require(np.allclose(result[target].to_numpy(float), values.astype(float), equal_nan=True),
                        f"Bound native/vendor {target} disagree: {export['scope_id']}")
            else:
                result[target] = values
            observed_fields[target] = "bound exact vendor cell table"
        if "x_um" in result and "y_um" in result:
            coordinate_source = "bound native/vendor exact-key geometry"

    for field in ("x_um", "y_um", "nucleus_area_um2", "cell_area_um2",
                  "genomic_control_counts", "control_probe_counts", "control_codeword_counts"):
        if field not in result:
            result[field] = np.nan
    result["coordinates_available"] = np.isfinite(result[["x_um", "y_um"]].to_numpy(float)).all(axis=1)
    result["scope_id"] = export["scope_id"]
    result["donor_id"] = export["donor_id"]
    result["source_serial"] = export["source_serial"]
    result["section_id"] = export["section_id"]
    result["lane"] = export["lane"]
    return result, {
        "scope_id": export["scope_id"], "native_matrix_rows": len(cells),
        "native_features": len(genes), "biological_rows": len(result),
        "available_genes": available, "unavailable_genes": missing,
        "raw_counts_integer_verified": True, "new_identity_gate": False,
        "coordinate_source": coordinate_source,
        "rows_with_finite_coordinates": int(result.coordinates_available.sum()),
        "measured_covariate_sources": observed_fields,
    }
