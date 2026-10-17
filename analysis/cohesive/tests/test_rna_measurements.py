import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from analysis.cohesive.rna_measurements import exact_index, project_expression, read_vector


def make_h5ad(path, *, encoding="csr_matrix", noninteger=False):
    values = np.array([[2, 0, 1], [0, 3, 0], [4, 1, 2]], dtype=float)
    if noninteger:
        values[0, 0] = 0.5
    matrix = {"csr_matrix": sparse.csr_matrix, "csc_matrix": sparse.csc_matrix}[encoding](values)
    strings = h5py.string_dtype("utf-8")
    with h5py.File(path, "w") as handle:
        obs = handle.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        cell_id = obs.create_group("cell_id")
        cell_id.create_dataset("values", data=["c1", "c2", "c3"], dtype=strings)
        cell_id.create_dataset("mask", data=[False, False, False])
        obs.create_dataset("x_centroid", data=[10., 20., 30.])
        obs.create_dataset("y_centroid", data=[1., 2., 3.])
        var = handle.create_group("var")
        var.attrs["_index"] = "gene"
        var.create_dataset("gene", data=["INS", "CCL21", "LTB"], dtype=strings)
        node = handle.create_group("X")
        node.attrs["encoding-type"] = encoding
        node.attrs["shape"] = matrix.shape
        for field in ("data", "indices", "indptr"):
            node.create_dataset(field, data=getattr(matrix, field))


@pytest.mark.parametrize("encoding", ["csr_matrix", "csc_matrix"])
def test_exact_reordered_raw_expression_preserves_zeros_and_unavailable(tmp_path, encoding):
    path = tmp_path / "matrix.h5ad"
    make_h5ad(path, encoding=encoding)
    export = {"scope_id": "pancreas_RNA|serial|section", "donor_id": "donor",
              "source_serial": "serial", "section_id": "section", "lane": "pancreas_RNA",
              "native_expression_matrix_binding": {"path": str(path)}, "native_geometry_bindings": {}}
    result, audit = project_expression(export, pd.DataFrame({"cell_id": ["c3", "c1"]}))
    assert result.cell_id.tolist() == ["c3", "c1"]
    assert result.RNA_INS.tolist() == [4, 2]
    assert result.RNA_CCL21.tolist() == [1, 0]
    assert result.raw_library_counts.tolist() == [7, 3]
    assert result.RNA_CCL19.isna().all()
    assert result.genomic_control_counts.isna().all()
    assert result.x_um.tolist() == [30, 10]
    assert audit["native_matrix_rows"] == 3
    assert audit["new_identity_gate"] is False


def test_normalized_values_are_not_misread_as_raw_counts(tmp_path):
    path = tmp_path / "matrix.h5ad"
    make_h5ad(path, noninteger=True)
    with pytest.raises(ValueError, match="raw integer counts"):
        project_expression({"scope_id": "test", "native_expression_matrix_binding": {"path": str(path)}},
                           pd.DataFrame({"cell_id": ["c1"]}))


def test_joins_reject_duplicate_missing_and_unmatched_keys():
    assert exact_index(["b", "a"], ["a", "b"], "test").tolist() == [1, 0]
    for source, requested in ((["a", "a"], ["a"]), (["a"], ["a", "a"]),
                              ([None, "a"], ["a"]), (["a"], ["missing"])):
        with pytest.raises(ValueError):
            exact_index(source, requested, "test")


def test_nullable_arrays_count_values_not_encoding_members_and_preserve_missingness(tmp_path):
    path = tmp_path / "encoding.h5"
    with h5py.File(path, "w") as handle:
        node = handle.create_group("nullable")
        node.create_dataset("values", data=[1., 2., 3., 4.])
        node.create_dataset("mask", data=[False, True, False, False])
        with pytest.raises(ValueError, match="Null identity"):
            read_vector(node)
        values = read_vector(node, allow_null=True)
        assert len(values) == 4 and np.isnan(values[1])
        assert values[[0, 2, 3]].tolist() == [1., 3., 4.]
