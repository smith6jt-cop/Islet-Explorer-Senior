"""The portable Xenium detection predicate keeps unavailable distinct from zero."""

import numpy as np
import pytest
from scipy import sparse

from multiplex5000.rna_evidence import RNAEvidenceError, detect_genes


def test_fixed_depth_detection_and_missing_gene_are_tristate():
    counts = sparse.csr_matrix(np.array([
        [1, 9999],    # one count meets the 1e4 target exactly
        [1, 19999],   # one count at twice that depth does not
        [2, 19998],   # two counts do
        [4, 6],       # measured, but no nucleus
        [0, 0],       # zero depth is unavailable
        [0, 5],       # assayed and eligible: zero is negative evidence
    ], dtype=np.int64))
    evidence = detect_genes(counts, ("CD3E", "OTHER"),
                            ("CD3E", "TRDC", "OTHER"),
                            np.array([True, True, True, False, True, True]))
    assert evidence.states["CD3E"].tolist() == [1, 0, 1, -1, -1, 0]
    assert evidence.states["OTHER"].tolist() == [1, 1, 1, -1, -1, 1]
    assert evidence.states["TRDC"].tolist() == [-1] * 6
    assert evidence.missing_genes == ("TRDC",)
    assert evidence.eligible.tolist() == [True, True, True, False, False, True]


def test_short_comparison_matches_recorded_lognorm_predicate():
    raw = np.array([0, 1, 2, 3, 5, 100], dtype=np.int64)
    other = np.array([0, 9_999, 19_998, 29_997, 49_995, 999_900], dtype=np.int64)
    counts = sparse.csr_matrix(np.column_stack([raw, other]))
    result = detect_genes(counts, ("GENE", "OTHER"), ("GENE",),
                          np.ones(len(raw), dtype=bool))
    depth = raw + other
    expected = np.full(len(raw), -1, dtype=np.int8)
    usable = depth > 0
    expected[usable] = (np.log1p(raw[usable] * 1e4 / depth[usable])
                        >= np.log1p(1.0)).astype(np.int8)
    np.testing.assert_array_equal(result.states["GENE"], expected)


def test_selected_columns_keep_requested_order_and_all_missing_is_unavailable():
    counts = sparse.csr_matrix([[2, 0, 3], [0, 5, 4]])
    genes = ("A", "B", "C")
    result = detect_genes(counts, genes, ("C", "A"), np.array([True, True]))
    assert list(result.states) == ["C", "A"]
    assert result.states["C"].tolist() == [1, 1]
    assert result.states["A"].tolist() == [1, 0]
    missing = detect_genes(counts, genes, ("NOT_ASSAYED",), np.array([True, True]))
    assert missing.states["NOT_ASSAYED"].tolist() == [-1, -1]


@pytest.mark.parametrize("counts,genes,requested,nuclei", [
    (sparse.csr_matrix([[1.0]]), ("CD3E",), ("CD3E",), np.array([True])),
    (sparse.csr_matrix([[-1]]), ("CD3E",), ("CD3E",), np.array([True])),
    (sparse.csr_matrix([[1, 2]]), ("CD3E", "CD3E"), ("CD3E",), np.array([True])),
    (sparse.csr_matrix([[1]]), ("CD3E",), ("CD3E", "CD3E"), np.array([True])),
    (sparse.csr_matrix([[1]]), ("CD3E",), ("CD3E",), np.array([1])),
])
def test_refuses_inputs_that_would_make_evidence_ambiguous(counts, genes,
                                                            requested, nuclei):
    with pytest.raises(RNAEvidenceError):
        detect_genes(counts, genes, requested, nuclei)
