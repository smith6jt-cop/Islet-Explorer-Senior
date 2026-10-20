"""nullgate_fast must be the same gate, not a similar one."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))                 # hpg/
sys.path.insert(0, str(HERE.parents[1]))             # xenium_celltype/

import nullgate                                       # noqa: E402
import nullgate_fast                                  # noqa: E402


def _fixture(n_cells=3000, n_genes=400, seed=3):
    rng = np.random.default_rng(seed)
    X = sp.random(n_cells, n_genes, density=0.05, format="csr", random_state=seed)
    genes = [f"G{i:04d}" for i in range(n_genes)]
    groups = rng.integers(0, 9, n_cells)
    # one group is deliberately tiny so the too_few_cells branch, and its effect on the RNG
    # stream, is exercised
    groups[:20] = 99
    cores = {"A": genes[:30], "B": genes[100:180], "C": genes[300:317]}
    return X, genes, groups, cores


def test_contrast_decomposition_is_exact():
    X, genes, groups, _ = _fixture()
    labels, sizes, delta = nullgate_fast.group_gene_deltas(X, groups)
    rng = np.random.default_rng(0)
    for g in labels[:4]:
        cols = sorted(rng.choice(X.shape[1], 41, replace=False).tolist())
        ref = nullgate._contrast(X, cols, groups == g)
        fast = float(delta[list(labels).index(g)][cols].mean())
        assert ref == pytest.approx(fast, rel=0, abs=1e-12)


def test_same_frame_as_reference():
    X, genes, groups, cores = _fixture()
    a = nullgate.annotate_groups(X, genes, groups, cores, n_draws=25, seed=7)
    b = nullgate_fast.annotate_groups(X, genes, groups, cores, n_draws=25, seed=7)
    pd.testing.assert_frame_equal(a, b, check_exact=False, rtol=1e-10, atol=1e-12)


def test_group_covering_every_cell_is_refused():
    X, genes, _, cores = _fixture()
    with pytest.raises(nullgate.NullGateError):
        nullgate_fast.annotate_groups(X, genes, np.zeros(X.shape[0], int), cores, n_draws=3)


def test_lineage_with_no_panel_gene_is_refused():
    X, genes, groups, _ = _fixture()
    with pytest.raises(nullgate.NullGateError):
        nullgate_fast.annotate_groups(X, genes, groups, {"Ghost": ["NOT_ON_PANEL"]}, n_draws=3)
