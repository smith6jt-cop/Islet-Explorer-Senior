"""The phenotyped base, and the two silent traps that nearly made it meaningless.

Both were caught by checking the round-trip rather than the write, and neither raised at the
point of damage. They are pinned here because the failure mode is a typing run that completes,
reports plausible numbers, and means nothing.
"""
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import build_phenotyped_base as B                       # noqa: E402


def test_string_columns_are_object_not_pandas3_string_dtype():
    """anndata < 0.13 REFUSES to write a nullable StringArray.

    Under pandas 3.0, `pd.Index(np.array([...], dtype=object))` returns an
    ArrowStringArray -- so the obvious construction produces a file the typer's own
    environment (scvi-env, anndata 0.12.4) cannot write. The fix is to pin the dtype, not to
    toggle `settings.allow_write_nullable_strings`: toggling it would make the file
    unwritable *there*, which is the one place it has to work.
    """
    idx = B._object_index(np.array(["aaa-1", "bbb-1"], dtype=object))
    assert idx.dtype == object, f"index dtype is {idx.dtype!r}, not object"
    assert isinstance(idx._values, np.ndarray)
    assert idx.name is None, "a NAMED index is the other trap -- see the round-trip test"

    cat = B._object_categorical(["panc", "panc", "pln"])
    assert cat.categories.dtype == object, f"categories are {cat.categories.dtype!r}"


def test_the_indices_survive_a_second_writer():
    """The trap that would have been invisible.

    `build_phenotyped_base` writes the h5ad; `finalize_panel.py` then reads it, stamps the
    panel masks and writes it again, in a DIFFERENT anndata. With a *named* index, that
    second write reset both indices to positional strings -- `var_names` became
    '0','1','2'... Every gene lookup then missed silently: the identity panel came back
    0/5101 and donor_sex 'unknown', with no error anywhere. scANVI would have trained on a
    feature set of nothing.

    This reproduces the two-writer round-trip on a tiny synthetic object rather than a real
    donor, so it runs in milliseconds and needs no bundle.
    """
    ad = pytest.importorskip("anndata")
    import scipy.sparse as sp

    genes = ["INS", "GCG", "KDM5D"]
    cells = ["aaaa-1", "bbbb-1"]
    X = sp.csr_matrix(np.array([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]], dtype=np.float32))
    a = ad.AnnData(X=X,
                   obs=pd.DataFrame(index=B._object_index(cells)),
                   var=pd.DataFrame(index=B._object_index(genes)))
    a.layers["counts"] = X
    a.obs["roi"] = B._object_categorical(["panc", "panc"])

    tmp = pathlib.Path(pytest.__file__).parent  # placeholder, replaced below
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = pathlib.Path(d) / "x.h5ad"
        a.write_h5ad(p)                     # writer 1 -- this script
        b = ad.read_h5ad(p)
        b.var["panel_for_embedding"] = b.var_names.isin(genes)
        b.write_h5ad(p)                     # writer 2 -- finalize_panel.py
        c = ad.read_h5ad(p)

    assert list(c.var_names) == genes, (
        f"gene names did not survive two writers: {list(c.var_names)}. "
        "This is the failure that returns identity_panel = 0/N with no error."
    )
    assert list(c.obs_names) == cells
    assert int(c.var["panel_for_embedding"].sum()) == len(genes)
    del tmp


def test_the_qc_filter_is_density_based_on_the_upper_bound():
    """Decision #10: absolute LOWER bounds, DENSITY upper bound.

    An absolute counts cap predominantly removes large *healthy* cells, because cell_area
    varies 5-10x in FFPE Xenium and counts track area. The test builds one cell that is
    large and bright at normal density and one that is small and equally bright: only the
    second is a density outlier, and an absolute cap would have taken the wrong one.
    """
    import scipy.sparse as sp
    n_genes = 40
    rows = []
    for total in (400, 400, 60, 60, 60, 60, 60, 60, 60, 60):
        r = np.zeros(n_genes); r[: n_genes // 2] = total / (n_genes // 2)
        rows.append(r)
    counts = sp.csr_matrix(np.array(rows))
    # cell 0 is large (normal density); cell 1 is small (outlier density)
    area = np.array([400.0, 20.0] + [60.0] * 8)

    keep, qc = B.qc_mask(counts, area)
    assert keep[0], "a large bright cell at normal density must be kept"
    assert not keep[1], "a small equally-bright cell is the density outlier"
    assert qc["min_genes"] == B.MIN_GENES and qc["min_counts"] == B.MIN_COUNTS
    assert qc["density_quantile"] == B.DENSITY_QUANTILE


def test_the_serial_is_read_from_the_bundle_name_and_refuses_a_path_without_one():
    """The serial is the sample id every downstream consumer keys on, so it is parsed, never
    guessed, and an unparseable path is an error rather than a fabricated id."""
    assert B.serial_of(
        "HuPrime5kPanTissue/Y9QRJ6/MM044_Xenium_nPODBatch10/"
        "output-XETG00298__0059865__Panc__20250328__202338") == "0059865"
    with pytest.raises(B.BuildError):
        B.serial_of("some/path/without/a/serial")
