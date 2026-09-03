"""
Tests for utils.xenium_explorer_export.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

anndata = pytest.importorskip("anndata")
zarr = pytest.importorskip("zarr")

# ``utils/__init__.py`` lazy-imports, so this works even when scanpy is not
# installed (which is the case in GitHub Actions CI).
from utils.xenium_explorer_export import (  # noqa: E402
    export_for_xenium_explorer,
    export_groups_to_csv,
    export_groups_to_zarr,
    generate_color_palette,
)


@pytest.fixture
def adata():
    rng = np.random.default_rng(0)
    n = 20
    obs = pd.DataFrame({
        "cell_id": [f"cell-{i}" for i in range(n)],
        "leiden": pd.Categorical(rng.integers(0, 3, size=n).astype(str)),
        "phenotype": pd.Categorical(
            rng.choice(["T_cell", "B_cell", "Fibroblast"], size=n)
        ),
        "int_cluster": rng.integers(0, 4, size=n),
    }, index=[f"cell-{i}" for i in range(n)])
    X = rng.random((n, 5)).astype("float32")
    return anndata.AnnData(X=X, obs=obs)


def test_generate_color_palette_returns_hex():
    colors = generate_color_palette(["a", "b", "c"])
    assert set(colors) == {"a", "b", "c"}
    for v in colors.values():
        assert v.startswith("#") and len(v) == 7


def test_export_groups_to_csv_roundtrip(adata, tmp_path):
    out = tmp_path / "cell_groups.csv"
    export_groups_to_csv(
        adata,
        group_keys=["leiden", "phenotype"],
        output_path=out,
        cell_id_key="cell_id",
    )
    df = pd.read_csv(out)
    assert list(df.columns) == ["cell_id", "leiden", "phenotype"]
    assert len(df) == adata.n_obs
    assert df["cell_id"].tolist() == adata.obs["cell_id"].tolist()


def test_export_groups_to_csv_rename(adata, tmp_path):
    out = tmp_path / "cell_groups.csv"
    export_groups_to_csv(
        adata,
        group_keys=["leiden"],
        output_path=out,
        rename={"leiden": "Clusters"},
    )
    df = pd.read_csv(out)
    assert "Clusters" in df.columns
    assert "leiden" not in df.columns


def test_export_groups_to_csv_numeric_column(adata, tmp_path):
    out = tmp_path / "cell_groups.csv"
    export_groups_to_csv(adata, group_keys="int_cluster", output_path=out)
    df = pd.read_csv(out)
    assert df["int_cluster"].dtype == object or pd.api.types.is_integer_dtype(df["int_cluster"])


def test_export_groups_to_csv_missing_key(adata, tmp_path):
    with pytest.raises(KeyError):
        export_groups_to_csv(adata, group_keys="no_such", output_path=tmp_path / "x.csv")


def _open_zarr_zip(path):
    # zarr.ZipStore lives in zarr.storage in v3.
    ZipStore = getattr(zarr, "ZipStore", None)
    if ZipStore is None:
        from zarr.storage import ZipStore as ZS
        return ZS(str(path), mode="r")
    return ZipStore(str(path), mode="r")


def test_export_groups_to_zarr_structure(adata, tmp_path):
    out = tmp_path / "analysis.zarr.zip"
    export_groups_to_zarr(
        adata,
        group_keys=["leiden", "phenotype"],
        output_path=out,
    )
    assert out.exists()

    store = _open_zarr_zip(out)
    try:
        root = zarr.open(store, mode="r")
        cg = root["cell_groups"]
        assert list(cg.attrs["grouping_names"]) == ["leiden", "phenotype"]
        assert len(cg.attrs["group_names"]) == 2

        for idx in range(2):
            sub = cg[str(idx)]
            indices = np.asarray(sub["indices"])
            indptr = np.asarray(sub["indptr"])
            # indices cover every cell exactly once
            assert indices.size == adata.n_obs
            assert sorted(indices.tolist()) == list(range(adata.n_obs))
            # indptr sums to len(indices)
            assert indptr[-1] == indices.size
            assert len(sub.attrs["color_palette"]) == len(cg.attrs["group_names"][idx])
    finally:
        store.close()


def test_export_groups_to_zarr_appends_suffix(adata, tmp_path):
    out = tmp_path / "analysis"
    written = export_groups_to_zarr(adata, group_keys="leiden", output_path=out)
    assert written.name.endswith(".zarr.zip")
    assert written.exists()


def test_export_for_xenium_explorer_writes_all_artifacts(adata, tmp_path):
    artifacts = export_for_xenium_explorer(
        adata,
        group_keys=["leiden", "phenotype"],
        output_dir=tmp_path,
        sample_name="sampleA",
    )
    assert set(artifacts) == {"csv", "zarr", "palette"}
    for p in artifacts.values():
        assert p.exists()

    palette = json.loads(artifacts["palette"].read_text())
    assert set(palette) == {"leiden", "phenotype"}
    # every value is a mapping of category -> hex color
    for cats in palette.values():
        for color in cats.values():
            assert color.startswith("#") and len(color) == 7
