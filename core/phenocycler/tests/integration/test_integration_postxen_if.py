"""S1c after the 2026-09-05 repoint: the re-stain lane's per-cell IF values reach the Xenium table as
protein evidence, aligned by cell_id through the bundle order, NaN preserved; the endocrine_subtype
stamp stays off (Q4)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.config import PipelineConfig
from phenocycler.integration import postxen
from phenocycler.integration.contract import partition_path, read_cell_table, write_cell_table


@pytest.fixture
def cfg(tmp_path):
    return PipelineConfig(data_dir=tmp_path)


@pytest.fixture
def xen_geometry(cfg, make_contract_table):
    df = make_contract_table(n=6, modality="xenium", donor_id="6436", roi="panc", profile="geometry")
    write_cell_table(df, partition_path(cfg.cells_xen_dir, "6436", "panc"))
    return df


def _bundle_and_arrays(tmp_path, ids, values_ins, values_gcg):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    pd.DataFrame({"cell_id": ids}).to_parquet(bundle / "cells.parquet", index=False)
    if_dir = tmp_path / "protein_rna"
    if_dir.mkdir()
    np.save(if_dir / "6436_0059907_cell_INS.npy", np.asarray(values_ins, dtype=float))
    np.save(if_dir / "6436_0059907_cell_GCG.npy", np.asarray(values_gcg, dtype=float))
    return bundle / "cells.parquet", if_dir


def test_if_values_are_stamped_by_cell_id_through_bundle_order_and_nan_is_kept(cfg, xen_geometry, tmp_path):
    ids = xen_geometry["cell_id"].astype(str).tolist()
    order = ids[::-1]                                   # the bundle lists cells in another order
    ins = [10.0, 20.0, np.nan, 40.0, 50.0, 60.0]        # one cell outside the imaged area
    gcg = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    bundle_cells, if_dir = _bundle_and_arrays(tmp_path, order, ins, gcg)
    res = postxen.stamp_if_values(cfg, "6436", "panc", if_dir=if_dir, bundle_cells=bundle_cells)
    assert res["markers"] == ["INS", "GCG"] and "error" not in res
    got = read_cell_table(cfg.cells_xen_dir, "6436", "panc", modality="xenium").df
    want_ins = dict(zip(order, ins))
    for cid, v in zip(got["cell_id"].astype(str), got["feat__if_INS"]):
        if np.isnan(want_ins[cid]):
            assert np.isnan(v)
        else:
            assert v == pytest.approx(want_ins[cid])
    assert (got["feat__if_GCG"].to_numpy() == np.array([dict(zip(order, gcg))[c] for c in got["cell_id"].astype(str)])).all()
    assert "endocrine_subtype" not in got.columns


def test_a_length_mismatch_is_refused_not_padded(cfg, xen_geometry, tmp_path):
    ids = xen_geometry["cell_id"].astype(str).tolist()
    bundle_cells, if_dir = _bundle_and_arrays(tmp_path, ids, [1.0] * 5, [1.0] * 5)   # 5 values, 6 cells
    with pytest.raises(ValueError, match="refusing to align by position"):
        postxen.stamp_if_values(cfg, "6436", "panc", if_dir=if_dir, bundle_cells=bundle_cells)


def test_endocrine_subtype_stamp_is_off_by_default():
    assert postxen.STAMP_ENDOCRINE_SUBTYPE is False
