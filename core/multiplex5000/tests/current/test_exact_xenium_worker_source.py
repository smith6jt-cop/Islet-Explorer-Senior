"""RNA workers may receive one accepted bundle instead of searching by donor."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from multiplex5000.tissues import lymph_node_xenium_worker as pln
from multiplex5000.tissues import pancreas_xenium_worker as panc
from multiplex5000.xenium_typing_input import XeniumTypingInputError


@pytest.mark.parametrize("worker", [panc, pln])
def test_exact_bundle_wins_when_native_ids_collide(tmp_path, worker, monkeypatch):
    ids = np.array(["shared-1", "shared-2"])
    a, b = tmp_path / "first.parquet", tmp_path / "second.parquet"
    pq.write_table(pa.table({"cell_id": ids, "nucleus_count": [1, 0]}), a)
    pq.write_table(pa.table({"cell_id": ids, "nucleus_count": [0, 2]}), b)
    # An exact-source invocation must not read the historical donor-level index.
    monkeypatch.setattr(worker, "XENIUM_PATHS", tmp_path / "absent-index.csv")
    got, source = worker.nucleus_for("6539", ids, b)
    assert got.tolist() == [0, 2]
    assert str(b) in source
    with pytest.raises(XeniumTypingInputError, match="cell IDs differ"):
        worker.nucleus_for("6539", np.array(["shared-1", "another-cell"]), b)


@pytest.mark.parametrize("worker", [panc, pln])
def test_missing_nucleus_is_not_converted_to_zero(tmp_path, worker):
    cells = tmp_path / "cells.parquet"
    pq.write_table(pa.table({"cell_id": ["a", "b"], "nucleus_count": [1, None]}), cells)
    with pytest.raises(XeniumTypingInputError, match="missing counts; no zero substituted"):
        worker.nucleus_for("6539", np.array(["a", "b"]), cells)


def test_pln_section_output_keeps_acquisition_identity(tmp_path, monkeypatch):
    cells = tmp_path / "cells.parquet"
    h5ad = tmp_path / "raw.h5ad"
    work = tmp_path / "out"
    seen = []
    monkeypatch.setattr(pln, "cluster_module", lambda: "same-package-partition")

    def partition(donor, section, path, mod, **kwargs):
        seen.append((donor, section, path, mod, kwargs))
        return {"n_cells_partitioned": 2}

    def call(donor, section, path, panels, **kwargs):
        seen.append((donor, section, path, panels, kwargs))
        return {"calls": pd.DataFrame({"cell_id": ["a", "b"], "lineage": ["B", "T"]}),
                "clusters": pd.DataFrame({"cluster_id": ["c1"]}),
                "gene_stats": pd.DataFrame({"gene": ["CD3E"]}),
                "report": {"n_cells_in_h5ad": 2}}

    monkeypatch.setattr(pln, "partition_section", partition)
    monkeypatch.setattr(pln, "call_section", call)
    report = pln.process_exact_section("6539", "0040001__PanT", h5ad, cells, {}, work,
                                       acquisition_id="xenium:0040001:PanT", roi="PanT")
    out = pd.read_parquet(work / "0040001__PanT.calls.parquet")
    assert out.cell_uid.tolist() == ["xenium:0040001:PanT:cells:a",
                                     "xenium:0040001:PanT:cells:b"]
    assert out.acquisition_id.eq("xenium:0040001:PanT").all()
    assert all(record[-1]["exact_cells"] == cells and record[-1]["out_root"] == work
               for record in seen)
    assert all(record[-1]["roi"] == "PanT" for record in seen)
    assert report["partition"]["n_cells_partitioned"] == 2
