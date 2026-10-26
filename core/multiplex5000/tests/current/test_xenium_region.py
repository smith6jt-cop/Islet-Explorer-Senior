"""Raw Xenium region reads keep cell identity and nucleus eligibility explicit."""

from __future__ import annotations

from pathlib import Path

import h5py
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from multiplex5000.xenium_region import XeniumRegionError, load_region
from multiplex5000.xenium_sources import XeniumSource


def _region(tmp_path: Path, serial: str, roi: str, *, table_ids=None,
            nuclei=(1, 0)) -> XeniumSource:
    folder = tmp_path / f"output-XETG00298__{serial}__{roi}__20261026__120000"
    folder.mkdir()
    with h5py.File(folder / "cell_feature_matrix.h5", "w") as file:
        matrix = file.create_group("matrix")
        matrix.create_dataset("shape", data=[3, 2])
        matrix.create_dataset("barcodes", data=[b"shared-1", b"other-1"])
        features = matrix.create_group("features")
        features.create_dataset("name", data=[b"CD3E", b"control", b"TRDC"])
        features.create_dataset("feature_type", data=[b"Gene Expression",
                                                       b"Negative Control Probe",
                                                       b"Gene Expression"])
        # CSC(features x cells) is CSR(cells x features): [2,9,1], [3].
        matrix.create_dataset("data", data=[2, 9, 1, 3])
        matrix.create_dataset("indices", data=[0, 1, 2, 2])
        matrix.create_dataset("indptr", data=[0, 3, 4])
    ids = table_ids if table_ids is not None else ["other-1", "shared-1"]
    pq.write_table(pa.table({"cell_id": ids, "x_centroid": [20.0, 10.0],
                             "y_centroid": [40.0, 30.0], "nucleus_count": list(nuclei)}),
                   folder / "cells.parquet")
    return XeniumSource(f"xenium:{serial}:{roi}", f"{serial}__{roi}", "6539", roi,
                        "pancreas" if roi == "panc" else "pancreatic_lymph_node",
                        folder.name, str(folder), str(folder), "ready", (), "")


def test_region_join_keeps_raw_counts_and_excludes_no_nucleus(tmp_path):
    source = _region(tmp_path, "0059865", "panc")
    region = load_region(source)
    assert region.genes == ("CD3E", "TRDC")
    assert region.counts.toarray().tolist() == [[2, 1], [0, 3]]
    assert region.cell_id.tolist() == ["shared-1", "other-1"]
    assert region.cell_uid.tolist() == ["xenium:0059865:panc:cells:shared-1",
                                         "xenium:0059865:panc:cells:other-1"]
    assert region.x_um.tolist() == [10.0, 20.0]
    assert region.nucleus_count.tolist() == [0, 1]
    assert region.eligible.tolist() == [False, True]


def test_same_native_id_on_second_region_stays_distinct(tmp_path):
    pancreas = load_region(_region(tmp_path, "0059865", "panc"))
    pln = load_region(_region(tmp_path, "0059865", "pln"))
    assert pancreas.cell_id[0] == pln.cell_id[0]
    assert pancreas.cell_uid[0] != pln.cell_uid[0]


def test_missing_or_repeated_cell_ids_cannot_join(tmp_path):
    source = _region(tmp_path, "0059865", "panc", table_ids=["other-1", "lost-1"])
    with pytest.raises(XeniumRegionError, match="cell IDs differ"):
        load_region(source)
    pq.write_table(pa.table({"cell_id": ["shared-1", "shared-1"],
                             "x_centroid": [1.0, 2.0], "y_centroid": [1.0, 2.0],
                             "nucleus_count": [1, 1]}),
                   Path(source.source_path) / "cells.parquet")
    with pytest.raises(XeniumRegionError, match="duplicate cells.parquet"):
        load_region(source)


def test_unready_region_is_refused(tmp_path):
    from dataclasses import replace

    source = _region(tmp_path, "0059865", "panc")
    with pytest.raises(XeniumRegionError, match="not ready"):
        load_region(replace(source, status="remote_only"))
