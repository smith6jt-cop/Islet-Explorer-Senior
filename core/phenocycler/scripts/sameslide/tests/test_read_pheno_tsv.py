"""read_pheno_tsv.py accepts TSV, single-file Parquet, and split cells+context Parquet.

All three input shapes must produce the identical contract-valid table for the same
underlying cells -- that is the whole point of adding Parquet support without touching the
downstream contract. Synthetic fixtures only; no data/ dependency.

    ../../../../.venv/bin/python -m pytest -q core/phenocycler/scripts/sameslide/tests
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve()
LANE = HERE.parents[1]
sys.path.insert(0, str(LANE))

import read_pheno_tsv as R  # noqa: E402


#: One synthetic QuPath export, in the exact column spelling this module matches by pattern.
#: Six cells: two inside Islet_3, one in unannotated "Root object" territory (dropped by
#: default), three plain tissue cells. INS/GCG/SST are carried like any other marker -- the
#: whole point of keeping them (see module docstring) is that this reader must not special-case
#: them away.
def _raw_rows() -> list[dict]:
    return [
        {
            "Object ID": f"obj-{i}",
            "Image": "6533.qptiff",
            "Object type": "Cell",
            "Centroid X µm": 100.0 + i,
            "Centroid Y µm": 200.0 + i,
            "Cell: Area µm^2": 50.0 + i,
            "Parent": parent,
            "Signed distance to annotation with Islet µm": dist_islet,
            "Signed distance to annotation with Tissue µm": dist_tissue,
            "Nucleus: DAPI: Mean": 10.0 + i,
            "Cell: DAPI: Mean": 11.0 + i,
            "Membrane: DAPI: Mean": 12.0 + i,
            "Cytoplasm: DAPI: Mean": 13.0 + i,
            "Cell: INS: Mean": 0.1 * i,
            "Cell: GCG: Mean": 0.2 * i,
            "Cell: SST: Mean": 0.3 * i,
        }
        for i, (parent, dist_islet, dist_tissue) in enumerate([
            ("Islet_3", -5.0, 120.0),
            ("Islet_3", -3.0, 118.0),
            ("Tissue", 30.0, 90.0),
            ("Tissue", 45.0, 80.0),
            ("Tissue", 60.0, 70.0),
            ("Root object (Image)", 500.0, 500.0),
        ])
    ]


def _raw_frame() -> pd.DataFrame:
    return pd.DataFrame(_raw_rows())


#: Columns the v3 cells export carries: measurements only, no geometry/parent context.
_CELLS_COLUMNS = [
    "Object ID", "Cell: Area µm^2",
    "Nucleus: DAPI: Mean", "Cell: DAPI: Mean", "Membrane: DAPI: Mean", "Cytoplasm: DAPI: Mean",
    "Cell: INS: Mean", "Cell: GCG: Mean", "Cell: SST: Mean",
]

#: Columns the v3 context export carries: Object ID, Parent, centroids, signed distances.
_CONTEXT_COLUMNS = [
    "Object ID", "Image", "Object type", "Centroid X µm", "Centroid Y µm", "Parent",
    "Signed distance to annotation with Islet µm", "Signed distance to annotation with Tissue µm",
]


def _write_tsv(path: Path) -> None:
    _raw_frame().to_csv(path, sep="\t", index=False)


def _write_parquet(path: Path, columns: list[str] | None = None) -> None:
    frame = _raw_frame()
    if columns is not None:
        frame = frame[columns]
    frame.to_parquet(path, index=False)


def _sorted(df: pd.DataFrame) -> pd.DataFrame:
    """Row- and column-order-independent view for comparison.

    Column order legitimately differs between the single-file and split-file readers (the
    split reader appends cells-file columns before context-file ones), which is not a
    contract violation -- ``coerce_cell_table`` fixes the order downstream. Only the raw
    ``read_table*`` helpers are compared here, before that canonicalisation runs.
    """
    return df.sort_values("cell_id").reset_index(drop=True).sort_index(axis=1)


def test_tsv_and_single_file_parquet_give_identical_output(tmp_path):
    tsv_path = tmp_path / "6533.tsv"
    parquet_path = tmp_path / "6533.parquet"
    _write_tsv(tsv_path)
    _write_parquet(parquet_path)

    header = R.read_header(tsv_path)
    meta, features, _ = R.plan_columns(header)
    tsv_df, tsv_stats = R.read_table(tsv_path, meta, features)

    p_header = R.read_header_parquet(parquet_path)
    p_meta, p_features, _ = R.plan_columns(p_header)
    parquet_df, parquet_stats = R.read_table_parquet(parquet_path, p_meta, p_features)

    assert tsv_stats["n_cells"] == parquet_stats["n_cells"] == 5   # one Root-object cell dropped
    assert tsv_stats["n_features"] == parquet_stats["n_features"]
    pd.testing.assert_frame_equal(_sorted(tsv_df), _sorted(parquet_df))


def test_tsv_and_split_context_parquet_give_identical_output(tmp_path):
    tsv_path = tmp_path / "6533.tsv"
    cells_path = tmp_path / "cells" / "6533.parquet"
    context_path = tmp_path / "context" / "6533.parquet"
    cells_path.parent.mkdir()
    context_path.parent.mkdir()
    _write_tsv(tsv_path)
    _write_parquet(cells_path, _CELLS_COLUMNS)
    _write_parquet(context_path, _CONTEXT_COLUMNS)

    header = R.read_header(tsv_path)
    meta, features, _ = R.plan_columns(header)
    tsv_df, tsv_stats = R.read_table(tsv_path, meta, features)

    cells_header = R.read_header_parquet(cells_path)
    context_header = R.read_header_parquet(context_path)
    split_meta = R.find_meta_columns(context_header, cells_header)
    split_features = R.find_feature_columns(cells_header)
    joined_df, joined_stats = R.read_table_parquet_joined(cells_path, context_path,
                                                          split_meta, split_features)

    assert tsv_stats["n_cells"] == joined_stats["n_cells"] == 5
    assert tsv_stats["n_features"] == joined_stats["n_features"]
    pd.testing.assert_frame_equal(_sorted(tsv_df), _sorted(joined_df))


def test_object_id_mismatch_fails_loudly(tmp_path):
    cells_path = tmp_path / "cells.parquet"
    context_path = tmp_path / "context.parquet"
    _write_parquet(cells_path, _CELLS_COLUMNS)

    context_frame = _raw_frame()[_CONTEXT_COLUMNS].iloc[:-1]   # drop the last cell's context row
    context_frame.to_parquet(context_path, index=False)

    cells_header = R.read_header_parquet(cells_path)
    context_header = R.read_header_parquet(context_path)
    meta = R.find_meta_columns(context_header, cells_header)
    features = R.find_feature_columns(cells_header)

    with pytest.raises(SystemExit, match="Object ID mismatch"):
        R.read_table_parquet_joined(cells_path, context_path, meta, features)


def test_context_parquet_without_parquet_is_rejected(tmp_path):
    context_path = tmp_path / "context.parquet"
    _write_parquet(context_path, _CONTEXT_COLUMNS)

    with pytest.raises(SystemExit):
        R.main(["--context-parquet", str(context_path), "--donor", "6533", "--serial", "0080870"])


def test_main_writes_the_same_cell_table_from_tsv_and_split_parquet(tmp_path):
    tsv_path = tmp_path / "6533.tsv"
    cells_path = tmp_path / "cells" / "6533.parquet"
    context_path = tmp_path / "context" / "6533.parquet"
    cells_path.parent.mkdir()
    context_path.parent.mkdir()
    _write_tsv(tsv_path)
    _write_parquet(cells_path, _CELLS_COLUMNS)
    _write_parquet(context_path, _CONTEXT_COLUMNS)

    out_tsv = tmp_path / "out_tsv"
    out_parquet = tmp_path / "out_parquet"
    rc_tsv = R.main(["--tsv", str(tsv_path), "--donor", "6533", "--roi", "panc",
                     "--serial", "0080870", "--out", str(out_tsv)])
    rc_parquet = R.main(["--parquet", str(cells_path), "--context-parquet", str(context_path),
                         "--donor", "6533", "--roi", "panc", "--serial", "0080870",
                         "--out", str(out_parquet)])
    assert rc_tsv == 0 and rc_parquet == 0

    written_tsv = pd.read_parquet(out_tsv / "donor_id=6533" / "roi=panc" / "serial=0080870"
                                  / "cells_pheno.parquet")
    written_parquet = pd.read_parquet(out_parquet / "donor_id=6533" / "roi=panc"
                                      / "serial=0080870" / "cells_pheno.parquet")
    pd.testing.assert_frame_equal(_sorted(written_tsv), _sorted(written_parquet))

    # INS/GCG/SST are kept, not stripped -- the whole reason this pair exists.
    kept = {R.as_feature("Cell__INS"), R.as_feature("Cell__GCG"), R.as_feature("Cell__SST")}
    assert kept <= set(written_tsv.columns)
