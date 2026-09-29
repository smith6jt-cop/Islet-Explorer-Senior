"""Both spellings of QuPath's signed distance, and the v3 per-image context Parquet.

Synthetic, built in ``tmp_path``; reads nothing under ``data/``.
"""
import pathlib
import sys

import duckdb
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import distance  # noqa: E402


def _context(tmp_path, with_word: bool, suffix: str) -> pathlib.Path:
    word = "with " if with_word else ""
    frame = pd.DataFrame({
        "Image": ["img-1", "img-1"],
        "Object ID": ["a", "b"],
        f"Signed distance to annotation {word}Islet µm": [-2.5, 7.0],
        f"Signed distance to annotation {word}Tissue µm": [4.0, None],
    })
    path = tmp_path / f"export{suffix}"
    if suffix == ".parquet":
        frame.drop(columns=["Image"]).to_parquet(path, index=False)
    else:
        frame.to_csv(path, index=False)
    return path


@pytest.mark.parametrize("suffix", [".csv", ".parquet"])
@pytest.mark.parametrize("with_word", [True, False])
def test_both_spellings_are_found_and_extracted(tmp_path, with_word, suffix):
    path = _context(tmp_path, with_word, suffix)
    con = duckdb.connect()
    islet, tissue = distance.assert_columns_present(con, path)
    word = "with " if with_word else ""
    assert islet == f"Signed distance to annotation {word}Islet µm"
    assert tissue == f"Signed distance to annotation {word}Tissue µm"

    out = tmp_path / "out"
    con.execute(distance.extract_sql(path, {"img-1": "6374"}, out,
                                     islet_col=islet, tissue_col=tissue))
    got = pd.read_parquet(out / "donor_id=6374").sort_values("object_id")
    assert got["object_id"].tolist() == ["a", "b"]
    assert got["signed_dist_islet_um"].tolist() == [-2.5, 7.0]
    # A missing measurement stays unavailable, never 0.
    assert got["signed_dist_tissue_um"].isna().tolist() == [False, True]


def test_an_export_with_neither_spelling_is_refused(tmp_path):
    path = tmp_path / "x.csv"
    pd.DataFrame({"Image": ["i"], "Object ID": ["a"]}).to_csv(path, index=False)
    with pytest.raises(RuntimeError, match="annotation with Islet"):
        distance.assert_columns_present(duckdb.connect(), path)


def test_a_per_image_context_parquet_is_refused_for_two_images(tmp_path):
    path = _context(tmp_path, True, ".parquet")
    with pytest.raises(ValueError, match="per-image context Parquet"):
        distance.extract_sql(path, {"i1": "1", "i2": "2"}, tmp_path / "out")
