"""`read_window` must equal `page.asarray()[slice]`, on real TIFFs.

`tiff_window` reaches into tifffile internals -- `dataoffsets`, `databytecounts`, `decode`,
`chunked` -- because a level-0 qptiff plane is ~1.2 GB and the pilot needs 2048-px windows out
of 51 of them. Internals move between tifffile releases, and the failure mode is not an
exception: it is a window that decodes to the wrong pixels and types cells on them.

So these tests write real tiled and stripped TIFFs to `tmp_path` and compare against
`asarray()`. They touch no cohort data and no qptiff; the files are a few KB. A mock page
would only test the mock.
"""
import pathlib
import sys

import numpy as np
import pytest
import tifffile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from tiff_window import read_window, selftest_bit_exact  # noqa: E402

TILE, STRIP = 32, 2          # the qptiff/mask chunk geometries, scaled down
H, W = 200, 168              # not a whole number of chunks in either axis


def _data(dtype=np.uint16):
    """Content where every pixel is distinguishable, so a misplaced chunk cannot alias."""
    y, x = np.mgrid[0:H, 0:W]
    return ((y * W + x) % np.iinfo(dtype).max).astype(dtype)


@pytest.fixture(scope="module")
def tiled(tmp_path_factory):
    """Tiled + LZW, like a qptiff level-0 plane."""
    p = tmp_path_factory.mktemp("tw") / "tiled.tif"
    tifffile.imwrite(p, _data(), tile=(TILE, TILE), compression="lzw")
    with tifffile.TiffFile(p) as tf:
        page = tf.pages[0]
        assert page.is_tiled and page.tilelength == TILE
        yield page


@pytest.fixture(scope="module")
def stripped(tmp_path_factory):
    """Stripped + Deflate, like the rasterized cell mask."""
    p = tmp_path_factory.mktemp("tw") / "stripped.tif"
    tifffile.imwrite(p, _data(np.uint32), rowsperstrip=STRIP, compression="deflate")
    with tifffile.TiffFile(p) as tf:
        page = tf.pages[0]
        assert not page.is_tiled and int(page.tags["RowsPerStrip"].value) == STRIP
        yield page


@pytest.fixture(params=["tiled", "stripped"])
def page(request):
    return request.getfixturevalue(request.param)


# ------------------------------------------------------------ the core contract


def test_selftest_bit_exact_passes(page):
    assert selftest_bit_exact(page, size=48, verbose=False) is True


@pytest.mark.parametrize(
    "box",
    [
        (0, H, 0, W),                       # whole page
        (TILE, TILE * 3, TILE, TILE * 2),   # chunk-aligned
        (7, 61, 13, 97),                    # unaligned on every edge
        (0, 1, 0, 1),                       # single pixel, origin
        (H - 1, H, W - 1, W),               # single pixel, far corner
        (H - 3, H, W - 5, W),               # partial last chunk in both axes
        (5, 6, 0, W),                       # one full row
        (0, H, 11, 12),                     # one full column
    ],
)
def test_matches_asarray(page, box):
    y0, y1, x0, x1 = box
    got = read_window(page, y0, y1, x0, x1)
    ref = page.asarray()[y0:y1, x0:x1]
    assert got.dtype == ref.dtype
    assert got.shape == ref.shape
    assert np.array_equal(got, ref)


def test_every_chunk_boundary_is_exact(page):
    """Sweep windows that straddle each chunk seam -- where an off-by-one would live."""
    ch = page.tilelength if page.is_tiled else int(page.tags["RowsPerStrip"].value)
    full = page.asarray()
    for cy in range(0, H, ch):
        y0, y1 = max(0, cy - 3), min(H, cy + 3)
        assert np.array_equal(read_window(page, y0, y1, 0, W), full[y0:y1, 0:W]), f"row seam {cy}"


# ------------------------------------------------------------ clamping contract


def test_out_of_range_is_clamped_not_wrapped(page):
    """Negative and past-the-end coordinates clamp; numpy would wrap a negative index."""
    got = read_window(page, -10, 5, -10, 5)
    assert np.array_equal(got, page.asarray()[0:5, 0:5])
    got = read_window(page, H - 4, H + 100, W - 4, W + 100)
    assert np.array_equal(got, page.asarray()[H - 4:H, W - 4:W])


@pytest.mark.parametrize("box", [(50, 50, 0, 10), (0, 10, 60, 60), (H + 5, H + 9, 0, 10)])
def test_empty_window_keeps_the_shape_contract(page, box):
    got = read_window(page, *box)
    assert got.dtype == page.dtype
    assert 0 in got.shape


# ------------------------------------------------- sparse chunks read as zeros


class _SparsePage:
    """A page whose chunks are all absent. tifffile writes these for unwritten regions."""

    is_tiled = True
    tilelength = tilewidth = TILE
    shape = (TILE * 2, TILE * 2)
    dtype = np.dtype(np.uint16)
    chunked = (2, 2)
    dataoffsets = (0, 0, 0, 0)
    databytecounts = (0, 0, 0, 0)
    parent = None                       # never reached: the reader opens the handle lazily

    def decode(self, *a, **k):          # pragma: no cover - must not be called
        raise AssertionError("decode() called for a chunk with zero bytes")


def test_sparse_chunks_are_zeros_and_read_nothing():
    got = read_window(_SparsePage(), 0, TILE * 2, 0, TILE * 2)
    assert got.shape == (TILE * 2, TILE * 2)
    assert got.dtype == np.uint16
    assert not got.any()
