"""Windowed reads from a tiled or stripped TIFF, without zarr.

``tifffile`` can only hand back a whole page without the optional ``zarr``
dependency (``aszarr``/``imread(selection=...)`` both raise ``ModuleNotFoundError``
here, and ``aszarr`` is deprecated). A whole level-0 qptiff plane is ~1.2 GB, and
the pilot needs 2048-px windows out of 51 of them, so we decode only the chunks a
window touches -- 16 tiles per window on this data instead of 2,303.

This leans on tifffile internals (``dataoffsets``, ``databytecounts``,
``decode``, ``chunked``), so ``selftest_bit_exact`` must pass before the reader is
trusted: it checks ``read_window`` against ``page.asarray()[slice]`` at interior,
unaligned, top-left, bottom-right and partial-last-chunk positions.

Each worker process must open its OWN ``TiffFile``. Under the default ``fork``
start method a shared handle means a shared file offset, and concurrent seeks
interleave into corrupt reads.
"""
from __future__ import annotations
import numpy as np


def read_window(page, y0: int, y1: int, x0: int, x1: int) -> np.ndarray:
    """Return ``page[y0:y1, x0:x1]`` decoding only the chunks it touches.

    Handles tiled pages (qptiff, LZW 512x512) and stripped pages (the rasterized
    mask, Deflate, RowsPerStrip=2) with the same code path: a strip is just a
    chunk whose width is the full image.
    """
    H, W = page.shape[-2], page.shape[-1]
    y0, x0 = max(0, y0), max(0, x0)
    y1, x1 = min(H, y1), min(W, x1)
    if y1 <= y0 or x1 <= x0:
        # Keep the caller's shape contract even when the window is fully clipped.
        return np.zeros((max(0, y1 - y0), max(0, x1 - x0)), dtype=page.dtype)

    if page.is_tiled:
        ch, cw = page.tilelength, page.tilewidth
    else:
        ch, cw = int(page.tags["RowsPerStrip"].value), W
    n_down, n_across = page.chunked[0], page.chunked[1]

    out = np.zeros((y1 - y0, x1 - x0), dtype=page.dtype)
    fh = None                                  # opened lazily: a fully sparse window reads nothing
    for ci in range(y0 // ch, min((y1 - 1) // ch + 1, n_down)):
        for cj in range(x0 // cw, min((x1 - 1) // cw + 1, n_across)):
            idx = ci * n_across + cj
            count = page.databytecounts[idx]
            if not count:                      # sparse chunk: leave as zeros
                continue
            if fh is None:
                fh = page.parent.filehandle
            fh.seek(page.dataoffsets[idx])
            chunk, _, shape = page.decode(fh.read(count), idx)
            chunk = np.asarray(chunk).reshape(shape[-3], shape[-2])[:ch, :cw]
            cy, cx = ci * ch, cj * cw
            sy0, sx0 = max(y0, cy), max(x0, cx)
            sy1, sx1 = min(y1, cy + ch), min(x1, cx + cw)
            out[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = chunk[sy0 - cy:sy1 - cy, sx0 - cx:sx1 - cx]
    return out


def selftest_bit_exact(page, size: int = 300, verbose: bool = True) -> bool:
    """Assert read_window == page.asarray()[slice] at five awkward positions."""
    H, W = page.shape[-2], page.shape[-1]
    ch = page.tilelength if page.is_tiled else int(page.tags["RowsPerStrip"].value)
    cw = page.tilewidth if page.is_tiled else W
    full = page.asarray()
    # Interior anchors must land inside the image. A stripped page has
    # chunk width == W, so "2 chunks in" is off the right edge; clamp to a
    # chunk boundary that exists.
    ay = min(2 * ch, max(0, ((H - size) // ch) * ch))
    ax = min(2 * cw, max(0, ((W - size) // cw) * cw))
    cases = {
        "chunk-aligned interior": (ay, ay + size, ax, ax + size),
        "unaligned interior":     (ay + 7, ay + 7 + size, ax + 13, ax + 13 + size),
        "top-left corner":        (0, size, 0, size),
        "bottom-right corner":    (H - size, H, W - size, W),
        "spans partial last chunk": (H - ch - 5, H, max(0, W - cw - 5), W),
    }
    # Fail closed: an empty comparison is not evidence of correctness, and a
    # self-test with no cases must never report success.
    if not cases:
        if verbose:
            print("  FAIL  no test cases could be placed in this page")
        return False
    ok = True
    for name, (a, b, c, d) in cases.items():
        got, ref = read_window(page, a, b, c, d), full[a:b, c:d]
        if ref.size == 0:
            same, why = False, "vacuous: reference slice is empty"
        else:
            same = got.shape == ref.shape and np.array_equal(got, ref)
            why = "" if same else "mismatch"
        ok &= same
        if verbose:
            print(f"  {'PASS' if same else 'FAIL'}  {name:<26} [{a}:{b}, {c}:{d}] {got.shape} {why}")
    return bool(ok)
