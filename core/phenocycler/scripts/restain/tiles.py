"""Per-islet crops of the IF, on a grid that shares the image's own pixel lattice.

Two things here are easy to get wrong and neither raises.

**The reader clamps.**  ``tiff_window.read_window`` returns the INTERSECTION of the requested
window with the page, so a window that hangs off an edge comes back SMALLER than asked for --
not zero-padded.  Assigning that into a fixed-size array, or trusting ``arr.shape``, shifts
every pixel in the crop.  :func:`read_window_on_grid` therefore allocates the grid-shaped
buffer itself and places the returned block at its true offset.

**A sparse chunk decodes to black.**  A tile whose ``databytecounts`` entry is 0 is not
stored, and ``read_window`` leaves it as zeros.  Zero is a legal intensity here, so nothing
downstream can tell "no signal" from "never written".  The returned ``valid`` mask marks those
pixels, read from the tag alone with no decode, and it has to reach the FIGURES as well as the
statistics -- an unwritten tile painted at the display floor reads as measured background.

The grid origin is snapped so ``origin_um / pixel_um`` is an integer.  Label rasters and image
crops then land on the same lattice as the qptiff itself, so nothing is ever resampled and a
half-pixel shift cannot creep in between the geometry and the pixels.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "dct"))

from tiff_window import read_window                               # noqa: E402

from phenocycler.integration.rasterize import GridSpec            # noqa: E402


class TileError(RuntimeError):
    """A window cannot be read as the geometry it claims to cover."""


def window_grid(bounds_um, pixel_um: float) -> tuple[GridSpec, tuple[int, int, int, int]]:
    """Bounds ``(x0, y0, x1, y1)`` in microns -> (grid, ``(y0, y1, x0, x1)`` in pixels).

    Snapped OUTWARD to whole pixels of the image's own lattice, so the grid covers at least
    the requested bounds and ``origin_um / pixel_um`` is exactly an integer.
    """
    x0_um, y0_um, x1_um, y1_um = (float(v) for v in bounds_um)
    if not (x1_um > x0_um and y1_um > y0_um):
        raise TileError(f"empty bounds {bounds_um!r}")
    if not pixel_um > 0:
        raise TileError(f"pixel_um must be positive, got {pixel_um!r}")

    x0 = int(np.floor(x0_um / pixel_um))
    y0 = int(np.floor(y0_um / pixel_um))
    x1 = int(np.ceil(x1_um / pixel_um))
    y1 = int(np.ceil(y1_um / pixel_um))
    grid = GridSpec(origin_x_um=x0 * pixel_um, origin_y_um=y0 * pixel_um,
                    pixel_um=pixel_um, height=y1 - y0, width=x1 - x0)
    return grid, (y0, y1, x0, x1)


def _sparse_chunk_mask(page, y0: int, y1: int, x0: int, x1: int) -> np.ndarray:
    """True where a pixel falls in a chunk that was never written to the file.

    Read from ``databytecounts`` only -- no decode -- so it costs nothing and can be asked
    for on every window.
    """
    H, W = page.shape[-2], page.shape[-1]
    cy0, cx0 = max(0, y0), max(0, x0)
    cy1, cx1 = min(H, y1), min(W, x1)
    out = np.zeros((max(0, cy1 - cy0), max(0, cx1 - cx0)), dtype=bool)
    if out.size == 0:
        return out
    if page.is_tiled:
        ch, cw = page.tilelength, page.tilewidth
    else:
        ch, cw = int(page.tags["RowsPerStrip"].value), W
    n_down, n_across = page.chunked[0], page.chunked[1]
    for ci in range(cy0 // ch, min((cy1 - 1) // ch + 1, n_down)):
        for cj in range(cx0 // cw, min((cx1 - 1) // cw + 1, n_across)):
            if page.databytecounts[ci * n_across + cj]:
                continue
            sy0, sx0 = max(cy0, ci * ch), max(cx0, cj * cw)
            sy1, sx1 = min(cy1, ci * ch + ch), min(cx1, cj * cw + cw)
            out[sy0 - cy0:sy1 - cy0, sx0 - cx0:sx1 - cx0] = True
    return out


def read_window_on_grid(page, grid: GridSpec, pixel_um: float) -> tuple[np.ndarray, np.ndarray]:
    """``(image, valid)``, both exactly ``grid.shape``.

    ``valid`` is False where the grid falls outside the page and where it lands on a chunk the
    file does not store.  The image itself is 0 there, which is why the mask exists.
    """
    if abs(grid.pixel_um - pixel_um) > 1e-12:
        raise TileError(f"grid.pixel_um={grid.pixel_um} but the read was asked for {pixel_um}")
    y0 = int(round(grid.origin_y_um / pixel_um))
    x0 = int(round(grid.origin_x_um / pixel_um))
    y1, x1 = y0 + grid.height, x0 + grid.width

    H, W = page.shape[-2], page.shape[-1]
    cy0, cx0 = max(0, y0), max(0, x0)
    cy1, cx1 = min(H, y1), min(W, x1)

    image = np.zeros(grid.shape, dtype=page.dtype)
    valid = np.zeros(grid.shape, dtype=bool)
    if cy1 <= cy0 or cx1 <= cx0:
        return image, valid                       # entirely off the page: nothing is available

    block = read_window(page, cy0, cy1, cx0, cx1)
    ry0, rx0 = cy0 - y0, cx0 - x0
    if block.shape != (cy1 - cy0, cx1 - cx0):
        raise TileError(
            f"read_window returned {block.shape} for a {(cy1 - cy0, cx1 - cx0)} window; "
            "the reader's clamping contract changed")
    image[ry0:ry0 + block.shape[0], rx0:rx0 + block.shape[1]] = block
    valid[ry0:ry0 + block.shape[0], rx0:rx0 + block.shape[1]] = True
    valid[ry0:ry0 + block.shape[0], rx0:rx0 + block.shape[1]] &= \
        ~_sparse_chunk_mask(page, cy0, cy1, cx0, cx1)
    return image, valid


def level_for(levels_shapes, base_um_per_px: float, target_um: float) -> int:
    """Index of the pyramid level closest to ``target_um``, never finer than it needs.

    Separate from ``rasterize._pick_pyramid_level`` only in that the base scale is REQUIRED:
    this lane always knows its recorded pixel size and must never fall back to a tag.
    """
    h0 = levels_shapes[0][0]
    best, best_err = 0, np.inf
    for i, (h, _w) in enumerate(levels_shapes):
        if h <= 0:
            continue
        um = base_um_per_px * (h0 / h)
        if um > target_um * 1.6:
            continue
        err = abs(um - target_um)
        if err < best_err:
            best, best_err = i, err
    return best
