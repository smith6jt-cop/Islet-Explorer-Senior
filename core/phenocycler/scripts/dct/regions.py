"""Region selection: 8 anchor points per donor, drawn by a rule that cannot see the calls.

Everything executed before 2026-08-24 rested on one donor and one region, and it was the
wrong donor -- peri-islet cells are 0.03 % of donor 6539 against a cohort median of 1.95 %,
a ~65x outlier at the bottom of the one quantity this study measures.  This module fixes
the sampling frame in advance so no later draw can be steered by a result.

Four regions per donor come from a hash over the donor id and the origin, which is
deterministic, reproducible from the donor id alone and blind to content.  Four more come
from the same hash order restricted to cores that contain an islet, because the study
question is the peri-islet immune environment and a blind draw over whole-pancreas area
misses islets entirely on the low-peri donors.

**Regions are anchor points, not tiles.**  A window-size arm reads concentric windows around
the same anchor and scores only the common core, so it changes the model's normalization
context and not the cell population -- which is what the superseded `--median-tile` arm got
wrong, having changed the tile and the tissue at once.

Selection needs no mask and no image: centroids and QuPath's own islet parent are already in
the run's expression parquet.

**Trap, verified 2026-08-24 against the archived mask:** ``X_centroid``/``Y_centroid`` in that
parquet are **microns**, not pixels.  Divided by ``pixel_size_um_x`` from the manifest, 98.2 %
of a window's cells land inside it; taken as pixels, 0 % do.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

__all__ = ["Region", "select_regions", "region_hash", "CORE_PX", "MIN_CELLS"]

CORE_PX = 1024      # the scored core, common to every window size
MIN_CELLS = 500     # below this an anchor cannot power the gate's anchors
MAX_WINDOW = 4096   # the largest concentric window an arm may read
N_BLIND = 4
N_ISLET = 4


@dataclass(frozen=True)
class Region:
    donor_id: str
    y0: int              # level-0 pixel origin of the scored core
    x0: int
    core_px: int
    n_cells: int
    n_islet_core_cells: int
    draw: str            # "blind" or "islet"
    order: int           # position in the hash order it was drawn from
    digest: str


def region_hash(donor_id: str, y0: int, x0: int) -> str:
    """Deterministic, content-blind, and reproducible from the donor id alone."""
    return hashlib.sha256(f"{donor_id}:{y0}:{x0}".encode()).hexdigest()


def _centroids_px(expression_parquet: pathlib.Path, pixel_size_um: float):
    df = pd.read_parquet(expression_parquet, columns=["X_centroid", "Y_centroid", "cell_region"])
    if not np.isfinite(pixel_size_um) or pixel_size_um <= 0:
        raise ValueError(f"pixel size {pixel_size_um!r} is not usable")
    x = df["X_centroid"].to_numpy(float) / pixel_size_um
    y = df["Y_centroid"].to_numpy(float) / pixel_size_um
    islet = df["cell_region"].astype(str).to_numpy() == "core"
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok], islet[ok]


def select_regions(donor_id: str, expression_parquet, pixel_size_um: float,
                   image_shape, *, core_px: int = CORE_PX, min_cells: int = MIN_CELLS,
                   n_blind: int = N_BLIND, n_islet: int = N_ISLET,
                   max_window: int = MAX_WINDOW) -> list:
    """The eight anchors for one donor.  Raises rather than returning a short list."""
    H, W = image_shape
    x, y, islet = _centroids_px(pathlib.Path(expression_parquet), pixel_size_um)

    # A core must sit far enough inside the slide for the largest concentric window to fit,
    # so a window-size arm never silently reads a clipped field.
    margin = (max_window - core_px) // 2
    gy, gx = (y // core_px).astype(int), (x // core_px).astype(int)
    counts, islets = {}, {}
    for j, i, is_islet in zip(gy, gx, islet):
        counts[(j, i)] = counts.get((j, i), 0) + 1
        islets[(j, i)] = islets.get((j, i), 0) + int(is_islet)

    candidates = []
    for (j, i), n in counts.items():
        # plain ints: numpy scalars leak through into the JSON artifact otherwise
        y0, x0 = int(j) * core_px, int(i) * core_px
        if n < min_cells:
            continue
        if y0 < margin or x0 < margin or y0 + core_px + margin > H or x0 + core_px + margin > W:
            continue
        candidates.append((region_hash(donor_id, y0, x0), y0, x0, int(n), int(islets[(j, i)])))
    candidates.sort()

    chosen, taken = [], set()
    for rank, (dg, y0, x0, n, ni) in enumerate(candidates):
        if len(chosen) >= n_blind:
            break
        chosen.append(Region(donor_id, y0, x0, core_px, n, ni, "blind", rank, dg))
        taken.add((y0, x0))
    n_islet_taken = 0
    for rank, (dg, y0, x0, n, ni) in enumerate(candidates):
        if n_islet_taken >= n_islet:
            break
        if ni < 1 or (y0, x0) in taken:
            continue
        chosen.append(Region(donor_id, y0, x0, core_px, n, ni, "islet", rank, dg))
        taken.add((y0, x0))
        n_islet_taken += 1

    if len(chosen) < n_blind + n_islet:
        raise ValueError(
            f"donor {donor_id}: only {len(chosen)} of {n_blind + n_islet} regions available "
            f"({len(candidates)} cores with >= {min_cells} cells and room for a {max_window}-px "
            f"window; {sum(1 for c in candidates if c[4] >= 1)} of them contain an islet). "
            "A short draw is a finding about the donor, not a smaller sample to proceed with."
        )
    return chosen


def build_all(manifest_path, expression_root, donors, out_path=None) -> dict:
    """Write every donor's origins as one tracked artifact, before any inference runs."""
    import tifffile

    man = json.loads(pathlib.Path(manifest_path).read_text())
    by_donor = {im["donor_id"]: im for im in man["images"]}
    out = {"schema_version": 1, "method_version": "dct-regions-2026-08-24",
           "core_px": CORE_PX, "min_cells": MIN_CELLS, "max_window": MAX_WINDOW,
           "rule": ("sha256(donor:y0:x0) ascending over 1024-px cores with >= 500 cells and room "
                    "for a 4096-px concentric window; first 4, then the first 4 of the same order "
                    "that contain a cell whose QuPath parent is Islet_N. Content-blind and "
                    "reproducible from the donor id alone."),
           "centroid_units": "X_centroid/Y_centroid in the expression parquet are MICRONS; "
                             "divide by pixel_size_um_x for level-0 pixels (verified against the "
                             "archived mask, 2026-08-24)",
           "regions": {}}
    for d in donors:
        im = by_donor[d]
        with tifffile.TiffFile(im["qptiff"]["path"]) as tf:
            page = tf.series[0].pages[0]
            shape = (int(page.imagelength), int(page.imagewidth))
        regs = select_regions(d, pathlib.Path(expression_root) / f"donor_id={d}" / "data.parquet",
                              im["pixel_size_um_x"], shape)
        out["regions"][d] = [asdict(r) for r in regs]
    if out_path:
        pathlib.Path(out_path).write_text(json.dumps(out, indent=2) + "\n")
    return out
