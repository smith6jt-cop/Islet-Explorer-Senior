"""Empty-slide nulls: the background a channel reads on glass, the floor of an edge or gate level (levels.level_for).

Ported from Islet-Explorer scripts/uniform/typing/redsea_cells.py at 578767e (``positivity.SOURCE_COMMIT``):
``EMPTY``, ``empty_windows``, ``window_means`` and the per-channel null record of ``process``. The method is unchanged:
means of 19 x 19-pixel windows (about a cell's measurement area) wholly in the slide area at least 20 um outside the
tissue, 20,000 seeded windows, their p50 / p99 / p99.9, on the despeckled channel (and p99 before the despeckle).

What the CODEX code assumed, and what is explicit here:
- One CODEX image is one section, and everything in the image outside that image's Tissue outline was slide area.
  An IES qptiff can hold several sections (a 67-plex whole slide holds four), and another section's tissue is not
  glass. The empty area is therefore computed once per slide, outside the union of the tissue outlines of EVERY
  section on that slide (``slide_nulls(outlines=...)``), and each section on the slide takes that slide's nulls.
- The CODEX images were stitched mosaics with no unscanned area. A qptiff keeps the slide's unscanned area as exact
  zeros (13 % of slide 0110618, 43 % of HDL011's qptiff at 1/16), which would enter as glass reading 0. The empty area
  is also restricted to scanned pixels, those above zero in ``scanned_channel`` (DAPI): ``empty_windows(scanned=)``.
Outlines are GeoJSON features in the qptiff's level-0 pixel frame, the frame Protein Measurement's geometry is in.
"""
from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

import numpy as np
import pandas as pd

from .coverage import despeckle

EMPTY = {'window_px': 19, 'windows': 20000, 'margin_um': 20.0, 'seed': 0}
NULL_COLUMNS = ('windows', 'null_p50', 'null_p99', 'null_p999', 'null_p99_before_despeckle')


class NullsError(ValueError):
    """A slide's empty area cannot be measured as declared."""


def empty_windows(outlines: Iterable[dict], H: int, W: int, pixel_um: float,
                  scanned: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Centers of EMPTY['window_px']-square windows wholly in the slide area at least EMPTY['margin_um'] outside every
    tissue outline in `outlines` (GeoJSON features, level-0 pixels), a seeded sample of EMPTY['windows']. `scanned`:
    the slide's scanned pixels (True); unscanned pixels are not slide area."""
    import cv2
    from shapely.geometry import shape

    mask = np.ones((H, W), np.uint8)
    if scanned is not None:
        if scanned.shape != (H, W):
            raise NullsError(f'scanned mask {scanned.shape} is not the slide {(H, W)}')
        mask &= scanned.astype(np.uint8)
    for feature in outlines:
        geom = shape(feature['geometry'])
        for poly in getattr(geom, 'geoms', [geom]):
            exterior = np.round(np.asarray(poly.exterior.coords)).astype(np.int32)
            x0, y0 = np.maximum(exterior.min(axis=0), 0)
            x1, y1 = np.minimum(exterior.max(axis=0) + 1, (W, H))
            if x1 <= x0 or y1 <= y0:
                continue
            # A hole belongs only to its own polygon. Painting holes into the
            # shared empty mask can erase a Tissue island already drawn there.
            local = np.zeros((y1 - y0, x1 - x0), np.uint8)
            origin = np.asarray((x0, y0), dtype=np.int32)
            cv2.fillPoly(local, [exterior - origin], 1)
            for hole in poly.interiors:
                ring = np.round(np.asarray(hole.coords)).astype(np.int32)
                cv2.fillPoly(local, [ring - origin], 0)
            mask[y0:y1, x0:x1] &= 1 - local
    r = int(round(EMPTY['margin_um'] / pixel_um))
    mask = cv2.erode(mask, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1)))
    w = EMPTY['window_px']
    ok = cv2.erode(mask, np.ones((w, w), np.uint8)).astype(bool)
    e = w // 2 + 1
    ok[:e], ok[-e:], ok[:, :e], ok[:, -e:] = False, False, False, False
    ys, xs = np.nonzero(ok)
    if len(ys):
        pick = np.random.default_rng(EMPTY['seed']).choice(len(ys), min(EMPTY['windows'], len(ys)), replace=False)
        ys, xs = ys[pick], xs[pick]
    return ys, xs


def window_means(a: np.ndarray, ys: np.ndarray, xs: np.ndarray) -> np.ndarray:
    import cv2

    ii = cv2.integral(np.ascontiguousarray(a, dtype=a.dtype.newbyteorder('=')), sdepth=cv2.CV_64F)
    h = EMPTY['window_px'] // 2
    y0, y1, x0, x1 = ys - h, ys + h + 1, xs - h, xs + h + 1
    return (ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0]) / EMPTY['window_px'] ** 2


def channel_null(raw: np.ndarray, ys: np.ndarray, xs: np.ndarray) -> dict:
    """One channel's null record (redsea_cells.process): on the despeckled channel, and p99 before the despeckle."""
    wm, wr = window_means(despeckle(raw), ys, xs), window_means(raw, ys, xs)
    return {'windows': len(wm), 'null_p50': float(np.percentile(wm, 50)), 'null_p99': float(np.percentile(wm, 99)),
            'null_p999': float(np.percentile(wm, 99.9)), 'null_p99_before_despeckle': float(np.percentile(wr, 99))}


def slide_nulls(read_channel: Callable[[str], np.ndarray], channels: Mapping[str, str], shape: tuple[int, int],
                outlines: Iterable[dict], pixel_um: float, scanned_channel: str = 'DAPI') -> dict[str, dict]:
    """{marker: null record} for one slide. `read_channel(name)` returns a level-0 channel as stored; `channels` maps
    each marker to its channel name; `outlines` are the tissue outlines of every section on the slide."""
    H, W = shape
    if scanned_channel not in set(channels.values()):
        raise NullsError(f'no {scanned_channel} channel to tell scanned slide from unscanned')
    scanned = read_channel(scanned_channel) > 0
    ys, xs = empty_windows(list(outlines), H, W, pixel_um, scanned)
    if not len(ys):
        raise NullsError('no empty window: the slide has no glass 20 um clear of every section')
    return {marker: channel_null(read_channel(name), ys, xs) for marker, name in channels.items()}


def read_outlines(path: Path) -> list[dict]:
    import json

    raw = json.loads(Path(path).read_text())
    return raw['features'] if isinstance(raw, dict) else raw


def cohort_nulls(manifest_path: Path, tissue_dir: Path, sections: list[str] | None = None) -> pd.DataFrame:
    """The empty-slide nulls of every section of a QuPath cohort manifest: (image, marker, windows, null_*).

    Sections are grouped by the qptiff the manifest names. A slide's empty area is outside the outlines of every
    manifest section on it (``tissue_dir/<section>.geojson``) and of any tissue on it that is no cohort section
    (``tissue_dir/_slide__<qptiff stem>.geojson``, optional). A section on the slide without an outline stops the run:
    its tissue would otherwise be read as glass. `sections` limits which sections' rows are returned; every section on
    their slides still bounds the glass."""
    import tifffile

    from ..artifacts import QuPathCohortManifest

    manifest = QuPathCohortManifest.read_json(manifest_path)
    slides: dict[str, list] = {}
    for image in manifest.images:
        slides.setdefault(str(Path(image.qptiff.path).resolve()), []).append(image)
    wanted = set(sections) if sections is not None else {i.donor_id for i in manifest.images}
    rows = []
    for qptiff, images in sorted(slides.items()):
        if not wanted & {i.donor_id for i in images}:
            continue
        missing = [i.donor_id for i in images if not (Path(tissue_dir) / f'{i.donor_id}.geojson').is_file()]
        if missing:
            raise NullsError(f'{qptiff}: no tissue outline for section(s) {missing} on this slide')
        outlines = [f for i in images for f in read_outlines(Path(tissue_dir) / f'{i.donor_id}.geojson')]
        extra = Path(tissue_dir) / f'_slide__{Path(qptiff).name.split(".")[0]}.geojson'
        if extra.is_file():
            outlines += read_outlines(extra)
        pixel = {float(i.pixel_size_um_x) for i in images}
        if len(pixel) != 1:
            raise NullsError(f'{qptiff}: sections disagree on the pixel size {sorted(pixel)}')
        channel_maps = {tuple((c.marker, c.channel_name, c.channel_index) for c in i.channel_map) for i in images}
        if len(channel_maps) != 1:
            raise NullsError(f'{qptiff}: sections on one slide record different channel maps')
        cmap = images[0].channel_map
        index = {c.channel_name: c.channel_index for c in cmap}
        with tifffile.TiffFile(qptiff) as tf:
            level = tf.series[0].levels[0]
            shape = level.shape[-2:]
            def read(name, pages=level.pages, index=index):
                return pages[index[name]].asarray()

            per = slide_nulls(read, {c.marker: c.channel_name for c in cmap}, shape, outlines, pixel.pop())
        for image in images:
            if image.donor_id in wanted:
                rows += [{'image': image.donor_id, 'marker': m, 'slide': qptiff, **v} for m, v in per.items()]
    return pd.DataFrame(rows, columns=['image', 'marker', 'slide', *NULL_COLUMNS])
