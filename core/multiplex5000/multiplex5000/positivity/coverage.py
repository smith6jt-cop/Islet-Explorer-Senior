"""Area coverage: per cell and marker, the proportion of the cell's pixels at or above the marker's level.

Ported from Islet-Explorer scripts/uniform/typing/area_coverage.py at 578767e (``positivity.SOURCE_COMMIT``), with the
raster and despeckle of redsea_cells.py it uses. The computation is unchanged. Changed in the move: the image and the
geometry are passed in (in IES, the qptiff and the cell GeoJSON the QuPath cohort manifest names, already in the
qptiff's level-0 frame) instead of read from a CODEX run, and the run-directory fingerprints are left to the stage's
receipt. The despeckle is the same 3x3 median with replicated borders, computed with scipy instead of OpenCV.

Owner 2026-10-29: "There are two markers for a T Cell CD8a and both are brighter and cover more cell area than any
myeloid marker so this should not be a conflicting call." These marginal fractions describe each marker separately.
Identity conflicts use exclusive native-interior pixel votes in interior_signal.py.

Pixels: the despeckled channel, uncompensated. Level: the image's level of the marker (levels level_raw, on the
marker's measurement). Measured over the measurement's pixels: the whole cell, or REDSEA's 1-pixel inner boundary band
for a marker the level selection reads on the membrane.

Interior coverage (covint__<marker>, owner 2026-10-29 "Interior coverage, clear lead"): the same proportion over the
cell's interior, the pixels at least INTERIOR_PX from its outline, for every marker whatever its measurement. A
neighbor's signal left at the contact rim does not reach the interior; on 127 tie-settled cells checked against the
image (output/typing_accuracy_20261029/), the difference in interior coverage separated "the placed group fills the
cell" from "both or neither fill it" at AUC 0.87, whole-cell coverage at 0.76.

``compute`` returns (object_id, n_interior_pixels, cov__<marker>, covint__<marker>) for every typing marker with a
finite level.
"""
from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np
import pandas as pd
from scipy.ndimage import maximum_filter, median_filter, minimum_filter
from skimage.draw import polygon as draw_polygon
from skimage.segmentation import find_boundaries

from . import registry as RG

PREFIX = 'cov__'
INTERIOR_PREFIX = 'covint__'
INTERIOR_PX = 2
INTERIOR_COUNT = 'n_interior_pixels'
SCHEMA_VERSION = 1
DESPECKLE = {'filter': 'median', 'size': 3, 'borders': 'replicated', 'channels': 'all'}


def despeckle(a: np.ndarray) -> np.ndarray:
    """ImageJ's Despeckle: the 3x3 median of each pixel (borders replicated), on the channel's own integer values.

    OpenCV's medianBlur where it takes the type (uint8, uint16, float32 at this size; the source's own call), SciPy's
    median filter otherwise: equal outputs (tests), and on a 1.15 Gpx slide channel about 250x apart in time (2 min
    against well under a second)."""
    native = np.ascontiguousarray(a, dtype=a.dtype.newbyteorder('='))       # some TIFFs are big-endian
    if native.ndim == 2 and native.dtype in (np.uint8, np.uint16, np.float32):
        try:
            import cv2
        except ImportError:
            pass
        else:
            return cv2.medianBlur(native, DESPECKLE['size'])
    return median_filter(native, size=DESPECKLE['size'], mode='nearest')


# --------------------------------------------------------------------------- raster (redsea_cells.py, unchanged)
def draw_geometry(target, geom, label, shape):
    polys = geom['coordinates'] if geom['type'] == 'MultiPolygon' else [geom['coordinates']]
    for poly in polys:
        ext = np.asarray(poly[0], dtype=np.float64)
        if ext.shape[0] < 3:
            continue
        pr, pc = draw_polygon(ext[:, 1], ext[:, 0], shape=shape)
        target[pr, pc] = label
        for hole in poly[1:]:
            h = np.asarray(hole, dtype=np.float64)
            if h.shape[0] < 3:
                continue
            hr, hc = draw_polygon(h[:, 1], h[:, 0], shape=shape)
            target[hr, hc] = 0


def rasterize(features, shape):
    mask = np.zeros(shape, dtype=np.int32)
    nuc = np.zeros(shape, dtype=np.int32)
    oids, has_nucleus = [], []
    for i, ft in enumerate(features):
        oids.append(ft['id'])
        draw_geometry(mask, ft['geometry'], i + 1, shape)
        ng = ft.get('nucleusGeometry')
        present = bool(ng and ng.get('coordinates'))
        has_nucleus.append(present)
        if present:
            draw_geometry(nuc, ng, i + 1, shape)
    np.putmask(nuc, nuc != mask, 0)       # nucleus a strict subset of its own cell's final pixels
    return mask, nuc, oids, np.asarray(has_nucleus)


# --------------------------------------------------------------------------- coverage
def validate_frame(frame: pd.DataFrame, object_ids: pd.Series, markers: list[str], tag: str) -> None:
    required = ['object_id', INTERIOR_COUNT] + [INTERIOR_PREFIX + m for m in markers]
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f'{tag}: coverage missing columns: {missing}; recompute coverage')
    ids = frame['object_id']
    if ids.isna().any() or ids.duplicated().any() or object_ids.isna().any() or object_ids.duplicated().any():
        raise ValueError(f'{tag}: coverage/source object IDs must be non-null and unique')
    missing_ids, extra_ids = set(object_ids) - set(ids), set(ids) - set(object_ids)
    if missing_ids or extra_ids:
        raise ValueError(f'{tag}: coverage object IDs differ from REDSEA: {len(missing_ids)} missing, '
                         f'{len(extra_ids)} extra')
    count = frame[INTERIOR_COUNT].to_numpy(float)
    if not np.isfinite(count).all() or (count < 0).any() or (count != np.floor(count)).any():
        raise ValueError(f'{tag}: invalid interior pixel count')
    values = frame[[INTERIOR_PREFIX + m for m in markers]].to_numpy(float)
    has = count > 0
    if not np.isfinite(values[has]).all() or ((values[has] < 0) | (values[has] > 1)).any():
        raise ValueError(f'{tag}: invalid interior coverage for cells with interior pixels')
    if not np.isnan(values[~has]).all():
        raise ValueError(f'{tag}: cells with no interior must have NaN interior coverage')


def read_validated(frame: pd.DataFrame, all_ids: pd.Series, tag: str, object_ids: pd.Series,
                   markers: list[str]) -> pd.DataFrame:
    """The coverage of `object_ids` (identity's cells), in their order, after the source frame is checked against every
    cell the section's REDSEA table holds (`all_ids`). area_coverage.read_validated without the run fingerprints."""
    validate_frame(frame, all_ids, markers, tag)
    if object_ids.isna().any() or object_ids.duplicated().any() or not set(object_ids) <= set(all_ids):
        raise ValueError(f'{tag}: identity input object IDs must be unique and present in coverage')
    return frame.set_index('object_id').loc[object_ids].reset_index()


def interior_positions(mask):
    """Flat positions owned by the same cell throughout the interior erosion window."""
    width = 2 * INTERIOR_PX + 1
    return np.flatnonzero(((minimum_filter(mask, size=width) == mask) &
                           (maximum_filter(mask, size=width) == mask) & (mask > 0)).ravel())


def compute(tag: str, features: Iterable[dict], shape: tuple[int, int], read_channel: Callable[[str], np.ndarray],
            channel_names: list[str], levels: pd.DataFrame, reg: RG.Registry) -> pd.DataFrame:
    """area_coverage.one(): `features` the section's cell GeoJSON features, `shape` the image (H, W), `read_channel`
    returns one channel as stored (its own integer type) by name, `channel_names` the image's channels, `levels` the
    levels table (image, marker, measurement, level_raw)."""
    lv = levels[(levels['image'] == tag) & np.isfinite(levels['level_raw'])].set_index('marker')
    H, W = shape
    mask, _, oids, _ = rasterize(features, (H, W))
    n = len(oids)
    flat = mask.ravel()
    size = np.bincount(flat, minlength=n + 1)[1:].astype(float)
    band = np.flatnonzero(find_boundaries(mask, mode='inner', connectivity=2, background=0).ravel())
    interior = interior_positions(mask)
    del mask
    iflat = flat[interior]
    isize = np.bincount(iflat, minlength=n + 1)[1:].astype(float)
    bflat = flat[band]
    bsize = np.bincount(bflat, minlength=n + 1)[1:].astype(float)
    frame = {'object_id': oids, INTERIOR_COUNT: isize.astype(np.int64)}
    for m in reg.typing_markers():
        if m not in lv.index:
            continue
        ch = reg.table.loc[m, 'channel']
        ch = ch if ch in channel_names else m       # the registry's PGP9_5 is PGP9.5 in the images
        if ch not in channel_names:
            raise ValueError(f'{tag}: no image channel for required coverage marker {m}')
        above = despeckle(read_channel(ch)).ravel() >= lv.loc[m, 'level_raw']
        with np.errstate(invalid='ignore', divide='ignore'):
            if lv.loc[m, 'measurement'] == 'membrane':
                frame[PREFIX + m] = np.bincount(bflat, weights=above[band], minlength=n + 1)[1:] / \
                    np.where(bsize > 0, bsize, np.nan)
            else:
                frame[PREFIX + m] = np.bincount(flat, weights=above, minlength=n + 1)[1:] / \
                    np.where(size > 0, size, np.nan)
            frame[INTERIOR_PREFIX + m] = np.bincount(iflat, weights=above[interior], minlength=n + 1)[1:] / \
                np.where(isize > 0, isize, np.nan)
        del above
    out = pd.DataFrame(frame)
    validate_frame(out, out['object_id'], [c.removeprefix(INTERIOR_PREFIX) for c in out
                                           if c.startswith(INTERIOR_PREFIX)], tag)
    return out
