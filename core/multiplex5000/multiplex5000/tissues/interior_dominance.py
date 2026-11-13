"""Resolve positive, tied protein identities from the native segmentation and image.

Area is a difference in fractions of ALL owned interior pixels. Co-positive
pixels have one vote at most. Intensity calibration is marker-specific and
never clipped at the cell-table display ceiling.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_erosion, median_filter

METHOD = "native-interior-area-brightness-v1"
TYPED = "interior_area_brightness"


def _record(status, **fields):
    return dict(status=status, winner=None, interior_pixels=0, votes={},
                medians={}, area_lead=None) | fields


def score_interior(signals: dict, gates: dict) -> dict:
    """Score calibrated per-marker vectors from the SAME interior, in native order."""
    markers = list(signals)
    if len(markers) < 2:
        raise ValueError("a conflict needs at least two distinct contenders")
    vectors = [np.asarray(signals[m], dtype=np.float64) for m in markers]
    if any(v.ndim != 1 or v.shape != vectors[0].shape for v in vectors):
        raise ValueError("signals must be matching one-dimensional interior vectors")
    x = np.stack(vectors)
    n = x.shape[1]
    if not n:
        return _record("empty_interior")
    threshold = np.asarray([gates.get(m, np.nan) for m in markers], dtype=float)
    if not np.isfinite(x).all() or not np.isfinite(threshold).all():
        return _record("unavailable_signal_or_gate", interior_pixels=n)
    eligible = x >= threshold[:, None]
    highest = np.where(eligible, x, -np.inf).max(axis=0)
    top = eligible & (x == highest)
    counts = (top & (top.sum(axis=0) == 1)).sum(axis=1)
    med = np.median(x, axis=1)
    idx = int(counts.argmax())
    lead = int(counts[idx] - np.max(np.delete(counts, idx)))
    evidence = _record("insufficient_area_lead", interior_pixels=n,
                       votes=dict(zip(markers, map(int, counts))),
                       medians=dict(zip(markers, map(float, med))), area_lead=lead / n)
    # An exact 15 percentage points qualifies. Never round a fraction to test the cutoff.
    if 100 * lead >= 15 * n:
        if med[idx] > np.max(np.delete(med, idx)):
            evidence.update(status=TYPED, winner=markers[idx])
        else:
            evidence.update(status="brightness_does_not_confirm_area")
    return evidence


def _bounds(geom):
    if geom["type"] not in {"Polygon", "MultiPolygon"}:
        raise ValueError("expected QuPath cell polygons")
    polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
    pts = [p for poly in polys for p in poly[0]]
    xs, ys = zip(*[(p[0], p[1]) for p in pts])
    return min(xs), min(ys), max(xs), max(ys)


def _bins(bounds):
    x0, y0, x1, y1 = bounds
    for y in range(int(y0 // 128), int(y1 // 128) + 1):
        for x in range(int(x0 // 128), int(x1 // 128) + 1):
            yield x, y


def _target_masks(path, ids, height, width):
    """Replay all intersecting polygons in export order, including holes and overwrites."""
    from ..artifacts import iter_qupath_geojson_features
    from ..redsea import _draw_geometry

    windows, index = {}, defaultdict(list)
    for f in iter_qupath_geojson_features(path):
        oid = str(f.get("id"))
        if oid not in ids:
            continue
        if oid in windows:
            raise ValueError(f"duplicate target object id: {oid}")
        bx, by, rx, ry = _bounds(f["geometry"])
        x0, y0 = max(0, int(np.floor(bx))-4), max(0, int(np.floor(by))-4)
        x1, y1 = min(width, int(np.ceil(rx))+5), min(height, int(np.ceil(ry))+5)
        if x0 >= x1 or y0 >= y1:
            raise ValueError(f"target {oid} lies outside its contracted image")
        windows[oid] = (x0, y0, np.zeros((y1-y0, x1-x0), dtype=bool))
        for key in _bins((x0, y0, x1-1, y1-1)):
            index[key].append(oid)
    if windows:
        for f in iter_qupath_geojson_features(path):
            geom = f["geometry"]
            bx, by, rx, ry = bounds = _bounds(geom)
            for oid in {oid for key in _bins(bounds) for oid in index.get(key, ())}:
                x, y, mask = windows[oid]
                h, w = mask.shape
                if bx <= x+w-1 and by <= y+h-1 and rx >= x and ry >= y:
                    _draw_geometry(mask, geom, str(f.get("id")) == oid,
                                   mask.shape, 1.0, origin=(x, y))
    return windows


def measure_conflicts(image, contenders: dict, calibration: dict, gate) -> dict:
    """Measure exactly these Object IDs in a QuPathImageManifest at native level 0.

    calibration maps marker -> (raw alignment multiplier, aligned floor, aligned span);
    gate is either one aligned threshold or a marker-specific threshold mapping.
    Missing evidence is returned explicitly; channel/geometry contract errors raise.
    """
    from ..redsea import qptiff_channels, validated_channel_markers
    from ..tiff_window import read_window

    results, ready = {}, {}
    for oid, markers in contenders.items():
        if len(markers) < 2 or len(set(markers)) != len(markers):
            raise ValueError(f"{oid}: contenders must be two or more distinct markers")
        gates = gate if isinstance(gate, dict) else dict.fromkeys(markers, gate)
        def usable(m):
            c = calibration.get(m)
            return (c is not None and len(c) == 3 and np.isfinite(c).all() and c[0] > 0
                    and c[2] > 0 and gates.get(m) is not None and np.isfinite(gates[m]))
        if all(usable(m) for m in markers):
            ready[oid] = markers
        else:
            results[oid] = _record("unavailable_calibration_or_gate")
    if not ready:
        return results
    for file in (image.cell_geojson, image.qptiff):
        path = Path(file.path)
        if not path.is_file():
            return results | {oid: _record("missing_native_source") for oid in ready}
        stat = path.stat()
        if stat.st_size != file.size_bytes:
            raise ValueError(f"native source differs from its manifest: {path}")
        if stat.st_mtime_ns != file.mtime_ns:
            from ..artifacts import FileFingerprint
            current = FileFingerprint.capture(path, hash_mode=file.hash_mode,
                                              sample_size_bytes=file.sample_size_bytes,
                                              sample_count=file.sample_count)
            if current.sha256 != file.sha256:
                raise ValueError(f"native source digest differs from its manifest: {path}")
    tif, pages, names, (height, width), _ = qptiff_channels(Path(image.qptiff.path), 0)
    try:
        channels = validated_channel_markers(image, names)
        for oid, markers in list(ready.items()):
            if any(channels.count(m) != 1 for m in markers):
                results[oid] = _record("unavailable_channel")
                del ready[oid]
        windows = _target_masks(image.cell_geojson.path, ready, height, width) if ready else {}
        for oid, markers in ready.items():
            if oid not in windows:
                results[oid] = _record("object_missing_in_segmentation")
                continue
            x, y, mask = windows[oid]
            h, w = mask.shape
            interior = binary_erosion(mask, structure=np.ones((5, 5), dtype=bool))
            vectors = {}
            for m in markers:
                raw = read_window(pages[channels.index(m)], y, y+h, x, x+w)
                smooth = median_filter(raw, size=3, mode="nearest")
                factor, floor, span = calibration[m]
                vectors[m] = (smooth[interior].astype(np.float64) * factor - floor) / span
            gates = gate if isinstance(gate, dict) else dict.fromkeys(markers, gate)
            results[oid] = score_interior(vectors, gates)
    finally:
        tif.close()
    return results
