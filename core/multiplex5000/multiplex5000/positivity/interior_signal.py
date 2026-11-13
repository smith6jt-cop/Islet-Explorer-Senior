"""Exclusive native-pixel evidence for a two-way protein identity conflict."""
from __future__ import annotations

import json
import numpy as np

from . import coverage as AC


def compare_signals(left, right):
    """Return exclusive pixel counts and median brightness for two calibrated groups.

    Input matrices are markers by interior pixels, on the uncapped scale
    (native signal - marker level) / marker span. A group is detectable at a
    pixel if any of its active, group-specific markers clears its level.
    """
    left, right = np.asarray(left, float), np.asarray(right, float)
    if (left.ndim != 2 or right.ndim != 2 or left.shape[1] != right.shape[1]
            or not left.size or not right.size or not np.isfinite(left).all()
            or not np.isfinite(right).all()):
        return None
    a, b = left.mean(axis=0), right.mean(axis=0)
    a_on, b_on = (left >= 0).any(axis=0), (right >= 0).any(axis=0)
    a_wins = a_on & (~b_on | (a > b))
    b_wins = b_on & (~a_on | (b > a))
    return a.size, int(a_wins.sum()), int(b_wins.sum()), float(np.median(a)), float(np.median(b))


class InteriorSignal:
    """Cache only the requested cell/channel samples from one section's native image."""

    def __init__(self, open_window, object_ids, counts, reg, markers, levels, tops):
        self.open_window = open_window
        self.object_ids = np.asarray(object_ids)
        self.counts = np.asarray(counts)
        self.markers = markers
        self.reg = reg
        self.levels = np.asarray(levels)
        self.spans = np.asarray(tops) - self.levels
        self.window = None
        self.positions = {}
        self.samples = {}
        self.evidence = {}

    def _open(self):
        features, shape, read, channels = self.open_window()
        if not features:
            raise ValueError("no native cell geometries for spatial conflict resolution")
        mask, nucleus, ids, _ = AC.rasterize(features, shape)
        del nucleus
        self.labels = {oid: label for label, oid in enumerate(ids, 1)}
        if len(self.labels) != len(ids) or None in self.labels:
            raise ValueError("native cell geometry IDs must be unique and non-null")
        self.window = read, channels
        positions = AC.interior_positions(mask)
        self.owners = np.zeros(mask.size, dtype=mask.dtype)
        self.owners[positions] = mask.ravel()[positions]

    def _load(self, needed):
        if self.window is None:
            self._open()
        missing = {row for rows in needed.values() for row in rows if row not in self.positions}
        if missing:
            labels = {self.labels[self.object_ids[row]]: row for row in missing}
            selected = np.flatnonzero(np.isin(self.owners, list(labels), kind="table"))
            selected = selected[np.argsort(self.owners[selected], kind="stable")]
            owners = self.owners[selected]
            for label, row in labels.items():
                start, end = np.searchsorted(owners, (label, label + 1))
                if end - start != self.counts[row]:
                    raise ValueError(f"{self.object_ids[row]}: native interior pixel count differs from coverage")
                self.positions[row] = selected[start:end]
        read, channels = self.window
        for col, rows in needed.items():
            marker = self.markers[col]
            name = self.reg.table.loc[marker, "channel"]
            plane = AC.despeckle(read(name if name in channels else marker)).ravel()
            for row in rows:
                raw = plane[self.positions[row]].astype(float)
                self.samples[row, col] = (raw - self.levels[col]) / self.spans[col]

    def measure(self, requests, stage):
        active, needed = [], {}
        for row, left, right in requests:
            shared = set(left) & set(right)
            left = tuple(m for m in left if m not in shared)
            right = tuple(m for m in right if m not in shared)
            cols = list(left + right)
            if not left or not right or self.counts[row] == 0:
                continue
            if not (np.isfinite(self.spans[cols]) & (self.spans[cols] > 0) & np.isfinite(self.levels[cols])).all():
                continue
            active.append((row, left, right))
            for col in cols:
                if (row, col) not in self.samples:
                    needed.setdefault(col, set()).add(row)
        if needed:
            self._load(needed)
        result = {}
        for row, left, right in active:
            measurement = compare_signals([self.samples[row, m] for m in left], [self.samples[row, m] for m in right])
            if measurement is None:
                continue
            result[row] = measurement
            n, a, b, ma, mb = measurement
            self.evidence.setdefault(row, []).append({
                "stage": stage, "markers": [[self.markers[m] for m in left], [self.markers[m] for m in right]],
                "interior_pixels": n, "winner_pixels": [a, b], "median_scaled": [ma, mb]})
        return result

    def audit(self, n):
        return [json.dumps(self.evidence[i], separators=(",", ":")) if i in self.evidence else "" for i in range(n)]
