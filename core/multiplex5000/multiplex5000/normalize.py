"""Shared UniFORM background-mode shifts for protein measurements.

All samples for a marker enter one histogram basis. Each marker receives a
single log shift per sample; no positive fraction or disease label enters the
fit. Positivity cutoffs and tissue calls are downstream, per donor decisions.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import numpy as np

UNIFORM_COMMIT = "7c33396d12e5a19dc0b36e46a94c3293de0cd814"
NUM_BINS = 1024


class NormalizeError(ValueError):
    """A marker or pinned UniFORM checkout cannot be normalized as declared."""


@dataclass(frozen=True)
class ShiftResult:
    reference: str
    scale: dict[str, float]
    shift_log: dict[str, float]
    mode_sd_before: float
    mode_sd_after: float
    mode_fold_before: float
    mode_fold_after: float
    n_cells: dict[str, int]


def load_uniform(checkout: Path) -> tuple[ModuleType, ModuleType]:
    """Load the pinned separate UniFORM source without changing ``sys.path``."""
    checkout = Path(checkout).resolve()
    try:
        head = subprocess.run(["git", "-C", str(checkout), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True,
                              timeout=30).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as err:
        raise NormalizeError(f"cannot read UniFORM commit at {checkout}") from err
    if head != UNIFORM_COMMIT:
        raise NormalizeError(f"UniFORM is at {head}, expected {UNIFORM_COMMIT}")
    package = checkout / "UniFORM"
    location = package / "__init__.py"
    if not location.is_file():
        raise NormalizeError(f"UniFORM package missing at {package}")
    loaded = sys.modules.get("UniFORM")
    if loaded is not None:
        if Path(loaded.__file__).resolve() != location:
            raise NormalizeError("another UniFORM source is already imported")
    else:
        spec = importlib.util.spec_from_file_location(
            "UniFORM", location, submodule_search_locations=[str(package)])
        if spec is None or spec.loader is None:
            raise NormalizeError(f"cannot import UniFORM at {package}")
        loaded = importlib.util.module_from_spec(spec)
        sys.modules["UniFORM"] = loaded
        try:
            spec.loader.exec_module(loaded)
        except Exception:
            del sys.modules["UniFORM"]
            raise
    from UniFORM import normalization, registration

    return registration, normalization


def marker_histograms(values: dict[str, np.ndarray], num_bins: int = NUM_BINS):
    """UniFORM's log-floor histogram, with one global range for every sample."""
    if len(values) < 2:
        raise NormalizeError("a marker needs at least two measurable samples")
    logs = {}
    for sample, data in values.items():
        array = np.asarray(data, dtype=float)
        array = array[np.isfinite(array) & (array >= 1.0)]
        if not len(array):
            raise NormalizeError(f"{sample}: no finite marker values at or above log floor 1")
        logs[sample] = np.log(array)
    lo = min(float(v.min()) for v in logs.values())
    hi = max(float(v.max()) for v in logs.values())
    if lo >= hi:
        raise NormalizeError("marker has no intensity range to align")
    counts = {sample: np.histogram(log, bins=num_bins, range=(lo, hi))[0]
              for sample, log in logs.items()}
    return counts, {"global_min": lo, "global_max": hi,
                    "min_list": [float(v.min()) for v in logs.values()],
                    "max_list": [float(v.max()) for v in logs.values()]}


def background_mode_spread(counts: dict, ranges: dict, shifts: dict | None = None,
                           num_bins: int = NUM_BINS) -> dict:
    """Per-sample log modes, their SD and fold span, skipping the floor bin."""
    lo, hi = ranges["global_min"], ranges["global_max"]
    width = (hi - lo) / (num_bins - 1)
    peaks = {}
    for sample, count in counts.items():
        peak = lo + (int(np.argmax(count[1:])) + 1) * width
        if shifts is not None:
            peak -= shifts[sample] * width
        peaks[sample] = peak
    positions = np.asarray(list(peaks.values()))
    return {"peaks": peaks, "sd_log": float(positions.std()),
            "fold": float(np.exp(positions.max() - positions.min()))}


def fit_marker(marker: str, values: dict[str, np.ndarray], registration: ModuleType,
               normalization: ModuleType, *, reference: str | None = None,
               num_bins: int = NUM_BINS) -> ShiftResult:
    """Fit UniFORM's rigid shift for one marker on the same sample basis."""
    counts, ranges = marker_histograms(values, num_bins)
    samples = list(values)
    if reference is not None and reference not in samples:
        raise NormalizeError(f"{marker}: reference sample {reference} is unavailable")
    hist = {marker: {"counts": [counts[s] for s in samples],
                     "bin_edges": [None] * len(samples)}}
    shifts, references, _ = registration.automatic_registration(
        hist, [marker], [marker], samples, [reference] if reference else None, None, num_bins)
    bins = dict(zip(samples, shifts[marker], strict=True))
    logs = dict(zip(samples, normalization.calculate_shift_in_log_pixels(
        {marker: ranges}, [marker], num_bins, shifts)[marker], strict=True))
    before = background_mode_spread(counts, ranges, num_bins=num_bins)
    after = background_mode_spread(counts, ranges, bins, num_bins=num_bins)
    if not all(np.isfinite(value) for value in logs.values()):
        raise NormalizeError(f"{marker}: UniFORM returned a nonfinite shift")
    scales = {sample: float(np.exp(-logs[sample])) for sample in samples}
    if not all(np.isfinite(value) and value > 0 for value in scales.values()):
        raise NormalizeError(f"{marker}: UniFORM returned an invalid scale")
    return ShiftResult(references[0],
                       scales,
                       {sample: float(logs[sample]) for sample in samples},
                       before["sd_log"], after["sd_log"], before["fold"], after["fold"],
                       {sample: int(np.count_nonzero(np.isfinite(data) & (data >= 1.0)))
                        for sample, data in values.items()})
