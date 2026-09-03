"""Display windows for image panels, with the black point nailed to the measured floor.

The rule this implements is stated in ``scripts/celltype/figures.py:3-13`` and has, until now,
had no implementation anywhere in the repo -- it was prose, re-applied by hand:

    window on the workflow's own measured floor, because a percentile black point hides the
    zero-inflation, which is itself a finding; no gamma < 1.

So ``DisplayRange.vmin`` is a hardcoded 0.0 and there is **no parameter that can change it**.
That is deliberate: a black-point argument with a default of 0 is one keyword away from being
set, and the one non-compliant display path in this repo
(``scripts/review/boundary_contact_audit.py:789`` -- a 1st-percentile black point, then a clip)
got there by exactly that route. ``integration/rasterize.py:407 normalize_image`` is likewise
a percentile stretch; it is registration input, never a figure, and is not imported here.

The zeros are also COUNTED. Pinning the black point is only worth anything if the thing it
protects is reported, so every range carries ``frac_at_zero`` and every panel prints it.

A percentile WHITE point is allowed, because a handful of hot pixels otherwise compresses the
whole ramp -- but it must be stated on the panel, and the fraction it burns must be counted.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class DisplayError(RuntimeError):
    """A display window cannot be built as specified."""


@dataclass(frozen=True)
class DisplayRange:
    """An unsigned intensity window. ``vmin`` is 0.0 and cannot be anything else."""

    vmax: float
    high_percentile: float
    ceiling: float
    n_finite: int
    n_at_zero: int
    n_at_or_above_vmax: int
    n_at_ceiling: int
    n_nonfinite: int
    fallback_reason: str | None = None
    vmin: float = 0.0

    def __post_init__(self):
        if self.vmin != 0.0:
            raise DisplayError(
                "the black point is the measured floor, 0. A percentile black point hides the "
                "zero-inflation, which is the thing these panels exist to show.")
        if not (self.vmax > 0.0):
            raise DisplayError(f"vmax must be positive, got {self.vmax!r}")

    @property
    def frac_at_zero(self) -> float:
        return self.n_at_zero / self.n_finite if self.n_finite else 0.0

    @property
    def frac_at_or_above_vmax(self) -> float:
        return self.n_at_or_above_vmax / self.n_finite if self.n_finite else 0.0

    @property
    def frac_at_ceiling(self) -> float:
        return self.n_at_ceiling / self.n_finite if self.n_finite else 0.0

    def caption(self) -> str:
        """What has to appear on the panel for the window to be checkable."""
        return (f"0-{self.vmax:.0f} (p{self.high_percentile:g}), "
                f"n={self.n_finite:,}, =0 {self.frac_at_zero:.1%}, "
                f">=max {self.frac_at_or_above_vmax:.2%}, "
                f"sat {self.frac_at_ceiling:.2%}, gamma 1")

    def as_dict(self) -> dict:
        d = {k: v for k, v in self.__dict__.items()}
        d.update(frac_at_zero=self.frac_at_zero,
                 frac_at_or_above_vmax=self.frac_at_or_above_vmax,
                 frac_at_ceiling=self.frac_at_ceiling, gamma=None)
        return d


def unsigned_display_range(values, *, high_percentile: float = 99.5,
                           ceiling: float = 255.0) -> DisplayRange:
    """Window for a channel that is clamped at 0 by the scanner.

    ``values`` should be the SELECTION-INDEPENDENT population -- every crop of that channel,
    not just the ones being drawn -- so that changing how many islets are on the page cannot
    move the ramp.
    """
    a = np.asarray(values, dtype=np.float64).ravel()
    finite = np.isfinite(a)
    n_nonfinite = int((~finite).sum())
    a = a[finite]
    if a.size == 0:
        return DisplayRange(vmax=ceiling, high_percentile=high_percentile, ceiling=ceiling,
                            n_finite=0, n_at_zero=0, n_at_or_above_vmax=0, n_at_ceiling=0,
                            n_nonfinite=n_nonfinite,
                            fallback_reason="no finite values; window is the sensor ceiling")
    vmax = float(np.percentile(a, high_percentile))
    reason = None
    if not (vmax > 0.0):
        # A channel that is >= 99.5% zero. Falling back to the ceiling keeps the panel honest
        # (everything really is at the floor) instead of dividing by a zero range.
        vmax, reason = float(ceiling), (
            f"p{high_percentile:g} is 0; more than that share of this channel is at the floor")
    return DisplayRange(vmax=vmax, high_percentile=high_percentile, ceiling=ceiling,
                        n_finite=int(a.size), n_at_zero=int((a == 0).sum()),
                        n_at_or_above_vmax=int((a >= vmax).sum()),
                        n_at_ceiling=int((a >= ceiling).sum()),
                        n_nonfinite=n_nonfinite, fallback_reason=reason)


@dataclass(frozen=True)
class SymmetricRange:
    """A window for a SIGNED score, centred on zero because zero is what it means.

    The RNA beta-alpha contrast is a difference, so its neutral value is 0 and the colour scale
    must put 0 at the middle of the map. ``n_below``/``n_above`` are separate because the two
    arms are not the same size -- donor 6442's islet cells run 5,194 beta-like to 1,018
    alpha-like -- and a pooled "outside the window" count would hide saturation of the minority
    arm entirely.
    """

    half_width: float
    high_percentile: float
    n_finite: int
    n_below: int
    n_above: int
    n_at_zero: int
    n_nonfinite: int

    @property
    def vmin(self) -> float:
        return -self.half_width

    @property
    def vmax(self) -> float:
        return self.half_width

    def caption(self) -> str:
        return (f"+/-{self.half_width:.2f} (p{self.high_percentile:g} of |value|), "
                f"n={self.n_finite:,}, below {self.n_below / max(1, self.n_finite):.2%}, "
                f"above {self.n_above / max(1, self.n_finite):.2%}, 0 at centre, gamma 1")

    def as_dict(self) -> dict:
        d = {k: v for k, v in self.__dict__.items()}
        d.update(vmin=self.vmin, vmax=self.vmax, gamma=None)
        return d


def symmetric_display_range(values, *, high_percentile: float = 99.5) -> SymmetricRange:
    a = np.asarray(values, dtype=np.float64).ravel()
    finite = np.isfinite(a)
    n_nonfinite = int((~finite).sum())
    a = a[finite]
    if a.size == 0:
        return SymmetricRange(half_width=1.0, high_percentile=high_percentile, n_finite=0,
                              n_below=0, n_above=0, n_at_zero=0, n_nonfinite=n_nonfinite)
    hw = float(np.percentile(np.abs(a), high_percentile))
    if not (hw > 0.0):
        hw = 1.0
    return SymmetricRange(half_width=hw, high_percentile=high_percentile, n_finite=int(a.size),
                          n_below=int((a < -hw).sum()), n_above=int((a > hw).sum()),
                          n_at_zero=int((a == 0).sum()), n_nonfinite=n_nonfinite)


#: Grey for "measured, but not available" -- borrowed from `integration/figures.py`'s
#: Unassigned swatch. The legend always names the REASON; there is no catch-all bucket.
UNAVAILABLE = "#BBBBBB"


def mapped_colormap(name: str, *, under: str | None = None, over: str | None = None):
    """A colormap whose bad/under/over colours are set explicitly rather than defaulted.

    A NaN drawn at the colormap's bottom colour is indistinguishable from a real minimum, which
    is the same failure as writing 0 for a missing measurement.
    """
    import matplotlib.pyplot as plt

    extremes = {"bad": UNAVAILABLE}
    if under is not None:
        extremes["under"] = under
    if over is not None:
        extremes["over"] = over
    return plt.get_cmap(name).with_extremes(**extremes)
