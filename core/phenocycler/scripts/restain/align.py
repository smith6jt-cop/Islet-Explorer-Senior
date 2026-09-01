"""Is the post-Xenium IF aligned to its Xenium section at CELL scale?

Nothing had ever checked what a registration leaves at the scale of one cell -- and if the two
images are out by a cell diameter, then sampling the image over a cell's real outline instead of
a disk around its centre cannot help, and the result alone would never say so.  So this runs
FIRST and its verdict gates the rest.  ``--transform`` selects which registration is under test:
the tissue-outline ECC fit, VALIS's DAPI-to-DAPI fit, or the nucleus-matched transform.

**How independent this is depends on which transform it is measuring.**  For the tissue-outline
and VALIS arms it was genuinely independent -- neither had ever seen a nucleus.  The
nucleus-matched transform is fitted on Xenium nuclei against IF DAPI detections, which is closely
related to what this measures, so for that arm this is a consistency check.  What is not
circular there is the spatial hold-out in ``run_refine.py`` and the placement null.

**The measurement.**  Xenium's nucleus segmentation comes from Xenium's own DAPI on the same
physical section; the IF re-stained that section and its own DAPI is recorded
``channel_dapi = NUCLEUS``.  So the two should correspond nucleus for nucleus.  Per islet:
warp the Xenium nucleus polygons forward into the IF frame, rasterise them, and find the shift
that best lines that mask up with the IF DAPI, by normalised cross-correlation over a search
disk.  The answer is a residual offset in MICRONS, reported beside the donor's own median cell
diameter, so "is this within a cell?" is answerable without trusting a threshold.

**No correction fitted here may reach the hormone channels.**  Fitting a shift on DAPI and
then applying it to the hormone measurement would tune the alignment on the very thing this
exists to test.  The gate is handed a `functools.partial` over a COPY of the affine, it writes
only ``alignment_panc.json``, and nothing in `run_hormone` reads that file.  That is friction,
not proof: a determined caller can still defeat it, and claiming otherwise would be a stronger
statement than the code supports.

**Every unavailability reason preferentially removes badly-aligned islets** -- a flat
correlation surface, too few nuclei, a blank window.  So the headline is reported twice: over
the islets that were measured, and with every unmeasured islet entered at the search ceiling.
A survivors-only median can only ever flatter.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import tiles as TL                                                # noqa: E402
from rasterize_labels import rasterize_binary                     # noqa: E402

from phenocycler.integration.rasterize import GridSpec            # noqa: E402

# --------------------------------------------------------------------------- #
# Declared BEFORE the run. Every one of these is a knob, and a knob chosen after
# seeing the answer is the tuning this project forbids.
# --------------------------------------------------------------------------- #

#: Correlate at 0.5 um/px (pyramid level 1). Nuclei are ~6 um across, so this is ~12 px per
#: nucleus -- enough to localise one, coarse enough that a 30 um search is cheap. The 0.25
#: um/px measurement crops are a separate lane; this does not share them.
ALIGN_TARGET_UM = 0.5

#: A cell is ~9 um, so 30 um would already be three cell diameters -- but on donor 6442 the
#: per-islet magnitudes piled up against a 30 um rim (p90 29.9 of a 30.0 limit) while the
#: WHOLE-SECTION peak sat at 14.1 um and did not move between a 30 um and a 1000 um search.
#: That is a real per-islet spread around a real global residual, with its tail clipped. 60 um
#: covers the observed mean plus ~3 sd. Widening can only make a residual look LARGER, never
#: smaller, so it cannot flatter the verdict -- and a measurement pinned to its own search
#: boundary is not a measurement.
SEARCH_RADIUS_UM = 60.0

#: Removes the slow intensity roll-off across a window so the peak is driven by nuclei rather
#: than by overall brightness. 15 um is larger than a nucleus and smaller than an islet.
HIGHPASS_SIGMA_UM = 15.0

#: Below this the correlation surface carries no nucleus-scale structure at all.
MIN_NCC_PEAK = 0.15
#: The peak must stand above the rest of the surface, or the argmax is noise with a location.
MIN_PEAK_Z = 4.0
#: Fewer nuclei than this and the peak is one or two objects, not a registration estimate.
MIN_NUCLEI_PER_ISLET = 40
#: The window is padded by at least the search radius, so the search DISK fits inside it. A
#: disk clipped by the array bound would silently search a square-ish region on one side and a
#: sliver on the other, biasing the argmax toward the padded direction.
ISLET_BUFFER_UM = 60.0

#: A peak this close to the rim is reported as truncated: the surface has not been allowed to
#: turn over, so the value is a lower bound on that islet's residual.
TRUNCATION_FRACTION_OF_R = 0.9

#: A verdict on a minority of islets is not a verdict on the section.
MIN_MEASURED_FRACTION = 0.60
#: DAPI must actually be brighter inside the warped nuclei than in the gaps between them,
#: on this many islets, or the two DAPI signals do not correspond and the offset is moot.
MIN_CONTRAST = 1.05
MIN_CONTRAST_ISLET_FRACTION = 0.60


class AlignError(RuntimeError):
    """The gate cannot be computed as specified -- never silently downgraded."""


@dataclass(frozen=True)
class Offset:
    """One islet's residual. ``dx_um``/``dy_um`` is the shift to ADD to the warped Xenium mask
    to land it on the IF DAPI. ``None`` everywhere means unavailable, and ``reason`` says why."""

    islet_index: int
    n_nuclei: int
    cx_um: float = float("nan")
    cy_um: float = float("nan")
    dx_um: float | None = None
    dy_um: float | None = None
    magnitude_um: float | None = None
    ncc_peak: float | None = None
    peak_z: float | None = None
    contrast: float | None = None
    contrast_shifted: float | None = None
    half_split_um: float | None = None
    truncated: bool = False
    reason: str = "measured"

    def as_dict(self) -> dict:
        out = {}
        for k, v in self.__dict__.items():
            if v is None or isinstance(v, (bool, str)):
                out[k] = v
            elif isinstance(v, (int, np.integer)):
                out[k] = int(v)
            else:
                out[k] = float(v)
        return out


def self_overlap_offset_um(diameter_um: float, retained_fraction: float) -> float:
    """Displacement at which a disc of ``diameter_um`` still covers ``retained_fraction`` of
    its own footprint.

    This is where the verdict's thresholds come from, and it is geometry rather than a
    preference: two circles of radius r whose centres are t apart overlap on
    ``(2/pi)(acos(u) - u*sqrt(1-u^2))`` of one circle, with ``u = t / 2r``.  At the half point
    a sampler is taking as much of its neighbours as of its own cell, whatever aperture it
    uses.
    """
    if not 0.0 < retained_fraction <= 1.0:
        raise AlignError(f"retained_fraction must be in (0, 1], got {retained_fraction!r}")

    def f(u):
        u = min(max(u, 0.0), 1.0)
        return (2.0 / np.pi) * (np.arccos(u) - u * np.sqrt(max(0.0, 1.0 - u * u)))

    lo, hi = 0.0, 1.0
    for _ in range(80):                     # f is monotone decreasing in u
        mid = 0.5 * (lo + hi)
        if f(mid) > retained_fraction:
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi) * diameter_um)


def _parabolic(vm1: float, v0: float, vp1: float) -> float:
    """Sub-sample peak position from three samples, 0 when the fit is degenerate."""
    denom = vm1 - 2.0 * v0 + vp1
    if denom == 0.0 or not np.isfinite(denom):
        return 0.0
    d = 0.5 * (vm1 - vp1) / denom
    return float(d) if abs(d) <= 1.0 else 0.0


def ncc_surface(dapi: np.ndarray, valid: np.ndarray, template: np.ndarray, *,
                pixel_um: float,
                highpass_sigma_um: float = HIGHPASS_SIGMA_UM):
    """``(surface, None)`` or ``(None, reason)`` -- the normalised correlation surface,
    ``fftshift``ed so zero shift sits at the centre.

    Shared by :func:`offset_by_ncc` and the figures, so the picture a reader checks is the
    same array the verdict was taken from rather than a second computation that could differ.
    """
    from scipy.ndimage import gaussian_filter

    a = np.asarray(dapi, dtype=np.float64)
    a = np.where(valid, a, np.nan)
    finite = np.isfinite(a)
    if finite.sum() < 0.5 * a.size:
        return None, "DAPI window mostly unavailable"
    a = np.where(finite, a, np.nanmean(a))
    if float(np.std(a)) == 0.0:
        return None, "DAPI window blank or constant"

    a = a - gaussian_filter(a, max(1.0, highpass_sigma_um / pixel_um))
    b = np.asarray(template, dtype=np.float64)
    if b.sum() == 0:
        return None, "no nuclei rasterised in this window"

    a = a - a.mean()
    b = b - b.mean()
    sa, sb = float(a.std()), float(b.std())
    if sa == 0.0 or sb == 0.0:
        return None, "no variance to correlate"

    cc = np.fft.irfft2(np.fft.rfft2(a) * np.conj(np.fft.rfft2(b)), s=a.shape)
    return np.fft.fftshift(cc) / (a.size * sa * sb), None


def offset_by_ncc(dapi: np.ndarray, valid: np.ndarray, template: np.ndarray, *,
                  pixel_um: float, islet_index: int = -1, n_nuclei: int = 0,
                  search_radius_um: float = SEARCH_RADIUS_UM,
                  highpass_sigma_um: float = HIGHPASS_SIGMA_UM,
                  min_peak: float = MIN_NCC_PEAK, min_peak_z: float = MIN_PEAK_Z) -> Offset:
    """Shift that lines ``template`` (the warped nucleus mask) up with ``dapi``.

    Sign convention, recorded because it is the thing most easily inverted: the returned
    ``(dx_um, dy_um)`` is what must be ADDED to the mask's coordinates to land it on the DAPI.
    """
    if dapi.shape != template.shape or dapi.shape != valid.shape:
        raise AlignError(f"shapes disagree: dapi {dapi.shape}, template {template.shape}, "
                         f"valid {valid.shape}")
    kw = dict(islet_index=islet_index, n_nuclei=n_nuclei)
    if n_nuclei < MIN_NUCLEI_PER_ISLET:
        return Offset(**kw, reason="too few nuclei to localise")
    if not valid.any():
        return Offset(**kw, reason="DAPI window unavailable")

    cc, reason = ncc_surface(dapi, valid, template, pixel_um=pixel_um,
                             highpass_sigma_um=highpass_sigma_um)
    if cc is None:
        return Offset(**kw, reason=reason)
    cy, cx = cc.shape[0] // 2, cc.shape[1] // 2
    a = cc                                    # only its shape is used from here on

    r_px = int(round(search_radius_um / pixel_um))
    if r_px >= min(cy, cx):
        return Offset(**kw, reason="search disk does not fit the window")
    yy, xx = np.ogrid[:a.shape[0], :a.shape[1]]
    disk = ((yy - cy) ** 2 + (xx - cx) ** 2) <= r_px ** 2
    masked = np.where(disk, cc, -np.inf)
    py, px_ = np.unravel_index(int(np.argmax(masked)), masked.shape)
    peak = float(cc[py, px_])

    off_disk = cc[disk]
    bg_sd = float(off_disk.std())
    peak_z = float((peak - off_disk.mean()) / bg_sd) if bg_sd > 0 else np.inf

    if peak < min_peak:
        return Offset(**kw, ncc_peak=peak, peak_z=peak_z, reason="no nucleus-scale peak")
    if peak_z < min_peak_z:
        return Offset(**kw, ncc_peak=peak, peak_z=peak_z, reason="peak is ambiguous")
    truncated = ((py - cy) ** 2 + (px_ - cx) ** 2) >= (TRUNCATION_FRACTION_OF_R * r_px) ** 2

    sy = _parabolic(cc[py - 1, px_], peak, cc[py + 1, px_]) if 0 < py < a.shape[0] - 1 else 0.0
    sx = _parabolic(cc[py, px_ - 1], peak, cc[py, px_ + 1]) if 0 < px_ < a.shape[1] - 1 else 0.0
    dy = (py - cy + sy) * pixel_um
    dx = (px_ - cx + sx) * pixel_um
    return Offset(**kw, dx_um=dx, dy_um=dy, magnitude_um=float(np.hypot(dx, dy)),
                  ncc_peak=peak, peak_z=peak_z, truncated=bool(truncated),
                  reason="measured")


def contrast_inside_vs_gaps(dapi: np.ndarray, valid: np.ndarray, template: np.ndarray,
                            *, shift_px: tuple[int, int] = (0, 0)) -> float | None:
    """Mean DAPI inside the warped nuclei over mean DAPI in the gaps, or ``None``.

    ``None`` rather than ``inf`` when the gaps are not strictly positive: an ``inf`` breaks
    ``json.dumps(allow_nan=False)`` exactly the way a NaN does, and this record must stay
    writable.
    """
    mask = np.roll(np.asarray(template, dtype=bool), shift_px, axis=(0, 1))
    a = np.asarray(dapi, dtype=np.float64)
    inside = mask & valid
    gaps = (~mask) & valid
    if not inside.any() or not gaps.any():
        return None
    gap_mean = float(a[gaps].mean())
    if not (gap_mean > 0.0):
        return None
    return float(a[inside].mean() / gap_mean)


def verdict(median_um: float | None, *, diameter_um: float, n_measured: int, n_total: int,
            contrast_fraction: float, mode: str = "nucleus",
            precision_um: float | None = None) -> tuple[str, str]:
    """``(PASS | INCONCLUSIVE | FAIL, why)`` -- a pure function of the declared constants."""
    t_pass = self_overlap_offset_um(diameter_um, 0.75)
    t_fail = self_overlap_offset_um(diameter_um, 0.50)
    frac = (n_measured / n_total) if n_total else 0.0

    if median_um is None or n_measured == 0:
        return "INCONCLUSIVE", "no islet produced a measurable offset"
    if frac < MIN_MEASURED_FRACTION:
        return ("INCONCLUSIVE",
                f"only {n_measured}/{n_total} islets ({frac:.0%}) were measurable; "
                f"a verdict needs {MIN_MEASURED_FRACTION:.0%}")
    if contrast_fraction < MIN_CONTRAST_ISLET_FRACTION:
        return ("FAIL",
                f"DAPI is not brighter inside the warped nuclei than between them on "
                f"{1 - contrast_fraction:.0%} of islets; the two DAPI signals do not correspond")
    if median_um >= t_fail:
        return ("FAIL",
                f"median residual {median_um:.2f} um >= {t_fail:.2f} um, at which a "
                f"{diameter_um:.2f} um cell covers less than half its own footprint")
    if precision_um is not None and precision_um > t_pass:
        return ("INCONCLUSIVE",
                f"the estimator's own split-half precision is {precision_um:.2f} um, wider "
                f"than the {t_pass:.2f} um it would have to resolve")
    if mode != "nucleus":
        return ("INCONCLUSIVE",
                f"measured against the {mode} mask, which cannot separate a half-cell offset "
                "from a segmentation difference; only the nucleus arm can certify")
    if median_um <= t_pass:
        return ("PASS",
                f"median residual {median_um:.2f} um <= {t_pass:.2f} um, at which a "
                f"{diameter_um:.2f} um cell still covers 75% of its own footprint")
    return ("INCONCLUSIVE",
            f"median residual {median_um:.2f} um is between {t_pass:.2f} and {t_fail:.2f} um "
            f"for a {diameter_um:.2f} um cell")
