"""Step 2: put every marker on a common scale with UniFORM's rigid log-shift.

UniFORM (Wang et al., doi:10.1101/2024.12.06.626879) aligns each marker's **background mode**
across samples and applies the alignment as a single constant in log space -- equivalently a
per-(donor, marker) multiplicative rescaling.  Two properties are why this lane uses it:

* the shift is estimated from the population that *should* be identical across donors, so it
  does not confound staining brightness with how many positive cells a donor has.  A
  per-donor median or p99 does confound them: a T1D donor with no beta cells has a low INS
  median *because the biology changed*, and dividing by it would inflate the survivors;
* the shift is **monotone**, so distribution shape and within-donor rank are preserved
  exactly and the positive fraction is untouched.  That is what
  ``docs/INTENSITY_SCALE.md``'s "standardize the unit, never the distribution" asks for, and
  it means every rank statistic in Step 1 is invariant under it -- so those statistics can
  check this step without circularity.

**Where the shift comes from.**  ``registration.automatic_registration`` and
``normalization.calculate_shift_in_log_pixels`` are UniFORM's own functions, called
unmodified.  The one thing done here rather than there is building the histograms: UniFORM's
``preprocessing.process_sample_distributions`` bundles histogramming with per-sample GMM
fitting and a 56-panel figure, and a GMM per (marker, donor) over 25M cells is hours of work
for a number the registration does not use.  ``marker_histograms`` reproduces its
histogramming exactly -- same ``log_transform_intensities(min_value=1)``, same global range,
same ``np.histogram(..., bins=num_bins, range=(global_min, global_max))`` -- at full
resolution, and the GMM cutoffs are fitted separately on a capped subsample.

**Three attempts to improve on UniFORM's default, all refuted by measurement.**  Eleven
markers stay above the 0.10 background-mode bar after the automatic shift, and none of these
helps -- they are recorded so they are not retried:

* *normalise each donor's histogram to sum 1 before cross-correlating* (donors differ 8x in
  cell count, and reference selection is an L2 distance to the mean histogram).  It changes
  the chosen reference on 8 of 14 markers and the outcome on none: median residual SD 0.1507
  against 0.1501.
* *keep the exact zeros* by histogramming ``log1p`` instead of dropping everything below 1.
  Motivated by the observation that residual SD tracks ``zero_frac`` at rho +0.60 while the
  abundance hypothesis fails (+0.36, wrong sign).  Substantially **worse**: median 0.4437.
  The zeros are all clamped to one value, so they add a fixed-position spike that anchors the
  cross-correlation.
* *place the landmark at the negative component's peak* via ``fit_two_component_gmm``, as a
  reproducible stand-in for the interactive ``landmark.landmark_refinement``.  **Worse on
  every marker tested**, including the four that were already fine (CD3e 0.015 -> 0.194).
  ``docs/METHOD_HISTORY.md`` already records why: on these distributions a two-component GMM
  splits distribution *shape* rather than negative from positive, assigning ~95% of control
  to the "contaminated" component for 18 of 24 markers.

None of them helps because none of them addresses the actual cause, which is measurable:
**the eleven differ between donors in distribution SHAPE, not in location.**  After the
shift, each donor's histogram correlates with the reference's at median 0.9704 for the eleven
against 0.9915 for the other 44 -- Spearman(residual SD, shape correlation) = -0.607 over all
55 markers, Mann-Whitney p = 1.6e-05.  Vimentin is the clearest case: its histogram is bimodal
and the two modes sit at different relative positions in different donors, because it marks
stroma, endothelium and immune cells in donor-varying proportions, so the mixture weights and
therefore the shape differ.

A rigid shift is a location model.  It cannot align distributions whose shape differs, and
neither can landmark refinement -- moving a distribution does not reshape it.  So the eleven
are carried with their measured residual, and the consequence for anything downstream is
stated rather than left implicit: **for these markers a single cohort-wide cutoff is not
valid even after alignment**, and only within-donor rank statistics are safe.

**What a shift cannot repair.**  A rigid shift models a per-donor scale difference.  It does
not model a different measurement, which differs in specificity and background as well as
gain.  Step 1 flags markers whose level steps with the acquisition group, and Step 2 leaves
such a marker unaligned when UniFORM's shift does not bring its background modes within the
bar either (``run_step2.REFUSED``; adjusted 2026-10-02 -- until then the flag alone refused
the marker, which with 25 donors in one cohort refused CD8 although UniFORM aligns it to
0.065).
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                        # noqa: E402

#: The UniFORM checkout this lane was built against.  It is a gitignored clone under
#: /external/ because the repo carries no packaging metadata, so pip cannot pin it; the
#: assertion below is what stops a drifted checkout from silently changing the shift.
UNIFORM_COMMIT = "7c33396d12e5a19dc0b36e46a94c3293de0cd814"
UNIFORM_PATH = panel.REPO / "external" / "UniFORM"

NUM_BINS = 1024
GMM_MAX_CELLS = 200_000          # the GMM cutoff is a Step-3 input, not a shift input


def _import_uniform():
    """Import UniFORM from the pinned checkout, refusing a drifted one."""
    if not UNIFORM_PATH.exists():
        raise RuntimeError(
            f"UniFORM checkout missing at {UNIFORM_PATH}.\n"
            f"  git clone https://github.com/smith6jt-cop/UniFORM.git external/UniFORM\n"
            f"  git -C external/UniFORM checkout {UNIFORM_COMMIT}")
    head = subprocess.run(["git", "-C", str(UNIFORM_PATH), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    if head != UNIFORM_COMMIT:
        raise RuntimeError(
            f"UniFORM checkout is at {head}, this lane is pinned to {UNIFORM_COMMIT}. "
            "A different commit can change the registration without anything in the output "
            "being able to show it.")
    sys.path.insert(0, str(UNIFORM_PATH))
    from UniFORM import registration, normalization, preprocessing   # noqa: F401
    return registration, normalization, preprocessing


def marker_histograms(values_by_donor: dict, num_bins: int = NUM_BINS) -> tuple:
    """UniFORM's histogramming for one marker: log, common global range, fixed bins.

    Returns ``(counts_by_donor, {'global_min','global_max','min_list','max_list'})`` in the
    shapes ``automatic_registration`` and ``calculate_shift_in_log_pixels`` expect.
    """
    logs, mins, maxs = {}, [], []
    for d, v in values_by_donor.items():
        v = np.asarray(v, dtype=float)
        v = v[np.isfinite(v) & (v >= 1.0)]          # log_transform_intensities(min_value=1)
        lv = np.log(v) if len(v) else np.array([0.0])
        logs[d] = lv
        mins.append(float(lv.min()))
        maxs.append(float(lv.max()))
    lo, hi = min(mins), max(maxs)
    counts = {d: np.histogram(lv, bins=num_bins, range=(lo, hi))[0]
              for d, lv in logs.items()}
    return counts, {"global_min": lo, "global_max": hi,
                    "min_list": mins, "max_list": maxs}


def gmm_cutoff(values: np.ndarray, ranges: dict, num_bins: int = NUM_BINS,
               max_cells: int = GMM_MAX_CELLS, seed: int = 0) -> float:
    """UniFORM's two-component negative cutoff, on a capped subsample.

    Not used by the shift.  Carried because Step 3 needs a per-marker "above background"
    level, and after alignment one cutoff is valid across donors.
    """
    _, _, preprocessing = _import_uniform()
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v) & (v >= 1.0)]
    if len(v) < 1000:
        return float("nan")
    rng = np.random.default_rng(seed)
    if len(v) > max_cells:
        v = v[rng.choice(len(v), max_cells, replace=False)]
    cut, *_ = preprocessing.fit_two_component_gmm(
        np.log(v), num_bins, ranges["global_min"], ranges["global_max"])
    return float(np.exp(cut))


def gmm_background_landmark(values_by_donor: dict, ranges: dict,
                            num_bins: int = NUM_BINS, max_cells: int = GMM_MAX_CELLS,
                            seed: int = 0) -> list:
    """Per-donor bin index of the NEGATIVE component's peak, as a landmark.

    UniFORM's automatic registration cross-correlates whole histograms, so it aligns whatever
    peak is tallest.  For a marker present on most cells that is the *positive* population,
    and the alignment then follows biology rather than background -- which is the case its
    interactive ``landmark.landmark_refinement`` exists to fix, by having a person point at
    the negative mode.

    ``fit_two_component_gmm`` already separates the two components, so the landmark can be
    placed by rule instead: take the peak of the lower-mean component.  This is UniFORM's own
    model and UniFORM's own landmark pathway (``compute_landmark_shifts``); the only thing
    replaced is the human's mouse, which makes the choice reproducible and auditable rather
    than a click nobody can re-derive.
    """
    _, _, preprocessing = _import_uniform()
    rng = np.random.default_rng(seed)
    lo, hi = ranges["global_min"], ranges["global_max"]
    binw = (hi - lo) / (num_bins - 1)
    out = []
    for d, v in values_by_donor.items():
        v = np.asarray(v, dtype=float)
        v = v[np.isfinite(v) & (v >= 1.0)]
        if len(v) > max_cells:
            v = v[rng.choice(len(v), max_cells, replace=False)]
        _, pts, p0, p1 = preprocessing.fit_two_component_gmm(np.log(v), num_bins, lo, hi)
        neg = p0 if pts[int(np.argmax(p0))] < pts[int(np.argmax(p1))] else p1
        out.append(float((pts[int(np.argmax(neg))] - lo) / binw))
    return out


def background_mode_spread(counts: dict, ranges: dict, shifts: dict | None = None,
                           num_bins: int = NUM_BINS) -> dict:
    """Where each donor's histogram peaks, before and (optionally) after the shift.

    This is the quantity the shift is trying to align, so it is the quantity the acceptance
    bar is stated on -- not the positive fraction, which must not be tuned toward uniformity.
    """
    lo, hi = ranges["global_min"], ranges["global_max"]
    binw = (hi - lo) / (num_bins - 1)
    peaks = {}
    for d, c in counts.items():
        c = np.asarray(c, dtype=float)
        peak = lo + (int(np.argmax(c[1:])) + 1) * binw       # skip bin 0: the log floor
        if shifts is not None:
            peak -= shifts[d] * binw
        peaks[d] = peak
    v = np.array(list(peaks.values()))
    return {"peaks": peaks, "sd_log": float(v.std()),
            "fold": float(np.exp(v.max() - v.min()))}
