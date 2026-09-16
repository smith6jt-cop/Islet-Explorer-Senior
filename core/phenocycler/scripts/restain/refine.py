#!/usr/bin/env python
"""Candidate registrations of the two DAPI nucleus fields, and one way of scoring them all.

The IF re-stained the same physical section Xenium measured, so this is a point-set problem:
the same nuclei, twice. Every candidate here produces a ``Transform`` mapping Xenium microns to
IF microns, and every candidate is scored by the same measurement -- the distance between a
Xenium nucleus and the IF nucleus it lands on.

**Why that score and not the correlation the fitters optimise.** A per-cell distance is the
thing the lane actually needs ("read this cell in the IF"), it is in microns rather than in
arbitrary correlation units, and it is comparable across candidates fitted by different
packages on different objectives. It also carries its own null: nuclei are ~9 um apart, so an
unregistered field would sit near half that, and a candidate cannot look good by accident.

**What is NOT a score.** ``run_align.py`` stays out of this file. It is the independent gate --
rasterised masks and windowed NCC, a different representation of the same section -- and it
only sees a transform after one has been chosen here.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from phenocycler.integration.rasterize import GridSpec              # noqa: E402
from phenocycler.integration.register import icp                    # noqa: E402
from phenocycler.integration.transform import (Transform,           # noqa: E402
                                               fit_affine_lstsq,
                                               fit_rigid_umeyama)

#: A pair further apart than this is not the same nucleus. Nuclei are ~9 um apart here, so half
#: that spacing is the point past which "nearest" stops meaning "corresponding".
PAIR_MAX_UM = 4.5
#: Radii at which the matched fraction is reported. 2 um is well inside a nucleus, 5 um is a
#: nucleus away; quoting the curve rather than one number keeps the choice from mattering.
MATCH_RADII_UM = (1.0, 2.0, 3.0, 5.0)
#: Offset applied to build the "no correspondence" null. Larger than any residual under test and
#: larger than a nucleus, so a null pair is a coincidence rather than a near miss.
NULL_OFFSET_UM = 25.0


class RefineError(RuntimeError):
    """A candidate cannot be fitted or scored as asked."""


@dataclass
class Candidate:
    name: str
    transform: Transform
    stages: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    score: dict = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# The score
# --------------------------------------------------------------------------- #

def mutual_pairs(moving_xy: np.ndarray, fixed_xy: np.ndarray, *,
                 max_dist_um: float = PAIR_MAX_UM):
    """Indices of mutual nearest neighbours within ``max_dist_um``, and their distances.

    Mutual rather than one-sided because the two fields do not have the same members: Xenium
    segments nuclei the IF detector misses and vice versa, and plain nearest-neighbour would
    pair every orphan with whatever happens to be closest and report that as a residual.
    """
    from scipy.spatial import cKDTree

    if len(moving_xy) == 0 or len(fixed_xy) == 0:
        return np.empty(0, int), np.empty(0, int), np.empty(0)
    tf_tree = cKDTree(fixed_xy)
    mv_tree = cKDTree(moving_xy)
    d_mf, j = tf_tree.query(moving_xy, k=1, workers=-1)
    _d_fm, i = mv_tree.query(fixed_xy, k=1, workers=-1)
    mutual = (i[j] == np.arange(len(moving_xy))) & (d_mf <= max_dist_um)
    idx = np.flatnonzero(mutual)
    return idx, j[idx], d_mf[idx]


def score_transform(tf: Transform, moving_xy: np.ndarray, fixed_xy: np.ndarray, *,
                    max_dist_um: float = PAIR_MAX_UM,
                    radii_um=MATCH_RADII_UM,
                    null_offset_um: float = NULL_OFFSET_UM) -> dict:
    """Per-cell agreement of ``tf``, with the null that says what chance would give.

    Reports the pair-distance quantiles over mutual pairs AND the fraction of all Xenium nuclei
    with an IF nucleus inside each radius. The first is the residual on nuclei that correspond;
    the second is what stops a candidate looking good by matching only the easy ones.
    """
    from scipy.spatial import cKDTree

    warped = tf.apply(moving_xy)
    _i, _j, d = mutual_pairs(warped, fixed_xy, max_dist_um=max_dist_um)

    tree = cKDTree(fixed_xy)
    d_all, _ = tree.query(warped, k=1, workers=-1)
    d_null, _ = tree.query(warped + np.array([null_offset_um, 0.0]), k=1, workers=-1)

    def q(a):
        return {} if a.size == 0 else {
            "p50": float(np.percentile(a, 50)), "p90": float(np.percentile(a, 90)),
            "p99": float(np.percentile(a, 99)), "mean": float(a.mean())}

    out = {
        "n_moving": int(len(moving_xy)), "n_fixed": int(len(fixed_xy)),
        "n_pairs": int(d.size), "paired_fraction": float(d.size / max(1, len(moving_xy))),
        "pair_distance_um": q(d),
        "nn_distance_um": q(d_all),
        "nn_distance_um_null": q(d_null),
        "matched_within": {f"{r:g}um": float((d_all <= r).mean()) for r in radii_um},
        "matched_within_null": {f"{r:g}um": float((d_null <= r).mean()) for r in radii_um},
        "max_dist_um": float(max_dist_um),
        "null_offset_um": float(null_offset_um),
    }
    if d.size:
        # Direction of what is left, which is what says whether more DOF would help.
        res = fixed_xy[_j] - warped[_i]
        out["residual_median_xy_um"] = [float(np.median(res[:, 0])), float(np.median(res[:, 1]))]
        out["residual_after_constant_um"] = float(np.median(
            np.linalg.norm(res - np.median(res, axis=0), axis=1)))
    return out


# --------------------------------------------------------------------------- #
# Candidates
# --------------------------------------------------------------------------- #

def seed_from_points(src_um: np.ndarray, dst_um: np.ndarray) -> Transform:
    """Affine through a point mapping -- used to turn VALIS's warped centroids into a seed.

    Only a seed: it need only land inside half the nucleus spacing (~4.5 um) for the mutual
    nearest-neighbour matching below to be unambiguous, and VALIS is already at 1.4-2.8 um.
    """
    ok = np.isfinite(src_um).all(axis=1) & np.isfinite(dst_um).all(axis=1)
    return fit_affine_lstsq(src_um[ok], dst_um[ok])


def refine_icp(seed: Transform, moving_xy: np.ndarray, fixed_xy: np.ndarray, *,
               max_dist_um: float = PAIR_MAX_UM, max_iter: int = 30,
               affine: bool = False) -> Transform:
    """``register.icp`` on the two nucleus fields, similarity first and then optionally affine.

    ``register.icp`` (register.py:252) already restricts to mutual pairs inside a radius, which
    is the whole difficulty here: the two fields differ in membership, and plain closest-point
    ICP would drag every unmatched nucleus onto its nearest stranger and bias the fit.
    """
    return icp(moving_xy, fixed_xy, initial=seed, max_dist_um=max_dist_um,
               max_iter=max_iter, similarity_only=not affine)


def refine_translation(seed: Transform, moving_xy: np.ndarray, fixed_xy: np.ndarray, *,
                       max_dist_um: float = PAIR_MAX_UM) -> Transform:
    """Seed plus one rigid translation: the median residual of its own mutual pairs.

    The smallest possible correction, and the one the measured residual asks for -- it is a
    near-constant vector. Kept as its own candidate so that "a translation was enough" is a
    result rather than an assumption buried inside a six-parameter fit.
    """
    warped = seed.apply(moving_xy)
    i, j, _d = mutual_pairs(warped, fixed_xy, max_dist_um=max_dist_um)
    if i.size < 3:
        raise RefineError(f"only {i.size} mutual pairs -- seed is too far out to translate from")
    shift = np.median(fixed_xy[j] - warped[i], axis=0)
    aff = np.array(seed.affine, dtype=np.float64).copy()
    aff[:, 2] += shift
    return Transform(affine=aff, displacement=seed.displacement, grid=seed.grid,
                     meta={**dict(seed.meta), "translation_um": shift.tolist()})


def refine_local(seed: Transform, moving_xy: np.ndarray, fixed_xy: np.ndarray, *,
                 pixel_um: float = 250.0, max_dist_um: float = PAIR_MAX_UM,
                 min_pairs_per_bin: int = 20, smooth_bins: float = 1.0) -> Transform:
    """Seed plus a smooth displacement field fitted to its own matched-pair residuals.

    Only worth applying if the residual after a global fit still has spatial structure; the
    caller decides that from ``score_transform``'s ``residual_after_constant_um``.

    A bin with too few pairs is UNKNOWN, not zero. Leaving it at zero would make the sampler
    ramp from a real displacement down to nothing across the tissue edge -- a fabricated
    gradient exactly where the geometry is least constrained -- so sparse bins are filled from
    the nearest well-populated one and then the whole field is smoothed.
    """
    from scipy.ndimage import distance_transform_edt, gaussian_filter
    from scipy.stats import binned_statistic_2d

    warped = seed.apply(moving_xy)
    i, j, _d = mutual_pairs(warped, fixed_xy, max_dist_um=max_dist_um)
    if i.size < 3:
        raise RefineError("too few mutual pairs to fit a field")
    pts, res = warped[i], fixed_xy[j] - warped[i]

    x0, y0 = pts.min(axis=0) - pixel_um
    x1, y1 = pts.max(axis=0) + pixel_um
    w = int(np.ceil((x1 - x0) / pixel_um)) + 1
    h = int(np.ceil((y1 - y0) / pixel_um)) + 1
    # GridSpec maps micron `origin + k * pixel_um` to index k, so the bin EDGES must start half
    # a bin before the origin for a bin's median to sit at the node the sampler reads.
    xedges = x0 - 0.5 * pixel_um + np.arange(w + 1) * pixel_um
    yedges = y0 - 0.5 * pixel_um + np.arange(h + 1) * pixel_um
    grid = GridSpec(float(x0), float(y0), float(pixel_um), h, w)

    counts, *_ = binned_statistic_2d(pts[:, 1], pts[:, 0], None, statistic="count",
                                     bins=[yedges, xedges])
    populated = counts >= min_pairs_per_bin
    if not populated.any():
        raise RefineError(f"no bin holds {min_pairs_per_bin} pairs at {pixel_um} um")
    _d2, idx = distance_transform_edt(~populated, return_distances=True, return_indices=True)

    field = np.zeros((2, h, w), dtype=np.float32)
    for k in (0, 1):
        stat, *_ = binned_statistic_2d(pts[:, 1], pts[:, 0], res[:, k], statistic="median",
                                       bins=[yedges, xedges])
        stat = np.where(np.isfinite(stat), stat, 0.0)
        filled = np.where(populated, stat, stat[tuple(idx)])
        field[k] = gaussian_filter(filled, smooth_bins) if smooth_bins else filled

    return Transform(affine=np.array(seed.affine, dtype=np.float64), displacement=field,
                     grid=grid, meta={**dict(seed.meta), "field_pixel_um": float(pixel_um),
                                      "field_bins_populated": int(populated.sum()),
                                      "field_bins_total": int(populated.size)})


def global_translation_by_image(seed: Transform, moving_xy: np.ndarray, if_image: np.ndarray,
                                *, um_per_px: float, origin_um: float,
                                search_um: float = 30.0,
                                highpass_sigma_um: float = 15.0) -> tuple[float, float, float]:
    """Residual translation from INTENSITY, not points: ``(dx_um, dy_um, ncc_peak)``.

    Rasterises the warped Xenium nuclei over the whole section and correlates that against the
    high-passed IF DAPI. This is a genuinely different mechanism from the point matching above
    -- it never forms a correspondence -- so agreeing with the ICP translation is evidence, and
    disagreeing is a warning worth acting on.
    """
    from scipy.ndimage import gaussian_filter

    img = np.asarray(if_image, dtype=np.float32)
    tmpl = _nucleus_raster(seed, moving_xy, img.shape, um_per_px=um_per_px, origin_um=origin_um)

    a = img - gaussian_filter(img, max(1.0, highpass_sigma_um / um_per_px))
    a -= a.mean()
    b = tmpl - tmpl.mean()
    sa, sb = float(a.std()), float(b.std())
    if sa == 0 or sb == 0:
        raise RefineError("no variance to correlate")
    cc = np.fft.irfft2(np.fft.rfft2(a) * np.conj(np.fft.rfft2(b)), s=a.shape)
    cc = np.fft.fftshift(cc) / (a.size * sa * sb)

    cy, cx = cc.shape[0] // 2, cc.shape[1] // 2
    r_px = int(round(search_um / um_per_px))
    sub = cc[cy - r_px:cy + r_px + 1, cx - r_px:cx + r_px + 1]
    py, px = np.unravel_index(int(np.argmax(sub)), sub.shape)
    peak = float(sub[py, px])
    fy = _parabolic(sub[py - 1, px], peak, sub[py + 1, px]) if 0 < py < sub.shape[0] - 1 else 0.0
    fx = _parabolic(sub[py, px - 1], peak, sub[py, px + 1]) if 0 < px < sub.shape[1] - 1 else 0.0
    return ((px - r_px + fx) * um_per_px, (py - r_px + fy) * um_per_px, peak)


def _nucleus_raster(seed: Transform, moving_xy: np.ndarray, shape, *,
                    um_per_px: float, origin_um: float) -> np.ndarray:
    """Warped Xenium nuclei as a count raster on the IF image's own lattice."""
    h, w = shape
    warped = seed.apply(moving_xy)
    col = np.rint((warped[:, 0] - origin_um) / um_per_px).astype(np.int64)
    row = np.rint((warped[:, 1] - origin_um) / um_per_px).astype(np.int64)
    keep = (col >= 0) & (col < w) & (row >= 0) & (row < h)
    if not keep.any():
        raise RefineError("no warped nucleus lands inside the IF image")
    out = np.zeros((h, w), dtype=np.float32)
    np.add.at(out, (row[keep], col[keep]), 1.0)
    return out


def global_translation_by_phase(seed: Transform, moving_xy: np.ndarray, if_image: np.ndarray,
                                *, um_per_px: float, origin_um: float,
                                upsample: int = 100,
                                highpass_sigma_um: float = 15.0) -> tuple[float, float]:
    """The same question as :func:`global_translation_by_image`, asked by a third mechanism.

    ``skimage.registration.phase_cross_correlation`` normalises the cross-power spectrum, so it
    weights every spatial frequency equally instead of letting the brightest structures carry
    the peak, and it refines to a fraction of a pixel by a local DFT rather than by fitting a
    parabola. Three estimates of one translation -- nearest-neighbour nuclei, plain NCC, and
    phase correlation -- from three different packages: where they agree, the number is a
    property of the sections rather than of a method.
    """
    from scipy.ndimage import gaussian_filter
    from skimage.registration import phase_cross_correlation

    img = np.asarray(if_image, dtype=np.float32)
    a = img - gaussian_filter(img, max(1.0, highpass_sigma_um / um_per_px))
    b = _nucleus_raster(seed, moving_xy, img.shape, um_per_px=um_per_px, origin_um=origin_um)
    b = b - gaussian_filter(b, max(1.0, highpass_sigma_um / um_per_px))
    shift, _err, _phase = phase_cross_correlation(a, b, upsample_factor=upsample,
                                                  normalization="phase")
    # `shift` is what to ADD to `b`'s indices to line it up with `a`, in (row, col).
    return float(shift[1] * um_per_px), float(shift[0] * um_per_px)


def _parabolic(vm1, v0, vp1) -> float:
    denom = float(vm1) - 2.0 * float(v0) + float(vp1)
    if denom == 0.0 or not np.isfinite(denom):
        return 0.0
    d = 0.5 * (float(vm1) - float(vp1)) / denom
    return d if abs(d) <= 1.0 else 0.0


def translated(tf: Transform, dx_um: float, dy_um: float) -> Transform:
    """``tf`` with a constant added to its translation, field and metadata preserved."""
    aff = np.array(tf.affine, dtype=np.float64).copy()
    aff[:, 2] += np.array([dx_um, dy_um], dtype=np.float64)
    return Transform(affine=aff, displacement=tf.displacement, grid=tf.grid,
                     meta={**dict(tf.meta), "translation_um": [float(dx_um), float(dy_um)]})


def describe(tf: Transform) -> dict:
    """The transform in units a reader can check: microns, degrees, parts per million."""
    return {"affine": np.asarray(tf.affine).tolist(),
            "scale": float(tf.scale), "scale_ppm": float((tf.scale - 1.0) * 1e6),
            "rotation_deg": float(tf.rotation_deg),
            "translation_um": [float(v) for v in tf.translation_um],
            "mirrored": bool(tf.is_mirrored),
            "has_field": tf.displacement is not None,
            "max_field_um": float(tf.max_displacement_um)}
