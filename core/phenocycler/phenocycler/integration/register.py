"""
S4 — registering two serial sections.

Two independent tracks, both run, best-by-QC wins and the choice is recorded:

**Track A — point set.** RANSAC over candidate islet correspondences, then ICP refinement.
Needs no images at all, so it works on any donor whose cell tables exist. Islets are natural
fiducials in pancreas: there are tens to hundreds of them, they are well separated, and they
carry size and shape features that constrain matching far better than position alone.

**Track B — image.** Coarse orientation search, then opencv ECC (Euclidean then affine),
then an optional SimpleITK B-spline. Runs on DAPI-to-DAPI where the morphology images are
mounted, and falls back to tissue/lineage-density rasters derived from the cell tables where
they are not.

Three things this design takes seriously that the prototype it replaces did not.

*Rotation and mirroring are not optional.* Serial sections are picked up and mounted by hand;
they arrive at arbitrary rotation and are sometimes flipped. A translation-only registration
(``phase_cross_correlation`` and nothing else) cannot represent either, so it will happily
return a confident-looking shift that is simply wrong. The coarse stage searches rotation
exhaustively and evaluates the mirrored hypothesis explicitly.

*Non-rigid deformation is real but bounded.* Tissue tears, folds and stretches during
sectioning and mounting. A B-spline absorbs that. But given enough freedom a B-spline will
also happily warp any two images into agreement, so displacement is capped at
``reg_max_disp_um`` and the cap is enforced, not merely configured.

*A registration that fails must say so.* Every result carries a QC block, and
``qc.py`` gates on it. A donor whose sections genuinely do not correspond — different blocks,
too large a gap, tissue destroyed — should drop out of structure-level analysis rather than
contribute a plausible-looking alignment to it.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import CellTable, has_identity, read_cell_table
from .manifest import _pheno_paths, load_manifest, paired_rows
from .rasterize import (GridSpec, grid_for, lineage_density, normalize_image,
                        rasterize_points, rasterize_structures, tissue_mask)
from .structures import structures_path
from .transform import (Transform, affine_from_pixels, affine_to_pixels, fit_affine_lstsq,
                        fit_rigid_umeyama, from_rigid)

#: Lineages used as density channels when registering without morphology images. Endocrine
#: is effectively an islet map and Vascular a vessel map; both are structural and largely
#: disease-independent.
DENSITY_LINEAGES = ("Endocrine", "Vascular", "Exocrine")


@dataclass
class RegistrationResult:
    """A transform plus everything needed to judge whether to believe it."""

    transform: Transform
    track: str                       # 'point_set' | 'image'
    grid: Optional[GridSpec] = None
    n_inliers: int = 0
    n_candidates: int = 0
    rmse_um: float = float("nan")
    tissue_dice: float = float("nan")
    islet_dice: float = float("nan")
    # Out-of-sample evidence. `footprint_*` is image structure inside the moving footprint,
    # scored against other placements on the same section; `islet_peak_*` is where the joint
    # rotation+translation correspondence landscape peaks relative to this transform.
    footprint_ncc: float = float("nan")
    footprint_ncc_z: float = float("nan")
    footprint_containment: float = float("nan")
    # The null itself, not just the z derived from it. Across 20 donors the null's spread
    # varies 14-fold (sd 0.006 to 0.083), so the same z means different things on different
    # sections and a reader cannot check that without these.
    footprint_null_mean: float = float("nan")
    footprint_null_sd: float = float("nan")
    footprint_null_n: int = 0
    footprint_null_floor: float = float("nan")
    islet_peak_ratio: float = float("nan")
    islet_peak_offset_um: float = float("nan")
    islet_peak_angle_deg: float = float("nan")
    mirrored: bool = False
    stages: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def inlier_fraction(self) -> float:
        return self.n_inliers / self.n_candidates if self.n_candidates else 0.0

    def score(self) -> float:
        """Single comparable number for picking between tracks.

        **Tissue Dice deliberately carries no weight.** It used to carry half of this score, on
        the reasoning that two outlines cannot be made to overlap unless the transform is close
        to right. That reasoning fails when one section is a sub-field of the other, and it
        failed in exactly the way that matters: on donor 6436, correcting a transform that was
        wrong by 12.5 deg and 0.8 mm -- taking islet correspondence from 12 of 136 (chance is
        5) to 69 -- *lowered* tissue Dice from 0.549 to 0.530 and footprint containment from
        0.910 to 0.890. A silhouette can be brought into agreement by many placements, and the
        one that agrees best is not the one that matches anatomy. Dice is still measured and
        still reported; it is not allowed to choose.

        What is weighted instead is evidence the transform was not fitted against: image
        structure inside the moving footprint measured against other placements of that same
        footprint, and whether the joint rotation+translation correspondence landscape peaks
        *here* rather than somewhere else. Terms that could not be computed are dropped and the
        weights renormalised, so a track that structurally cannot supply one -- the point-set
        track has no images -- is not scored as though it supplied a zero.
        """
        terms = [
            (0.45, float(np.clip(self.footprint_ncc_z / 10.0, 0.0, 1.0))
             if np.isfinite(self.footprint_ncc_z) else None),
            (0.25, float(np.clip(self.islet_peak_ratio / 3.0, 0.0, 1.0))
             if np.isfinite(self.islet_peak_ratio) else None),
            (0.10, float(np.clip(1.0 - self.islet_peak_offset_um / 500.0, 0.0, 1.0))
             if np.isfinite(self.islet_peak_offset_um) else None),
            (0.12, float(np.clip(self.islet_dice / 0.2, 0.0, 1.0))
             if np.isfinite(self.islet_dice) else None),
            (0.08, float(np.exp(-self.rmse_um / 150.0))
             if np.isfinite(self.rmse_um) else None),
        ]
        available = [(w, v) for w, v in terms if v is not None]
        if not available:
            return 0.0
        return float(sum(w * v for w, v in available) / sum(w for w, _v in available))

    def to_dict(self) -> dict:
        return {
            "track": self.track,
            # Per-axis, because the single isotropic figure under "transform" cannot show a
            # one-axis stretch -- the error that silently grows across a section.
            "scale_report": scale_report(self.transform),
            # Whether the non-rigid part is a map that can be inverted at all. `max_disp_um`
            # alone cannot say: a field can sit inside its cap and still fold.
            "displacement_report": displacement_report(self.transform),
            "n_inliers": self.n_inliers,
            "n_candidates": self.n_candidates,
            "inlier_fraction": self.inlier_fraction,
            "rmse_um": self.rmse_um,
            # Measured and reported, but deliberately unweighted in score() -- see the note
            # there. On donor 6436 it moved the wrong way at the correct transform.
            "tissue_dice": self.tissue_dice,
            "islet_dice": self.islet_dice,
            "footprint_ncc": self.footprint_ncc,
            "footprint_ncc_z": self.footprint_ncc_z,
            "footprint_containment": self.footprint_containment,
            "footprint_null_mean": self.footprint_null_mean,
            "footprint_null_sd": self.footprint_null_sd,
            "footprint_null_n": self.footprint_null_n,
            "footprint_null_floor": self.footprint_null_floor,
            "islet_peak_ratio": self.islet_peak_ratio,
            "islet_peak_offset_um": self.islet_peak_offset_um,
            "islet_peak_angle_deg": self.islet_peak_angle_deg,
            "mirrored": self.mirrored,
            "score": self.score(),
            "stages": self.stages,
            "notes": self.notes,
            "transform": self.transform.to_dict(),
        }


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #

def dice(a: np.ndarray, b: np.ndarray) -> float:
    """Sorensen-Dice between two boolean masks."""
    a = np.asarray(a, dtype=bool)
    b = np.asarray(b, dtype=bool)
    denom = a.sum() + b.sum()
    if denom == 0:
        return float("nan")
    return float(2.0 * np.logical_and(a, b).sum() / denom)


def nearest_neighbour_rmse(src: np.ndarray, dst: np.ndarray,
                           max_dist_um: float = np.inf) -> tuple[float, int]:
    """RMSE of mutual-nearest-neighbour pairs, and how many there were.

    Mutual nearest rather than one-way: a one-way nearest-neighbour distance can be driven
    arbitrarily low by collapsing everything onto one dense cluster, which is exactly the
    degenerate solution an unconstrained optimiser will find.
    """
    if len(src) == 0 or len(dst) == 0:
        return float("nan"), 0
    from scipy.spatial import cKDTree

    ts, td = cKDTree(src), cKDTree(dst)
    d_sd, i_sd = td.query(src, k=1)
    _d_ds, i_ds = ts.query(dst, k=1)
    mutual = i_ds[i_sd] == np.arange(len(src))
    keep = mutual & (d_sd <= max_dist_um)
    if not keep.any():
        return float("nan"), 0
    return float(np.sqrt(np.mean(d_sd[keep] ** 2))), int(keep.sum())


def residual_vote(
    moving_xy: np.ndarray,
    fixed_xy: np.ndarray,
    *,
    angles_deg: Optional[Sequence[float]] = None,
    try_mirror: bool = True,
    bin_um: float = 50.0,
    span_um: float = 2500.0,
    smooth_um: float = 75.0,
) -> dict:
    """Joint (rotation, translation) correspondence landscape between two structure sets.

    The question this answers is whether a stored transform sits at the right place at all,
    and it exists because the diagnostic that came before it could not ask that. Sweeping a
    residual *rotation* about the warped moving centroid traces a one-parameter subgroup of
    SE(2) that fixes that centroid by construction: no angle in it can express a
    translation. A registration that is wrong in both angle and offset therefore shows a flat
    curve at every angle, which reads as "no correspondence exists" when it means "the search
    never looked".

    So translation is made free at every angle. For each hypothesis the moving centroids are
    rotated about their own centroid and *every* moving-to-fixed offset vector is binned; the
    best translation for that angle is wherever those offsets pile up. Two structure sets that
    correspond vote for one offset and produce a peak; two that do not spread their votes over
    the whole search area and produce noise.

    The reported statistic is ``peak / p95``, both measured on the same donor's own landscape,
    so it is a within-donor null -- the same structures under other placements -- and carries
    no expectation about how many structures *should* match. The volume is normalised by its
    own median, which is why a featureless landscape reads ~1.0 everywhere.

    Costs one 2-D histogram per hypothesis over ``len(moving) * len(fixed)`` offsets, which is
    a few seconds for a whole 360 deg x mirror sweep on a section's islets, and needs no
    image data.

    Returns a dict carrying ``peak``, ``p95``, ``contrast`` (peak/p95), the winning
    ``angle_deg`` / ``mirror`` / ``dx_um`` / ``dy_um``, and the per-hypothesis ``landscape``
    (n_angles, 2) of peak heights for the un-mirrored and mirrored halves.
    """
    from scipy import ndimage

    mv = np.asarray(moving_xy, dtype=np.float64)
    fx = np.asarray(fixed_xy, dtype=np.float64)
    if len(mv) == 0 or len(fx) == 0:
        return {"peak": float("nan"), "p95": float("nan"), "contrast": float("nan"),
                "angle_deg": float("nan"), "mirror": False,
                "dx_um": float("nan"), "dy_um": float("nan"),
                "n_moving": len(mv), "n_fixed": len(fx),
                "angles_deg": np.zeros(0), "landscape": np.zeros((0, 2))}

    angles = np.asarray(list(angles_deg) if angles_deg is not None
                        else np.arange(-180.0, 180.0, 2.0), dtype=np.float64)
    edges = np.arange(-span_um, span_um + bin_um / 2, bin_um)
    sigma = max(smooth_um / bin_um, 1e-6)
    centre = mv.mean(axis=0)
    mirrors = (False, True) if try_mirror else (False,)

    volume = np.zeros((len(angles), len(mirrors), len(edges) - 1, len(edges) - 1))
    for mi, mirror in enumerate(mirrors):
        base = mv.copy()
        if mirror:
            base[:, 0] = 2.0 * centre[0] - base[:, 0]
        for ai, ang in enumerate(angles):
            th = np.radians(ang)
            c, s = np.cos(th), np.sin(th)
            rot = (base - centre) @ np.array([[c, s], [-s, c]]) + centre
            # Every moving-to-fixed offset. The translation that best explains this angle is
            # wherever they agree, so it never has to be optimised for -- it is read off.
            off = (fx[None, :, :] - rot[:, None, :]).reshape(-1, 2)
            hist, _xe, _ye = np.histogram2d(off[:, 0], off[:, 1], bins=(edges, edges))
            volume[ai, mi] = ndimage.gaussian_filter(hist, sigma=sigma, mode="constant")

    # Normalise by the volume's own median so a featureless landscape reads ~1 everywhere and
    # the contrast below is dimensionless. A median of zero means nothing landed in range.
    median = float(np.median(volume))
    if median > 0:
        volume = volume / median

    flat = int(np.argmax(volume))
    ai, mi, xi, yi = np.unravel_index(flat, volume.shape)
    centres = (edges[:-1] + edges[1:]) / 2.0
    peak = float(volume[ai, mi, xi, yi])
    p95 = float(np.percentile(volume, 95))

    return {
        "peak": peak,
        "p95": p95,
        "median": float(np.median(volume)),
        "contrast": float(peak / p95) if p95 > 0 else float("inf"),
        "angle_deg": float(angles[ai]),
        "mirror": bool(mirrors[mi]),
        "dx_um": float(centres[xi]),
        "dy_um": float(centres[yi]),
        "n_moving": int(len(mv)),
        "n_fixed": int(len(fx)),
        "bin_um": float(bin_um),
        "span_um": float(span_um),
        "smooth_um": float(smooth_um),
        "angles_deg": angles,
        "landscape": volume.max(axis=(2, 3)),
        "best_map": volume[ai, mi],
        "bin_centres_um": centres,
    }


# --------------------------------------------------------------------------- #
# Track A — point-set registration on structure centroids
# --------------------------------------------------------------------------- #

def _candidate_pairs(
    moving: pd.DataFrame,
    fixed: pd.DataFrame,
    *,
    area_ratio: float = 3.0,
    size_ratio: float = 3.0,
) -> list[tuple[int, int]]:
    """Plausible islet correspondences, filtered on size before geometry is considered.

    An islet's cell count and hull area are roughly conserved across a few microns of
    sectioning, so a 30-cell islet is not a candidate partner for a 300-cell one no matter
    how close they land. Filtering on that first shrinks the RANSAC search space by an order
    of magnitude and, more importantly, removes the correspondences most likely to produce a
    confidently wrong fit.
    """
    pairs = []
    m_n = moving["n_cells"].to_numpy(np.float64)
    f_n = fixed["n_cells"].to_numpy(np.float64)
    m_a = moving.get("area_um2", pd.Series(np.zeros(len(moving)))).to_numpy(np.float64)
    f_a = fixed.get("area_um2", pd.Series(np.zeros(len(fixed)))).to_numpy(np.float64)

    for i in range(len(moving)):
        for j in range(len(fixed)):
            rn = m_n[i] / max(f_n[j], 1e-9)
            if not (1 / size_ratio <= rn <= size_ratio):
                continue
            if m_a[i] > 0 and f_a[j] > 0:
                ra = m_a[i] / f_a[j]
                if not (1 / area_ratio <= ra <= area_ratio):
                    continue
            pairs.append((i, j))
    return pairs


def ransac_affine(
    moving_xy: np.ndarray,
    fixed_xy: np.ndarray,
    candidates: Sequence[tuple[int, int]],
    *,
    n_iter: int = 4000,
    inlier_tol_um: float = 150.0,
    rng: Optional[np.random.Generator] = None,
    similarity_only: bool = True,
) -> tuple[Optional[Transform], np.ndarray]:
    """RANSAC a transform from noisy candidate correspondences.

    Fits a *similarity* (rotation, isotropic scale, translation, optional mirror) rather than
    a free affine by default. With 4 degrees of freedom instead of 6 it cannot shear a bad
    correspondence set into apparent agreement, and two serial sections of the same block
    genuinely differ by little more than a rigid motion at this stage — shear belongs in the
    later refinement, if anywhere.
    """
    rng = rng or np.random.default_rng(0)
    if len(candidates) < 3:
        return None, np.zeros(0, dtype=int)

    cand = np.asarray(candidates)
    src_all = moving_xy[cand[:, 0]]
    dst_all = fixed_xy[cand[:, 1]]

    best_tf, best_inliers = None, np.zeros(0, dtype=int)
    n_sample = 2 if similarity_only else 3

    for _ in range(n_iter):
        idx = rng.choice(len(cand), size=n_sample, replace=False)
        s, d = src_all[idx], dst_all[idx]
        if np.linalg.norm(s[0] - s[1]) < 1e-6 or np.linalg.norm(d[0] - d[1]) < 1e-6:
            continue
        try:
            tf = (fit_rigid_umeyama(s, d, allow_reflection=True) if similarity_only
                  else fit_affine_lstsq(s, d))
        except (ValueError, np.linalg.LinAlgError):
            continue

        # Reject physically implausible fits before scoring: both modalities report calibrated
        # microns, so a scale far from 1 means a bad sample, not a real observation.
        if not (0.8 <= tf.scale <= 1.25):
            continue

        resid = np.linalg.norm(tf.apply(src_all) - dst_all, axis=1)
        inliers = np.flatnonzero(resid <= inlier_tol_um)
        if len(inliers) > len(best_inliers):
            best_tf, best_inliers = tf, inliers

    if best_tf is None or len(best_inliers) < 3:
        return best_tf, best_inliers

    # Refit on all inliers; a 2-point sample fixes the model but is a poor estimate of it.
    try:
        refined = (fit_rigid_umeyama(src_all[best_inliers], dst_all[best_inliers],
                                     allow_reflection=True)
                   if similarity_only else
                   fit_affine_lstsq(src_all[best_inliers], dst_all[best_inliers]))
        resid = np.linalg.norm(refined.apply(src_all) - dst_all, axis=1)
        best_inliers = np.flatnonzero(resid <= inlier_tol_um)
        best_tf = refined
    except (ValueError, np.linalg.LinAlgError):
        pass
    return best_tf, best_inliers


def icp(
    moving_xy: np.ndarray,
    fixed_xy: np.ndarray,
    initial: Optional[Transform] = None,
    *,
    max_iter: int = 60,
    max_dist_um: float = 250.0,
    tol_um: float = 0.05,
    similarity_only: bool = True,
) -> Transform:
    """Iterative closest point refinement over mutual nearest neighbours.

    Restricted to mutual pairs within ``max_dist_um``: plain closest-point ICP on two islet
    fields that differ in count (islets appear and vanish between sections) will drag
    unmatched points onto whatever is nearest and bias the fit.
    """
    from scipy.spatial import cKDTree

    tf = initial or Transform()
    tree_fixed = cKDTree(fixed_xy)
    prev = np.inf

    for _ in range(max_iter):
        cur = tf.apply(moving_xy)
        d, j = tree_fixed.query(cur, k=1)
        tree_cur = cKDTree(cur)
        _d2, i2 = tree_cur.query(fixed_xy, k=1)
        mutual = i2[j] == np.arange(len(cur))
        keep = mutual & (d <= max_dist_um)
        if keep.sum() < 3:
            break

        try:
            step = (fit_rigid_umeyama(cur[keep], fixed_xy[j[keep]]) if similarity_only
                    else fit_affine_lstsq(cur[keep], fixed_xy[j[keep]]))
            # NB: no reflection here. The seed already fixed handedness; letting each ICP
            # step flip would let it oscillate between the two solutions.
        except (ValueError, np.linalg.LinAlgError):
            break
        tf = step.compose(tf)

        rmse = float(np.sqrt(np.mean(d[keep] ** 2)))
        if abs(prev - rmse) < tol_um:
            break
        prev = rmse
    return tf


def register_point_set(
    moving_structures: pd.DataFrame,
    fixed_structures: pd.DataFrame,
    *,
    inlier_tol_um: float = 150.0,
    seed: int = 0,
) -> RegistrationResult:
    """Track A: align two islet fields."""
    res = RegistrationResult(transform=Transform(), track="point_set")
    if moving_structures is None or fixed_structures is None:
        res.notes.append("no structures available")
        return res
    if len(moving_structures) < 3 or len(fixed_structures) < 3:
        res.notes.append(
            f"too few structures to register ({len(moving_structures)} moving, "
            f"{len(fixed_structures)} fixed; need >= 3 each)")
        return res

    mv = moving_structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    fx = fixed_structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)

    cands = _candidate_pairs(moving_structures, fixed_structures)
    res.n_candidates = len(cands)
    if not cands:
        res.notes.append("no size-compatible islet pairs")
        return res

    rng = np.random.default_rng(seed)
    tf, inliers = ransac_affine(mv, fx, cands, inlier_tol_um=inlier_tol_um, rng=rng)
    if tf is None:
        res.notes.append("RANSAC found no consistent transform")
        return res
    res.stages.append("ransac")

    tf = icp(mv, fx, initial=tf, max_dist_um=inlier_tol_um * 1.5)
    res.stages.append("icp")

    rmse, n_pairs = nearest_neighbour_rmse(tf.apply(mv), fx, max_dist_um=inlier_tol_um * 2)
    res.transform = tf
    res.n_inliers = int(n_pairs)
    res.rmse_um = rmse
    res.mirrored = tf.is_mirrored
    return res


# --------------------------------------------------------------------------- #
# Track B — image / raster registration
# --------------------------------------------------------------------------- #

def coarse_orientation(
    moving: np.ndarray,
    fixed: np.ndarray,
    grid: GridSpec,
    *,
    angles_deg: Optional[Sequence[float]] = None,
    try_mirror: bool = True,
    coarse_pixel_um: float = 25.0,
) -> list[Transform]:
    """Candidate coarse transforms, ranked by correlation after a translation search.

    Exhaustive over rotation because there is no cheap way to recover a large rotation from
    a translation-only correlation, and serial sections routinely differ by tens of degrees.
    Each hypothesis gets its optimal translation from phase correlation, which is exact and
    fast, so the search is only over angle (and the mirror flag).

    Runs at ``coarse_pixel_um`` rather than the registration grid's own resolution. A real
    section is 10-20 mm across, so at the 2 um working resolution the raster is 5000-10000 px
    on a side and rotating it ~144 times (72 angles x 2 mirror hypotheses) is intractable.
    Orientation is a low-frequency property — a 25 um pixel resolves it perfectly well — and
    because the transform is expressed in **microns** it drops straight into the fine stage at
    full resolution with no conversion. That is precisely why ``Transform`` is defined in
    micron space rather than pixels.

    Returns several candidates rather than one: the best-correlating angle is not reliably
    the right one on low-contrast masks, so the fine stage refines the top few and picks by
    its own metric.
    """
    from scipy import ndimage
    from skimage.registration import phase_cross_correlation

    angles = list(angles_deg) if angles_deg is not None else list(np.arange(0, 360, 5.0))

    factor = max(1, int(round(coarse_pixel_um / grid.pixel_um)))
    if factor > 1:
        before = fixed.shape
        moving = ndimage.zoom(moving, 1.0 / factor, order=1)
        fixed = ndimage.zoom(fixed, 1.0 / factor, order=1)
        # Take the pixel size the zoom actually achieved, not the one it was asked for.
        # ndimage.zoom rounds the output shape, so at the defaults (2 um grid, 25 um coarse)
        # `pixel_um * 12` asserts 24.000 while the real pixel is 24.016 -- a 0.07% scale error,
        # ~10 um across a section, injected into every coarse candidate before the fine stage
        # ever sees it.
        achieved = 0.5 * (grid.pixel_um * before[0] / fixed.shape[0]
                          + grid.pixel_um * before[1] / fixed.shape[1])
        grid = GridSpec(grid.origin_x_um, grid.origin_y_um, achieved,
                        fixed.shape[0], fixed.shape[1])

    # ndimage.rotate(reshape=False) rotates about (shape - 1) / 2, not shape / 2. Half a
    # coarse pixel is 12 um of systematic offset in every candidate this function returns.
    cx = (grid.width - 1) / 2.0
    cy = (grid.height - 1) / 2.0
    fixed_n = normalize_image(fixed)

    scored: list[tuple[float, Transform]] = []
    for mirror in ([False, True] if try_mirror else [False]):
        base = moving[:, ::-1] if mirror else moving
        for ang in angles:
            rot = ndimage.rotate(base, ang, reshape=False, order=1, mode="constant", cval=0.0)
            rot_n = normalize_image(rot)
            try:
                shift, _err, _phase = phase_cross_correlation(
                    fixed_n, rot_n, upsample_factor=4, normalization=None)
            except TypeError:      # older scikit-image signature
                shift, _err, _phase = phase_cross_correlation(fixed_n, rot_n, upsample_factor=4)
            shifted = ndimage.shift(rot_n, shift, order=1, mode="constant", cval=0.0)

            denom = np.linalg.norm(shifted) * np.linalg.norm(fixed_n)
            corr = float((shifted * fixed_n).sum() / denom) if denom > 0 else 0.0

            # Build the equivalent micron transform. scipy rotates about the array centre and
            # `shift` is (row, col); both conversions happen here so nothing downstream has to
            # know about either convention.
            th = np.radians(ang)
            c, s = np.cos(th), np.sin(th)
            lin = np.array([[c, s], [-s, c]])         # ndimage.rotate is clockwise-positive
            if mirror:
                lin = lin @ np.array([[-1.0, 0.0], [0.0, 1.0]])
            centre = np.array([cx, cy])
            t_px = centre - lin @ centre + np.array([shift[1], shift[0]])
            tf = affine_from_pixels(np.column_stack([lin, t_px]), grid)
            tf.meta.update({"coarse_angle_deg": ang, "coarse_mirror": mirror,
                            "coarse_corr": corr})
            scored.append((corr, tf))

    scored.sort(key=lambda t: -t[0])
    return [tf for _corr, tf in scored]


def refine_ecc(
    moving: np.ndarray,
    fixed: np.ndarray,
    grid: GridSpec,
    initial: Transform,
    *,
    motion: str = "affine",
    iterations: int = 300,
    eps: float = 1e-6,
    max_pixels: int = 4_000_000,
    moving_mask: Optional[np.ndarray] = None,
) -> tuple[Transform, float]:
    """Refine with opencv's ECC (enhanced correlation coefficient).

    ECC optimises correlation directly and is robust to the intensity differences that make
    cross-modal registration hard — a DAPI image and a protein channel do not share a linear
    intensity relationship, but they do share structure. Returns ``(transform, ecc_value)``;
    on failure to converge it returns the input transform and ``nan`` rather than raising, so
    one bad donor does not abort a cohort run.

    ``moving_mask`` marks where the moving image actually has data. It matters here because
    the Xenium capture area is a sub-field of the PhenoCycler section — 39% and 54% of it on
    the two donors registered so far — so an unmasked ECC spends most of its correlation
    budget on fixed-frame tissue that has no counterpart at all. OpenCV takes it in the
    ``inputMask`` slot, documented as "valid values of inputImage", and re-warps it with the
    current estimate at every iteration, so the mask tracks the transform rather than being a
    fixed stencil. Mask the *moving* side only, never the intersection: a rigid warp preserves
    the moving mask's area, whereas an intersection mask lets the optimiser buy correlation by
    shrinking the overlap — which is the same partial-overlap failure one level down.

    **A masked ECC value is not comparable to an unmasked one** — measured at 0.31 against
    0.49 for the same transform on donor 6436 — because rho is normalised over whatever the
    mask admits. Callers that compare values across invocations must mask them all or none.
    """
    import cv2
    from scipy import ndimage

    warp_mode = {"euclidean": cv2.MOTION_EUCLIDEAN, "affine": cv2.MOTION_AFFINE}[motion]
    f = normalize_image(fixed).astype(np.float32)
    m = normalize_image(moving).astype(np.float32)

    # Cap the working resolution. ECC is O(pixels) per iteration and runs hundreds of them,
    # so a full-resolution whole-section raster is not viable. Downsampling is safe because
    # the affine lives in micron space: `grid` is rescaled alongside the images, so the
    # recovered transform is identical whatever resolution it was estimated at.
    if f.size > max_pixels:
        shrink = float(np.sqrt(max_pixels / f.size))
        f = ndimage.zoom(f, shrink, order=1)
        m = ndimage.zoom(m, shrink, order=1)
        grid = GridSpec(grid.origin_x_um, grid.origin_y_um, grid.pixel_um / shrink,
                        f.shape[0], f.shape[1])

    # ECC estimates the fixed->moving warp; we hold moving->fixed, so seed with the inverse.
    warp = affine_to_pixels(initial.invert(), grid).astype(np.float32)
    if motion == "euclidean":
        # Force an exact rotation matrix; ECC rejects a seed that is not one.
        u, _s, vt = np.linalg.svd(warp[:, :2])
        warp[:, :2] = u @ vt

    # Resample the mask from `m`'s achieved shape, not from `shrink`: ndimage.zoom rounds, and
    # a one-pixel disagreement is a hard cv2.error rather than a warning.
    mask_u8 = None
    if moving_mask is not None:
        mm = np.asarray(moving_mask)
        if mm.shape != m.shape:
            mm = ndimage.zoom(mm.astype(np.float32),
                              (m.shape[0] / mm.shape[0], m.shape[1] / mm.shape[1]),
                              order=1)[:m.shape[0], :m.shape[1]]
        mask_u8 = np.ascontiguousarray(((np.asarray(mm) > 0.5).astype(np.uint8) * 255))

    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, iterations, eps)
    try:
        value, warp = cv2.findTransformECC(f, m, warp, warp_mode, criteria, mask_u8, 5)
    except cv2.error:
        return initial, float("nan")

    fixed_to_moving = affine_from_pixels(np.asarray(warp, dtype=np.float64), grid)
    return fixed_to_moving.invert(), float(value)


def refine_bspline(
    moving: np.ndarray,
    fixed: np.ndarray,
    grid: GridSpec,
    *,
    max_disp_um: float = 200.0,
    grid_spacing_um: float = 400.0,
    iterations: int = 200,
    max_pixels: int = 2_000_000,
) -> Optional[np.ndarray]:
    """Non-rigid B-spline refinement; returns a ``(2, H, W)`` micron displacement field.

    Deliberately coarse: the control-point spacing is hundreds of microns, so the field can
    absorb the smooth stretch and shear of sectioning without being able to shrink-wrap one
    section's islets onto another's. The result is then hard-capped at ``max_disp_um``.

    That cap is the important part. A B-spline with enough freedom can align any two images,
    including two that have nothing to do with each other, and a "successful" registration of
    unrelated sections is far more damaging than an obvious failure.
    """
    try:
        import SimpleITK as sitk
    except ImportError:
        print("[register] SimpleITK not installed — skipping the non-rigid stage "
              "(pip install -e '.[integration]')", flush=True)
        return None

    from scipy import ndimage

    fa = normalize_image(fixed).astype(np.float32)
    ma = normalize_image(moving).astype(np.float32)
    if fa.size > max_pixels:
        shrink = float(np.sqrt(max_pixels / fa.size))
        fa = ndimage.zoom(fa, shrink, order=1)
        ma = ndimage.zoom(ma, shrink, order=1)
        grid = GridSpec(grid.origin_x_um, grid.origin_y_um, grid.pixel_um / shrink,
                        fa.shape[0], fa.shape[1])
    f = sitk.GetImageFromArray(fa)
    m = sitk.GetImageFromArray(ma)

    mesh = [max(1, int(round(s * grid.pixel_um / grid_spacing_um)))
            for s in (grid.width, grid.height)]
    try:
        tx = sitk.BSplineTransformInitializer(f, mesh)
        reg = sitk.ImageRegistrationMethod()
        reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=32)
        reg.SetMetricSamplingStrategy(reg.RANDOM)
        reg.SetMetricSamplingPercentage(0.2, seed=1)
        reg.SetInterpolator(sitk.sitkLinear)
        reg.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=iterations,
                                          convergenceMinimumValue=1e-6, convergenceWindowSize=10)
        reg.SetOptimizerScalesFromPhysicalShift()
        reg.SetInitialTransform(tx, inPlace=True)
        reg.SetShrinkFactorsPerLevel([4, 2, 1])
        reg.SetSmoothingSigmasPerLevel([2, 1, 0])
        reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
        out_tx = reg.Execute(f, m)
    except Exception as exc:  # noqa: BLE001 - ITK raises a wide variety on bad input
        print(f"[register] B-spline stage failed ({exc}); keeping the affine result",
              flush=True)
        return None

    disp_filter = sitk.TransformToDisplacementFieldFilter()
    disp_filter.SetReferenceImage(f)
    disp_filter.SetOutputPixelType(sitk.sitkVectorFloat64)
    disp = sitk.GetArrayFromImage(disp_filter.Execute(out_tx))   # (H, W, 2), index units

    # Sign convention. SimpleITK's registration returns the transform used to *resample* the
    # moving image into the fixed frame, i.e. it maps fixed -> moving. `Transform.displacement`
    # is the opposite: moving -> fixed, added after the affine. So the field is negated here.
    # Getting this backwards does not fail loudly — it produces a transform that doubles the
    # misalignment instead of removing it — so it is asserted by tests/test_integration_register.
    # Vector components are (x, y); pixel size converts index units to microns.
    field = -np.stack([disp[..., 0], disp[..., 1]]).astype(np.float32) * grid.pixel_um
    bspline_grid = grid

    mag = np.linalg.norm(field, axis=0)
    over = mag > max_disp_um
    if over.any():
        scale = np.ones_like(mag)
        scale[over] = max_disp_um / mag[over]
        field = field * scale
        print(f"[register] capped {100*over.mean():.1f}% of the displacement field at "
              f"{max_disp_um:.0f} um", flush=True)
    return field, bspline_grid


def displacement_report(tf: Transform, *, n_samples: int = 20_000, seed: int = 0) -> dict:
    """Whether a non-rigid field is a map that can actually be inverted.

    ``max_displacement_um`` cannot answer this, and it was the only thing recorded. A field can
    sit inside its cap everywhere and still fold -- ``x -> x + D(x)`` stops being injective the
    moment ``det(I + grad D)`` goes non-positive, and where it has folded there is no inverse to
    find, so ``Transform.inverse_apply``'s fixed point does not merely converge slowly, it does
    not exist. Every annotation vertex landing in such a region is lost.

    The cap is also implicated in causing it: ``refine_bspline`` clips magnitude **per pixel**,
    which is a non-smooth operation applied to whatever fraction of the field is over the
    bound, and it introduces gradient discontinuities exactly at the cap boundary. Measured on
    donor 6436's first field, 14.7% of pixels were pinned, 0.85% of the field was folded, and
    the round trip missed by up to 400 um.

    Reported rather than assumed, because it is pure geometry with no biology in it and so is
    one of the few things here that is safe to gate on hard.
    """
    if tf.displacement is None or tf.grid is None:
        return {"has_displacement": False}

    g = tf.grid
    d = np.asarray(tf.displacement, dtype=np.float64)
    mag = np.linalg.norm(d, axis=0)
    dxdy, dxdx = np.gradient(d[0], g.pixel_um)
    dydy, dydx = np.gradient(d[1], g.pixel_um)
    det = (1.0 + dxdx) * (1.0 + dydy) - dxdy * dydx

    rng = np.random.default_rng(seed)
    pts = np.column_stack([
        rng.uniform(g.origin_x_um, g.origin_x_um + g.width * g.pixel_um, n_samples),
        rng.uniform(g.origin_y_um, g.origin_y_um + g.height * g.pixel_um, n_samples)])
    rt = np.linalg.norm(np.asarray(tf.apply(tf.inverse_apply(pts, n_iter=10))) - pts, axis=1)

    return {
        "has_displacement": True,
        "max_displacement_um": float(mag.max()),
        "fraction_at_cap": float(np.mean(mag >= 0.999 * mag.max())),
        "jacobian_det_min": float(det.min()),
        "fraction_folded": float(np.mean(det <= 0.0)),
        "roundtrip_p99_um": float(np.percentile(rt, 99)),
        "roundtrip_max_um": float(rt.max()),
    }


def displacement_is_invertible(
    tf: Transform,
    *,
    min_jacobian_det: float = 0.2,
    max_roundtrip_um: float = 1.0,
) -> tuple[bool, str]:
    """Whether the non-rigid part can be inverted well enough to carry geometry back."""
    report = displacement_report(tf)
    if not report.get("has_displacement"):
        return True, "affine only"
    if report["jacobian_det_min"] < min_jacobian_det:
        return False, (
            f"displacement folds: min Jacobian det {report['jacobian_det_min']:.3f} "
            f"below {min_jacobian_det:.2f}, {100 * report['fraction_folded']:.2f}% of the "
            f"field non-injective"
        )
    if report["roundtrip_max_um"] > max_roundtrip_um:
        return False, (
            f"displacement does not round trip: max {report['roundtrip_max_um']:.2f} um "
            f"over the {max_roundtrip_um:.2f} um bound "
            f"(p99 {report['roundtrip_p99_um']:.2f} um)"
        )
    return True, (
        f"displacement invertible (det min {report['jacobian_det_min']:.3f}, "
        f"round trip max {report['roundtrip_max_um']:.3f} um)"
    )


def scale_report(tf: Transform) -> dict:
    """Decompose a transform's linear part into the quantities that must stay physical.

    ``Transform.scale`` reports one isotropic number (the square root of |det|), which hides
    exactly the failure this lane is prone to: a stretch along one axis only. Donor 6436's
    affine came back at scale 1.0509 -- unremarkable on its own -- while actually holding
    sx=1.1057 against sy=0.9989, an 11% x-only stretch worth 509 um of drift across a 10 mm
    section. So report both axes, their ratio, and the shear.
    """
    a = np.asarray(tf.affine, dtype=np.float64)[:, :2]
    sx = float(np.linalg.norm(a[:, 0]))
    sy = float(np.linalg.norm(a[:, 1]))
    det = float(np.linalg.det(a))
    lo, hi = (sx, sy) if sx <= sy else (sy, sx)
    cosine = float(a[:, 0] @ a[:, 1] / (sx * sy)) if sx > 0 and sy > 0 else 0.0
    return {
        "scale_x": sx,
        "scale_y": sy,
        "anisotropy": (hi / lo) if lo > 0 else np.inf,
        "linear_scale": float(np.sqrt(abs(det))),
        "shear_cosine": cosine,
        "shear_deg": float(90.0 - np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))),
    }


def scale_is_physical(
    tf: Transform,
    *,
    max_anisotropy: float = 1.03,
    max_scale_deviation: float = 0.05,
) -> tuple[bool, str]:
    """Whether a transform's scale is one two calibrated micron frames could actually produce.

    Serial sections do stretch, but they stretch a little and roughly isotropically. A large
    single-axis scale is a fitting artefact, not tissue, and it is the one error that is
    invisible in a whole-section overlay yet ruins every downstream micron measurement.
    """
    report = scale_report(tf)
    if not np.isfinite(report["anisotropy"]) or report["anisotropy"] > max_anisotropy:
        return False, (
            f"anisotropy {report['anisotropy']:.4f} exceeds {max_anisotropy:.4f} "
            f"(sx={report['scale_x']:.4f}, sy={report['scale_y']:.4f})"
        )
    deviation = abs(report["linear_scale"] - 1.0)
    if deviation > max_scale_deviation:
        return False, (
            f"linear scale {report['linear_scale']:.4f} departs from 1 by {deviation:.4f}, "
            f"over the {max_scale_deviation:.4f} bound"
        )
    return True, (
        f"scale ok (sx={report['scale_x']:.4f}, sy={report['scale_y']:.4f}, "
        f"anisotropy {report['anisotropy']:.4f}, linear {report['linear_scale']:.4f})"
    )


def select_coarse_candidates(
    candidates: Sequence[Transform],
    *,
    keep_frac: float = 0.97,
    min_candidates: int = 8,
    max_candidates: int = 24,
    max_rotation_deg: float = 45.0,
) -> list[Transform]:
    """Which coarse hypotheses are worth refining.

    A fixed rank cut was the single defect that made donor 6436's registration wrong. Its
    correct hypothesis ranked 4th of 144 with a coarse score of 0.4352 against a best of
    0.4433 -- a 1.8% spread across four hypotheses 175 deg apart -- and the cut was 3.

    That bunching is structural rather than bad luck. The coarse score is a whole-frame
    cosine, and when the moving section covers only a fraction of the fixed one (39% and 54%
    on the two donors registered so far), most placements overlap a similar amount of tissue
    whichever way the moving image is turned. The ranking still carries real signal -- it put
    the truth 1st for one donor and 4th for the other -- but its *margins* are small, so the
    honest rule is a plateau, not a rank: refine everything within ``keep_frac`` of the best,
    with a floor so a sharp winner still gets company and a ceiling so a flat landscape does
    not run ECC 144 times.
    """
    if not candidates:
        return []

    # Drop hypotheses the sections cannot physically be in. The maintainer records that no
    # section is mounted past 45 degrees from its partner (2026-09-11), corroborated by their
    # own capture-box inspection across 17 donors. This bound is what makes widening the
    # candidate set safe: the coarse score is nearly angle-blind under partial overlap -- a
    # compact region inside a larger tissue mask overlaps it about equally whichever way it is
    # turned -- so without it, refining more candidates simply gives a masked ECC more wrong
    # angles to polish into something locally plausible. Measured across the cohort it did
    # exactly that, returning -100, -99, +89, -73 and +58 degrees on near-upright sections.
    if np.isfinite(max_rotation_deg):
        within = [c for c in candidates
                  if abs((float(c.rotation_deg) + 180.0) % 360.0 - 180.0) <= max_rotation_deg]
        if within:
            candidates = within
    scores = np.array([float(c.meta.get("coarse_corr", 0.0)) for c in candidates])
    best = float(scores[0])
    n = int((scores >= keep_frac * best).sum()) if best > 0 else min_candidates
    n = int(np.clip(n, min(min_candidates, len(candidates)), max_candidates))
    return list(candidates[:n])


def register_image(
    moving: np.ndarray,
    fixed: np.ndarray,
    grid: GridSpec,
    *,
    nonrigid: bool = True,
    max_disp_um: float = 200.0,
    coarse_step_deg: float = 5.0,
    coarse_keep_frac: float = 0.97,
    coarse_min_candidates: int = 8,
    coarse_max_candidates: int = 24,
    max_rotation_deg: float = 45.0,
    max_anisotropy: float = 1.03,
    max_scale_deviation: float = 0.05,
    moving_mask: Optional[np.ndarray] = None,
) -> RegistrationResult:
    """Track B: coarse orientation -> ECC euclidean -> ECC affine -> optional B-spline."""
    res = RegistrationResult(transform=Transform(), track="image", grid=grid)

    candidates = coarse_orientation(
        moving, fixed, grid, angles_deg=np.arange(0, 360, coarse_step_deg))
    if not candidates:
        res.notes.append("coarse orientation search produced no candidates")
        return res
    res.stages.append("coarse")

    refine_set = select_coarse_candidates(
        candidates, keep_frac=coarse_keep_frac, min_candidates=coarse_min_candidates,
        max_candidates=coarse_max_candidates, max_rotation_deg=max_rotation_deg)
    res.notes.append(f"coarse: rotation constrained to +/-{max_rotation_deg:.0f} deg "
                     f"(recorded mounting bound)")
    res.notes.append(
        f"coarse: refined {len(refine_set)} of {len(candidates)} hypotheses "
        f"(corr {float(candidates[0].meta.get('coarse_corr', float('nan'))):.4f} .. "
        f"{float(refine_set[-1].meta.get('coarse_corr', float('nan'))):.4f})")
    if moving_mask is not None:
        res.stages.append("masked_ecc")

    # Refine a RESIDUAL, never the whole transform. ECC linearises around its current estimate,
    # so a seed far from identity is not merely slow to converge -- it diverges outright. Donor
    # 6450's sections differ by ~170 deg and every one of its eight candidates failed with
    # "NaN encountered", masked and unmasked alike; 6436 and 6442 only ever worked because
    # their coarse winners happened to sit at -20 and +5 deg.
    #
    # So the candidate is applied to the moving image first and ECC is handed a near-identity
    # problem, then the two are composed back in micron space. This costs one resampling of the
    # moving raster per candidate and it is also what makes the mirror search real: a reflected
    # seed cannot be expressed as a MOTION_EUCLIDEAN warp at all -- cv2 rejects the non-proper
    # rotation the SVD returns -- so before this, every mirrored candidate raised cv2.error and
    # was silently dropped. Pre-warping puts the reflection in the image, where it costs
    # nothing. On 6450 all eight candidates converge, mirrored ones included.
    best, best_val = None, -np.inf
    for cand in refine_set:
        pre = cand.warp_mask(np.asarray(moving, dtype=np.float32), grid)
        pre_mask = (cand.warp_mask(np.asarray(moving_mask, dtype=bool), grid)
                    if moving_mask is not None else None)
        residual, val = refine_ecc(pre, fixed, grid, Transform(), motion="euclidean",
                                   moving_mask=pre_mask)
        if not np.isfinite(val):
            continue
        if val > best_val:
            best, best_val = residual.compose(cand), val
    if best is None:
        # Every ECC seed diverged; keep the best correlating coarse candidate so the caller
        # still gets a usable starting point rather than the identity.
        res.transform = candidates[0]
        res.mirrored = candidates[0].is_mirrored
        res.notes.append("ECC did not converge from any coarse seed; using coarse only")
        return res
    res.stages.append("ecc_euclidean")

    # The affine stage is the only one that can change scale, and a six-parameter fit driven by
    # a partial-overlap silhouette will happily buy correlation with an unphysical stretch. Take
    # it only when it is both better AND a scale two calibrated micron frames could produce;
    # otherwise keep the rigid fit, which has scale exactly 1 by construction.
    # Same treatment for the affine stage, and for the same reason: it is seeded with `best`,
    # which after a large-rotation coarse win is itself far from identity.
    pre_best = best.warp_mask(np.asarray(moving, dtype=np.float32), grid)
    pre_best_mask = (best.warp_mask(np.asarray(moving_mask, dtype=bool), grid)
                     if moving_mask is not None else None)
    residual, val = refine_ecc(pre_best, fixed, grid, Transform(), motion="affine",
                               moving_mask=pre_best_mask)
    tf = residual.compose(best)
    if np.isfinite(val) and val >= best_val:
        physical, why = scale_is_physical(
            tf, max_anisotropy=max_anisotropy, max_scale_deviation=max_scale_deviation)
        if physical:
            best, best_val = tf, val
            res.stages.append("ecc_affine")
            res.notes.append(f"affine accepted: {why}")
        else:
            res.notes.append(f"affine refused, kept rigid: {why}")

    res.transform = best
    res.mirrored = best.is_mirrored
    res.notes.append(f"ecc={best_val:.4f}")

    if nonrigid:
        warped = best.warp_mask(moving, grid)
        out = refine_bspline(warped, fixed, grid, max_disp_um=max_disp_um)
        if out is not None:
            field, field_grid = out
            res.transform = Transform(affine=best.affine, displacement=field,
                                      grid=field_grid, meta=dict(best.meta))
            res.stages.append("bspline")
    return res


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #

def _load_structures(cfg: PipelineConfig, modality: str, donor: str, roi: str,
                     kind: str = "islet") -> Optional[pd.DataFrame]:
    path = structures_path(cfg, modality, donor, roi, kind)
    if not path.exists():
        return None
    return pd.read_parquet(path)


def footprint_ncc(
    moving_img: np.ndarray,
    fixed_img: np.ndarray,
    moving_mask: np.ndarray,
    grid: GridSpec,
    tf: Transform,
    warped_mask: Optional[np.ndarray] = None,
) -> tuple[float, float]:
    """Correlation between the two images *inside the warped moving footprint only*.

    Returns ``(r, n_pixels)``, where ``n_pixels`` is how many grid pixels the correlation was
    actually measured over, so a correlation computed on almost nothing is visible as such
    rather than looking like a good score. Containment -- the share of the footprint landing on
    fixed tissue -- is a different quantity and is owned by :func:`evaluate`, which has the
    fixed mask to compute it against.

    Restricting to the footprint is what makes this comparable across placements: an
    unrestricted correlation over a canvas that is mostly background is dominated by
    background agreeing with background, which every wrong placement also achieves.
    """
    warped = tf.warp_mask(np.asarray(moving_img, dtype=np.float32), grid)
    mask = (tf.warp_mask(np.asarray(moving_mask, dtype=bool), grid)
            if warped_mask is None else np.asarray(warped_mask, dtype=bool))
    if mask.sum() < 100:
        return float("nan"), 0
    a = np.asarray(warped, dtype=np.float64)[mask]
    b = np.asarray(fixed_img, dtype=np.float64)[mask]
    if a.std() <= 0 or b.std() <= 0:
        return float("nan"), int(mask.sum())
    r = float(np.corrcoef(a, b)[0, 1])
    return r, int(mask.sum())


def footprint_ncc_null(
    moving_img: np.ndarray,
    fixed_img: np.ndarray,
    moving_mask: np.ndarray,
    fixed_mask: np.ndarray,
    grid: GridSpec,
    tf: Transform,
    *,
    n: int = 120,
    max_angle_deg: float = 180.0,
    max_shift_um: float = 3000.0,
    min_containment: float = 0.85,
    min_floor: float = 0.5,
    seed: int = 0,
) -> tuple[np.ndarray, float]:
    """The same correlation at ``n`` other placements of the same footprint.

    A **within-donor, within-section** null: staining, cell density, tissue quality and the
    donor's disease state are all held fixed, and only where the footprint sits varies. That
    matters because an absolute correlation cannot be gated on here -- donor 6436's correct
    placement reads higher than donor 6442's correct placement, so a bar set from one would
    reject the other. What generalises is how far a placement stands above the alternatives
    available on its own section.

    Placements are constrained to keep at least ``min_containment`` of the footprint on fixed
    tissue, so the null is drawn from placements that are *plausible*, not from ones that hang
    the footprint off the slide and trivially score nothing.

    That floor is **relaxed until enough placements exist**, down to ``min_floor``, and the
    floor actually reached is returned alongside the samples. A fixed floor is not safe: donor
    6450's own placement sits at 0.823 containment and almost nothing else on its section
    reaches even 0.78, so a fixed 0.85 produced an empty null and silently removed the primary
    bar for that donor. Relaxing is not free and is not the safe direction -- a less contained
    placement puts more of the footprint over background, where the moving raster has signal
    and the fixed one has none, which drives the null's correlation *down* and makes the
    observation easier to beat. So the floor reached is reported and travels with the z.
    """
    rng = np.random.default_rng(seed)
    # Rotate about the WARPED FOOTPRINT'S OWN CENTROID, not the grid centre. Rotating a
    # footprint that does not sit at the centre of the canvas about that centre couples every
    # rotation to a translation the size of its offset, which throws it off the section --
    # donor 6450 produced zero usable placements out of 2,400 attempts that way, silently
    # removing its primary bar. The null is supposed to vary orientation and position while
    # holding coverage roughly fixed; turning the footprint in place is what does that.
    rows0, cols0 = np.nonzero(np.asarray(moving_mask, dtype=bool))
    if len(rows0) == 0:
        return np.zeros(0)
    seed_centre = grid.to_microns(
        np.array([[float(cols0.mean()), float(rows0.mean())]]))[0]
    warped_centre = tf.apply(seed_centre.reshape(1, 2))[0]
    centre = (float(warped_centre[0]), float(warped_centre[1]))

    # Probe points for the containment test: a subsample of the moving footprint's own pixels,
    # carried forward and looked up in the fixed mask. Warping the whole raster just to reject
    # a placement costs 160 ms and most placements are rejected, so the probe is what makes a
    # null of this size affordable at all.
    rows, cols = rows0, cols0
    take = rng.choice(len(rows), size=min(4000, len(rows)), replace=False)
    probe_um = grid.to_microns(np.column_stack([cols[take], rows[take]]).astype(np.float64))

    def _containment(cand: Transform) -> float:
        px = grid.to_pixels(cand.apply(probe_um))
        c = np.rint(px[:, 0]).astype(int)
        r = np.rint(px[:, 1]).astype(int)
        ok = (c >= 0) & (c < grid.width) & (r >= 0) & (r < grid.height)
        if not ok.any():
            return 0.0
        hit = np.zeros(len(px), dtype=bool)
        hit[ok] = fixed_mask[r[ok], c[ok]]
        return float(hit.mean())

    # Draw the candidate placements once, then decide how far the floor has to drop to admit
    # enough of them. Sampling once keeps the relaxation honest: the same placements are being
    # judged at every floor, so a lower floor only ever adds, never reshuffles.
    proposals = []
    for _ in range(n * 20):
        cand = from_rigid(float(rng.uniform(-max_angle_deg, max_angle_deg)),
                          (float(rng.uniform(-max_shift_um, max_shift_um)),
                           float(rng.uniform(-max_shift_um, max_shift_um))),
                          center_um=centre).compose(tf)
        proposals.append((_containment(cand), cand))
    proposals.sort(key=lambda t: -t[0])

    floor = float(min_containment)
    for trial in (min_containment, 0.9 * min_containment, 0.8 * min_containment,
                  0.7 * min_containment, min_floor):
        trial = max(float(trial), float(min_floor))
        if sum(1 for c, _ in proposals if c >= trial) >= n:
            floor = trial
            break
        floor = trial

    out: list[float] = []
    for containment, cand in proposals:
        if len(out) >= n:
            break
        if containment < floor:
            break
        warped = cand.warp_mask(np.asarray(moving_mask, dtype=bool), grid)
        if warped.sum() < 100:
            continue
        r, _c = footprint_ncc(moving_img, fixed_img, moving_mask, grid, cand,
                              warped_mask=warped)
        if np.isfinite(r):
            out.append(r)
    return np.asarray(out, dtype=np.float64), floor


def evaluate(
    result: RegistrationResult,
    moving_table: CellTable,
    fixed_table: CellTable,
    grid: GridSpec,
    moving_structures: Optional[pd.DataFrame] = None,
    fixed_structures: Optional[pd.DataFrame] = None,
    *,
    qc_pixel_um: float = 10.0,
    moving_image: Optional[np.ndarray] = None,
    fixed_image: Optional[np.ndarray] = None,
    null_placements: int = 120,
) -> RegistrationResult:
    """Fill in the QC block for a candidate transform.

    Deliberately measured on quantities the transform was *not* optimised against wherever
    possible: the point-set track never sees the tissue mask, and the image track never sees
    islet centroids, so each track's headline metric is an out-of-sample check on the other.

    Evaluated on a deliberately coarse grid. Dice between two tissue outlines is a
    low-frequency measurement — at 10 um a whole-section mask is a few hundred pixels on a
    side, the metric is unchanged to three decimals, and the warp is two orders of magnitude
    cheaper than at the 2 um working resolution. Cost matters because this runs once per
    track per donor.
    """
    if qc_pixel_um > grid.pixel_um:
        factor = qc_pixel_um / grid.pixel_um
        grid = GridSpec(grid.origin_x_um, grid.origin_y_um, qc_pixel_um,
                        max(1, int(round(grid.height / factor))),
                        max(1, int(round(grid.width / factor))))

    mv_mask = tissue_mask(moving_table.xy(), grid)
    fx_mask = tissue_mask(fixed_table.xy(), grid)
    warped = result.transform.warp_mask(mv_mask, grid)
    result.tissue_dice = dice(warped, fx_mask)
    result.footprint_containment = float(
        np.logical_and(warped, fx_mask).sum() / max(warped.sum(), 1))

    # Image structure inside the moving footprint, against other placements of that same
    # footprint on this same section. The images arrive on the registration grid, so they are
    # coarsened here the same way the masks were.
    if moving_image is not None and fixed_image is not None:
        from scipy import ndimage

        def _to_qc(img: np.ndarray) -> np.ndarray:
            a = np.asarray(img, dtype=np.float32)
            if a.shape == grid.shape:
                return a
            return ndimage.zoom(a, (grid.height / a.shape[0], grid.width / a.shape[1]),
                                order=1)[:grid.height, :grid.width]

        mv_i, fx_i = _to_qc(moving_image), _to_qc(fixed_image)
        if mv_i.shape != grid.shape or fx_i.shape != grid.shape:
            # Never skip the primary bar in silence: without a note here the verdict would
            # read "no image pair was registered" for a donor that had one.
            result.notes.append(
                f"footprint NCC not computed: rasters coarsened to {mv_i.shape} / "
                f"{fx_i.shape} against a QC grid of {grid.shape}")
        else:
            r, n_px = footprint_ncc(mv_i, fx_i, mv_mask, grid, result.transform)
            result.footprint_ncc = r
            # Draw the null from placements at least as well contained as THIS one, not from an
            # absolute standard the donor may not reach. Donor 6450 sits at 0.822 containment,
            # below the 0.85 floor, so every sampled placement was rejected and the null came
            # back empty -- which silently removed the primary bar for that donor. The null's
            # job is to vary placement while holding coverage roughly fixed; tying the floor to
            # the observation is what actually does that.
            floor = min(0.85, 0.95 * result.footprint_containment) \
                if np.isfinite(result.footprint_containment) else 0.85
            null, floor = footprint_ncc_null(
                mv_i, fx_i, mv_mask, fx_mask, grid, result.transform,
                n=null_placements, min_containment=floor)
            result.footprint_null_n = int(len(null))
            result.footprint_null_floor = float(floor)
            if len(null):
                result.footprint_null_mean = float(null.mean())
                result.footprint_null_sd = float(null.std())
            if np.isfinite(r) and len(null) >= 20 and null.std() > 0:
                result.footprint_ncc_z = float((r - null.mean()) / null.std())
                result.notes.append(
                    f"footprint NCC {r:.3f} over {n_px:,} px vs null "
                    f"{null.mean():.3f}+/-{null.std():.3f} over {len(null)} other placements "
                    f"at containment >= {floor:.2f} -> z={result.footprint_ncc_z:.1f}")
            elif np.isfinite(r):
                result.notes.append(
                    f"footprint NCC {r:.3f} over {n_px:,} px; the null produced only "
                    f"{len(null)} usable placements at containment >= {floor:.2f}, so no z is "
                    f"reported and the primary bar is UNAVAILABLE for this donor")

    if moving_structures is not None and fixed_structures is not None \
            and len(moving_structures) and len(fixed_structures):
        mv = moving_structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
        fx = fixed_structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
        rmse, n = nearest_neighbour_rmse(result.transform.apply(mv), fx, max_dist_um=300.0)
        if result.track == "image" or not np.isfinite(result.rmse_um):
            result.rmse_um = rmse
            result.n_inliers = n
            result.n_candidates = result.n_candidates or min(len(mv), len(fx))

        # Islet Dice: dilate each centroid to a ~75 um disc so the metric measures "did the
        # islets land on each other" rather than "did two zero-area points coincide", which
        # would be 0 for any transform. `iterations` on a 3x3 element rather than one large
        # structuring element — same result, far cheaper, but only with the 8-connected
        # element: scipy's default is the 4-connected cross, which iterates to an L1 diamond
        # whose diagonal reach is 1/sqrt(2) of its axial reach, so the "75 um disc" was
        # really 75 um across the axes and 53 um across the diagonals.
        from scipy import ndimage

        disc = ndimage.generate_binary_structure(2, 2)
        r = max(1, int(round(75.0 / grid.pixel_um)))
        mv_islets = ndimage.binary_dilation(rasterize_structures(moving_structures, grid) > 0,
                                            structure=disc, iterations=r)
        fx_islets = ndimage.binary_dilation(rasterize_structures(fixed_structures, grid) > 0,
                                            structure=disc,
                                            iterations=r)
        result.islet_dice = dice(result.transform.warp_mask(mv_islets, grid), fx_islets)

        # Does the joint rotation+translation correspondence landscape peak *here*? An
        # absolute islet Dice or RMSE cannot answer that -- both saturate, so a transform in
        # the wrong basin returns the null value rather than a bad one, which is
        # indistinguishable from a donor whose sections genuinely do not correspond.
        vote = residual_vote(result.transform.apply(mv), fx)
        result.islet_peak_ratio = vote["contrast"]
        result.islet_peak_angle_deg = vote["angle_deg"]
        result.islet_peak_offset_um = float(np.hypot(vote["dx_um"], vote["dy_um"]))
    return result


def register_donor(
    cfg: PipelineConfig,
    donor: str,
    roi: str = "panc",
    *,
    use_images: bool = True,
    nonrigid: Optional[bool] = None,
) -> Optional[RegistrationResult]:
    """Register one donor's two sections, running both tracks and keeping the better one."""
    nonrigid = cfg.reg_nonrigid if nonrigid is None else nonrigid

    fixed_is_pheno = cfg.fixed_modality == "phenocycler"
    fixed_root = cfg.cells_pheno_dir if fixed_is_pheno else cfg.cells_xen_dir
    moving_root = cfg.cells_xen_dir if fixed_is_pheno else cfg.cells_pheno_dir
    fixed_mod = "phenocycler" if fixed_is_pheno else "xenium"
    moving_mod = "xenium" if fixed_is_pheno else "phenocycler"

    try:
        fixed_table = read_cell_table(fixed_root, donor, roi, modality=fixed_mod)
        moving_table = read_cell_table(moving_root, donor, roi, modality=moving_mod)
    except FileNotFoundError as exc:
        print(f"[register] donor {donor}: {exc}", flush=True)
        return None

    grid = grid_for(fixed_table.xy(), cfg.reg_pixel_um, other=moving_table.xy())

    fx_struct = _load_structures(cfg, fixed_mod, donor, roi)
    mv_struct = _load_structures(cfg, moving_mod, donor, roi)

    results: list[RegistrationResult] = []

    # Track A — point set.
    ps = register_point_set(mv_struct, fx_struct,
                            inlier_tol_um=cfg.match_max_dist_um * 0.75)
    if ps.transform is not None and (ps.n_inliers or ps.stages):
        results.append(evaluate(ps, moving_table, fixed_table, grid, mv_struct, fx_struct))

    # Track B — image/raster.
    if use_images:
        mv_img, fx_img = _registration_rasters(cfg, moving_table, fixed_table, grid,
                                               moving_mod, fixed_mod, donor)
        # Where the moving modality actually looked. The same primitive `evaluate` scores
        # against below, so the mask ECC optimises under and the mask QC measures against
        # cannot drift apart.
        mv_mask = tissue_mask(moving_table.xy(), grid)
        img = register_image(mv_img, fx_img, grid, nonrigid=nonrigid,
                             max_disp_um=cfg.reg_max_disp_um,
                             max_anisotropy=cfg.reg_max_anisotropy,
                             max_scale_deviation=cfg.reg_max_scale_deviation,
                             coarse_keep_frac=cfg.reg_coarse_keep_frac,
                             coarse_min_candidates=cfg.reg_coarse_min_candidates,
                             coarse_max_candidates=cfg.reg_coarse_max_candidates,
                             max_rotation_deg=cfg.reg_max_rotation_deg,
                             moving_mask=mv_mask)
        results.append(evaluate(img, moving_table, fixed_table, grid, mv_struct, fx_struct,
                                moving_image=mv_img, fixed_image=fx_img))

    if not results:
        print(f"[register] donor {donor}: no track produced a transform", flush=True)
        return None

    best = max(results, key=lambda r: r.score())
    best.grid = grid
    if len(results) > 1:
        other = [r for r in results if r is not best][0]
        best.notes.append(
            f"chosen over {other.track} (score {best.score():.3f} vs {other.score():.3f})")
    return best


def _registration_rasters(
    cfg: PipelineConfig,
    moving_table: CellTable,
    fixed_table: CellTable,
    grid: GridSpec,
    moving_mod: str,
    fixed_mod: str,
    donor: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Best available raster pair: DAPI-to-DAPI if the images are mounted, else densities.

    Falls back rather than failing. Most donors in a cohort will not have both morphology
    trees mounted at the same time, and a density-based registration is a genuinely useful
    result, not a placeholder.
    """
    from .rasterize import image_to_grid, read_qptiff_channel, read_xenium_morphology

    images: dict[str, Optional[np.ndarray]] = {"phenocycler": None, "xenium": None}

    # Resolve the qptiff for THIS SECTION, not just this donor: a donor's pancreas and its
    # lymph nodes are separate images, and for a lymph node this raster is the only
    # registration track available (the point-set track needs islets to seed on).
    section_roi = getattr(moving_table, "roi", "") or getattr(fixed_table, "roi", "")
    qptiff_path, _ = _pheno_paths(cfg, donor, section_roi)
    qptiff = Path(qptiff_path) if qptiff_path else None
    if qptiff is not None:
        try:
            img, um = read_qptiff_channel(qptiff, channel=0, target_um=cfg.reg_pixel_um)
            images["phenocycler"] = image_to_grid(img, um, grid)
            print(f"[register] {donor}/{section_roi or '?'}: PhenoCycler DAPI from "
                  f"{qptiff.name} ({um:.2f} um/px)", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[register] {donor}/{section_roi or '?'}: could not read {qptiff.name}: "
                  f"{exc}", flush=True)
    elif section_roi:
        print(f"[register] {donor}/{section_roi}: no qptiff for this section under "
              f"{cfg.images_dir} — falling back to a cell-density raster for the PhenoCycler "
              f"side. Another section's image is deliberately NOT substituted: it would "
              f"align this section against a different piece of tissue.", flush=True)

    man = load_manifest(cfg)
    row = man[(man["donor_id"] == donor) & (man["roi"] == moving_table.roi)]
    bundle = str(row.iloc[0]["xenium_bundle"]) if len(row) else ""
    if bundle and Path(bundle).exists():
        try:
            img, um = read_xenium_morphology(Path(bundle), channel=0,
                                             target_um=cfg.reg_pixel_um)
            images["xenium"] = image_to_grid(img, um, grid)
            print(f"[register] donor {donor}: Xenium DAPI from {Path(bundle).name} "
                  f"({um:.2f} um/px)", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[register] donor {donor}: could not read Xenium morphology: {exc}",
                  flush=True)

    def fallback(table: CellTable) -> np.ndarray:
        if has_identity(table.df):
            stack = lineage_density(table, grid, DENSITY_LINEAGES)
            return normalize_image(stack.sum(axis=0))
        # True geometry fallback: all-cell density. A geometry profile cannot select
        # Endocrine/Vascular/Exocrine without inventing identity, but its complete point field
        # still supplies a tissue-density image for coarse/affine registration.
        return normalize_image(rasterize_points(table.xy(), grid))

    mv = images[moving_mod] if images[moving_mod] is not None else fallback(moving_table)
    fx = images[fixed_mod] if images[fixed_mod] is not None else fallback(fixed_table)
    if images[moving_mod] is None or images[fixed_mod] is None:
        print(f"[register] donor {donor}: registering on cell-density rasters "
              f"(DAPI unavailable for at least one modality)", flush=True)
    return normalize_image(mv), normalize_image(fx)


def registration_dir(cfg: PipelineConfig, donor: str, roi: str) -> Path:
    base = cfg.registration_dir / f"donor_id={donor}"
    if roi:
        base = base / f"roi={roi}"
    return base


def run_register(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    roi: Optional[str] = None,
    tissue: Optional[str] = None,
    use_images: bool = True,
) -> pd.DataFrame:
    """Stage entry point: register every paired section and persist the transforms.

    One call covers the whole selection. Registration is per (donor, roi) — a donor's
    pancreas and each of its lymph nodes are separate sections with separate transforms.
    """
    # Geometry only: coordinates and rasters. Registration does not read a lineage label.
    import json

    man = paired_rows(load_manifest(cfg), roi=roi, tissue=tissue)
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    if man.empty:
        print("[register] no paired donors to register", flush=True)
        return pd.DataFrame()

    rows = []
    for _, r in man.iterrows():
        donor, roi = r["donor_id"], r["roi"]
        res = register_donor(cfg, donor, roi, use_images=use_images)
        if res is None:
            continue
        out = registration_dir(cfg, donor, roi)
        out.mkdir(parents=True, exist_ok=True)
        res.transform.save(out)
        (out / "qc.json").write_text(json.dumps(res.to_dict(), indent=2))

        rows.append({
            "donor_id": donor, "roi": roi, "tissue": r.get("tissue", ""),
            "track": res.track,
            "n_inliers": res.n_inliers, "n_candidates": res.n_candidates,
            "scale_x": scale_report(res.transform)["scale_x"],
            "scale_y": scale_report(res.transform)["scale_y"],
            "anisotropy": scale_report(res.transform)["anisotropy"],
            "rmse_um": res.rmse_um, "tissue_dice": res.tissue_dice,
            "islet_dice": res.islet_dice, "mirrored": res.mirrored,
            "rotation_deg": res.transform.rotation_deg, "scale": res.transform.scale,
            "max_disp_um": res.transform.max_displacement_um,
            "score": res.score(), "stages": "|".join(res.stages),
        })
        print(f"[register] {donor}/{roi}: track={res.track} stages={'|'.join(res.stages)} "
              f"rot={res.transform.rotation_deg:+.1f}deg "
              f"tissue_dice={res.tissue_dice:.3f} islet_rmse={res.rmse_um:.1f}um",
              flush=True)

    return pd.DataFrame(rows)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Register PhenoCycler <-> Xenium serial sections")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default=None,
                    help="one ROI (default: every ROI of --tissue, or all tissues)")
    ap.add_argument("--tissue", default=None, help="restrict to one tissue's ROIs")
    ap.add_argument("--no-images", action="store_true",
                    help="point-set track only (skip all image reads)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    summary = run_register(cfg, args.donors, roi=args.roi, tissue=args.tissue,
                           use_images=not args.no_images)
    if len(summary):
        print()
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
