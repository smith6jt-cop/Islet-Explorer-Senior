"""Registration: can a known transform be recovered, and are the conventions right?

The headline tests apply a *known* rigid / affine / non-rigid warp to synthetic data and
assert that registration recovers it. Everything else in the spatial pipeline is downstream
of that being true, and several of the failure modes here are silent — a wrong displacement
sign doubles the misalignment instead of removing it, and a reflection the fitter cannot
represent comes back as a confident, wrong rotation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.integration.rasterize import (GridSpec, grid_for, normalize_image,
                                               rasterize_points)
from phenocycler.integration import register  # noqa: F401
from phenocycler.integration.register import (dice, icp, nearest_neighbour_rmse,
                                              ransac_affine, refine_bspline, refine_ecc,
                                              register_image, register_point_set,
                                              scale_is_physical, scale_report)
from phenocycler.integration.transform import (Transform, affine_from_pixels, affine_to_pixels,
                                               fit_affine_lstsq, fit_rigid_umeyama, from_rigid)


def _islets(n, seed, lo=300, hi=3700):
    r = np.random.default_rng(seed)
    xy = r.uniform(lo, hi, (n, 2))
    return pd.DataFrame({
        "structure_id": [f"s{i}" for i in range(n)],
        "x_um": xy[:, 0], "y_um": xy[:, 1],
        "n_cells": r.integers(25, 180, n).astype(float),
        "area_um2": r.uniform(3e3, 1.2e4, n),
        "eccentricity": r.uniform(0.1, 0.6, n),
    })


# --------------------------------------------------------------------------- #
# Transform algebra
# --------------------------------------------------------------------------- #

def test_umeyama_recovers_rotation_and_translation():
    rng = np.random.default_rng(0)
    truth = from_rigid(17.0, (250.0, -120.0), center_um=(2000, 2000))
    src = rng.uniform(0, 4000, (60, 2))
    est = fit_rigid_umeyama(src, truth.apply(src))
    assert est.rotation_deg == pytest.approx(17.0, abs=1e-6)
    assert est.scale == pytest.approx(1.0, abs=1e-9)
    assert np.abs(est.apply(src) - truth.apply(src)).max() < 1e-8


def test_umeyama_reflection_must_be_representable():
    """Serial sections are sometimes mounted face-down.

    Textbook Umeyama forces det=+1, which cannot represent that and silently returns a wrong
    proper rotation instead. `allow_reflection` is what makes flipped donors registrable.
    """
    rng = np.random.default_rng(1)
    truth = from_rigid(-73.0, (90.0, 410.0), mirror=True, center_um=(2000, 2000))
    src = rng.uniform(0, 4000, (40, 2))
    dst = truth.apply(src)

    forced = fit_rigid_umeyama(src, dst, allow_reflection=False)
    assert not forced.is_mirrored
    assert np.abs(forced.apply(src) - dst).max() > 100.0      # genuinely wrong

    ok = fit_rigid_umeyama(src, dst, allow_reflection=True)
    assert ok.is_mirrored
    assert np.abs(ok.apply(src) - dst).max() < 1e-6


def test_affine_lstsq_recovers_shear():
    rng = np.random.default_rng(2)
    A = np.array([[1.02, 0.06, 133.0], [-0.04, 0.98, -87.0]])
    src = rng.uniform(0, 4000, (50, 2))
    dst = src @ A[:, :2].T + A[:, 2]
    assert np.abs(fit_affine_lstsq(src, dst).apply(src) - dst).max() < 1e-8


def test_pixel_micron_affine_roundtrip():
    """The single conversion point between the two coordinate systems."""
    grid = GridSpec(-100.0, -50.0, 2.5, 800, 900)
    tf = from_rigid(23.0, (140.0, -75.0), center_um=(1000, 1000))
    back = affine_from_pixels(affine_to_pixels(tf, grid), grid)
    assert np.abs(back.affine - tf.affine).max() < 1e-9


def test_compose_and_invert():
    rng = np.random.default_rng(3)
    a = from_rigid(15.0, (100.0, 40.0))
    b = from_rigid(-40.0, (-70.0, 210.0), scale=1.05)
    pts = rng.uniform(0, 3000, (30, 2))
    assert np.abs(a.compose(b).apply(pts) - a.apply(b.apply(pts))).max() < 1e-8
    assert np.abs(a.invert().apply(a.apply(pts)) - pts).max() < 1e-8


def test_displacement_inverse_roundtrip():
    """`inverse_apply` must undo `apply` including the non-rigid part.

    `invert()` deliberately drops the displacement; using it for warping instead of
    `inverse_apply` makes the whole non-rigid stage a silent no-op.
    """
    grid = GridSpec(0.0, 0.0, 10.0, 400, 400)
    rows, cols = np.mgrid[0:400, 0:400]
    disp = np.stack([30 * np.sin(rows / 60.0), 20 * np.cos(cols / 70.0)]).astype(np.float32)
    tf = Transform(affine=from_rigid(9.0, (55.0, -20.0)).affine, displacement=disp, grid=grid)

    rng = np.random.default_rng(4)
    pts = rng.uniform(400, 3600, (300, 2))
    assert np.abs(tf.inverse_apply(tf.apply(pts)) - pts).max() < 1.0
    assert tf.max_displacement_um > 20.0


def test_transform_io_roundtrip(tmp_path):
    grid = GridSpec(0.0, 0.0, 10.0, 60, 60)
    disp = np.ones((2, 60, 60), dtype=np.float32) * 3.0
    tf = Transform(affine=from_rigid(11.0, (5.0, 6.0)).affine, displacement=disp, grid=grid,
                   meta={"track": "image"})
    tf.save(tmp_path)
    back = Transform.load(tmp_path)
    assert np.allclose(back.affine, tf.affine)
    assert np.allclose(back.displacement, tf.displacement)
    assert back.grid.pixel_um == grid.pixel_um
    assert back.meta["track"] == "image"


# --------------------------------------------------------------------------- #
# Track A — point set
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("angle,tx,ty,mirror,drop,jitter", [
    (12.0, 220.0, -140.0, False, 0.00, 5.0),
    (47.0, -310.0, 260.0, False, 0.25, 15.0),
    (-73.0, 90.0, 410.0, True, 0.15, 10.0),      # mounted flipped
    (155.0, -40.0, -330.0, True, 0.30, 20.0),    # flipped + a third of islets missing
])
def test_point_set_recovers_known_transform(angle, tx, ty, mirror, drop, jitter):
    """Islets appear and vanish between sections, so the fit must survive missing points."""
    rng = np.random.default_rng(11)
    fixed = _islets(40, 5)
    truth = from_rigid(angle, (tx, ty), mirror=mirror, center_um=(2000, 2000))

    mv_xy = truth.invert().apply(fixed.loc[:, ["x_um", "y_um"]].to_numpy())
    mv_xy = mv_xy + rng.normal(0, jitter, mv_xy.shape)
    moving = fixed.copy()
    moving["x_um"], moving["y_um"] = mv_xy[:, 0], mv_xy[:, 1]
    moving["n_cells"] *= rng.uniform(0.8, 1.25, len(moving))
    moving = moving[rng.random(len(moving)) >= drop].reset_index(drop=True)
    moving = pd.concat([moving, _islets(6, 99)], ignore_index=True)   # section-specific islets

    res = register_point_set(moving, fixed, inlier_tol_um=150.0)
    got = res.transform
    resid = np.linalg.norm(got.apply(mv_xy) - fixed.loc[:, ["x_um", "y_um"]].to_numpy(), axis=1)

    assert got.is_mirrored == mirror
    assert got.scale == pytest.approx(1.0, abs=0.05)
    assert np.median(resid) < 3 * max(jitter, 5.0)
    assert res.n_inliers >= 10


def test_point_set_declines_when_there_is_nothing_to_match():
    res = register_point_set(_islets(2, 1), _islets(2, 2))
    assert res.n_inliers == 0
    assert any("too few structures" in n for n in res.notes)


def test_ransac_rejects_implausible_scale():
    """Both instruments report calibrated microns, so a large scale change is a bad sample."""
    rng = np.random.default_rng(7)
    mv = rng.uniform(0, 1000, (20, 2))
    fx = mv * 3.0                                   # 3x scale: physically impossible here
    cands = [(i, i) for i in range(20)]
    tf, inliers = ransac_affine(mv, fx, cands, n_iter=500, rng=rng)
    assert tf is None or tf.scale <= 1.25


def test_icp_refines_a_rough_seed():
    rng = np.random.default_rng(9)
    fixed = _islets(35, 3).loc[:, ["x_um", "y_um"]].to_numpy()
    truth = from_rigid(8.0, (60.0, -35.0), center_um=(2000, 2000))
    moving = truth.invert().apply(fixed) + rng.normal(0, 3.0, fixed.shape)

    seed = from_rigid(5.0, (40.0, -20.0), center_um=(2000, 2000))
    before, _ = nearest_neighbour_rmse(seed.apply(moving), fixed, max_dist_um=400)
    after_tf = icp(moving, fixed, initial=seed, max_dist_um=400)
    after, _ = nearest_neighbour_rmse(after_tf.apply(moving), fixed, max_dist_um=400)
    assert after < before


# --------------------------------------------------------------------------- #
# Track B — image
# --------------------------------------------------------------------------- #

def _blob_image(seed=21, n_blobs=14, pixel_um=8.0):
    from scipy import ndimage

    r = np.random.default_rng(seed)
    pts = np.vstack([r.normal(c, 60, (250, 2)) for c in r.uniform(400, 3600, (n_blobs, 2))])
    grid = grid_for(pts, pixel_um, margin_um=400)
    img = ndimage.gaussian_filter(normalize_image(rasterize_points(pts, grid)), 2.0)
    return pts, grid, img


@pytest.mark.parametrize("angle,tx,ty", [(0.0, 150.0, -90.0), (23.0, -200.0, 120.0),
                                         (-41.0, 80.0, 260.0)])
def test_image_track_recovers_known_rigid_transform(angle, tx, ty):
    from scipy import ndimage

    pts, grid, fixed = _blob_image()
    truth = from_rigid(angle, (tx, ty), center_um=(2000, 2000))
    mv_pts = truth.invert().apply(pts)
    moving = ndimage.gaussian_filter(normalize_image(rasterize_points(mv_pts, grid)), 2.0)

    res = register_image(moving, fixed, grid, nonrigid=False, coarse_step_deg=3.0)
    got = res.transform
    err = abs(((got.rotation_deg - truth.rotation_deg + 180) % 360) - 180)
    resid = np.linalg.norm(got.apply(mv_pts) - pts, axis=1)

    assert err < 1.5
    assert np.median(resid) < 30.0
    assert "coarse" in res.stages


def test_bspline_improves_a_nonrigid_warp_and_respects_the_cap():
    """The displacement sign convention and the safety cap, together.

    SimpleITK returns a fixed->moving field while `Transform.displacement` is moving->fixed.
    Getting that backwards does not fail loudly — it doubles the misalignment — so the test
    asserts that alignment actually *improves*.
    """
    from scipy import ndimage

    _pts, grid, fixed = _blob_image(seed=31, n_blobs=16, pixel_um=10.0)
    H, W = grid.shape
    rows, cols = np.mgrid[0:H, 0:W]
    dx = 12.0 * np.sin(2 * np.pi * rows / H)
    dy = 12.0 * np.cos(2 * np.pi * cols / W)
    moving = ndimage.map_coordinates(fixed, [rows + dy, cols + dx], order=1, mode="constant")

    before = dice(moving > 0.05, fixed > 0.05)
    out = refine_bspline(moving, fixed, grid, max_disp_um=250.0, grid_spacing_um=400.0,
                         iterations=150)
    if out is None:
        pytest.skip("SimpleITK unavailable")
    field, field_grid = out
    tf = Transform(affine=np.array([[1.0, 0, 0], [0, 1.0, 0]]), displacement=field,
                   grid=field_grid)
    after = dice(tf.warp_mask(moving, field_grid) > 0.05, fixed > 0.05)
    assert after > before + 0.05, f"non-rigid stage did not improve alignment ({before}->{after})"

    capped = refine_bspline(moving, fixed, grid, max_disp_um=30.0, grid_spacing_um=400.0,
                            iterations=150)
    assert np.linalg.norm(capped[0], axis=0).max() <= 30.0 + 1e-4


def test_ecc_failure_is_not_fatal():
    """A donor whose ECC diverges must not abort a cohort run."""
    grid = GridSpec(0.0, 0.0, 10.0, 64, 64)
    blank = np.zeros(grid.shape, dtype=np.float32)
    tf, val = refine_ecc(blank, blank, grid, Transform(), motion="affine")
    assert np.isnan(val) or np.isfinite(val)


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #

def test_dice_bounds():
    a = np.zeros((10, 10), bool)
    a[2:8, 2:8] = True
    assert dice(a, a) == pytest.approx(1.0)
    assert dice(a, ~a) == pytest.approx(0.0)
    assert np.isnan(dice(np.zeros((4, 4), bool), np.zeros((4, 4), bool)))


def test_nearest_neighbour_rmse_uses_mutual_pairs():
    """One-way NN distance can be driven to zero by collapsing everything onto one point."""
    src = np.array([[0.0, 0.0], [1000.0, 1000.0]])
    dst = np.array([[0.0, 0.0]])
    rmse, n = nearest_neighbour_rmse(src, dst, max_dist_um=50.0)
    assert n == 1 and rmse == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# Scale plausibility
#
# PhenoCycler and Xenium both report calibrated microns, so the transform between two serial
# sections is rigid up to real tissue stretch. ECC's affine stage has six free parameters and,
# driven by a partial-overlap silhouette, will trade an unphysical stretch for correlation.
# Donor 6436 came back holding sx=1.1057 against sy=0.9989 -- an 11% x-only stretch worth
# 509 um of drift across a 10 mm section, and invisible in a whole-section overlay.
# --------------------------------------------------------------------------- #

def _linear(sx, sy, shear=0.0, angle_deg=0.0):
    th = np.radians(angle_deg)
    rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    lin = rot @ np.array([[sx, shear * sy], [0.0, sy]])
    return Transform(affine=np.column_stack([lin, np.zeros(2)]))


def test_scale_report_separates_the_two_axes_isotropic_scale_hides():
    """One isotropic number cannot show a single-axis stretch; that is the whole point."""
    tf = _linear(1.1057, 0.9989)
    report = scale_report(tf)
    assert report["scale_x"] == pytest.approx(1.1057, abs=1e-4)
    assert report["scale_y"] == pytest.approx(0.9989, abs=1e-4)
    assert report["anisotropy"] == pytest.approx(1.1069, abs=1e-3)
    # The isotropic figure looks nearly innocuous while the x axis is 11% out.
    assert report["linear_scale"] == pytest.approx(1.0509, abs=1e-3)
    assert tf.scale == pytest.approx(1.0509, abs=1e-3)


def test_a_single_axis_stretch_is_refused_however_ordinary_its_isotropic_scale_looks():
    ok, why = scale_is_physical(_linear(1.1057, 0.9989))
    assert not ok
    assert "anisotropy" in why


def test_a_rigid_transform_at_any_angle_is_accepted():
    for angle in (-30.0, -9.19, 0.0, 5.11, 120.0):
        ok, why = scale_is_physical(_linear(1.0, 1.0, angle_deg=angle))
        assert ok, (angle, why)


def test_a_plausible_serial_section_stretch_is_accepted():
    """6442's real transform: 1.2% anisotropy, 0.4% off unit scale. That is tissue, not error."""
    ok, why = scale_is_physical(_linear(1.0020, 0.9898))
    assert ok, why


def test_a_large_isotropic_scale_is_refused_even_with_no_anisotropy():
    ok, why = scale_is_physical(_linear(1.20, 1.20))
    assert not ok
    assert "linear scale" in why


def test_the_bounds_are_configurable_in_both_directions():
    tf = _linear(1.05, 1.0)
    assert not scale_is_physical(tf, max_anisotropy=1.01)[0]
    assert scale_is_physical(tf, max_anisotropy=1.10)[0]


def test_a_degenerate_linear_part_is_refused_rather_than_dividing_by_zero():
    ok, _why = scale_is_physical(Transform(affine=np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])))
    assert not ok


def test_image_track_keeps_the_rigid_fit_when_the_affine_stage_goes_unphysical():
    """End to end, on data that genuinely rewards an unphysical fit.

    The moving side really is stretched 15% along one axis, so the affine stage *will* beat
    the rigid one -- a rigid transform cannot represent that stretch at all. That is exactly
    the situation the gate exists for: preferring the fit that scores better would bake a
    stretch the instruments' calibration says cannot be real into every downstream micron.
    """
    from scipy import ndimage

    # A smooth, structured raster -- ECC optimises an intensity correlation and needs real
    # gradients, which a sparse rasterised point cloud does not provide.
    rng = np.random.default_rng(4)
    grid = GridSpec(0.0, 0.0, 8.0, 320, 320)
    field = ndimage.gaussian_filter(rng.random((grid.height, grid.width)).astype(np.float32), 9)
    fixed = np.zeros_like(field)
    fixed[30:290, 30:290] = field[30:290, 30:290]
    fixed = normalize_image(fixed)

    th = np.radians(6.0)
    rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    stretched = Transform(affine=np.column_stack([rot @ np.diag([1.15, 1.0]),
                                                  np.array([90.0, -40.0])]), grid=grid)
    # moving is the image that, warped by `stretched`, reproduces `fixed`.
    moving = normalize_image(stretched.invert().warp_mask(fixed, grid))

    loose = register_image(moving, fixed, grid, nonrigid=False,
                           max_anisotropy=1.30, max_scale_deviation=0.30)
    assert "ecc_affine" in loose.stages, "the stretch should make the affine stage win"
    assert scale_report(loose.transform)["anisotropy"] > 1.05

    strict = register_image(moving, fixed, grid, nonrigid=False,
                            max_anisotropy=1.03, max_scale_deviation=0.05)
    assert "ecc_affine" not in strict.stages
    assert any("affine refused" in n for n in strict.notes)
    # What survives is the rigid fit: scale 1 on both axes by construction.
    assert scale_is_physical(strict.transform, max_anisotropy=1.001,
                             max_scale_deviation=0.001)[0]


# --------------------------------------------------------------------------- #
# Coarse candidate selection — the defect that made donor 6436's transform wrong
# --------------------------------------------------------------------------- #

def _candidate(corr: float) -> Transform:
    tf = Transform()
    tf.meta["coarse_corr"] = corr
    return tf


def test_a_plateau_of_coarse_candidates_is_refined_not_a_fixed_rank():
    """Donor 6436's correct hypothesis ranked 4th of 144 and the cut was 3.

    Its top four scored 0.4433 / 0.4367 / 0.4358 / 0.4352 — a 1.8% spread across hypotheses
    175 deg apart. When one section is a sub-field of the other the coarse score bunches like
    this, so the rule has to be "how close to the best", not "how many".
    """
    scores = [0.4433, 0.4367, 0.4358, 0.4352, 0.3100, 0.3050]
    kept = register.select_coarse_candidates(
        [_candidate(c) for c in scores], keep_frac=0.97, min_candidates=1, max_candidates=24)
    assert len(kept) == 4, "the 0.4352 hypothesis is inside the plateau and must be refined"
    assert [c.meta["coarse_corr"] for c in kept] == scores[:4]


def test_a_sharp_coarse_winner_still_gets_company():
    """A floor, so one lucky hypothesis is never refined alone and unchecked."""
    kept = register.select_coarse_candidates(
        [_candidate(c) for c in (0.90, 0.20, 0.19, 0.18, 0.17)],
        keep_frac=0.97, min_candidates=3, max_candidates=24)
    assert len(kept) == 3


def test_a_flat_coarse_landscape_is_capped():
    """And a ceiling, so a featureless landscape does not run ECC on every hypothesis."""
    kept = register.select_coarse_candidates(
        [_candidate(0.5) for _ in range(144)],
        keep_frac=0.97, min_candidates=8, max_candidates=24)
    assert len(kept) == 24


# --------------------------------------------------------------------------- #
# Masked ECC
# --------------------------------------------------------------------------- #

def _partial_overlap_pair(shift_px: int = 4):
    """A fixed 'section' with a moving 'sub-field' covering part of it, plus its mask.

    Smoothed rather than raw noise: ECC optimises a correlation through image gradients, and
    white noise has no gradient structure to follow, so it diverges on it for reasons that
    have nothing to do with what these tests are about.
    """
    from scipy import ndimage

    rng = np.random.default_rng(0)
    fixed = ndimage.gaussian_filter(rng.random((160, 160)).astype(np.float32), 3.0)
    fixed = (fixed - fixed.min()) / max(float(np.ptp(fixed)), 1e-9)
    moving = np.zeros_like(fixed)
    mask = np.zeros(fixed.shape, dtype=bool)
    moving[20:90, 20:90] = fixed[20 - shift_px:90 - shift_px, 20:90]
    mask[20:90, 20:90] = True
    return moving, fixed, mask


def test_a_masked_ecc_value_is_not_comparable_to_an_unmasked_one():
    """Both are correlations, but over different supports.

    Measured at 0.31 against 0.49 for the same transform on donor 6436. register_image
    compares these values across candidates, so they must be masked uniformly or not at all —
    this pins the fact so nobody later assumes they interchange.
    """
    moving, fixed, mask = _partial_overlap_pair(shift_px=2)
    grid = GridSpec(0.0, 0.0, 2.0, fixed.shape[0], fixed.shape[1])
    seed = Transform()

    _tf_u, val_unmasked = register.refine_ecc(moving, fixed, grid, seed, motion="euclidean")
    _tf_m, val_masked = register.refine_ecc(moving, fixed, grid, seed, motion="euclidean",
                                            moving_mask=mask)
    assert np.isfinite(val_unmasked) and np.isfinite(val_masked)
    assert not np.isclose(val_unmasked, val_masked, atol=1e-3)


def test_masking_widens_the_basin_ecc_can_converge_from():
    """The reason to mask, beyond the score: the empty part of the canvas fights the fit.

    At a 6 px offset the unmasked fit diverges outright — every pixel of fixed-frame tissue
    with no moving counterpart is pulling on a correlation that is supposed to be about the
    part that does correspond. Masked, the same pair converges to 0.997.
    """
    moving, fixed, mask = _partial_overlap_pair(shift_px=6)
    grid = GridSpec(0.0, 0.0, 2.0, fixed.shape[0], fixed.shape[1])

    _tf_u, val_unmasked = register.refine_ecc(moving, fixed, grid, Transform(),
                                              motion="euclidean")
    _tf_m, val_masked = register.refine_ecc(moving, fixed, grid, Transform(),
                                            motion="euclidean", moving_mask=mask)
    assert not np.isfinite(val_unmasked)
    assert np.isfinite(val_masked) and val_masked > 0.9


def test_the_mask_is_resampled_to_whatever_resolution_ecc_actually_ran_at():
    """ndimage.zoom rounds, and a one-pixel disagreement is a hard cv2.error, not a warning."""
    moving, fixed, mask = _partial_overlap_pair()
    grid = GridSpec(0.0, 0.0, 2.0, fixed.shape[0], fixed.shape[1])
    # Force the internal downsample by capping well below the raster's own size.
    _tf, val = register.refine_ecc(moving, fixed, grid, Transform(), motion="euclidean",
                                   moving_mask=mask, max_pixels=4_000)
    assert np.isfinite(val), "a shape mismatch would have surfaced as nan via cv2.error"


# --------------------------------------------------------------------------- #
# The joint residual search
# --------------------------------------------------------------------------- #

def test_the_residual_vote_finds_a_joint_rotation_and_translation_offset():
    """The failure a rotation-only sweep cannot see, by construction.

    Rotating about the warped moving centroid fixes that centroid, so no angle in that sweep
    carries a translation. Here the answer is both, and the vote has to recover both.
    """
    rng = np.random.default_rng(3)
    fixed = rng.uniform(0, 8000, size=(150, 2))
    angle, offset = -14.0, np.array([-700.0, 400.0])
    th = np.radians(-angle)
    c, s = np.cos(th), np.sin(th)
    centre = fixed.mean(axis=0)
    moving = (fixed[:80] - offset - centre) @ np.array([[c, s], [-s, c]]) + centre

    out = register.residual_vote(moving, fixed, angles_deg=np.arange(-30, 30, 2.0))
    assert out["contrast"] > 3.0
    assert abs(out["angle_deg"] - angle) <= 2.0
    assert out["mirror"] is False

    # The recovered residual is expressed about the MOVING set's own centroid, which is not
    # the centre the offset above was constructed around — so the parameters are checked by
    # applying them, the way pheno_residual_search.py does, not by comparing them literally.
    recovered = from_rigid(out["angle_deg"], (out["dx_um"], out["dy_um"]),
                           center_um=tuple(float(v) for v in moving.mean(axis=0)))
    placed = recovered.apply(moving)
    from scipy.spatial import cKDTree

    d, _i = cKDTree(fixed).query(placed)
    assert np.median(d) < 60.0, "the recovered residual must actually bring them together"


def test_two_unrelated_point_sets_produce_no_peak():
    """The statistic has to be able to say 'these do not correspond', or it says nothing."""
    rng = np.random.default_rng(5)
    out = register.residual_vote(rng.uniform(0, 8000, size=(120, 2)),
                                 rng.uniform(0, 8000, size=(300, 2)),
                                 angles_deg=np.arange(-30, 30, 2.0))
    assert out["contrast"] < 3.0


# --------------------------------------------------------------------------- #
# Displacement invertibility
# --------------------------------------------------------------------------- #

def test_a_folded_displacement_field_is_refused_by_its_jacobian_not_its_magnitude():
    """max|D| cannot see folding. Donor 6436's first field sat at its cap and still folded."""
    grid = GridSpec(0.0, 0.0, 10.0, 60, 60)
    yy, xx = np.mgrid[0:60, 0:60]
    # A gradient steeper than the grid: x + D(x) stops being monotonic, so it has no inverse.
    field = np.zeros((2, 60, 60), np.float32)
    field[0] = (-30.0 * xx).astype(np.float32)
    folded = Transform(affine=np.array([[1.0, 0, 0], [0, 1.0, 0]]),
                       displacement=field, grid=grid)
    ok, why = register.displacement_is_invertible(folded)
    assert not ok and "folds" in why

    gentle = Transform(affine=np.array([[1.0, 0, 0], [0, 1.0, 0]]),
                       displacement=np.zeros((2, 60, 60), np.float32), grid=grid)
    ok, _why = register.displacement_is_invertible(gentle)
    assert ok


def test_an_affine_only_transform_is_invertible_by_definition():
    ok, why = register.displacement_is_invertible(Transform())
    assert ok and why == "affine only"


def _structured(shape=(220, 220), seed=11):
    """An image with asymmetric large-scale structure.

    Smoothed noise will not do: it is close enough to rotationally self-similar that a coarse
    orientation search cannot tell 165 deg from 10 deg, which makes it a bad model of a
    section and a worse test fixture.
    """
    from scipy import ndimage

    rng = np.random.default_rng(seed)
    img = np.zeros(shape, np.float32)
    for cy, cx, amp in ((60, 70, 1.0), (70, 150, 0.7), (150, 90, 0.85), (170, 160, 0.5),
                        (110, 115, 0.35)):
        yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
        img += amp * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 22.0 ** 2)))
    img += 0.05 * rng.random(shape)
    img = ndimage.gaussian_filter(img, 1.5)
    return ((img - img.min()) / max(float(np.ptp(img)), 1e-9)).astype(np.float32)


def test_a_large_rotation_seed_is_refined_as_a_residual_not_as_a_whole_transform():
    """ECC linearises around its estimate, so a far-from-identity seed diverges.

    Donor 6450's sections differ by ~170 deg and every one of its eight candidates returned
    "NaN encountered" -- masked and unmasked alike -- until the candidate was applied to the
    moving image first. 6436 and 6442 only ever worked because their coarse winners happened
    to land at -20 and +5 deg.

    The seed is built directly rather than taken from the coarse search, so what is under test
    is the refinement formulation and nothing else.
    """
    fixed = _structured()
    grid = GridSpec(0.0, 0.0, 2.0, fixed.shape[0], fixed.shape[1])
    centre = (grid.width * grid.pixel_um / 2.0, grid.height * grid.pixel_um / 2.0)
    truth = from_rigid(165.0, (0.0, 0.0), center_um=centre)
    moving = truth.invert().warp_mask(fixed, grid)

    _tf_direct, direct = register.refine_ecc(moving, fixed, grid, truth, motion="euclidean")
    pre = truth.warp_mask(moving.astype(np.float32), grid)
    residual, val = register.refine_ecc(pre, fixed, grid, Transform(), motion="euclidean")

    assert np.isfinite(val), "the residual formulation must converge on a 165 deg seed"
    assert abs(residual.compose(truth).rotation_deg - 165.0) < 5.0
    # The direct call is what fails on real data; if it happens to survive here, the residual
    # must still be at least as good, which is the property the pipeline relies on.
    if np.isfinite(direct):
        assert val >= direct - 1e-6


def test_the_mirror_search_survives_the_fine_stage():
    """MOTION_EUCLIDEAN cannot represent a reflection, so a mirrored seed used to be dropped.

    cv2 rejects the non-proper rotation the SVD orthogonalisation returns, `refine_ecc` caught
    the cv2.error and returned nan, and every mirrored candidate vanished -- silently, because
    the mirror hypothesis is advertised at the top of this module. Pre-warping puts the
    reflection in the image, where cv2 never sees it as a warp.
    """
    fixed = _structured(seed=12)
    grid = GridSpec(0.0, 0.0, 2.0, fixed.shape[0], fixed.shape[1])
    centre = (grid.width * grid.pixel_um / 2.0, grid.height * grid.pixel_um / 2.0)
    truth = from_rigid(20.0, (0.0, 0.0), mirror=True, center_um=centre)
    moving = truth.invert().warp_mask(fixed, grid)
    assert truth.is_mirrored

    _tf, direct = register.refine_ecc(moving, fixed, grid, truth, motion="euclidean")
    assert not np.isfinite(direct), "a reflected seed is not a MOTION_EUCLIDEAN warp"

    pre = truth.warp_mask(moving.astype(np.float32), grid)
    residual, val = register.refine_ecc(pre, fixed, grid, Transform(), motion="euclidean")
    assert np.isfinite(val), "a mirrored candidate must reach a finite ECC value"
    assert residual.compose(truth).is_mirrored, "the reflection must survive composition"
