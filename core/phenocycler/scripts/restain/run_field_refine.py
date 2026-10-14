#!/usr/bin/env python
"""Fit a displacement field to run_align's own per-window offsets, and iterate on it.

``run_intensity_refine.py`` fits ONE low-dimensional affine residual, which is the right model
when the two sections differ by an instrument-scale error.  It is the wrong model when the
section itself changed shape between the two acquisitions: donor 6533 had a kimwipe dragged
vertically across it after the PhenoCycler run and before Xenium, so the tissue is physically
sheared and no global affine -- six parameters for the whole section -- can register it.  Its
``refined`` gate reads 38.80 um.

What is available instead is ``run_align.py``'s own measurement grid.  Every occupied 500 um
window reports the shift that lands the warped Xenium nuclei on the IF DAPI, so the gate is
already a sampled displacement field; this stage turns those samples into a dense one on a
250 um grid, composes it with the source transform, and writes a NEW stem::

    python scripts/restain/run_align.py --donors 6533 --transform refined --windows occupied
    python scripts/restain/run_field_refine.py --donors 6533
    python scripts/restain/run_align.py --donors 6533 --transform refined_field --windows occupied
    python scripts/restain/run_field_refine.py --donors 6533 \
        --source-stem refined_field --output-stem refined_field2
    python scripts/restain/run_align.py --donors 6533 --transform refined_field2 --windows occupied

Each pass measures the residual the previous pass left behind, so the stage is deliberately
re-runnable on its own output: the source may already carry a field, in which case the new one
is composed onto it on the grid it already uses and ``field_iterations`` counts the chain.
That is the whole provenance of 6533's ``refined_field2``: 38.80 um FAIL -> 2.03 INCONCLUSIVE
-> 0.202 PASS.

**The field is defined on the whole fixed frame, anchored at micron zero.**  A transform is a
property of the frame, not of the measurement that estimated it: two runs that happened to
measure different windows must still write the same grid, and a cell outside the measured
region must read a defined displacement rather than fall off the end of the field and be
silently clamped to its edge.  So the grid spans the IF frame the gate measured against --
``record.pyramid_shapes[0] x record.pixel_um``, the same image ``run_align`` re-verifies
against the inventory before it measures anything -- from (0, 0), at ``FIELD_PIXEL_UM``, and
it is NOT the bounding box of whichever windows reported a peak.  For 6533 that is 97 x 41
nodes at origin (0, 0), pixel 250 um: the grid its released ``refined_field`` and
``refined_field2`` already carry, and the contract a reproduction has to hit before any of the
numbers below mean anything.  A second and later pass inherits the source field's own grid
unchanged, because the new field is added to it node for node.

**No window fold is withheld here, and that is a deliberate difference from the sibling.**
``run_intensity_refine.py`` can fit its affine on one checkerboard fold and reserve the other,
because six parameters are recoverable from half the windows.  A dense field has no parameters
to spare: dropping every second window doubles the spacing of the only measurements that
describe the tear.  The generalisation claim is made instead by leave-one-window-out -- each
window is predicted from a field fitted without it -- and the release claim is still made by a
fresh ``run_align.py`` pass over all occupied windows against the new stem, which re-measures
DAPI rather than re-reading the offsets this stage consumed.

The output is a new, explicitly named transform.  It never overwrites its source, and it is not
releasable until that fresh gate says PASS.

**Reproducing the released field.**  6533's ``refined_field``/``refined_field2`` predate this
module, and the estimator that made them was recovered by measuring the released npz rather
than by trusting a description of it.  It is ``--estimator griddata --grid-anchor half-node
--window-selection measured``, which are the defaults, and it rebuilds them like this::

    # pass 1: the refined gate (252 measured windows, median 38.80 um)
    python scripts/restain/run_field_refine.py --donors 6533 \
        --source-stem refined --output-stem refined_field3a \
        --inventory data/development/sameslide_6533/registration_inputs.csv \
        --images-dir data/development/sameslide_6533/images \
        --out data/development/sameslide_6533/out
    python scripts/restain/run_align.py --donors 6533 --transform refined_field3a \
        --windows occupied --window-fold all ...        # measures the residual pass 1 leaves
    # pass 2: that gate (318 measured windows), composed onto pass 1's own grid
    python scripts/restain/run_field_refine.py --donors 6533 \
        --source-stem refined_field3a --output-stem refined_field3 ...

Pass 1 rebuilt this way matches the released ``refined_field`` to 0.0000 um on all 97 x 41
nodes; pass 2 matches ``refined_field2`` to 1e-5 um, which is float32 storage.  Two details
had to be measured rather than assumed, and both are recorded in the constants below: the
node values are the estimator evaluated HALF A NODE (125 um) inside the corner the JSON
declares (``GRID_ANCHORS``), and pass 2 is fitted to all 318 measured windows, not to the 317
that also improve DAPI contrast (``WINDOW_SELECTIONS``).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import align as AL  # noqa: E402
import registration_inputs as RI  # noqa: E402
import run_align as RA  # noqa: E402

# The source-gate contract (identity against the inventory row, transform and displacement
# sha256, "all occupied windows") and the atomic write order are imported rather than restated:
# two refiners that disagree about what a verified source is would be a provenance bug, and a
# copy of those checks is exactly how that disagreement starts.
#
# ``IR._field_grid`` is deliberately NOT imported.  The sibling fits six affine parameters and
# only rasterises them, so a grid around the windows it happened to fit is enough for it.  A
# dense field is the transform, and the grid it lives on is part of what is released -- see
# ``_frame_grid`` and the module docstring.
import run_intensity_refine as IR  # noqa: E402
from run_register import OUT  # noqa: E402

from phenocycler.integration.rasterize import GridSpec  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402

#: Field node spacing.  Shared with ``run_intensity_refine`` and with ``refine.refine_local``;
#: half the 500 um probe step.  ``run_align`` anchors its occupied windows at micron zero and
#: this grid is anchored there too, so every window centre lands exactly on a node instead of
#: between four of them.
FIELD_PIXEL_UM = IR.FIELD_PIXEL_UM
#: The fixed frame is padded up to a whole block of this many microns before its node count is
#: taken.  Not decoration: 6533's frame is 19200 px x 0.5078 um = 9749.76 um, a quarter of a
#: micron below a node boundary, so without a block a 0.003% change in the fitted pixel size
#: would add a whole column and put every later field on a different grid from the one before
#: it.  A 2 mm block -- eight nodes, four probe windows -- makes the grid a property of the
#: section rather than of the last digit of a calibration, and it is the rule under which 6533's
#: released field is 97 x 41 at (0, 0).
FRAME_BLOCK_UM = 8.0 * FIELD_PIXEL_UM
#: Gaussian bandwidth of the field fit, in microns.  A window reports one shift for a 500 um
#: square, so the field cannot honestly resolve structure finer than a window; 0.75 of the
#: probe step blends each node with the few windows that physically overlap it and leaves the
#: field smooth between them, rather than ringing through measurements it must pass exactly.
FIELD_SMOOTH_UM = 0.75 * RA.WINDOW_STEP_UM
#: Below this total kernel weight from the OTHER windows, a window has no neighbour close
#: enough to be predicted from, so it is left out of the leave-one-window-out claim (it still
#: anchors the field around itself).  At the default bandwidth that is about 2 mm -- four window
#: widths -- from every other window.  This is a statement about the generalisation claim only:
#: it is NOT the rule for filling a node, because a node is not asked to be predicted from
#: anything, it is asked for the best estimate available.  The two used to be one constant, and
#: over the whole frame that silently turned 43% of the grid into a piecewise-constant
#: nearest-window fill with a 2 mm discontinuity ring around the tissue.
MIN_NEIGHBOUR_WEIGHT = 1e-6
#: A dense field needs the windows to blanket the section, not merely to exist.  Twice the
#: sibling's affine minimum, which is itself five times that model's six parameters.
MIN_FIELD_WINDOWS = 60
#: Offsets pinned at the rim of the search disk are lower bounds on a shift that ran off the
#: end of the measurement, so a field fitted mostly to them is fitted to the rim.  A minority
#: of them is expected on a torn section and is kept (see ``_measured_windows``).
MAX_TRUNCATED_FRACTION = 0.5
#: Reported for context: the largest shift a window is able to measure at all.
SEARCH_RADIUS_UM = AL.SEARCH_RADIUS_UM
#: The two estimators this stage can fit, and which one the released field used.
#:
#: ``griddata`` is the estimator that produced 6533's released ``refined_field`` and
#: ``refined_field2``, recovered by measurement rather than by memory: piecewise-linear
#: barycentric interpolation of the window offsets on their Delaunay triangulation
#: (``scipy.interpolate.griddata(..., method="linear")``), the NaN hull exterior filled with
#: the nearest window's own offset (``method="nearest"``), then one pass of
#: ``scipy.ndimage.gaussian_filter`` at ``GRIDDATA_SMOOTH_NODES`` on the node grid.  Rebuilt
#: that way, pass 1 reproduces the released ``refined_field`` field to 0.0000 um on every one
#: of its 97 x 41 nodes and pass 2 reproduces ``refined_field2`` to 1e-5 um, which is float32
#: storage.  It is therefore the default: this stage exists to make that field reproducible.
#:
#: ``kernel`` is the Gaussian-kernel regression this module was written around.  It is a
#: convex combination of the measurements, so it cannot invent a shift larger than the largest
#: one measured and it does not pass exactly through offsets that carry the estimator's own
#: precision -- properties the interpolant does not have.  It is kept, and it is a different
#: transform: it does not reproduce the released field and would have to earn its own gate.
ESTIMATORS = ("griddata", "kernel")
#: Gaussian smoothing applied to the interpolated field, in NODES (1 node = ``FIELD_PIXEL_UM``
#: = 250 um).  Part of the recovered estimator, not a free parameter.
GRIDDATA_SMOOTH_NODES = 1.0
#: How the node grid's corner is anchored when the field is SAMPLED.
#:
#: ``GridSpec`` declares ``pixel = (micron - origin) / pixel_um``, so node ``(row, col)`` of a
#: field whose JSON says ``origin_x_um = origin_y_um = 0`` is read back by ``Transform.apply``
#: at micron ``(col * 250, row * 250)`` -- the ``corner`` convention.  6533's released
#: ``refined_field`` and ``refined_field2`` were not built that way.  Measured against the
#: released npz, their node values are the estimator evaluated at ``((col + 0.5) * 250,
#: (row + 0.5) * 250)``: half a node -- 125 um -- inside the corner their own JSON declares.
#: As released, the field is therefore read 125 um away from the point it was estimated at, in
#: both axes.
#:
#: The offset is reproduced by default rather than quietly corrected, because it is not free.
#: Pass 2 measured the residual that pass 1 left behind AS APPLIED, i.e. with the half-node
#: lookup already in it, so ``refined_field2``'s 0.202 um gate is a gate on the whole chain
#: including the offset, and removing it from one pass without re-running the other would
#: produce a field that has never been measured.  ``--grid-anchor corner`` is the corrected
#: convention: a different transform, which has to earn its own gate.
GRID_ANCHORS = ("half-node", "corner")
#: Whether a window must demonstrably improve DAPI contrast to be fitted to.
#:
#: ``measured`` -- a finite centre and offset from a window ``run_align`` called measured -- is
#: what the released field used: pass 2 is reproduced exactly by all 318 measured windows and
#: not by the 317 that also gain contrast.  ``contrast`` additionally requires
#: ``contrast_shifted > contrast``, which is the sibling refiner's evidence test; on pass 1 the
#: two agree (all 252 measured windows gain contrast) so only pass 2 distinguishes them.
WINDOW_SELECTIONS = ("measured", "contrast")
#: Where ``align.py`` actually calls a peak truncated -- 0.9 of the search radius, not the
#: radius itself.  Quoting the radius in the refusal would overstate by 6 um the shift a
#: rim-pinned window is a lower bound on.
TRUNCATION_RADIUS_UM = AL.TRUNCATION_FRACTION_OF_R * AL.SEARCH_RADIUS_UM


class FieldRefineError(RuntimeError):
    """The field cannot be fitted from this gate under the declared contract."""


def _verified_source(donor_dir: Path, stem: str, *,
                     record: RI.RegistrationInput) -> tuple[Transform, dict, Path]:
    """The sibling's source contract, re-raised under this stage's error type."""
    try:
        return IR._verified_source(donor_dir, stem, record=record)
    except IR.IntensityRefineError as exc:
        raise FieldRefineError(str(exc)) from exc


def _measured_windows(gate: dict, *, selection: str = "measured") -> list[dict]:
    """The windows whose offsets the field is allowed to be fitted to.

    ``selection="measured"`` takes every window ``run_align`` called measured, which is what
    the released field was fitted to.  ``selection="contrast"`` adds the sibling's evidence
    test -- the peak must be demonstrably useful, i.e. moving by it increases DAPI contrast.
    Either way a TRUNCATED window is kept here, which is a deliberate difference from the
    sibling.  The sibling rejects truncation because a magnitude clipped at the
    search rim biases an affine gradient over the whole section.  A field is local: a clipped
    window still carries the right direction at the right place, dropping it would punch a hole
    in the field exactly where the tissue moved most, and the under-correction it causes is
    what the next iteration measures and removes.

    A window's CENTRE is checked for finiteness as well as its offset.  ``AL.Offset`` defaults
    ``cx_um``/``cy_um`` to NaN and only a measured window sets them, so a malformed row that
    still says ``"measured"`` can carry a NaN centre behind a perfectly good-looking offset.
    That is not one bad window: a single NaN coordinate makes every kernel weight NaN, so every
    node fails its support test, every node takes the fallback, and ``argmin`` over a row of
    NaNs returns that same window -- the entire field silently collapses onto one offset.  The
    offset check alone never sees it.
    """
    if selection not in WINDOW_SELECTIONS:
        raise FieldRefineError(f"unknown window selection {selection!r}")
    out = []
    for w in gate["alignment"].get("windows", []):
        if w.get("reason") != "measured":
            continue
        if selection == "contrast":
            if w.get("contrast") is None or w.get("contrast_shifted") is None:
                continue
            if not w["contrast_shifted"] > w["contrast"]:
                continue
        cx, cy = w.get("cx_um"), w.get("cy_um")
        if cx is None or cy is None or not (np.isfinite(cx) and np.isfinite(cy)):
            continue
        dx, dy = w.get("dx_um"), w.get("dy_um")
        if dx is None or dy is None or not (np.isfinite(dx) and np.isfinite(dy)):
            continue
        out.append(w)
    if len(out) < MIN_FIELD_WINDOWS:
        raise FieldRefineError(
            f"only {len(out)} occupied windows measured a DAPI peak that improves contrast; "
            f"a {FIELD_PIXEL_UM:.0f} um field needs {MIN_FIELD_WINDOWS}. Re-gate with "
            f"--windows occupied, or accept that this section cannot be refined non-rigidly")
    truncated = sum(1 for w in out if w.get("truncated"))
    if truncated > MAX_TRUNCATED_FRACTION * len(out):
        raise FieldRefineError(
            f"{truncated} of {len(out)} windows are pinned at the "
            f"{TRUNCATION_RADIUS_UM:.0f} um rim of the {SEARCH_RADIUS_UM:.0f} um search disk; "
            f"their magnitudes are lower bounds, so a field fitted to them is fitted to the "
            f"rim of the search disk rather than to the tissue")
    return out


def _linear_nearest(points_um: np.ndarray, offsets_um: np.ndarray):
    """The two interpolants ``griddata`` builds, kept so they are built once, not per query."""
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    return (LinearNDInterpolator(points_um, offsets_um),
            NearestNDInterpolator(points_um, offsets_um))


def _interpolate(fit: dict, xy_um: np.ndarray) -> np.ndarray:
    """Barycentric-linear inside the hull of the windows, nearest window outside it.

    This is ``scipy.interpolate.griddata(..., method="linear")`` with its NaN exterior filled
    by ``method="nearest"``, which is what the released field did, written as two persistent
    interpolators rather than two triangulations per call.  UNSMOOTHED: the Gaussian pass is
    a property of the node grid, not of a query point, so it lives in
    :func:`_predict_on_grid`.
    """
    xy = np.asarray(xy_um, dtype=np.float64).reshape(-1, 2)
    linear, nearest = fit["interpolators"]
    out = np.asarray(linear(xy), dtype=np.float64).reshape(-1, 2)
    outside = ~np.all(np.isfinite(out), axis=1)
    if np.any(outside):
        out[outside] = np.asarray(nearest(xy[outside]), dtype=np.float64).reshape(-1, 2)
    return out


def _predict_on_grid(fit: dict, xy_um: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """The fitted field at ``shape``-ordered query points, smoothed if the estimator says so.

    ``xy_um`` is the node grid's query points in row-major order -- for a second pass they are
    the nodes displaced by the source field, not the nodes themselves -- so the Gaussian is
    over nodes, in node units, exactly as the released field was built.
    """
    if fit["estimator"] == "kernel":
        return _predict(fit, xy_um)
    from scipy.ndimage import gaussian_filter
    raw = _interpolate(fit, xy_um)
    out = np.empty_like(raw)
    for k in range(2):
        out[:, k] = gaussian_filter(raw[:, k].reshape(shape),
                                    sigma=GRIDDATA_SMOOTH_NODES).ravel()
    return out


def _heldout_griddata(p: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Leave-one-window-out prediction from the interpolant, window by window.

    The kernel estimator gets this for free -- zero the diagonal of the weight matrix -- but an
    interpolant has no weights to zero: dropping a window changes the triangulation, so the
    only honest held-out prediction is to rebuild it without that window.  252-318 Delaunay
    triangulations of ~300 points take about a second.
    """
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    n = len(p)
    out = np.empty((n, 2), dtype=np.float64)
    for i in range(n):
        keep = np.ones(n, dtype=bool)
        keep[i] = False
        value = np.asarray(LinearNDInterpolator(p[keep], d[keep])(p[i:i + 1]),
                           dtype=np.float64).reshape(1, 2)
        if not np.all(np.isfinite(value)):
            value = np.asarray(NearestNDInterpolator(p[keep], d[keep])(p[i:i + 1]),
                               dtype=np.float64).reshape(1, 2)
        out[i] = value[0]
    return out


def _node_microns(grid: GridSpec, anchor: str) -> np.ndarray:
    """The fixed-frame microns each node of ``grid`` is ESTIMATED at, row-major.

    ``corner`` is what ``GridSpec`` declares and what ``Transform.apply`` reads back.
    ``half-node`` is what the released field was built with -- see ``GRID_ANCHORS``.
    """
    if anchor not in GRID_ANCHORS:
        raise FieldRefineError(f"unknown grid anchor {anchor!r}")
    rows, cols = np.mgrid[:grid.height, :grid.width]
    u = grid.to_microns(np.column_stack([cols.ravel(), rows.ravel()]))
    if anchor == "half-node":
        u = u + 0.5 * grid.pixel_um
    return u


def _fit_field(windows: list[dict], *, estimator: str = "griddata") -> dict:
    """Gaussian-kernel regression of the measured offsets, plus its own generalisation.

    A weighted mean, not an interpolant: the field is a convex combination of measurements, so
    it cannot invent a shift larger than the largest one measured, and it does not have to pass
    exactly through offsets that carry the estimator's own ~0.03 um precision plus whatever the
    window's nuclei could not resolve.  ``heldout_*`` predicts each window from a fit that
    never saw it, which is the only claim this stage can make about a node between windows.
    """
    if estimator not in ESTIMATORS:
        raise FieldRefineError(f"unknown estimator {estimator!r}")
    p = np.asarray([[w["cx_um"], w["cy_um"]] for w in windows], dtype=np.float64)
    d = np.asarray([[w["dx_um"], w["dy_um"]] for w in windows], dtype=np.float64)
    baseline = np.linalg.norm(d, axis=1)

    interpolators = None
    if estimator == "kernel":
        weight = _kernel(p, p)
        np.fill_diagonal(weight, 0.0)
        total = weight.sum(axis=1)
        # A window with no neighbour cannot be predicted from the others, so it is left out of
        # the generalisation claim rather than refused outright -- it still anchors the field
        # around itself.  Too many of them means the set is shattered, not smooth.
        supported = total > MIN_NEIGHBOUR_WEIGHT
        if int(supported.sum()) < MIN_FIELD_WINDOWS:
            raise FieldRefineError(
                f"only {int(supported.sum())} of {len(windows)} windows have another window "
                f"within {FIELD_SMOOTH_UM:.0f} um; the offsets are too scattered to fit a "
                f"field")
        heldout = (weight[supported] @ d) / total[supported, None]
    else:
        # Every window is predictable from an interpolant rebuilt without it: outside the
        # remaining hull it falls back to the nearest window, so none is dropped.
        supported = np.ones(len(p), dtype=bool)
        interpolators = _linear_nearest(p, d)
        heldout = _heldout_griddata(p, d)
    residual = np.linalg.norm(d[supported] - heldout, axis=1)
    if float(np.median(residual)) >= float(np.median(baseline[supported])):
        raise FieldRefineError(
            f"leave-one-window-out residual {np.median(residual):.3f} um does not improve on "
            f"the measured {np.median(baseline):.3f} um: these offsets have no spatial "
            f"coherence at {FIELD_SMOOTH_UM:.0f} um, so a field fitted to them would only "
            f"memorise the windows it was fitted to")
    return {
        "points_um": p,
        "offsets_um": d,
        "estimator": estimator,
        "interpolators": interpolators,
        "n_windows": len(windows),
        "n_truncated": int(sum(1 for w in windows if w.get("truncated"))),
        "baseline_median_um": float(np.median(baseline)),
        "baseline_p90_um": float(np.percentile(baseline, 90)),
        "baseline_max_um": float(np.max(baseline)),
        "heldout_median_um": float(np.median(residual)),
        "heldout_p90_um": float(np.percentile(residual, 90)),
        "n_heldout": int(supported.sum()),
    }


def _square_distances(xy_um: np.ndarray, points_um: np.ndarray) -> np.ndarray:
    dx = xy_um[:, 0][:, None] - points_um[None, :, 0]
    dy = xy_um[:, 1][:, None] - points_um[None, :, 1]
    return dx * dx + dy * dy


def _kernel(xy_um: np.ndarray, points_um: np.ndarray) -> np.ndarray:
    return np.exp(-_square_distances(xy_um, points_um)
                  / (2.0 * FIELD_SMOOTH_UM * FIELD_SMOOTH_UM))


def _predict(fit: dict, xy_um: np.ndarray, *, block: int = 4096) -> np.ndarray:
    """The fitted field at arbitrary fixed-frame microns, in blocks over the query points.

    The nearest-window fallback fires on exactly one condition: the weighted mean is not
    computable, i.e. the kernel weights underflowed to zero outright and the ratio would be
    0/0.  At the default bandwidth that is about 14 mm from every window, and the point of the
    branch is to return the nearest window's own offset there instead of a NaN -- by DISTANCE,
    because ``argmin`` over a row of zeros would silently return the first window.

    It is deliberately not a cutoff at some weight the author thinks is small.  A weighted mean
    of 250-odd offsets is defined and finite long after any single weight stops mattering, and
    it tends smoothly to the nearest window's offset as the query point recedes.  Replacing it
    early with the nearest value instead puts a discontinuity into the field at that radius --
    over the whole frame, with a 1e-6 cutoff, a ring at 2 mm around the tissue and a
    piecewise-constant plateau beyond it covering 43% of 6533's grid, worth up to 4.1 um at the
    step.  Against the released field the change is worth 4.07 um -> 1.85 um on the worst node.
    """
    xy = np.asarray(xy_um, dtype=np.float64).reshape(-1, 2)
    p, d = fit["points_um"], fit["offsets_um"]
    out = np.empty_like(xy)
    for start in range(0, len(xy), block):
        chunk = xy[start:start + block]
        square = _square_distances(chunk, p)
        weight = np.exp(-square / (2.0 * FIELD_SMOOTH_UM * FIELD_SMOOTH_UM))
        total = weight.sum(axis=1)
        computable = total > 0.0
        out[start:start + block][computable] = (
            (weight[computable] @ d) / total[computable, None])
        if not np.all(computable):
            far = np.flatnonzero(~computable)
            out[start:start + block][far] = d[np.argmin(square[far], axis=1)]
    return out


def _frame_blocks(extent_um: float) -> int:
    """Node intervals covering ``extent_um``, rounded up to a whole ``FRAME_BLOCK_UM`` block."""
    per_block = int(round(FRAME_BLOCK_UM / FIELD_PIXEL_UM))
    blocks = int(np.ceil(float(extent_um) / FRAME_BLOCK_UM))
    return max(1, blocks) * per_block


def _frame_grid(record: RI.RegistrationInput, windows: list[dict]) -> GridSpec:
    """The grid the field is defined on: micron zero out to the far corner of the fixed frame.

    The frame is the IF image the gate measured against -- ``pyramid_shapes[0]`` at the
    inventory's ``pixel_um`` -- which ``run_align`` re-reads from the qptiff and refuses to
    measure against unless it still matches the inventory row, and which ``_verified_source``
    has already bound to this gate.  So the grid is a property of the section's frame, fixed
    before a single window is measured.

    It is emphatically NOT the bounding box of the measured windows.  That box moves whenever a
    window stops reporting a peak, so two passes over the same section would write two
    differently shaped transforms; it stops at the tissue, so a nucleus in the margin samples
    off the end of the field and is silently clamped to its edge value; and it is why a chain
    run end to end produced a 57 x 31 field at (250, 8750) where 6533's released
    ``refined_field``/``refined_field2`` are 97 x 41 at (0, 0).  A field on a different grid is
    not a reproduction of one.

    The windows are still consulted, but only to widen the frame -- never to originate it.  A
    measured window outside the declared image is a contradiction worth refusing, not rounding.
    """
    rows_px, cols_px = (int(v) for v in record.pyramid_shapes[0])
    pixel_um = float(record.pixel_um)
    frame_x_um, frame_y_um = cols_px * pixel_um, rows_px * pixel_um

    q = np.asarray([[w["cx_um"], w["cy_um"]] for w in windows], dtype=np.float64) \
        if windows else np.empty((0, 2))
    if len(q):
        half = RA.WINDOW_UM / 2.0
        if float(q[:, 0].min()) < -half or float(q[:, 1].min()) < -half:
            raise FieldRefineError(
                f"a probe window centred at ({q[:, 0].min():.0f}, {q[:, 1].min():.0f}) um lies "
                f"before the origin of the {frame_x_um:.0f} x {frame_y_um:.0f} um fixed frame; "
                f"the gate and the transform disagree about where the frame starts")
        frame_x_um = max(frame_x_um, float(q[:, 0].max()) + half)
        frame_y_um = max(frame_y_um, float(q[:, 1].max()) + half)

    return GridSpec(0.0, 0.0, FIELD_PIXEL_UM,
                    _frame_blocks(frame_y_um) + 1, _frame_blocks(frame_x_um) + 1)


def _corrected_transform(source: Transform, fit: dict, all_windows: list[dict],
                         *, record: RI.RegistrationInput,
                         iterations: int, source_stem: str,
                         anchor: str = "half-node") -> tuple[Transform, float]:
    # A source that already carries a field fixes the grid: the new field is ADDED to it node
    # for node, so re-gridding here would resample the source's own displacement and lose it.
    if source.displacement is not None:
        grid = source.grid
        old_field = np.asarray(source.displacement, dtype=np.float64)
        if not np.all(np.isfinite(old_field)):
            raise FieldRefineError(
                f"the source field on {grid.height} x {grid.width} nodes is not finite; "
                f"composing onto it would write NaN into the released transform")
    else:
        grid = _frame_grid(record, all_windows)
        old_field = np.zeros((2, grid.height, grid.width), dtype=np.float64)
    u = _node_microns(grid, anchor)
    old = np.column_stack([old_field[0].ravel(), old_field[1].ravel()])
    # Transform.apply computes q = u + D(u), and the gate measured its offsets at q -- the
    # place the source transform currently puts the tissue.  The new field is therefore
    # D_new(u) = D(u) + field(q), which preserves the source exactly and adds the shift that
    # was measured where the tissue actually landed.
    new = old + _predict_on_grid(fit, u + old, grid.shape)
    added = float(np.max(np.linalg.norm(new - old, axis=1)))
    # A weighted mean of the measurements cannot exceed their largest magnitude.  If it does,
    # the fill path or the weights are wrong, not the tissue.
    if added > fit["baseline_max_um"] * (1.0 + 1e-6):
        raise FieldRefineError(
            f"fitted correction reaches {added:.2f} um but the largest offset any window "
            f"measured is {fit['baseline_max_um']:.2f} um; the field is extrapolating")
    field = np.stack([new[:, 0].reshape(grid.shape),
                      new[:, 1].reshape(grid.shape)]).astype(np.float32)
    if not np.all(np.isfinite(field)):
        # Transform.max_displacement_um reads the field with nanmax, so a NaN node would be
        # written, declared and gated without anything saying so.
        raise FieldRefineError(
            f"{int((~np.isfinite(field)).sum())} of {field.size} field values are not finite")
    # The transform's meta carries what the transform IS -- where the field came from, how many
    # windows, its node spacing, which pass -- and nothing else.  The fit's own diagnostics
    # (bandwidth, largest added correction, held-out residual) live in the refinement report
    # beside the numbers they have to be read with; duplicating them here only produces a
    # second copy to disagree with the first, and the released refined_field2 carries neither.
    meta = dict(source.meta)
    meta.update({
        "field_source": f"run_align occupied windows ({source_stem})",
        "n_windows": int(fit["n_windows"]),
        "field_pixel_um": float(grid.pixel_um),
        "field_iterations": int(iterations),
        "field_estimator": fit["estimator"],
        "field_grid_anchor": anchor,
    })
    return Transform(affine=np.array(source.affine), displacement=field, grid=grid,
                     meta=meta), added


def _iteration_count(source: Transform) -> int:
    """Which pass this is over the same affine.

    A source written before ``field_iterations`` was recorded still declares ``field_source``,
    so a chain started by an earlier version of this stage keeps counting rather than resetting.
    """
    prior = source.meta.get("field_iterations")
    if prior is None:
        prior = 1 if source.meta.get("field_source") else 0
    return int(prior) + 1


def refine_one(record: RI.RegistrationInput, donor_dir: Path, *,
               source_stem: str, output_stem: str, overwrite: bool = False,
               estimator: str = "griddata", grid_anchor: str = "half-node",
               selection: str = "measured") -> dict:
    donor_dir = Path(donor_dir)
    source_stem = IR._safe_stem(source_stem, role="source")
    output_stem = IR._safe_stem(output_stem, role="output")
    if source_stem == output_stem:
        raise FieldRefineError("source and output stems must differ")
    output_path = donor_dir / f"{output_stem}.json"
    displacement_out = donor_dir / f"{output_stem}_displacement.npz"
    for path in (output_path, displacement_out):
        if path.exists() and not overwrite:
            raise FieldRefineError(f"{path} exists; pass --overwrite to replace it")

    source, gate, gate_path = _verified_source(donor_dir, source_stem, record=record)
    windows = _measured_windows(gate, selection=selection)
    fit = _fit_field(windows, estimator=estimator)
    iterations = _iteration_count(source)
    corrected, added = _corrected_transform(
        source, fit, gate["alignment"]["windows"], record=record,
        iterations=iterations, source_stem=source_stem, anchor=grid_anchor)
    transform_path, displacement_path = IR._save_transform_atomic(
        corrected, donor_dir, output_stem)

    grid = corrected.grid
    report = {
        "donor_id": gate["donor_id"], "serial": gate["serial"], "roi": gate["roi"],
        "source": {
            "stem": source_stem, "gate_file": str(gate_path),
            "gate_sha256": RA._sha256(gate_path),
            "transform_sha256": gate["alignment"]["transform_sha256"],
            "displacement_sha256": gate["alignment"].get("displacement_sha256"),
            "verdict": gate["alignment"]["verdict"],
            "median_offset_um": gate["alignment"]["median_offset_um"],
            "n_windows_measured": gate["alignment"].get("n_measured"),
            "had_displacement": source.displacement is not None,
        },
        "fit": {
            "family": ("linear barycentric interpolation of per-window DAPI offsets, nearest "
                       "window outside the hull, Gaussian-smoothed on the node grid"
                       if estimator == "griddata"
                       else "Gaussian-kernel regression of per-window DAPI offsets"),
            "estimator": estimator,
            "grid_anchor": grid_anchor,
            "smoothing_um": (GRIDDATA_SMOOTH_NODES * FIELD_PIXEL_UM if estimator == "griddata"
                             else FIELD_SMOOTH_UM),
            "selection": (f"{selection}; search-rim truncation kept"),
            "n_windows": fit["n_windows"],
            "n_windows_measured_in_gate": sum(
                1 for w in gate["alignment"].get("windows", [])
                if w.get("reason") == "measured"),
            "n_truncated": fit["n_truncated"],
            "measured_median_um": fit["baseline_median_um"],
            "measured_p90_um": fit["baseline_p90_um"],
            "measured_max_um": fit["baseline_max_um"],
            "heldout_residual_median_um": fit["heldout_median_um"],
            "heldout_residual_p90_um": fit["heldout_p90_um"],
            "heldout_scheme": "leave-one-window-out",
            "n_heldout": fit["n_heldout"],
            "field_iterations": iterations,
            "grid": grid.as_dict(),
            "grid_from_source": source.displacement is not None,
            "grid_frame": ("inherited from the source field" if source.displacement is not None
                           else "fixed frame from micron zero"),
            "frame_um": [float(record.pyramid_shapes[0][1]) * float(record.pixel_um),
                         float(record.pyramid_shapes[0][0]) * float(record.pixel_um)],
            "frame_block_um": FRAME_BLOCK_UM,
            "max_added_correction_um": added,
        },
        "output": {
            "stem": output_stem, "transform_file": str(transform_path),
            "transform_sha256": RA._sha256(transform_path),
            "displacement_file": str(displacement_path),
            "displacement_sha256": RA._sha256(displacement_path),
            "max_displacement_um": corrected.max_displacement_um,
            "required_gate_command": (
                f"run_align.py --donors {gate['donor_id']} --transform {output_stem} "
                f"--windows occupied"),
        },
        "constants": {
            "ESTIMATOR": estimator,
            "GRID_ANCHOR": grid_anchor,
            "WINDOW_SELECTION": selection,
            "GRIDDATA_SMOOTH_NODES": GRIDDATA_SMOOTH_NODES,
            "FIELD_PIXEL_UM": FIELD_PIXEL_UM,
            "FRAME_BLOCK_UM": FRAME_BLOCK_UM,
            "FIELD_SMOOTH_UM": FIELD_SMOOTH_UM,
            "MIN_NEIGHBOUR_WEIGHT": MIN_NEIGHBOUR_WEIGHT,
            "MIN_FIELD_WINDOWS": MIN_FIELD_WINDOWS,
            "MAX_TRUNCATED_FRACTION": MAX_TRUNCATED_FRACTION,
            "SEARCH_RADIUS_UM": SEARCH_RADIUS_UM,
            "TRUNCATION_RADIUS_UM": TRUNCATION_RADIUS_UM,
            "WINDOW_UM": RA.WINDOW_UM,
            "WINDOW_STEP_UM": RA.WINDOW_STEP_UM,
        },
    }
    report_path = donor_dir / f"field_refinement_{output_stem}.json"
    IR._atomic_json(report_path, report)
    report["report_file"] = str(report_path)
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--source-stem", default="refined")
    ap.add_argument("--output-stem", default="refined_field")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--estimator", choices=ESTIMATORS, default="griddata",
                    help="griddata (the estimator the released 6533 field used) or kernel")
    ap.add_argument("--grid-anchor", choices=GRID_ANCHORS, default="half-node",
                    help="half-node reproduces the released field's 125 um corner offset; "
                         "corner is the convention GridSpec declares")
    ap.add_argument("--window-selection", choices=WINDOW_SELECTIONS, default="measured",
                    help="measured (what the released field used) or measured-and-"
                         "contrast-improving")
    ap.add_argument("--inventory", type=Path, default=RI.REGISTRATION_INPUTS)
    ap.add_argument("--images-dir", type=Path, default=RI.IMAGES)
    ap.add_argument("--geojson-dir", type=Path, default=RI.GEOJSON_DIR)
    ap.add_argument("--xenium-paths", type=Path, default=RI.XENIUM_PATHS)
    ap.add_argument("--xenium-root", type=Path, default=RI.XENIUM_ROOT)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    records = RI.load_registration_inputs(
        args.inventory, args.xenium_paths, args.images_dir, args.geojson_dir,
        args.xenium_root, roi=args.roi)
    selected = RI.select_registration_inputs(records, args.donors,
                                             context="run_field_refine")
    for record in selected:
        RI.require_registration_files(record)
    rows = []
    for record in selected:
        row = refine_one(record, args.out / "valis" / f"donor_id={record.donor_id}",
                         source_stem=args.source_stem, output_stem=args.output_stem,
                         overwrite=args.overwrite, estimator=args.estimator,
                         grid_anchor=args.grid_anchor, selection=args.window_selection)
        rows.append(row)
        f, o = row["fit"], row["output"]
        print(f"[{record.donor_id}] {f['estimator']}/{f['grid_anchor']} "
              f"pass {f['field_iterations']} over "
              f"{f['n_windows']} windows ({f['n_truncated']} at the search rim): measured "
              f"{f['measured_median_um']:.2f} um, held-out "
              f"{f['heldout_residual_median_um']:.2f} um; added up to "
              f"{f['max_added_correction_um']:.2f} um -> {o['stem']} "
              f"(max field {o['max_displacement_um']:.2f} um). Gate it with "
              f"{o['required_gate_command']}", flush=True)
    if not rows:
        print("no donor matched", file=sys.stderr)
        return 1
    summary = args.out / f"field_refinement_{args.roi}_{args.output_stem}.json"
    IR._atomic_json(summary, {"rows": rows})
    print(f"-> {summary}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
