#!/usr/bin/env python
"""The registered re-stain, in the file formats Xenium Explorer imports.

This lane fitted a transform that maps **Xenium microns to IF microns** and measured it at a
median per-cell distance of 0.546 um (6436) and 0.630 um (6442). Xenium Explorer wants the
same correspondence expressed the other way round and in a different unit -- so the one thing
this module must get exactly right is the change of frame, and every step of it is derived
here rather than assumed.

Three sources fix the target convention, and they agree:

* 10x's file-format documentation: the morphology images are **0.2125 um/px** at full
  resolution, "the origin of the coordinate system is the upper left of the TIFF image", and
  micron coordinates "can be converted to pixel coordinates by dividing by the pixel size".
  So Xenium pixel = Xenium micron / 0.2125, with no half-pixel offset.
* 10x's image-alignment documentation: the alignment file is a **3x3 transformation matrix
  (row-major), last row [0,0,1], no header**, and "the transformation value units are in
  pixels".
* ``spatialdata_io.xenium_aligned_image`` reads that CSV with ``header=None`` and applies it
  as ``Affine(alignment, input_axes=("x","y"), output_axes=("x","y"))`` into the same global
  frame it puts Xenium geometry in by scaling microns by ``1/pixel_size``. Which settles the
  direction that prose alone leaves open: the matrix maps **imported-image pixels to Xenium
  morphology pixels**, not the reverse.

The IF side is where a half pixel can hide. ``export_dapi.py`` measured this qptiff's pyramid
to be a **box mean** (r=0.99992 against 0.9891 for the best decimation), so level-``L`` pixel
``k`` is the mean of level-0 pixels ``f*k .. f*k+f-1`` and its centre sits at ``(k + 0.5) * s``
microns for ``s = f * 0.25``. That ``s/2`` is the frame the affine was fitted in, so it belongs
in the matrix. What is *not* settled by any document is whether Explorer indexes an imported
image by pixel corner or pixel centre; the gap that leaves is ``s/2`` = **0.125 um at full
resolution**, which is a fifth of the residual this lane reports and is recorded in the
bundle's manifest rather than quietly assumed away.

Nothing here re-fits anything. It reads the transform ``run_refine.py`` chose and re-expresses
it, so a defect in this module shows up as a constant offset in Explorer rather than as a
number that looks plausible.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

#: 10x Genomics, Xenium file-format documentation: the morphology images are 0.2125 um/px at
#: full resolution and micron coordinates convert to pixels by dividing by it. Not a fit and
#: not negotiable -- the alignment file lands in Explorer's frame or it lands nowhere.
XENIUM_UM_PER_PX = 0.2125

#: 10x's conversion guidance: "tile size should be 1024 px (important for performance in
#: Xenium Explorer; smaller tile sizes may cause laggy performance)". The source qptiff is
#: tiled at 512, so re-tiling is the point of rewriting it rather than an incidental change.
EXPLORER_TILE = 1024

#: The alignment file Xenium Explorer itself writes, byte for byte. 10x's published
#: ``Xenium_V1_hPancreas_Cancer_Add_on_FFPE_he_imagealignment.csv`` (73 bytes, hexdumped) is
#: comma-delimited, has no header, CRLF between rows, NO newline after the last row, a last row
#: that is literally ``0,0,1``, and numbers of at most nine significant digits with trailing
#: zeros stripped -- every one of its nine tokens is what ``'%.9g'`` produces. Two other
#: 10x-published matrices carry 16-17 digit floats, but nothing says Explorer wrote those, so
#: this writer matches the exemplar instead of arguing that a parser ought to cope: on these
#: sections nine digits move the far corner of the Xenium field by under 3e-4 px (6e-5 um).
ALIGNMENT_CSV_EXEMPLAR = "Xenium_V1_hPancreas_Cancer_Add_on_FFPE_he_imagealignment.csv"
ALIGNMENT_SIGNIFICANT_DIGITS = 9

#: Explorer's keypoints file (v1.3+): ``fixed`` is the Xenium morphology (DAPI) image in pixels,
#: ``alignment`` is the imported image in pixels. 10x's documented example and its published
#: keypoints file both carry full-precision floats in every column, so those are written as-is.
KEYPOINTS_HEADER = ("fixedX", "fixedY", "alignmentX", "alignmentY")

#: How far the similarity fitted through the keypoints may sit from the transform's own
#: rotation and scale before the keypoints path is declared inconsistent. A direction error is
#: 0.24 degrees / a factor 1.39 on these sections; the shear a similarity absorbs is ~0.01
#: degrees / ~1e-4.
CONSISTENCY_ROTATION_DEG = 0.05
CONSISTENCY_SCALE_REL = 1e-3

#: A similarity fitted through keypoints that encode an affine misses them by the affine's
#: anisotropy times the field (measured 4e-4 of the extent); through a reflected or otherwise
#: incoherent correspondence it misses by the extent itself. One in a hundred of the extent
#: separates the two by orders of magnitude on both sides.
SIMILARITY_MAX_RESIDUAL_FRACTION = 1e-2

#: The bar ``hormone.py`` was pre-registered at, before any of it was run. Quoted here because
#: this module decides what may be called a hormone-positive cell, and that decision is not
#: allowed to move the bar in either direction.
PREREGISTERED_AUROC_BAR = 0.75


class ExplorerExportError(RuntimeError):
    """A bundle cannot be written as the thing it claims to be -- never silently degraded."""


# --------------------------------------------------------------------------- #
# The change of frame
# --------------------------------------------------------------------------- #

def alignment_matrix(affine: np.ndarray, *, if_um_per_px: float, if_origin_um: float,
                     xenium_um_per_px: float = XENIUM_UM_PER_PX) -> np.ndarray:
    """Xenium Explorer's 3x3 alignment matrix: exported IF pixels -> Xenium morphology pixels.

    ``affine`` is this lane's own 2x3, mapping **Xenium microns to IF microns**. Explorer
    needs the inverse of that, sandwiched between the two pixel/micron scalings::

        S_if   = [[s, 0, o], [0, s, o], [0, 0, 1]]   exported IF pixel -> IF micron
        A      = [[a, b, tx], [c, d, ty], [0, 0, 1]] Xenium micron     -> IF micron
        S_xen  = [[k, 0, 0], [0, k, 0], [0, 0, 1]]   Xenium pixel      -> Xenium micron

        M = S_xen^-1 @ A^-1 @ S_if

    ``o`` (``if_origin_um``) is not cosmetic. It is where exported pixel coordinate 0 sits in
    the micron frame the affine was fitted in, and ``export_dapi.origin_um`` derives it from
    the *measured* pyramid rule. Dropping it puts every cell a constant ``o`` off in both axes.
    """
    A = np.asarray(affine, dtype=np.float64)
    if A.shape != (2, 3):
        raise ExplorerExportError(f"expected a 2x3 micron affine, got shape {A.shape}")
    if not np.isfinite(A).all():
        raise ExplorerExportError("the affine carries non-finite entries")
    if abs(float(np.linalg.det(A[:, :2]))) < 1e-12:
        raise ExplorerExportError("the affine's linear part is singular and cannot be inverted")
    for name, v in (("if_um_per_px", if_um_per_px), ("xenium_um_per_px", xenium_um_per_px)):
        if not float(v) > 0.0:
            raise ExplorerExportError(f"{name} must be positive, got {v!r}")

    A_h = np.vstack([A, [0.0, 0.0, 1.0]])
    s, o, k = float(if_um_per_px), float(if_origin_um), float(xenium_um_per_px)
    S_if = np.array([[s, 0.0, o], [0.0, s, o], [0.0, 0.0, 1.0]])
    S_xen_inv = np.array([[1.0 / k, 0.0, 0.0], [0.0, 1.0 / k, 0.0], [0.0, 0.0, 1.0]])

    M = S_xen_inv @ np.linalg.inv(A_h) @ S_if
    M[2] = [0.0, 0.0, 1.0]          # exact, rather than 1e-17 away from it
    return M


def matrix_round_trip_um(M: np.ndarray, affine: np.ndarray, *, if_um_per_px: float,
                         if_origin_um: float, xenium_um_per_px: float = XENIUM_UM_PER_PX,
                         xy_xenium_um: np.ndarray) -> float:
    """Worst-case microns by which ``M`` disagrees with ``affine`` on the given points.

    The matrix is a re-expression, so this has no fitting freedom in it: it is a check that
    the composition above was assembled correctly, and it is written into the manifest so a
    reader does not have to take the algebra on trust. Points are pushed Xenium micron ->
    IF micron -> IF pixel -> (via M) Xenium pixel -> Xenium micron and compared with where
    they started.
    """
    xy = np.asarray(xy_xenium_um, dtype=np.float64).reshape(-1, 2)
    A = np.asarray(affine, dtype=np.float64)
    if_um = xy @ A[:, :2].T + A[:, 2]
    if_px = (if_um - float(if_origin_um)) / float(if_um_per_px)
    xen_px = np.column_stack([if_px, np.ones(len(if_px))]) @ np.asarray(M).T
    back = xen_px[:, :2] * float(xenium_um_per_px)
    return float(np.max(np.linalg.norm(back - xy, axis=1)))


def format_alignment_number(v: float) -> str:
    """One matrix entry the way the exemplar writes it: <= 9 significant digits, positional.

    ``'%.9g'`` reproduces every token of the exemplar, but switches to an exponent below 1e-4,
    which the exemplar never shows; ``format_float_positional`` gives the same digits without
    ever doing that. Negative zero is written as ``0``.
    """
    x = float(v)
    if not np.isfinite(x):
        raise ExplorerExportError(f"cannot write a non-finite matrix entry ({x!r})")
    if x == 0.0:
        return "0"
    s = np.format_float_positional(x, precision=ALIGNMENT_SIGNIFICANT_DIGITS, unique=False,
                                   fractional=False, trim="-")
    if s.endswith("."):
        s = s[:-1]
    return s


def rounded_alignment_matrix(M: np.ndarray) -> np.ndarray:
    """The matrix as it reads back from the file this module writes -- what actually ships."""
    M = np.asarray(M, dtype=np.float64)
    out = np.array([[float(format_alignment_number(v)) for v in M[r]] for r in (0, 1)]
                   + [[0.0, 0.0, 1.0]])
    return out


def write_alignment_csv(path: Path, M: np.ndarray) -> Path:
    """The matrix in the byte layout of Explorer's own export (see ``ALIGNMENT_CSV_EXEMPLAR``).

    Three rows, no header, CRLF between rows and **no newline after the last row**, the last
    row literally ``0,0,1``, the other entries at most nine significant digits. Nothing is left
    to a library default: the previous version of this function used ``csv.writer`` and
    ``repr``, which added a trailing CRLF and wrote ``0.0,0.0,1.0`` with 17-digit floats, and
    Explorer 4.x refused that file with "Failed to load file." Whether any of those three
    differences was the reason is unknown -- 10x documents no validation rule -- so all three
    are removed. The rounding is worth under 3e-4 Xenium px at the far corner of the field.
    """
    M = np.asarray(M, dtype=np.float64)
    if M.shape != (3, 3):
        raise ExplorerExportError(f"an alignment matrix is 3x3, got {M.shape}")
    if not np.allclose(M[2], [0.0, 0.0, 1.0]):
        raise ExplorerExportError(f"the last row must be [0,0,1], got {M[2].tolist()}")
    rows = [",".join(format_alignment_number(v) for v in M[r]) for r in (0, 1)] + ["0,0,1"]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes("\r\n".join(rows).encode("ascii"))
    return path


def read_alignment_csv(path: Path) -> np.ndarray:
    """The file back as a 3x3, tolerating either line ending -- how a reader sees it."""
    raw = Path(path).read_bytes().decode("ascii")
    rows = [r for r in raw.replace("\r\n", "\n").split("\n") if r.strip()]
    if len(rows) != 3:
        raise ExplorerExportError(f"{path}: expected 3 rows, found {len(rows)}")
    cells = [r.split(",") for r in rows]
    if any(len(c) != 3 for c in cells):
        raise ExplorerExportError(f"{path}: expected 3 comma-separated numbers per row, got "
                                  f"{[len(c) for c in cells]}")
    try:
        return np.array([[float(v) for v in c] for c in cells], dtype=np.float64)
    except ValueError as e:
        raise ExplorerExportError(f"{path}: a cell is not a number ({e})") from None


def alignment_csv_format() -> dict:
    """What the bytes are, for the manifest, so a report of a refused file can be matched."""
    return {"delimiter": ",", "header": False, "row_terminator": "CRLF",
            "trailing_newline": False, "last_row": "0,0,1",
            "number_format": f"positional, <= {ALIGNMENT_SIGNIFICANT_DIGITS} significant digits, "
                             "trailing zeros stripped (the '%.9g' rule)",
            "exemplar": ALIGNMENT_CSV_EXEMPLAR}


# --------------------------------------------------------------------------- #
# The keypoints, and the similarity Explorer would fit through them
# --------------------------------------------------------------------------- #

def keypoint_pairs(cells_xy_um: np.ndarray, M: np.ndarray, *, morph_hw: tuple[int, int],
                   if_hw: tuple[int, int], grid: tuple[int, int] = (4, 5),
                   cell_ids=None, xenium_um_per_px: float = XENIUM_UM_PER_PX) -> dict:
    """Exact correspondences under ``M``, anchored on real Xenium nuclei.

    Explorer's second import vehicle is a keypoints file it fits its own similarity through.
    These are not a fit: each ``fixed`` point is a Xenium cell centroid (in morphology pixels)
    and its ``alignment`` point is exactly ``M^-1`` of it, so the pairs encode the shipped
    matrix and nothing else. They sit on nuclei rather than on a bare grid so that, in
    Explorer's keypoint view, each one can be seen landing on the same nucleus in both images.
    A pair outside either image is refused -- that is the one thing Explorer can legitimately
    reject.
    """
    from scipy.spatial import cKDTree

    xy = np.asarray(cells_xy_um, dtype=np.float64).reshape(-1, 2)
    M = np.asarray(M, dtype=np.float64)
    if M.shape != (3, 3) or not np.allclose(M[2], [0.0, 0.0, 1.0]):
        raise ExplorerExportError("keypoints need a 3x3 alignment matrix with last row [0,0,1]")
    nx, ny = int(grid[0]), int(grid[1])
    if nx * ny < 3 or len(xy) < 3:
        raise ExplorerExportError("Explorer needs at least 3 keypoints")
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    gx = lo[0] + (hi[0] - lo[0]) * (np.arange(nx) + 0.5) / nx
    gy = lo[1] + (hi[1] - lo[1]) * (np.arange(ny) + 0.5) / ny
    nodes = np.array([(x, y) for y in gy for x in gx], dtype=np.float64)
    _, idx = cKDTree(xy).query(nodes)
    idx = np.array(sorted(set(int(i) for i in idx)))
    if len(idx) < 3:
        raise ExplorerExportError(f"only {len(idx)} distinct nuclei under the {nx}x{ny} grid")

    fixed = xy[idx] / float(xenium_um_per_px)
    alignment = (np.column_stack([fixed, np.ones(len(fixed))]) @ np.linalg.inv(M).T)[:, :2]
    H_m, W_m = int(morph_hw[0]), int(morph_hw[1])
    H_i, W_i = int(if_hw[0]), int(if_hw[1])
    off_m = ((fixed[:, 0] < 0) | (fixed[:, 0] >= W_m) | (fixed[:, 1] < 0) | (fixed[:, 1] >= H_m))
    # M's IF input is centre-indexed (pixel p is centred at coordinate p), so the image spans
    # [-0.5, W-0.5] x [-0.5, H-0.5] in that frame; the morphology side is corner-origin.
    off_i = ((alignment[:, 0] < -0.5) | (alignment[:, 0] > W_i - 0.5)
             | (alignment[:, 1] < -0.5) | (alignment[:, 1] > H_i - 0.5))
    if off_m.any() or off_i.any():
        raise ExplorerExportError(
            f"{int(off_m.sum())} keypoint(s) fall outside the morphology image and "
            f"{int(off_i.sum())} outside the imported image; refusing to write them")
    kp = np.column_stack([fixed, alignment])
    return {"keypoints": kp, "index": idx, "grid": [nx, ny],
            "cell_ids": ([str(cell_ids[i]) for i in idx] if cell_ids is not None else None)}


def write_keypoints_csv(path: Path, kp: np.ndarray) -> Path:
    """``fixedX,fixedY,alignmentX,alignmentY`` with full-precision floats, CRLF, no trailing newline.

    The header and the float columns follow 10x's own keypoints exemplars; whether Explorer
    minds a trailing newline on THIS file is unverified, so the matrix file's convention (none)
    is used and the choice is recorded in the manifest.
    """
    kp = np.asarray(kp, dtype=np.float64)
    if kp.ndim != 2 or kp.shape[1] != 4 or kp.shape[0] < 3:
        raise ExplorerExportError(f"keypoints must be (n>=3, 4), got {kp.shape}")
    if not np.isfinite(kp).all():
        raise ExplorerExportError("keypoints carry non-finite entries")
    rows = [",".join(KEYPOINTS_HEADER)]
    rows += [",".join(repr(float(v)) for v in row) for row in kp]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes("\r\n".join(rows).encode("ascii"))
    return path


def read_keypoints_csv(path: Path) -> np.ndarray:
    raw = Path(path).read_bytes().decode("ascii")
    rows = [r for r in raw.replace("\r\n", "\n").split("\n") if r.strip()]
    if not rows or tuple(rows[0].split(",")) != KEYPOINTS_HEADER:
        raise ExplorerExportError(f"{path}: header is not {','.join(KEYPOINTS_HEADER)}")
    cells = [r.split(",") for r in rows[1:]]
    if not cells or any(len(c) != 4 for c in cells):
        raise ExplorerExportError(f"{path}: expected 4 comma-separated numbers per row, got "
                                  f"{[len(c) for c in cells]}")
    try:
        return np.array([[float(v) for v in c] for c in cells], dtype=np.float64)
    except ValueError as e:
        raise ExplorerExportError(f"{path}: a cell is not a number ({e})") from None


def keypoints_csv_format() -> dict:
    return {"delimiter": ",", "header": list(KEYPOINTS_HEADER), "row_terminator": "CRLF",
            "trailing_newline": False, "number_format": "repr (shortest round-trip float)",
            "fixed": "Xenium morphology image pixels (cell centroid um / 0.2125)",
            "alignment": "imported IF image pixels = M^-1 @ fixed"}


def keypoints_residual_px(kp: np.ndarray, M: np.ndarray) -> float:
    """Max |M @ alignment - fixed| in Xenium pixels: ~1e-9 by construction, or something is wrong."""
    kp = np.asarray(kp, dtype=np.float64)
    pushed = (np.column_stack([kp[:, 2:4], np.ones(len(kp))]) @ np.asarray(M).T)[:, :2]
    return float(np.max(np.linalg.norm(pushed - kp[:, 0:2], axis=1)))


def is_exact_similarity(M: np.ndarray) -> bool:
    M = np.asarray(M, dtype=np.float64)
    return bool(M[0, 0] == M[1, 1] and M[0, 1] == -M[1, 0])


def similarity_from_keypoints(kp: np.ndarray) -> np.ndarray:
    """The least-squares similarity through the keypoints, as Explorer's own aligner would fit.

    10x documents Explorer's alignment as a similarity (Horn 1987: translation, uniform scale,
    rotation, reflection); ``fit_rigid_umeyama`` is the same closed form. The linear part is
    rebuilt from the fitted (scale, angle) so ``a == d`` and ``b == -c`` hold **bit-exactly**,
    which is what every Explorer-exported matrix shows, and the translation is recomputed
    against that canonical rotation. This is a projection of the shipped matrix at chosen
    points -- a diagnostic with a recorded cost, never a registration.
    """
    from phenocycler.integration.transform import fit_rigid_umeyama

    kp = np.asarray(kp, dtype=np.float64)
    fixed, alignment = kp[:, 0:2], kp[:, 2:4]
    T = fit_rigid_umeyama(alignment, fixed, allow_scale=True, allow_reflection=False)
    A = T.affine
    s = float(np.sqrt(abs(np.linalg.det(A[:, :2]))))
    theta = float(np.arctan2(A[1, 0], A[0, 0]))
    c, si = s * np.cos(theta), s * np.sin(theta)
    L = np.array([[c, -si], [si, c]], dtype=np.float64)
    t = fixed.mean(axis=0) - L @ alignment.mean(axis=0)
    S = np.eye(3)
    S[:2, :2] = L
    S[:2, 2] = t
    if not s > 0.0:
        raise ExplorerExportError("the similarity through the keypoints has zero scale")
    # ``det(L) = s^2 > 0`` and ``a == d, b == -c`` hold by construction, so they cannot be the
    # check. What can fail is the fit itself: a proper-rotation fit through a reflected or
    # incoherent correspondence still returns a matrix, one that misses the keypoints by the
    # size of the field. Refuse that here rather than write it as a "diagnostic".
    extent = float(np.ptp(fixed, axis=0).max())
    resid = keypoints_residual_px(kp, S)
    if extent > 0 and resid > SIMILARITY_MAX_RESIDUAL_FRACTION * extent:
        raise ExplorerExportError(
            f"the similarity misses the keypoints by {resid:.1f} px over a {extent:.0f} px "
            f"field; the correspondence is not a similarity plus a small affine residual "
            f"(reflected? swapped columns?)")
    return S


def similarity_deviation_um(M: np.ndarray, S: np.ndarray, pts_if_px: np.ndarray,
                            *, xenium_um_per_px: float = XENIUM_UM_PER_PX) -> np.ndarray:
    """Per-point |M p - S p| in microns, for imported-image pixel points ``p``."""
    p = np.asarray(pts_if_px, dtype=np.float64).reshape(-1, 2)
    h = np.column_stack([p, np.ones(len(p))])
    d = (h @ np.asarray(M).T - h @ np.asarray(S).T)[:, :2]
    return np.linalg.norm(d, axis=1) * float(xenium_um_per_px)


def matrix_anisotropy(M: np.ndarray) -> dict:
    """Why ``M`` is not a similarity, in numbers a reader can check."""
    L = np.asarray(M, dtype=np.float64)[:2, :2]
    sv = np.linalg.svd(L, compute_uv=False)
    return {"a_minus_d": float(L[0, 0] - L[1, 1]), "b_plus_c": float(L[0, 1] + L[1, 0]),
            "singular_values": [float(v) for v in sv],
            "axis_ratio": float(sv.max() / sv.min()),
            "is_exact_similarity": is_exact_similarity(M),
            "mirrored": bool(np.linalg.det(L) < 0)}


def similarity_consistency(S: np.ndarray, affine: np.ndarray, *, if_um_per_px: float,
                           xenium_um_per_px: float = XENIUM_UM_PER_PX) -> dict:
    """The one free check the exact keypoints residual cannot make.

    ``affine`` maps Xenium um -> IF um; ``S`` maps IF px -> Xenium px. So the rotation of ``S``
    must be the negative of the affine's and its scale must be ``(if_um_per_px / k) / scale(A)``.
    A sign or direction error in the keypoints path shows up here as 0.24 degrees or a factor
    1.39, where the keypoints residual -- exact by construction -- shows nothing. The tolerances
    are wide enough for what a least-squares similarity through a sheared affine legitimately
    absorbs (measured 0.008 degrees and 7e-5 on donor 6436) and a few times too narrow for
    either failure.
    """
    A = np.asarray(affine, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)
    rot_A = float(np.degrees(np.arctan2(A[1, 0], A[0, 0])))
    scale_A = float(np.sqrt(abs(np.linalg.det(A[:, :2]))))
    rot_S = float(np.degrees(np.arctan2(S[1, 0], S[0, 0])))
    scale_S = float(np.sqrt(abs(np.linalg.det(S[:2, :2]))))
    expected_rot = -rot_A
    expected_scale = (float(if_um_per_px) / float(xenium_um_per_px)) / scale_A
    return {"rotation_deg": rot_S, "expected_rotation_deg": expected_rot,
            "rotation_diff_deg": rot_S - expected_rot,
            "scale": scale_S, "expected_scale": expected_scale,
            "scale_rel_diff": (scale_S - expected_scale) / expected_scale,
            "consistent": bool(abs(rot_S - expected_rot) < CONSISTENCY_ROTATION_DEG
                               and abs(scale_S - expected_scale) / expected_scale
                               < CONSISTENCY_SCALE_REL)}


# --------------------------------------------------------------------------- #
# The morphology grid Explorer places everything on
# --------------------------------------------------------------------------- #

def morphology_grid(bundle: Path, *, xenium_um_per_px: float = XENIUM_UM_PER_PX) -> dict:
    """Level-0 size and pyramid of the Xenium morphology image, read from the bundle itself.

    Opened with ``is_ome=False``: the OME-XML declares four files of which only ``_0000`` is
    on this machine, and following it logs three failures for nothing. Page 0 and its SubIFDs
    are the whole pyramid. Refused unless the pixel size is Xenium's and every level is exactly
    floor-halved -- the warped image is written onto this grid, level for level.
    """
    import re
    import tifffile

    from registration_inputs import morphology_dapi_path
    p = morphology_dapi_path(bundle)
    if not p.exists():
        raise ExplorerExportError(f"{p} is absent")
    with tifffile.TiffFile(str(p), is_ome=False) as tif:
        page = tif.pages[0]
        shapes = [(int(lv.shape[-2]), int(lv.shape[-1])) for lv in tif.series[0].levels]
        m = re.search(r'PhysicalSizeX="([^"]+)"', page.description or "")
    if m is None:
        raise ExplorerExportError(f"{p} carries no PhysicalSizeX")
    um = float(m.group(1))
    if abs(um - float(xenium_um_per_px)) > 1e-12:
        raise ExplorerExportError(f"{p} is {um} um/px, not Xenium's {xenium_um_per_px}")
    H, W = shapes[0]
    expect = [(H >> L, W >> L) for L in range(len(shapes))]
    if shapes != expect:
        raise ExplorerExportError(f"{p} pyramid {shapes} is not floor-halved {expect}")
    return {"file": str(p), "height": H, "width": W, "n_levels": len(shapes),
            "level_shapes": [list(s) for s in shapes], "um_per_px": um}


# --------------------------------------------------------------------------- #
# The image
# --------------------------------------------------------------------------- #

def qptiff_level_shapes(qptiff: Path) -> list[tuple[int, int]]:
    """``(height, width)`` of every pyramid level of the qptiff's first series."""
    import tifffile

    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        return [(int(lv.shape[-2]), int(lv.shape[-1]))
                for lv in getattr(series, "levels", [series])]


def integer_factors(shapes: list[tuple[int, int]]) -> list[int]:
    """Downsample factor of each level, refused unless every one is an exact integer.

    ``shapes[0][0] / shapes[lvl][0]`` is not the factor -- the levels are floor-halved, so the
    ratio can be 4.000176 where the factor is 4, a scale error worth 0.6 um at the far corner
    of a section. ``export_dapi.pyramid_factor`` makes the same point for the DAPI export and
    this is the same rule applied to every level at once, because the exported OME-TIFF
    carries all of them and Explorer will navigate between them.
    """
    import export_dapi as ED

    return [ED.pyramid_factor(shapes, lvl) for lvl in range(len(shapes))]


def _tile_stream(zstore, channels: list[int], height: int, width: int, dtype,
                 tile: int = EXPLORER_TILE):
    """Tiles in the order ``TiffWriter`` consumes them: per channel, per row, per column.

    Streamed rather than materialised. A full-resolution channel of these sections is 2.9 GB,
    and the point of writing tiles is that the peak is one tile rather than the slide.
    """
    for c in channels:
        for y0 in range(0, height, tile):
            for x0 in range(0, width, tile):
                patch = np.asarray(zstore[c, y0:y0 + tile, x0:x0 + tile])
                if patch.shape != (tile, tile):
                    pad = np.zeros((tile, tile), dtype=dtype)
                    pad[:patch.shape[0], :patch.shape[1]] = patch
                    patch = pad
                yield patch


def write_explorer_ome_tiff(qptiff: Path, out_path: Path, channels: list[int],
                            channel_names: list[str], *, base_um_per_px: float,
                            compression: str = "zlib", tile: int = EXPLORER_TILE,
                            max_levels: int | None = None, workers: int = 8) -> dict:
    """The selected qptiff channels as a pyramidal, tiled OME-TIFF Explorer can import.

    10x's requirement is "pyramidal, tiled, OME-TIFF"; this writes exactly that, and takes
    each pyramid level **from the qptiff's own levels** rather than downsampling. That is not
    an optimisation. ``export_dapi.py`` measured this pyramid to be a box mean, and the
    alignment matrix is derived from that same rule -- so re-deriving the levels with some
    other filter would put the matrix and the image on subtly different grids.

    ``compression`` defaults to ``zlib``: lossless, universally readable, and it does not put
    a lossy codec between the colleague and a hormone channel they are being asked to judge.
    10x's own morphology files are JPEG-2000, which ``imagecodecs`` can also write; it is a
    flag rather than the default for that reason.
    """
    import tifffile

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if len(channels) != len(channel_names):
        raise ExplorerExportError("every exported channel needs a name")
    if not channels:
        raise ExplorerExportError("no channels selected; refusing to write an empty image")

    shapes = qptiff_level_shapes(qptiff)
    factors = integer_factors(shapes)
    n_levels = len(shapes) if max_levels is None else min(len(shapes), int(max_levels))

    written = []
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        dtype = series.dtype
        with tifffile.TiffWriter(str(out_path), bigtiff=True) as tw:
            for lvl in range(n_levels):
                h, w = shapes[lvl]
                store = series.aszarr(level=lvl)
                try:
                    import zarr
                    z = zarr.open(store, mode="r")
                    shape = (len(channels), h, w)
                    # Compression is the wall clock here -- a full-resolution channel is
                    # 2.9 GB and there are three of them -- and it is per-tile, so it
                    # parallelises exactly. The tile ORDER is still serial, so the file is
                    # byte-identical whatever this is set to.
                    opts = dict(tile=(tile, tile), compression=compression,
                                photometric="minisblack", maxworkers=int(workers))
                    if lvl == 0:
                        um = float(base_um_per_px)
                        tw.write(_tile_stream(z, channels, h, w, dtype, tile=tile),
                                 shape=shape, dtype=dtype, subifds=n_levels - 1,
                                 metadata={"axes": "CYX",
                                           "Channel": {"Name": list(channel_names)},
                                           "PhysicalSizeX": um, "PhysicalSizeY": um,
                                           "PhysicalSizeXUnit": "µm",
                                           "PhysicalSizeYUnit": "µm"},
                                 **opts)
                    else:
                        tw.write(_tile_stream(z, channels, h, w, dtype, tile=tile),
                                 shape=shape, dtype=dtype, subfiletype=1, **opts)
                finally:
                    store.close()
                written.append({"level": lvl, "height": h, "width": w,
                                "factor": int(factors[lvl]),
                                "um_per_px": float(base_um_per_px) * int(factors[lvl])})
                print(f"[explorer]   level {lvl}: {h} x {w} x {len(channels)}", flush=True)

    return {"path": str(out_path), "dtype": str(dtype), "tile": int(tile),
            "compression": compression, "channels": list(channel_names),
            "channel_indices": [int(c) for c in channels], "levels": written,
            "bytes": int(out_path.stat().st_size)}


# --------------------------------------------------------------------------- #
# The image on the Xenium grid, for "Yes, skip alignment"
# --------------------------------------------------------------------------- #
#
# 10x documents "Yes, skip alignment" for an image "aligned and transformed into the Xenium
# coordinate system with community-developed software", and its reference file for that option
# has exactly the morphology image's dimensions at 0.2125 um/px: Explorer places such an image
# with the identity in pixel space and ignores OME PhysicalSize. So the IF is resampled onto the
# morphology's own level-0 grid through the same matrix, and needs no CSV at all.
#
# Conventions, each pinned by a test. Xenium pixels are corner-origin (10x: px = um / 0.2125),
# so pixel (k, l) is centred at M-frame coordinate (k + 0.5, l + 0.5). The matrix's IF input is
# centre-indexed (S_if carries the measured 0.125 um), and cv2 treats integer coordinates as
# pixel centres, so no correction is needed on the IF side: output pixel (k, l) samples the IF
# at M^-1 (k + 0.5, l + 0.5).

def warp_tile_matrix(Minv: np.ndarray, x0: int, y0: int, cx0: int, cy0: int) -> np.ndarray:
    """cv2's 2x3 for one output tile: output-local (i, j) -> crop-local IF coordinate."""
    Minv = np.asarray(Minv, dtype=np.float64)
    L, t = Minv[:2, :2], Minv[:2, 2]
    A = np.empty((2, 3), dtype=np.float64)
    A[:, :2] = L
    A[:, 2] = L @ np.array([x0 + 0.5, y0 + 0.5]) + t - np.array([cx0, cy0], dtype=np.float64)
    return A


def warp_plane_to_xenium_grid(plane: np.ndarray, M: np.ndarray, height: int, width: int, *,
                              tile: int = EXPLORER_TILE, margin: int = 3,
                              out: np.ndarray | None = None) -> np.ndarray:
    """One channel resampled onto a ``height x width`` Xenium grid, tile by tile.

    Per output tile the corner extents are pushed through ``M^-1``, that crop of the source
    (plus a margin) is warped with a crop-local matrix, and a tile whose source falls entirely
    off the image is left at zero -- outside the IF is 0, never an edge replicated. Bilinear on
    the stored values, no intensity rescaling. cv2 parallelises inside each call, which is why
    the loop is a plain serial one.
    """
    import cv2

    plane = np.asarray(plane)
    if plane.ndim != 2:
        raise ExplorerExportError(f"a plane is 2-D, got shape {plane.shape}")
    M = np.asarray(M, dtype=np.float64)
    Minv = np.linalg.inv(M)
    H_if, W_if = plane.shape
    if out is None:
        out = np.zeros((int(height), int(width)), dtype=plane.dtype)
    elif out.shape != (int(height), int(width)) or out.dtype != plane.dtype:
        raise ExplorerExportError("`out` must match the requested grid and the plane's dtype")
    flags = cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP
    for y0 in range(0, int(height), tile):
        th = min(tile, int(height) - y0)
        for x0 in range(0, int(width), tile):
            tw = min(tile, int(width) - x0)
            corners = np.array([[x0, y0, 1.0], [x0 + tw, y0, 1.0],
                                [x0, y0 + th, 1.0], [x0 + tw, y0 + th, 1.0]])
            s = corners @ Minv.T
            cx0 = max(int(np.floor(s[:, 0].min())) - margin, 0)
            cy0 = max(int(np.floor(s[:, 1].min())) - margin, 0)
            cx1 = min(int(np.ceil(s[:, 0].max())) + margin + 1, W_if)
            cy1 = min(int(np.ceil(s[:, 1].max())) + margin + 1, H_if)
            if cx1 <= cx0 or cy1 <= cy0:
                out[y0:y0 + th, x0:x0 + tw] = 0          # written, not skipped: `out` may be reused
                continue
            crop = np.ascontiguousarray(plane[cy0:cy1, cx0:cx1])
            A = warp_tile_matrix(Minv, x0, y0, cx0, cy0)
            out[y0:y0 + th, x0:x0 + tw] = cv2.warpAffine(
                crop, A, (tw, th), flags=flags, borderMode=cv2.BORDER_CONSTANT, borderValue=0)
    return out


def halve_box_mean(stack: np.ndarray) -> np.ndarray:
    """``(C, h, w) -> (C, h//2, w//2)`` by 2x2 box mean, floor-halved, round half up for ints."""
    stack = np.asarray(stack)
    if stack.ndim != 3:
        raise ExplorerExportError(f"a stack is (C, h, w), got shape {stack.shape}")
    C, h, w = stack.shape
    h1, w1 = h >> 1, w >> 1
    if h1 < 1 or w1 < 1:
        raise ExplorerExportError(f"cannot halve a {h}x{w} level")
    blocks = stack[:, :2 * h1, :2 * w1].reshape(C, h1, 2, w1, 2)
    if np.issubdtype(stack.dtype, np.integer):
        # Four values of itemsize n fit in itemsize 2n: uint16 for uint8, so the level-1
        # transient of a (3, 68119, 37187) stack is 3.8 GB rather than the 46 GB of uint64.
        acc = {1: np.uint16, 2: np.uint32}.get(stack.dtype.itemsize, np.uint64)
        s = blocks.sum(axis=(2, 4), dtype=acc)
        s += 2
        s >>= 2
        return s.astype(stack.dtype)
    return blocks.mean(axis=(2, 4), dtype=np.float64).astype(stack.dtype)


def xenium_pyramid(level0: np.ndarray, n_levels: int) -> list[np.ndarray]:
    """``n_levels`` levels from level 0, each an exact floor-halved box mean of the previous.

    Level-L pixel k covers level-0 pixels ``2^L k .. 2^L (k+1) - 1``, so its centre sits at
    ``0.2125 * 2^L * (k + 0.5)`` um. (10x's own morphology pyramid is decimated rather than
    averaged, which differs from this by under half a level-L pixel and only at coarse levels;
    level 0, the one that matters, is exact.)
    """
    level0 = np.asarray(level0)
    if level0.ndim != 3:
        raise ExplorerExportError(f"level 0 is (C, H, W), got shape {level0.shape}")
    C, H, W = level0.shape
    levels = [level0]
    for L in range(1, int(n_levels)):
        levels.append(halve_box_mean(levels[-1]))
        if levels[-1].shape != (C, H >> L, W >> L):
            raise ExplorerExportError(f"level {L} is {levels[-1].shape}, not {(C, H >> L, W >> L)}")
    return levels


def write_warped_ome_tiff(levels: list[np.ndarray], out_path: Path, channel_names: list[str], *,
                          um_per_px: float = XENIUM_UM_PER_PX, compression: str = "zlib",
                          tile: int = EXPLORER_TILE, workers: int = 8) -> dict:
    """The in-RAM pyramid as a tiled, SubIFD-pyramidal, BigTIFF OME-TIFF on the Xenium grid.

    ``ome=True`` is passed explicitly rather than inferred from the file name. PhysicalSize is
    Xenium's, and the TIFF resolution tags are written as well (10x's Explorer-exported aligned
    images carry both). Reuses ``_tile_stream`` over the in-RAM arrays, so the tile order and
    padding are the same as for the native file.
    """
    import tifffile

    if not levels:
        raise ExplorerExportError("no pyramid levels to write")
    C, H, W = levels[0].shape
    if len(channel_names) != C:
        raise ExplorerExportError(f"{C} channels but {len(channel_names)} names")
    for L, arr in enumerate(levels):
        if arr.ndim != 3 or arr.shape != (C, H >> L, W >> L):
            raise ExplorerExportError(f"level {L} has shape {arr.shape}, expected {(C, H >> L, W >> L)}")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    written = []
    with tifffile.TiffWriter(str(out_path), bigtiff=True, ome=True) as tw:
        for L, arr in enumerate(levels):
            _, h, w = arr.shape
            um = float(um_per_px) * (2 ** L)
            opts = dict(tile=(tile, tile), compression=compression, photometric="minisblack",
                        maxworkers=int(workers), resolution=(1e4 / um, 1e4 / um),
                        resolutionunit="CENTIMETER")
            stream = _tile_stream(arr, list(range(C)), h, w, arr.dtype, tile=tile)
            if L == 0:
                tw.write(stream, shape=(C, h, w), dtype=arr.dtype, subifds=len(levels) - 1,
                         metadata={"axes": "CYX", "Channel": {"Name": list(channel_names)},
                                   "PhysicalSizeX": float(um_per_px),
                                   "PhysicalSizeY": float(um_per_px),
                                   "PhysicalSizeXUnit": "µm", "PhysicalSizeYUnit": "µm"},
                         **opts)
            else:
                tw.write(stream, shape=(C, h, w), dtype=arr.dtype, subfiletype=1, **opts)
            written.append({"level": L, "height": h, "width": w, "factor": 2 ** L,
                            "um_per_px": um})
    return {"path": str(out_path), "dtype": str(levels[0].dtype), "tile": int(tile),
            "compression": compression, "channels": list(channel_names), "levels": written,
            "resolution_unit": "CENTIMETER", "pixels_per_cm": 1e4 / float(um_per_px),
            "bytes": int(out_path.stat().st_size)}


# --------------------------------------------------------------------------- #
# The annotations
# --------------------------------------------------------------------------- #

#: QuPath's own packed-ARGB integers, so the classes arrive in Explorer coloured the way the
#: maintainer drew them rather than in whatever order an importer assigns.
QUPATH_COLOR = {"Tissue": -16744448, "Islet": -3670016}


def warp_geoms_to_xenium(geoms, transform):
    """IF-micron shapely geometries -> Xenium-micron shapely geometries.

    ``inverse_apply`` rather than ``apply``: this lane's transform is fitted Xenium -> IF, and
    the annotations were drawn on the IF. Getting this backwards does not error -- it produces
    polygons roughly 2 * (714, 282) um from where they belong on 6436, which looks like a
    registration failure rather than like a sign error.
    """
    import shapely

    return [shapely.transform(g, lambda a: transform.inverse_apply(a)) for g in geoms]


def _rings(geom) -> list[list[list[float]]]:
    """GeoJSON coordinate rings for a Polygon, exterior first."""
    import shapely

    if geom.geom_type != "Polygon":
        raise ExplorerExportError(f"expected a Polygon, got {geom.geom_type}")
    out = [[[float(x), float(y)] for x, y in shapely.get_coordinates(geom.exterior)]]
    for hole in geom.interiors:
        out.append([[float(x), float(y)] for x, y in shapely.get_coordinates(hole)])
    return out


def annotations_geojson(named_geoms: dict[str, list], path: Path) -> dict:
    """Class name -> geometries, written as a QuPath-shaped GeoJSON FeatureCollection.

    Coordinates are **Xenium microns**. Explorer's annotation import asks which of microns or
    pixels a file is in, so the frame is stated in the filename as well as here; a file that
    is silently in the wrong one imports without complaint and lands 4.7x off.

    ``MultiPolygon`` is split into one feature per part. Explorer's importer takes polygons,
    and the ``Tissue`` class arrives from QuPath as a single MultiPolygon of every fragment on
    the slide, which as one feature would be one annotation covering the whole section.
    """
    features = []
    for cls, geoms in named_geoms.items():
        parts = []
        for g in geoms:
            if g.is_empty:
                continue
            parts.extend(list(g.geoms) if g.geom_type == "MultiPolygon" else [g])
        for i, g in enumerate(parts):
            features.append({
                "type": "Feature",
                "id": f"{cls}-{i + 1}",
                "geometry": {"type": "Polygon", "coordinates": _rings(g)},
                "properties": {"objectType": "annotation",
                               "name": f"{cls} {i + 1}",
                               "classification": {"name": cls,
                                                  "colorRGB": QUPATH_COLOR.get(cls, -1)}},
            })
    doc = {"type": "FeatureCollection", "features": features}
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc))
    return {"path": str(path), "n_features": len(features),
            "per_class": {c: sum(1 for f in features
                                 if f["properties"]["classification"]["name"] == c)
                          for c in named_geoms}}


def annotations_csv(named_geoms: dict[str, list], path: Path) -> dict:
    """The same polygons as ``Selection,X,Y``, for Explorer versions older than the GeoJSON path.

    One row per vertex, exteriors only -- the CSV form has no way to express a hole, and
    silently dropping one is better stated than discovered. Coordinates are Xenium microns.
    """
    import shapely

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = n_sel = n_holes = 0
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Selection", "X", "Y"])
        for cls, geoms in named_geoms.items():
            parts = []
            for g in geoms:
                if g.is_empty:
                    continue
                parts.extend(list(g.geoms) if g.geom_type == "MultiPolygon" else [g])
            for i, g in enumerate(parts):
                n_sel += 1
                n_holes += len(g.interiors)
                for x, y in shapely.get_coordinates(g.exterior):
                    w.writerow([f"{cls} {i + 1}", repr(float(x)), repr(float(y))])
                    n_rows += 1
    return {"path": str(path), "n_selections": n_sel, "n_vertices": n_rows,
            "n_holes_dropped": n_holes}


# --------------------------------------------------------------------------- #
# The cell groups
# --------------------------------------------------------------------------- #

def write_cell_groups(path: Path, cell_ids, groups, colors: dict[str, str]) -> dict:
    """Explorer's cell-group CSV: ``cell_id,group,color``, header in that exact order.

    ``color`` is optional from Explorer v4.1 and is written here because the alternative is a
    palette assigned by import order, which would give the same label a different colour in
    two donors. A cell whose group has no colour is refused rather than written blank.
    """
    import pandas as pd

    cell_ids = np.asarray(cell_ids)
    groups = np.asarray(groups)
    if len(cell_ids) != len(groups):
        raise ExplorerExportError(f"{len(cell_ids)} cells against {len(groups)} labels")
    missing = sorted(set(map(str, np.unique(groups))) - set(colors))
    if missing:
        raise ExplorerExportError(f"no colour recorded for group(s) {missing}")

    df = pd.DataFrame({"cell_id": cell_ids, "group": groups,
                       "color": [colors[str(g)] for g in groups]})
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    counts = df.group.value_counts().to_dict()
    return {"path": str(path), "n_cells": int(len(df)),
            "groups": {str(k): int(v) for k, v in sorted(counts.items())}}


def two_component_threshold(values: np.ndarray, *, max_iter: int = 200,
                            tol: float = 1e-8, seed: int = 0) -> dict:
    """Where one donor's own log intensity stops looking like its lower component.

    A two-component Gaussian mixture fitted by EM on ``log1p(value)``, with the threshold at
    the crossing where the posterior for the upper component reaches 0.5.

    **What makes this admissible under the repository's own rule.** The rule forbids tuning
    toward an expected marker frequency, cross-donor uniformity or disease ordering. Nothing
    here is told what fraction should come out positive: the mixing weight is *estimated from
    this donor's own distribution* and then reported, so the positive fraction is an output to
    be read, never an input to be hit. Every fitted parameter is written into the manifest for
    the same reason -- a threshold nobody can see the derivation of is a threshold nobody can
    disagree with.

    It is still a model, and it is only ever applied to a channel whose AUROC has already
    cleared the pre-registered bar. On a channel that has not, this number would be a
    confident split of noise.
    """
    v = np.asarray(values, dtype=np.float64)
    x = np.log1p(v[np.isfinite(v)])
    if x.size < 100:
        raise ExplorerExportError(f"only {x.size} finite values; refusing to fit a mixture")

    lo, hi = np.percentile(x, [25.0, 75.0])
    mu = np.array([lo, hi], dtype=np.float64)
    if not mu[1] > mu[0]:
        raise ExplorerExportError("the channel has no spread between its quartiles")
    sd = np.full(2, max(float(x.std()), 1e-6))
    w = np.array([0.5, 0.5])

    ll_prev = -np.inf
    for it in range(int(max_iter)):
        p = w * np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))
        tot = p.sum(axis=1)
        tot[tot <= 0] = np.finfo(float).tiny
        r = p / tot[:, None]
        nk = r.sum(axis=0)
        nk[nk <= 0] = np.finfo(float).tiny
        w = nk / len(x)
        mu = (r * x[:, None]).sum(axis=0) / nk
        sd = np.sqrt((r * (x[:, None] - mu) ** 2).sum(axis=0) / nk)
        sd = np.maximum(sd, 1e-6)
        ll = float(np.log(tot).sum())
        if abs(ll - ll_prev) < tol * max(1.0, abs(ll)):
            break
        ll_prev = ll

    order = np.argsort(mu)
    mu, sd, w = mu[order], sd[order], w[order]

    grid = np.linspace(mu[0], mu[1], 20_001)
    pu = w[1] * np.exp(-0.5 * ((grid - mu[1]) / sd[1]) ** 2) / sd[1]
    pl = w[0] * np.exp(-0.5 * ((grid - mu[0]) / sd[0]) ** 2) / sd[0]
    cross = np.nonzero(pu >= pl)[0]
    if cross.size == 0:
        raise ExplorerExportError(
            "the upper component never outweighs the lower one between their means; "
            "this channel is not two-component and no threshold is defined for it")
    thr_log = float(grid[cross[0]])

    return {"threshold": float(np.expm1(thr_log)), "threshold_log1p": thr_log,
            "means_log1p": mu.tolist(), "sds_log1p": sd.tolist(), "weights": w.tolist(),
            "n_fitted": int(x.size), "iterations": int(it + 1),
            "log_likelihood": float(ll),
            "separation_sds": float((mu[1] - mu[0]) / np.sqrt((sd ** 2).mean()))}


def decile_bands(values: np.ndarray) -> tuple[np.ndarray, dict]:
    """Per-donor deciles of a channel, as ordered labels -- a display, not a call.

    Cells whose value is unavailable (off the imported image) are labelled ``Unavailable``
    with the reason recorded beside them, rather than pooled into a bottom band where they
    would read as "measured, and negative".
    """
    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(v)
    labels = np.full(v.shape, "Unavailable", dtype=object)
    if ok.sum() < 10:
        return labels, {"n_available": int(ok.sum()), "edges": []}
    edges = np.percentile(v[ok], np.arange(10, 100, 10))
    idx = np.searchsorted(edges, v[ok], side="right")
    labels[ok] = [f"d{i + 1:02d}" for i in idx]
    return labels, {"n_available": int(ok.sum()),
                    "edges": [float(e) for e in edges]}
