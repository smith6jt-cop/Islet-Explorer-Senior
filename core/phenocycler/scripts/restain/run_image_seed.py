#!/usr/bin/env python
"""Recover a fail-closed DAPI image seed when VALIS cannot produce one.

This stage is deliberately narrower than VALIS.  It fits one orientation-preserving
similarity from the two exported DAPI planes, refines that similarity with masked local ECC,
and asks a separately constructed, acquisition-field-aware nucleus detector whether the
placement contains more correspondence than a 25-um shifted null.  The transform convention
is always Xenium microns -> post-Xenium IF microns.

Passing writes ``image_seed.json`` beside ``export.json``.  Failing writes only
``image_seed_report.json``; it never creates a transform that downstream code could mistake
for releasable geometry, and it never reads or overwrites ``refined.json``::

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/run_image_seed.py --donors 112
    ../../.venv/bin/python scripts/restain/run_align.py --donors 112 \
        --transform image_seed --window-fold all
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import refine as RF  # noqa: E402
import registration_inputs as RI  # noqa: E402
from run_register import OUT  # noqa: E402

from phenocycler.integration.transform import Transform  # noqa: E402

# Coarse image fit.  Eight-micron pixels retain ducts and tissue compartments while keeping a
# whole slide small enough for deterministic feature extraction.
FEATURE_UM = 8.0
FEATURE_MODES = ("plain", "clahe", "highpass")
MAX_FEATURES = 50_000
FEATURE_RATIO = 0.78
FEATURE_RANSAC_UM = 12.0
MIN_FEATURE_INLIERS = 8
MIN_FEATURE_SPAN_FRACTION = 0.20

# Same-slide microscopes both report calibrated microns.  A large scale change is therefore a
# unit/provenance failure, not biology.  Reflection is impossible in this fallback by model
# construction and checked again before persistence.
MIN_TOTAL_SCALE = 0.95
MAX_TOTAL_SCALE = 1.05
MAX_SIMILARITY_ANISOTROPY = 1.001

# Fine image fit.  ECC is local so blank acquisition fields cannot pull the whole section.
ECC_PATCH_UM = 512.0
ECC_STEP_UM = 384.0
ECC_MASK_MARGIN_UM = 30.0
ECC_MIN_VALID_FRACTION = 0.25
ECC_MIN_CC = 0.30
ECC_MAX_SHIFT_UM = 40.0
ECC_RANSAC_UM = 3.0
ECC_MIN_PATCHES = 20
ECC_MIN_INLIERS = 12
ECC_MIN_INLIER_FRACTION = 0.25
ECC_MAX_RESIDUAL_MEDIAN_UM = 3.0
ECC_MIN_CORRECTION_SCALE = 0.98
ECC_MAX_CORRECTION_SCALE = 1.02
ECC_MAX_CORRECTION_ROTATION_DEG = 2.0

# Independent evidence.  The camera field dimensions come from qptiff XML; these constants
# describe nuclei, not donor-specific tiling.  z>=2 is intentionally permissive because the
# evidence is an enrichment against an identically sampled shifted null, not a claim that every
# peak is a nucleus.
DETECT_LOW_SIGMA_UM = 1.5
DETECT_HIGH_SIGMA_UM = 4.5
DETECT_MIN_SEPARATION_UM = 4.0
DETECT_MIN_Z = 2.0
DETECT_FIELD_MARGIN_UM = 14.0
MIN_XENIUM_NUCLEI_PER_FIELD = 20
EVIDENCE_RADIUS_UM = 2.0
EVIDENCE_SECOND_RADIUS_UM = 5.0
EVIDENCE_NULL_UM = 25.0

# Release requires both effect size and replication.  Five hundred excess correspondences are
# >100 observations per similarity degree of freedom; the ratio floors prevent a large dense
# field from turning a negligible effect into a nominally significant one.  The field-level
# sign test makes acquisition fields, rather than hundreds of thousands of correlated nuclei,
# the units of replication.
MIN_EVIDENCE_NUCLEI = 10_000
MIN_DETECTED_PEAKS = 5_000
# A same-slide qptiff acquisition must cover essentially the whole Xenium cell field.  These
# are physical support checks, not intensity-registration scores.  Existing independently
# registered cohort controls were all in bounds and placed 99.63--100% of nuclei on nonzero IF
# support.  The deliberately looser floors below admit edge cropping and dark pixels while
# preventing a locally attractive fragment from certifying a whole-slide transform.
MIN_FINITE_NUCLEUS_FRACTION = 1.0
MIN_IN_IF_BOUNDS_FRACTION = 0.98
MIN_IF_SUPPORT_FRACTION = 0.95
MIN_EXCESS_MATCHES_2UM = 500
MIN_EXCESS_MUTUAL_PAIRS = 500
MIN_ENRICHMENT_2UM = 1.20
MIN_ENRICHMENT_5UM = 1.10
MIN_BLOCKS = 20
MIN_POSITIVE_BLOCK_FRACTION = 0.60
MAX_BLOCK_SIGN_P = 0.01

OUTPUT_STEM = "image_seed"
REPORT_NAME = "image_seed_report.json"


class ImageSeedError(RuntimeError):
    """The image-only fallback cannot establish releasable geometry."""


def _atomic_json(path: Path, value: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False))
    tmp.replace(path)


def _persist_failure(donor_dir: Path, report: dict) -> dict:
    """Durably fail the stage before removing any stale candidate stem."""
    report_path = donor_dir / REPORT_NAME
    transform_path = donor_dir / f"{OUTPUT_STEM}.json"
    report["output"] = {
        "releasable": False,
        "image_seed_qc_pass": False,
        "release_status": "FAILED_IMAGE_SEED_QC",
        "transform_file": None,
        "required_gate_command": None,
    }
    _atomic_json(report_path, report)
    invalidated = bool(transform_path.exists())
    transform_path.unlink(missing_ok=True)
    (donor_dir / f"{OUTPUT_STEM}_displacement.npz").unlink(missing_ok=True)
    (donor_dir / f"{OUTPUT_STEM}_failure_report.json").unlink(missing_ok=True)
    report["output"]["invalidated_prior_transform"] = invalidated
    # The FAIL verdict was already durable before its transform stem vanished.  Re-write only
    # to record whether this attempt invalidated a stale candidate.
    _atomic_json(report_path, report)
    report["report_file"] = str(report_path)
    return report


def _display_u8(image: np.ndarray, mode: str) -> np.ndarray:
    import cv2

    a = np.asarray(image, dtype=np.float32)
    nz = a[a > 0]
    if nz.size < 100:
        raise ImageSeedError("DAPI plane has too few nonzero pixels for image registration")
    lo, hi = np.percentile(nz, [1.0, 99.8])
    if not hi > lo:
        raise ImageSeedError("DAPI plane has no usable intensity range")
    out = np.clip((a - lo) * (255.0 / (hi - lo)), 0, 255).astype(np.uint8)
    if mode == "clahe":
        out = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(32, 32)).apply(out)
    elif mode == "highpass":
        background = cv2.GaussianBlur(out, (0, 0), 5.0)
        out = cv2.normalize(cv2.subtract(out, background), None, 0, 255, cv2.NORM_MINMAX)
    elif mode != "plain":
        raise ValueError(f"unknown feature mode {mode!r}")
    return out


def _resize_for_features(image: np.ndarray, native_um: float, mode: str):
    import cv2

    scale = float(native_um) / FEATURE_UM
    if not np.isfinite(scale) or scale <= 0:
        raise ImageSeedError(f"invalid native pixel size {native_um}")
    shape = (max(1, int(round(image.shape[1] * scale))),
             max(1, int(round(image.shape[0] * scale))))
    small = cv2.resize(image, shape, interpolation=cv2.INTER_AREA)
    mask = cv2.resize((np.asarray(image) > 0).astype(np.uint8), shape,
                      interpolation=cv2.INTER_NEAREST)
    mask = cv2.erode(mask, np.ones((3, 3), np.uint8)) * 255
    return _display_u8(small, mode), mask, scale


def keypoints_um(keypoints, native_um: float, origin_um: float,
                 resize_scale: float) -> np.ndarray:
    """OpenCV resized-pixel keypoints -> exact exported-image micron centres."""
    p = np.asarray([point.pt for point in keypoints], dtype=np.float64)
    native_px = (p + 0.5) / float(resize_scale) - 0.5
    return float(origin_um) + float(native_um) * native_px


def _similarity_metrics(affine: np.ndarray) -> dict:
    a = np.asarray(affine, dtype=np.float64).reshape(2, 3)
    if not np.isfinite(a).all():
        return {"scale": None, "determinant": None, "anisotropy": None,
                "rotation_deg": None, "mirrored": None}
    singular = np.linalg.svd(a[:, :2], compute_uv=False)
    det = float(np.linalg.det(a[:, :2]))
    scale = float(np.sqrt(abs(det)))
    anisotropy = (float(singular.max() / singular.min())
                  if singular.min() > 0 else None)
    rotation = float(np.degrees(np.arctan2(a[1, 0], a[0, 0])))
    return {"scale": scale, "determinant": det, "anisotropy": anisotropy,
            "rotation_deg": rotation, "mirrored": bool(det < 0)}


def _validate_similarity(affine: np.ndarray, *, correction: bool = False) -> dict:
    metrics = _similarity_metrics(affine)
    lo, hi = ((ECC_MIN_CORRECTION_SCALE, ECC_MAX_CORRECTION_SCALE) if correction
              else (MIN_TOTAL_SCALE, MAX_TOTAL_SCALE))
    if metrics["scale"] is None:
        raise ImageSeedError("similarity contains non-finite coefficients")
    if metrics["mirrored"]:
        raise ImageSeedError("image fallback refuses a mirrored transform")
    if not lo <= metrics["scale"] <= hi:
        raise ImageSeedError(
            f"similarity scale {metrics['scale']:.6f} is outside [{lo}, {hi}]")
    if (metrics["anisotropy"] is None
            or metrics["anisotropy"] > MAX_SIMILARITY_ANISOTROPY):
        raise ImageSeedError(
            f"linear anisotropy {metrics['anisotropy']} is not a similarity")
    if correction and abs(metrics["rotation_deg"]) > ECC_MAX_CORRECTION_ROTATION_DEG:
        raise ImageSeedError(
            f"ECC correction rotation {metrics['rotation_deg']:.3f} exceeds "
            f"{ECC_MAX_CORRECTION_ROTATION_DEG} degrees")
    return metrics


def coarse_similarity(xenium: np.ndarray, if_image: np.ndarray, meta: dict):
    """SIFT/RANSAC similarity on an exact common physical resolution."""
    import cv2

    cv2.setRNGSeed(0)
    attempts, accepted = [], []
    for mode in FEATURE_MODES:
        x, xm, sx = _resize_for_features(
            xenium, float(meta["xenium_um_per_px"]), mode)
        y, ym, sy = _resize_for_features(if_image, float(meta["if_um_per_px"]), mode)
        sift = cv2.SIFT_create(nfeatures=MAX_FEATURES, contrastThreshold=0.005,
                               edgeThreshold=20, sigma=1.2)
        kx, dx = sift.detectAndCompute(x, xm)
        ky, dy = sift.detectAndCompute(y, ym)
        rec = {"mode": mode, "n_xenium_features": len(kx), "n_if_features": len(ky)}
        if dx is None or dy is None or len(dx) < 2 or len(dy) < 2:
            rec["reason"] = "too few descriptors"
            attempts.append(rec)
            continue
        knn = cv2.BFMatcher(cv2.NORM_L2).knnMatch(dx, dy, k=2)
        matches = [a for a, b in knn if a.distance < FEATURE_RATIO * b.distance]
        rec["n_ratio_matches"] = len(matches)
        if len(matches) < 3:
            rec["reason"] = "too few ratio-test matches"
            attempts.append(rec)
            continue
        xum = keypoints_um(kx, float(meta["xenium_um_per_px"]),
                           float(meta["xenium_origin_um"]), sx)
        yum = keypoints_um(ky, float(meta["if_um_per_px"]),
                           float(meta["if_origin_um"]), sy)
        src = np.ascontiguousarray([xum[m.queryIdx] for m in matches], dtype=np.float64)
        dst = np.ascontiguousarray([yum[m.trainIdx] for m in matches], dtype=np.float64)
        affine, mask = cv2.estimateAffinePartial2D(
            src, dst, method=cv2.RANSAC, ransacReprojThreshold=FEATURE_RANSAC_UM,
            maxIters=100_000, confidence=0.999, refineIters=50)
        if affine is None or mask is None:
            rec["reason"] = "RANSAC produced no similarity"
            attempts.append(rec)
            continue
        inside = mask.ravel().astype(bool)
        residual = np.linalg.norm(dst - (src @ affine[:, :2].T + affine[:, 2]), axis=1)
        full_span = np.array([
            xenium.shape[1] * float(meta["xenium_um_per_px"]),
            xenium.shape[0] * float(meta["xenium_um_per_px"]),
        ])
        span = np.ptp(src[inside], axis=0) / full_span if inside.any() else np.zeros(2)
        rec.update({
            "affine": affine.tolist(), "n_inliers": int(inside.sum()),
            "inlier_fraction": float(inside.mean()),
            "inlier_span_fraction_xy": span.tolist(),
            "inlier_residual_median_um": (float(np.median(residual[inside]))
                                           if inside.any() else None),
            "metrics": _similarity_metrics(affine),
        })
        reasons = []
        try:
            _validate_similarity(affine)
        except ImageSeedError as exc:
            reasons.append(str(exc))
        if inside.sum() < MIN_FEATURE_INLIERS:
            reasons.append(f"{inside.sum()} inliers < {MIN_FEATURE_INLIERS}")
        if np.any(span < MIN_FEATURE_SPAN_FRACTION):
            reasons.append(
                f"inlier span {span.tolist()} is below {MIN_FEATURE_SPAN_FRACTION}")
        rec["reason"] = "; ".join(reasons) if reasons else None
        attempts.append(rec)
        if not reasons:
            accepted.append((int(inside.sum()), float(np.prod(span)),
                             -float(np.median(residual[inside])), affine, rec))
    if not accepted:
        raise ImageSeedError("no coarse feature similarity cleared physical/coverage QC")
    _nin, _coverage, _residual, affine, winner = max(accepted, key=lambda item: item[:3])
    return np.asarray(affine, dtype=np.float64), {
        "common_um_per_px": FEATURE_UM, "winner_mode": winner["mode"],
        "attempts": attempts,
    }


def _normalised_float(image: np.ndarray) -> np.ndarray | None:
    a = np.asarray(image, dtype=np.float32)
    nz = a[a > 0]
    if nz.size < 1_000:
        return None
    lo, hi = np.percentile(nz, [1.0, 99.8])
    if not hi > lo:
        return None
    return np.clip((a - lo) / (hi - lo), 0, 1).astype(np.float32)


def _sample_xenium_patch(xenium: np.ndarray, affine: np.ndarray, meta: dict,
                          x0: int, y0: int, width: int, height: int) -> np.ndarray:
    """Sample Xenium DAPI on an IF-export pixel window through Xenium-um -> IF-um."""
    import cv2

    inverse = np.linalg.inv(np.vstack([affine, [0.0, 0.0, 1.0]]))[:2]
    cols, rows = np.meshgrid(np.arange(x0, x0 + width), np.arange(y0, y0 + height))
    dst_um = np.stack([
        float(meta["if_origin_um"]) + cols * float(meta["if_um_per_px"]),
        float(meta["if_origin_um"]) + rows * float(meta["if_um_per_px"]),
    ], axis=-1)
    src_um = dst_um @ inverse[:, :2].T + inverse[:, 2]
    map_x = ((src_um[..., 0] - float(meta["xenium_origin_um"]))
             / float(meta["xenium_um_per_px"])).astype(np.float32)
    map_y = ((src_um[..., 1] - float(meta["xenium_origin_um"]))
             / float(meta["xenium_um_per_px"])).astype(np.float32)
    return cv2.remap(xenium, map_x, map_y, cv2.INTER_LINEAR,
                     borderMode=cv2.BORDER_CONSTANT)


def _compose_affine(after: np.ndarray, before: np.ndarray) -> np.ndarray:
    return (np.vstack([after, [0.0, 0.0, 1.0]])
            @ np.vstack([before, [0.0, 0.0, 1.0]]))[:2]


def refine_with_ecc(coarse: np.ndarray, xenium: np.ndarray, if_image: np.ndarray,
                    xenium_nuclei: np.ndarray, meta: dict):
    """Fit one small similarity correction from replicated masked local translations."""
    import cv2

    if_um = float(meta["if_um_per_px"])
    patch_px = max(64, int(round(ECC_PATCH_UM / if_um)))
    step_px = max(32, int(round(ECC_STEP_UM / if_um)))
    margin_px = max(1, int(round(ECC_MASK_MARGIN_UM / if_um)))
    if patch_px >= min(if_image.shape):
        raise ImageSeedError("IF export is smaller than one ECC patch")
    warped = xenium_nuclei @ coarse[:, :2].T + coarse[:, 2]
    finite = np.isfinite(warped).all(axis=1)
    if not finite.any():
        raise ImageSeedError("Xenium nucleus coordinates contain no finite points")
    warped_px = np.column_stack([
        (warped[:, 0] - float(meta["if_origin_um"])) / if_um,
        (warped[:, 1] - float(meta["if_origin_um"])) / if_um,
    ])
    lo = np.floor(np.nanmin(warped_px[finite], axis=0)).astype(int)
    hi = np.ceil(np.nanmax(warped_px[finite], axis=0)).astype(int)
    x_start = max(patch_px // 2, lo[0] + patch_px // 2)
    x_stop = min(if_image.shape[1] - patch_px // 2, hi[0] - patch_px // 2)
    y_start = max(patch_px // 2, lo[1] + patch_px // 2)
    y_stop = min(if_image.shape[0] - patch_px // 2, hi[1] - patch_px // 2)
    if x_stop < x_start or y_stop < y_start:
        raise ImageSeedError("coarse similarity has no patch-sized IF overlap")
    centers_x = np.arange(x_start, x_stop + 1, step_px, dtype=int)
    centers_y = np.arange(y_start, y_stop + 1, step_px, dtype=int)
    eroder = np.ones((2 * margin_px + 1, 2 * margin_px + 1), np.uint8)
    gauss = max(1, int(round(15.0 / if_um)))
    gauss = gauss if gauss % 2 else gauss + 1
    rows = []
    for cy in centers_y:
        for cx in centers_x:
            x0, y0 = cx - patch_px // 2, cy - patch_px // 2
            src = _sample_xenium_patch(
                xenium, coarse, meta, x0, y0, patch_px, patch_px)
            dst = if_image[y0:y0 + patch_px, x0:x0 + patch_px]
            valid = cv2.erode(((src > 0) & (dst > 0)).astype(np.uint8), eroder)
            valid_fraction = float(valid.mean())
            if valid_fraction < ECC_MIN_VALID_FRACTION:
                continue
            src_f, dst_f = _normalised_float(src), _normalised_float(dst)
            if src_f is None or dst_f is None:
                continue
            warp = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float32)
            try:
                cc, warp = cv2.findTransformECC(
                    dst_f, src_f, warp, cv2.MOTION_TRANSLATION,
                    (cv2.TERM_CRITERIA_COUNT | cv2.TERM_CRITERIA_EPS, 300, 1e-6),
                    valid, gauss)
            except cv2.error:
                continue
            shift_um = warp[:, 2].astype(np.float64) * if_um
            if not np.isfinite(cc) or not np.isfinite(shift_um).all():
                continue
            if np.linalg.norm(shift_um) > ECC_MAX_SHIFT_UM:
                continue
            center_um = np.array([
                float(meta["if_origin_um"]) + cx * if_um,
                float(meta["if_origin_um"]) + cy * if_um,
            ])
            rows.append({"center_um": center_um.tolist(), "cc": float(cc),
                         "inverse_shift_um": shift_um.tolist(),
                         "valid_fraction": valid_fraction})
    reliable = [row for row in rows if row["cc"] >= ECC_MIN_CC]
    if len(reliable) < ECC_MIN_PATCHES:
        raise ImageSeedError(
            f"only {len(reliable)} ECC patches clear cc>={ECC_MIN_CC}; "
            f"need {ECC_MIN_PATCHES}")
    src = np.ascontiguousarray([row["center_um"] for row in reliable], dtype=np.float64)
    # ECC's warp is used with WARP_INVERSE_MAP: it maps an IF/template position back into the
    # prewarped Xenium input.  The forward correction therefore has the opposite sign.
    inv_shift = np.asarray([row["inverse_shift_um"] for row in reliable], dtype=np.float64)
    dst = np.ascontiguousarray(src - inv_shift)
    cv2.setRNGSeed(0)
    correction, mask = cv2.estimateAffinePartial2D(
        src, dst, method=cv2.RANSAC, ransacReprojThreshold=ECC_RANSAC_UM,
        maxIters=100_000, confidence=0.999, refineIters=50)
    if correction is None or mask is None:
        raise ImageSeedError("ECC patch translations produced no similarity correction")
    inside = mask.ravel().astype(bool)
    residual = np.linalg.norm(dst - (src @ correction[:, :2].T + correction[:, 2]), axis=1)
    if inside.sum() < ECC_MIN_INLIERS:
        raise ImageSeedError(f"only {inside.sum()} ECC RANSAC inliers; need {ECC_MIN_INLIERS}")
    if inside.mean() < ECC_MIN_INLIER_FRACTION:
        raise ImageSeedError(
            f"ECC inlier fraction {inside.mean():.3f} < {ECC_MIN_INLIER_FRACTION}")
    median = float(np.median(residual[inside]))
    if median > ECC_MAX_RESIDUAL_MEDIAN_UM:
        raise ImageSeedError(
            f"ECC inlier residual median {median:.3f} um > "
            f"{ECC_MAX_RESIDUAL_MEDIAN_UM}")
    correction_metrics = _validate_similarity(correction, correction=True)
    affine = _compose_affine(correction, coarse)
    total_metrics = _validate_similarity(affine)
    return affine, {
        "patch_um": ECC_PATCH_UM, "step_um": ECC_STEP_UM,
        "n_attempted_grid": int(len(centers_x) * len(centers_y)),
        "n_usable": len(rows), "n_reliable": len(reliable),
        "n_ransac_inliers": int(inside.sum()),
        "inlier_fraction": float(inside.mean()),
        "inlier_residual_median_um": median,
        "correction_affine": correction.tolist(),
        "correction_metrics": correction_metrics, "total_metrics": total_metrics,
        "patches": rows,
    }


def camera_field_shape(meta: dict) -> dict:
    """Acquisition-camera field size on the exported IF lattice, from qptiff XML."""
    import tifffile

    required = ("qptiff", "if_pixel_um", "if_dapi_page", "if_um_per_px",
                "if_channel_names", "if_pyramid_shapes")
    missing = [key for key in required if key not in meta]
    if missing:
        raise ImageSeedError(f"export.json is missing exact IF grid metadata {missing}")
    qptiff = Path(meta["qptiff"])
    if not qptiff.is_file():
        raise ImageSeedError(f"qptiff recorded in export.json is absent: {qptiff}")
    with tifffile.TiffFile(str(qptiff)) as tif:
        levels = list(getattr(tif.series[0], "levels", [tif.series[0]]))
        current_shapes = [tuple(level.shape[-2:]) for level in levels]
        pages = levels[0].pages
        dapi_page = int(meta["if_dapi_page"])
        if not 0 <= dapi_page < len(pages):
            raise ImageSeedError(f"if_dapi_page={dapi_page} is outside the qptiff")
        description = pages[dapi_page].description or ""
    recorded_shapes = [tuple(shape) for shape in meta["if_pyramid_shapes"]]
    if current_shapes != recorded_shapes:
        raise ImageSeedError(
            f"qptiff pyramid changed since export: {current_shapes} != {recorded_shapes}")
    names = re.findall(r"<Name>([^<]*)</Name>", description)
    recorded_names = list(meta["if_channel_names"])
    if (not names or names[0].strip().casefold() != "dapi"
            or not 0 <= dapi_page < len(recorded_names)
            or str(recorded_names[dapi_page]).strip().casefold() != "dapi"):
        raise ImageSeedError(
            f"export.json if_dapi_page={dapi_page} is not literal DAPI in current qptiff")
    try:
        root = ET.fromstring(description)
        camera = root.find(".//CameraSettings")
        if camera is None:
            raise ValueError("no CameraSettings element")
        width0 = int(camera.findtext(".//Width", ""))
        height0 = int(camera.findtext(".//Height", ""))
    except (ET.ParseError, TypeError, ValueError) as exc:
        # Some qptiff writers prepend non-XML text; keep a strict, scoped fallback rather than
        # accepting the first arbitrary Width in the description.
        block = re.search(r"<CameraSettings(?:\s[^>]*)?>.*?</CameraSettings>",
                          description, flags=re.DOTALL)
        widths = [] if block is None else re.findall(r"<Width>\s*(\d+)\s*</Width>", block.group())
        heights = [] if block is None else re.findall(
            r"<Height>\s*(\d+)\s*</Height>", block.group())
        if len(widths) != 1 or len(heights) != 1:
            raise ImageSeedError(
                f"cannot read one camera field shape from qptiff XML: {exc}") from exc
        width0, height0 = int(widths[0]), int(heights[0])
    pixel_values = {
        float(value) for value in re.findall(
            r"<PixelSizeMicrons>\s*([^<]+?)\s*</PixelSizeMicrons>", description)
    }
    base_um = float(meta["if_pixel_um"])
    # The XML tag is the preferred witness and is still checked whenever it exists. Some
    # exports carry none at all (this 70-plex Fusion qptiff does), in which case export.json
    # declares where its value came from and that declaration is what gets verified. Silently
    # accepting a mismatch is still refused.
    if pixel_values:
        if (len(pixel_values) != 1
                or not np.isclose(next(iter(pixel_values)), base_um, rtol=0, atol=1e-12)):
            raise ImageSeedError(
                f"qptiff DAPI-page PixelSizeMicrons {sorted(pixel_values)} does not equal "
                f"export.json if_pixel_um={base_um}")
    elif meta.get("pixel_size_source") == "qptiff XML PixelSizeMicrons":
        raise ImageSeedError(
            "export.json claims the pixel size came from qptiff XML PixelSizeMicrons, but the "
            "DAPI page carries no such tag")
    export_um = float(meta["if_um_per_px"])
    factor_float = export_um / base_um
    factor = int(round(factor_float))
    if factor <= 0 or not np.isclose(factor_float, factor, rtol=0, atol=1e-9):
        raise ImageSeedError(
            f"IF export factor {factor_float} from {export_um}/{base_um} is not an integer")
    if width0 % factor or height0 % factor:
        raise ImageSeedError(
            f"camera field {(width0, height0)} is not divisible by pyramid factor {factor}")
    return {"level0_width_px": width0, "level0_height_px": height0,
            "pyramid_factor": factor, "export_width_px": width0 // factor,
            "export_height_px": height0 // factor,
            "export_width_um": width0 // factor * export_um,
            "export_height_um": height0 // factor * export_um,
            "source": "qptiff DAPI-page CameraSettings XML + export.json pixel scales"}


def _parabolic(vm1: np.ndarray, v0: np.ndarray, vp1: np.ndarray) -> np.ndarray:
    denom = vm1 - 2.0 * v0 + vp1
    out = np.zeros_like(v0, dtype=np.float64)
    ok = np.isfinite(denom) & (denom != 0)
    out[ok] = 0.5 * (vm1[ok] - vp1[ok]) / denom[ok]
    out[np.abs(out) > 1.0] = 0.0
    return out


def tile_aware_peaks(if_image: np.ndarray, seed: Transform, xenium_nuclei: np.ndarray,
                     meta: dict, camera: dict):
    """Permissive DAPI maxima with local background and camera-field edge exclusion."""
    import cv2

    um = float(meta["if_um_per_px"])
    origin = float(meta["if_origin_um"])
    fw, fh = int(camera["export_width_px"]), int(camera["export_height_px"])
    margin = max(2, int(round(DETECT_FIELD_MARGIN_UM / um)))
    sep = max(1, int(round(DETECT_MIN_SEPARATION_UM / um)))
    low = max(0.8, DETECT_LOW_SIGMA_UM / um)
    high = max(2.0 * low, DETECT_HIGH_SIGMA_UM / um)
    warped = seed.apply(xenium_nuclei)
    px = np.column_stack([(warped[:, 0] - origin) / um, (warped[:, 1] - origin) / um])
    finite = np.isfinite(px).all(axis=1)
    keys = np.zeros((len(px), 2), dtype=np.int64)
    keys[finite] = np.floor(px[finite] / np.array([fw, fh])).astype(np.int64)
    unique, counts = np.unique(keys[finite], axis=0, return_counts=True)
    selected = {
        (int(key[0]), int(key[1])) for key, count in zip(unique, counts, strict=True)
        if count >= MIN_XENIUM_NUCLEI_PER_FIELD and key[0] >= 0 and key[1] >= 0
        and (key[0] + 1) * fw <= if_image.shape[1]
        and (key[1] + 1) * fh <= if_image.shape[0]
    }
    if not selected:
        raise ImageSeedError("coarse similarity overlaps no complete acquisition fields")
    eroder = np.ones((2 * margin + 1, 2 * margin + 1), np.uint8)
    maxima_kernel = np.ones((2 * sep + 1, 2 * sep + 1), np.uint8)
    points, responses, field_ids = [], [], []
    field_records = []
    for ix, iy in sorted(selected):
        x0, y0 = ix * fw, iy * fh
        image = np.asarray(if_image[y0:y0 + fh, x0:x0 + fw], dtype=np.float32)
        valid = cv2.erode((image > 0).astype(np.uint8), eroder).astype(bool)
        # OpenCV's erosion uses a foreground-valued border for this dtype.  Camera-field
        # boundaries are invalid regardless of pixel value, so exclude them explicitly.
        valid[:margin, :] = valid[-margin:, :] = False
        valid[:, :margin] = valid[:, -margin:] = False
        if float(valid.mean()) < 0.10:
            field_records.append({"ix": ix, "iy": iy, "status": "insufficient_valid"})
            continue
        dog = cv2.GaussianBlur(image, (0, 0), low) - cv2.GaussianBlur(image, (0, 0), high)
        sample = dog[valid]
        median = float(np.median(sample))
        sigma = float(1.4826 * np.median(np.abs(sample - median)))
        if sigma <= 0:
            sigma = float(sample.std())
        if not sigma > 0:
            field_records.append({"ix": ix, "iy": iy, "status": "no_spread"})
            continue
        maxima = (dog == cv2.dilate(dog, maxima_kernel)) & valid
        row, col = np.nonzero(maxima)
        z = (dog[row, col] - median) / sigma
        keep = z >= DETECT_MIN_Z
        row, col, z = row[keep], col[keep], z[keep]
        if len(row):
            v0 = dog[row, col].astype(np.float64)
            dr = _parabolic(dog[row - 1, col].astype(np.float64), v0,
                            dog[row + 1, col].astype(np.float64))
            dc = _parabolic(dog[row, col - 1].astype(np.float64), v0,
                            dog[row, col + 1].astype(np.float64))
            points.append(np.column_stack([
                origin + (x0 + col + dc) * um,
                origin + (y0 + row + dr) * um,
            ]))
            responses.append(z)
            field_ids.extend([(ix, iy)] * len(row))
        field_records.append({"ix": ix, "iy": iy, "status": "measured",
                              "valid_fraction": float(valid.mean()),
                              "n_peaks": int(len(row)), "median": median, "sigma": sigma})
    if not points:
        raise ImageSeedError("tile-aware detector found no DAPI peaks")
    return np.concatenate(points), np.concatenate(responses), np.asarray(field_ids), {
        "n_fields_selected": len(selected),
        "n_fields_measured": int(sum(row["status"] == "measured" for row in field_records)),
        "field_records": field_records,
        "parameters": {"low_sigma_um": DETECT_LOW_SIGMA_UM,
                       "high_sigma_um": DETECT_HIGH_SIGMA_UM,
                       "min_separation_um": DETECT_MIN_SEPARATION_UM,
                       "min_response_z": DETECT_MIN_Z,
                       "field_margin_um": DETECT_FIELD_MARGIN_UM},
    }


def enrichment_evidence(seed: Transform, xenium_nuclei: np.ndarray, peaks: np.ndarray,
                        if_image: np.ndarray, meta: dict, camera: dict) -> dict:
    """Target-vs-shifted-null nucleus agreement with acquisition fields as blocks."""
    from scipy.spatial import cKDTree
    from scipy.stats import binomtest

    um = float(meta["if_um_per_px"])
    origin = float(meta["if_origin_um"])
    fw, fh = int(camera["export_width_px"]), int(camera["export_height_px"])
    warped = seed.apply(xenium_nuclei)
    px = np.column_stack([(warped[:, 0] - origin) / um, (warped[:, 1] - origin) / um])
    finite = np.isfinite(px).all(axis=1)
    rounded = np.zeros((len(px), 2), dtype=np.int64)
    rounded[finite] = np.rint(px[finite]).astype(np.int64)
    in_bounds = (finite
                 & (rounded[:, 0] >= 0) & (rounded[:, 0] < if_image.shape[1])
                 & (rounded[:, 1] >= 0) & (rounded[:, 1] < if_image.shape[0]))
    on_support = np.zeros(len(px), dtype=bool)
    in_idx = np.flatnonzero(in_bounds)
    on_support[in_idx] = if_image[rounded[in_idx, 1], rounded[in_idx, 0]] > 0
    n_total = len(xenium_nuclei)
    finite_fraction = float(finite.sum() / n_total) if n_total else 0.0
    in_bounds_fraction = float(in_bounds.sum() / n_total) if n_total else 0.0
    support_fraction = float(on_support.sum() / n_total) if n_total else 0.0
    keys = np.zeros((len(px), 2), dtype=np.int64)
    keys[finite] = np.floor(px[finite] / np.array([fw, fh])).astype(np.int64)
    local = px - keys * np.array([fw, fh])
    margin_px = (EVIDENCE_NULL_UM + EVIDENCE_RADIUS_UM + DETECT_FIELD_MARGIN_UM) / um
    keep = (finite
            & (keys[:, 0] >= 0) & (keys[:, 1] >= 0)
            & ((keys[:, 0] + 1) * fw <= if_image.shape[1])
            & ((keys[:, 1] + 1) * fh <= if_image.shape[0])
            & (local[:, 0] >= margin_px)
            & (local[:, 0] < fw - margin_px)
            & (local[:, 1] >= margin_px)
            & (local[:, 1] < fh - margin_px))
    idx = np.flatnonzero(keep)
    if idx.size:
        col = np.rint(px[idx, 0]).astype(int)
        row = np.rint(px[idx, 1]).astype(int)
        null_col = np.rint(px[idx, 0] + EVIDENCE_NULL_UM / um).astype(int)
        valid = ((if_image[row, col] > 0) & (if_image[row, null_col] > 0))
        idx = idx[valid]
    moving = xenium_nuclei[idx]
    score = RF.score_transform(seed, moving, peaks, null_offset_um=EVIDENCE_NULL_UM)
    placed = seed.apply(moving)
    target_mutual = len(RF.mutual_pairs(placed, peaks, max_dist_um=RF.PAIR_MAX_UM)[2])
    null_mutual = len(RF.mutual_pairs(
        placed + np.array([EVIDENCE_NULL_UM, 0.0]), peaks,
        max_dist_um=RF.PAIR_MAX_UM)[2])
    tree = cKDTree(peaks)
    target_distance = tree.query(placed, k=1)[0]
    null_distance = tree.query(placed + np.array([EVIDENCE_NULL_UM, 0.0]), k=1)[0]
    target2 = target_distance <= EVIDENCE_RADIUS_UM
    null2 = null_distance <= EVIDENCE_RADIUS_UM
    target5 = target_distance <= EVIDENCE_SECOND_RADIUS_UM
    null5 = null_distance <= EVIDENCE_SECOND_RADIUS_UM
    blocks = []
    if len(idx):
        evidence_keys = keys[idx]
        for key in np.unique(evidence_keys, axis=0):
            block = np.all(evidence_keys == key, axis=1)
            if block.sum() < MIN_XENIUM_NUCLEI_PER_FIELD:
                continue
            blocks.append({"ix": int(key[0]), "iy": int(key[1]),
                           "n_xenium": int(block.sum()),
                           "target_2um": int(target2[block].sum()),
                           "null_2um": int(null2[block].sum()),
                           "delta_2um": int(target2[block].sum() - null2[block].sum())})
    deltas = np.asarray([block["delta_2um"] for block in blocks], dtype=int)
    non_tied = deltas[deltas != 0]
    positive = int((non_tied > 0).sum())
    sign_p = (float(binomtest(positive, len(non_tied), 0.5, alternative="greater").pvalue)
              if len(non_tied) else 1.0)
    n_target2, n_null2 = int(target2.sum()), int(null2.sum())
    n_target5, n_null5 = int(target5.sum()), int(null5.sum())
    # Keep reports strict-JSON even in a sparse synthetic field with zero null hits.  The raw
    # counts remain beside these conservative denominator-floored ratios.
    ratio2 = float(n_target2 / max(1, n_null2))
    ratio5 = float(n_target5 / max(1, n_null5))
    positive_fraction = float(positive / len(non_tied)) if len(non_tied) else 0.0
    evidence = {
        "n_xenium_evidence": int(len(moving)), "n_if_peaks": int(len(peaks)),
        "whole_slide_support": {
            "n_xenium_total": int(n_total),
            "n_finite": int(finite.sum()), "finite_fraction": finite_fraction,
            "n_in_if_bounds": int(in_bounds.sum()),
            "in_if_bounds_fraction": in_bounds_fraction,
            "n_on_nonzero_if_pixel": int(on_support.sum()),
            "on_nonzero_if_pixel_fraction": support_fraction,
            "meaning": "whole-slide physical acquisition support; not a DAPI match score",
        },
        "score": score,
        "target_mutual_pairs": int(target_mutual), "null_mutual_pairs": int(null_mutual),
        "excess_mutual_pairs": int(target_mutual - null_mutual),
        "target_matches_2um": n_target2, "null_matches_2um": n_null2,
        "excess_matches_2um": n_target2 - n_null2,
        "enrichment_2um": ratio2,
        "target_matches_5um": n_target5, "null_matches_5um": n_null5,
        "excess_matches_5um": n_target5 - n_null5,
        "enrichment_5um": ratio5,
        "blocks": blocks, "n_non_tied_blocks": int(len(non_tied)),
        "positive_blocks": positive, "positive_block_fraction": positive_fraction,
        "one_sided_sign_test_p": sign_p,
        "null": {"translation_um": [EVIDENCE_NULL_UM, 0.0],
                 "sampling": "same Xenium nuclei, same acquisition-field interiors"},
    }
    reasons = []
    checks = [
        (finite_fraction >= MIN_FINITE_NUCLEUS_FRACTION,
         f"finite-nucleus fraction {finite_fraction:.3f} < "
         f"{MIN_FINITE_NUCLEUS_FRACTION}"),
        (in_bounds_fraction >= MIN_IN_IF_BOUNDS_FRACTION,
         f"nuclei in IF bounds {in_bounds_fraction:.3f} < "
         f"{MIN_IN_IF_BOUNDS_FRACTION}"),
        (support_fraction >= MIN_IF_SUPPORT_FRACTION,
         f"nuclei on nonzero IF acquisition support {support_fraction:.3f} < "
         f"{MIN_IF_SUPPORT_FRACTION}"),
        (len(moving) >= MIN_EVIDENCE_NUCLEI,
         f"evidence nuclei {len(moving)} < {MIN_EVIDENCE_NUCLEI}"),
        (len(peaks) >= MIN_DETECTED_PEAKS,
         f"detected peaks {len(peaks)} < {MIN_DETECTED_PEAKS}"),
        (target_mutual - null_mutual >= MIN_EXCESS_MUTUAL_PAIRS,
         f"excess mutual pairs {target_mutual-null_mutual} < {MIN_EXCESS_MUTUAL_PAIRS}"),
        (n_target2 - n_null2 >= MIN_EXCESS_MATCHES_2UM,
         f"excess 2-um matches {n_target2-n_null2} < {MIN_EXCESS_MATCHES_2UM}"),
        (ratio2 >= MIN_ENRICHMENT_2UM,
         f"2-um enrichment {ratio2:.3f} < {MIN_ENRICHMENT_2UM}"),
        (ratio5 >= MIN_ENRICHMENT_5UM,
         f"5-um enrichment {ratio5:.3f} < {MIN_ENRICHMENT_5UM}"),
        (len(non_tied) >= MIN_BLOCKS,
         f"non-tied acquisition fields {len(non_tied)} < {MIN_BLOCKS}"),
        (positive_fraction >= MIN_POSITIVE_BLOCK_FRACTION,
         f"positive-field fraction {positive_fraction:.3f} < "
         f"{MIN_POSITIVE_BLOCK_FRACTION}"),
        (sign_p <= MAX_BLOCK_SIGN_P,
         f"field sign-test p {sign_p:.4g} > {MAX_BLOCK_SIGN_P}"),
    ]
    reasons.extend(reason for passed, reason in checks if not passed)
    evidence["qc"] = {"verdict": "PASS" if not reasons else "FAIL", "reasons": reasons}
    return evidence


def _required_metadata(meta: dict, xenium: np.ndarray, if_image: np.ndarray) -> dict:
    required = ("donor_id", "serial", "roi", "xenium_um_per_px", "if_um_per_px",
                "xenium_origin_um", "if_origin_um", "qptiff", "if_pixel_um",
                "if_dapi_page", "if_channel_names", "if_pyramid_shapes",
                "pixel_size_source")
    missing = [key for key in required if key not in meta]
    if missing:
        raise ImageSeedError(f"export.json is missing {missing}")
    for key in ("xenium_um_per_px", "if_um_per_px", "xenium_origin_um", "if_origin_um"):
        if not np.isfinite(float(meta[key])):
            raise ImageSeedError(f"export.json records non-finite {key}={meta[key]}")
    if float(meta["xenium_um_per_px"]) <= 0 or float(meta["if_um_per_px"]) <= 0:
        raise ImageSeedError("exported pixel sizes must be positive")
    trusted_sources = ("qptiff XML PixelSizeMicrons",
                       "registration inventory (fitted, cross-checked)")
    if meta["pixel_size_source"] not in trusted_sources:
        raise ImageSeedError(
            f"untrusted IF pixel-size source {meta['pixel_size_source']!r}")
    nested_required = {
        "xenium": ("shape", "um_per_px", "origin_um", "factor", "level", "pyramid",
                   "levels"),
        "if": ("shape", "um_per_px", "origin_um", "factor", "level", "pyramid", "qptiff"),
    }
    for plane, keys in nested_required.items():
        if not isinstance(meta.get(plane), dict):
            raise ImageSeedError(f"export.json is missing the exact {plane!r} grid record")
        absent = [key for key in keys if key not in meta[plane]]
        if absent:
            raise ImageSeedError(f"export.json {plane!r} grid is missing {absent}")
    recorded_x = tuple(meta["xenium"]["shape"])
    recorded_if = tuple(meta["if"]["shape"])
    if recorded_x != tuple(xenium.shape):
        raise ImageSeedError(f"Xenium DAPI shape changed: {xenium.shape} != {recorded_x}")
    if recorded_if != tuple(if_image.shape):
        raise ImageSeedError(f"IF DAPI shape changed: {if_image.shape} != {recorded_if}")
    consistency = (
        ("xenium_um_per_px", "xenium", "um_per_px"),
        ("xenium_origin_um", "xenium", "origin_um"),
        ("if_um_per_px", "if", "um_per_px"),
        ("if_origin_um", "if", "origin_um"),
    )
    for top, plane, nested in consistency:
        if not np.isclose(float(meta[top]), float(meta[plane][nested]), rtol=0, atol=1e-12):
            raise ImageSeedError(
                f"export.json contradicts itself: {top}={meta[top]} but "
                f"{plane}.{nested}={meta[plane][nested]}")
    if Path(meta["qptiff"]).resolve() != Path(meta["if"]["qptiff"]).resolve():
        raise ImageSeedError("export.json top-level and IF-grid qptiff paths disagree")
    origin_check = meta.get("xenium_origin_check", {})
    if not bool(origin_check.get("brackets_cell_table", False)):
        raise ImageSeedError("export.json lacks a passing Xenium image-origin extent check")
    if "base_um_per_px" not in origin_check:
        raise ImageSeedError("export.json lacks the Xenium base pixel size used for its grid")
    for plane, expected_pyramid in (("xenium", "decimated"), ("if", "box_mean")):
        level = int(meta[plane]["level"])
        factor = int(meta[plane]["factor"])
        if level < 0 or factor != 2 ** level:
            raise ImageSeedError(
                f"export.json {plane} level/factor is not dyadic: level={level}, factor={factor}")
        if meta[plane]["pyramid"] != expected_pyramid:
            raise ImageSeedError(
                f"export.json {plane} pyramid={meta[plane]['pyramid']!r}; expected "
                f"{expected_pyramid!r} for exact pixel centres")
    xenium_base = float(origin_check["base_um_per_px"])
    xenium_factor = int(meta["xenium"]["factor"])
    if not np.isclose(float(meta["xenium_um_per_px"]), xenium_base * xenium_factor,
                      rtol=0, atol=1e-12):
        raise ImageSeedError("export.json Xenium scale does not equal base scale times factor")
    if not np.isclose(float(meta["xenium_origin_um"]), xenium_base / 2.0,
                      rtol=0, atol=1e-12):
        raise ImageSeedError("export.json Xenium decimated-grid origin is not base_um/2")
    if_factor = int(meta["if"]["factor"])
    if not np.isclose(float(meta["if_um_per_px"]), float(meta["if_pixel_um"]) * if_factor,
                      rtol=0, atol=1e-12):
        raise ImageSeedError("export.json IF scale does not equal qptiff scale times factor")
    if not np.isclose(float(meta["if_origin_um"]), float(meta["if_pixel_um"]) * if_factor / 2.0,
                      rtol=0, atol=1e-12):
        raise ImageSeedError("export.json IF box-mean-grid origin is not half an output pixel")
    if_shapes = meta.get("if_pyramid_shapes")
    if (not isinstance(if_shapes, list)
            or tuple(if_shapes[int(meta["if"]["level"])]) != recorded_if):
        raise ImageSeedError("export.json IF exported shape disagrees with its pyramid record")
    xenium_levels = meta["xenium"]["levels"]
    if (not isinstance(xenium_levels, list)
            or tuple(xenium_levels[int(meta["xenium"]["level"])]) != recorded_x):
        raise ImageSeedError("export.json Xenium exported shape disagrees with its pyramid record")
    return {"xenium": {"shape": list(xenium.shape),
                       "um_per_px": float(meta["xenium_um_per_px"]),
                       "origin_um": float(meta["xenium_origin_um"])},
            "if": {"shape": list(if_image.shape),
                   "um_per_px": float(meta["if_um_per_px"]),
                   "origin_um": float(meta["if_origin_um"]),
                   "base_um_per_px": float(meta["if_pixel_um"]),
                   "dapi_page": int(meta["if_dapi_page"]),
                   "qptiff": str(meta["qptiff"])}}


def fit_image_seed(donor_dir: Path) -> tuple[Transform, dict]:
    """Compute a candidate and its independent release evidence; do not write it."""
    import tifffile

    donor_dir = Path(donor_dir)
    need = [donor_dir / "export.json",
            donor_dir / "slides" / "01_xenium_dapi.ome.tiff",
            donor_dir / "slides" / "02_if_dapi.ome.tiff",
            donor_dir / "xenium_nuclei_um.npy"]
    missing = [str(path) for path in need if not path.is_file()]
    if missing:
        raise ImageSeedError(f"missing image-seed inputs: {missing}")
    meta = json.loads((donor_dir / "export.json").read_text())
    xenium = tifffile.imread(str(need[1]))
    if_image = tifffile.imread(str(need[2]))
    xenium_nuclei = np.load(need[3])
    physical_grid = _required_metadata(meta, xenium, if_image)
    camera = camera_field_shape(meta)
    coarse, coarse_report = coarse_similarity(xenium, if_image, meta)
    affine, ecc_report = refine_with_ecc(
        coarse, xenium, if_image, xenium_nuclei, meta)
    transform = Transform(affine=affine, meta={
        "source": "run_image_seed.py constrained common-grid SIFT + masked local ECC",
        "direction": "Xenium microns -> IF microns",
        "release_requires": "image_seed_report.json qc PASS and run_align.py image gate PASS",
    })
    peaks, responses, _field_ids, detector = tile_aware_peaks(
        if_image, transform, xenium_nuclei, meta, camera)
    detector.update({"n_peaks": int(len(peaks)),
                     "response_z_quantiles": {
                         str(q): float(np.percentile(responses, q))
                         for q in (0, 25, 50, 75, 90, 99, 100)}})
    evidence = enrichment_evidence(
        transform, xenium_nuclei, peaks, if_image, meta, camera)
    report = {
        "stage": OUTPUT_STEM, "donor_id": str(meta["donor_id"]),
        "serial": str(meta["serial"]), "roi": str(meta["roi"]),
        "direction": "Xenium microns -> IF microns",
        "physical_grid": physical_grid, "camera_field": camera,
        "coarse": coarse_report, "ecc": ecc_report,
        "candidate": transform.to_dict(), "detector": detector,
        "evidence": evidence,
        "qc": evidence["qc"],
        "constants": {
            "MIN_TOTAL_SCALE": MIN_TOTAL_SCALE, "MAX_TOTAL_SCALE": MAX_TOTAL_SCALE,
            "MIN_FEATURE_INLIERS": MIN_FEATURE_INLIERS,
            "MIN_FEATURE_SPAN_FRACTION": MIN_FEATURE_SPAN_FRACTION,
            "ECC_MIN_CC": ECC_MIN_CC, "ECC_MIN_PATCHES": ECC_MIN_PATCHES,
            "ECC_MIN_INLIERS": ECC_MIN_INLIERS,
            "EVIDENCE_MUTUAL_RADIUS_UM": RF.PAIR_MAX_UM,
            "EVIDENCE_NULL_UM": EVIDENCE_NULL_UM,
            "MIN_EVIDENCE_NUCLEI": MIN_EVIDENCE_NUCLEI,
            "MIN_DETECTED_PEAKS": MIN_DETECTED_PEAKS,
            "MIN_FINITE_NUCLEUS_FRACTION": MIN_FINITE_NUCLEUS_FRACTION,
            "MIN_IN_IF_BOUNDS_FRACTION": MIN_IN_IF_BOUNDS_FRACTION,
            "MIN_IF_SUPPORT_FRACTION": MIN_IF_SUPPORT_FRACTION,
            "MIN_EXCESS_MATCHES_2UM": MIN_EXCESS_MATCHES_2UM,
            "MIN_EXCESS_MUTUAL_PAIRS": MIN_EXCESS_MUTUAL_PAIRS,
            "MIN_ENRICHMENT_2UM": MIN_ENRICHMENT_2UM,
            "MIN_ENRICHMENT_5UM": MIN_ENRICHMENT_5UM,
            "MIN_BLOCKS": MIN_BLOCKS,
            "MIN_POSITIVE_BLOCK_FRACTION": MIN_POSITIVE_BLOCK_FRACTION,
            "MAX_BLOCK_SIGN_P": MAX_BLOCK_SIGN_P,
        },
        "qc_contract": {
            "unit_of_replication": "qptiff acquisition camera field",
            "null": "same nuclei and valid field interiors translated +25 um in IF x",
            "whole_slide_support_basis": (
                "same physical section must be covered by both acquisitions; independently "
                "registered cohort controls place 99.63--100% of nuclei on nonzero IF support"
            ),
            "basis": (
                "predeclared operational effect-size floors plus a one-sided exact sign test; "
                "500 excess mutual pairs is the minimum seed support supplied to refinement"
            ),
        },
    }
    return transform, report


def image_seed_one(donor_dir: Path, *, overwrite: bool = False) -> dict:
    """Fit, report, and persist only a candidate that clears the enrichment contract."""
    donor_dir = Path(donor_dir)
    transform_path = donor_dir / f"{OUTPUT_STEM}.json"
    report_path = donor_dir / REPORT_NAME
    if transform_path.exists() and not overwrite:
        raise ImageSeedError(f"{transform_path} exists; pass --overwrite to replace it")
    try:
        transform, report = fit_image_seed(donor_dir)
    except Exception as exc:  # noqa: BLE001 - turn every computation failure into durable QC
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        report = {
            "stage": OUTPUT_STEM, "donor_dir": str(donor_dir),
            "direction": "Xenium microns -> IF microns",
            "qc": {"verdict": "FAIL",
                   "reasons": [f"{type(exc).__name__}: {exc}"]},
        }
        return _persist_failure(donor_dir, report)

    passed = report["qc"]["verdict"] == "PASS"
    report["output"] = {
        # image_seed is always a candidate.  Registration release requires a separate PASS
        # from run_align, whose report is deliberately not forged or updated by this stage.
        "releasable": False,
        "image_seed_qc_pass": bool(passed),
        "release_status": "PENDING_ALIGNMENT_GATE" if passed else "FAILED_IMAGE_SEED_QC",
        "transform_file": str(transform_path) if passed else None,
        "required_gate_command": (
            f"run_align.py --donors {report['donor_id']} --transform {OUTPUT_STEM} "
            "--window-fold all") if passed else None,
    }
    if passed:
        # Transform.save writes one JSON for an affine-only transform.  Keep it temporary until
        # the PASS report is durable; readers can therefore never observe a candidate stem
        # whose evidence report was not successfully serialized.
        tmp_stem = f".{OUTPUT_STEM}.tmp-{os.getpid()}"
        tmp = transform.save(donor_dir, stem=tmp_stem)
        report["output"]["transform_sha256"] = hashlib.sha256(tmp.read_bytes()).hexdigest()
        _atomic_json(report_path, report)
        tmp.replace(transform_path)
        (donor_dir / f"{OUTPUT_STEM}_failure_report.json").unlink(missing_ok=True)
    else:
        return _persist_failure(donor_dir, report)
    report["report_file"] = str(report_path)
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--overwrite", action="store_true")
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
    selected = RI.select_registration_inputs(records, args.donors, context="run_image_seed")
    for record in selected:
        RI.require_registration_files(record)
    rows, failed = [], False
    for record in selected:
        row = image_seed_one(
            args.out / "valis" / f"donor_id={record.donor_id}", overwrite=args.overwrite)
        rows.append(row)
        qc = row["qc"]
        failed |= qc["verdict"] != "PASS"
        print(f"[{record.donor_id}] image seed {qc['verdict']}: "
              f"{'; '.join(qc.get('reasons', [])) or 'cleared physical and enrichment QC'}",
              flush=True)
        if row.get("output", {}).get("required_gate_command"):
            print(f"    next: {row['output']['required_gate_command']}", flush=True)
    if not rows:
        print("no donor matched", file=sys.stderr)
        return 1
    summary = args.out / f"image_seed_{args.roi}.json"
    _atomic_json(summary, rows)
    print(f"-> {summary}", flush=True)
    return 4 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
