"""Match internal tissue structure when a serial section images only part of its partner."""

from __future__ import annotations

import numpy as np

from .transform import affine_from_pixels


def _texture(image):
    import cv2

    smooth = cv2.GaussianBlur(image, (0, 0), 1.4)
    mean = cv2.GaussianBlur(smooth, (0, 0), 9)
    variance = cv2.GaussianBlur(smooth * smooth, (0, 0), 9) - mean * mean
    return np.clip((smooth - mean) / np.sqrt(np.maximum(variance, 0.0004)), -3, 3)


def _working_images(moving, fixed, mask, grid):
    import cv2

    # Retain fine regional structure without allocating another full-resolution slide.
    stride = max(1, 16 / grid.pixel_um, max(grid.shape) / 1600)
    h, w = grid.shape
    size = (max(1, round(w / stride)), max(1, round(h / stride)))
    sx, sy = size[0] / w, size[1] / h
    # OpenCV resize maps pixel centres with a half-pixel offset.
    resize = np.array([[sx, 0, (sx - 1) / 2], [0, sy, (sy - 1) / 2], [0, 0, 1]])

    def down(a):
        return cv2.resize(a, size, interpolation=cv2.INTER_AREA).astype(np.float32)

    support = down(np.asarray(mask, dtype=np.float32)) > 0.5
    return _texture(down(moving)), _texture(down(fixed)), support, resize


def _coarse_texture(moving, fixed, mask, step_deg, max_rotation_deg):
    import cv2
    from scipy.signal import fftconvolve

    # Trim the empty moving canvas. Rotate on a bounding square and correlate
    # linearly over all translations. A circular shift on the shared canvas can
    # wrap the cropped section onto the wrong lobe of its partner.
    ys, xs = np.nonzero(mask)
    left, top = xs.min(), ys.min()
    moving = (moving * mask)[top : ys.max() + 1, left : xs.max() + 1]
    source = cv2.pyrDown(moving)
    target = cv2.pyrDown(fixed)
    h, w = source.shape
    side = int(np.ceil(np.hypot(h, w))) + 2
    center = np.array([(w - 1) / 2, (h - 1) / 2])
    best_score = -np.inf
    best = None
    for mirrored in (False, True):
        reflected = np.diag([-1.0, 1.0] if mirrored else [1.0, 1.0])
        for angle in np.arange(0, 360, step_deg):
            rotation = cv2.getRotationMatrix2D((0.0, 0.0), float(angle), 1)[:, :2] @ reflected
            orientation = np.degrees(np.arctan2(rotation[1, 0], rotation[0, 0]))
            if abs(orientation) > max_rotation_deg + 1e-8:
                continue
            offset = (side - 1) / 2 - rotation @ center
            template = cv2.warpAffine(source, np.column_stack((rotation, offset)), (side, side))
            energy = np.sum(template**2)
            if energy <= 0:
                continue
            response = fftconvolve(target, template[::-1, ::-1], mode="full")
            y, x = np.unravel_index(response.argmax(), response.shape)
            score = response[y, x] / np.sqrt(energy)
            if score > best_score:
                best_score = score
                shift = 2 * (offset + [x - side + 1, y - side + 1]) - rotation @ [left, top]
                best = np.column_stack((rotation, shift))
    return best


def regional_texture_seed(moving, fixed, mask, grid, *, step_deg, max_rotation_deg):
    """Propose one interior-structure seed for the same masked ECC comparison."""
    if mask is None or min(grid.shape) < 32:
        return None
    moving, fixed, support, resize = _working_images(moving, fixed, mask, grid)
    if support.sum() < 64 or not moving[support].any() or not fixed.any():
        return None
    matrix = _coarse_texture(moving, fixed, support, step_deg, max_rotation_deg)
    if matrix is None:
        return None
    native = np.linalg.inv(resize) @ np.vstack((matrix, [0, 0, 1])) @ resize
    seed = affine_from_pixels(native[:2], grid)
    seed.meta["coarse_method"] = "regional_texture"
    return seed
