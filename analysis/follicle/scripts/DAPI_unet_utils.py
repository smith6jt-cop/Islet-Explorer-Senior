"""Shared utilities for the DAPI / Phenocycler follicle U-Net.

Companion to `HE_unet_utils.py`. The two pipelines solve the same problem —
segment splenic follicles from a single-channel *nuclear density* image — and
deliberately share the same on-disk data contract, the same tiled-inference
core, and the same donor-level CV, so their training sets can be pooled.

What is different here:

  * **Input.** The nuclear channel is DAPI read straight out of a multiplex
    OME-TIFF, not a stain deconvolution. `dapi_to_uint8` normalises and
    **inverts** it so high nuclear density is DARK — matching the polarity of
    `HE_unet_utils.rgb_to_h`. Without that inversion the two cohorts could not
    share weights.
  * **Labels.** They come from the `SpleenRegions*` QuPath pixel classifier,
    which used 25-33 marker channels. This model sees only DAPI, so it is
    distilling a multi-channel teacher into a single-channel student. Some
    error is irreducible — most of all Follicle vs PALS, which the teacher
    separated with CD20 vs CD3e and which look near-identical in DAPI alone.
  * **Head.** One decoder, `classes = 1 + K`: channel 0 is the binary Follicle
    logit (the deliverable), channels 1.. are an auxiliary region softmax
    (RedPulp / Follicle / PALS / Trabeculae / LargeVessel) used only as extra
    training signal.

PARITY REQUIREMENT: `otsu_from_hist256`, `hist256` and `dapi_to_uint8` below MUST
stay numerically identical to the corresponding blocks in
`export_DAPI_training_pairs.groovy`. Every step is integer / rational arithmetic
over a 256-bin histogram precisely so the two languages agree to the byte. Any
change here has to be made in both files in the same commit, and
`scripts/check_dapi_parity.py` re-verifies it against real exports.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from PIL import Image

import HE_unet_utils as U

PROJECT = U.PROJECT
DATA_DIR = U.DATA_ROOT / "data_dapi"
MODELS_DIR = U.OUTPUT_ROOT / "analysis" / "DAPI_unet_models"
PRED_DIR = U.OUTPUT_ROOT / "analysis" / "DAPI_unet_predictions"

TARGET_UM_PER_PX = U.TARGET_UM_PER_PX  # 2.0 — same grid as the H&E pipeline
PATCH = U.PATCH
POS_FRAC = U.POS_FRAC

# Normalisation percentiles, applied within the tissue mask.
LO_PCT = 1.0
HI_PCT = 99.5

# ---------------------------------------------------------------------------
# Region label space
# ---------------------------------------------------------------------------
# Index 0 is "background or unlabeled". QuPath's `Ignore*` class produces no
# annotation (Workflow.groovy does not pass INCLUDE_IGNORED), so an unpainted
# pixel means "the classifier declined to label this", not "confidently
# background" — worth remembering when reading the aux head's numbers.
#
# The exporter paints in THIS order, later winning, so Follicle is painted last
# and `regions == FOLLICLE_INDEX` agrees pixel-for-pixel with the binary `masks`
# PNG (which the exporter derives from the label map for exactly that reason).
REGION_CLASSES = ["Background", "RedPulp", "Trabeculae", "LargeVessel", "PALS", "Follicle"]
FOLLICLE_INDEX = REGION_CLASSES.index("Follicle")
PALS_INDEX = REGION_CLASSES.index("PALS")

# SmallVessel is deliberately absent: its labels were made at native
# 0.508 µm/px with a 10 µm² floor (~3.5 µm across), so at the 2 µm/px training
# scale it is sub-pixel. `--region-classes` can add it back for anyone who wants
# to train a dedicated vessel model at native resolution.

# Collapsing Follicle and PALS into one white-pulp class. Useful as an honest
# ablation: if the model scores far better this way, the residual error really
# was the DAPI-invisible Follicle/PALS boundary rather than segmentation quality.
MERGED_REGION_CLASSES = ["Background", "RedPulp", "Trabeculae", "LargeVessel", "WhitePulp"]


def merge_white_pulp(region: np.ndarray) -> np.ndarray:
    """Map a REGION_CLASSES label map onto MERGED_REGION_CLASSES."""
    out = np.asarray(region).copy()
    out[out == FOLLICLE_INDEX] = PALS_INDEX          # Follicle -> WhitePulp slot
    return out


# ---------------------------------------------------------------------------
# DAPI normalisation — the Groovy/Python parity contract
# ---------------------------------------------------------------------------
def hist256(raw: np.ndarray) -> np.ndarray:
    """256-bin histogram of a uint16 image, binned by `value >> 8`.

    Binning by a shift (not by min/max) keeps the bin edges slide-independent,
    so `lo`/`hi`/Otsu are functions of the pixel values alone and reproduce
    exactly in Groovy.
    """
    raw = np.asarray(raw)
    if raw.dtype != np.uint16:
        raw = np.clip(raw, 0, 65535).astype(np.uint16)
    return np.bincount((raw >> 8).ravel(), minlength=256)[:256].astype(np.int64)


def otsu_from_hist256(hist: np.ndarray) -> int:
    """Otsu threshold as a BIN index in [0, 255] (integer arithmetic only).

    Returns the bin `t` maximising between-class variance for the split
    `bins <= t` (background) vs `bins > t` (foreground). Note the asymmetry:
    bin `t` itself belongs to the BACKGROUND side, so tissue is `bin > t`, which
    is why `tissue_from_raw` cuts at `(t + 1) << 8` and not `t << 8`.

    Degenerate histograms (all mass in one bin) return 0, which makes almost the
    whole image tissue — the safe direction, since a too-generous tissue mask
    only widens negative-patch sampling.
    """
    hist = np.asarray(hist, dtype=np.float64)
    total = hist.sum()
    if total <= 0:
        return 0
    idx = np.arange(256, dtype=np.float64)
    sum_all = float((hist * idx).sum())
    w_b = np.cumsum(hist)                    # weight of bins <= t
    sum_b = np.cumsum(hist * idx)
    w_f = total - w_b
    valid = (w_b > 0) & (w_f > 0)
    if not valid.any():
        return 0
    mean_b = np.divide(sum_b, w_b, out=np.zeros_like(sum_b), where=w_b > 0)
    mean_f = np.divide(sum_all - sum_b, w_f, out=np.zeros_like(sum_b), where=w_f > 0)
    between = w_b * w_f * (mean_b - mean_f) ** 2
    between[~valid] = -1.0
    return int(np.argmax(between))


def otsu3_from_hist256(hist: np.ndarray) -> tuple[int, int]:
    """Three-class Otsu over the 256-bin histogram: returns `(t1, t2)`.

    Classes are bins `[0..t1]`, `[t1+1..t2]`, `[t2+1..255]`, chosen to maximise
    `sum_c(S_c^2 / W_c)` (the standard between-class-variance form). First
    maximum wins, so the Groovy transcription only has to use a strict `>`.

    Two classes are not enough here. A spleen DAPI histogram has THREE modes —
    glass, red pulp, white pulp — and 2-class Otsu on it puts the threshold
    between red pulp and white pulp, not between glass and tissue. That yields a
    "tissue" mask containing only white pulp, which would then supply every
    negative training patch and every percentile. `t1` is the glass/tissue split
    we actually want; `t2` (roughly red pulp vs white pulp) is returned for
    provenance only.

    Falls back to the 2-class threshold when the histogram cannot support three
    non-empty classes.
    """
    hist = np.asarray(hist, dtype=np.float64)
    total = float(hist.sum())
    if total <= 0:
        return 0, 0
    idx = np.arange(256, dtype=np.float64)
    # W[i] / S[i] cover bins strictly below i, so class (a..b] is W[b+1]-W[a+1].
    W = np.concatenate([[0.0], np.cumsum(hist)])
    S = np.concatenate([[0.0], np.cumsum(hist * idx)])
    best = -1.0
    bt = (0, 0)
    for t1 in range(255):
        w0, s0 = W[t1 + 1], S[t1 + 1]
        if w0 <= 0.0:
            continue
        for t2 in range(t1 + 1, 255):
            w1 = W[t2 + 1] - W[t1 + 1]
            if w1 <= 0.0:
                continue
            w2 = total - W[t2 + 1]
            if w2 <= 0.0:
                continue
            s1 = S[t2 + 1] - S[t1 + 1]
            s2 = S[256] - S[t2 + 1]
            v = s0 * s0 / w0 + s1 * s1 / w1 + s2 * s2 / w2
            if v > best:
                best = v
                bt = (t1, t2)
    if best < 0.0:
        t = otsu_from_hist256(hist)
        return t, t
    return bt


def percentile_from_hist256(hist: np.ndarray, pct: float) -> int:
    """Percentile as a raw uint16 value, read off the 256-bin histogram.

    Returns the *low edge* of the containing bin (`bin << 8`), so the answer is
    an exact integer in both languages. The 256-bin quantisation costs at most
    255 raw counts of precision, which is immaterial against a 16-bit range and
    is what buys byte-exact agreement with the Groovy export.
    """
    hist = np.asarray(hist, dtype=np.int64)
    total = int(hist.sum())
    if total <= 0:
        return 0
    target = pct / 100.0 * total
    cum = np.cumsum(hist)
    b = int(np.searchsorted(cum, target, side="left"))
    return min(b, 255) << 8


def dapi_transform_params(raw: np.ndarray) -> dict:
    """Derive `{otsu_bin, lo, hi}` from a raw uint16 DAPI image.

    Order matters and is mirrored in the exporter: histogram -> Otsu -> tissue
    mask -> percentiles *within tissue* -> normalise. Computing the percentiles
    over the whole image instead would let a large empty background dominate the
    low end and crush the tissue into a few grey levels.
    """
    raw = np.asarray(raw)
    full = hist256(raw)
    t1, t2 = otsu3_from_hist256(full)
    otsu_bin = int(t1)
    tissue_hist = full.copy()
    tissue_hist[:otsu_bin + 1] = 0      # keep bins > otsu_bin, i.e. tissue
    if tissue_hist.sum() <= 0:          # everything below threshold — use all
        tissue_hist = full
    lo = percentile_from_hist256(tissue_hist, LO_PCT)
    hi = percentile_from_hist256(tissue_hist, HI_PCT)
    if hi <= lo:                        # flat slide; avoid a zero-width window
        hi = lo + 256
    tissue_frac = float(tissue_hist.sum() / max(full.sum(), 1))
    return {"otsu_bin": otsu_bin, "otsu_upper_bin": int(t2),
            "lo": int(lo), "hi": int(hi),
            "lo_pct": LO_PCT, "hi_pct": HI_PCT,
            # Recorded so an implausible mask is visible in the sidecar rather
            # than only showing up as poor training later. A spleen section
            # normally covers a good fraction of the slide; a few percent means
            # the glass/tissue split went wrong.
            "tissue_frac": round(tissue_frac, 6),
            "invert": True, "polarity": "high_DAPI_dark"}


def dapi_to_uint8(raw: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """Normalise raw uint16 DAPI to uint8 and INVERT (high DAPI -> dark).

    `inten = floor(255 * (1 - clip((v - lo) / (hi - lo), 0, 1)))`

    The inversion is not cosmetic: `HE_unet_utils.rgb_to_h` writes high
    hematoxylin as dark, and a shared-weights model needs both cohorts to agree
    on which way "more nuclei" points.
    """
    raw = np.asarray(raw, dtype=np.float64)
    if hi <= lo:
        raise ValueError(f"dapi_to_uint8 needs hi > lo, got lo={lo} hi={hi}")
    x = np.clip((raw - lo) / (hi - lo), 0.0, 1.0)
    return np.floor(255.0 * (1.0 - x)).astype(np.uint8)


def tissue_from_raw(raw: np.ndarray, otsu_bin: int) -> np.ndarray:
    """Boolean tissue mask: raw DAPI strictly above the Otsu threshold bin.

    The `+ 1` is load-bearing. Otsu's `t` is the last BACKGROUND bin, so cutting
    at `t << 8` leaves the entire background bin inside the mask — on a slide
    that is mostly empty, that is enough to drag the 1st percentile down into the
    background and crush the tissue into a handful of grey levels.
    """
    return np.asarray(raw) >= ((int(otsu_bin) + 1) << 8)


# ---------------------------------------------------------------------------
# Multiplex OME-TIFF reading
# ---------------------------------------------------------------------------
def find_dapi_index(path: Path) -> int:
    """Index of the DAPI channel in an OME-TIFF, by name, falling back to 0.

    Channel names have been lost on some of these images before — that is why
    `scripts/RestoreChannelNames.groovy` exists — so trusting the name alone is
    not safe, and trusting index 0 alone is not either. Resolve by name, warn
    and fall back when the name is missing.
    """
    import tifffile
    try:
        with tifffile.TiffFile(path) as tif:
            xml = tif.ome_metadata or ""
    except Exception:  # noqa: BLE001
        xml = ""
    names = re.findall(r'<Channel\b[^>]*\bName="([^"]*)"', xml)
    for i, n in enumerate(names):
        if n.strip().upper().startswith("DAPI"):
            if i != 0:
                print(f"  note: {Path(path).name}: DAPI is channel {i}, not 0")
            return i
    print(f"  WARNING: {Path(path).name}: no channel named DAPI in OME-XML "
          f"({len(names)} names found); falling back to channel 0")
    return 0


def _ome_pixel_size_um(path: Path) -> float:
    """PhysicalSizeX from OME-XML, or NaN. Same parse as predict_HE_unet."""
    import tifffile
    try:
        with tifffile.TiffFile(path) as tif:
            desc = tif.ome_metadata or ""
        m = re.search(r'PhysicalSizeX="([0-9eE.+\-]+)"', desc)
        if m:
            return float(m.group(1))
    except Exception:  # noqa: BLE001
        pass
    return float("nan")


def read_dapi_channel(path: Path, target_um_per_px: float = TARGET_UM_PER_PX
                      ) -> tuple[np.ndarray, float]:
    """Read the DAPI plane of a multiplex OME-TIFF at ~`target_um_per_px`.

    Returns `(raw_uint16, effective_px_size_um)` — raw, NOT normalised; callers
    pair this with `dapi_transform_params` + `dapi_to_uint8`.

    OpenSlide (which `predict_HE_unet.py` uses) cannot select one channel of a
    33-channel uint16 OME-TIFF, so this goes through tifffile and the pyramid
    levels directly. It picks the finest level that is still at least as coarse
    as the target where one exists, then resamples to exactly the target —
    downsampling from a finer level rather than upsampling from a coarser one.
    """
    import tifffile

    path = Path(path)
    px = _ome_pixel_size_um(path)
    if not np.isfinite(px) or px <= 0:
        raise RuntimeError(f"{path.name}: missing pixel calibration (PhysicalSizeX).")
    ci = find_dapi_index(path)

    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        levels = list(getattr(series, "levels", [series])) or [series]
        base_w = _axis_len(levels[0], "X")
        # Pick the finest level whose pixel size does not exceed the target.
        best, best_px = levels[0], px
        for lvl in levels:
            w = _axis_len(lvl, "X")
            if w <= 0:
                continue
            lvl_px = px * (base_w / w)
            if lvl_px <= target_um_per_px * 1.001 and lvl_px >= best_px:
                best, best_px = lvl, lvl_px
        arr = _read_channel_plane(best, ci)

    raw = np.asarray(arr)
    if raw.ndim != 2:
        raise RuntimeError(f"{path.name}: expected a 2-D DAPI plane, got {raw.shape}")
    if raw.dtype != np.uint16:
        raw = np.clip(raw, 0, 65535).astype(np.uint16)

    scale = best_px / target_um_per_px
    if abs(scale - 1.0) > 1e-3:
        new_w = max(1, int(round(raw.shape[1] * scale)))
        new_h = max(1, int(round(raw.shape[0] * scale)))
        # BILINEAR on 'I;16' is unsupported in Pillow — resample as float32 then
        # round back, which also avoids the uint16 wraparound of a naive cast.
        img = Image.fromarray(raw.astype(np.float32), mode="F")
        raw = np.asarray(img.resize((new_w, new_h), Image.BILINEAR))
        raw = np.clip(np.rint(raw), 0, 65535).astype(np.uint16)
    return raw, float(target_um_per_px)


def _axis_len(level, axis: str) -> int:
    axes = getattr(level, "axes", "")
    shape = getattr(level, "shape", ())
    i = axes.find(axis)
    return int(shape[i]) if 0 <= i < len(shape) else 0


def _read_channel_plane(level, ci: int) -> np.ndarray:
    """Read plane `ci` of a tifffile series level along its channel axis.

    Reads the single IFD page for that channel where possible. `level.asarray()`
    would decode EVERY channel first: these files are CYX with one page per
    channel, so on a 33-channel slide at 2 µm/px that is ~7 GB decoded to reach
    one ~200 MB plane. The whole reason the Groovy export uses
    `extractChannels` is to avoid exactly that, and the Python side has to match.
    """
    axes = getattr(level, "axes", "")
    c = axes.find("C")
    if c < 0:
        return level.asarray()
    n = int(level.shape[c])
    idx = min(max(ci, 0), max(n - 1, 0))

    # Fast path: one page per channel, channel is the outermost axis.
    pages = getattr(level, "pages", None)
    if c == 0 and pages is not None and len(pages) == n:
        page = pages[idx]
        if page is not None:
            plane = np.asarray(page.asarray())
            if plane.ndim == 3 and plane.shape[0] == 1:
                plane = plane[0]
            if plane.ndim == 2:
                return plane

    return np.take(level.asarray(), idx, axis=c)


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------
def attach_regions(samples: list[U.DonorSample], data_dir: Path | None = None
                   ) -> dict[str, Path]:
    """Map sample name -> `regions/<name>.png`, for samples that have one.

    Samples without a region PNG (every H&E donor, in a pooled run) are simply
    absent from the map; `DapiRegionDataset` then zeroes their auxiliary loss
    rather than inventing labels.
    """
    out: dict[str, Path] = {}
    for s in samples:
        # Derive from the sample's own image path, which discover_samples found
        # under `<base>/images/`, so this is right even when a caller passes a
        # data_dir that does not match where the sample actually came from (the
        # pooled case, where samples arrive from several dirs at once).
        candidates = [s.image_path.parent.parent / "regions" / s.image_path.name]
        if data_dir is not None:
            candidates.append(Path(data_dir) / "regions" / s.image_path.name)
        for p in candidates:
            if p.exists():
                out[s.name] = p
                break
    return out


def discover_samples(data_dir: Path | None = None) -> list[U.DonorSample]:
    """`HE_unet_utils.discover_samples` rooted at the DAPI data dir by default."""
    return U.discover_samples(data_dir if data_dir is not None else DATA_DIR)


class DapiRegionDataset(U.FolliclePatchDataset):
    """Patch sampler returning `(img, follicle_mask, region_label, aux_valid)`.

    Extends the H&E sampler rather than replacing it, so patch geometry,
    positive/negative balance, augmentation and donor round-robin stay
    bit-identical between the two pipelines.

    `aux_valid` is 1.0 for donors with a `regions/` PNG and 0.0 otherwise, which
    is exactly what lets a pooled DAPI+H&E run train the shared Follicle head on
    everything while the auxiliary region head only sees DAPI donors.
    """

    def __init__(self, samples: list[U.DonorSample], region_paths: dict[str, Path],
                 patch: int = PATCH, length: int = 2000, augment: bool = True,
                 pos_frac: float = POS_FRAC, seed: int = 0,
                 n_region_classes: int = len(REGION_CLASSES),
                 merge_wp: bool = False):
        super().__init__(samples, patch=patch, length=length, augment=augment,
                         pos_frac=pos_frac, seed=seed)
        self.n_region_classes = n_region_classes
        self.merge_wp = merge_wp
        self._regions: dict[str, np.ndarray] = {}
        for s in samples:
            p = region_paths.get(s.name)
            if p is None:
                continue
            reg = np.array(Image.open(p).convert("L"), dtype=np.uint8)
            img = self._cache[s.name][0]
            if reg.shape != img.shape:
                raise ValueError(
                    f"{s.name}: region map {reg.shape} does not match image "
                    f"{img.shape} — stale/mismatched export. Re-run "
                    "export_DAPI_training_pairs.groovy.")
            if merge_wp:
                reg = merge_white_pulp(reg)
            if reg.max() >= n_region_classes:
                raise ValueError(
                    f"{s.name}: region map has label {int(reg.max())} but only "
                    f"{n_region_classes} classes are configured.")
            self._regions[s.name] = reg

    def __getitem__(self, idx: int):
        s = self.samples[idx % len(self.samples)]
        img, mask, tissue = self._cache[s.name]
        reg = self._regions.get(s.name)
        positive = self.rng.random() < self.pos_frac and mask.any()
        y, x = self._sample_center(mask, tissue, positive)
        # One origin for all three arrays — they cannot drift apart, and the
        # shared helper pads donors smaller than one patch up to full size.
        y0, x0 = self._window_origin(y, x, *mask.shape)
        img_patch = self._crop_to_patch(img, y0, x0)
        mask_patch = self._crop_to_patch(mask, y0, x0)
        # Zeros are a valid label ("Background"); aux_valid, not the content,
        # is what tells the loss to ignore them.
        reg_patch = (self._crop_to_patch(reg, y0, x0) if reg is not None
                     else np.zeros_like(mask_patch))
        aux_valid = 1.0 if reg is not None else 0.0

        if self.augment:
            # Apply the SAME geometric ops to the region map. Stack it with the
            # mask so one call cannot desynchronise the two label maps, then
            # take the image-only intensity jitter from the parent.
            k = int(self.rng.integers(4))
            img_patch = np.rot90(img_patch, k=k).copy()
            mask_patch = np.rot90(mask_patch, k=k).copy()
            reg_patch = np.rot90(reg_patch, k=k).copy()
            if self.rng.random() < 0.5:
                img_patch = np.fliplr(img_patch).copy()
                mask_patch = np.fliplr(mask_patch).copy()
                reg_patch = np.fliplr(reg_patch).copy()
            scale = float(self.rng.uniform(0.9, 1.1))
            img_patch = np.clip(img_patch.astype(np.float32) * scale, 0, 255).astype(np.uint8)
            gamma = self.rng.uniform(0.85, 1.15)
            img_patch = np.clip(255.0 * (img_patch / 255.0) ** gamma, 0, 255).astype(np.uint8)

        img_t = torch.from_numpy(img_patch).float().unsqueeze(0) / 255.0
        mask_t = torch.from_numpy(mask_patch).float().unsqueeze(0)
        reg_t = torch.from_numpy(reg_patch.astype(np.int64))
        return img_t, mask_t, reg_t, torch.tensor(aux_valid, dtype=torch.float32)


# ---------------------------------------------------------------------------
# Model + loss
# ---------------------------------------------------------------------------
def build_multitask_unet(encoder: str = "efficientnet-b3", pretrained: bool = True,
                         n_region_classes: int = len(REGION_CLASSES)):
    """U-Net with `1 + n_region_classes` output channels.

    Channel 0 is the binary Follicle logit; channels 1.. are the auxiliary
    region logits. `n_region_classes=0` collapses to exactly the H&E model, so
    the same checkpoint machinery covers both.
    """
    return U.build_unet(encoder, pretrained=pretrained, in_channels=1,
                        classes=1 + max(0, int(n_region_classes)))


def multitask_activation(logits: torch.Tensor) -> torch.Tensor:
    """sigmoid on channel 0, softmax across channels 1.. — for `predict_wsi`."""
    if logits.shape[1] == 1:
        return torch.sigmoid(logits)
    return torch.cat([torch.sigmoid(logits[:, :1]),
                      torch.softmax(logits[:, 1:], dim=1)], dim=1)


def follicle_activation(logits: torch.Tensor) -> torch.Tensor:
    """Sigmoid of channel 0 only — pair with `predict_wsi(..., n_out=1)`.

    The model still computes every channel; this just stops `predict_wsi` from
    stitching the auxiliary ones. That matters at slide scale: a 2 µm/px raster
    of the largest section is ~10^8 px, so accumulating 1+6 float32 channels is
    ~3 GB, versus ~430 MB for the one channel anything downstream reads. Used by
    the per-epoch validation and by whole-slide inference.
    """
    return torch.sigmoid(logits[:, :1])


def region_class_weights(region_paths: dict[str, Path], n_classes: int,
                         merge_wp: bool = False) -> list[float]:
    """Inverse-sqrt-frequency class weights over the training region maps.

    RedPulp covers most of the tissue; unweighted cross-entropy would let the
    aux head score well by predicting it almost everywhere. Weights are computed
    once from the *training* donors and stored in the checkpoint so evaluation
    and any later fine-tune reuse the same numbers.
    """
    counts = np.zeros(n_classes, dtype=np.float64)
    for p in region_paths.values():
        reg = np.array(Image.open(p).convert("L"), dtype=np.uint8)
        if merge_wp:
            reg = merge_white_pulp(reg)
        counts += np.bincount(reg.ravel(), minlength=n_classes)[:n_classes]
    if counts.sum() <= 0:
        return [1.0] * n_classes
    freq = counts / counts.sum()
    present = counts > 0
    w = np.zeros(n_classes, dtype=np.float64)
    # Classes absent from the training set get weight 0, not a huge one. Some
    # donors have no Trabeculae or LargeVessel at all; 1/sqrt(clipped-zero) would
    # hand those a ~1000x weight and, once normalised, crush the weights of the
    # classes that actually occur.
    w[present] = 1.0 / np.sqrt(freq[present])
    w[present] = w[present] / w[present].mean()   # comparable loss scale
    return [float(v) for v in w]


def region_loss(logits: torch.Tensor, target: torch.Tensor,
                class_weights: torch.Tensor | None = None,
                eps: float = 1e-6) -> torch.Tensor:
    """Per-sample weighted CE + macro soft-dice over the region channels.

    Returns shape `(B,)`, NOT a scalar — `multitask_loss` needs per-sample values
    so it can zero out donors that have no region labels.
    """
    ce = F.cross_entropy(logits, target, weight=class_weights, reduction="none")
    ce = ce.mean(dim=(1, 2))
    probs = torch.softmax(logits, dim=1)
    onehot = F.one_hot(target, num_classes=logits.shape[1]).permute(0, 3, 1, 2).float()
    num = 2 * (probs * onehot).sum(dim=(2, 3))
    den = (probs + onehot).sum(dim=(2, 3)) + eps
    per_class = 1 - num / den
    # Average over classes PRESENT in this patch's target only. A 512 px patch
    # rarely contains all six regions, and counting an absent class as a full
    # 1.0 dice penalty would make the loss depend mostly on which classes
    # happened to be in frame rather than on how well they were segmented.
    present = (onehot.sum(dim=(2, 3)) > 0).float()
    dice = (per_class * present).sum(dim=1) / present.sum(dim=1).clamp(min=1.0)
    return 0.5 * ce + 0.5 * dice


def multitask_loss(logits: torch.Tensor, mask: torch.Tensor,
                   region: torch.Tensor | None = None,
                   aux_valid: torch.Tensor | None = None,
                   lam: float = 0.3,
                   class_weights: torch.Tensor | None = None) -> torch.Tensor:
    """`follicle_loss + lam * mean(aux_valid * region_loss)`.

    The auxiliary term is averaged over the *whole* batch, not over the valid
    samples only, so a batch that happens to contain mostly H&E donors
    contributes proportionally less auxiliary gradient instead of having a
    handful of DAPI patches speak for it at full weight.
    """
    primary = U.bce_dice_loss(logits[:, :1], mask)
    if region is None or logits.shape[1] <= 1 or lam == 0:
        return primary
    per_sample = region_loss(logits[:, 1:], region, class_weights=class_weights)
    if aux_valid is not None:
        per_sample = per_sample * aux_valid.to(per_sample.dtype)
    return primary + lam * per_sample.mean()


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def region_confusion(pred_labels: np.ndarray, gt_labels: np.ndarray,
                     class_names: list[str] | None = None) -> dict:
    """Per-class IoU plus the Follicle<->PALS confusion, from label maps.

    The Follicle/PALS pair is called out on purpose: DAPI cannot see CD20 vs
    CD3e, so a model can post a respectable follicle count while systematically
    swapping the two. Reporting only aggregate IoU would hide that.
    """
    names = class_names if class_names is not None else REGION_CLASSES
    cm = confusion_counts(pred_labels, gt_labels, len(names))
    return region_metrics_from_cm(cm, names)


def confusion_counts(pred_labels: np.ndarray, gt_labels: np.ndarray,
                     n_classes: int) -> np.ndarray:
    """`(n, n)` confusion counts; rows = ground truth, cols = prediction.

    Split out from `region_confusion` because confusion matrices ADD, which is
    what lets region metrics be accumulated over image bands on a slide too
    large to hold all K probability channels at once.
    """
    pred = np.asarray(pred_labels).ravel().astype(np.int64)
    gt = np.asarray(gt_labels).ravel().astype(np.int64)
    return np.bincount(gt * n_classes + pred,
                       minlength=n_classes * n_classes
                       ).reshape(n_classes, n_classes).astype(np.int64)


def region_metrics_from_cm(cm: np.ndarray, class_names: list[str]) -> dict:
    """Per-class IoU + Follicle<->PALS confusion from an accumulated matrix."""
    names = list(class_names)
    cm = np.asarray(cm, dtype=np.float64)

    # A class absent from BOTH the prediction and the truth reports None, not
    # 1.0. Several donors have no Trabeculae or LargeVessel at all, and scoring
    # those as a perfect 1.0 would inflate the headline mean by a class the model
    # was never asked to find. `iou_score` keeps its both-empty==1.0 convention
    # for the binary Follicle head; this is a different question.
    per_class: dict[str, float | None] = {}
    for i, name in enumerate(names):
        inter = cm[i, i]
        union = cm[i, :].sum() + cm[:, i].sum() - inter
        per_class[name] = float(inter / union) if union > 0 else None

    scored = [v for v in per_class.values() if v is not None]
    out = {"per_class_iou": per_class,
           "mean_iou": float(np.mean(scored)) if scored else 0.0,
           "n_classes_scored": len(scored),
           "confusion": cm.astype(np.int64).tolist(),
           "class_names": list(names)}
    if "Follicle" in names and "PALS" in names:
        f, p = names.index("Follicle"), names.index("PALS")
        gt_f, gt_p = cm[f, :].sum(), cm[p, :].sum()
        out["follicle_as_pals"] = float(cm[f, p] / gt_f) if gt_f > 0 else 0.0
        out["pals_as_follicle"] = float(cm[p, f] / gt_p) if gt_p > 0 else 0.0
    return out


def region_confusion_banded(model, img: np.ndarray, gt_regions: np.ndarray,
                            device, n_region_classes: int,
                            class_names: list[str], patch: int = PATCH,
                            overlap: int = 64, band_rows: int = 4096,
                            tta: bool = True) -> dict:
    """Region metrics computed in horizontal bands, accumulating the confusion.

    Holding all K probability channels for a whole 2 µm/px slide is ~3 GB;
    banding caps it at `band_rows` at a time.

    Each band OWNS a disjoint range of rows but is predicted with a `patch`-row
    margin on either side, so interior rows still see their neighbours. Only the
    owned rows are counted, so every pixel contributes exactly once and the
    accumulated matrix totals `gt_regions.size` — the naive version, where bands
    simply overlapped, double-counted the seams.
    """
    img = np.asarray(img)
    h = img.shape[0]
    n = len(class_names)
    cm = np.zeros((n, n), dtype=np.int64)
    band_rows = max(1, int(band_rows))
    for own0 in range(0, h, band_rows):
        own1 = min(h, own0 + band_rows)
        y0 = max(0, own0 - patch)
        y1 = min(h, own1 + patch)
        prob = U.predict_wsi(model, img[y0:y1], device, patch=patch,
                             overlap=overlap, tta=tta,
                             n_out=1 + n_region_classes,
                             activation=multitask_activation)
        pred = np.argmax(prob[1:], axis=0).astype(np.uint8)
        del prob
        cm += confusion_counts(pred[own0 - y0:own1 - y0],
                               gt_regions[own0:own1], n)
    return region_metrics_from_cm(cm, class_names)


def check_region_integrity(samples: list[U.DonorSample],
                           region_paths: dict[str, Path],
                           min_follicle_frac: float = 1e-4) -> list[dict]:
    """Per-donor region-map sanity check, run before training.

    Deliberately NOT `check_mask_integrity` with different numbers: that one
    flags a high positive fraction as suspicious, which is wrong here — RedPulp
    legitimately covers most of the tissue. What actually goes wrong with a
    region export is a missing class or a map that is one flat label, so those
    are what this checks.
    """
    reports = []
    for s in samples:
        p = region_paths.get(s.name)
        if p is None:
            reports.append({"name": s.name, "status": "no_regions",
                            "present_classes": [], "follicle_frac": 0.0})
            continue
        reg = np.array(Image.open(p).convert("L"), dtype=np.uint8)
        counts = np.bincount(reg.ravel(), minlength=len(REGION_CLASSES))
        present = [REGION_CLASSES[i] for i in range(len(REGION_CLASSES))
                   if i < len(counts) and counts[i] > 0]
        foll_frac = float(counts[FOLLICLE_INDEX] / reg.size) if len(counts) > FOLLICLE_INDEX else 0.0
        if len(present) <= 1:
            status = "single_label"
        elif foll_frac < min_follicle_frac:
            status = "no_follicle"
        else:
            status = "ok"
        reports.append({"name": s.name, "status": status,
                        "present_classes": present,
                        "follicle_frac": foll_frac})
    return reports


# ---------------------------------------------------------------------------
# Sidecar helpers
# ---------------------------------------------------------------------------
def read_meta(sample: U.DonorSample) -> dict:
    return json.loads(Path(sample.meta_path).read_text())


def sidecar_transform(meta: dict) -> dict:
    """Pull `dapi_transform` out of a sidecar, erroring loudly if it is absent.

    Inference must reuse the exact `lo`/`hi` the export baked into the training
    PNG; silently recomputing them would put the model on a different intensity
    scale than it trained on.
    """
    t = meta.get("dapi_transform")
    if not t or "lo" not in t or "hi" not in t:
        raise KeyError("sidecar has no usable `dapi_transform` — re-export with "
                       "the current export_DAPI_training_pairs.groovy")
    return t
