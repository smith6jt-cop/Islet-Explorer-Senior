"""Shared utilities for the H&E follicle U-Net pipeline.

Dataset, model factory, augmentation, loss, and post-processing. Imported by
both `train_HE_unet.py` and `predict_HE_unet.py`.
"""
from __future__ import annotations

import json
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from PIL import Image
from torch.utils.data import Dataset

CODE_ROOT = Path(__file__).resolve().parent.parent
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from identity import SampleManifestError, resolve_sample_id  # noqa: E402
from paths import DATA_ROOT, OUTPUT_ROOT  # noqa: E402

# Compatibility for callers that use U.PROJECT for external data paths.
PROJECT = DATA_ROOT
DATA_DIR = DATA_ROOT / "data"
IMAGES_DIR = DATA_DIR / "images"
MASKS_DIR = DATA_DIR / "masks"
META_DIR = DATA_DIR / "meta"
MODELS_DIR = OUTPUT_ROOT / "analysis" / "HE_unet_models"
PRED_DIR = OUTPUT_ROOT / "analysis" / "HE_unet_predictions"

TARGET_UM_PER_PX = 2.0
PATCH = 512
POS_FRAC = 0.5

# ---------------------------------------------------------------------------
# RGB → hematoxylin (H) deconvolution
# ---------------------------------------------------------------------------
# Reference H&E optical-density unit vectors. These MUST stay numerically
# identical to the deconvolution in `export_HE_training_pairs.groovy` — any
# change here has to be made in both files in the same commit (see the
# export↔Python parity check in HE_unet_README.md).
_H_REF = np.array([0.5626, 0.7201, 0.4062])
_E_REF = np.array([0.2159, 0.8012, 0.5581])
_OD_WHITE = 255.0   # white reference for OD = -log10(max(channel, 1) / white)
_CH_MAX = 2.0       # H concentration clip / display normalisation ceiling


def _h_projection() -> np.ndarray:
    """First row of (MMᵀ)⁻¹M for M = [Ĥ; Ê]: the OD→H-concentration vector."""
    Hn = _H_REF / np.linalg.norm(_H_REF)
    En = _E_REF / np.linalg.norm(_E_REF)
    hh, he, ee = Hn @ Hn, Hn @ En, En @ En
    det = hh * ee - he * he
    a00, a01 = ee / det, -he / det
    return a00 * Hn + a01 * En


def rgb_to_h(rgb: np.ndarray) -> np.ndarray:
    """Deconvolve an RGB H&E image to a single-channel hematoxylin (H) image.

    Inference-side counterpart of the deconvolution baked into
    `export_HE_training_pairs.groovy`; the two transforms MUST agree to the byte
    (high hematoxylin → dark). Returns a uint8 `(H, W)` array.

    OD = -log10(max(channel, 1) / 255); cH = w·OD clipped to [0, _CH_MAX];
    intensity = floor(255·(1 − cH/_CH_MAX)).
    """
    rgb = np.asarray(rgb)
    if rgb.ndim != 3 or rgb.shape[-1] < 3:
        raise ValueError(f"rgb_to_h expects an (H, W, 3) image, got {rgb.shape}")
    w = _h_projection()
    # OD lookup over the 8-bit range — matches the export's per-pixel LUT.
    lut = -np.log10(np.maximum(np.arange(256), 1) / _OD_WHITE)
    r, g, b = rgb[..., 0].astype(np.intp), rgb[..., 1].astype(np.intp), rgb[..., 2].astype(np.intp)
    cH = w[0] * lut[r] + w[1] * lut[g] + w[2] * lut[b]
    cH = np.clip(cH, 0.0, _CH_MAX)
    inten = np.clip(255.0 * (1.0 - cH / _CH_MAX), 0.0, 255.0)
    return np.floor(inten).astype(np.uint8)


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------
@dataclass
class DonorSample:
    name: str
    image_path: Path
    mask_path: Path
    meta_path: Path
    tissue_path: Path | None = None
    sample_id: str | None = None

    @property
    def meta(self) -> dict:
        return json.loads(self.meta_path.read_text())

    @property
    def identity(self) -> str:
        """Canonical sample ID declared by metadata or the exact manifest."""

        if self.sample_id:
            return self.sample_id
        return resolve_sample_id(self.name)


def _identity_from_metadata(name: str, meta_path: Path) -> str:
    """Resolve an export without interpreting any part of its filename."""

    metadata = json.loads(meta_path.read_text())
    declared = str(metadata.get("sample_id", "")).strip()
    if declared:
        return declared

    resolved = set()
    attempted = []
    for value in (metadata.get("source_name"), metadata.get("image"), name):
        if value is None:
            continue
        alias = str(value).strip()
        if not alias or alias in attempted:
            continue
        attempted.append(alias)
        try:
            resolved.add(resolve_sample_id(alias))
        except SampleManifestError:
            pass
    if len(resolved) == 1:
        return resolved.pop()
    if len(resolved) > 1:
        raise SampleManifestError(
            f"{meta_path} aliases resolve to conflicting samples: {sorted(resolved)}")
    raise SampleManifestError(
        f"{meta_path} has no sample_id and none of its exact aliases are in the "
        f"sample manifest: {attempted}")


def discover_samples(data_dir: Path | None = None) -> list[DonorSample]:
    """Find (image, mask, meta) triples under `data_dir` (defaults to DATA_DIR).

    Passing an explicit dir lets the smoke test / CI point at a synthetic
    dataset without touching the repo's `data/`. The RGB-derived tissue mask
    (`data/tissue/<name>.png`) is attached when present.
    """
    base = Path(data_dir) if data_dir is not None else DATA_DIR
    images_dir, masks_dir, meta_dir = base / "images", base / "masks", base / "meta"
    tissue_dir = base / "tissue"
    samples = []
    for img in sorted(images_dir.glob("*.png")):
        mask = masks_dir / img.name
        meta = meta_dir / (img.stem + ".json")
        if mask.exists() and meta.exists():
            tissue = tissue_dir / img.name
            samples.append(DonorSample(
                img.stem, img, mask, meta,
                tissue if tissue.exists() else None,
                _identity_from_metadata(img.stem, meta)))
    return samples


def _load_pair(sample: DonorSample) -> tuple[np.ndarray, np.ndarray]:
    """Load the single-channel H image (2-D) and the binary follicle mask."""
    img = np.array(Image.open(sample.image_path).convert("L"), dtype=np.uint8)
    mask = np.array(Image.open(sample.mask_path).convert("L"), dtype=np.uint8)
    mask = (mask > 127).astype(np.uint8)
    if img.shape[:2] != mask.shape:
        raise ValueError(f"{sample.name}: image {img.shape[:2]} vs mask {mask.shape}")
    return img, mask


def _tissue_mask(img: np.ndarray) -> np.ndarray:
    """Cheap whitespace mask: anything not near-white is tissue.

    Fallback only — used when the export's RGB-derived tissue PNG is absent.
    Accepts a 3-D RGB image (thresholds mean(R,G,B)) or the 2-D H image (high
    H = dark = tissue). CAVEAT: deriving tissue from H alone under-counts red
    pulp (low H, near-white), which is exactly why the export writes a proper
    RGB-derived tissue mask.
    """
    gray = img.mean(axis=2) if img.ndim == 3 else img
    return gray < 220


def _load_tissue(sample: DonorSample, h_img: np.ndarray) -> np.ndarray:
    """Tissue mask from the export's RGB-derived PNG when present, else the H
    fallback (with the red-pulp caveat documented in `_tissue_mask`).

    A present-but-shape-mismatched tissue PNG is a data-contract violation
    (stale / mismatched export) and raises rather than silently falling back —
    quietly deriving tissue from H would hide the mismatch and bias negative
    sampling.
    """
    if sample.tissue_path is not None and sample.tissue_path.exists():
        t = np.array(Image.open(sample.tissue_path).convert("L"), dtype=np.uint8)
        if t.shape != h_img.shape:
            raise ValueError(
                f"{sample.name}: tissue mask {t.shape} does not match H image "
                f"{h_img.shape} — stale/mismatched export. Re-run "
                "export_HE_training_pairs.groovy to regenerate data/.")
        return t > 127
    return _tissue_mask(h_img)


class FolliclePatchDataset(Dataset):
    """Patch-based sampler.

    Each `__getitem__` returns one (image, mask) tile, biased so `POS_FRAC` of
    tiles are centered on a follicle pixel and the rest are random tissue tiles.
    """

    def __init__(self, samples: list[DonorSample], patch: int = PATCH,
                 length: int = 2000, augment: bool = True,
                 pos_frac: float = POS_FRAC, seed: int = 0):
        self.samples = samples
        self.patch = patch
        self.length = length
        self.augment = augment
        self.pos_frac = pos_frac
        self.rng = np.random.default_rng(seed)

        # Cache the on-disk single-channel H image, the follicle mask, and the
        # disk tissue mask. No rgb_to_h call in the training path — the export
        # already wrote H images.
        self._cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for s in samples:
            img, mask = _load_pair(s)
            tissue = _load_tissue(s, img)
            self._cache[s.name] = (img, mask, tissue)

    def __len__(self) -> int:
        return self.length

    def _sample_center(self, mask: np.ndarray, tissue: np.ndarray,
                       positive: bool) -> tuple[int, int]:
        h, w = mask.shape
        half = self.patch // 2
        if positive:
            candidate = mask
        else:
            # True negatives: tissue away from any follicle. Fall back to plain
            # tissue if the slide is densely follicular, then to the centre.
            candidate = tissue & ~mask.astype(bool)
        ys, xs = np.where(candidate)
        if len(ys) == 0:
            ys, xs = np.where(tissue)
        if len(ys) == 0:
            return half, half
        i = self.rng.integers(len(ys))
        # Return the raw centre. Clamping to [half, h-half] here is wrong when
        # the image is SMALLER than the patch: np.clip with a_min > a_max
        # silently returns a_max, which then slices from a negative index and
        # yields an undersized tile. `_crop_to_patch` does the clamping (and
        # padding) correctly instead.
        return int(ys[i]), int(xs[i])

    def _window_origin(self, y: int, x: int, h: int, w: int) -> tuple[int, int]:
        """Top-left corner of the patch window, clamped into the image.

        Shared by every array cropped for one sample so the image, the follicle
        mask and (in the DAPI subclass) the region map cannot drift apart.
        """
        half = self.patch // 2
        y0 = min(max(y - half, 0), max(0, h - self.patch))
        x0 = min(max(x - half, 0), max(0, w - self.patch))
        return y0, x0

    def _crop_to_patch(self, arr: np.ndarray, y0: int, x0: int) -> np.ndarray:
        """Crop at (y0, x0) and pad to EXACTLY (patch, patch).

        Padding is what makes donors whose export is smaller than one patch
        usable at all. Without it the crop silently comes back undersized and
        the failure surfaces much later as a DataLoader collate error — and only
        when such a donor is batched alongside a normal one, which is exactly
        what a pooled DAPI + H&E run does.
        """
        p = self.patch
        out = arr[y0:y0 + p, x0:x0 + p]
        ph, pw = p - out.shape[0], p - out.shape[1]
        if ph <= 0 and pw <= 0:
            return out
        # `reflect` needs pad < dim on that axis; fall back to edge replication
        # for donors smaller than half a patch rather than raising.
        mode = ("reflect" if (out.shape[0] > ph and out.shape[1] > pw
                              and min(out.shape[:2]) > 1) else "edge")
        return np.pad(out, ((0, max(0, ph)), (0, max(0, pw))), mode=mode)

    def _augment(self, img: np.ndarray, mask: np.ndarray
                 ) -> tuple[np.ndarray, np.ndarray]:
        # Flips + 90° rotations
        k = int(self.rng.integers(4))
        img = np.rot90(img, k=k).copy()
        mask = np.rot90(mask, k=k).copy()
        if self.rng.random() < 0.5:
            img = np.fliplr(img).copy()
            mask = np.fliplr(mask).copy()
        # Intensity jitter on the single H channel — cheap surrogate for the
        # day-to-day staining variation that survives Macenko normalization.
        scale = float(self.rng.uniform(0.9, 1.1))
        img = np.clip(img.astype(np.float32) * scale, 0, 255).astype(np.uint8)
        # Mild gamma
        gamma = self.rng.uniform(0.85, 1.15)
        img = np.clip(255.0 * (img / 255.0) ** gamma, 0, 255).astype(np.uint8)
        return img, mask

    def __getitem__(self, idx: int) -> tuple[torch.Tensor, torch.Tensor]:
        s = self.samples[idx % len(self.samples)]
        img, mask, tissue = self._cache[s.name]
        positive = self.rng.random() < self.pos_frac and mask.any()
        y, x = self._sample_center(mask, tissue, positive)
        y0, x0 = self._window_origin(y, x, *mask.shape)
        img_patch = self._crop_to_patch(img, y0, x0)
        mask_patch = self._crop_to_patch(mask, y0, x0)
        if self.augment:
            img_patch, mask_patch = self._augment(img_patch, mask_patch)
        # Single-channel H image: (H, W) → (1, H, W) via unsqueeze, not permute.
        img_t = torch.from_numpy(img_patch).float().unsqueeze(0) / 255.0
        mask_t = torch.from_numpy(mask_patch).float().unsqueeze(0)
        return img_t, mask_t


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
def build_unet(encoder: str = "efficientnet-b3", pretrained: bool = True,
               in_channels: int = 1, classes: int = 1) -> nn.Module:
    """Return a U-Net with an ImageNet-pretrained encoder.

    Uses segmentation_models_pytorch if available; otherwise falls back to a
    plain-PyTorch U-Net so the pipeline still runs (worse, but trainable).
    `pretrained=False` skips the ImageNet weight download — used by the smoke
    test / CI so a crash-check needs no network access.

    `classes` > 1 gives a multi-head output: the DAPI pipeline uses channel 0 as
    the binary Follicle logit and channels 1.. as an auxiliary region softmax
    (see DAPI_unet_utils). Defaults reproduce the single-channel H&E model
    exactly, so existing checkpoints keep loading.
    """
    try:
        import segmentation_models_pytorch as smp
        # Single-channel H input. smp adapts the ImageNet-pretrained first conv
        # to one input channel by summing the RGB filter weights.
        return smp.Unet(
            encoder_name=encoder,
            encoder_weights="imagenet" if pretrained else None,
            in_channels=in_channels,
            classes=classes,
            activation=None,
        )
    except ImportError:
        return _FallbackUNet(in_ch=in_channels, n_classes=classes)


class _ConvBlock(nn.Module):
    def __init__(self, c_in: int, c_out: int):
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv2d(c_in, c_out, 3, padding=1, bias=False),
            nn.BatchNorm2d(c_out),
            nn.ReLU(inplace=True),
            nn.Conv2d(c_out, c_out, 3, padding=1, bias=False),
            nn.BatchNorm2d(c_out),
            nn.ReLU(inplace=True),
        )

    def forward(self, x):
        return self.block(x)


class _FallbackUNet(nn.Module):
    def __init__(self, in_ch: int = 1, n_classes: int = 1):
        super().__init__()
        chs = [32, 64, 128, 256, 512]
        self.enc = nn.ModuleList()
        c_in = in_ch  # single-channel nuclear image (H or inverted DAPI)
        for c in chs:
            self.enc.append(_ConvBlock(c_in, c))
            c_in = c
        self.pool = nn.MaxPool2d(2)
        self.up = nn.ModuleList(
            nn.ConvTranspose2d(chs[i + 1], chs[i], 2, stride=2) for i in range(len(chs) - 1)
        )
        self.dec = nn.ModuleList(
            _ConvBlock(chs[i + 1], chs[i]) for i in range(len(chs) - 1)
        )
        self.head = nn.Conv2d(chs[0], n_classes, 1)

    def forward(self, x):
        skips = []
        for i, block in enumerate(self.enc):
            x = block(x)
            if i < len(self.enc) - 1:
                skips.append(x)
                x = self.pool(x)
        for i in reversed(range(len(self.up))):
            x = self.up[i](x)
            x = torch.cat([x, skips[i]], dim=1)
            x = self.dec[i](x)
        return self.head(x)


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------
def dice_loss(logits: torch.Tensor, targets: torch.Tensor, eps: float = 1e-6
              ) -> torch.Tensor:
    probs = torch.sigmoid(logits)
    num = 2 * (probs * targets).sum(dim=(2, 3))
    den = (probs + targets).sum(dim=(2, 3)) + eps
    return (1 - num / den).mean()


def bce_dice_loss(logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
    return 0.5 * F.binary_cross_entropy_with_logits(logits, targets) + 0.5 * dice_loss(logits, targets)


# ---------------------------------------------------------------------------
# Tiled inference with cosine-weighted blending
# ---------------------------------------------------------------------------
def _cosine_window(size: int) -> np.ndarray:
    x = np.linspace(-np.pi / 2, np.pi / 2, size, dtype=np.float32)
    w = np.cos(x) ** 2
    return np.outer(w, w)


@torch.no_grad()
def predict_wsi(model: nn.Module, himg: np.ndarray, device: torch.device,
                patch: int = PATCH, overlap: int = 64, tta: bool = True,
                n_out: int = 1, activation=None) -> np.ndarray:
    """Tile, infer, stitch. Returns a float32 probability map.

    Shape is `(H, W)` when `n_out == 1` (the H&E case) and `(n_out, H, W)`
    otherwise — the DAPI multitask model stitches its Follicle channel and its
    auxiliary region channels in one pass.

    `activation` maps a `(1, n_out, P, P)` logit tensor to probabilities; it
    defaults to `torch.sigmoid`. The DAPI model passes a mixed
    sigmoid(channel 0) + softmax(channels 1..) activation. Activations are
    applied per TTA op and the *probabilities* averaged, not the logits.

    Consumes a single-channel nuclear-density image (2-D, or `(H, W, 1)`).
    Callers convert at the boundary (`predict_HE_unet.py` via `rgb_to_h`,
    `predict_DAPI_unet.py` via `dapi_to_uint8`); this function does NOT convert.
    """
    if activation is None:
        activation = torch.sigmoid
    model.eval()
    himg = np.asarray(himg)
    if himg.ndim == 3 and himg.shape[-1] == 1:
        himg = himg[..., 0]  # tolerate (H, W, 1)
    elif himg.ndim != 2:
        # Reject RGB / multi-channel: silently taking channel 0 would run the
        # model on the red channel and look superficially valid. Callers must
        # convert RGB->H via rgb_to_h before calling predict_wsi.
        raise ValueError(
            f"predict_wsi expects a single-channel H image (2-D or (H, W, 1)), "
            f"got shape {himg.shape}. Convert RGB->H via rgb_to_h first.")
    h, w = himg.shape[:2]
    stride = patch - overlap
    if stride <= 0:
        # Otherwise this surfaces much later as a bare ZeroDivisionError from the
        # padding arithmetic, which says nothing about the actual mistake.
        raise ValueError(
            f"predict_wsi needs overlap < patch, got patch={patch} "
            f"overlap={overlap}.")
    pad_h = (-(h - patch) % stride) if h > patch else (patch - h)
    pad_w = (-(w - patch) % stride) if w > patch else (patch - w)
    padded = np.pad(himg, ((0, pad_h), (0, pad_w)), mode="reflect")
    H, W = padded.shape[:2]

    prob = np.zeros((n_out, H, W), dtype=np.float32)
    weight = np.zeros((H, W), dtype=np.float32)
    # cos^2 is exactly 0 at the tile edge. A pixel on the outer border of the
    # padded image is covered only by that zero, so its accumulated weight is 0
    # and the `maximum(weight, 1e-6)` guard below turns it into probability 0
    # rather than a prediction. Harmless for a binary map (the slide border is
    # not follicle), but it makes the multi-class softmax channels sum to 0
    # instead of 1 along the border. Flooring the window keeps the blend a true
    # convex combination everywhere at a negligible cost to the taper.
    win = np.maximum(_cosine_window(patch), 1e-3)

    ys = list(range(0, H - patch + 1, stride))
    xs = list(range(0, W - patch + 1, stride))
    if ys[-1] != H - patch:
        ys.append(H - patch)
    if xs[-1] != W - patch:
        xs.append(W - patch)

    tta_ops = [
        (lambda x: x, lambda x: x),
        (lambda x: torch.flip(x, [-1]), lambda x: torch.flip(x, [-1])),
        (lambda x: torch.flip(x, [-2]), lambda x: torch.flip(x, [-2])),
        (lambda x: torch.rot90(x, 1, [-2, -1]), lambda x: torch.rot90(x, -1, [-2, -1])),
    ] if tta else [(lambda x: x, lambda x: x)]

    for y in ys:
        for x in xs:
            tile = padded[y:y + patch, x:x + patch]
            # (patch, patch) → (1, 1, patch, patch)
            t = torch.from_numpy(tile).float().unsqueeze(0).unsqueeze(0) / 255.0
            t = t.to(device)
            acc = torch.zeros((1, n_out, patch, patch), device=device)
            for fwd, inv in tta_ops:
                logits = model(fwd(t))
                acc = acc + inv(activation(logits))
            # (1, n_out, P, P) -> (n_out, P, P); win broadcasts over the channel axis.
            p = (acc / len(tta_ops))[0].cpu().numpy()
            prob[:, y:y + patch, x:x + patch] += p * win
            weight[y:y + patch, x:x + patch] += win

    prob = prob / np.maximum(weight, 1e-6)
    prob = prob[:, :h, :w]
    # Preserve the historical 2-D return for the single-output H&E model.
    return prob[0] if n_out == 1 else prob


# ---------------------------------------------------------------------------
# Post-processing: probability map → instance polygons → records
# ---------------------------------------------------------------------------
def extract_follicle_records(prob: np.ndarray, px_size_um: float,
                             threshold: float = 0.5,
                             min_area_um2: float = 10_000,
                             split_touching: bool = True,
                             image_name: str = "",
                             classification: str = "Follicle",
                             ) -> tuple[list[dict], list[list[tuple[float, float]]]]:
    """Convert a probability map into per-instance records + polygons.

    Records use the same column schema as 27_donor_spleen_measurements1.csv so
    they can be concatenated and fed to H18 unchanged. `classification` sets the
    class string written into the records (the DAPI pipeline reuses this for
    region classes other than Follicle).

    Polygon coordinates are in µm. `write_geojson(scale=...)` converts back to
    pixels when the output has to re-import into QuPath.
    """
    from skimage.measure import find_contours, label, regionprops
    from skimage.morphology import closing, disk, remove_small_objects

    mask = prob > threshold
    mask = closing(mask, disk(3))
    min_area_px = max(1, int(min_area_um2 / (px_size_um ** 2)))
    # max_size in skimage ≥0.26 removes objects with size <= max_size
    # (old min_size removed size < min_size); subtract 1 to match.
    mask = remove_small_objects(mask, max_size=max(0, min_area_px - 1))

    if split_touching:
        instances = _watershed_split(mask, px_size_um)
    else:
        instances = label(mask, connectivity=2)

    records, polygons = [], []
    for r in regionprops(instances):
        area_um2 = r.area * px_size_um ** 2
        if area_um2 < min_area_um2:
            continue
        # Centroid in µm
        cy_um = r.centroid[0] * px_size_um
        cx_um = r.centroid[1] * px_size_um
        # Contour for GeoJSON / perimeter — use the full-image binary so we get
        # one outer contour per instance in image coordinates.
        comp = (instances == r.label)
        contours = find_contours(comp.astype(np.uint8), 0.5)
        if not contours:
            continue
        contour = max(contours, key=len)
        polygon_um = [(float(c[1] * px_size_um), float(c[0] * px_size_um))
                      for c in contour]
        perimeter_um = _polygon_perimeter(polygon_um)
        records.append({
            "Image": image_name,
            "Object ID": f"unet-{image_name}-{r.label}",
            "Object type": "Annotation",
            "Name": "",
            "Classification": classification,
            "Parent": "Root object (Image)",
            "ROI": "Polygon",
            "Centroid X µm": round(cx_um, 2),
            "Centroid Y µm": round(cy_um, 2),
            "Area µm^2": round(area_um2, 2),
            "Perimeter µm": round(perimeter_um, 2),
        })
        polygons.append(polygon_um)
    return records, polygons


def _watershed_split(mask: np.ndarray, px_size_um: float) -> np.ndarray:
    """Watershed-on-distance-transform to separate touching follicles."""
    from scipy.ndimage import distance_transform_edt, label as ndi_label
    from skimage.feature import peak_local_max
    from skimage.measure import label as sk_label
    from skimage.segmentation import watershed

    if not mask.any():
        return np.zeros(mask.shape, dtype=np.int32)
    dist = distance_transform_edt(mask)
    # Minimum centre-to-centre separation: ~50 µm.
    min_dist_px = max(3, int(50 / px_size_um))
    coords = peak_local_max(dist, min_distance=min_dist_px, labels=mask)
    if len(coords) == 0:
        # Flat distance plateau or single tiny blob — skip watershed, just label.
        return sk_label(mask, connectivity=2)
    markers = np.zeros(mask.shape, dtype=np.int32)
    for i, (yy, xx) in enumerate(coords, start=1):
        markers[yy, xx] = i
    markers, n_markers = ndi_label(markers > 0)
    if n_markers == 0:
        return sk_label(mask, connectivity=2)
    return watershed(-dist, markers, mask=mask)


def _polygon_perimeter(poly: list[tuple[float, float]]) -> float:
    if len(poly) < 2:
        return 0.0
    arr = np.asarray(poly)
    d = np.diff(arr, axis=0, append=arr[:1])
    return float(np.linalg.norm(d, axis=1).sum())


def write_geojson(polygons: list[list[tuple[float, float]]], path: Path,
                  image_name: str = "", classification: str = "Follicle",
                  scale: float = 1.0) -> None:
    """Write instance polygons as a GeoJSON FeatureCollection.

    Input polygons are in µm (as produced by `extract_follicle_records`).
    `scale` multiplies every coordinate on the way out: leave it at 1.0 to keep
    µm, or pass `1 / px_size_um` to emit **full-resolution pixel** coordinates,
    which is what QuPath's *Objects → Import objects* expects. Existing callers
    keep the historical µm output.
    """
    features = []
    for i, poly in enumerate(polygons):
        ring = poly + [poly[0]] if poly[0] != poly[-1] else poly
        if scale != 1.0:
            ring = [(x * scale, y * scale) for x, y in ring]
        features.append({
            "type": "Feature",
            "id": f"unet-{image_name}-{i}",
            "geometry": {"type": "Polygon", "coordinates": [list(map(list, ring))]},
            "properties": {"classification": {"name": classification}},
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


# ---------------------------------------------------------------------------
# Donor-level k-fold split (stratified by genotype)
# ---------------------------------------------------------------------------
def donor_id_from_name(name: str) -> str:
    """Resolve a legacy caller's exact alias through the sample manifest."""

    return resolve_sample_id(name)


def donor_kfold(samples: list[DonorSample], k: int = 5,
                genotype_xlsx: Path | None = None,
                seed: int = 0) -> list[tuple[list[DonorSample], list[DonorSample]]]:
    """Genotype-stratified k-fold over donors. Falls back to plain shuffle
    when no genotype file is provided.

    Splitting is by *donor* (`donor_id_from_name`), not by file — several files
    of the same donor always land in the same fold.
    """
    rng = np.random.default_rng(seed)
    geno_map: dict[str, str] = {}
    if genotype_xlsx and Path(genotype_xlsx).exists():
        import pandas as pd
        g = pd.read_excel(genotype_xlsx)
        col = next((c for c in g.columns if "rs3184504" in c), None)
        if col is not None:
            for _, row in g.iterrows():
                sid = str(row.get("Sample ID", "")).strip()
                if sid:
                    geno_map[sid] = str(row[col])

    # Group files by donor FIRST, then stratify and deal whole donors. Dealing
    # individual files round-robin would split a multi-file donor across folds
    # even with a correct donor id, which is the leakage this split exists to
    # prevent. `dict` preserves insertion order, so the grouping is
    # deterministic given sorted `samples` and only `rng.shuffle` varies it.
    by_donor: dict[str, list[DonorSample]] = {}
    for s in samples:
        by_donor.setdefault(s.identity, []).append(s)

    bins: dict[str, list[str]] = {}
    for donor in by_donor:
        bins.setdefault(geno_map.get(donor, "unknown"), []).append(donor)

    folds: list[list[DonorSample]] = [[] for _ in range(k)]
    for stratum in bins.values():
        rng.shuffle(stratum)
        for i, donor in enumerate(stratum):
            folds[i % k].extend(by_donor[donor])

    splits = []
    for i in range(k):
        val = folds[i]
        train = [s for j, f in enumerate(folds) if j != i for s in f]
        splits.append((train, val))
    return splits


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def iou_score(pred: np.ndarray, target: np.ndarray) -> float:
    pred_b = pred > 0.5
    target_b = target > 0.5
    inter = (pred_b & target_b).sum()
    union = (pred_b | target_b).sum()
    return float(inter / union) if union else 1.0


def instance_metrics(prob: np.ndarray, gt_mask: np.ndarray, px_size_um: float,
                     threshold: float = 0.5, min_area_um2: float = 10_000.0
                     ) -> dict:
    """Deployment-aligned follicle metrics: predicted vs ground-truth instance
    count and mean area (µm²).

    Predictions go through the same watershed + min-area post-processing as
    `predict_HE_unet.py` (`extract_follicle_records`). Ground-truth instances
    are the connected components of the mask with the SAME min-area floor
    applied — otherwise the GT count is inflated by sub-threshold specks and the
    count concordance H18 cares about is biased. Pure (arrays in, dict out) so
    it is unit-testable without a model.
    """
    from skimage.measure import label as sk_label, regionprops

    gt_b = (gt_mask > 0).astype(np.uint8)
    gt_areas = [r.area * px_size_um ** 2 for r in regionprops(sk_label(gt_b, connectivity=2))]
    gt_areas = [a for a in gt_areas if a >= min_area_um2]

    records, _ = extract_follicle_records(
        prob, px_size_um=px_size_um, threshold=threshold,
        min_area_um2=min_area_um2, image_name="")
    pred_areas = [r["Area µm^2"] for r in records]

    return {
        "gt_count": len(gt_areas),
        "pred_count": len(pred_areas),
        "gt_mean_area_um2": float(np.mean(gt_areas)) if gt_areas else 0.0,
        "pred_mean_area_um2": float(np.mean(pred_areas)) if pred_areas else 0.0,
    }


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
def set_seed(seed: int, deterministic: bool = True) -> None:
    """Seed python/numpy/torch (+CUDA) for reproducible runs.

    With `num_workers=0` (the default) this makes CPU runs bit-exact across
    reruns at a fixed seed. On GPU, cuDNN/AMP leave minor nondeterminism;
    `deterministic=True` sets the cuDNN flags, but we intentionally do NOT call
    `torch.use_deterministic_algorithms` — it needs CUBLAS_WORKSPACE_CONFIG set
    before the process starts and can hard-error on SMP upsampling kernels.
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if deterministic:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def seed_worker(worker_id: int) -> None:
    """DataLoader `worker_init_fn`: reseed numpy/random per worker from torch's
    per-worker base seed. Only needed when `num_workers > 0` — the global
    `set_seed` does not reach forked workers."""
    worker_seed = torch.initial_seed() % 2 ** 32
    np.random.seed(worker_seed)
    random.seed(worker_seed)


# ---------------------------------------------------------------------------
# Data integrity
# ---------------------------------------------------------------------------
def check_mask_integrity(samples: list[DonorSample],
                         max_pos_frac_warn: float = 0.6,
                         single_blob_area_frac: float = 0.5) -> list[dict]:
    """Per-donor mask sanity check, run before training.

    Returns one report dict per donor with `pos_frac`, `n_components`, and a
    `status`. Flags `empty` masks (no Follicle pixels) and shapes that suggest a
    coarse region outline rather than per-follicle annotations — a single
    connected component covering a large fraction of tissue
    (`suspicious_single_blob`) or an implausibly high positive fraction
    (`suspicious_high_pos`). Reuses `_load_pair` (which also validates the
    image/mask shape match).
    """
    from skimage.measure import label as sk_label

    reports = []
    for s in samples:
        _, mask = _load_pair(s)
        pos = int(mask.sum())
        pos_frac = pos / mask.size
        n_cc = int(sk_label(mask, connectivity=2).max())
        if pos == 0:
            status = "empty"
        elif n_cc == 1 and pos_frac > single_blob_area_frac:
            status = "suspicious_single_blob"
        elif pos_frac > max_pos_frac_warn:
            status = "suspicious_high_pos"
        else:
            status = "ok"
        reports.append({"name": s.name, "pos_frac": pos_frac,
                        "n_components": n_cc, "status": status})
    return reports
