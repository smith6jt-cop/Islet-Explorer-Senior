"""Unit tests for the DAPI/Phenocycler follicle U-Net utilities.

Companion to `test_HE_unet_utils.py`. The emphasis here is on the two things the
DAPI pipeline has that the H&E one does not:

  * a **Groovy/Python parity contract** for the intensity transform, which the
    H&E pipeline only documents as a manual snippet. `test_otsu_*`,
    `test_percentile_*` and `test_groovy_reference_*` pin down the exact
    arithmetic the export mirrors.
  * a **multi-head model** (binary Follicle + auxiliary region softmax) whose
    auxiliary loss must vanish for donors that have no region labels — the
    mechanism the pooled DAPI+H&E run depends on.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
import torch
from PIL import Image

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

import DAPI_unet_utils as D  # noqa: E402
import HE_unet_utils as U  # noqa: E402
from make_synthetic_dapi_data import make_synthetic_dapi_donor  # noqa: E402


# ---------------------------------------------------------------------------
# Parity contract: histogram / Otsu / percentile / normalise
# ---------------------------------------------------------------------------
def test_hist256_bins_by_high_byte():
    raw = np.array([[0, 255, 256, 511, 65535]], dtype=np.uint16)
    h = D.hist256(raw)
    assert h.sum() == 5
    assert h[0] == 2       # 0 and 255 both land in bin 0
    assert h[1] == 2       # 256 and 511 in bin 1
    assert h[255] == 1


def test_otsu_separates_two_clear_modes():
    hist = np.zeros(256, dtype=np.int64)
    hist[10] = 1000
    hist[200] = 1000
    t = D.otsu_from_hist256(hist)
    assert 10 <= t < 200, f"threshold {t} should fall between the two modes"


def test_otsu_degenerate_histograms_return_zero():
    assert D.otsu_from_hist256(np.zeros(256, dtype=np.int64)) == 0
    single = np.zeros(256, dtype=np.int64)
    single[42] = 500
    # One populated bin has no valid two-class split; 0 makes everything tissue,
    # which is the safe direction.
    assert D.otsu_from_hist256(single) == 0


def test_otsu_returns_first_maximum():
    """numpy argmax takes the FIRST max; the Groovy loop uses strict '>' to match.

    A symmetric histogram is the case where 'first' vs 'last' actually differ, so
    it is what pins the two implementations together.
    """
    hist = np.zeros(256, dtype=np.int64)
    hist[100] = 500
    hist[101] = 500
    t = D.otsu_from_hist256(hist)
    manual = _reference_otsu(hist)
    assert t == manual


def test_otsu3_separates_glass_from_tissue_not_tissue_from_tissue():
    """The reason 2-class Otsu is not enough for a spleen DAPI histogram.

    Three modes — glass, red pulp, white pulp. 2-class Otsu puts its threshold
    between red pulp and white pulp, which would make "tissue" mean white pulp
    only. The 3-class lower threshold has to land between glass and red pulp.
    """
    hist = np.zeros(256, dtype=np.int64)
    hist[1] = 30_000       # glass
    hist[45] = 62_000      # red pulp
    hist[195] = 7_000      # white pulp

    t1, t2 = D.otsu3_from_hist256(hist)
    assert 1 <= t1 < 45, f"3-class lower threshold {t1} should split glass/tissue"
    assert 45 <= t2 < 195, f"3-class upper threshold {t2} should split the pulps"

    # Tissue is `bin > t`. Under 2-class Otsu the red-pulp mode at bin 45 falls
    # on the background side, so the mask would contain white pulp only — the
    # exact failure 3-class Otsu exists to avoid. Asserted so this test still
    # means something if the 2-class implementation ever changes.
    t2class = D.otsu_from_hist256(hist)
    assert t2class >= 45, ("this histogram is the case 2-class Otsu gets wrong; "
                           "if it no longer does, the test lost its point")


def test_otsu3_degenerate_histograms():
    assert D.otsu3_from_hist256(np.zeros(256, dtype=np.int64)) == (0, 0)
    single = np.zeros(256, dtype=np.int64)
    single[7] = 100
    t1, t2 = D.otsu3_from_hist256(single)      # falls back to the 2-class value
    assert t1 == t2 == D.otsu_from_hist256(single)


def test_tissue_mask_excludes_the_glass_bin():
    raw = np.array([[100, 300, 20000, 40000]], dtype=np.uint16)
    # bin 0 is the threshold's background side; tissue must start above it.
    assert list(D.tissue_from_raw(raw, otsu_bin=0)[0]) == [False, True, True, True]
    assert list(D.tissue_from_raw(raw, otsu_bin=1)[0]) == [False, False, True, True]


def test_percentile_from_hist_is_bin_low_edge():
    hist = np.zeros(256, dtype=np.int64)
    hist[4] = 50
    hist[8] = 50
    assert D.percentile_from_hist256(hist, 1.0) == 4 << 8
    assert D.percentile_from_hist256(hist, 99.5) == 8 << 8
    assert D.percentile_from_hist256(np.zeros(256, dtype=np.int64), 50.0) == 0


def test_dapi_to_uint8_inverts_and_clips():
    raw = np.array([[0, 1000, 2000, 3000, 65535]], dtype=np.uint16)
    out = D.dapi_to_uint8(raw, lo=1000, hi=3000)
    assert out.dtype == np.uint8
    assert out[0, 0] == 255 and out[0, 1] == 255      # at/below lo -> white
    assert out[0, 3] == 0 and out[0, 4] == 0          # at/above hi -> black
    assert 100 < int(out[0, 2]) < 155                 # midpoint -> mid grey
    # Monotonically decreasing: high DAPI must be dark, matching rgb_to_h.
    assert list(out[0]) == sorted(out[0], reverse=True)


def test_dapi_to_uint8_rejects_degenerate_window():
    with pytest.raises(ValueError):
        D.dapi_to_uint8(np.zeros((4, 4), dtype=np.uint16), lo=100, hi=100)


def test_polarity_matches_hematoxylin():
    """Dense nuclei must be DARK in both modalities, or pooling is meaningless."""
    dark_he = U.rgb_to_h(np.array([[[40, 30, 90]]], dtype=np.uint8))    # H-dense
    light_he = U.rgb_to_h(np.array([[[245, 245, 245]]], dtype=np.uint8))  # white
    assert dark_he[0, 0] < light_he[0, 0]

    dapi = D.dapi_to_uint8(np.array([[60000, 500]], dtype=np.uint16), 1000, 50000)
    assert dapi[0, 0] < dapi[0, 1]     # bright DAPI -> dark, dim DAPI -> light


def _reference_otsu(hist: np.ndarray) -> int:
    """Straight transcription of the Groovy `otsuFromHist` loop.

    Deliberately written as the imperative loop the Groovy uses rather than
    reusing the vectorised Python, so agreement is evidence and not a tautology.
    """
    total = float(hist.sum())
    sum_all = float(sum(float(hist[i]) * i for i in range(256)))
    if total <= 0:
        return 0
    w_b = 0.0
    sum_b = 0.0
    best = -1.0
    best_t = 0
    for t in range(256):
        w_b += hist[t]
        sum_b += float(hist[t]) * t
        w_f = total - w_b
        if w_b <= 0 or w_f <= 0:
            continue
        mean_b = sum_b / w_b
        mean_f = (sum_all - sum_b) / w_f
        between = w_b * w_f * (mean_b - mean_f) ** 2
        if between > best:
            best = between
            best_t = t
    return best_t


def _reference_otsu3(hist: np.ndarray) -> tuple[int, int]:
    """Straight transcription of the Groovy `otsu3FromHist` loop."""
    total = 0.0
    W = [0.0] * 257
    S = [0.0] * 257
    for i in range(256):
        total += float(hist[i])
        W[i + 1] = W[i] + float(hist[i])
        S[i + 1] = S[i] + float(hist[i]) * i
    if total <= 0.0:
        return 0, 0
    best = -1.0
    bt1 = bt2 = 0
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
                best, bt1, bt2 = v, t1, t2
    if best < 0.0:
        t = _reference_otsu(hist)
        return t, t
    return bt1, bt2


def _reference_percentile(hist: np.ndarray, pct: float) -> int:
    """Transcription of the Groovy `percentileFromHist` loop."""
    total = int(hist.sum())
    if total <= 0:
        return 0
    target = pct / 100.0 * total
    cum = 0
    for b in range(256):
        cum += int(hist[b])
        if cum >= target:
            return b << 8
    return 255 << 8


def _reference_normalise(raw: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """Transcription of the Groovy per-pixel normalise+invert."""
    out = np.empty(raw.shape, dtype=np.uint8)
    span = float(hi - lo)
    flat = raw.ravel()
    res = np.empty(flat.size, dtype=np.uint8)
    for i, v in enumerate(flat.tolist()):
        x = (v - lo) / span
        x = 0.0 if x < 0.0 else (1.0 if x > 1.0 else x)
        res[i] = int(math.floor(255.0 * (1.0 - x)))
    return res.reshape(out.shape)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_groovy_reference_transform_matches_vectorised(seed):
    """The whole chain, imperative reference vs the shipped vectorised code.

    This is the automated form of the parity check the H&E README only describes
    in prose. It cannot run Groovy, but it does pin the exact arithmetic the
    Groovy is a transcription of, so a drift in either direction fails here.
    """
    rng = np.random.default_rng(seed)
    raw = np.concatenate([
        rng.integers(200, 3000, size=4000),
        rng.integers(20000, 60000, size=1500),
    ]).astype(np.uint16).reshape(110, 50)

    h = D.hist256(raw)
    assert D.otsu_from_hist256(h) == _reference_otsu(h)
    assert D.otsu3_from_hist256(h) == _reference_otsu3(h)

    otsu = D.otsu3_from_hist256(h)[0]
    tissue_hist = h.copy()
    tissue_hist[:otsu + 1] = 0     # tissue is bin > otsu (Otsu's t is background)
    for pct in (1.0, 99.5):
        assert D.percentile_from_hist256(tissue_hist, pct) == \
            _reference_percentile(tissue_hist, pct)

    params = D.dapi_transform_params(raw)
    assert params["otsu_bin"] == otsu
    assert params["lo"] == _reference_percentile(tissue_hist, 1.0)
    assert params["hi"] == _reference_percentile(tissue_hist, 99.5)

    got = D.dapi_to_uint8(raw, params["lo"], params["hi"])
    assert np.array_equal(got, _reference_normalise(raw, params["lo"], params["hi"]))


def test_transform_params_use_tissue_not_whole_image():
    """Percentiles are taken within the Otsu tissue mask, not over everything.

    A slide that is mostly empty background would otherwise put `lo` down in the
    background and crush the tissue into a handful of grey levels.
    """
    raw = np.full((200, 200), 300, dtype=np.uint16)          # 97.5% background
    raw[:32, :] = np.random.default_rng(0).integers(
        30000, 50000, size=(32, 200)).astype(np.uint16)       # tissue strip
    params = D.dapi_transform_params(raw)
    assert params["lo"] > 10000, "lo was pulled into the background"


# ---------------------------------------------------------------------------
# Region label space
# ---------------------------------------------------------------------------
def test_follicle_is_painted_last():
    """Follicle must be the highest index so the exporter's paint order makes
    `regions == FOLLICLE_INDEX` agree with the binary mask."""
    assert D.REGION_CLASSES[-1] == "Follicle"
    assert D.FOLLICLE_INDEX == len(D.REGION_CLASSES) - 1
    assert D.PALS_INDEX == D.FOLLICLE_INDEX - 1


def test_merge_white_pulp_maps_into_the_merged_space():
    reg = np.array([[0, 1, 2, 3, D.PALS_INDEX, D.FOLLICLE_INDEX]], dtype=np.uint8)
    merged = D.merge_white_pulp(reg)
    assert merged.max() < len(D.MERGED_REGION_CLASSES)
    assert D.MERGED_REGION_CLASSES[merged[0, -1]] == "WhitePulp"
    assert D.MERGED_REGION_CLASSES[merged[0, -2]] == "WhitePulp"
    assert merged[0, 1] == 1     # RedPulp untouched


def test_region_confusion_absent_class_is_none_not_one():
    """A class in neither prediction nor truth must not score a free 1.0.

    `iou_score` keeps its both-empty==1.0 convention for the binary head; scoring
    an absent region class the same way would inflate the reported mean by a
    class the model was never asked to find.
    """
    gt = np.zeros((8, 8), dtype=np.uint8)
    gt[:4] = D.FOLLICLE_INDEX
    pred = gt.copy()
    r = D.region_confusion(pred, gt)
    assert r["per_class_iou"]["Follicle"] == pytest.approx(1.0)
    assert r["per_class_iou"]["Trabeculae"] is None
    assert r["n_classes_scored"] == 2         # Background + Follicle
    assert r["mean_iou"] == pytest.approx(1.0)


def test_region_confusion_detects_follicle_pals_swap():
    """The failure mode this cohort is most prone to has to be visible."""
    gt = np.full((10, 10), D.FOLLICLE_INDEX, dtype=np.uint8)
    pred = np.full((10, 10), D.PALS_INDEX, dtype=np.uint8)
    r = D.region_confusion(pred, gt)
    assert r["follicle_as_pals"] == pytest.approx(1.0)
    assert r["per_class_iou"]["Follicle"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Donor identity / leakage
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name,expected", [
    ("HDL011_PC33", "HDL011"),
    ("HDL073_PC29", "HDL073"),
    ("1901HBMP004_PC29", "1901HBMP004"),
    ("HBMP006", "HBMP006"),
])
def test_donor_id_from_name(name, expected):
    assert U.donor_id_from_name(name) == expected


def test_donor_id_from_name_rejects_unlisted_alias():
    with pytest.raises(U.SampleManifestError, match="unlisted sample alias"):
        U.donor_id_from_name("HDL900_synthetic")


def test_multi_file_donor_stays_in_one_fold():
    """HDL073_PC19 and HDL073_PC29 are one donor; splitting them leaks patches.

    This is the concrete bug the widened regex plus donor-grouping fixes — under
    the old file-level dealing these two landed in different folds.
    """
    names = ["HDL073_PC19", "HDL073_PC29", "HDL011_PC33", "HDL043_PC33",
             "HDL055_PC33", "1901HBMP004_PC29"]
    identities = {
        "HDL073_PC19": "HDL073", "HDL073_PC29": "HDL073",
        "HDL011_PC33": "HDL011", "HDL043_PC33": "HDL043",
        "HDL055_PC33": "HDL055", "1901HBMP004_PC29": "1901HBMP004",
    }
    samples = [U.DonorSample(n, Path(f"{n}.png"), Path(f"{n}.png"),
                             Path(f"{n}.json"), sample_id=identities[n])
               for n in names]
    for k in (2, 3):
        for train, val in U.donor_kfold(samples, k=k, seed=0):
            tr = {s.name for s in train}
            va = {s.name for s in val}
            assert not (tr & va)
            assert ("HDL073_PC19" in tr) == ("HDL073_PC29" in tr)
            assert ("HDL073_PC19" in va) == ("HDL073_PC29" in va)


# ---------------------------------------------------------------------------
# Model, loss, inference
# ---------------------------------------------------------------------------
def test_multitask_model_head_width():
    model = D.build_multitask_unet("efficientnet-b0", pretrained=False,
                                   n_region_classes=len(D.REGION_CLASSES))
    out = model(torch.randn(1, 1, 64, 64))
    assert out.shape == (1, 1 + len(D.REGION_CLASSES), 64, 64)


def test_no_aux_model_is_shape_identical_to_the_he_model():
    """K=0 must give a model an H&E checkpoint can be warm-started from."""
    dapi = D.build_multitask_unet("efficientnet-b0", pretrained=False,
                                  n_region_classes=0)
    he = U.build_unet("efficientnet-b0", pretrained=False)
    assert {k: v.shape for k, v in dapi.state_dict().items()} == \
           {k: v.shape for k, v in he.state_dict().items()}


def test_multitask_activation_normalises_region_channels():
    logits = torch.randn(2, 1 + 4, 8, 8)
    p = D.multitask_activation(logits)
    assert torch.all((p[:, 0] >= 0) & (p[:, 0] <= 1))
    assert torch.allclose(p[:, 1:].sum(dim=1), torch.ones(2, 8, 8), atol=1e-5)
    # K=0 must still work — that is the H&E-shaped path.
    assert D.multitask_activation(torch.randn(1, 1, 4, 4)).shape == (1, 1, 4, 4)


def test_multitask_loss_is_finite_and_differentiable():
    logits = torch.randn(2, 1 + len(D.REGION_CLASSES), 32, 32, requires_grad=True)
    mask = (torch.rand(2, 1, 32, 32) > 0.5).float()
    region = torch.randint(0, len(D.REGION_CLASSES), (2, 32, 32))
    valid = torch.ones(2)
    loss = D.multitask_loss(logits, mask, region, valid, lam=0.3)
    assert torch.isfinite(loss)
    loss.backward()
    assert torch.isfinite(logits.grad).all()


def test_aux_loss_vanishes_when_no_region_labels():
    """The pooling mechanism: H&E donors must contribute Follicle gradient only.

    With aux_valid=0 the loss has to equal the pure Follicle loss exactly, no
    matter what garbage sits in the region target.
    """
    torch.manual_seed(0)
    logits = torch.randn(3, 1 + len(D.REGION_CLASSES), 16, 16)
    mask = (torch.rand(3, 1, 16, 16) > 0.5).float()
    region = torch.randint(0, len(D.REGION_CLASSES), (3, 16, 16))
    primary = U.bce_dice_loss(logits[:, :1], mask)

    off = D.multitask_loss(logits, mask, region, torch.zeros(3), lam=0.3)
    assert torch.allclose(off, primary, atol=1e-6)

    on = D.multitask_loss(logits, mask, region, torch.ones(3), lam=0.3)
    assert on > primary

    # A mixed batch sits strictly between the two.
    mixed = D.multitask_loss(logits, mask, region,
                             torch.tensor([1.0, 0.0, 0.0]), lam=0.3)
    assert primary < mixed < on


def test_region_loss_ignores_classes_absent_from_the_patch():
    """A patch containing only Background+Follicle must not be penalised for the
    four regions that are not in frame."""
    torch.manual_seed(0)
    n = len(D.REGION_CLASSES)
    logits = torch.randn(1, n, 16, 16)
    target = torch.full((1, 16, 16), D.FOLLICLE_INDEX, dtype=torch.long)
    perfect = torch.full((1, n, 16, 16), -10.0)
    perfect[:, D.FOLLICLE_INDEX] = 10.0
    loss = D.region_loss(perfect, target)
    assert float(loss) < 0.05, "absent classes are leaking a fixed dice penalty"
    assert torch.isfinite(D.region_loss(logits, target)).all()


def test_region_class_weights_zero_for_absent_classes(tmp_path):
    reg = np.zeros((32, 32), dtype=np.uint8)
    reg[:16] = D.FOLLICLE_INDEX
    p = tmp_path / "d.png"
    Image.fromarray(reg).save(p)
    w = D.region_class_weights({"d": p}, len(D.REGION_CLASSES))
    assert w[D.REGION_CLASSES.index("Trabeculae")] == 0.0
    assert w[D.FOLLICLE_INDEX] > 0
    assert w[0] > 0


def test_follicle_activation_matches_channel_zero_of_the_full_map():
    """Stitching only channel 0 must give the same answer as stitching all of
    them and slicing — that equivalence is what makes the memory saving free."""
    model = D.build_multitask_unet("efficientnet-b0", pretrained=False,
                                   n_region_classes=3)
    img = np.random.default_rng(1).integers(0, 256, size=(80, 96)).astype(np.uint8)
    dev = torch.device("cpu")
    full = U.predict_wsi(model, img, dev, patch=64, overlap=16, tta=False,
                         n_out=4, activation=D.multitask_activation)
    only = U.predict_wsi(model, img, dev, patch=64, overlap=16, tta=False,
                         n_out=1, activation=D.follicle_activation)
    assert only.shape == (80, 96)
    assert np.allclose(full[0], only, atol=1e-6)


def test_region_confusion_banded_matches_single_pass():
    """Banding must not change the metrics: confusion matrices add."""
    model = D.build_multitask_unet("efficientnet-b0", pretrained=False,
                                   n_region_classes=3)
    names = ["Background", "A", "B"]
    img = np.random.default_rng(2).integers(0, 256, size=(200, 96)).astype(np.uint8)
    gt = np.random.default_rng(3).integers(0, 3, size=(200, 96)).astype(np.uint8)
    dev = torch.device("cpu")

    prob = U.predict_wsi(model, img, dev, patch=64, overlap=16, tta=False,
                         n_out=4, activation=D.multitask_activation)
    single = D.region_confusion(np.argmax(prob[1:], axis=0).astype(np.uint8),
                                gt, names)
    banded = D.region_confusion_banded(model, img, gt, dev, 3, names,
                                       patch=64, overlap=16, band_rows=200,
                                       tta=False)
    assert banded["per_class_iou"] == single["per_class_iou"]

    # A genuinely multi-band run stays close (only the seam rows blend
    # differently), and must cover every pixel exactly once.
    multi = D.region_confusion_banded(model, img, gt, dev, 3, names,
                                      patch=64, overlap=16, band_rows=96,
                                      tta=False)
    assert np.array_equal(np.array(multi["confusion"]).sum(), gt.size)


def test_predict_wsi_rejects_overlap_at_or_above_patch():
    """Otherwise this surfaces as a bare ZeroDivisionError from the padding
    arithmetic, far from the actual mistake."""
    model = U.build_unet("efficientnet-b0", pretrained=False)
    img = np.zeros((80, 80), dtype=np.uint8)
    with pytest.raises(ValueError, match="overlap < patch"):
        U.predict_wsi(model, img, torch.device("cpu"), patch=64, overlap=64)


def test_confusion_counts_are_additive():
    gt = np.array([[0, 1], [2, 1]], dtype=np.uint8)
    pred = np.array([[0, 0], [2, 1]], dtype=np.uint8)
    whole = D.confusion_counts(pred, gt, 3)
    halves = (D.confusion_counts(pred[:1], gt[:1], 3)
              + D.confusion_counts(pred[1:], gt[1:], 3))
    assert np.array_equal(whole, halves)
    assert whole.sum() == gt.size


def test_predict_wsi_multichannel_output_shape():
    """`n_out > 1` returns (C, H, W); n_out == 1 keeps the historical 2-D shape."""
    model = D.build_multitask_unet("efficientnet-b0", pretrained=False,
                                   n_region_classes=3)
    img = np.random.default_rng(0).integers(0, 256, size=(96, 128)).astype(np.uint8)
    out = U.predict_wsi(model, img, torch.device("cpu"), patch=64, overlap=16,
                        tta=False, n_out=4, activation=D.multitask_activation)
    assert out.shape == (4, 96, 128)
    assert np.isfinite(out).all()
    assert np.allclose(out[1:].sum(axis=0), 1.0, atol=1e-4)

    single = U.build_unet("efficientnet-b0", pretrained=False)
    out1 = U.predict_wsi(single, img, torch.device("cpu"), patch=64, overlap=16,
                         tta=False)
    assert out1.shape == (96, 128)


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------
def _synth(tmp_path: Path, name: str = "HDL900_PC33_synthetic", **kw):
    make_synthetic_dapi_donor(tmp_path, name, size=256, n_follicles=2,
                              radius=30, seed=0, **kw)
    samples = U.discover_samples(tmp_path)
    return samples, D.attach_regions(samples, tmp_path)


def test_synthetic_donor_roundtrip(tmp_path):
    samples, regions = _synth(tmp_path)
    assert len(samples) == 1
    s = samples[0]
    assert s.name in regions
    meta = json.loads(s.meta_path.read_text())
    assert meta["target_px_size_um"] == D.TARGET_UM_PER_PX
    assert meta["image_type"] == "dapi"
    D.sidecar_transform(meta)     # must not raise

    reg = np.array(Image.open(regions[s.name]).convert("L"))
    msk = np.array(Image.open(s.mask_path).convert("L"))
    # The exporter derives masks from regions; the generator must match.
    assert np.array_equal(msk > 127, reg == D.FOLLICLE_INDEX)


def test_dataset_returns_four_tensors(tmp_path):
    samples, regions = _synth(tmp_path)
    ds = D.DapiRegionDataset(samples, regions, patch=128, length=4, augment=True,
                             seed=0)
    img, mask, region, valid = ds[0]
    assert img.shape == (1, 128, 128) and img.dtype == torch.float32
    assert mask.shape == (1, 128, 128)
    assert region.shape == (128, 128) and region.dtype == torch.int64
    assert float(valid) == 1.0
    assert 0.0 <= float(img.min()) and float(img.max()) <= 1.0
    assert set(torch.unique(mask).tolist()) <= {0.0, 1.0}
    assert int(region.max()) < len(D.REGION_CLASSES)


def test_dataset_pads_donors_smaller_than_one_patch(tmp_path):
    """A donor whose export is smaller than the patch must still yield full
    patch-sized tiles.

    Regression test. `_sample_center` used to clamp with
    `np.clip(y, half, h - half)`; when `h < patch` that has a_min > a_max, numpy
    returns a_max, and the slice then started at a negative index and came back
    undersized. Nothing raised — the failure surfaced far away as a DataLoader
    collate error, and only when such a donor was batched next to a normal one.
    """
    samples, regions = _synth(tmp_path, name="HDL900_small")
    ds = D.DapiRegionDataset(samples, regions, patch=512, length=6,
                             augment=True, seed=0)
    for i in range(6):
        img, mask, region, _ = ds[i]
        assert img.shape == (1, 512, 512)
        assert mask.shape == (1, 512, 512)
        assert region.shape == (512, 512)


def test_pooled_mixed_size_donors_collate(tmp_path):
    """The exact shape CI hit: pooled donors of different export sizes.

    A 256 px donor and a 640 px donor in one batch must stack, which they only
    do if every tile is padded/cropped to exactly the patch size.
    """
    from torch.utils.data import DataLoader

    make_synthetic_dapi_donor(tmp_path, "HDL900_small", size=256,
                              n_follicles=2, radius=30, seed=0)
    make_synthetic_dapi_donor(tmp_path, "HDL901_big", size=640,
                              n_follicles=3, radius=40, seed=1)
    samples = U.discover_samples(tmp_path)
    assert len(samples) == 2
    ds = D.DapiRegionDataset(samples, D.attach_regions(samples, tmp_path),
                             patch=512, length=8, augment=True, seed=0)
    img, mask, region, valid = next(iter(DataLoader(ds, batch_size=4, shuffle=True)))
    assert img.shape == (4, 1, 512, 512)
    assert mask.shape == (4, 1, 512, 512)
    assert region.shape == (4, 512, 512)
    assert valid.shape == (4,)


def test_he_dataset_also_pads_small_donors(tmp_path):
    """The fix lives in the shared base class, so the H&E pipeline gets it too."""
    from make_synthetic_he_data import make_synthetic_donor

    make_synthetic_donor(tmp_path, "HDL900_synthetic", size=256, n_follicles=2,
                         radius=30, seed=0)
    samples = U.discover_samples(tmp_path)
    ds = U.FolliclePatchDataset(samples, patch=512, length=4, augment=True, seed=0)
    for i in range(4):
        img, mask = ds[i]
        assert img.shape == (1, 512, 512) and mask.shape == (1, 512, 512)


def test_dataset_without_regions_reports_aux_invalid(tmp_path):
    """The pooled-run path: a donor with no region map trains the Follicle head
    only, and must not crash or fabricate labels."""
    samples, _ = _synth(tmp_path)
    ds = D.DapiRegionDataset(samples, {}, patch=128, length=2, augment=False, seed=0)
    _, _, region, valid = ds[0]
    assert float(valid) == 0.0
    assert int(region.max()) == 0


def test_dataset_augmentation_keeps_mask_and_region_aligned(tmp_path):
    """Rotating the image but not the region map would silently destroy the aux
    labels, and the Follicle IoU would still look fine."""
    samples, regions = _synth(tmp_path)
    ds = D.DapiRegionDataset(samples, regions, patch=128, length=8, augment=True,
                             seed=3)
    for i in range(8):
        _, mask, region, _ = ds[i]
        m = mask[0].numpy() > 0.5
        r = region.numpy() == D.FOLLICLE_INDEX
        assert np.array_equal(m, r), "mask and region map desynchronised"


def test_dataset_rejects_mismatched_region_shape(tmp_path):
    samples, regions = _synth(tmp_path)
    Image.fromarray(np.zeros((64, 64), dtype=np.uint8)).save(regions[samples[0].name])
    with pytest.raises(ValueError, match="does not match"):
        D.DapiRegionDataset(samples, regions, patch=64, length=2, seed=0)


def test_check_region_integrity_statuses(tmp_path):
    samples, regions = _synth(tmp_path)
    ok = D.check_region_integrity(samples, regions)
    assert ok[0]["status"] == "ok"
    assert "Follicle" in ok[0]["present_classes"]

    assert D.check_region_integrity(samples, {})[0]["status"] == "no_regions"

    flat = np.zeros_like(np.array(Image.open(regions[samples[0].name]).convert("L")))
    Image.fromarray(flat).save(regions[samples[0].name])
    assert D.check_region_integrity(samples, regions)[0]["status"] == "single_label"


# ---------------------------------------------------------------------------
# Multiplex reader
# ---------------------------------------------------------------------------
def test_read_dapi_channel_picks_the_target_level(tmp_path):
    """The synthetic OME-TIFF is a 2-level pyramid whose coarse level holds the
    authored raster verbatim, so a correct level choice returns it byte-exact."""
    _synth(tmp_path, write_ome=True, n_channels=4)
    ome = next((tmp_path / "ome").glob("*.ome.tiff"))
    raw, px = D.read_dapi_channel(ome, D.TARGET_UM_PER_PX)
    assert px == D.TARGET_UM_PER_PX
    assert raw.shape == (256, 256)
    assert raw.dtype == np.uint16


def test_ome_roundtrip_reproduces_the_exported_png(tmp_path):
    """End-to-end parity: OME-TIFF -> read -> normalise == the training PNG.

    This is the check that would catch the export and the inference path drifting
    apart, which is otherwise invisible until model quality quietly degrades.
    """
    samples, _ = _synth(tmp_path, write_ome=True)
    s = samples[0]
    ome = next((tmp_path / "ome").glob("*.ome.tiff"))
    raw, _ = D.read_dapi_channel(ome, D.TARGET_UM_PER_PX)

    meta = json.loads(s.meta_path.read_text())
    t = D.sidecar_transform(meta)
    assert np.array_equal(
        D.dapi_to_uint8(raw, int(t["lo"]), int(t["hi"])),
        np.array(Image.open(s.image_path).convert("L"), dtype=np.uint8))

    # And re-deriving lo/hi from the raw pixels lands on the same window, which
    # is what inference does when there is no sidecar.
    p = D.dapi_transform_params(raw)
    assert (p["otsu_bin"], p["lo"], p["hi"]) == (t["otsu_bin"], t["lo"], t["hi"])


def test_channel_read_takes_one_page_not_the_whole_stack(tmp_path):
    """Reading one channel must not decode all of them.

    These files are CYX with one IFD per channel; `level.asarray()` would decode
    every channel to reach one plane — ~7 GB on a 33-channel slide at 2 µm/px.
    The single-page path has to return byte-identical pixels to the slow path.
    """
    import tifffile
    _synth(tmp_path, write_ome=True, n_channels=5)
    ome = next((tmp_path / "ome").glob("*.ome.tiff"))
    with tifffile.TiffFile(ome) as tif:
        level = tif.series[0].levels[-1]
        assert level.axes.startswith("C") and len(level.pages) == 5
        fast = D._read_channel_plane(level, 0)
        slow = np.take(level.asarray(), 0, axis=0)
    assert fast.shape == slow.shape
    assert np.array_equal(fast, slow)


def test_find_dapi_index_is_zero_for_the_synthetic_stack(tmp_path):
    _synth(tmp_path, write_ome=True, n_channels=3)
    ome = next((tmp_path / "ome").glob("*.ome.tiff"))
    assert D.find_dapi_index(ome) == 0


# ---------------------------------------------------------------------------
# GeoJSON units
# ---------------------------------------------------------------------------
def test_write_geojson_scale_converts_um_to_px(tmp_path):
    """QuPath's object import wants pixels; the H&E predictions are in µm. The
    scale factor is the only thing standing between the two."""
    poly = [(0.0, 0.0), (100.0, 0.0), (100.0, 50.0)]
    um_path, px_path = tmp_path / "um.geojson", tmp_path / "px.geojson"
    U.write_geojson([poly], um_path, image_name="x")
    U.write_geojson([poly], px_path, image_name="x", scale=1.0 / 2.0)

    um = json.loads(um_path.read_text())["features"][0]["geometry"]["coordinates"][0]
    px = json.loads(px_path.read_text())["features"][0]["geometry"]["coordinates"][0]
    assert um[1] == [100.0, 0.0]
    assert px[1] == [50.0, 0.0]
    assert len(um) == len(px) == 4          # ring closed in both


def test_write_geojson_honours_classification():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "g.geojson"
        U.write_geojson([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]], p,
                        classification="PALS")
        feat = json.loads(p.read_text())["features"][0]
        assert feat["properties"]["classification"]["name"] == "PALS"


def test_extract_records_honours_classification():
    prob = np.zeros((80, 80), dtype=np.float32)
    yy, xx = np.ogrid[:80, :80]
    prob[(yy - 40) ** 2 + (xx - 40) ** 2 < 30 ** 2] = 1.0
    recs, _ = U.extract_follicle_records(prob, px_size_um=2.0, min_area_um2=1000,
                                         image_name="im", classification="PALS")
    assert recs and all(r["Classification"] == "PALS" for r in recs)
    # Default is unchanged for every existing caller.
    recs2, _ = U.extract_follicle_records(prob, px_size_um=2.0, min_area_um2=1000)
    assert recs2[0]["Classification"] == "Follicle"
