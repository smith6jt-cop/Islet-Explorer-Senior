"""Unit tests for the H&E follicle U-Net utilities.

These run on CPU, use no real H&E data, and skip anything that requires
downloading ImageNet weights. They cover the bits of HE_unet_utils.py that
are easy to get wrong: post-processing geometry, donor-level k-fold, the
patch dataset, loss/metric numerics, and a tiny end-to-end forward pass.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

import HE_unet_utils as U  # noqa: E402


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def test_cosine_window_properties():
    w = U._cosine_window(64)
    assert w.shape == (64, 64)
    assert w.dtype == np.float32
    # Cosine² peaks at the centre, zero at the edges.
    # Peak near centre (cos² peaks at 0; linspace doesn't include exactly 0
    # for even N, so accept anything within ~1% of unity).
    assert w[32, 32] == pytest.approx(1.0, abs=5e-3)
    assert w[0, 0] < 1e-3
    assert w[-1, -1] < 1e-3
    # Strictly non-negative.
    assert (w >= 0).all()


def test_polygon_perimeter_unit_square():
    square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    # Closed perimeter (last → first edge included by np.diff append).
    assert U._polygon_perimeter(square) == pytest.approx(4.0)


def test_polygon_perimeter_degenerate():
    assert U._polygon_perimeter([]) == 0.0
    assert U._polygon_perimeter([(0.0, 0.0)]) == 0.0


# ---------------------------------------------------------------------------
# IoU metric
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("pred,target,expected", [
    (np.ones((10, 10)), np.ones((10, 10)), 1.0),
    (np.zeros((10, 10)), np.zeros((10, 10)), 1.0),           # both empty
    (np.ones((10, 10)), np.zeros((10, 10)), 0.0),            # disjoint
])
def test_iou_score_edge_cases(pred, target, expected):
    assert U.iou_score(pred, target) == pytest.approx(expected)


def test_iou_score_partial_overlap():
    pred = np.zeros((10, 10))
    pred[:5, :] = 1
    target = np.zeros((10, 10))
    target[3:8, :] = 1
    # |∩| = 2 rows × 10 = 20, |∪| = 8 rows × 10 = 80
    assert U.iou_score(pred, target) == pytest.approx(0.25)


# ---------------------------------------------------------------------------
# Watershed split — these are the cases Copilot flagged.
# ---------------------------------------------------------------------------
def test_watershed_empty_mask():
    mask = np.zeros((50, 50), dtype=bool)
    labels = U._watershed_split(mask, px_size_um=2.0)
    assert labels.shape == mask.shape
    assert labels.max() == 0


def test_watershed_single_small_blob_falls_back():
    # A blob too small for peak_local_max to find any peaks at the 50 µm
    # min-distance setting — must still produce one labelled instance, not zero.
    mask = np.zeros((50, 50), dtype=bool)
    mask[20:24, 20:24] = True
    labels = U._watershed_split(mask, px_size_um=2.0)
    assert (labels > 0).any()
    assert labels.max() == 1


def test_watershed_two_touching_blobs_split():
    # Two circles touching at the equator — watershed should split them
    # because their peaks are >50 µm apart.
    mask = np.zeros((200, 400), dtype=bool)
    yy, xx = np.ogrid[:200, :400]
    mask |= (yy - 100) ** 2 + (xx - 100) ** 2 < 60 ** 2
    mask |= (yy - 100) ** 2 + (xx - 300) ** 2 < 60 ** 2
    labels = U._watershed_split(mask, px_size_um=2.0)
    assert labels.max() == 2


# ---------------------------------------------------------------------------
# extract_follicle_records — schema + thresholding
# ---------------------------------------------------------------------------
EXPECTED_COLUMNS = [
    "Image", "Object ID", "Object type", "Name", "Classification", "Parent",
    "ROI", "Centroid X µm", "Centroid Y µm", "Area µm^2", "Perimeter µm",
]


def _synthetic_prob(disk_centres_radii):
    prob = np.zeros((300, 300), dtype=np.float32)
    yy, xx = np.ogrid[:300, :300]
    for (cy, cx, r) in disk_centres_radii:
        prob[(yy - cy) ** 2 + (xx - cx) ** 2 < r ** 2] = 0.95
    return prob


def test_extract_follicle_records_schema():
    prob = _synthetic_prob([(80, 80, 30), (220, 220, 35)])
    records, polygons = U.extract_follicle_records(
        prob, px_size_um=2.0, image_name="dummy.tif",
    )
    assert len(records) == 2
    assert len(polygons) == 2
    for rec in records:
        assert list(rec.keys()) == EXPECTED_COLUMNS
        assert rec["Image"] == "dummy.tif"
        assert rec["Classification"] == "Follicle"
        assert rec["Object type"] == "Annotation"
        assert rec["Area µm^2"] > 0
        assert rec["Perimeter µm"] > 0


def test_extract_follicle_records_min_area_filter():
    # 12 px radius @ 2 µm/px = 24 µm radius → area ≈ π·24² ≈ 1810 µm²,
    # well under the 10 000 µm² floor, so the record must be filtered out.
    prob = _synthetic_prob([(150, 150, 12)])
    records, _ = U.extract_follicle_records(
        prob, px_size_um=2.0, min_area_um2=10_000, image_name="dummy.tif",
    )
    assert records == []


def test_extract_follicle_records_polygons_in_um():
    prob = _synthetic_prob([(150, 150, 40)])
    records, polygons = U.extract_follicle_records(
        prob, px_size_um=2.0, image_name="dummy.tif",
    )
    assert len(records) == 1
    # Polygon vertices should be in µm and bracket the centroid.
    xs = [p[0] for p in polygons[0]]
    ys = [p[1] for p in polygons[0]]
    cx = records[0]["Centroid X µm"]
    cy = records[0]["Centroid Y µm"]
    assert min(xs) <= cx <= max(xs)
    assert min(ys) <= cy <= max(ys)


# ---------------------------------------------------------------------------
# GeoJSON writer
# ---------------------------------------------------------------------------
def test_write_geojson_roundtrip(tmp_path):
    polygons = [[(0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0)]]
    out = tmp_path / "follicles.geojson"
    U.write_geojson(polygons, out, image_name="dummy.tif")
    payload = json.loads(out.read_text())
    assert payload["type"] == "FeatureCollection"
    assert len(payload["features"]) == 1
    feat = payload["features"][0]
    assert feat["geometry"]["type"] == "Polygon"
    assert feat["properties"]["classification"]["name"] == "Follicle"
    # Ring is closed.
    ring = feat["geometry"]["coordinates"][0]
    assert ring[0] == ring[-1]


# ---------------------------------------------------------------------------
# Donor k-fold — every donor seen exactly once across val folds
# ---------------------------------------------------------------------------
def _fake_samples(n: int) -> list[U.DonorSample]:
    return [U.DonorSample(
        name=f"HDL{i:03d}_HxE",
        image_path=Path(f"/tmp/im_{i}.png"),
        mask_path=Path(f"/tmp/mk_{i}.png"),
        meta_path=Path(f"/tmp/me_{i}.json"),
        sample_id=f"SYNTH_{i:03d}",
    ) for i in range(n)]


def test_donor_kfold_partitions_uniquely():
    samples = _fake_samples(20)
    folds = U.donor_kfold(samples, k=5, genotype_xlsx=None, seed=42)
    assert len(folds) == 5
    val_names: list[str] = []
    for train, val in folds:
        names = {s.name for s in train} | {s.name for s in val}
        assert names == {s.name for s in samples}
        assert {s.name for s in train}.isdisjoint({s.name for s in val})
        val_names.extend(s.name for s in val)
    assert sorted(val_names) == sorted(s.name for s in samples)


def test_donor_kfold_single_fold_for_smoke_runs():
    samples = _fake_samples(5)
    folds = U.donor_kfold(samples, k=1, genotype_xlsx=None, seed=0)
    assert len(folds) == 1
    train, val = folds[0]
    # k=1 puts everything in val and nothing in train; the trainer handles
    # the small-n smoke path via --folds 1 (random 80/20).
    assert len(train) + len(val) == len(samples)


# ---------------------------------------------------------------------------
# Loss + dataset + a tiny end-to-end forward pass with the fallback U-Net
# ---------------------------------------------------------------------------
def test_bce_dice_loss_is_finite_and_differentiable():
    logits = torch.zeros(2, 1, 32, 32, requires_grad=True)
    target = torch.zeros(2, 1, 32, 32)
    target[:, :, 8:24, 8:24] = 1.0
    loss = U.bce_dice_loss(logits, target)
    assert torch.isfinite(loss)
    loss.backward()
    assert logits.grad is not None
    assert torch.isfinite(logits.grad).all()


def _make_synthetic_donor(tmp_path: Path, name: str = "HDL999_HxE",
                          size: int = 256):
    """Write a (image, mask, meta) triple under tmp_path.

    Thin wrapper over the shared generator in `scripts/make_synthetic_he_data.py`
    (single source of truth, also used by CI) that returns a DonorSample.
    """
    from make_synthetic_he_data import make_synthetic_donor
    make_synthetic_donor(tmp_path, name=name, size=size, n_follicles=2,
                         radius=min(30, size // 6), seed=0)
    return U.DonorSample(
        name=name,
        image_path=tmp_path / "images" / f"{name}.png",
        mask_path=tmp_path / "masks" / f"{name}.png",
        meta_path=tmp_path / "meta" / f"{name}.json",
        sample_id=name,
    )


def test_patch_dataset_tile_shapes(tmp_path):
    sample = _make_synthetic_donor(tmp_path, size=256)
    ds = U.FolliclePatchDataset([sample], patch=128, length=4, augment=True)
    for i in range(len(ds)):
        img, mask = ds[i]
        assert img.shape == (1, 128, 128)
        assert mask.shape == (1, 128, 128)
        assert img.dtype == torch.float32
        assert mask.dtype == torch.float32
        assert 0.0 <= float(img.min()) and float(img.max()) <= 1.0
        assert set(np.unique(mask.numpy())).issubset({0.0, 1.0})


def test_fallback_unet_forward_pass():
    # Force the fallback by removing smp from sys.modules — keeps CI fast and
    # avoids downloading ImageNet weights.
    import builtins
    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.startswith("segmentation_models_pytorch"):
            raise ImportError("blocked for test")
        return real_import(name, *args, **kwargs)

    builtins.__import__ = blocked
    try:
        model = U.build_unet()
    finally:
        builtins.__import__ = real_import
    assert isinstance(model, U._FallbackUNet)
    x = torch.randn(1, 1, 64, 64)
    y = model(x)
    assert y.shape == (1, 1, 64, 64)
    assert torch.isfinite(y).all()


def test_predict_wsi_runs_on_fallback_model():
    model = U._FallbackUNet().eval()
    himg = (np.random.default_rng(1).integers(0, 255,
            size=(96, 128))).astype(np.uint8)
    prob = U.predict_wsi(model, himg, device=torch.device("cpu"),
                         patch=64, overlap=16, tta=False)
    assert prob.shape == himg.shape
    assert prob.dtype == np.float32
    assert (prob >= 0).all() and (prob <= 1).all()


# ---------------------------------------------------------------------------
# Reproducibility — seeded dataset + DataLoader generator
# ---------------------------------------------------------------------------
def test_seeded_training_inputs_are_reproducible(tmp_path):
    from torch.utils.data import DataLoader
    from make_synthetic_he_data import make_synthetic_donor
    for i in range(2):
        make_synthetic_donor(tmp_path, name=f"HDL90{i}_synthetic", size=256,
                             n_follicles=2, radius=30, seed=i)
    samples = U.discover_samples(tmp_path)
    assert len(samples) == 2

    def first_batches(seed, n=3):
        U.set_seed(seed)
        ds = U.FolliclePatchDataset(samples, patch=128, length=12,
                                    augment=True, seed=seed)
        g = torch.Generator()
        g.manual_seed(seed)
        loader = DataLoader(ds, batch_size=2, shuffle=True, num_workers=0,
                            generator=g)
        it = iter(loader)
        return [next(it)[0] for _ in range(n)]

    a, b, c = first_batches(0), first_batches(0), first_batches(1)
    assert all(torch.equal(x, y) for x, y in zip(a, b))      # same seed → identical
    assert not all(torch.equal(x, y) for x, y in zip(a, c))  # different seed → differs


# ---------------------------------------------------------------------------
# Data-integrity preflight
# ---------------------------------------------------------------------------
def test_check_mask_integrity_statuses(tmp_path):
    from PIL import Image
    from make_synthetic_he_data import make_synthetic_donor

    ok_dir = tmp_path / "ok"
    make_synthetic_donor(ok_dir, name="HDL900_synthetic", size=512,
                         n_follicles=3, radius=45, seed=0)
    ok = U.check_mask_integrity(U.discover_samples(ok_dir))[0]
    assert ok["status"] == "ok"
    assert ok["n_components"] >= 2

    empty_dir = tmp_path / "empty"
    make_synthetic_donor(empty_dir, name="HDL901_synthetic", size=256,
                         n_follicles=2, radius=30, seed=0)
    Image.fromarray(np.zeros((256, 256), dtype=np.uint8)).save(
        empty_dir / "masks" / "HDL901_synthetic.png")
    assert U.check_mask_integrity(U.discover_samples(empty_dir))[0]["status"] == "empty"

    blob_dir = tmp_path / "blob"
    make_synthetic_donor(blob_dir, name="HDL902_synthetic", size=256,
                         n_follicles=2, radius=30, seed=0)
    big = np.zeros((256, 256), dtype=np.uint8)
    big[20:236, 20:236] = 255  # one component covering >50% of the image
    Image.fromarray(big).save(blob_dir / "masks" / "HDL902_synthetic.png")
    rep = U.check_mask_integrity(U.discover_samples(blob_dir))[0]
    assert rep["status"] == "suspicious_single_blob"


# ---------------------------------------------------------------------------
# Deployment-aligned instance metrics
# ---------------------------------------------------------------------------
def test_instance_metrics_counts_and_handles_empty_pred():
    size = 400
    prob = np.zeros((size, size), dtype=np.float32)
    gt = np.zeros((size, size), dtype=np.uint8)
    yy, xx = np.ogrid[:size, :size]
    for cy, cx in [(100, 100), (300, 300)]:
        disk = (yy - cy) ** 2 + (xx - cx) ** 2 < 60 ** 2
        prob[disk] = 0.95
        gt[disk] = 255
    m = U.instance_metrics(prob, gt, px_size_um=2.0)
    assert m["gt_count"] == 2
    assert m["pred_count"] == 2
    assert m["gt_mean_area_um2"] > 10_000
    assert m["pred_mean_area_um2"] > 10_000

    # All-zero prediction must not crash; counts/areas collapse to zero.
    m0 = U.instance_metrics(np.zeros((size, size), np.float32), gt, px_size_um=2.0)
    assert m0["pred_count"] == 0
    assert m0["pred_mean_area_um2"] == 0.0
    assert m0["gt_count"] == 2


# ---------------------------------------------------------------------------
# Synthetic-data generator (shared by tests + CI)
# ---------------------------------------------------------------------------
def test_make_synthetic_donor_roundtrip(tmp_path):
    from make_synthetic_he_data import make_synthetic_donor
    name, n = make_synthetic_donor(tmp_path, name="HDL900_synthetic", size=384,
                                   n_follicles=3, radius=45, seed=0)
    assert n >= 1
    samples = U.discover_samples(tmp_path)
    assert len(samples) == 1
    s = samples[0]
    assert s.name == "HDL900_synthetic"
    assert s.meta["target_px_size_um"] == 2.0
    from skimage.measure import label as sk_label
    _, mask = U._load_pair(s)
    assert int(sk_label(mask, connectivity=2).max()) == n
