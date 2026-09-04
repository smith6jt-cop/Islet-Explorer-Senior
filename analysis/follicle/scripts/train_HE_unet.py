"""Train the H&E follicle U-Net with donor-stratified k-fold CV.

Assumes `export_HE_training_pairs.groovy` has populated data/
with (images, masks, meta) PNG/JSON triples.

  # Full genotype-stratified 5-fold CV (the real, reportable run):
  python scripts/train_HE_unet.py --epochs 50 --folds 5

  # Single-fold quality estimate (random 80/20; train to a plateau):
  python scripts/train_HE_unet.py --folds 1 --epochs 60 --patience 12

  # Fast smoke test — synthetic data, CPU, no pretrained weights, no network:
  python scripts/make_synthetic_he_data.py --out-dir /tmp/he --n-donors 3
  python scripts/train_HE_unet.py --smoke --data-dir /tmp/he

Outputs one checkpoint per fold to <out-dir>/fold_{i}.pt and a cv_summary.json
that records, per fold, the best validation IoU, the FULL per-epoch history
(loss / val IoU / lr), and best-checkpoint follicle count + mean area vs. ground
truth. With >=3 validation donors it also writes pooled count/area Spearman
(the H18 success-bar metric) to pooled_metrics.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch
from PIL import Image
from torch.utils.data import DataLoader
from tqdm import tqdm

import HE_unet_utils as U

HEARTBEAT_EVERY = 50

# Args that --smoke overrides. argparse leaves them None when the user doesn't
# pass them, so resolve_defaults() can fill the mode-appropriate value while
# still honouring anything the user set explicitly.
NORMAL_DEFAULTS = dict(epochs=50, folds=5, batch_size=16, patches_per_epoch=2000,
                       patience=0, pretrained=True)
SMOKE_DEFAULTS = dict(epochs=2, folds=1, batch_size=2, patches_per_epoch=8,
                      patience=1, pretrained=False)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    # Mode-resolved args default to None (sentinel = "user did not set this").
    p.add_argument("--epochs", type=int, default=None, help="(default 50; 2 under --smoke)")
    p.add_argument("--folds", type=int, default=None, help="(default 5; 1 under --smoke)")
    p.add_argument("--batch-size", type=int, default=None, help="(default 16; 2 under --smoke)")
    p.add_argument("--patches-per-epoch", type=int, default=None,
                   help="(default 2000; 8 under --smoke)")
    p.add_argument("--patience", type=int, default=None,
                   help="epochs without val-IoU improvement before stopping; "
                        "0 disables (default 0; 1 under --smoke)")
    p.add_argument("--pretrained", action=argparse.BooleanOptionalAction, default=None,
                   help="ImageNet-pretrained encoder (default True; False under --smoke)")
    p.add_argument("--lr", type=float, default=1e-4)
    p.add_argument("--lr-tmax", type=int, default=None,
                   help="cosine-annealing period in epochs (default: --epochs). "
                        "Decoupled from --epochs so early stopping doesn't depend "
                        "on the schedule bottoming out.")
    p.add_argument("--encoder", default="efficientnet-b3")
    p.add_argument("--genotype-xlsx",
                   default=str(U.CODE_ROOT / "Spleen_rs3184504_Genotypes.xlsx"))
    p.add_argument("--data-dir", default=None,
                   help="directory holding {images,masks,meta} (default: repo data/)")
    p.add_argument("--extra-data-dir", action="append", default=[],
                   help="additional data dir(s) to pool in, e.g. the DAPI "
                        "`data_dapi/`. Repeatable. Both cohorts write "
                        "single-channel images with high nuclear density = dark, "
                        "so they can share one Follicle head.")
    p.add_argument("--init-from", default=None,
                   help="checkpoint to warm-start from; only shape-matching "
                        "weights are loaded (e.g. a DAPI fold_0.pt, whose head "
                        "is wider and is therefore skipped)")
    p.add_argument("--out-dir", default=None,
                   help="where to write fold_{i}.pt + cv_summary.json "
                        "(default: analysis/HE_unet_models/; <data-dir>/HE_unet_models "
                        "under --smoke so the git-tracked dir stays clean)")
    p.add_argument("--smoke", action="store_true",
                   help="fast crash-check: tiny epochs/patches/batch, CPU, no "
                        "pretrained weights")
    # 0 by default: the dataset caches whole-donor RGB/mask arrays, so each
    # forked worker would duplicate that cache in RAM. Bump only if you've
    # confirmed memory headroom.
    p.add_argument("--num-workers", type=int, default=0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--device", default=None,
                   help="(default cuda if available else cpu; cpu under --smoke)")
    return resolve_defaults(p.parse_args())


def resolve_defaults(args: argparse.Namespace) -> argparse.Namespace:
    """Fill mode-appropriate values for any arg the user left unset (None)."""
    base = SMOKE_DEFAULTS if args.smoke else NORMAL_DEFAULTS
    for k, v in base.items():
        if getattr(args, k) is None:
            setattr(args, k, v)
    if args.device is None:
        args.device = ("cpu" if args.smoke
                       else ("cuda" if torch.cuda.is_available() else "cpu"))
    if args.lr_tmax is None:
        args.lr_tmax = args.epochs
    if args.out_dir is None:
        base_dir = Path(args.data_dir) if args.data_dir else U.DATA_DIR
        args.out_dir = str((base_dir / "HE_unet_models") if args.smoke else U.MODELS_DIR)
    return args


def evaluate(model: torch.nn.Module, samples: list[U.DonorSample],
             device: torch.device) -> dict:
    """Whole-image semantic IoU on the validation donors (fast; no TTA)."""
    model.eval()
    per_donor = {}
    for s in samples:
        # On-disk training image is already single-channel H — read it as such.
        img = np.array(Image.open(s.image_path).convert("L"))
        mask = np.array(Image.open(s.mask_path).convert("L"))
        gt = (mask > 127).astype(np.uint8)
        prob = U.predict_wsi(model, img, device, tta=False)
        per_donor[s.name] = U.iou_score(prob, gt)
    return per_donor


def evaluate_instances(model: torch.nn.Module, samples: list[U.DonorSample],
                       device: torch.device) -> dict:
    """Deployment-faithful follicle count + mean area vs. GT (TTA on, watershed).

    Run once on the best checkpoint — this is the metric H18 depends on, unlike
    the per-epoch semantic IoU used for model selection.
    """
    model.eval()
    per_donor = {}
    for s in samples:
        px = float(s.meta["target_px_size_um"])
        # On-disk training image is already single-channel H — read it as such.
        img = np.array(Image.open(s.image_path).convert("L"))
        gt = np.array(Image.open(s.mask_path).convert("L"))
        prob = U.predict_wsi(model, img, device, tta=True)
        metrics = U.instance_metrics(prob, gt, px_size_um=px)
        metrics["sample_id"] = s.identity
        per_donor[s.name] = metrics
    return per_donor


def load_all_samples(args: argparse.Namespace) -> list[U.DonorSample]:
    """Discover donors from --data-dir plus any --extra-data-dir.

    Names must be unique across dirs: the dataset caches by name and the CV
    splits by donor, so a silent collision would corrupt the cache and hide a
    train/val leak at the same time.
    """
    dirs = [Path(args.data_dir) if args.data_dir else U.DATA_DIR]
    dirs += [Path(d) for d in args.extra_data_dir]
    samples: list[U.DonorSample] = []
    seen: dict[str, Path] = {}
    for d in dirs:
        found = U.discover_samples(d)
        for s in found:
            if s.name in seen:
                raise SystemExit(
                    f"Duplicate donor name {s.name!r} in {d} and {seen[s.name]}. "
                    "Rename one export — names must be unique across pooled dirs.")
            seen[s.name] = d
        samples.extend(found)
        if len(dirs) > 1:
            print(f"  {d}: {len(found)} donors")
    return samples


def warm_start(model: torch.nn.Module, ckpt_path: str) -> dict:
    """Load only shape-matching weights from another checkpoint.

    A DAPI checkpoint has a `1 + K` channel head against this model's single
    channel, so a strict load is impossible by construction. Report what was
    skipped rather than dropping it silently — a warm start that quietly loaded
    almost nothing would look like a training failure much later.
    """
    try:
        ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=True)
    except TypeError:
        ckpt = torch.load(ckpt_path, map_location="cpu")
    src = ckpt.get("state_dict", ckpt)
    own = model.state_dict()
    loaded, skipped = [], []
    for k, v in src.items():
        if k in own and own[k].shape == v.shape:
            own[k] = v
            loaded.append(k)
        else:
            skipped.append(k)
    model.load_state_dict(own)
    frac = len(loaded) / max(1, len(loaded) + len(skipped))
    print(f"  warm start from {ckpt_path}: loaded {len(loaded)} tensors, "
          f"skipped {len(skipped)} ({frac:.0%} transferred)")
    if skipped:
        print(f"    skipped e.g. {skipped[:4]}")
    if frac < 0.5:
        # Almost always an encoder mismatch (e.g. a b0 checkpoint into a b3
        # model). Training still runs and merely looks disappointing rather than
        # broken, so say it loudly here.
        print(f"    WARNING: only {frac:.0%} of the checkpoint transferred. "
              "Check that --encoder matches the checkpoint's; a mismatched "
              "encoder skips nearly everything and the run is effectively "
              "training from scratch.")
    return {"init_from": str(ckpt_path), "loaded": len(loaded),
            "skipped": len(skipped), "fraction_transferred": round(frac, 4)}


def train_one_fold(fold: int, train_samples: list[U.DonorSample],
                   val_samples: list[U.DonorSample],
                   args: argparse.Namespace) -> dict:
    device = torch.device(args.device)
    U.set_seed(args.seed + fold)  # reproducible decoder init + shuffle, per fold
    out_dir = Path(args.out_dir)
    print(f"\n=== Fold {fold} | train={len(train_samples)} val={len(val_samples)} | "
          f"device={device} ===")

    train_ds = U.FolliclePatchDataset(
        train_samples, length=args.patches_per_epoch,
        augment=True, seed=args.seed + fold,
    )
    loader_gen = torch.Generator()
    loader_gen.manual_seed(args.seed + fold)
    train_loader = DataLoader(
        train_ds, batch_size=args.batch_size, shuffle=True,
        num_workers=args.num_workers, pin_memory=(device.type == "cuda"),
        generator=loader_gen,
        worker_init_fn=(U.seed_worker if args.num_workers > 0 else None),
    )

    model = U.build_unet(args.encoder, pretrained=args.pretrained).to(device)
    warm = warm_start(model, args.init_from) if args.init_from else None
    optim = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(
        optim, T_max=args.lr_tmax, eta_min=1e-6)
    scaler = torch.amp.GradScaler(device.type, enabled=(device.type == "cuda"))

    best = {"epoch": -1, "mean_iou": -1.0, "per_donor": {}}
    history: list[dict] = []
    epochs_since_improve = 0
    stopped_early = False
    last_epoch = -1
    n_batches = len(train_loader)
    for epoch in range(args.epochs):
        last_epoch = epoch
        model.train()
        t0 = time.time()
        running = 0.0
        pbar = tqdm(train_loader, desc=f"fold {fold} epoch {epoch:03d}",
                    leave=False, mininterval=2.0)
        for i, (img, mask) in enumerate(pbar, start=1):
            img = img.to(device, non_blocking=True)
            mask = mask.to(device, non_blocking=True)
            optim.zero_grad(set_to_none=True)
            with torch.amp.autocast(device.type, enabled=(device.type == "cuda")):
                logits = model(img)
                loss = U.bce_dice_loss(logits, mask)
            scaler.scale(loss).backward()
            scaler.step(optim)
            scaler.update()
            running += loss.item()
            pbar.set_postfix(loss=f"{running / i:.4f}")
            if i % HEARTBEAT_EVERY == 0:
                tqdm.write(f"    [fold {fold} epoch {epoch:03d}] "
                           f"batch {i}/{n_batches}  "
                           f"loss_avg={running / i:.4f}  "
                           f"elapsed={time.time() - t0:.0f}s")
        epoch_lr = optim.param_groups[0]["lr"]  # rate used THIS epoch
        sched.step()
        train_loss = running / max(1, n_batches)

        ious = evaluate(model, val_samples, device)
        mean_iou = float(np.mean(list(ious.values()))) if ious else 0.0
        print(f"  epoch {epoch:03d}  loss={train_loss:.4f}  val IoU={mean_iou:.3f}"
              f"  lr={epoch_lr:.2e}  ({time.time() - t0:.1f}s)")
        history.append({"epoch": epoch, "train_loss": round(train_loss, 6),
                        "val_iou": round(mean_iou, 6), "lr": epoch_lr})

        if mean_iou > best["mean_iou"]:
            best = {"epoch": epoch, "mean_iou": mean_iou, "per_donor": ious}
            out_dir.mkdir(parents=True, exist_ok=True)
            torch.save({
                "state_dict": model.state_dict(),
                "encoder": args.encoder,
                "fold": fold,
                "epoch": epoch,
                "val_iou": mean_iou,
                "val_per_donor": ious,
            }, out_dir / f"fold_{fold}.pt")
            epochs_since_improve = 0
        else:
            epochs_since_improve += 1
            if args.patience and epochs_since_improve >= args.patience:
                stopped_early = True
                print(f"  early stop: no val-IoU improvement for "
                      f"{args.patience} epoch(s) (stopped at epoch {epoch})")
                break

    best["history"] = history
    best["stopped_early"] = stopped_early
    best["stopped_at_epoch"] = last_epoch
    if warm:
        best["warm_start"] = warm

    # Deployment-faithful instance metrics on the BEST checkpoint (not last).
    if best["epoch"] >= 0:
        ckpt_path = out_dir / f"fold_{fold}.pt"
        try:
            ckpt = torch.load(ckpt_path, map_location=device, weights_only=True)
        except TypeError:
            ckpt = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ckpt["state_dict"])
        best["instance_metrics"] = evaluate_instances(model, val_samples, device)

    print(f"  fold {fold} best: epoch={best['epoch']} mean IoU={best['mean_iou']:.3f}"
          f"{' (early-stopped)' if stopped_early else ''}")
    return best


def main() -> None:
    args = parse_args()
    if args.smoke:
        print("=" * 70)
        print("SMOKE MODE: tiny epochs/patches/batch, CPU, no pretrained weights.")
        print("This is a crash-check / quick estimate, NOT a quality run.")
        print("=" * 70)
    U.set_seed(args.seed)
    samples = load_all_samples(args)
    if not samples:
        raise SystemExit(
            f"No training pairs in {args.data_dir or U.DATA_DIR}. Run "
            "export_HE_training_pairs.groovy (or make_synthetic_he_data.py) first.")

    # Data-integrity preflight: catch coarse-outline / empty masks before training.
    reports = U.check_mask_integrity(samples)
    for r in reports:
        if r["status"].startswith("suspicious"):
            print(f"WARNING: {r['name']} {r['status']} "
                  f"(pos_frac={r['pos_frac']:.3f}, n_components={r['n_components']})")
    empty = [r["name"] for r in reports if r["status"] == "empty"]
    if empty:
        raise SystemExit(f"Empty masks (no Follicle pixels): {empty}. "
                         "Re-export or remove these donors.")

    if args.folds == 1:
        rng = np.random.default_rng(args.seed)
        # Split by DONOR, not by file. With --extra-data-dir a donor can
        # contribute several files, and splitting files would put two exports of
        # the same donor on both sides of the split.
        by_donor: dict[str, list[U.DonorSample]] = {}
        for s in samples:
            by_donor.setdefault(s.identity, []).append(s)
        keys = list(by_donor)
        order = rng.permutation(len(keys))
        n_val = max(1, len(keys) // 5)
        val = [s for i in order[:n_val] for s in by_donor[keys[i]]]
        train = [s for i in order[n_val:] for s in by_donor[keys[i]]]
        splits = [(train, val)]
    else:
        splits = U.donor_kfold(
            samples, k=args.folds,
            genotype_xlsx=Path(args.genotype_xlsx) if args.genotype_xlsx else None,
            seed=args.seed,
        )

    results = []
    for i, (train, val) in enumerate(splits):
        # Files and donor IDs are recorded separately: they differ whenever one
        # donor contributes several files (a pooled run with --extra-data-dir),
        # and a reader that counted `*_donors` would otherwise treat one donor
        # as several independent ones.
        results.append({
            "fold": i,
            "train_files": [s.name for s in train],
            "val_files": [s.name for s in val],
            "train_donors": sorted({s.identity for s in train}),
            "val_donors": sorted({s.identity for s in val}),
            "best": train_one_fold(i, train, val, args),
        })

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "cv_summary.json").write_text(json.dumps(results, indent=2))
    mean = float(np.mean([r["best"]["mean_iou"] for r in results]))
    print(f"\nCV mean IoU across {len(results)} folds: {mean:.3f}")

    # Pooled deployment metric: per-donor predicted vs GT follicle count/area.
    # Spearman needs >=3 points; the single-fold (n_val=1) smoke path can't
    # produce it, so we report the per-donor pairs and skip the correlation.
    pooled = [m for r in results
              for m in r["best"].get("instance_metrics", {}).values()]
    if len(pooled) >= 3:
        from scipy.stats import spearmanr
        rho_c = float(spearmanr([m["gt_count"] for m in pooled],
                                [m["pred_count"] for m in pooled]).correlation)
        rho_a = float(spearmanr([m["gt_mean_area_um2"] for m in pooled],
                                [m["pred_mean_area_um2"] for m in pooled]).correlation)
        print(f"Pooled follicle-count Spearman rho={rho_c:.3f}  "
              f"mean-area rho={rho_a:.3f}  (n={len(pooled)} donors)")
        (out_dir / "pooled_metrics.json").write_text(json.dumps({
            "n_donors": len(pooled),
            "count_spearman": rho_c,
            "area_spearman": rho_a,
            "per_donor": pooled,
        }, indent=2))
    else:
        print(f"Pooled count/area Spearman skipped "
              f"(need >=3 val donors, have {len(pooled)}).")


if __name__ == "__main__":
    main()
