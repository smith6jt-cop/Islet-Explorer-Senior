"""Train the DAPI/Phenocycler follicle U-Net with donor-stratified k-fold CV.

Assumes `export_DAPI_training_pairs.groovy` has populated data_dapi/ with
(images, masks, tissue, regions, meta) tuples.

  # Single-fold quality estimate first — find the plateau before committing:
  python scripts/train_DAPI_unet.py --folds 1 --epochs 60 --patience 12

  # Reportable run (14 donors; 4 folds keeps >=3 val donors per fold, which is
  # what the pooled Spearman needs):
  python scripts/train_DAPI_unet.py --folds 4 --epochs 50

  # Pooled with the H&E cohort — H&E donors have no regions/, so the auxiliary
  # region head simply sees no gradient from them:
  python scripts/train_DAPI_unet.py --data-dir data_dapi --extra-data-dir data

  # Fast smoke test — synthetic data, CPU, no pretrained weights, no network:
  python scripts/make_synthetic_dapi_data.py --out-dir /tmp/dapi --n-donors 3
  python scripts/train_DAPI_unet.py --smoke --data-dir /tmp/dapi

The model has one decoder and `1 + K` output channels: channel 0 is the binary
Follicle logit (the deliverable, and what model selection is scored on) and
channels 1.. are the auxiliary region softmax. Selecting on the primary head
means a checkpoint can never win on region accuracy while the follicle
segmentation degrades.
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

import DAPI_unet_utils as D
import HE_unet_utils as U

HEARTBEAT_EVERY = 50

NORMAL_DEFAULTS = dict(epochs=50, folds=4, batch_size=16, patches_per_epoch=2000,
                       patience=0, pretrained=True)
SMOKE_DEFAULTS = dict(epochs=2, folds=1, batch_size=2, patches_per_epoch=8,
                      patience=1, pretrained=False)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--epochs", type=int, default=None, help="(default 50; 2 under --smoke)")
    p.add_argument("--folds", type=int, default=None, help="(default 4; 1 under --smoke)")
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
                   help="cosine-annealing period in epochs (default: --epochs)")
    p.add_argument("--encoder", default="efficientnet-b3")
    p.add_argument("--genotype-xlsx",
                   default=str(U.CODE_ROOT / "Spleen_rs3184504_Genotypes.xlsx"))
    p.add_argument("--data-dir", default=None,
                   help="directory holding {images,masks,tissue,regions,meta} "
                        "(default: repo data_dapi/)")
    p.add_argument("--extra-data-dir", action="append", default=[],
                   help="additional data dir(s) to pool in, e.g. the H&E `data/`. "
                        "Repeatable. Donors without a regions/ map contribute to "
                        "the Follicle head only.")
    p.add_argument("--out-dir", default=None)
    p.add_argument("--aux-weight", type=float, default=0.3,
                   help="lambda on the auxiliary region loss (0 disables it)")
    p.add_argument("--no-aux", action="store_true",
                   help="train the binary Follicle head only (K=0), giving a "
                        "model architecturally identical to the H&E one")
    p.add_argument("--merge-white-pulp", action="store_true",
                   help="collapse Follicle+PALS into one WhitePulp region class. "
                        "DAPI cannot see CD20 vs CD3e, so this is the honest "
                        "ablation for how much of the region error is that split.")
    p.add_argument("--init-from", default=None,
                   help="checkpoint to warm-start from; only shape-matching "
                        "weights are loaded (e.g. an H&E fold_0.pt)")
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--num-workers", type=int, default=0,
                   help="0 by default: the dataset caches whole-donor arrays in "
                        "RAM, which every forked worker would duplicate")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--device", default=None)
    return resolve_defaults(p.parse_args())


def resolve_defaults(args: argparse.Namespace) -> argparse.Namespace:
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
        base_dir = Path(args.data_dir) if args.data_dir else D.DATA_DIR
        args.out_dir = str((base_dir / "DAPI_unet_models") if args.smoke else D.MODELS_DIR)
    if args.no_aux:
        args.aux_weight = 0.0
    return args


def region_class_names(args: argparse.Namespace) -> list[str]:
    if args.no_aux:
        return []
    return D.MERGED_REGION_CLASSES if args.merge_white_pulp else D.REGION_CLASSES


def load_all_samples(args: argparse.Namespace) -> tuple[list[U.DonorSample], dict[str, Path]]:
    """Discover the DAPI donors plus any pooled extra dirs, and map region maps.

    Names must be unique across dirs — the dataset caches by name and
    `donor_kfold` splits by donor, so a silent collision would both corrupt the
    cache and hide a leak. Collisions abort rather than get renamed.
    """
    dirs = [Path(args.data_dir) if args.data_dir else D.DATA_DIR]
    dirs += [Path(d) for d in args.extra_data_dir]
    samples: list[U.DonorSample] = []
    region_paths: dict[str, Path] = {}
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
        region_paths.update(D.attach_regions(found, d))
        print(f"  {d}: {len(found)} donors, "
              f"{sum(1 for s in found if s.name in region_paths)} with region maps")
    return samples, region_paths


def evaluate(model: torch.nn.Module, samples: list[U.DonorSample],
             device: torch.device) -> dict:
    """Whole-image Follicle IoU on the validation donors (fast; no TTA)."""
    model.eval()
    per_donor = {}
    for s in samples:
        img = np.array(Image.open(s.image_path).convert("L"))
        mask = np.array(Image.open(s.mask_path).convert("L"))
        gt = (mask > 127).astype(np.uint8)
        # Stitch channel 0 only. The model still computes the region channels,
        # but accumulating them for a whole slide every epoch would cost ~3 GB
        # per donor for a metric that reads one channel.
        prob = U.predict_wsi(model, img, device, tta=False, n_out=1,
                             activation=D.follicle_activation)
        per_donor[s.name] = U.iou_score(prob, gt)
    return per_donor


def evaluate_instances(model: torch.nn.Module, samples: list[U.DonorSample],
                       device: torch.device, n_region_classes: int,
                       region_paths: dict[str, Path], merge_wp: bool,
                       class_names: list[str]) -> tuple[dict, dict]:
    """Follicle count + mean area vs GT (TTA on, watershed), and region confusion.

    Run once on the best checkpoint, not per epoch. Region metrics are computed
    here rather than during whole-slide inference because they need the full K
    probability channels, which `predict_DAPI_unet.py` deliberately does not
    stitch; `region_confusion_banded` keeps even this pass bounded in memory.
    """
    model.eval()
    per_donor, per_region = {}, {}
    for s in samples:
        px = float(s.meta["target_px_size_um"])
        img = np.array(Image.open(s.image_path).convert("L"))
        gt = np.array(Image.open(s.mask_path).convert("L"))
        foll = U.predict_wsi(model, img, device, tta=True, n_out=1,
                             activation=D.follicle_activation)
        metrics = U.instance_metrics(foll, gt, px_size_um=px)
        metrics["sample_id"] = s.identity
        per_donor[s.name] = metrics
        del foll
        rp = region_paths.get(s.name)
        if n_region_classes and rp is not None:
            gt_reg = np.array(Image.open(rp).convert("L"), dtype=np.uint8)
            if merge_wp:
                gt_reg = D.merge_white_pulp(gt_reg)
            # Banded so the K-channel map is never held for the whole slide.
            per_region[s.name] = D.region_confusion_banded(
                model, img, gt_reg, device, n_region_classes, class_names)
    return per_donor, per_region


def warm_start(model: torch.nn.Module, ckpt_path: str) -> dict:
    """Load only shape-matching weights from another checkpoint.

    The output head differs between the H&E model (1 channel) and this one
    (1 + K), so a strict load is impossible by construction. Report what was
    skipped instead of silently dropping it — a warm start that quietly loaded
    almost nothing would look like a training failure later.
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
        # model). Training will still run and will look merely disappointing
        # rather than broken, so say it loudly here.
        print(f"    WARNING: only {frac:.0%} of the checkpoint transferred. "
              "Check that --encoder matches the checkpoint's; a mismatched "
              "encoder skips nearly everything and the run is effectively "
              "training from scratch.")
    return {"init_from": str(ckpt_path), "loaded": len(loaded),
            "skipped": len(skipped), "fraction_transferred": round(frac, 4)}


def train_one_fold(fold: int, train_samples: list[U.DonorSample],
                   val_samples: list[U.DonorSample], region_paths: dict[str, Path],
                   args: argparse.Namespace) -> dict:
    device = torch.device(args.device)
    U.set_seed(args.seed + fold)
    out_dir = Path(args.out_dir)
    class_names = region_class_names(args)
    k = len(class_names)
    n_aux = sum(1 for s in train_samples if s.name in region_paths)
    print(f"\n=== Fold {fold} | train={len(train_samples)} "
          f"({n_aux} with regions) val={len(val_samples)} | device={device} ===")

    # With --no-aux there is no region head, so don't load region maps at all —
    # the dataset would otherwise validate them against a 1-class label space.
    train_regions = ({s.name: region_paths[s.name]
                      for s in train_samples if s.name in region_paths}
                     if k else {})
    weights_t = None
    class_weights = None
    if k:
        # Computed from the TRAINING donors only and stored in the checkpoint, so
        # evaluation and any later fine-tune reuse the same numbers.
        class_weights = D.region_class_weights(train_regions, k,
                                               merge_wp=args.merge_white_pulp)
        weights_t = torch.tensor(class_weights, dtype=torch.float32, device=device)
        print("  region class weights: "
              + ", ".join(f"{n}={w:.2f}" for n, w in zip(class_names, class_weights)))

    train_ds = D.DapiRegionDataset(
        train_samples, train_regions, length=args.patches_per_epoch,
        augment=True, seed=args.seed + fold,
        n_region_classes=max(k, 1), merge_wp=args.merge_white_pulp,
    )
    loader_gen = torch.Generator()
    loader_gen.manual_seed(args.seed + fold)
    train_loader = DataLoader(
        train_ds, batch_size=args.batch_size, shuffle=True,
        num_workers=args.num_workers, pin_memory=(device.type == "cuda"),
        generator=loader_gen,
        worker_init_fn=(U.seed_worker if args.num_workers > 0 else None),
    )

    model = D.build_multitask_unet(args.encoder, pretrained=args.pretrained,
                                   n_region_classes=k).to(device)
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
        for i, (img, mask, region, aux_valid) in enumerate(pbar, start=1):
            img = img.to(device, non_blocking=True)
            mask = mask.to(device, non_blocking=True)
            region = region.to(device, non_blocking=True)
            aux_valid = aux_valid.to(device, non_blocking=True)
            optim.zero_grad(set_to_none=True)
            with torch.amp.autocast(device.type, enabled=(device.type == "cuda")):
                logits = model(img)
                loss = D.multitask_loss(logits, mask,
                                        region=region if k else None,
                                        aux_valid=aux_valid,
                                        lam=args.aux_weight,
                                        class_weights=weights_t)
            scaler.scale(loss).backward()
            scaler.step(optim)
            scaler.update()
            running += loss.item()
            pbar.set_postfix(loss=f"{running / i:.4f}")
            if i % HEARTBEAT_EVERY == 0:
                tqdm.write(f"    [fold {fold} epoch {epoch:03d}] batch {i}/{n_batches}  "
                           f"loss_avg={running / i:.4f}  elapsed={time.time() - t0:.0f}s")
        epoch_lr = optim.param_groups[0]["lr"]
        sched.step()
        train_loss = running / max(1, n_batches)

        ious = evaluate(model, val_samples, device)
        mean_iou = float(np.mean(list(ious.values()))) if ious else 0.0
        print(f"  epoch {epoch:03d}  loss={train_loss:.4f}  val Follicle IoU={mean_iou:.3f}"
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
                "n_region_classes": k,
                "region_classes": class_names,
                "region_class_weights": class_weights,
                "merge_white_pulp": bool(args.merge_white_pulp),
                "target_px_size_um": D.TARGET_UM_PER_PX,
                "modality": "dapi",
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

    if best["epoch"] >= 0:
        ckpt_path = out_dir / f"fold_{fold}.pt"
        try:
            ckpt = torch.load(ckpt_path, map_location=device, weights_only=True)
        except TypeError:
            ckpt = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ckpt["state_dict"])
        inst, regions = evaluate_instances(
            model, val_samples, device, k, region_paths,
            args.merge_white_pulp, class_names)
        best["instance_metrics"] = inst
        if regions:
            best["region_metrics"] = regions

    print(f"  fold {fold} best: epoch={best['epoch']} Follicle IoU={best['mean_iou']:.3f}"
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
    samples, region_paths = load_all_samples(args)
    if not samples:
        raise SystemExit(
            f"No training pairs in {args.data_dir or D.DATA_DIR}. Run "
            "export_DAPI_training_pairs.groovy (or make_synthetic_dapi_data.py) first.")

    reports = U.check_mask_integrity(samples)
    for r in reports:
        if r["status"].startswith("suspicious"):
            print(f"WARNING: {r['name']} {r['status']} "
                  f"(pos_frac={r['pos_frac']:.3f}, n_components={r['n_components']})")
    empty = [r["name"] for r in reports if r["status"] == "empty"]
    if empty:
        raise SystemExit(f"Empty masks (no Follicle pixels): {empty}. "
                         "Re-export or remove these donors.")

    if not args.no_aux:
        for r in D.check_region_integrity(samples, region_paths):
            if r["status"] not in ("ok", "no_regions"):
                print(f"WARNING: {r['name']} region map {r['status']} "
                      f"(classes={r['present_classes']}, "
                      f"follicle_frac={r['follicle_frac']:.5f})")

    donors = sorted({s.identity for s in samples})
    print(f"\n{len(samples)} files across {len(donors)} donors: {', '.join(donors)}")
    if args.folds > 1 and len(donors) // args.folds < 3:
        print(f"NOTE: {len(donors)} donors / {args.folds} folds leaves "
              f"{len(donors) // args.folds} validation donors per fold; the pooled "
              "count/area Spearman needs >=3 pooled val donors to be reported.")

    if args.folds == 1:
        rng = np.random.default_rng(args.seed)
        # Split by DONOR even in the single-fold path — splitting files would put
        # two exports of one donor on both sides of the split.
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
        # Record BOTH the files and the donor IDs. The two differ whenever a
        # donor contributes several files (HDL073_PC19 + HDL073_PC29, or a
        # pooled DAPI+H&E run), and a downstream reader that counted
        # `*_donors` would then treat one donor as several independent ones —
        # the same conflation the donor-level split exists to prevent.
        results.append({
            "fold": i,
            "train_files": [s.name for s in train],
            "val_files": [s.name for s in val],
            "train_donors": sorted({s.identity for s in train}),
            "val_donors": sorted({s.identity for s in val}),
            "best": train_one_fold(i, train, val, region_paths, args),
        })

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "cv_summary.json").write_text(json.dumps(results, indent=2))
    mean = float(np.mean([r["best"]["mean_iou"] for r in results]))
    print(f"\nCV mean Follicle IoU across {len(results)} folds: {mean:.3f}")

    pooled = [m for r in results
              for m in r["best"].get("instance_metrics", {}).values()]
    summary: dict = {"n_donors": len(pooled), "per_donor": pooled}
    if len(pooled) >= 3:
        from scipy.stats import spearmanr
        rho_c = float(spearmanr([m["gt_count"] for m in pooled],
                                [m["pred_count"] for m in pooled]).correlation)
        rho_a = float(spearmanr([m["gt_mean_area_um2"] for m in pooled],
                                [m["pred_mean_area_um2"] for m in pooled]).correlation)
        print(f"Pooled follicle-count Spearman rho={rho_c:.3f}  "
              f"mean-area rho={rho_a:.3f}  (n={len(pooled)} donors)")
        summary["count_spearman"] = rho_c
        summary["area_spearman"] = rho_a
    else:
        print(f"Pooled count/area Spearman skipped "
              f"(need >=3 val donors, have {len(pooled)}).")

    # Region metrics are the honest read on the DAPI-only ceiling: a model can
    # post a good follicle count while systematically swapping Follicle and PALS,
    # which aggregate IoU alone would not reveal.
    reg = [m for r in results for m in r["best"].get("region_metrics", {}).values()]
    if reg:
        names = reg[0]["class_names"]
        # Average each class over the donors where it is actually present; a
        # class absent everywhere stays None rather than becoming a number.
        per_class: dict[str, float | None] = {}
        for n in names:
            vals = [m["per_class_iou"][n] for m in reg
                    if m["per_class_iou"].get(n) is not None]
            per_class[n] = float(np.mean(vals)) if vals else None
        summary["region_mean_iou"] = float(np.mean([m["mean_iou"] for m in reg]))
        summary["region_per_class_iou"] = per_class
        print("Region IoU: " + ", ".join(
            f"{n}=" + ("absent" if v is None else f"{v:.3f}")
            for n, v in per_class.items()))
        if "follicle_as_pals" in reg[0]:
            fap = float(np.mean([m["follicle_as_pals"] for m in reg]))
            paf = float(np.mean([m["pals_as_follicle"] for m in reg]))
            summary["follicle_as_pals"] = fap
            summary["pals_as_follicle"] = paf
            print(f"Follicle->PALS confusion={fap:.3f}  PALS->Follicle={paf:.3f}  "
                  "(DAPI cannot see CD20 vs CD3e; this is the expected ceiling)")
    (out_dir / "pooled_metrics.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
