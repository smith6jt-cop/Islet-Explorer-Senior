"""Run the trained DAPI follicle U-Net over multiplex Phenocycler slides.

  python scripts/predict_DAPI_unet.py --wsi-dir Images \
      --models 'analysis/DAPI_unet_models/fold_*.pt'

For each slide it:
  1. Reads the DAPI channel at TARGET_UM_PER_PX (tifffile + pyramid levels —
     OpenSlide cannot select one channel of a 33-channel uint16 OME-TIFF).
  2. Normalises + inverts with the same transform the QuPath export applied.
  3. Runs ensemble inference with TTA across the supplied checkpoint(s).
  4. Post-processes into instance polygons with watershed splitting.
  5. Writes per-image GeoJSON and appends rows to
     Measurements/DAPI_follicles_unet.csv.

GeoJSON defaults to **pixel** coordinates here (`--geojson-units px`) so the
output re-imports into QuPath via Objects -> Import objects. The H&E predictions
are in µm, which is why that flag exists rather than a hardcoded convention.

Only the Follicle channel is stitched across the slide. The auxiliary region
channels are a training signal, not a deliverable, and a whole slide at 2 µm/px
(~10^8 px) times six float32 region maps is several GB. Region agreement is
reported by train_DAPI_unet.py on the much smaller 2 µm/px training rasters.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

import DAPI_unet_utils as D
import HE_unet_utils as U


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--wsi-dir", required=True,
                   help="Directory of multiplex OME-TIFFs to process.")
    p.add_argument("--models", nargs="+", required=True,
                   help="Checkpoint paths or glob patterns "
                        "(e.g. analysis/DAPI_unet_models/fold_*.pt).")
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--min-area-um2", type=float, default=10_000)
    p.add_argument("--no-tta", action="store_true")
    p.add_argument("--no-watershed", action="store_true")
    p.add_argument("--geojson-units", choices=["px", "um"], default="px",
                   help="px (default) writes full-resolution pixel coordinates, "
                        "which is what QuPath's object import expects; um matches "
                        "the H&E predictions.")
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument("--csv-out", default=str(U.PROJECT / "Measurements"
                                            / "DAPI_follicles_unet.csv"))
    p.add_argument("--pred-dir", default=str(D.PRED_DIR))
    p.add_argument("--extensions", nargs="+",
                   default=[".ome.tif", ".ome.tiff", ".tif", ".tiff", ".qptiff"])
    return p.parse_args()


def load_models(patterns: list[str], device: torch.device
                ) -> tuple[list[torch.nn.Module], int]:
    """Load checkpoints and return `(models, n_region_classes)`.

    All checkpoints must agree on the head width — averaging a 1-channel model
    with a 1+K-channel one would silently mix the Follicle logit with a region
    logit, so a mismatch is an error rather than a warning.
    """
    paths: list[str] = []
    for pat in patterns:
        matches = sorted(glob.glob(pat))
        paths.extend(matches if matches else [pat])
    models: list[torch.nn.Module] = []
    n_region = None
    for p in paths:
        try:
            ckpt = torch.load(p, map_location=device, weights_only=True)
        except TypeError:
            ckpt = torch.load(p, map_location=device)
        k = int(ckpt.get("n_region_classes", 0))
        if n_region is None:
            n_region = k
        elif k != n_region:
            raise SystemExit(
                f"{p} has {k} region classes but an earlier checkpoint has "
                f"{n_region}. Ensemble members must share a head width.")
        # pretrained=False: the weights are about to be overwritten anyway, and
        # this keeps inference from reaching for an ImageNet download.
        model = D.build_multitask_unet(ckpt.get("encoder", "efficientnet-b3"),
                                       pretrained=False, n_region_classes=k)
        model.load_state_dict(ckpt["state_dict"])
        model.to(device).eval()
        models.append(model)
        print(f"loaded {p}  (val IoU={ckpt.get('val_iou', float('nan')):.3f}, "
              f"{k} region classes)")
    if not models:
        raise SystemExit("No checkpoints matched --models patterns.")
    return models, int(n_region or 0)


@torch.no_grad()
def ensemble_follicle_prob(models: list[torch.nn.Module], img: np.ndarray,
                           device: torch.device, tta: bool) -> np.ndarray:
    """Mean Follicle probability map across the ensemble, channel 0 only.

    Head width deliberately does not appear here: `follicle_activation` slices
    channel 0 whatever K is, so there is nothing for a region-count argument to
    change.
    """
    prob = None
    for m in models:
        # n_out=1 + follicle_activation: the model computes every channel, but
        # only channel 0 is stitched. A slide at 2 µm/px is ~10^8 px, so keeping
        # the region channels too would cost several GB for maps nothing here
        # reads.
        p = U.predict_wsi(m, img, device, tta=tta, n_out=1,
                          activation=D.follicle_activation)
        prob = p if prob is None else prob + p
    return prob / len(models)


def main() -> None:
    args = parse_args()
    device = torch.device(args.device)
    models, n_region = load_models(args.models, device)
    wsi_dir = Path(args.wsi_dir)

    wsis: list[Path] = []
    for ext in args.extensions:
        wsis.extend(sorted(wsi_dir.rglob(f"*{ext}")))
    wsis = sorted(set(wsis))
    if not wsis:
        raise SystemExit(f"No slides under {wsi_dir} matching {args.extensions}")

    pred_dir = Path(args.pred_dir)
    pred_dir.mkdir(parents=True, exist_ok=True)
    all_records: list[dict] = []
    for wsi in wsis:
        try:
            raw, eff_px = D.read_dapi_channel(wsi, D.TARGET_UM_PER_PX)
            # Recompute the transform from this slide's own pixels — the same
            # per-slide rule the export used. Nothing upstream normalises these
            # images, so a fixed global window would put each slide on a
            # different effective scale.
            params = D.dapi_transform_params(raw)
            img = D.dapi_to_uint8(raw, params["lo"], params["hi"])
            del raw
        except Exception as e:  # noqa: BLE001
            print(f"  skip {wsi.name}: {e}")
            continue

        prob = ensemble_follicle_prob(models, img, device, tta=not args.no_tta)
        records, polygons = U.extract_follicle_records(
            prob, px_size_um=eff_px,
            threshold=args.threshold,
            min_area_um2=args.min_area_um2,
            split_touching=not args.no_watershed,
            image_name=wsi.name,
        )
        all_records.extend(records)

        # extract_follicle_records returns µm; px output divides back out by the
        # effective pixel size so coordinates land in full-resolution image space.
        scale = 1.0 if args.geojson_units == "um" else 1.0 / eff_px
        U.write_geojson(polygons, pred_dir / f"{wsi.stem}.geojson",
                        image_name=wsi.name, classification="Follicle",
                        scale=scale)
        (pred_dir / f"{wsi.stem}.meta.json").write_text(json.dumps({
            "image": wsi.name,
            "effective_px_size_um": eff_px,
            "n_follicles": len(records),
            "threshold": args.threshold,
            "min_area_um2": args.min_area_um2,
            "n_models": len(models),
            "n_region_classes": n_region,
            "tta": not args.no_tta,
            "watershed": not args.no_watershed,
            "geojson_units": args.geojson_units,
            "dapi_transform": params,
        }, indent=2))
        print(f"{wsi.name}: {len(records)} follicles "
              f"(lo={params['lo']} hi={params['hi']})")

    out_csv = Path(args.csv_out)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if all_records:
        pd.DataFrame(all_records).to_csv(out_csv, index=False)
    else:
        # Write the header anyway — a headerless empty CSV breaks pd.read_csv
        # downstream, and H18/H19 should see "no rows", not a parse error.
        pd.DataFrame(columns=[
            "Image", "Object ID", "Object type", "Name", "Classification",
            "Parent", "ROI", "Centroid X µm", "Centroid Y µm",
            "Area µm^2", "Perimeter µm"]).to_csv(out_csv, index=False)
    print(f"\nWrote {len(all_records)} follicles across {len(wsis)} slides → {out_csv}")


if __name__ == "__main__":
    main()
