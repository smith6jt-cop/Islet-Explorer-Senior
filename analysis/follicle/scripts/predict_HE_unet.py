"""Run the trained H&E follicle U-Net over whole-slide images.

  python scripts/predict_HE_unet.py --wsi-dir /path/to/HE_wsis \
      --models analysis/HE_unet_models/fold_*.pt

For each WSI (OME-TIFF, SVS, NDPI — anything OpenSlide can open) it:
  1. Reads the slide at TARGET_UM_PER_PX.
  2. Runs ensemble inference with TTA across the supplied checkpoint(s).
  3. Post-processes into instance polygons with watershed splitting.
  4. Writes per-image GeoJSON (µm coords) into analysis/HE_unet_predictions/
     and appends rows to Measurements/HE_follicles_unet.csv.

The CSV uses the same schema as 27_donor_spleen_measurements1.csv so
H18_HE_follicles.ipynb can read it by changing one filename.
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import torch

import HE_unet_utils as U


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--wsi-dir", required=True,
                   help="Directory of H&E WSIs to process.")
    p.add_argument("--models", nargs="+", required=True,
                   help="Checkpoint paths or glob patterns "
                        "(e.g. analysis/HE_unet_models/fold_*.pt).")
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--min-area-um2", type=float, default=10_000)
    p.add_argument("--no-tta", action="store_true")
    p.add_argument("--no-watershed", action="store_true")
    p.add_argument("--geojson-units", choices=["um", "px"], default="um",
                   help="um (default, unchanged) or px. QuPath's Objects -> "
                        "Import objects expects full-resolution PIXEL "
                        "coordinates, so the existing µm files do not import "
                        "correctly; pass px when that is what you need.")
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument("--csv-out", default=str(U.PROJECT / "Measurements"
                                            / "HE_follicles_unet.csv"))
    p.add_argument("--extensions", nargs="+",
                   default=[".ome.tif", ".ome.tiff", ".tif", ".tiff",
                            ".svs", ".ndpi"])
    return p.parse_args()


def load_models(patterns: list[str], device: torch.device) -> list[torch.nn.Module]:
    paths = []
    for pat in patterns:
        matches = sorted(glob.glob(pat))
        paths.extend(matches if matches else [pat])
    models = []
    for p in paths:
        # weights_only=True refuses arbitrary pickle code on load. Falls back
        # to the unsafe path for older torch — only do this on trusted files.
        try:
            ckpt = torch.load(p, map_location=device, weights_only=True)
        except TypeError:
            ckpt = torch.load(p, map_location=device)
        model = U.build_unet(ckpt.get("encoder", "efficientnet-b3"))
        model.load_state_dict(ckpt["state_dict"])
        model.to(device).eval()
        models.append(model)
        print(f"loaded {p}  (val IoU={ckpt.get('val_iou', float('nan')):.3f})")
    if not models:
        raise SystemExit("No checkpoints matched --models patterns.")
    return models


def read_wsi_at_target(path: Path) -> tuple[np.ndarray, float]:
    """Return (rgb_array, effective_px_size_um) at U.TARGET_UM_PER_PX."""
    from PIL import Image

    # Prefer OpenSlide when available — covers SVS/NDPI and most OME-TIFFs.
    try:
        import openslide
    except ImportError:
        openslide = None

    if openslide is not None:
        slide = None
        try:
            slide = openslide.OpenSlide(str(path))
            mpp_x = float(slide.properties.get(openslide.PROPERTY_NAME_MPP_X, "nan"))
            if not np.isfinite(mpp_x) or mpp_x <= 0:
                mpp_x = _ome_pixel_size_um(path)
            if not np.isfinite(mpp_x) or mpp_x <= 0:
                raise RuntimeError(f"{path.name}: missing pixel calibration.")
            target_level = slide.get_best_level_for_downsample(
                U.TARGET_UM_PER_PX / mpp_x)
            level_ds = slide.level_downsamples[target_level]
            level_w, level_h = slide.level_dimensions[target_level]
            img = slide.read_region(
                (0, 0), target_level, (level_w, level_h)).convert("RGB")
            eff_px = mpp_x * level_ds
            if abs(eff_px - U.TARGET_UM_PER_PX) / U.TARGET_UM_PER_PX > 0.05:
                scale = eff_px / U.TARGET_UM_PER_PX
                new_w = int(round(level_w * scale))
                new_h = int(round(level_h * scale))
                img = img.resize((new_w, new_h), Image.BILINEAR)
                eff_px = U.TARGET_UM_PER_PX
            return np.array(img, dtype=np.uint8), eff_px
        except Exception as e:  # noqa: BLE001
            print(f"  OpenSlide failed on {path.name} ({e}); falling back to tifffile")
        finally:
            if slide is not None:
                slide.close()

    # Fallback for plain OME-TIFFs / single-resolution TIFFs.
    import tifffile
    with tifffile.TiffFile(path) as tif:
        page = tif.pages[0]
        arr = page.asarray()
    px = _ome_pixel_size_um(path)
    if arr.ndim == 2:
        arr = np.stack([arr] * 3, axis=-1)
    if arr.shape[-1] > 3:
        arr = arr[..., :3]
    if arr.dtype != np.uint8:
        arr = (arr.astype(np.float32) / arr.max() * 255).astype(np.uint8)
    if not np.isfinite(px) or px <= 0:
        raise RuntimeError(f"{path.name}: missing pixel calibration.")
    scale = px / U.TARGET_UM_PER_PX
    new_w = int(round(arr.shape[1] * scale))
    new_h = int(round(arr.shape[0] * scale))
    rgb = np.array(
        Image.fromarray(arr).resize((new_w, new_h), Image.BILINEAR),
        dtype=np.uint8,
    )
    return rgb, U.TARGET_UM_PER_PX


def _ome_pixel_size_um(path: Path) -> float:
    """Parse PhysicalSizeX from OME-XML; return NaN if absent."""
    try:
        import tifffile
        with tifffile.TiffFile(path) as tif:
            desc = tif.ome_metadata or ""
        m = re.search(r'PhysicalSizeX="([0-9eE.+\-]+)"', desc)
        if m:
            return float(m.group(1))
    except Exception:  # noqa: BLE001
        pass
    return float("nan")


@torch.no_grad()
def ensemble_predict(models: list[torch.nn.Module], himg: np.ndarray,
                     device: torch.device, tta: bool) -> np.ndarray:
    prob = None
    for m in models:
        p = U.predict_wsi(m, himg, device, tta=tta)
        prob = p if prob is None else prob + p
    return prob / len(models)


def main() -> None:
    args = parse_args()
    device = torch.device(args.device)
    models = load_models(args.models, device)
    wsi_dir = Path(args.wsi_dir)

    wsis = []
    for ext in args.extensions:
        wsis.extend(sorted(wsi_dir.rglob(f"*{ext}")))
    if not wsis:
        raise SystemExit(f"No WSIs under {wsi_dir} matching {args.extensions}")

    U.PRED_DIR.mkdir(parents=True, exist_ok=True)
    all_records: list[dict] = []
    for wsi in wsis:
        try:
            rgb, eff_px = read_wsi_at_target(wsi)
        except Exception as e:  # noqa: BLE001
            print(f"  skip {wsi.name}: {e}")
            continue
        # Convert the normalized whole-slide RGB to single-channel H once, at the
        # inference boundary — the model trained on H images from the export.
        himg = U.rgb_to_h(rgb)
        prob = ensemble_predict(models, himg, device, tta=not args.no_tta)
        records, polygons = U.extract_follicle_records(
            prob, px_size_um=eff_px,
            threshold=args.threshold,
            min_area_um2=args.min_area_um2,
            split_touching=not args.no_watershed,
            image_name=wsi.name,
        )
        all_records.extend(records)
        scale = 1.0 if args.geojson_units == "um" else 1.0 / eff_px
        U.write_geojson(polygons, U.PRED_DIR / f"{wsi.stem}.geojson",
                        image_name=wsi.name, scale=scale)
        (U.PRED_DIR / f"{wsi.stem}.meta.json").write_text(json.dumps({
            "image": wsi.name,
            "effective_px_size_um": eff_px,
            "geojson_units": args.geojson_units,
            "n_follicles": len(records),
            "threshold": args.threshold,
            "min_area_um2": args.min_area_um2,
            "n_models": len(models),
            "tta": not args.no_tta,
            "watershed": not args.no_watershed,
        }, indent=2))
        print(f"{wsi.name}: {len(records)} follicles")

    out_csv = Path(args.csv_out)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_records).to_csv(out_csv, index=False)
    print(f"\nWrote {len(all_records)} follicles across {len(wsis)} WSIs → {out_csv}")


if __name__ == "__main__":
    main()
