"""Generate tiny synthetic H&E follicle data for smoke tests and CI.

Writes (H image, mask, tissue, meta) PNG/JSON tuples in the same
`data/{images,masks,tissue,meta}` layout that
`export_HE_training_pairs.groovy` produces (single-channel hematoxylin
images, not RGB), so
`train_HE_unet.py --smoke --data-dir <out>` and the unit tests can exercise the
full pipeline with no real H&E data and no network access.

  python scripts/make_synthetic_he_data.py --out-dir /tmp/he_synth --n-donors 3

This is the single source of truth for synthetic donors — `tests/` imports
`make_synthetic_donor` from here rather than duplicating the logic.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image

TARGET_UM_PER_PX = 2.0


def make_synthetic_donor(out_dir: Path, name: str, size: int = 768,
                         n_follicles: int = 4, radius: int = 45,
                         seed: int = 0,
                         sample_id: str | None = None) -> tuple[str, int]:
    """Write one (H image, mask, tissue, meta) tuple under
    `out_dir/{images,masks,tissue,meta}`.

    The image is a single-channel hematoxylin (H) grayscale (high H = dark), as
    produced by `export_HE_training_pairs.groovy`: light tissue background with
    dark follicle disks. Follicles are disks of `radius` px. At the 2 µm/px
    target a 45 px disk is π·(45·2)² ≈ 25 000 µm², comfortably above the
    10 000 µm² post-processing floor. Returns `(name, n_follicles_placed)`.
    """
    out_dir = Path(out_dir)
    rng = np.random.default_rng(seed)
    # Light tissue background (low H), kept under the 220 whitespace threshold so
    # the patch sampler / tissue mask treat the whole slide as tissue.
    img = rng.integers(175, 211, size=(size, size)).astype(np.uint8)
    mask = np.zeros((size, size), dtype=np.uint8)
    yy, xx = np.ogrid[:size, :size]

    # Place roughly non-overlapping follicles so connected components == count.
    margin = radius + 5
    min_sep = (2 * radius + 10) ** 2
    centres: list[tuple[int, int]] = []
    attempts = 0
    while len(centres) < n_follicles and attempts < 200:
        attempts += 1
        cy = int(rng.integers(margin, size - margin))
        cx = int(rng.integers(margin, size - margin))
        if all((cy - py) ** 2 + (cx - px) ** 2 > min_sep for py, px in centres):
            centres.append((cy, cx))
    for cy, cx in centres:
        disk = (yy - cy) ** 2 + (xx - cx) ** 2 < radius ** 2
        mask[disk] = 255
        # Follicles are hematoxylin-dense → dark in the H channel (~70–110).
        img[disk] = rng.integers(70, 111, size=int(disk.sum())).astype(np.uint8)

    tissue = (img < 220).astype(np.uint8) * 255

    images_dir, masks_dir, tissue_dir, meta_dir = (
        out_dir / "images", out_dir / "masks",
        out_dir / "tissue", out_dir / "meta")
    for d in (images_dir, masks_dir, tissue_dir, meta_dir):
        d.mkdir(parents=True, exist_ok=True)
    Image.fromarray(img).save(images_dir / f"{name}.png")
    Image.fromarray(mask).save(masks_dir / f"{name}.png")
    Image.fromarray(tissue).save(tissue_dir / f"{name}.png")
    (meta_dir / f"{name}.json").write_text(json.dumps({
        "image": name,
        "sample_id": sample_id or name,
        "source_name": f"{name}.ome.tif - synthetic",
        "image_type": "hematoxylin",
        "channels": 1,
        "tissue_mask": f"tissue/{name}.png",
        "native_px_size_um": TARGET_UM_PER_PX / 4.0,
        "target_px_size_um": TARGET_UM_PER_PX,
        "downsample": 4.0,
        "width_native_px": size * 4,
        "height_native_px": size * 4,
        "width_export_px": size,
        "height_export_px": size,
        "n_follicle_annotations": len(centres),
    }, indent=2))
    return name, len(centres)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n-donors", type=int, default=3)
    p.add_argument("--size", type=int, default=768)
    p.add_argument("--n-follicles", type=int, default=4)
    p.add_argument("--radius", type=int, default=45)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    out = Path(args.out_dir)
    for i in range(args.n_donors):
        name = f"HDL{900 + i}_synthetic"
        _, n = make_synthetic_donor(out, name, size=args.size,
                                    n_follicles=args.n_follicles,
                                    radius=args.radius, seed=args.seed + i,
                                    sample_id=f"SYNTH_HE_{i:03d}")
        print(f"wrote {name} ({n} follicles) -> {out}")


if __name__ == "__main__":
    main()
