"""Generate tiny synthetic DAPI training data for smoke tests and CI.

Writes (DAPI image, Follicle mask, tissue, region map, meta) tuples in the same
`data_dapi/{images,masks,tissue,regions,meta}` layout that
`export_DAPI_training_pairs.groovy` produces, so
`train_DAPI_unet.py --smoke --data-dir <out>` and the unit tests exercise the
full pipeline with no real Phenocycler data and no network access.

  python scripts/make_synthetic_dapi_data.py --out-dir /tmp/dapi --n-donors 3

With `--write-ome-tiff` it also emits a small multi-channel pyramidal OME-TIFF
per donor under `<out>/ome/`, following the conventions in
`scripts/convert_to_ome_tiff.py` (CYX, uint16, DAPI first, PhysicalSizeX in µm).
That lets CI exercise `read_dapi_channel` and the whole `predict_DAPI_unet.py`
path — a gap the H&E pipeline never closed, since it has no synthetic WSI.

This is the single source of truth for synthetic DAPI donors — `tests/` imports
`make_synthetic_dapi_donor` from here rather than duplicating the logic.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image

TARGET_UM_PER_PX = 2.0
# The synthetic slide is authored at the TARGET resolution and the OME-TIFF is
# built around it as a 2-level pyramid: full resolution at 1 µm/px, SubIFD level
# at 2 µm/px holding exactly the authored raster. So `read_dapi_channel` has a
# real level choice to make, and when it makes the right one the bytes it
# returns are identical to what the PNGs were derived from — which is what makes
# the round-trip a parity test rather than an approximation.
NATIVE_UM_PER_PX = 1.0
PYRAMID_FACTOR = 2
REGION_CLASSES = ["Background", "RedPulp", "Trabeculae", "LargeVessel", "PALS", "Follicle"]
FOLLICLE_INDEX = REGION_CLASSES.index("Follicle")
PALS_INDEX = REGION_CLASSES.index("PALS")
REDPULP_INDEX = REGION_CLASSES.index("RedPulp")


def make_synthetic_dapi_donor(out_dir: Path, name: str, size: int = 768,
                               n_follicles: int = 4, radius: int = 45,
                               seed: int = 0, write_ome: bool = False,
                               n_channels: int = 4,
                               sample_id: str | None = None) -> tuple[str, int]:
    """Write one synthetic DAPI donor. Returns `(name, n_follicles_placed)`.

    The synthetic slide mimics the real structure the model has to learn:
    RedPulp background at moderate nuclear density, follicles as dense-nuclei
    disks, and a PALS annulus around each follicle. The PALS ring is deliberately
    given the SAME intensity distribution as the follicle core — in a real slide
    DAPI cannot tell white-pulp B-zone from T-zone either, so a synthetic set
    where it could would quietly validate a distinction the model will never
    make on real data.

    Images on disk are already normalised and INVERTED (high DAPI = dark), as
    the QuPath export writes them.
    """
    out_dir = Path(out_dir)
    rng = np.random.default_rng(seed)

    # Raw uint16 DAPI, before normalisation. Dim background, brighter tissue.
    raw = rng.integers(200, 900, size=(size, size)).astype(np.uint16)
    region = np.zeros((size, size), dtype=np.uint8)
    yy, xx = np.ogrid[:size, :size]

    # Tissue occupies a large central block; the rim stays background so the
    # Otsu tissue mask has something real to separate.
    margin_bg = size // 12
    tissue_box = (slice(margin_bg, size - margin_bg), slice(margin_bg, size - margin_bg))
    raw[tissue_box] = rng.integers(9000, 14000, size=raw[tissue_box].shape).astype(np.uint16)
    region[tissue_box] = REDPULP_INDEX

    pals_r = radius + 14
    margin = pals_r + 5
    min_sep = (2 * pals_r + 10) ** 2
    centres: list[tuple[int, int]] = []
    attempts = 0
    while len(centres) < n_follicles and attempts < 400:
        attempts += 1
        cy = int(rng.integers(margin_bg + margin, size - margin_bg - margin))
        cx = int(rng.integers(margin_bg + margin, size - margin_bg - margin))
        if all((cy - py) ** 2 + (cx - px) ** 2 > min_sep for py, px in centres):
            centres.append((cy, cx))

    for cy, cx in centres:
        d2 = (yy - cy) ** 2 + (xx - cx) ** 2
        ring = d2 < pals_r ** 2
        core = d2 < radius ** 2
        # PALS first, Follicle second — same "later wins" order as the exporter.
        region[ring] = PALS_INDEX
        region[core] = FOLLICLE_INDEX
        white_pulp = ring
        raw[white_pulp] = rng.integers(
            40000, 60000, size=int(white_pulp.sum())).astype(np.uint16)

    mask = (region == FOLLICLE_INDEX).astype(np.uint8) * 255

    # Normalise + invert with the SAME code path the exporter mirrors, so the
    # synthetic set is a real exercise of the parity contract, not a shortcut.
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import DAPI_unet_utils as D

    params = D.dapi_transform_params(raw)
    img = D.dapi_to_uint8(raw, params["lo"], params["hi"])
    tissue = D.tissue_from_raw(raw, params["otsu_bin"]).astype(np.uint8) * 255

    dirs = {k: out_dir / k for k in ("images", "masks", "tissue", "regions", "meta")}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    Image.fromarray(img).save(dirs["images"] / f"{name}.png")
    Image.fromarray(mask).save(dirs["masks"] / f"{name}.png")
    Image.fromarray(tissue).save(dirs["tissue"] / f"{name}.png")
    Image.fromarray(region).save(dirs["regions"] / f"{name}.png")

    px_area = TARGET_UM_PER_PX ** 2
    counts = np.bincount(region.ravel(), minlength=len(REGION_CLASSES))
    downsample = float(PYRAMID_FACTOR)
    (dirs["meta"] / f"{name}.json").write_text(json.dumps({
        "image": name,
        "sample_id": sample_id or name,
        "source_name": f"{name}.ome.tiff - synthetic",
        "image_type": "dapi",
        "channels": 1,
        "dapi_channel_index": 0,
        "dapi_channel_name": "DAPI",
        "tissue_mask": f"tissue/{name}.png",
        "region_map": f"regions/{name}.png",
        "region_classes": REGION_CLASSES,
        "dapi_transform": params,
        "native_px_size_um": NATIVE_UM_PER_PX,
        "target_px_size_um": TARGET_UM_PER_PX,
        "downsample": downsample,
        "width_native_px": int(size * downsample),
        "height_native_px": int(size * downsample),
        "width_export_px": size,
        "height_export_px": size,
        "n_follicle_annotations": len(centres),
        "class_area_um2": {c: int(counts[i] * px_area) if i < len(counts) else 0
                           for i, c in enumerate(REGION_CLASSES)},
    }, indent=2))

    if write_ome:
        _write_synthetic_ome(out_dir / "ome" / f"{name}.ome.tiff", raw,
                             n_channels=n_channels, rng=rng)
    return name, len(centres)


def _write_synthetic_ome(path: Path, dapi_target: np.ndarray, n_channels: int,
                         rng: np.random.Generator) -> None:
    """Write a small multi-channel pyramidal CYX uint16 OME-TIFF, DAPI first.

    Follows `scripts/convert_to_ome_tiff.py`: **one IFD page per channel** (not
    one multi-sample page), SubIFD pyramid levels, big-endian uint16, ADOBE_DEFLATE,
    and OME-XML passed as `bytes` so the raw UTF-8 µ survives tifffile's
    ASCII-only description check.

    `dapi_target` is the authored raster at TARGET_UM_PER_PX. The full-resolution
    page is a nearest-neighbour 2x expansion of it at NATIVE_UM_PER_PX, and the
    SubIFD level is `dapi_target` verbatim — so a reader that correctly picks the
    2 µm/px level gets the authored pixels back exactly.
    """
    import tifffile

    path.parent.mkdir(parents=True, exist_ok=True)
    h, w = dapi_target.shape
    fh, fw = h * PYRAMID_FACTOR, w * PYRAMID_FACTOR
    names = ["DAPI"] + [f"Marker{i}" for i in range(1, n_channels)]

    def full_res(level_img: np.ndarray) -> np.ndarray:
        return np.repeat(np.repeat(level_img, PYRAMID_FACTOR, axis=0),
                         PYRAMID_FACTOR, axis=1)

    chan_xml = "".join(
        f'<Channel ID="Channel:0:{c}" Name="{names[c]}" SamplesPerPixel="1">'
        f"<LightPath/></Channel>" for c in range(n_channels))
    tiffdata = "".join(
        f'<TiffData FirstC="{c}" FirstT="0" FirstZ="0" IFD="{c}" PlaneCount="1"/>'
        for c in range(n_channels))
    ome = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
        f'<Image ID="Image:0" Name="{path.stem}">'
        f'<Pixels ID="Pixels:0" DimensionOrder="XYCZT" Type="uint16" '
        f'BigEndian="true" Interleaved="false" SizeX="{fw}" SizeY="{fh}" '
        f'SizeC="{n_channels}" SizeZ="1" SizeT="1" '
        f'PhysicalSizeX="{NATIVE_UM_PER_PX}" PhysicalSizeXUnit="µm" '
        f'PhysicalSizeY="{NATIVE_UM_PER_PX}" PhysicalSizeYUnit="µm">'
        f"{chan_xml}{tiffdata}</Pixels></Image></OME>"
    ).encode("utf-8")

    compress = tifffile.COMPRESSION.ADOBE_DEFLATE
    with tifffile.TiffWriter(path, bigtiff=True, byteorder=">") as tw:
        for c in range(n_channels):
            level = (dapi_target if c == 0
                     else rng.integers(0, 5000, size=(h, w)).astype(np.uint16))
            page_kwargs = dict(tile=(256, 256), compression=compress,
                               subifds=1, metadata=None)
            if c == 0:
                page_kwargs["description"] = ome
            tw.write(full_res(level), **page_kwargs)
            tw.write(level, tile=(256, 256), compression=compress,
                     subfiletype=1, metadata=None)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n-donors", type=int, default=3)
    p.add_argument("--size", type=int, default=768)
    p.add_argument("--n-follicles", type=int, default=4)
    p.add_argument("--radius", type=int, default=45)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--write-ome-tiff", action="store_true",
                   help="also emit <out>/ome/<name>.ome.tiff so the multiplex "
                        "reader and predict_DAPI_unet.py can be tested")
    p.add_argument("--n-channels", type=int, default=4,
                   help="channels in the synthetic OME-TIFF (DAPI is channel 0)")
    args = p.parse_args()

    out = Path(args.out_dir)
    for i in range(args.n_donors):
        name = f"HDL{900 + i}_PC33_synthetic"
        _, n = make_synthetic_dapi_donor(
            out, name, size=args.size, n_follicles=args.n_follicles,
            radius=args.radius, seed=args.seed + i,
            write_ome=args.write_ome_tiff, n_channels=args.n_channels,
            sample_id=f"SYNTH_DAPI_{i:03d}")
        print(f"wrote {name} ({n} follicles) -> {out}")


if __name__ == "__main__":
    main()
