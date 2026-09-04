#!/usr/bin/env python3
"""Convert a directory of single-channel TIF files into a pyramidal OME-TIFF.

Produces a Bio-Formats/QuPath-compatible file:
- Pyramidal OME-TIFF with SubIFDs (6 resolution levels: 1x, 2x, 4x, 8x, 16x, 32x)
- 512x512 tiles, DEFLATE compression, uint16
- CYX axis order (one IFD page per channel)
- Manual OME-XML with per-channel TiffData blocks for correct channel mapping

Usage:
    python convert_to_ome_tiff.py <input_dir> <output_path> [--pixel-size 0.508]
"""

import argparse
import sys
from pathlib import Path
from uuid import uuid4

import numpy as np
import tifffile


NUM_SUBIFDS = 5  # 2x, 4x, 8x, 16x, 32x


def discover_channels(input_dir: Path) -> list[Path]:
    """Find channel TIF files, sort with DAPI first, rename DAPI-01 -> DAPI."""
    tifs = sorted(
        p for p in input_dir.glob("*.tif")
        if p.stem not in ("BLANK",) and not p.stem.startswith("BLANK")
    )
    if not tifs:
        sys.exit(f"No TIF files found in {input_dir}")

    # Sort: DAPI first, then alphabetical
    dapi = [p for p in tifs if p.stem in ("DAPI", "DAPI-01")]
    others = [p for p in tifs if p.stem not in ("DAPI", "DAPI-01")]
    return dapi + others


def channel_name(path: Path) -> str:
    """Extract channel name from filename, normalizing DAPI-01 to DAPI."""
    name = path.stem
    if name == "DAPI-01":
        return "DAPI"
    return name


def downsample_2x(img: np.ndarray) -> np.ndarray:
    """Block-mean 2x downsample of a 2D (Y, X) array, trimming odd dims."""
    h, w = img.shape
    h2, w2 = h - h % 2, w - w % 2
    return img[:h2, :w2].reshape(h2 // 2, 2, w2 // 2, 2).mean(axis=(1, 3)).astype(np.uint16)


def build_ome_xml(
    names: list[str],
    size_y: int,
    size_x: int,
    pixel_size: float,
    filename: str,
) -> bytes:
    """Build Bio-Formats-compatible OME-XML with per-channel TiffData blocks.

    Returns UTF-8 bytes (not str) so callers can pass directly to
    tifffile's ``description`` parameter, bypassing its 7-bit ASCII check.
    Bio-Formats writes UTF-8 µ in the description tag; we must match that.
    """
    n_channels = len(names)
    image_uuid = f"urn:uuid:{uuid4()}"

    # Channel elements (with LightPath to match Bio-Formats output)
    channel_elements = []
    for name in names:
        channel_elements.append(
            f'<Channel ID="Channel:0:{len(channel_elements)}" '
            f'Name="{name}" SamplesPerPixel="1">'
            f'<LightPath/></Channel>'
        )

    # TiffData elements — one per channel, each mapping to its IFD page
    tiffdata_elements = []
    for c in range(n_channels):
        tiffdata_elements.append(
            f'<TiffData FirstC="{c}" FirstT="0" FirstZ="0" '
            f'IFD="{c}" PlaneCount="1">'
            f'<UUID FileName="{filename}">{image_uuid}</UUID>'
            f'</TiffData>'
        )

    ome_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
        f'UUID="{image_uuid}" '
        'xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 '
        'http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">'
        f'<Image ID="Image:0" Name="{filename}">'
        f'<Pixels BigEndian="true" DimensionOrder="XYCZT" '
        f'ID="Pixels:0" Interleaved="false" '
        f'PhysicalSizeX="{pixel_size}" PhysicalSizeXUnit="\u00b5m" '
        f'PhysicalSizeY="{pixel_size}" PhysicalSizeYUnit="\u00b5m" '
        f'SizeC="{n_channels}" SizeT="1" SizeX="{size_x}" '
        f'SizeY="{size_y}" SizeZ="1" Type="uint16">'
        + "".join(channel_elements)
        + "".join(tiffdata_elements)
        + '</Pixels></Image></OME>'
    )
    # Return UTF-8 bytes to bypass tifffile's ASCII-only description check
    return ome_xml.encode('utf-8')


def convert(input_dir: Path, output_path: Path, pixel_size: float) -> None:
    """Convert single-channel TIFs to pyramidal OME-TIFF."""
    channel_paths = discover_channels(input_dir)
    names = [channel_name(p) for p in channel_paths]
    n_channels = len(names)

    print(f"Found {n_channels} channels: {', '.join(names)}")

    # Read first channel to get dimensions
    first = tifffile.imread(str(channel_paths[0]))
    size_y, size_x = first.shape
    del first
    print(f"Image dimensions: {size_y} x {size_x}")

    # Build OME-XML
    filename = output_path.name
    ome_xml = build_ome_xml(names, size_y, size_x, pixel_size, filename)

    # Resolution in pixels per centimeter (matching Bio-Formats convention)
    resolution_ppcm = 1e4 / pixel_size  # µm/pixel → pixels/cm
    resolution = (resolution_ppcm, resolution_ppcm)

    # Use ADOBE_DEFLATE (code 8) — Bio-Formats uses this, not DEFLATE (32946)
    compress = tifffile.COMPRESSION.ADOBE_DEFLATE

    subifd_options = dict(
        tile=(512, 512),
        compression=compress,
        subfiletype=1,
        metadata=None,
    )

    print(f"Writing {output_path} ...")
    with tifffile.TiffWriter(str(output_path), bigtiff=True, byteorder=">") as tw:
        for i, (path, name) in enumerate(zip(channel_paths, names)):
            print(f"  [{i + 1}/{n_channels}] {name}")

            # Load single 2D channel
            img = tifffile.imread(str(path)).astype(np.uint16)

            # Write full-resolution page
            page_kwargs = dict(
                tile=(512, 512),
                compression=compress,
                subifds=NUM_SUBIFDS,
                metadata=None,
                resolution=resolution,
                resolutionunit=3,  # CENTIMETER
                software="OME Bio-Formats 8.2.0",
            )
            if i == 0:
                page_kwargs["description"] = ome_xml
            tw.write(img, **page_kwargs)

            # Write 5 SubIFD pyramid levels (2x, 4x, 8x, 16x, 32x)
            sub = img
            del img
            for _ in range(NUM_SUBIFDS):
                sub = downsample_2x(sub)
                tw.write(sub, **subifd_options)
            del sub

    print(f"Wrote {output_path} ({output_path.stat().st_size / 1e9:.2f} GB)")


def main():
    parser = argparse.ArgumentParser(
        description="Convert single-channel TIFs to pyramidal OME-TIFF"
    )
    parser.add_argument("input_dir", type=Path, help="Directory of channel TIF files")
    parser.add_argument("output_path", type=Path, help="Output OME-TIFF path")
    parser.add_argument(
        "--pixel-size",
        type=float,
        default=0.5077663810243286,
        help="Pixel size in \u00b5m (default: 0.5077663810243286)",
    )
    args = parser.parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"Input directory not found: {args.input_dir}")

    args.output_path.parent.mkdir(parents=True, exist_ok=True)

    convert(args.input_dir, args.output_path, args.pixel_size)


if __name__ == "__main__":
    main()
