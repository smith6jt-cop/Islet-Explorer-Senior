#!/usr/bin/env python3
"""Build a TissUUmaps ``.tmap`` project for one donor's viewer OME-TIFF.

Companion to ``convert_viewer_ometiff.sh``: that script produces
``/data/islet_ome_tiff/<case_id>.ome.tiff`` (full-res base, LZW, 5-level
pyramid); this one turns a donor slide into something the TissUUmaps server can
open directly, without going through its ``/slide`` auto-detection route.

Why not just point TissUUmaps at the OME-TIFF
---------------------------------------------
The fork's ``/slide`` route *does* detect multi-channel TIFFs and split them,
but two things make it unusable for this cohort (both verified — see
``docs/tissuumaps_evaluation.md``):

1. **Layer names lose marker identity.** ``_extract_multichannel_layers()``
   names every layer ``<basename>_Channel_<i>``, discarding the channel names
   ``_is_multichannel_tiff()`` collected (and that collection is itself broken:
   its ``Name="([^"]+)"`` regex matches the OME-XML ``<Image Name=...>``
   attribute first, so channel 0 comes back as ``Image0``). Positional names
   are actively wrong here because the cohort **mixes panels** — batch-1 slides
   are 59-channel with SST at index 35, batch-2 put SST at index 1 and three of
   them are 58-channel. This is the same trap ``build_channel_config_b64()``
   avoids on the Avivator side by resolving markers BY NAME.
2. **Extraction runs inside the HTTP request.** First view of a 35 GB slide
   would block a Flask worker for as long as the conversion takes.

So we pre-generate the per-channel pyramids offline and emit a ``.tmap`` that
names each layer after its marker, resolved from the slide's own OME-XML.

Usage
-----
    python scripts/senior/build_tissuumaps_project.py \
        --image /data/islet_ome_tiff/6539.ome.tiff \
        --out-dir /data/tissuumaps/projects \
        --regions data/app_data/json/6539.geojson \
        --markers data/app_data/cells/6539_tissue.csv

Idempotent and resumable: an existing channel pyramid is reused unless it is
older than the source slide (or ``--force`` is given), so a killed run can be
restarted. Writing the ``.tmap`` itself is always cheap.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

# Marker -> colour for the four channels that are visible by default, matching
# the Avivator defaults in viewer_helpers.R::build_channel_config_b64().
# TissUUmaps' "Color" filter takes "R,G,B" on a 0-200 scale where 100 is
# neutral, so a plain hex colour maps to 2 * (channel / 255) * 100.
DEFAULT_ON = {
    "INS": (228, 26, 28),      # red
    "GCG": (55, 126, 184),     # blue
    "SST": (255, 204, 0),      # yellow
    "DAPI": (127, 127, 127),   # grey
}

# Slides carry OME-XML channel names verbatim from the qptiff <Biomarker> XML,
# which does not always spell markers the way the app does.
NAME_ALIASES = {
    "PAN-CYTOKERATIN": "Pan_Cytokeratin",
    "PAN CYTOKERATIN": "Pan_Cytokeratin",
    "E-CADHERIN": "ECAD",
    "HLA-DR": "HLADR",
    "CD8": "CD8a",
    "PODOPLANIN": "PDPN",
    "COLLAGEN IV": "CollagenIV",
    "GRANZYME B": "GranzymeB",
}

# Schema version of the emitted project; matches tissuumaps-schema ~=1.3.
SCHEMA_VERSION = "1.3"

# Full-res base of the viewer OME-TIFFs. convert_viewer_ometiff.sh runs
# bioformats2raw with `-s 0 --use-existing-resolutions`, so level 0 is the
# slide's native resolution. NOTE the app's PIXEL_SIZE_UM (0.5078) is the
# *resolution #1* value used for the QuPath GeoJSONs — half this one. Anything
# exported from QuPath at resolution #1 therefore needs --region-scale 2.
DEFAULT_UM_PER_PX = 0.2539


def log(msg: str) -> None:
    print(f"[tmap] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Channel names
# ---------------------------------------------------------------------------


def read_channel_names(image: Path) -> list[str]:
    """Return per-channel marker names from the slide's OME-XML.

    Parses the ``<Channel Name=...>`` elements in document order rather than
    regex-matching the first ``Name="`` in the XML (which is the ``<Image>``
    element's own name attribute). Falls back to ``Channel_<i>`` for any
    channel without a usable name.
    """
    import tifffile

    with tifffile.TiffFile(str(image)) as tif:
        desc = tif.pages[0].description or ""
        n_pages = len(tif.pages)

    names: list[str] = []
    if "<Channel" in desc:
        try:
            root = ET.fromstring(desc)
            # OME-XML is namespaced; match on the local tag name.
            for el in root.iter():
                if el.tag.rsplit("}", 1)[-1] != "Channel":
                    continue
                nm = (el.get("Name") or el.get("Fluor") or "").strip()
                names.append(nm)
        except ET.ParseError:
            log("WARNING: OME-XML did not parse; falling back to regex")
            names = [m.strip() for m in re.findall(r'<Channel[^>]*\sName="([^"]*)"', desc)]

    if len(names) != n_pages:
        if names:
            log(f"WARNING: {len(names)} channel names for {n_pages} pages — padding")
        names = (names + [""] * n_pages)[:n_pages]

    return [normalize_name(nm, i) for i, nm in enumerate(names)]


def normalize_name(name: str, index: int) -> str:
    name = (name or "").strip()
    if not name:
        return f"Channel_{index}"
    return NAME_ALIASES.get(name.upper(), name)


def safe_stem(name: str) -> str:
    """Filesystem-safe channel stem (marker names contain '/', spaces, '-')."""
    return re.sub(r"[^A-Za-z0-9._+-]", "_", name)


# ---------------------------------------------------------------------------
# Channel extraction
# ---------------------------------------------------------------------------


def extract_channel(
    image: Path,
    page: int,
    out_file: Path,
    compression: str,
    quality: int,
    lo_pct: float,
    hi_pct: float,
    force: bool,
) -> bool:
    """Write one channel as a pyramidal 8-bit TIFF. Returns True if it ran."""
    import pyvips

    if out_file.exists() and not force:
        if out_file.stat().st_mtime >= image.stat().st_mtime:
            return False

    src = pyvips.Image.new_from_file(str(image), page=page)

    # Percentile stretch measured on a shrunk copy — scanning a 35 GB slide's
    # full page twice is far too expensive, and percentiles are stable under
    # downsampling.
    shrink = max(1, int(min(src.width, src.height) // 2048))
    stats = src.shrink(shrink, shrink) if shrink > 1 else src
    lo, hi = stats.percent(lo_pct), stats.percent(hi_pct)
    if hi <= lo:
        lo, hi = 0, max(1, int(stats.max()))

    # Clip in float, then scaleimage() to uchar. Do NOT substitute `.cast()`
    # here: on libvips 8.15 the cast path makes JPEG-in-TIFF fail with
    # "VipsJpeg: Unsupported color conversion request" for single-band images,
    # while the scaleimage path (what TissUUmaps itself uses) writes fine.
    # Clipping anchors the range at 0/255 so scaleimage is effectively a no-op
    # and channels stay comparable to each other.
    out = (255.0 * (src - lo)) / (hi - lo)
    out = (out < 0).ifthenelse(0, out)
    out = (out > 255).ifthenelse(255, out)
    out = out.scaleimage()

    kwargs = dict(pyramid=True, tile=True, tile_width=256, tile_height=256,
                  compression=compression, bigtiff=True, properties=True)
    if compression == "jpeg":
        kwargs["Q"] = quality
    out_file.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_file.with_suffix(".tif.part")
    try:
        out.tiffsave(str(tmp), **kwargs)
    except Exception as exc:  # noqa: BLE001 - re-raised with an actionable hint
        tmp.unlink(missing_ok=True)
        if compression == "jpeg" and "color conversion" in str(exc):
            # Greyscale JPEG-in-TIFF is build-dependent: libvips 8.15.1 rejects
            # some single-band pipelines with "VipsJpeg: Unsupported color
            # conversion request". deflate/lzw always work (and are lossless).
            raise RuntimeError(
                f"libvips refused JPEG for channel {page}; re-run with "
                f"--compression deflate (lossless) or lzw"
            ) from exc
        raise
    tmp.replace(out_file)
    return True


# ---------------------------------------------------------------------------
# Regions / markers
# ---------------------------------------------------------------------------


def copy_regions(src: Path, dst: Path, scale: float, layer_index: int = 0) -> int:
    """Copy a QuPath GeoJSON, scaling coordinates into full-res pixel space.

    Also stamps ``properties.collectionIndex`` on every feature. Without it,
    ``regionUtils.geoJSON2regions()`` sees >1 image layer and pops a blocking
    "which layer index?" modal, which silently defeats ``autoLoad`` — verified:
    regions stay at 0 until the prompt is answered. Every channel layer is
    pixel-aligned, so any index works; 0 (the first channel) is the reference.

    QuPath's own ``properties.name`` and ``properties.classification.name`` are
    read natively by TissUUmaps, so islet names and classes carry over as-is.
    """
    with open(src) as fh:
        gj = json.load(fh)

    def walk(coords):
        if coords and isinstance(coords[0], (int, float)):
            return [coords[0] * scale, coords[1] * scale] + list(coords[2:])
        return [walk(c) for c in coords]

    feats = gj.get("features", gj if isinstance(gj, list) else [gj])
    for feat in feats:
        if not isinstance(feat, dict):
            continue
        geom = feat.get("geometry")
        if geom and geom.get("coordinates") is not None and scale != 1.0:
            geom["coordinates"] = walk(geom["coordinates"])
        props = feat.setdefault("properties", {})
        props["collectionIndex"] = layer_index
        if not props.get("name") and feat.get("id"):
            props["name"] = feat["id"]

    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w") as fh:
        json.dump(gj, fh)
    return len(feats)


# ---------------------------------------------------------------------------
# Project assembly
# ---------------------------------------------------------------------------


def build_project(
    case_id: str,
    channels: list[str],
    channel_files: list[str],
    um_per_px: float,
    marker_rel: str | None,
    marker_cols: tuple[str, str, str],
    marker_coord_factor: float,
    region_rel: str | None,
) -> dict:
    layers, opacities, visibilities, filters = [], {}, {}, {}

    # Resolve the default-on markers BY NAME, never by index — the cohort mixes
    # panels (batch-1 SST@35 vs batch-2 SST@1, 58- vs 59-channel).
    upper = {c.upper(): i for i, c in enumerate(channels)}
    for idx, (name, rel) in enumerate(zip(channels, channel_files)):
        layers.append({"name": name, "tileSource": rel})
        opacities[str(idx)] = 1.0
        visibilities[str(idx)] = False

    matched = []
    for marker, (r, g, b) in DEFAULT_ON.items():
        idx = upper.get(marker.upper())
        if idx is None:
            log(f"WARNING: default channel {marker} not found in this panel")
            continue
        matched.append(f"{marker}@{idx}")
        visibilities[str(idx)] = True
        filters[str(idx)] = [
            {"name": "Color", "value": f"{round(r/255*200)},{round(g/255*200)},{round(b/255*200)}"}
        ]
    log(f"default-on channels: {', '.join(matched) if matched else 'NONE'}")

    project: dict = {
        "schemaVersion": SCHEMA_VERSION,
        "filename": case_id,
        "description": f"Islet Explorer — donor {case_id} ({len(channels)}-plex)",
        "layers": layers,
        "layerOpacities": opacities,
        "layerVisibilities": visibilities,
        "layerFilters": filters,
        "filters": ["Color", "Brightness", "Contrast"],
        # Additive blending is what makes overlaid fluorescence channels read
        # correctly; the default "source-over" just paints the top layer.
        "compositeMode": "lighter",
        "mpp": um_per_px,
        "markerFiles": [],
        "regionFiles": [],
        "hideTabs": False,
        "hideNavigator": False,
    }

    if marker_rel:
        x_col, y_col, group_col = marker_cols
        project["markerFiles"].append({
            "title": "Cells",
            "comment": "Single-cell centroids, coloured by phenotype",
            "path": marker_rel,
            "autoLoad": False,
            "expectedHeader": {
                "X": x_col,
                "Y": y_col,
                "gb_col": group_col,
                "coord_factor": marker_coord_factor,
            },
            "expectedRadios": {"cb_gr": True, "cb_gr_rand": True},
        })

    if region_rel:
        project["regionFiles"].append({
            "title": "Islets (QuPath)",
            "path": region_rel,
            "autoLoad": True,
        })

    return project


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--image", required=True, type=Path,
                    help="Donor OME-TIFF, e.g. /data/islet_ome_tiff/6539.ome.tiff")
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="TissUUmaps project root (the server's slide dir)")
    ap.add_argument("--case-id", default=None, help="Defaults to the image stem")
    ap.add_argument("--regions", type=Path, default=None, help="QuPath GeoJSON to attach")
    ap.add_argument("--markers", type=Path, default=None, help="Per-cell CSV to attach")
    ap.add_argument("--marker-cols", default="X_centroid,Y_centroid,phenotype",
                    help="x,y,group columns in the marker CSV")
    ap.add_argument("--um-per-px", type=float, default=DEFAULT_UM_PER_PX,
                    help=f"Slide resolution, µm/px (default {DEFAULT_UM_PER_PX}, the full-res base)")
    ap.add_argument("--region-scale", type=float, default=2.0,
                    help="Multiply GeoJSON coords by this (QuPath resolution #1 -> full res)")
    ap.add_argument("--region-layer", type=int, default=0,
                    help="Layer index regions attach to (suppresses the import prompt)")
    ap.add_argument("--compression", default="deflate", choices=["jpeg", "deflate", "lzw", "none"],
                    help="Channel pyramid codec. deflate (default) is lossless and always works; "
                         "jpeg is much smaller but greyscale JPEG-in-TIFF is libvips-build-dependent")
    ap.add_argument("--quality", type=int, default=90, help="JPEG quality when --compression jpeg")
    ap.add_argument("--percentiles", default="0.5,99.5", help="lo,hi percentile stretch to 8-bit")
    ap.add_argument("--channels", default=None,
                    help="Comma-separated marker names to extract (default: all)")
    ap.add_argument("--force", action="store_true", help="Re-extract channels that already exist")
    ap.add_argument("--dry-run", action="store_true", help="Report the plan, extract nothing")
    args = ap.parse_args()

    image = args.image.resolve()
    if not image.is_file():
        log(f"ERROR: no such image: {image}")
        return 1

    case_id = args.case_id or re.sub(r"\.ome$", "", image.stem, flags=re.I)
    out_dir = args.out_dir.resolve()
    ch_dir = out_dir / "channels" / case_id
    lo_pct, hi_pct = (float(v) for v in args.percentiles.split(","))

    names = read_channel_names(image)
    log(f"{case_id}: {len(names)} channels")
    for m in DEFAULT_ON:
        idx = next((i for i, n in enumerate(names) if n.upper() == m), None)
        log(f"   {m:<5} -> index {idx if idx is not None else 'MISSING'}")

    wanted = None
    if args.channels:
        wanted = {c.strip().upper() for c in args.channels.split(",")}

    keep: list[tuple[int, str]] = [
        (i, n) for i, n in enumerate(names) if wanted is None or n.upper() in wanted
    ]
    if wanted:
        log(f"extracting {len(keep)} of {len(names)} channels (--channels)")

    if args.dry_run:
        log(f"DRY RUN: would write {len(keep)} pyramids to {ch_dir}")
        log(f"DRY RUN: would write {out_dir / (case_id + '.tmap')}")
        return 0

    t0 = time.time()
    channel_files, kept_names = [], []
    for n_done, (page, name) in enumerate(keep, 1):
        stem = f"{page:02d}_{safe_stem(name)}"
        out_file = ch_dir / f"{stem}.tif"
        ran = extract_channel(image, page, out_file, args.compression, args.quality,
                              lo_pct, hi_pct, args.force)
        log(f"  [{n_done}/{len(keep)}] {name:<20} {'extracted' if ran else 'cached'} "
            f"({out_file.stat().st_size / 1e6:.1f} MB)")
        channel_files.append(f"channels/{case_id}/{stem}.tif.dzi")
        kept_names.append(name)

    total = sum(f.stat().st_size for f in ch_dir.glob("*.tif"))
    log(f"channels: {total / 1e9:.2f} GB from a {image.stat().st_size / 1e9:.2f} GB slide "
        f"({total / max(1, image.stat().st_size):.0%}) in {time.time() - t0:.0f}s")

    region_rel = None
    if args.regions:
        region_rel = f"regions/{case_id}.geojson"
        n_feat = copy_regions(args.regions, out_dir / region_rel, args.region_scale,
                              args.region_layer)
        log(f"regions: {n_feat} features scaled x{args.region_scale} -> {region_rel}")

    marker_rel = None
    if args.markers:
        marker_rel = f"markers/{case_id}.csv"
        dst = out_dir / marker_rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.markers, dst)
        log(f"markers: {dst.stat().st_size / 1e6:.1f} MB -> {marker_rel}")

    x_col, y_col, group_col = (c.strip() for c in args.marker_cols.split(","))
    project = build_project(
        case_id, kept_names, channel_files, args.um_per_px,
        marker_rel, (x_col, y_col, group_col),
        # Marker CSVs are in µm; TissUUmaps wants image pixels.
        marker_coord_factor=round(1.0 / args.um_per_px, 6),
        region_rel=region_rel,
    )

    tmap_path = out_dir / f"{case_id}.tmap"
    with open(tmap_path, "w") as fh:
        json.dump(project, fh, indent=2)
    log(f"wrote {tmap_path}")

    try:
        from tissuumaps_schema.utils import CURRENT_SCHEMA_MODULE as schema

        if schema.VERSION != SCHEMA_VERSION:
            log(f"NOTE: installed schema is {schema.VERSION}, project declares {SCHEMA_VERSION}")
        schema.project.Project.model_validate(project)
        log(f"schema {schema.VERSION}: valid")
    except ImportError:
        log("tissuumaps-schema not installed — skipped validation")
    except Exception as exc:  # noqa: BLE001 - validation errors are informational
        log(f"WARNING: schema validation failed: {exc}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
