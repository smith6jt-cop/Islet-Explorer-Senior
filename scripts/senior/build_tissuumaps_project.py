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
        --regions data/Islets \
        --markers data/app_data/cells/6539_tissue.csv

``--regions`` takes either a GeoJSON file or the directory the exports live in
(``data/Islets``), in which case the donor's file is found by case id; ``.gz`` is
read transparently. Region coordinates are rescaled into full-resolution pixels
by ``--region-scale``, which defaults to deriving the factor from the two pixel
sizes and then checking it against the image bounds.

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

# Colours for the four channels visible by default, following the Avivator
# palette in viewer_helpers.R::build_channel_config_b64().
#
# TissUUmaps' "Color" filter values run **0-100 per channel, where 100 leaves
# that channel untouched** -- filterUtils applies `channels(value - 100)`, so
# every value is a subtraction and 0 means "remove this channel entirely".
# Getting this scale wrong is not a subtle mis-tint: the UI renders the value
# through `componentToHex(c) = floor(c * 255 / 100)`, so anything above 100
# overflows to a 3-digit hex component, the resulting string is not a valid
# colour, the picker falls back to #000000, and the layer reads back as
# "0,0,0" -> channels(-100) -> **solid black and invisible**.
#
# Pure primaries are used rather than approximations of specific hex colours:
# they isolate a fluorescence channel cleanly under additive blending, and they
# round-trip through the picker exactly (100 -> "ff" -> 100).
DEFAULT_ON = {
    "INS": "100,0,0",      # red
    "GCG": "0,0,100",      # blue
    "SST": "100,100,0",    # yellow
    "DAPI": "50,50,50",    # dim grey, so nuclei don't wash out the markers
}

# Every other channel gets this explicitly. 100,100,100 is a no-op, i.e. plain
# greyscale -- the point is that no layer is ever left without a Color entry,
# so switching one on can never land on a stale or black picker value.
NEUTRAL_COLOR = "100,100,100"

# Filter controls offered per layer, and their neutral starting values.
#
# NB "Min"/"Max" -- 0-255 black/white-point sliders, the closest thing
# TissUUmaps has to Avivator's contrast limits -- are deliberately NOT used.
# filterUtils.js implements them, but they are missing from the Filter enum in
# tissuumaps-schema 1.3, and views.py validates every .tmap on load: including
# them puts a pydantic traceback in a warning banner over the viewer every time
# the project is opened. They also cannot express a real window, because each
# reads the other's *global* default rather than the sibling layer value
# (filterUtils.getFilterParams returns _filters[name].params).
#
# Everything starts neutral on purpose: the contrast decision that matters is
# the percentile stretch applied during extraction on the full-resolution
# 16-bit page, where there is far more signal than in these 8-bit proxies.
FILTER_DEFAULTS = [
    ("Color", None),        # filled in per layer
    ("Brightness", 0),      # slider range -50..50, 0 = no change
    ("Contrast", 1),        # slider range 0..8,    1 = no change
]
FILTERS = [name for name, _ in FILTER_DEFAULTS]

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
# *resolution #1* value used for the QuPath GeoJSONs — twice this one, which is
# where the default region scale of x2 comes from (see DEFAULT_REGION_UM_PER_PX).
DEFAULT_UM_PER_PX = 0.2539

# Pixel size the QuPath structure GeoJSONs are exported at. Same constant as the
# app's PIXEL_SIZE_UM (00_globals.R), which was measured, not assumed: 0.5078
# µm/px, the median over 2,355 islets across 3 donors. Dividing it by
# DEFAULT_UM_PER_PX is where the default x2 region scale comes from, rather than
# a bare "2" nobody can check.
DEFAULT_REGION_UM_PER_PX = 0.5078


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


def read_image_size(image: Path) -> tuple[int, int]:
    """Width and height of the full-resolution page, in pixels."""
    import tifffile

    with tifffile.TiffFile(str(image)) as tif:
        page = tif.pages[0]
        return int(page.imagewidth), int(page.imagelength)


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


def _histogram(src):
    """Full-resolution histogram of an integer image, as int64 counts."""
    import numpy as np

    hist = src.hist_find()
    dtype = {"uchar": np.uint8, "ushort": np.uint16, "uint": np.uint32,
             "float": np.float32, "double": np.float64}.get(hist.format, np.uint32)
    return np.ndarray(buffer=hist.write_to_memory(), dtype=dtype,
                      shape=(hist.width,)).astype("int64")


def channel_stats(src, stretch: str, lo_pct: float, hi_pct: float,
                  bg_k: float, signal_q: float) -> dict:
    """Choose display cut points for one channel, from a single histogram pass.

    Marker intensity is not a biological readout in this cohort -- essentially
    all of its variation is technical (saturation, background level, and markers
    whose positive fraction ranges over orders of magnitude). The viewer
    therefore has to stay legible across those, which a fixed percentile pair
    cannot do: one pair cannot serve a marker present in 0.05% of pixels and one
    present in 40%.

    ``robust`` (default) anchors to the background peak instead:

      * black point = background mode + ``bg_k`` * sigma, sigma estimated from a
        one-sided MAD *below* the mode, since signal only contaminates the upper
        side. Background therefore renders black whatever level it sits at.
      * white point = the ``signal_q`` quantile of pixels above the black point,
        with the saturation bin excluded so clipped pixels cannot drag it. This
        adapts to however rare the marker is.
      * a channel with no population above background is left dark rather than
        having its noise stretched to full range, so an absent marker looks
        absent instead of looking like snow.

    Measured against ground-truth masks on synthetic channels (rendered 0-255,
    background mean -> signal mean, and how much signal saturates):

        channel      percentile 0.5/99.5        robust
        dense        bg 3  sig 160  clip  1.5%  bg 0  sig 153  clip 1.0%
        sparse       bg 127 sig 255 clip 100%   bg 0  sig 177  clip 3.7%
        high bg      bg 23 sig 214  clip 25.4%  bg 0  sig 154  clip 1.3%
        saturated    bg 1  sig 151  clip 25.2%  bg 0  sig 121  clip 8.9%
        absent       bg 127 (noise as signal)   bg 0  (stays black)

    ``percentile`` keeps the old fixed-percentile behaviour for comparison.
    """
    import numpy as np

    counts = _histogram(src)
    total = int(counts.sum())
    sat_bin = counts.size - 1
    if total <= 0:
        return {"lo": 0, "hi": 1, "clip_lo": 0.0, "clip_hi": 0.0,
                "sat": 0.0, "signal": 0.0, "mode": 0, "sigma": 0.0, "dark": True}

    if stretch == "percentile":
        cum = np.cumsum(counts)
        lo = int(np.searchsorted(cum, total * lo_pct / 100.0, side="left"))
        hi = int(np.searchsorted(cum, total * hi_pct / 100.0, side="left"))
        mode, sigma, dark = lo, 0.0, False
    else:
        # Detector saturation piles every clipped pixel into the top bin, where
        # it can out-count the background peak and be mistaken for it, so the
        # mode is searched strictly below saturation.
        mode = int(np.argmax(counts[:sat_bin])) if sat_bin > 0 else 0
        below = counts[: mode + 1]
        if below.sum() > 0:
            mad = int(np.searchsorted(np.cumsum(below[::-1]), below.sum() * 0.5))
        else:
            mad = 1
        sigma = max(1.0, 1.4826 * mad)
        lo = int(min(mode + bg_k * sigma, sat_bin))

        signal = counts[lo + 1: sat_bin]
        n_signal = int(signal.sum())
        dark = n_signal < total * 1e-4
        if dark:
            hi = int(min(lo + 6.0 * sigma, sat_bin))
        else:
            hi = lo + 1 + int(np.searchsorted(np.cumsum(signal),
                                              n_signal * signal_q / 100.0))
            hi = min(max(hi, lo + int(6.0 * sigma)), sat_bin)

    if hi <= lo:
        lo, hi = 0, max(1, sat_bin)

    return {
        "lo": lo, "hi": hi, "mode": mode, "sigma": sigma, "dark": dark,
        "clip_lo": float(counts[:lo].sum()) / total * 100.0,
        "clip_hi": float(counts[hi + 1:].sum()) / total * 100.0,
        "sat": float(counts[sat_bin]) / total * 100.0,
        "signal": float(counts[lo + 1: sat_bin].sum()) / total * 100.0,
    }


def linear_stretch_lut(lo: int, hi: int, fmt: str):
    """LUT mapping the source's integer range onto 0-255 over ``[lo, hi]``.

    Rounds to nearest. vips casts float->int by truncation, which darkens every
    pixel by up to a level; so does scikit-image's
    ``exposure.rescale_intensity(..., out_range=np.uint8)`` (measured mean error
    0.50 grey levels against 0.25 here), which is why neither is used directly.
    """
    import numpy as np
    import pyvips

    n_bins = {"uchar": 256, "ushort": 65536}.get(fmt)
    if n_bins is None:
        raise RuntimeError(f"maplut needs an integer source format, got {fmt!r}")

    ramp = np.arange(n_bins, dtype=np.float64)
    vals = np.clip(np.rint((255.0 * (ramp - lo)) / (hi - lo)), 0, 255).astype(np.uint8)
    return pyvips.Image.new_from_memory(vals.tobytes(), n_bins, 1, 1, "uchar")


def extract_channel(
    image: Path,
    page: int,
    out_file: Path,
    compression: str,
    quality: int,
    stretch: str,
    lo_pct: float,
    hi_pct: float,
    bg_k: float,
    signal_q: float,
    force: bool,
) -> dict | None:
    """Write one channel as a pyramidal 8-bit TIFF.

    Returns the stretch it applied, or None if a current file was reused.
    """
    import pyvips

    if out_file.exists() and not force:
        if out_file.stat().st_mtime >= image.stat().st_mtime:
            return None

    src = pyvips.Image.new_from_file(str(image), page=page)

    # Percentiles are measured on the FULL-RESOLUTION page.
    #
    # An earlier version measured them on a shrunk copy, justified as
    # "percentiles are stable under downsampling". They are not: shrink()
    # block-averages, which flattens precisely the sparse bright pixels the
    # upper percentile exists to find. Measured on a 4096x4096 channel with
    # sparse puncta, at the shrink factor a 40k-pixel slide would get:
    #
    #     shrink=1   99.5th pct = 9243   (true)
    #     shrink=4   99.5th pct = 4797
    #     shrink=14  99.5th pct = 1453   -- 6.4x too low
    #
    # Every feature above ~1453 then saturated to flat white with no internal
    # structure. The shrink also bought nothing: vips_shrink reads every input
    # pixel regardless, so it saved no I/O, only histogram bins. Measured cost
    # of percent() on the full page vs the shrunk copy: 10 ms vs 11 ms.
    # One histogram pass yields both cut points AND how much actually clips, so
    # the cost is lower than two percent() calls (each of which builds its own
    # histogram) and the quality of the stretch is visible per channel instead
    # of being taken on trust.
    st = channel_stats(src, stretch, lo_pct, hi_pct, bg_k, signal_q)
    lo, hi = st["lo"], st["hi"]

    if compression != "jpeg":
        # Precompute the whole tone curve as a lookup table and let vips apply
        # it with maplut: one table lookup per pixel instead of a subtract,
        # multiply, two compares and a cast. Measured 11 ms vs 40 ms on a
        # 2048x2048 channel, bit-identical output.
        #
        # The real reason to prefer it is that the mapping is now ordinary
        # numpy over a 65536-entry ramp, so any tone curve -- gamma, a
        # histogram-matched CDF, anything scikit-image can express -- can be
        # dropped in without the slide ever entering RAM. skimage's own
        # exposure.rescale_intensity() was measured against this and is NOT
        # used: its uint8 conversion truncates (mean error 0.50 grey levels vs
        # 0.25 for round-to-nearest below).
        out = src.maplut(linear_stretch_lut(lo, hi, src.format))
    else:
        out = (255.0 * (src - lo)) / (hi - lo)
        out = (out < 0).ifthenelse(0, out)
        out = (out > 255).ifthenelse(255, out)
        # scaleimage() is only used where libvips forces it: the JPEG-in-TIFF
        # writer rejects both the maplut and cast paths. It is not equivalent --
        # it RE-STRETCHES by the clipped image's own min/max, so a channel that
        # does not reach both ends of the range is silently remapped and is no
        # longer comparable to its neighbours. A dim test channel that should
        # have occupied 8-bit levels 10..14 came out as 0..255, an error of 241
        # levels. Another reason --compression jpeg is not the default.
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
    return st


# ---------------------------------------------------------------------------
# Regions / markers
# ---------------------------------------------------------------------------


def resolve_regions(spec: Path, case_id: str) -> Path | None:
    """Accept a GeoJSON file or a directory to search for the donor's file.

    Region exports live outside this repo and are named per-site, so pointing
    at the directory (``--regions data/Islets``) and letting the case id do the
    matching avoids having to encode somebody's filename convention here.
    Handles ``.gz`` as well; QuPath exports are often stored compressed.
    """
    if spec.is_file():
        return spec
    if not spec.is_dir():
        log(f"WARNING: --regions path does not exist: {spec}")
        return None

    hits = sorted(
        f for pattern in (f"*{case_id}*.geojson", f"*{case_id}*.geojson.gz",
                          f"*{case_id}*.json", f"*{case_id}*.json.gz")
        for f in spec.glob(pattern)
    )
    if not hits:
        log(f"WARNING: no region file matching '*{case_id}*' in {spec}")
        return None
    if len(hits) > 1:
        log(f"WARNING: {len(hits)} region files match '*{case_id}*'; using {hits[0].name}")
        for h in hits[1:]:
            log(f"           also matched: {h.name}")
    return hits[0]


def _load_geojson(src: Path):
    """Read a GeoJSON, transparently handling gzip."""
    if src.suffix == ".gz":
        import gzip

        with gzip.open(src, "rt") as fh:
            return json.load(fh)
    with open(src) as fh:
        return json.load(fh)


def _geometry_bounds(feats) -> tuple[float, float] | None:
    """Max x and y over every coordinate, or None if there is no geometry."""
    max_x = max_y = float("-inf")

    def walk(coords):
        nonlocal max_x, max_y
        if coords and isinstance(coords[0], (int, float)):
            max_x, max_y = max(max_x, coords[0]), max(max_y, coords[1])
            return
        for c in coords:
            walk(c)

    for feat in feats:
        geom = feat.get("geometry") if isinstance(feat, dict) else None
        if geom and geom.get("coordinates") is not None:
            walk(geom["coordinates"])
    return None if max_x == float("-inf") else (max_x, max_y)


def infer_region_scale(feats, img_w: int, img_h: int, nominal: float) -> float:
    """Settle the factor taking region coordinates into full-res image pixels.

    ``nominal`` comes from the two pixel sizes (the GeoJSON's µm/px over the
    slide's), which is the only *physical* statement available and so is what we
    trust. The geometry is used to falsify it, not to derive it: a region can
    never legitimately fall outside the image, so a nominal factor that pushes
    one past the edge is provably wrong and we back off to the largest smaller
    power of two that does fit.

    The reverse error -- coordinates that are already full-res being doubled --
    is only detectable when the tissue covers more than about half the slide,
    which it usually does here. The chosen factor and the resulting coverage are
    both logged so a survivor is at least visible; ``--region-scale`` overrides.
    """
    bounds = _geometry_bounds(feats)
    if bounds is None or img_w <= 0 or img_h <= 0:
        return nominal
    max_x, max_y = bounds

    def coverage(c: float) -> float:
        return max(max_x * c / img_w, max_y * c / img_h)

    # 2% slack: QuPath polygons trace the tissue edge and can sit a pixel or two
    # past the level boundary after rounding.
    if coverage(nominal) <= 1.02:
        cov = coverage(nominal)
        if cov < 0.25:
            log(f"WARNING: at x{nominal:g} the regions span only {cov:.0%} of the "
                f"slide. That is plausible for a small section but also what a "
                f"too-small factor looks like -- check the overlay.")
        return nominal

    feasible = [c for c in (8.0, 4.0, 2.0, 1.0, 0.5, 0.25)
                if c < nominal and coverage(c) <= 1.02]
    chosen = feasible[0] if feasible else 1.0
    log(f"WARNING: region extent {max_x:.0f}x{max_y:.0f} scaled by the expected "
        f"x{nominal:g} would fall outside the {img_w}x{img_h} image "
        f"({coverage(nominal):.2f}x its size), so that factor is wrong. "
        f"Using x{chosen:g} instead; pass --region-scale to override.")
    return chosen


def region_features(gj) -> list:
    """The feature list of a GeoJSON, whatever container QuPath wrapped it in.

    QuPath writes a FeatureCollection when several objects are exported and a
    bare list (or a single Feature) when one is.
    """
    if isinstance(gj, dict):
        return gj.get("features", [gj])
    return gj if isinstance(gj, list) else [gj]


def copy_regions(gj, dst: Path, scale: float, layer_index: int = 0) -> int:
    """Write a QuPath GeoJSON out, scaling coordinates into full-res pixel space.

    Also stamps ``properties.collectionIndex`` on every feature. Without it,
    ``regionUtils.geoJSON2regions()`` sees >1 image layer and pops a blocking
    "which layer index?" modal, which silently defeats ``autoLoad`` — verified:
    regions stay at 0 until the prompt is answered. Every channel layer is
    pixel-aligned, so any index works; 0 (the first channel) is the reference.

    QuPath's own ``properties.name`` and ``properties.classification.name`` are
    read natively by TissUUmaps, so islet names and classes carry over as-is.
    """
    def walk(coords):
        if coords and isinstance(coords[0], (int, float)):
            return [coords[0] * scale, coords[1] * scale] + list(coords[2:])
        return [walk(c) for c in coords]

    feats = region_features(gj)
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
        # Explicit neutral start for EVERY layer, not just the visible ones: a
        # channel switched on later must never come up black or carry a value
        # left over from another layer's picker.
        filters[str(idx)] = [
            {"name": name, "value": NEUTRAL_COLOR if value is None else value}
            for name, value in FILTER_DEFAULTS
        ]

    matched = []
    for marker, color in DEFAULT_ON.items():
        idx = upper.get(marker.upper())
        if idx is None:
            log(f"WARNING: default channel {marker} not found in this panel")
            continue
        matched.append(f"{marker}@{idx}")
        visibilities[str(idx)] = True
        filters[str(idx)][0]["value"] = color
    log(f"default-on channels: {', '.join(matched) if matched else 'NONE'}")

    project: dict = {
        "schemaVersion": SCHEMA_VERSION,
        "filename": case_id,
        "description": f"Islet Explorer — donor {case_id} ({len(channels)}-plex)",
        "layers": layers,
        "layerOpacities": opacities,
        "layerVisibilities": visibilities,
        "layerFilters": filters,
        "filters": FILTERS,
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
    ap.add_argument("--regions", type=Path, default=None,
                    help="QuPath GeoJSON to attach, or a directory to find the donor's "
                         "export in (e.g. data/Islets); .gz is read transparently")
    ap.add_argument("--markers", type=Path, default=None, help="Per-cell CSV to attach")
    ap.add_argument("--marker-cols", default="X_centroid,Y_centroid,phenotype",
                    help="x,y,group columns in the marker CSV")
    ap.add_argument("--um-per-px", type=float, default=DEFAULT_UM_PER_PX,
                    help=f"Slide resolution, µm/px (default {DEFAULT_UM_PER_PX}, the full-res base)")
    ap.add_argument("--region-scale", default="auto",
                    help="Multiply GeoJSON coords by this to reach full-res pixels. "
                         "'auto' (default) derives it from the two pixel sizes and "
                         "checks the result against the image bounds")
    ap.add_argument("--region-um-per-px", type=float, default=DEFAULT_REGION_UM_PER_PX,
                    help=f"Pixel size the GeoJSON was exported at, µm/px "
                         f"(default {DEFAULT_REGION_UM_PER_PX}, QuPath resolution #1)")
    ap.add_argument("--region-layer", type=int, default=0,
                    help="Layer index regions attach to (suppresses the import prompt)")
    ap.add_argument("--compression", default="deflate", choices=["jpeg", "deflate", "lzw", "none"],
                    help="Channel pyramid codec. deflate (default) is lossless and always works; "
                         "jpeg is much smaller but greyscale JPEG-in-TIFF is libvips-build-dependent")
    ap.add_argument("--quality", type=int, default=90, help="JPEG quality when --compression jpeg")
    ap.add_argument("--stretch", default="robust", choices=["robust", "percentile"],
                    help="robust (default) anchors black to the background peak and white to a "
                         "quantile of the above-background population, so the result survives "
                         "high background, saturation and markers of any rarity; percentile "
                         "applies fixed cut points to the whole image")
    ap.add_argument("--percentiles", default="0.5,99.5",
                    help="lo,hi cut points for --stretch percentile")
    ap.add_argument("--background-k", type=float, default=3.0,
                    help="black point = background mode + K sigma (--stretch robust)")
    ap.add_argument("--signal-quantile", type=float, default=99.0,
                    help="white point = this quantile of above-background pixels (--stretch robust)")
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
        stats = extract_channel(image, page, out_file, args.compression, args.quality,
                                args.stretch, lo_pct, hi_pct, args.background_k,
                                args.signal_quantile, args.force)
        size_mb = out_file.stat().st_size / 1e6
        if stats is None:
            log(f"  [{n_done}/{len(keep)}] {name:<20} cached    ({size_mb:.1f} MB)")
        else:
            flags = []
            if stats.get("dark"):
                flags.append("no signal above background - left dark")
            if stats.get("sat", 0) > 0.1:
                flags.append(f"{stats['sat']:.2f}% saturated")
            if stats["clip_hi"] > 2.0:
                flags.append(f"{stats['clip_hi']:.1f}% clipped white")
            note = ("  <-- " + "; ".join(flags)) if flags else ""
            log(f"  [{n_done}/{len(keep)}] {name:<20} extracted ({size_mb:.1f} MB) "
                f"bg={stats['mode']}+/-{stats['sigma']:.0f} -> 16-bit {stats['lo']}-{stats['hi']}, "
                f"signal {stats['signal']:.3f}% of px{note}")
        channel_files.append(f"channels/{case_id}/{stem}.tif.dzi")
        kept_names.append(name)

    # Report per-channel cost and project the full panel. Reporting only
    # "X GB from a Y GB slide" on a --channels run reads as a completeness
    # claim when it is a fraction of the panel, which is exactly the wrong
    # number to size storage from.
    on_disk = [f.stat().st_size for f in ch_dir.glob("*.tif")]
    total, n_files, slide = sum(on_disk), len(on_disk), image.stat().st_size
    per_channel = total / max(1, n_files)
    projected = per_channel * len(names)
    log(f"channels: {total / 1e9:.2f} GB across {n_files} of {len(names)} channels "
        f"({per_channel / 1e6:.0f} MB each) in {time.time() - t0:.0f}s")
    log(f"   full panel projects to {projected / 1e9:.1f} GB per donor, "
        f"{projected / max(1, slide):.0%} of the {slide / 1e9:.1f} GB slide"
        + ("" if n_files == len(names) else "  [extrapolated from a partial run]"))

    region_rel = None
    if args.regions:
        src = resolve_regions(args.regions.resolve(), case_id)
        if src is not None:
            gj = _load_geojson(src)
            feats = region_features(gj)
            if args.region_scale == "auto":
                img_w, img_h = read_image_size(image)
                nominal = args.region_um_per_px / args.um_per_px
                scale = infer_region_scale(feats, img_w, img_h, nominal)
                log(f"regions: {src.name} -> x{scale:g} "
                    f"({args.region_um_per_px} / {args.um_per_px} µm/px, "
                    f"checked against {img_w}x{img_h})")
            else:
                scale = float(args.region_scale)
            region_rel = f"regions/{case_id}.geojson"
            n_feat = copy_regions(gj, out_dir / region_rel, scale, args.region_layer)
            log(f"regions: {n_feat} features scaled x{scale:g} -> {region_rel}")

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
