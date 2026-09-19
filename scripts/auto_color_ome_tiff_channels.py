#!/usr/bin/env python3
"""
Auto-color OME-TIFF channels in-place using existing OME-XML channel names.

- Reads OME-XML from the TIFF ImageDescription and preserves channel names.
- If a Channel lacks a Color attribute (or is 0/white), assign a distinct color
  from a categorical palette and write it back as a 24-bit RGB decimal (RRGGBB).
- Pixel data is untouched.

Usage:
        python scripts/auto_color_ome_tiff_channels.py \
            --data-root local_images [--glob "*.ome.tif*"] [--dry-run] [--report-only] \
            [--fallback-root /path/to/original_tifs] [--qupath-root /path/to/qupath] \
            [--channel-names-file shiny/Channel_names]

NOTE 2026-09-18: there is no `shiny/` directory in this repository -- the app lives at
`apps/islet/` and the channel-name list is `apps/islet/Channel_names`. The `--channel-names-file`
DEFAULT below still reads `shiny/Channel_names`, so it silently resolves to nothing unless the
path is passed explicitly. Recorded as a finding; the default is not changed here.

Dependencies: tifftools, lxml
"""
import argparse
import glob
import os
import sys
from typing import List
import re
import json

try:
    import tifftools as tt
except Exception:
    print("ERROR: tifftools is required: pip install tifftools", file=sys.stderr)
    sys.exit(1)

try:
    from lxml import etree
except Exception:
    print("ERROR: lxml is required: pip install lxml", file=sys.stderr)
    sys.exit(1)


PALETTE_HEX = [
    "#FF6B6B",  # red
    "#4D96FF",  # blue
    "#6BCB77",  # green
    "#FFD93D",  # yellow
    "#B983FF",  # purple
    "#FF8E3C",  # orange
    "#2EC4B6",  # teal
    "#E76F51",  # burnt orange
    "#118AB2",  # ocean
    "#EF476F",  # pink/red
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data-root", required=True, help="Folder containing OME-TIFF files (searched recursively)")
    p.add_argument("--glob", default="*.ome.tif*", help="Filename glob to match (default: *.ome.tif*)")
    p.add_argument("--dry-run", action="store_true", help="Print planned changes without writing")
    p.add_argument("--fallback-root", help="Root path to search for original .tif/.tiff (to infer channel names)")
    p.add_argument("--qupath-root", help="Root path to search QuPath sidecar files (CSV/TSV/TXT/JSON/qpproj)")
    p.add_argument("--channel-names-file", default="shiny/Channel_names", help="Optional mapping file with lines like 'NAME (C#)'")
    p.add_argument("--report-only", action="store_true", help="Only report channel names/colors; no writes")
    return p.parse_args()


essential_ns = "http://www.openmicroscopy.org/Schemas/OME/2016-06"


def ensure_rgb_decimal(color_hex: str) -> str:
    color_hex = color_hex.strip()
    if not (color_hex.startswith("#") and len(color_hex) == 7):
        return ""
    try:
        r = int(color_hex[1:3], 16)
        g = int(color_hex[3:5], 16)
        b = int(color_hex[5:7], 16)
        return str((r << 16) | (g << 8) | b)
    except Exception:
        return ""


def make_palette(n: int) -> List[str]:
    # Repeat palette as needed; will cycle colors if channels > palette length
    out = []
    base = PALETTE_HEX
    for i in range(n):
        out.append(base[i % len(base)])
    return out

# Preferred colors for key channels
PREFERRED_HEX = {
    "DAPI": "#9CA3AF",  # grey
    "INS":  "#EF4444",  # red
    "GCG":  "#3B82F6",  # blue
    "SST":  "#F59E0B",  # yellow
}

def load_channel_names_file(path: str):
    mapping = {}
    if not path:
        return mapping
    # Allow relative to project root
    if not os.path.isabs(path):
        cand = os.path.join(os.getcwd(), path)
        if os.path.exists(cand):
            path = cand
    if not os.path.exists(path):
        return mapping
    rx = re.compile(r"^\s*(.*?)\s*\(C(\d+)\)\s*$")
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                m = rx.match(line)
                if not m:
                    continue
                name = m.group(1).strip()
                try:
                    cnum = int(m.group(2))
                except Exception:
                    continue
                idx0 = cnum - 1
                mapping[idx0] = name
    except Exception:
        return {}
    return mapping


def read_ome_xml(info):
    ifds = info.get("ifds", [])
    if not ifds:
        return None, None
    desc_tag = 270
    for ifd in ifds:
        tags = ifd.get("tags", {})
        if desc_tag in tags:
            val = tags[desc_tag].get("data")
            data = val if isinstance(val, (bytes, bytearray)) else str(val).encode("utf-8")
            if b"OME" in data and b"Pixels" in data:
                return ifd, tags[desc_tag]
    return None, None


def parse_channel_info(xml_bytes):
    # Ensure bytes for lxml when XML has an encoding declaration
    if isinstance(xml_bytes, str):
        xml_bytes = xml_bytes.encode("utf-8", errors="ignore")
    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.fromstring(xml_bytes, parser=parser)
    ns = {"ome": root.nsmap.get(None) or essential_ns}
    info = []
    # Iterate images and their channels
    for image in root.xpath("//ome:Image", namespaces=ns):
        pixels = image.find("{*}Pixels")
        if pixels is None:
            continue
        channels = [c for c in pixels.findall("{*}Channel")]
        for idx, ch in enumerate(channels):
            info.append({
                "image_name": image.get("Name") or image.get("ID") or "Image",
                "index": idx,
                "channel_id": ch.get("ID"),
                "channel_name": ch.get("Name"),
                "color": ch.get("Color"),
            })
    return root, info


def infer_name_from_id(ch_id: str) -> str:
    if not ch_id:
        return ""
    # Try patterns like "Channel:Insulin" or "Image:0:Channel:Insulin"
    # Heuristic: take last token after ':' or '/'
    token = ch_id.split(":")[-1].split("/")[-1]
    # Drop common prefixes
    for pref in ("Channel", "channel", "Ch", "ch"):
        if token.startswith(pref):
            token = token[len(pref):].lstrip(" _-#:")
    return token

def extract_annotations(root):
    """
    Extract channel names/colors from StructuredAnnotations.
    Returns: (names: dict[int,str], colors_hex: dict[int,str])
    Supported key patterns (case-insensitive):
      - Channel:0:Name, Channel:0:Color
      - Channel0Name, Channel0Color
      - Channel_0_Name, Channel_0_Color
    """
    ns = {"ome": root.nsmap.get(None) or essential_ns}
    names = {}
    colors = {}
    # MapAnnotation key-value pairs
    for m in root.xpath(".//ome:StructuredAnnotations//ome:MapAnnotation//ome:Value//*[@K or @k or @V or @v]", namespaces=ns):
        k = m.get("K") or m.get("k") or ""
        v = m.get("V") or m.get("v") or (m.text or "")
        if not k:
            continue
        k_low = k.lower()
        # Try to find channel index
        idx = None
        for pat in [r"channel[:_ ]?(\d+)[:_ ]?name", r"channel[:_ ]?(\d+)[:_ ]?color", r"channel(\d+)(name|color)"]:
            mobj = re.search(pat, k_low)
            if mobj:
                idx = int(mobj.group(1))
                break
        if idx is None:
            continue
        if "name" in k_low:
            names[idx] = v.strip()
        elif "color" in k_low:
            colors[idx] = v.strip()
    # XMLAnnotation OriginalMetadata style
    for omd in root.xpath(".//ome:StructuredAnnotations//ome:XMLAnnotation//ome:Value//*[@Key or *//Key]", namespaces=ns):
        # Expect structure like <OriginalMetadata><Key>Channel0Name</Key><Value>Insulin</Value></OriginalMetadata>
        for node in omd.findall(".//*"):
            if node.tag.endswith("OriginalMetadata"):
                key_el = node.find(".//*[(local-name(.)='Key')]")
                val_el = node.find(".//*[(local-name(.)='Value')]")
                if key_el is None or val_el is None:
                    continue
                k = (key_el.text or "").strip()
                v = (val_el.text or "").strip()
                if not k:
                    continue
                k_low = k.lower()
                idx = None
                for pat in [r"channel[:_ ]?(\d+)[:_ ]?name", r"channel[:_ ]?(\d+)[:_ ]?color", r"channel(\d+)(name|color)"]:
                    mobj = re.search(pat, k_low)
                    if mobj:
                        idx = int(mobj.group(1))
                        break
                if idx is None:
                    continue
                if "name" in k_low:
                    names[idx] = v
                elif "color" in k_low:
                    colors[idx] = v
    return names, colors

def _read_text_file(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return f.read()
    except Exception:
        try:
            with open(path, "r", encoding="latin-1") as f:
                return f.read()
        except Exception:
            return ""

def _try_parse_names_from_text(txt: str):
    names = []
    if not txt:
        return names
    # JSON object/array
    try:
        data = json.loads(txt)
        if isinstance(data, dict):
            for key in ("channels", "channel_names", "Channels", "ChannelNames"):
                if key in data and isinstance(data[key], list):
                    names = [str(x) for x in data[key]]
                    if names:
                        return names
        elif isinstance(data, list) and all(isinstance(x, str) for x in data) and 0 < len(data) <= 100:
            return [str(x) for x in data]
    except Exception:
        pass
    # key=value style like Channel0Name=Insulin
    m = re.findall(r"channel[_\s:]?(\d+)[_\s:]?name\s*[:=]\s*(.+)", txt, flags=re.I)
    if m:
        tmp = {}
        for idx, val in m:
            tmp[int(idx)] = val.strip()
        max_idx = max(tmp.keys())
        return [tmp.get(i, "") for i in range(max_idx + 1)]
    # one-per-line
    lines = [x.strip() for x in txt.splitlines() if x.strip()]
    if 0 < len(lines) <= 100 and all(len(x) < 100 for x in lines):
        return lines
    # csv or semicolon
    if "," in txt or ";" in txt:
        parts = [p.strip() for p in re.split(r"[,;]", txt) if p.strip()]
        if 0 < len(parts) <= 100 and all(len(x) < 100 for x in parts):
            return parts
    return []

def _find_original_tif(stem: str, root: str):
    if not root:
        return None
    candidates = []
    for pat in (f"**/{stem}.tif", f"**/{stem}.tiff", f"**/{stem}_*.tif", f"**/{stem}_*.tiff"):
        candidates.extend(glob.glob(os.path.join(root, pat), recursive=True))
    for c in candidates:
        base = os.path.basename(c)
        if os.path.splitext(base)[0] == stem:
            return c
    return candidates[0] if candidates else None

def _extract_names_from_tif(path: str):
    try:
        info = tt.read_tiff(path)
    except Exception:
        return []
    desc_tag = 270
    texts = []
    for ifd in info.get("ifds", []):
        tags = ifd.get("tags", {})
        if desc_tag in tags:
            val = tags[desc_tag].get("data")
            if isinstance(val, (bytes, bytearray)):
                texts.append(val.decode("utf-8", errors="ignore"))
            else:
                texts.append(str(val))
    for txt in texts:
        names = _try_parse_names_from_text(txt)
        if names:
            return names
    return []

def _find_qupath_sidecars(stem: str, root: str):
    if not root:
        return []
    cand = []
    for ext in ("*.csv", "*.tsv", "*.txt", "*.json", "*.qpproj"):
        cand.extend(glob.glob(os.path.join(root, f"**/*{stem}*{ext}"), recursive=True))
    return cand

def _extract_names_from_qupath(files: List[str]):
    for f in files:
        txt = _read_text_file(f)
        names = _try_parse_names_from_text(txt)
        if names:
            return names
    return []


def patch_names_and_colors(xml_bytes, ome_stem: str = "", fallback_root: str = None, qupath_root: str = None, names_map: dict | None = None) -> bytes:
    # Ensure bytes internally
    if isinstance(xml_bytes, str):
        xml_bytes = xml_bytes.encode("utf-8", errors="ignore")
    root, _ = parse_channel_info(xml_bytes)
    ns = {"ome": root.nsmap.get(None) or essential_ns}
    ann_names, ann_colors = extract_annotations(root)
    names_map = names_map or {}
    fallback_names = []
    if ome_stem:
        if fallback_root:
            tifp = _find_original_tif(ome_stem, fallback_root)
            if tifp:
                fallback_names = _extract_names_from_tif(tifp)
        if (not fallback_names) and qupath_root:
            sidecars = _find_qupath_sidecars(ome_stem, qupath_root)
            if sidecars:
                fallback_names = _extract_names_from_qupath(sidecars)
    for image in root.xpath("//ome:Image", namespaces=ns):
        pixels = image.find("{*}Pixels")
        if pixels is None:
            continue
        channels = [c for c in pixels.findall("{*}Channel")]
        hex_palette = make_palette(len(channels))
        for idx, ch in enumerate(channels):
            # Ensure Name
            name = ch.get("Name") or ""
            is_generic = bool(re.match(r"^Channel[ _]?\d+$", name.strip()))
            if not name.strip() or is_generic:
                # Try annotations first, then fallback lists, then ID inference
                inferred = ann_names.get(idx)
                if (not inferred) and fallback_names and idx < len(fallback_names):
                    inferred = fallback_names[idx]
                if not inferred:
                    inferred = infer_name_from_id(ch.get("ID") or "")
                # Avoid setting to pure digits
                if inferred.isdigit() or not inferred:
                    inferred = f"Channel {idx+1}"
                ch.set("Name", inferred)
            # Ensure Color
            cur = ch.get("Color")
            needs = True
            try:
                # Treat missing, zero, or negative as needing assignment
                needs = (cur is None) or (int(cur) <= 0)
            except Exception:
                needs = True
            # Always enforce preferred colors for these four names
            nm = (ch.get("Name") or names_map.get(idx) or "").strip().upper()
            enforce = PREFERRED_HEX.get(nm)
            if enforce:
                rgb_dec = ensure_rgb_decimal(enforce)
                if rgb_dec:
                    ch.set("Color", rgb_dec)
            elif needs:
                # Use annotation hex if provided, else palette
                hex_color = ann_colors.get(idx) or hex_palette[idx]
                rgb_dec = ensure_rgb_decimal(hex_color)
                if rgb_dec:
                    ch.set("Color", rgb_dec)
    return etree.tostring(root, pretty_print=True, xml_declaration=True, encoding="UTF-8")


def process_file(path: str, dry_run=False, fallback_root: str = None, qupath_root: str = None, names_map: dict | None = None):
    info = tt.read_tiff(path)
    ifd, tag = read_ome_xml(info)
    if ifd is None:
        print(f"SKIP (no OME-XML): {path}")
        return
    try:
        base = os.path.basename(path)
        stem = re.sub(r"\.ome\.tif(f)?$", "", base, flags=re.I)
        new_xml_bytes = patch_names_and_colors(tag["data"], ome_stem=stem, fallback_root=fallback_root, qupath_root=qupath_root, names_map=names_map)  # type: ignore[index]
    except Exception as e:
        print(f"ERROR patching {path}: {e}", file=sys.stderr)
        return
    if dry_run:
        print(f"DRY-RUN would set colors for {path}")
        return
    # Update all IFDs that contain OME-XML in ImageDescription
    new_xml_str = new_xml_bytes.decode("utf-8", errors="ignore")
    TAG_IMAGE_DESCRIPTION = 270
    for ifd in info.get("ifds", []):
        tg = ifd.get("tags", {}).get(TAG_IMAGE_DESCRIPTION)
        if not tg:
            continue
        val = tg.get("data")
        s = val if isinstance(val, str) else (val or b"").decode("utf-8", errors="ignore")
        if ("OME" in s) and ("Pixels" in s):
            tg["data"] = new_xml_str
    # tifftools.write_tiff expects a list of IFD dicts, not the whole info dict
    tt.write_tiff(path, info.get("ifds", []))
    print(f"UPDATED: {path}")


def main():
    args = parse_args()
    root = os.path.abspath(args.data_root)
    # Load optional channel names file
    names_map = load_channel_names_file(args.channel_names_file)
    patterns = args.glob.split(":") if ":" in args.glob else [args.glob]
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(root, "**", pat), recursive=True))
    files = [f for f in files if os.path.isfile(f)]
    if not files:
        print("No files found.")
        return
    if args.report_only:
        for f in sorted(files):
            ti = tt.read_tiff(f)
            ifd, tag = read_ome_xml(ti)
            if ifd is None:
                print(f"SKIP (no OME-XML): {f}")
                continue
            _, meta = parse_channel_info(tag["data"])  # type: ignore[index]
            print(f"FILE: {f}")
            for ch in meta:           
                print(f"  Image={ch['image_name']} idx={ch['index']} name={ch['channel_name']} color={ch['color']}")
    else:
        for f in sorted(files):
            process_file(f, dry_run=args.dry_run, fallback_root=args.fallback_root, qupath_root=args.qupath_root, names_map=names_map)


if __name__ == "__main__":
    main()
